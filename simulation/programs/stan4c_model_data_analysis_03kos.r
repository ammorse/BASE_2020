what.to.do <- (commandArgs(TRUE))
if(length(what.to.do)>0)
{
  data.file <- as.character(as.character(unlist(strsplit(what.to.do[1],"="))[2]))
  q.model.g1 <- as.numeric(as.character(unlist(strsplit(what.to.do[2],"="))[2]))
  q.model.g2 <- as.numeric(unlist(strsplit(what.to.do[3],"="))[2])
  stan.program <- as.character(unlist(strsplit(what.to.do[4],"="))[2])
  niter <- as.numeric(as.character(unlist(strsplit(what.to.do[5],"="))[2]))
  nburnin <- as.numeric(as.character(unlist(strsplit(what.to.do[6],"="))[2]))
  save.sample.freq <- as.numeric(as.character(unlist(strsplit(what.to.do[7],"="))[2]))
  reporting.freq <- as.numeric(as.character(unlist(strsplit(what.to.do[8],"="))[2]))
  stepsize <- as.numeric(as.character(unlist(strsplit(what.to.do[9],"="))[2]))
  nconditions<- as.numeric(as.character(unlist(strsplit(what.to.do[10],"="))[2]))
}

out.dir <- dirname(data.file)
print(out.dir)

base.out<-paste(gsub(".csv","",basename(data.file)), 
                                 stan.program,
                                 'qmodel-g1', q.model.g1, 'qmodel-g2', q.model.g2, 
                                 'iterations', niter, 'burnin', nburnin, 
                                 'thin', save.sample.freq, 'refresh', reporting.freq, 
                                 'stepsize', stepsize, sep="_")
fileout <- paste(paste(out.dir, base.out, sep = "/"),".csv",sep="")
print(fileout)

#########################
#Trying to generalize to an arbitrary number of conditions
#########################

seqcond<-paste("c",seq(1,nconditions),sep="")
fornames1<-paste("g1_",seqcond,sep="") #Tester given condition
fornames2 <- c("sampleprop", "theta", "q025", "q975", "Bayes_evidence", "AI_decision")
fornames<-paste(paste(rep(fornames1,each=length(fornames2)),fornames2,sep="_"),collapse=",")
#The three possible alignment states: better to tester (G1?), better to line (G2?), equally good (both).
#alnto<-c("g1","g2","Both")
counttester<-paste("counts_",seqcond,"_g1",sep="")
countline<-paste("counts_",seqcond,"_g2",sep="")
countboth<-paste("counts_",seqcond,"_both",sep="")
countheader<-paste(counttester,countline,countboth,sep=",",collapse=",")
priorstester<-paste("prior_",seqcond,"_g1",sep="")
priorsline<-paste("prior_",seqcond,"_g2",sep="")
allpriors<-sort(c(priorstester,priorsline))
seqreps<-paste(seqcond,"_num_reps",sep="")
activeflag<-paste(seqcond,"flag_analyze",sep="_")
firstheaders<-paste(c("comparison","FEATURE_ID", seqreps, countheader, allpriors, "H3_independence_Bayes_evidence", 
                      "probofsameAI"), collapse=",")

alpha_post_names<-paste(paste("alpha",seq(1,nconditions),"_postmean",sep=""),collapse=",")
cat("alpha names are:",alpha_post_names,"\n")
cat("stan4 firstheaders are:",firstheaders,"\n")
headers_out=paste(firstheaders,fornames,alpha_post_names,"sigma_alpha_mean","sigma_alpha_0.5","sigma_alpha_0.95",paste(activeflag,collapse=","),sep=",")
cat("stan4 headers_out is:",headers_out,"\n")
cat(headers_out,file=fileout,append=FALSE,sep="\n")


#End of arguments and names copied from the scripts used for the analysis of real data

#The special thing about this function is that it is able to accept a multiplicative factor for the estimation of the reads aligning on both genomes, to see how misspecification of that quantity affects results
options(warn=1)

library("rstan")
rstan_options(auto_write = FALSE)
gam.mles.data <- function(x){
  # CALCULATION OF THE MLES for gamma  GIVEN THE SAMPLE
  n <- length(x)
  xb <- mean(x)
  xd <- mean(log(x))
  s <- log(xb)-xd
  a0 <- (3.0-s+sqrt((s-3.0)^2+24.0*s))/12.0/s
  l <- 1
  repeat{
    ans <- (log(a0)-digamma(a0)-s)
    a1 <- a0-ans/(1.0/a0-trigamma(a0))
    if(abs(ans) <= 1.0e-7 | l >= 30){break}
    a0 <- a1
    l <- l+1}
  ah <- a1; bh <- xb/a1
  return(c(ah,bh))
}

#Computing prior hyperparameters for beta, the size effect
prior_empBayes_forbeta <- function(xs,ys,zs){
  #Function to compute the hyperparameters of the gamma distribution for
  #beta_1,...beta_K~beta(a_beta,b_beta) with b_beta~gamma(a_b_beta,b_b_beta)
  bbeta_est <- (xs+ys+zs)/2
  bbeta_est[which(bbeta_est==0)] <- 0.1
  tem <- gam.mles.data(bbeta_est) #MLE of the gamma function but it is parameterized so that E(x)=ab if x~gamma(a,b)
  tem[1] <- min(max(tem[1],10^(-3)),10^5)
  tem[2] <- min(max(tem[2],10^(-3)),10^5)  #Our gamma parameterization is E(x)=a/b
  a_beta <- tem[1]
  a_b_beta <- 2*tem[2]^(-1)#
  b_b_beta <- 2
  return(list(a_beta=a_beta, a_b_beta=a_b_beta, b_b_beta=b_b_beta))
}

library(tidyverse)
mydata <- as_tibble(read.csv(data.file))

for (gene in 1:nrow(mydata))
{
  
  nreps <- mydata[gene,] %>% dplyr::select(ends_with("num_reps"))
  nconditions <- length(nreps)            #Number of environments
  seqI<-nreps
  allI<-as.vector(seqI)  
  cat("seqI is")
  print(seqI)
 
  xs <- unlist(mydata[gene,] %>% dplyr::select(starts_with("counts")) %>% dplyr::select(contains("g1")))
  ys <- unlist(mydata[gene,] %>% dplyr::select(starts_with("counts")) %>% dplyr::select(contains("g2")))
  zs <- unlist(mydata[gene,] %>% dplyr::select(starts_with("counts")) %>% dplyr::select(contains("both")))
  xactiveflag <- unlist(mydata[gene,] %>% dplyr::select(contains("flag_analyze")))
  print(xs)
  print(ys)
  print(zs)
  #If at least one flaganalyze is 0 we fill the results wiht NA and go to the next
  if(prod(xactiveflag)==0)
  {
    howmanyNA<-length(unlist(strsplit(headers_out,",")))-2
	out=paste("line", "fusion_id", rep(NA,howmanyNA,collapse=","),sep=",")
    cat(out,file=fileout,append=TRUE,sep="\n")
    next
  }
  q_sim <- matrix(unlist((mydata[gene, ]%>% dplyr::select(starts_with("prior")))), nrow=length(nreps), byrow =TRUE, ncol=2)
  
  hyper_beta <- prior_empBayes_forbeta(xs, ys, zs)
  #Making the data ready for stan
  datastan <- list(K = sum(seqI),                 #Total number of bioresps
	    n_environment = nconditions,         #Number of environments, so far 2.
	    xenv = rep(seq(1,length(seqI)),seqI),       #Environment index 
	    xs = as.vector(xs),
	    ys = as.vector(ys),
	    zs = as.vector(zs),
	    r = q_sim,  #matrix of systematic bias corrections
	    a_beta = hyper_beta$a_beta,               #Set to MLE under the null model
	    a_b_beta = hyper_beta$a_b_beta,           #b_beta~gamma(a_b_beta,b_b_beta)
	    b_b_beta = hyper_beta$b_b_beta,
	    a_overdispersion = 2.01,#1,
	    b_overdispersion = 0.05 #100, #phi is apriori small, if inverse gamma is used prior mean b_phi/(a_phi-1)
	 )
  
  starting_values  <-  function(){ ###Not actually needed according to Luis
    out <- with(datastan, list(overdispersion=0.01, bbeta=xs+ys+zs, alpha=rep(1.0,n_environment)))
    return(out)
  }

	totalcounts=
	rbind(
	xs=tapply(datastan$xs,FUN=sum,INDEX=datastan$xenv),
	ys=tapply(datastan$ys,FUN=sum,INDEX=datastan$xenv),
	zs=tapply(datastan$zs,FUN=sum,INDEX=datastan$xenv)
	)
  
  datastan4c <- datastan
  datastan4c$c_sigma_alpha <- c(NA,0.01760,0.01472,0.01343,0.01265 )[datastan4c$n_environment]#Constants depends on number of environments
  cat("datastan is")
  print(datastan4c)

  fit1 <-rstan::stan(
    file = "environmentalmodel4c.stan", # Stan program
	  data = datastan4c,    # named list of data
	  chains = 1,             # number of Markov chains
	  warmup = nburnin,          # number of warmup iterations per chain
	  thin = save.sample.freq,
	  iter = niter,            # total number of iterations per chain
	  refresh = reporting.freq,             # no progress shown
	  init = "starting_values",
	  control = list(adapt_delta = stepsize)  # The bigger the value the smaller the step size
	  #pars=c("alpha","sigma_alpha")
	  #c("bbeta","alpha","overdispersion","theta","sigma_alpha") 
   )
	
   # fig.dir <- paste0(unlist(strsplit(out.dir, "/g3_sim_output"))[1], "/data_visualization")
   # jpeg(paste0(fig.dir,"/diagnostic_plots_stan4c_", gsub(".csv", "", tail(unlist(strsplit(fileout, "/")), 1)),".jpeg"))
   # rstan::pairs(fit1)
   # dev.off()
  
	palphasequal <- rstan::extract(fit1,pars=c("palphasequal"))$palphasequal
	theta <- rstan::extract(fit1,pars="theta")$theta #Estimated proportion of reads aligning to tester in mated after adjusting for systematic bias
	alpha <- rstan::extract(fit1,pars=c("alpha"))$alpha
	sigma_alpha <- rstan::extract(fit1,pars=c("sigma_alpha"))$sigma_alpha
	
	Stanresults <- matrix(NA,nrow=nconditions,ncol=4)
	Bayes_AI_pvalue <- rep(NA,nconditions)
	for(mystan in 1:nrow(Stanresults)) {
	  theta1 <- theta[,mystan]
	  alpha1 <- alpha[,mystan]
    Stanresults[mystan,]<-c(mean(theta1), quantile(theta1,c(0.025,0.975)),2*min(mean(theta1>1/2),mean(theta1<1/2))) 
	Bayes_AI_pvalue[mystan] <- 2*min(c(mean(alpha1>1), mean(alpha1<1)))
	} 
	
	colnames(Stanresults) <- c("mean", "q_025", "q_975", "AIbayesian-pval")
	cat("Stanresults",Stanresults,"\n")
	Stanresults[,"AIbayesian-pval"]

	sigma_alpha_out <- c(mean(sigma_alpha), quantile(sigma_alpha, c(0.5,0.95))) #POSSIBLE bug: c(0.5,0.95) or c(0.05,0.95)? 

	for(repcond in 1:nconditions) {
	  theta1 <- theta[,repcond]
	  smallStanres <- paste(round(totalcounts[row.names(totalcounts)=="xs", repcond]/(totalcounts[row.names(totalcounts)=="xs", repcond] + totalcounts[row.names(totalcounts)=="ys", repcond]), 4),
	                        round(mean(theta1), 4),
	                        paste(round(quantile(theta1, c(0.025,0.975)), 4), collapse=","),
	                        round(Bayes_AI_pvalue[repcond], 4), 
	                        ifelse(Bayes_AI_pvalue[repcond] < 0.05, 1, 0), sep=",")
	  if(repcond==1) fullStanres <- smallStanres else fullStanres <- paste(fullStanres, smallStanres, sep=",", collapse=",")
	}
	
	alpha1greateralpha2 <- NA
	#alpha1greateralpha2 is the bayesian independence test. Only meaningful for two conditions
	if(nconditions==2) {
	  alpha1greateralpha2 <- min(tem<-mean(alpha[,1]-alpha[,2]<0),1-tem)*2
	}

	probofsameAI <- mean(palphasequal)    
	
	totalcounts <- rbind(xs <- tapply(datastan$xs, FUN=sum, INDEX=datastan$xenv),
	                     ys <- tapply(datastan$ys, FUN=sum, INDEX=datastan$xenv),
	                     zs <- tapply(datastan$zs, FUN=sum, INDEX=datastan$xenv))
	
	out <- paste("line", "fusion_id", paste(allI,collapse=","),
	             paste(apply(totalcounts, 2, paste, collapse=","), collapse=","),
				 paste(t(datastan$r),collapse=","),
	             alpha1greateralpha2,
	             probofsameAI,
	             fullStanres,
	             paste(round(apply(alpha,2,mean),4),collapse=","),
	             paste(round(sigma_alpha_out,6),collapse=","),
				 paste(xactiveflag,collapse=","),
	             sep=",")
    print(xactiveflag)
	cat(out, file=fileout, append=TRUE, sep="\n")
    print(unlist(strsplit(headers_out,",")))
    print(unlist(strsplit(out,",")))

}
