what.to.do <- (commandArgs(TRUE))
if(length(what.to.do)>0)
{
  my.alpha <- as.numeric(as.character(unlist(strsplit(what.to.do[1],"="))[2]))
	# my.delta <- as.numeric(as.character(unlist(strsplit(what.to.do[2],"="))[2]))
	my.qtest <- as.numeric(as.character(unlist(strsplit(what.to.do[2],"="))[2]))
	my.qline <- as.numeric(as.character(unlist(strsplit(what.to.do[3],"="))[2]))
	mult <- as.numeric(as.character(unlist(strsplit(what.to.do[4],"="))[2]))
	myreads <- as.numeric(as.character(unlist(strsplit(what.to.do[5],"="))[2]))
	simruns <- as.numeric(as.character(unlist(strsplit(what.to.do[6],"="))[2]))
	out.file <- as.character(as.character(unlist(strsplit(what.to.do[7],"="))[2]))
  nconditions<- as.numeric(as.character(unlist(strsplit(what.to.do[8],"="))[2]))

}

#It is horrible I know!
#Temporarily we hardcode the number of replicates
nreps <- 3

#Removed set.seed(0) because always gave the same result!!!!!
#set.seed(0)

#Delta is like alpha, but in the second environment (i.e. delta=1 is no AI in ENV2)


###################################################
#
# Try to generalize for n environments
#
###################################################

print(nconditions)
seqcond <- paste("c", seq(1,nconditions), sep="")


repsuffix <- paste0("_rep", seq(1, nreps))
# fornames1 <- paste("g1_", seqcond, sep="") #Genotype 1 given condition number
# fornames2 <- c("sampleprop", "mean", "q025", "q975", "Bayes_evidence", "AI_decision")
# fornames <- paste(paste(rep(fornames1, each=length(fornames2)), fornames2, sep="_"), collapse=",")

#The three possible alignment states: better to tester (G1?), better to line (G2?), equally good (both).
#alnto<-c("g1","g2","Both")
counttester <- as.vector(sapply(paste("counts_", seqcond, "_g1", sep=""), paste0, repsuffix))
countline <- as.vector(sapply(paste("counts_", seqcond, "_g2", sep=""), paste0, repsuffix))
countboth <- as.vector(sapply(paste("counts_", seqcond, "_both", sep=""), paste0, repsuffix))
countheader <- paste(counttester, countline, countboth, sep=",", collapse=",")
priorstester <- paste("prior_", seqcond, "_g1", sep="")
priorsline <- paste("prior_", seqcond, "_g2", sep="")
allpriors <- sort(c(priorstester, priorsline))
seqreps <- paste(seqcond, "_num_reps", sep="")
activeflag <- paste(seqcond, "flag_analyze", sep="_")
headers <- paste(c("comparison", "FEATURE_ID", seqreps, countheader, allpriors, activeflag), collapse=",")
print(unlist(strsplit(headers,",")))
cat(headers, file=out.file, append=FALSE, sep="\n")


#End of arguments and names copied from the scripts used for the analysis of real data



##################################
#
# Simulation
#
##################################
print(simruns)
for (aaa in 1:simruns)
{
  #It simulates data according to a specific model (I need to check wich one)
  #and compares the inference for different versions of the model
  #
  #Parameters of true model
  allI <- rep(nreps, nconditions) #Number of replicates in each environment (we are currentyl stuck to having the same number of reps in the environments)
  # true_betas_row1 <- c(1, 1.2, 0.7)*myreads   #These are the biorep effets, the beta_i s
  # #At present we use true_betas that do not vary across conditions. This may be changed in the future
  # #When we were fixed to 2 environments we used different betas
  # true_betas_row2 <- c(1,1.3,0.8)*myreads
  true_betas_row1 <- c(1, 1, 1)*myreads
  true_alpha <- my.alpha            #if different from 1 AI at environment 1 (Mated)
  # true_delta <- my.delta            #if different from 1 AI at environment 1 (Virgin)
  true_gamma <- 1            #The current model assumes 1, but it used to be useful before
  true_tau <- 1              #The current model assumes 1, but at some point we considered it different from 1

  #If I well understand the first is q in tester and the second in line. q_row1 is enviornment1 and q_row2 environment 2
  
  allq <- matrix(rep(c(my.qtest, my.qline), nconditions), ncol=2, byrow=T)  
  #q_row1 <- c(my.qtest, my.qline)
  #q_row2 <- c(my.qtest, my.qline)
  
  true_phi <- 0.02           #Neg binomial dispersion parameter, the bigger it is the higher the variance=mu+phi mu^2 is

  #xs <- ys <- zs <- rep(NA, sum(allI))

  flaganalyze <- rep(1, nconditions)
  for(loopc in 1:nconditions)
	{
		cat("loopc is", loopc, "\n")
		means <- c(allq[loopc,1]/true_alpha, allq[loopc,2]*true_alpha, mult*((1-allq[loopc,1])/true_alpha+(1-allq[loopc,2])*true_alpha)*true_tau)
		for(i in 1:allI[loopc])
		{
			xs <- rnbinom(1, size=1/true_phi, mu=means[1]*true_betas_row1[i])
			ys <- rnbinom(1, size=1/true_phi, mu=means[2]*true_betas_row1[i])
			zs <- rnbinom(1, size=1/true_phi, mu=means[3]*true_betas_row1[i])
			# cat("xs of i",i,"is",xs[i],"\n")
			# cat("true_betas_row1[i] is",true_betas_row1[i],"\n")
			# cat("means[1] is",means[1],"\n")
		  if(loopc==1 & i==1) mycounts <- paste(xs, ys, zs, sep=",") else mycounts <- paste(mycounts, xs, ys, zs, sep=",")
		}	
	}

	cat("mycounts is", mycounts, "\n")
	
	out <- paste("line", "fusion_id", paste(allI, collapse=","),
	          mycounts, 
	          paste(as.vector(t(allq)), collapse=","), paste(flaganalyze, collapse=","), sep=",")
	cat(out, file=out.file, append=TRUE, sep="\n")
}
