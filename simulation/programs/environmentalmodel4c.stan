/*Program by Luis Leon Novelo
Sep 30 2020
The innovation is 
sigma_alpha~ p_{equa}\sim palphasequal*1(\sigma_\alpha=c_sigma_alpha)
            +(1-p_{equal})* 
             Inverse-gamma(sigma_alpha|shape_\sigma_\alpha,
                           scale_\sigma_\alpha,
                           truncated to be >c_sigma_alpha)

*/

functions{
    real ftem(real xx){
    real out;
    out=0.1;
    if(xx>0.1){out=xx;}
    return (out);}




//beta1=sapply(x1+y1+z1,ftem)/            sum(c(q1[1],q1[2]*alpha_est,((1-q1[1])+(1-q1[2])*alpha_est)*tau))
//gammatimesbeta2=sapply((x2+y2+z2),ftem)/sum(c(q2[1],q2[2]*delta_est,((1-q2[1])+(1-q2[2])*delta_est)*tau))

}
  



data{
   int <lower=2> K;             //Total number of bioresps
	 int n_environment;          //Number of environments. So far 2 mated and virgin
	 int xenv [K];               //Environment index mated or virgin
	 //int xbiorep[K];           // Biorep index
	 int xs [K];                 //tester counts     
	 int ys [K];                 //lin counts
	 int zs [K];                 //both counts
	 matrix [n_environment,2] r;             //Matrix of bias correction  
	 real <lower=0> a_beta;               //beta_1,...beta_K~gamma(a_beta,b_beta)
	 real a_b_beta;         //b_beta~gamma(a_b_beta,b_b_beta)
	 real b_b_beta;
	 real  a_overdispersion;
	 real  b_overdispersion;//phi is apriori small, if inverse gamma is used prior mean b_phi/(a_phi-1)	 
	 real c_sigma_alpha;     // value of c_sigma_alpha, it depends on the number of environments tested; 
	 }


parameters{
  real  <lower=0> overdispersion;
  real  <lower=0> b_beta;
  real mu_alpha;
  real <lower=c_sigma_alpha> sigma_alpha; // distribution of sigma alpha is truncated
  real <lower=0,upper=1> palphasequal;    // Probability of equal Allelic Imbalance
  vector <lower=0> [K] bbeta;             // Biorep specific
  vector <lower=0> [n_environment] alpha;
    }

transformed parameters{
}


model{
  
		for(k in 1:K) {
			   		xs[k]~ neg_binomial_2((1/alpha[xenv[k]])*r[xenv[k],1]*bbeta[k], 1/overdispersion);
			   		ys[k]~ neg_binomial_2(   alpha[xenv[k]] *r[xenv[k],2]*bbeta[k], 1/overdispersion);
			   		zs[k]~ neg_binomial_2(   (
			   		                         (1-r[xenv[k],1])*(1/alpha[xenv[k]])+
			   		                         (1-r[xenv[k],2])*    alpha[xenv[k]]
			   		                         )*bbeta[k], 1/overdispersion);
			   		}
		bbeta~    gamma(a_beta,b_beta);
		b_beta~   gamma(a_b_beta,2*b_b_beta);
		
		target += log_mix(palphasequal,
                     lognormal_lpdf(alpha | mu_alpha,  c_sigma_alpha),
                     lognormal_lpdf(alpha | mu_alpha, sigma_alpha));
                     

		mu_alpha~ normal(0,2);
		overdispersion ~  inv_gamma(a_overdispersion, b_overdispersion);
		sigma_alpha ~      normal(0,1);
    palphasequal ~      beta(19,1); 
}
	
	generated quantities{
	  vector  [n_environment] theta; //Proportion of reads coming from tester allele
	  for(i in 1:n_environment){
	    theta[i]=1/(alpha[i]^2+1);
	    }
	 	}
