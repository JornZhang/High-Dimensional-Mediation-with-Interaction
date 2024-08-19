simnum <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))




library(glmnet)
library(ncvreg)
library(freebird) # Dave method

source("MyMethod_functions.R")

set.seed(12301)



run=10

n=300
p=500
#p=1000

q=1


one.n=as.matrix(rep(1, times=n))

TRT=matrix( c(rep(1,times=(n/2)), rep(0, times=(n/2))) , nrow = n, ncol = q)
colnames(TRT) = "TRT"

theta0=0.5

c.T=0.75*4.2
theta.T=c.T*1


c.G = 1
s.G = 10
theta.G=c.G*matrix( c( c(1, 0.9, 0.8, 0.7, 0.6), 
                       c(1, 0.9, 0.8, 0.7, 0.6),
                       rep(0, p - s.G) ), 
                    ncol = 1 )



cstar.seq = c(-2,  
              -1.5, 
              -1.08, -1.06, -1.04, -1.02, -1, -0.98, -0.96, -0.94, -0.92,
              -0.5,
              0,
              0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.05, 1.1, 1.15,
              1.2, 1.5, 2)
s.D=10 # s.D is the sparsity of each column in Theta3




sigma_Y=0.5


# for bebta0 and A, I didn't put any sparsity in it. But Dave did it in the simulation.
c_beta0=1
#beta0=c_beta0 * matrix(runif(p, -1, 1), ncol=1)
#beta0=c_beta0 * matrix( c(0.2*c(1:5), rnorm(p-5, sd=0.1)) , ncol=1)
#beta0=c_beta0 * matrix( c(0.2*c(1:5), runif(p-5, -1,1)) , ncol=1)
#beta0=c_beta0 * matrix( c( c(-0.3, -0.2, 0.2, 0.3, 0.1), rnorm(p-5, sd=0.1)) , ncol=1)
beta0=c_beta0 * matrix( c( (-0.1)*c(1:10), 
                           rnorm(p-10, sd=0.1) ) , ncol=1)



c_B = 1
#B=c_B * matrix(runif(p*q, -1, 1), nrow=p, ncol = q)
#B=c_B * matrix( c( 0.2*c(1:5), 0.2*c(1:5), 0.2*c(1:5), 0.2*c(1:5), rnorm(p-20, sd=0.1)), nrow = p, ncol = q )
#B=c_B * matrix( c( 0.2*c(1:5), 0.2*c(1:5), 0.2*c(1:5), 0.2*c(1:5), runif(p-20, -1,1)), nrow = p, ncol = q )
#B=c_B * matrix( c( c(-0.3, -0.2, 0.2, 0.3, 0.1), rnorm(p-5, sd=0.1)), nrow = p, ncol = q )
B=c_B * matrix( c( 0.1*c(1:10) , 
                   rnorm(p-10, sd=0.1) ), 
                nrow = p, ncol = q )


rho_G=0.5 # if rho_G=0, Sigma_G=I_p * c_G
c_sigmaG=1 # c_G is the magnitude of the noise

# block diagonal
Gamma_G=matrix(0, p, p)
noc=1
blk.size=p/noc
diag.blk=matrix(rho_G, blk.size, blk.size) + (1-rho_G)*diag(blk.size)
for (i in 1:noc) {
  Gamma_G[((i-1)*blk.size+1):(i*blk.size), ((i-1)*blk.size+1):(i*blk.size)]=diag.blk
}
rm(diag.blk)


Sigma_G=c_sigmaG*Gamma_G
Q_G=eigen(Sigma_G)$vectors
Lambda0.5_G=diag(c(sqrt(eigen(Sigma_G)$values)))
Sigma0.5_G=Q_G%*%Lambda0.5_G%*%t(Q_G)
# check whether it is correct
print(max(abs(Sigma0.5_G%*%Sigma0.5_G-Sigma_G)))
rm(Q_G, Lambda0.5_G, Gamma_G)





p.fac<- rep(1, (q+p+q*p))
p.fac[1:(q)] <- 0


p.fac_runze = rep(1, (q+p))
p.fac_runze[1:(q)] = 0

p.fac_lasso = p.fac

##### SCAD
bias_theta.scad=matrix(0, nrow = run, ncol = length(cstar.seq))
sr_theta.scad=matrix(0, nrow = run, ncol = length(cstar.seq))

est_drct.10.scad=matrix(0, nrow = run, ncol = length(cstar.seq))
est_drct.11.scad=matrix(0, nrow = run, ncol = length(cstar.seq))
SE_drct.10.scad=matrix(0, nrow = run, ncol = length(cstar.seq))
SE_drct.11.scad=matrix(0, nrow = run, ncol = length(cstar.seq))
cp_drct.10.scad=matrix(0, nrow = run, ncol = length(cstar.seq))
cp_drct.11.scad=matrix(0, nrow = run, ncol = length(cstar.seq))
reject_drct.10.scad=matrix(0, nrow = run, ncol = length(cstar.seq))
reject_drct.11.scad=matrix(0, nrow = run, ncol = length(cstar.seq))

est_indrct.10.scad=matrix(0, nrow = run, ncol = length(cstar.seq))
est_indrct.11.scad=matrix(0, nrow = run, ncol = length(cstar.seq))
SE_indrct.10.scad=matrix(0, nrow = run, ncol = length(cstar.seq))
SE_indrct.11.scad=matrix(0, nrow = run, ncol = length(cstar.seq))
cp_indrct.10.scad=matrix(0, nrow = run, ncol = length(cstar.seq))
cp_indrct.11.scad=matrix(0, nrow = run, ncol = length(cstar.seq))
reject_indrct.10.scad=matrix(0, nrow = run, ncol = length(cstar.seq))
reject_indrct.11.scad=matrix(0, nrow = run, ncol = length(cstar.seq))


##### Oracle
est_drct.10.oracle=matrix(0, nrow = run, ncol = length(cstar.seq))
est_drct.11.oracle=matrix(0, nrow = run, ncol = length(cstar.seq))
SE_drct.10.oracle=matrix(0, nrow = run, ncol = length(cstar.seq))
SE_drct.11.oracle=matrix(0, nrow = run, ncol = length(cstar.seq))
cp_drct.10.oracle=matrix(0, nrow = run, ncol = length(cstar.seq))
cp_drct.11.oracle=matrix(0, nrow = run, ncol = length(cstar.seq))
reject_drct.10.oracle=matrix(0, nrow = run, ncol = length(cstar.seq))
reject_drct.11.oracle=matrix(0, nrow = run, ncol = length(cstar.seq))

est_indrct.10.oracle=matrix(0, nrow = run, ncol = length(cstar.seq))
est_indrct.11.oracle=matrix(0, nrow = run, ncol = length(cstar.seq))
SE_indrct.10.oracle=matrix(0, nrow = run, ncol = length(cstar.seq))
SE_indrct.11.oracle=matrix(0, nrow = run, ncol = length(cstar.seq))
cp_indrct.10.oracle=matrix(0, nrow = run, ncol = length(cstar.seq))
cp_indrct.11.oracle=matrix(0, nrow = run, ncol = length(cstar.seq))
reject_indrct.10.oracle=matrix(0, nrow = run, ncol = length(cstar.seq))
reject_indrct.11.oracle=matrix(0, nrow = run, ncol = length(cstar.seq))






##### Runze
est_drct.runze=matrix(0, nrow = run, ncol = length(cstar.seq))
SE_drct.runze=matrix(0, nrow = run, ncol = length(cstar.seq))
cp_drct10.runze=matrix(0, nrow = run, ncol = length(cstar.seq))
cp_drct11.runze=matrix(0, nrow = run, ncol = length(cstar.seq))
reject_drct.runze=matrix(0, nrow = run, ncol = length(cstar.seq))

est_indrct.runze=matrix(0, nrow = run, ncol = length(cstar.seq))
SE_indrct.runze=matrix(0, nrow = run, ncol = length(cstar.seq))
cp_indrct10.runze=matrix(0, nrow = run, ncol = length(cstar.seq))
cp_indrct11.runze=matrix(0, nrow = run, ncol = length(cstar.seq))
reject_indrct.runze=matrix(0, nrow = run, ncol = length(cstar.seq))







##### Lasso
bias_theta.lasso=matrix(0, nrow = run, ncol = length(cstar.seq))
sr_theta.lasso=matrix(0, nrow = run, ncol = length(cstar.seq))

est_drct.10.lasso=matrix(0, nrow = run, ncol = length(cstar.seq))
est_drct.11.lasso=matrix(0, nrow = run, ncol = length(cstar.seq))
bias_drct.10.lasso=matrix(0, nrow = run, ncol = length(cstar.seq))
bias_drct.11.lasso=matrix(0, nrow = run, ncol = length(cstar.seq))

est_indrct.10.lasso=matrix(0, nrow = run, ncol = length(cstar.seq))
est_indrct.11.lasso=matrix(0, nrow = run, ncol = length(cstar.seq))
bias_indrct.10.lasso=matrix(0, nrow = run, ncol = length(cstar.seq))
bias_indrct.11.lasso=matrix(0, nrow = run, ncol = length(cstar.seq))




##### Dave
est_drct.dave.cY=matrix(0, nrow = run, ncol = length(cstar.seq))
SE_drct.dave.cY=matrix(0, nrow = run, ncol = length(cstar.seq))
cp_drct10.dave.cY=matrix(0, nrow = run, ncol = length(cstar.seq))
cp_drct11.dave.cY=matrix(0, nrow = run, ncol = length(cstar.seq))
reject_drct.dave.cY=matrix(0, nrow = run, ncol = length(cstar.seq))

est_indrct.dave.cY=matrix(0, nrow = run, ncol = length(cstar.seq))
SE_indrct.dave.cY=matrix(0, nrow = run, ncol = length(cstar.seq))
cp_indrct10.dave.cY=matrix(0, nrow = run, ncol = length(cstar.seq))
cp_indrct11.dave.cY=matrix(0, nrow = run, ncol = length(cstar.seq))
reject_indrct.dave.cY=matrix(0, nrow = run, ncol = length(cstar.seq))



### Huang
est_indrct.10.Huang=matrix(0, nrow = run, ncol = length(cstar.seq))
est_indrct.11.Huang=matrix(0, nrow = run, ncol = length(cstar.seq))

SE_indrct.10.Huang=matrix(0, nrow = run, ncol = length(cstar.seq))
SE_indrct.11.Huang=matrix(0, nrow = run, ncol = length(cstar.seq))
cp_indrct.10.Huang=matrix(0, nrow = run, ncol = length(cstar.seq))
cp_indrct.11.Huang=matrix(0, nrow = run, ncol = length(cstar.seq))
reject_indrct.10.Huang=matrix(0, nrow = run, ncol = length(cstar.seq))
reject_indrct.11.Huang=matrix(0, nrow = run, ncol = length(cstar.seq))

SE_indrct.10.Huang.boot=matrix(0, nrow = run, ncol = length(cstar.seq))
SE_indrct.11.Huang.boot=matrix(0, nrow = run, ncol = length(cstar.seq))
cp_indrct.10.Huang.boot=matrix(0, nrow = run, ncol = length(cstar.seq))
cp_indrct.11.Huang.boot=matrix(0, nrow = run, ncol = length(cstar.seq))
reject_indrct.10.Huang.boot=matrix(0, nrow = run, ncol = length(cstar.seq))
reject_indrct.11.Huang.boot=matrix(0, nrow = run, ncol = length(cstar.seq))

PCAbias_indrct.10.Huang=matrix(0, nrow = run, ncol = length(cstar.seq))
PCAbias_indrct.11.Huang=matrix(0, nrow = run, ncol = length(cstar.seq))







#### set new seed for different arrays ####
set.seed(1234*simnum)




for (loop.c in 1:length(cstar.seq) ) {
  
  c.D = cstar.seq[loop.c]
  
  Theta.D=c.D*matrix( c( c(1, 0.9, 0.8, 0.7, 0.6), 
                         c(1, 0.9, 0.8, 0.7, 0.6),
                         rep(0, p - s.D) ), 
                      ncol = q )
  
  theta=rbind(theta0, theta.T, theta.G, Theta.D)
  
  NDE0 = as.numeric( theta.T + t(beta0)%*%Theta.D )
  NDE1 = as.numeric( theta.T + t(beta0 + B)%*%Theta.D  )
  
  for (loop.run in 1:run) {
    
    # generate mediation and response
    
    White_G=matrix(rnorm(n*p), nrow = n, ncol = p)
    #White_G=matrix(runif(n*p, min = -sqrt(3), max = sqrt(3)), n,p)
    #White_G=matrix((rpois(n*p, lambda = 1)-1), n, p)
    #White_G=matrix((2*rbinom(n*p,1,0.5)-1), n,p)
    #White_G=matrix( (sqrt(3/5)*rt(n*p, df=5)) , nrow = n, ncol = p)
    #White_G=matrix( (sqrt(6/8)*rt(n*p, df=8)) , nrow = n, ncol = p)
    #White_G=matrix( (sqrt(1/2)*rt(n*p, df=4)) , nrow = n, ncol = p)
    
    eps_G=White_G %*%Sigma0.5_G
    
    
    
    G=one.n%*%t(beta0)+TRT%*%t(B) +eps_G # G is n by p matrix whose ith row is the ith obs of the mediation
    
    colnames(G) = paste0("G", 1:ncol(G))
    
    eps_Y=as.matrix(rnorm(n, mean = 0, sd=sigma_Y)) # normally distributed
    #eps_Y= sigma_Y * as.matrix(runif(n, min = -sqrt(3), max = sqrt(3))) # uniformly distributed
    #eps_Y = sigma_Y * as.matrix(rpois(n, lambda = 1)-1) # poisson centered
    #eps_Y=sigma_Y * as.matrix(2*rbinom(n,1,0.5)-1) # symmetric bernoulli
    #eps_Y=sigma_Y * as.matrix( sqrt(3/5)*rt(n, df=5) )
    #eps_Y=sigma_Y * as.matrix( sqrt(6/8)*rt(n, df=8) )
    
    
    G.T=rbind(G[1:(n/2),], matrix(0, (n/2), p))
    
    colnames(G.T) = paste0("TRT.G", 1:ncol(G))
    
    Y=theta0*one.n + theta.T*TRT + G%*%theta.G + G.T%*%Theta.D + eps_Y  # Y is n by 1 matrix
    
    
    
    ################# Our Method ####################################
    
    # step 1 outcome model
    
    W_outcome=cbind(TRT, G, G.T)
    
    fitp.scad=cv.ncvreg(W_outcome,Y, penalty = "SCAD", penalty.factor = p.fac, lambda=exp(seq(10,-10, length.out=100)) )
    
    theta0.hat = coef(fitp.scad, s = "lambda.min")[1]
    
    theta.T.hat = coef(fitp.scad, s = "lambda.min")[1+q]
    
    theta.G.hat = coef(fitp.scad, s = "lambda.min")[(1+q+1):(1+q+p)]
    
    Theta.D.hat = coef(fitp.scad, s = "lambda.min")[(1+q+p+1):(1+q+p+p)]
    
    print(sum(abs(coef(fitp.scad, s = "lambda.min") - theta)))
    print(sum(sign(coef(fitp.scad, s = "lambda.min"))!=sign(theta)))
    
    names(theta0.hat)=NULL; names(theta.T.hat)=NULL; names(theta.G.hat)=NULL; names(Theta.D.hat)=NULL
    
    bias_theta.scad[loop.run,loop.c] = sum(abs(coef(fitp.scad, s = "lambda.min") - theta))
    sr_theta.scad[loop.run,loop.c] = sum(sign(coef(fitp.scad, s = "lambda.min"))!=sign(theta))
    
    # step 2 mediation model
    
    W_mediator=cbind(one.n, TRT)
    
    #beta.theta.G.hat=solve(t(W_mediator)%*%W_mediator) %*% t(W_mediator) %*% G %*% theta.G.hat
    
    #beta.theta.D.hat=solve(t(W_mediator)%*%W_mediator) %*% t(W_mediator) %*% G %*% Theta.D.hat
    
    # calculate the estimation error
    
    SE_smpl_scad = OraProjInt_SE_sample_noX(Y = Y, Trt= TRT, Mediator = G, Trt.Medator = G.T,
                                            theta0_hat = theta0.hat, 
                                            theta.T_hat = theta.T.hat, 
                                            theta.G_hat = theta.G.hat, 
                                            vtheta.D_hat = Theta.D.hat)
    # NDE
    est_drct.10.scad[loop.run,loop.c] = SE_smpl_scad$NDE_inference["NDE(T 1 = 0)","Estimate"] 
    est_drct.11.scad[loop.run,loop.c] = SE_smpl_scad$NDE_inference["NDE(T 1 = 1)","Estimate"]
    
    SE_drct.10.scad[loop.run, loop.c] = SE_smpl_scad$NDE_inference["NDE(T 1 = 0)", "SE"]
    SE_drct.11.scad[loop.run, loop.c] = SE_smpl_scad$NDE_inference["NDE(T 1 = 1)", "SE"]
    
    if ( ( ( NDE0 ) <= ( est_drct.10.scad[loop.run,loop.c] + qnorm(0.975)*SE_drct.10.scad[loop.run, loop.c] ) ) & ( ( NDE0 ) >= ( est_drct.10.scad[loop.run,loop.c] - qnorm(0.975)*SE_drct.10.scad[loop.run, loop.c] ) )  ){
      cp_drct.10.scad[loop.run, loop.c] = 1
    }
    
    if ( ( ( NDE1 ) <= ( est_drct.11.scad[loop.run,loop.c] + qnorm(0.975)*SE_drct.11.scad[loop.run, loop.c] ) ) & ( ( NDE1 ) >= ( est_drct.11.scad[loop.run,loop.c] - qnorm(0.975)*SE_drct.11.scad[loop.run, loop.c] ) )  ){
      cp_drct.11.scad[loop.run, loop.c] = 1
    }
    
    if ( SE_smpl_scad$NDE_inference["NDE(T 1 = 0)", "p value"]<0.05 ){
      reject_drct.10.scad[loop.run, loop.c] = 1
    }
    
    if ( SE_smpl_scad$NDE_inference["NDE(T 1 = 1)", "p value"]<0.05 ){
      reject_drct.11.scad[loop.run, loop.c] = 1
    }
    
    # NIE
    est_indrct.10.scad[loop.run,loop.c] = SE_smpl_scad$NIE_inference["NIE(T 1 = 0)","Estimate"] 
    est_indrct.11.scad[loop.run,loop.c] = SE_smpl_scad$NIE_inference["NIE(T 1 = 1)","Estimate"]
    
    SE_indrct.10.scad[loop.run, loop.c] = SE_smpl_scad$NIE_inference["NIE(T 1 = 0)", "SE"]
    SE_indrct.11.scad[loop.run, loop.c] = SE_smpl_scad$NIE_inference["NIE(T 1 = 1)", "SE"]
    
    if ( ( ( t(B)%*%theta.G ) <= ( est_indrct.10.scad[loop.run,loop.c] + qnorm(0.975)*SE_indrct.10.scad[loop.run, loop.c] ) ) & ( ( t(B)%*%theta.G ) >= ( est_indrct.10.scad[loop.run,loop.c] - qnorm(0.975)*SE_indrct.10.scad[loop.run, loop.c] ) )  ){
      cp_indrct.10.scad[loop.run, loop.c] = 1
    }
    
    if ( ( ( t(B)%*%theta.G+t(Theta.D)%*%B ) <= ( est_indrct.11.scad[loop.run,loop.c] + qnorm(0.975)*SE_indrct.11.scad[loop.run, loop.c] ) ) & ( ( t(B)%*%theta.G+t(Theta.D)%*%B ) >= ( est_indrct.11.scad[loop.run,loop.c] - qnorm(0.975)*SE_indrct.11.scad[loop.run, loop.c] ) )  ){
      cp_indrct.11.scad[loop.run, loop.c] = 1
    }
    
    if ( SE_smpl_scad$NIE_inference["NIE(T 1 = 0)", "p value"]<0.05 ){
      reject_indrct.10.scad[loop.run, loop.c] = 1
    }
    
    if ( SE_smpl_scad$NIE_inference["NIE(T 1 = 1)", "p value"]<0.05 ){
      reject_indrct.11.scad[loop.run, loop.c] = 1
    }
    #################################################################################################
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    ##################### oracle #########################################################
    A_G = which(theta.G != 0) # A_G is a vector, containing number location in 1:p with length s.G
    A_D = list() # A_D has q components. Each components contains number location in 1:p with length s.Dj
    for (j in 1:q) {
      A_D[[j]] = which(Theta.D[,j]!=0)
    }
    
    W_oracle = cbind(one.n, TRT, G[, A_G])
    for (j in 1:q) {
      tmpGT = G.T[, ((j-1)*p+1) : (j*p) ]
      W_oracle = cbind(W_oracle, G.T[ , A_D[[j]] ])
    }
    rm(tmpGT)
    
    theta_w_oracle = solve(t(W_oracle)%*%W_oracle)%*%t(W_oracle)%*%Y
    
    theta0.hat_oracle = theta_w_oracle[1]
    
    theta.T.hat_oracle = theta_w_oracle[1+q]
    
    if(length(A_G) == 0){
      theta.G.hat_oracle = rep( 0, p )  # I know where A_G are so I can do this (only in this special A_G)
    } else {
      theta.G.hat_oracleAG = theta_w_oracle[(1+q+1):(1+q+length(A_G))]
      theta.G.hat_oracle = c(theta.G.hat_oracleAG, rep(0, p-length(A_G) ) ) # I know where A_G are so I can do this (only in this special A_G)
    }
    
    if(length(A_D[[1]]) == 0){
      Theta.D.hat_oracle = rep( 0, p )  # I know where A_G are so I can do this (only in this special A_G)
    } else {
      Theta.D.hat_oracleAD = theta_w_oracle[(1+q+length(A_G)+1):(1+q+length(A_G)+length(A_D[[1]]))]
      Theta.D.hat_oracle = c(Theta.D.hat_oracleAD, rep(0, p-length(A_D[[1]]) ) ) # I know where A_D[[1]] are so I can do this. (only in this special A_D[[1]])
    }
    
    
    
    names(theta0.hat_oracle)=NULL; names(theta.T.hat_oracle)=NULL; names(theta.G.hat_oracle)=NULL; names(Theta.D.hat_oracle)=NULL
    
    beta.theta.G.hat_oracle=solve(t(W_mediator)%*%W_mediator) %*% t(W_mediator) %*% G %*% theta.G.hat_oracle
    
    beta.theta.D.hat_oracle=solve(t(W_mediator)%*%W_mediator) %*% t(W_mediator) %*% G %*% Theta.D.hat_oracle
    
    # calculate the estimation error
    
    SE_smpl_oracle = OraProjInt_SE_sample_noX(Y = Y, Trt = TRT, Mediator = G, Trt.Medator = G.T,
                                              theta0_hat = theta0.hat_oracle, 
                                              theta.T_hat = theta.T.hat_oracle, 
                                              theta.G_hat = theta.G.hat_oracle, 
                                              vtheta.D_hat = Theta.D.hat_oracle)
    # NDE
    est_drct.10.oracle[loop.run, loop.c] = SE_smpl_oracle$NDE_inference["NDE(T 1 = 0)", "Estimate"]
    est_drct.11.oracle[loop.run, loop.c] = SE_smpl_oracle$NDE_inference["NDE(T 1 = 1)", "Estimate"]
    
    SE_drct.10.oracle[loop.run, loop.c] = SE_smpl_oracle$NDE_inference["NDE(T 1 = 0)", "SE"]
    SE_drct.11.oracle[loop.run, loop.c] = SE_smpl_oracle$NDE_inference["NDE(T 1 = 1)", "SE"]
    
    if ( ( ( NDE0 ) <= ( est_drct.10.oracle[loop.run,loop.c] + qnorm(0.975)*SE_drct.10.oracle[loop.run, loop.c] ) ) & ( ( NDE0 ) >= ( est_drct.10.oracle[loop.run,loop.c] - qnorm(0.975)*SE_drct.10.oracle[loop.run, loop.c] ) )  ){
      cp_drct.10.oracle[loop.run, loop.c] = 1
    }
    
    if ( ( ( NDE1 ) <= ( est_drct.11.oracle[loop.run,loop.c] + qnorm(0.975)*SE_drct.11.oracle[loop.run, loop.c] ) ) & ( ( NDE1 ) >= ( est_drct.11.oracle[loop.run,loop.c] - qnorm(0.975)*SE_drct.11.oracle[loop.run, loop.c] ) )  ){
      cp_drct.11.oracle[loop.run, loop.c] = 1
    }
    
    if ( SE_smpl_oracle$NDE_inference["NDE(T 1 = 0)", "p value"]<0.05 ){
      reject_drct.10.oracle[loop.run, loop.c] = 1
    }
    
    if ( SE_smpl_oracle$NDE_inference["NDE(T 1 = 1)", "p value"]<0.05 ){
      reject_drct.11.oracle[loop.run, loop.c] = 1
    }
    
    # NIE
    est_indrct.10.oracle[loop.run, loop.c] = beta.theta.G.hat_oracle[2,1]
    est_indrct.11.oracle[loop.run, loop.c] = beta.theta.G.hat_oracle[2,1] + beta.theta.D.hat_oracle[2,1]
    SE_indrct.10.oracle[loop.run, loop.c] = SE_smpl_oracle$NIE_inference["NIE(T 1 = 0)", "SE"]
    SE_indrct.11.oracle[loop.run, loop.c] = SE_smpl_oracle$NIE_inference["NIE(T 1 = 1)", "SE"]
    
    if ( ( ( t(B)%*%theta.G ) <= ( est_indrct.10.oracle[loop.run,loop.c] + qnorm(0.975)*SE_indrct.10.oracle[loop.run, loop.c] ) ) & ( ( t(B)%*%theta.G ) >= ( est_indrct.10.oracle[loop.run,loop.c] - qnorm(0.975)*SE_indrct.10.oracle[loop.run, loop.c] ) )  ){
      cp_indrct.10.oracle[loop.run, loop.c] = 1
    }
    
    if ( ( ( t(B)%*%theta.G+t(Theta.D)%*%B ) <= ( est_indrct.11.oracle[loop.run,loop.c] + qnorm(0.975)*SE_indrct.11.oracle[loop.run, loop.c] ) ) & ( ( t(B)%*%theta.G+t(Theta.D)%*%B ) >= ( est_indrct.11.oracle[loop.run,loop.c] - qnorm(0.975)*SE_indrct.11.oracle[loop.run, loop.c] ) )  ){
      cp_indrct.11.oracle[loop.run, loop.c] = 1
    }
    
    if ( SE_smpl_oracle$NIE_inference["NIE(T 1 = 0)", "p value"]<0.05 ){
      reject_indrct.10.oracle[loop.run, loop.c] = 1
    }
    
    if ( SE_smpl_oracle$NIE_inference["NIE(T 1 = 1)", "p value"]<0.05 ){
      reject_indrct.11.oracle[loop.run, loop.c] = 1
    }
    #############################################################################################
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    ################## Runze ##############################################################
    W_outcome_runze=cbind(TRT, G)
    
    fitp.scad_runze = cv.ncvreg(W_outcome_runze, Y, penalty = "SCAD", penalty.factor = p.fac_runze, lambda=exp(seq(10,-10, length.out=100)) )
    
    theta0.hat_runze = coef(fitp.scad_runze, s = "lambda.min")[1]
    
    theta.T.hat_runze = coef(fitp.scad_runze, s = "lambda.min")[1+q]
    
    theta.G.hat_runze = coef(fitp.scad_runze, s = "lambda.min")[(1+q+1):(1+q+p)]
    
    names(theta0.hat_runze)=NULL; names(theta.T.hat_runze)=NULL; names(theta.G.hat_runze)=NULL
    
    # step 2 mediation model
    
    W_mediator=cbind(one.n, TRT)
    SE_smpl_runze = Runze_SE_sample_noX(Y = Y, Trt= TRT, Mediator = G, 
                                        theta0_hat = theta0.hat_runze, 
                                        theta.T_hat = theta.T.hat_runze, 
                                        theta.G_hat = theta.G.hat_runze)
    
    # NDE
    est_drct.runze[loop.run, loop.c] = SE_smpl_runze$NDE_est[2,1]
    SE_drct.runze[loop.run, loop.c] = SE_smpl_runze$NDE_SE[2,1]
    
    
    if ( ( ( NDE0 ) <= ( est_drct.runze[loop.run,loop.c] + qnorm(0.975)*SE_drct.runze[loop.run, loop.c] ) ) & ( ( NDE0 ) >= ( est_drct.runze[loop.run,loop.c] - qnorm(0.975)*SE_drct.runze[loop.run, loop.c] ) )  ){
      cp_drct10.runze[loop.run, loop.c] = 1
    }
    
    if ( ( ( NDE1 ) <= ( est_drct.runze[loop.run,loop.c] + qnorm(0.975)*SE_drct.runze[loop.run, loop.c] ) ) & ( ( NDE1 ) >= ( est_drct.runze[loop.run,loop.c] - qnorm(0.975)*SE_drct.runze[loop.run, loop.c] ) )  ){
      cp_drct11.runze[loop.run, loop.c] = 1
    }
    
    if ( ( 0 > ( est_drct.runze[loop.run,loop.c] + qnorm(0.975)*SE_drct.runze[loop.run, loop.c] ) ) | ( 0 < ( est_drct.runze[loop.run,loop.c] - qnorm(0.975)*SE_drct.runze[loop.run, loop.c] ) )  ){
      reject_drct.runze[loop.run, loop.c] = 1
    }
    
    # NIE
    if( length( which(theta.G.hat_runze!=0) ) == 0 ){
      est_indrct.runze[loop.run, loop.c] = 0
      SE_indrct.runze[loop.run, loop.c] = 0
    } else {
      est_indrct.runze[loop.run, loop.c] = SE_smpl_runze$NIE_est[2,1]
      SE_indrct.runze[loop.run, loop.c] = SE_smpl_runze$NIE_SE[2,1]
    }
    
    if ( ( ( t(B)%*%theta.G ) <= ( est_indrct.runze[loop.run,loop.c] + qnorm(0.975)*SE_indrct.runze[loop.run, loop.c] ) ) & ( ( t(B)%*%theta.G ) >= ( est_indrct.runze[loop.run,loop.c] - qnorm(0.975)*SE_indrct.runze[loop.run, loop.c] ) )  ){
      cp_indrct10.runze[loop.run, loop.c] = 1
    }
    
    if ( ( ( t(B)%*%theta.G+t(Theta.D)%*%B ) <= ( est_indrct.runze[loop.run,loop.c] + qnorm(0.975)*SE_indrct.runze[loop.run, loop.c] ) ) & ( ( t(B)%*%theta.G+t(Theta.D)%*%B ) >= ( est_indrct.runze[loop.run,loop.c] - qnorm(0.975)*SE_indrct.runze[loop.run, loop.c] ) )  ){
      cp_indrct11.runze[loop.run, loop.c] = 1
    }
    
    if ( ( 0 > ( est_indrct.runze[loop.run,loop.c] + qnorm(0.975)*SE_indrct.runze[loop.run, loop.c] ) ) | ( 0 < ( est_indrct.runze[loop.run,loop.c] - qnorm(0.975)*SE_indrct.runze[loop.run, loop.c] ) )  ){
      reject_indrct.runze[loop.run, loop.c] = 1
    }
    ############################################################################################
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    ######################### Lasso ###########################################################
    # step 1 outcome model
    
    W_outcome=cbind(TRT, G, G.T)
    
    fitp.lasso=cv.glmnet(W_outcome, Y,  penalty.factor=p.fac_lasso, lambda = exp(seq(-10,12,length.out=100 )))
    
    theta0.hat_lasso = coef(fitp.lasso, s = "lambda.min")[1]
    
    theta.T.hat_lasso = coef(fitp.lasso, s = "lambda.min")[1+q]
    
    theta.G.hat_lasso = coef(fitp.lasso, s = "lambda.min")[(1+q+1):(1+q+p)]
    
    Theta.D.hat_lasso = coef(fitp.lasso, s = "lambda.min")[(1+q+p+1):(1+q+p+p)]
    
    names(theta0.hat_lasso)=NULL; names(theta.T.hat_lasso)=NULL; names(theta.G.hat_lasso)=NULL; names(Theta.D.hat_lasso)=NULL
    
    bias_theta.lasso[loop.run,loop.c] = sum(abs(coef(fitp.lasso, s = "lambda.min") - theta))
    sr_theta.lasso[loop.run,loop.c] = sum(sign(coef(fitp.lasso, s = "lambda.min"))!=sign(theta))
    
    # step 2 mediation model
    
    W_mediator=cbind(one.n, TRT)
    
    beta.theta.G.hat_lasso=solve(t(W_mediator)%*%W_mediator) %*% t(W_mediator) %*% G %*% theta.G.hat_lasso
    beta.theta.D.hat_lasso=solve(t(W_mediator)%*%W_mediator) %*% t(W_mediator) %*% G %*% Theta.D.hat_lasso
    
    # calculate the estimation error
    # NDE
    est_drct.10.lasso[loop.run,loop.c] = theta.T.hat_lasso + beta.theta.D.hat_lasso[1,1]
    est_drct.11.lasso[loop.run,loop.c] = theta.T.hat_lasso + beta.theta.D.hat_lasso[1,1] + beta.theta.D.hat_lasso[2,1]
    
    bias_drct.10.lasso[loop.run, loop.c] = est_drct.10.lasso[loop.run,loop.c] - NDE0
    bias_drct.11.lasso[loop.run, loop.c] = est_drct.11.lasso[loop.run,loop.c] - NDE1
    
    # NIE
    est_indrct.10.lasso[loop.run,loop.c] = beta.theta.G.hat_lasso[2,1]
    est_indrct.11.lasso[loop.run,loop.c] = beta.theta.G.hat_lasso[2,1] + beta.theta.D.hat_lasso[2,1]
    
    bias_indrct.10.lasso[loop.run, loop.c] = est_indrct.10.lasso[loop.run,loop.c] - t(B)%*%theta.G
    bias_indrct.11.lasso[loop.run, loop.c] = est_indrct.11.lasso[loop.run,loop.c] - (t(B)%*%theta.G + t(Theta.D)%*%B)
    #######################################################################################################################
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    ################################## Dave debiased lasso ######################################################
    Inference_Dave.cY = Dave(Y=as.vector(Y), G=G, S=as.matrix(TRT), center.SG = T, center.Y = T, lam_list = ( sqrt(log(p)/n)/3 ) ) 
    
    # NDE
    est_drct.dave.cY[loop.run, loop.c] = as.numeric(Inference_Dave.cY$inference[2,"est"])
    SE_drct.dave.cY[loop.run, loop.c] = as.numeric(Inference_Dave.cY$inference[2,"SE"])
    if ( ( ( NDE0 ) <= ( est_drct.dave.cY[loop.run,loop.c] + qnorm(0.975)*SE_drct.dave.cY[loop.run, loop.c] ) ) & ( ( NDE0 ) >= ( est_drct.dave.cY[loop.run,loop.c] - qnorm(0.975)*SE_drct.dave.cY[loop.run, loop.c] ) )  ){
      cp_drct10.dave.cY[loop.run, loop.c] = 1
    }
    
    if ( ( ( NDE1 ) <= ( est_drct.dave.cY[loop.run,loop.c] + qnorm(0.975)*SE_drct.dave.cY[loop.run, loop.c] ) ) & ( ( NDE1 ) >= ( est_drct.dave.cY[loop.run,loop.c] - qnorm(0.975)*SE_drct.dave.cY[loop.run, loop.c] ) )  ){
      cp_drct11.dave.cY[loop.run, loop.c] = 1
    }
    
    if ( ( 0 > ( est_drct.dave.cY[loop.run,loop.c] + qnorm(0.975)*SE_drct.dave.cY[loop.run, loop.c] ) ) | ( 0 < ( est_drct.dave.cY[loop.run,loop.c] - qnorm(0.975)*SE_drct.dave.cY[loop.run, loop.c] ) )  ){
      reject_drct.dave.cY[loop.run, loop.c] = 1
    }
    
    # NIE
    est_indrct.dave.cY[loop.run, loop.c] = as.numeric(Inference_Dave.cY$inference[1,"est"])
    SE_indrct.dave.cY[loop.run, loop.c] = as.numeric(Inference_Dave.cY$inference[1,"SE"])
    if ( ( ( t(B)%*%theta.G ) <= ( est_indrct.dave.cY[loop.run,loop.c] + qnorm(0.975)*SE_indrct.dave.cY[loop.run, loop.c] ) ) & ( ( t(B)%*%theta.G ) >= ( est_indrct.dave.cY[loop.run,loop.c] - qnorm(0.975)*SE_indrct.dave.cY[loop.run, loop.c] ) )  ){
      cp_indrct10.dave.cY[loop.run, loop.c] = 1
    }
    
    if ( ( ( t(B)%*%theta.G+t(Theta.D)%*%B ) <= ( est_indrct.dave.cY[loop.run,loop.c] + qnorm(0.975)*SE_indrct.dave.cY[loop.run, loop.c] ) ) & ( ( t(B)%*%theta.G+t(Theta.D)%*%B ) >= ( est_indrct.dave.cY[loop.run,loop.c] - qnorm(0.975)*SE_indrct.dave.cY[loop.run, loop.c] ) )  ){
      cp_indrct11.dave.cY[loop.run, loop.c] = 1
    }
    
    if ( ( 0 > ( est_indrct.dave.cY[loop.run,loop.c] + qnorm(0.975)*SE_indrct.dave.cY[loop.run, loop.c] ) ) | ( 0 < ( est_indrct.dave.cY[loop.run,loop.c] - qnorm(0.975)*SE_indrct.dave.cY[loop.run, loop.c] ) )  ){
      reject_indrct.dave.cY[loop.run, loop.c] = 1
    }
    ####################################################################################################
    
    
    
    
    
    
    
    
    
    ######################################## Huang 2016 ###########################
    Inference_Huang = Huang2016_noX1(Y=Y, Trt = TRT, Mediator = G, Trt.Mediator = G.T, var_per = 0.8, boot = 500)
    
    est_indrct.10.Huang[loop.run,loop.c] = Inference_Huang$NIE_inference["NIE(T=0)", "Est"]
    est_indrct.11.Huang[loop.run,loop.c] = Inference_Huang$NIE_inference["NIE(T=1)", "Est"]
    
    SE_indrct.10.Huang[loop.run, loop.c] = Inference_Huang$NIE_inference["NIE(T=0)", "SE"]
    SE_indrct.11.Huang[loop.run, loop.c] = Inference_Huang$NIE_inference["NIE(T=1)", "SE"]
    
    if ( ( ( t(B)%*%theta.G ) <= ( est_indrct.10.Huang[loop.run,loop.c] + qnorm(0.975)*SE_indrct.10.Huang[loop.run, loop.c] ) ) & ( ( t(B)%*%theta.G ) >= ( est_indrct.10.Huang[loop.run,loop.c] - qnorm(0.975)*SE_indrct.10.Huang[loop.run, loop.c] ) )  ){
      cp_indrct.10.Huang[loop.run, loop.c] = 1
    }
    
    if ( ( ( t(B)%*%theta.G+t(Theta.D)%*%B ) <= ( est_indrct.11.Huang[loop.run,loop.c] + qnorm(0.975)*SE_indrct.11.Huang[loop.run, loop.c] ) ) & ( ( t(B)%*%theta.G+t(Theta.D)%*%B ) >= ( est_indrct.11.Huang[loop.run,loop.c] - qnorm(0.975)*SE_indrct.11.Huang[loop.run, loop.c] ) )  ){
      cp_indrct.11.Huang[loop.run, loop.c] = 1
    }
    
    if ( Inference_Huang$NIE_inference["NIE(T=0)", "p-value"]<0.05 ){
      reject_indrct.10.Huang[loop.run, loop.c] = 1
    }
    
    if ( Inference_Huang$NIE_inference["NIE(T=1)", "p-value"]<0.05 ){
      reject_indrct.11.Huang[loop.run, loop.c] = 1
    }
    
    
    
    # bootstrap SE
    SE_indrct.10.Huang.boot[loop.run, loop.c] = Inference_Huang$bootSE.NIE0
    SE_indrct.11.Huang.boot[loop.run, loop.c] = Inference_Huang$bootSE.NIE1
    
    if ( ( ( t(B)%*%theta.G ) <= ( est_indrct.10.Huang[loop.run,loop.c] + qnorm(0.975)*SE_indrct.10.Huang.boot[loop.run, loop.c] ) ) & ( ( t(B)%*%theta.G ) >= ( est_indrct.10.Huang[loop.run,loop.c] - qnorm(0.975)*SE_indrct.10.Huang.boot[loop.run, loop.c] ) )  ){
      cp_indrct.10.Huang.boot[loop.run, loop.c] = 1
    }
    
    if ( ( ( t(B)%*%theta.G+t(Theta.D)%*%B ) <= ( est_indrct.11.Huang[loop.run,loop.c] + qnorm(0.975)*SE_indrct.11.Huang.boot[loop.run, loop.c] ) ) & ( ( t(B)%*%theta.G+t(Theta.D)%*%B ) >= ( est_indrct.11.Huang[loop.run,loop.c] - qnorm(0.975)*SE_indrct.11.Huang.boot[loop.run, loop.c] ) )  ){
      cp_indrct.11.Huang.boot[loop.run, loop.c] = 1
    }
    
    if (  abs( est_indrct.10.Huang[loop.run,loop.c]/SE_indrct.10.Huang.boot[loop.run, loop.c] )>qnorm(0.975) ){
      reject_indrct.10.Huang.boot[loop.run, loop.c] = 1
    }
    
    if ( abs( est_indrct.11.Huang[loop.run,loop.c]/SE_indrct.11.Huang.boot[loop.run, loop.c] )>qnorm(0.975) ){
      reject_indrct.11.Huang.boot[loop.run, loop.c] = 1
    }
    
    PCAbias_indrct.10.Huang[loop.run, loop.c] = t(B)%*%theta.G - t(B)%*%(Inference_Huang$PCA_sigmaG$U)%*%t( Inference_Huang$PCA_sigmaG$U )%*%theta.G
    PCAbias_indrct.11.Huang[loop.run, loop.c] = t(B)%*%(theta.G+Theta.D) - t(B)%*%(Inference_Huang$PCA_sigmaG$U)%*%t( Inference_Huang$PCA_sigmaG$U )%*%(theta.G+Theta.D)
    ##########################################################################################################
    
    
    
    
    
    
    
    
    
    print(loop.run)
    
  }
  print(loop.c)
}





result = list( 
  run=run,
  n=n,
  p=p,
  q=q,
  one.n=one.n,
  TRT=TRT,
  theta0=theta0,
  c.T=c.T,
  theta.T=theta.T,
  c.G = c.G,
  s.G = s.G,
  theta.G=theta.G,
  cstar.seq = cstar.seq,
  s.D=s.D,
  sigma_Y=sigma_Y,
  c_beta0=c_beta0,
  beta0=beta0,
  c_B = c_B,
  B=B,
  rho_G=rho_G,
  c_sigmaG=c_sigmaG,
  noc=noc,
  blk.size=blk.size,
  Sigma_G=Sigma_G,
  p.fac=p.fac,
  p.fac_runze = p.fac_runze,
  p.fac_lasso = p.fac_lasso ,
  bias_theta.scad=bias_theta.scad,
  sr_theta.scad=sr_theta.scad,
  est_drct.10.scad=est_drct.10.scad,
  est_drct.11.scad=est_drct.11.scad,
  SE_drct.10.scad=SE_drct.10.scad,
  SE_drct.11.scad=SE_drct.11.scad,
  cp_drct.10.scad=cp_drct.10.scad,
  cp_drct.11.scad=cp_drct.11.scad,
  reject_drct.10.scad=reject_drct.10.scad,
  reject_drct.11.scad=reject_drct.11.scad,
  est_indrct.10.scad=est_indrct.10.scad,
  est_indrct.11.scad=est_indrct.11.scad,
  SE_indrct.10.scad=SE_indrct.10.scad,
  SE_indrct.11.scad=SE_indrct.11.scad,
  cp_indrct.10.scad=cp_indrct.10.scad,
  cp_indrct.11.scad=cp_indrct.11.scad,
  reject_indrct.10.scad=reject_indrct.10.scad,
  reject_indrct.11.scad=reject_indrct.11.scad,
  est_drct.10.oracle=est_drct.10.oracle,
  est_drct.11.oracle=est_drct.11.oracle,
  SE_drct.10.oracle=SE_drct.10.oracle,
  SE_drct.11.oracle=SE_drct.11.oracle,
  cp_drct.10.oracle=cp_drct.10.oracle,
  cp_drct.11.oracle=cp_drct.11.oracle,
  reject_drct.10.oracle=reject_drct.10.oracle,
  reject_drct.11.oracle=reject_drct.11.oracle,
  est_indrct.10.oracle=est_indrct.10.oracle,
  est_indrct.11.oracle=est_indrct.11.oracle,
  SE_indrct.10.oracle=SE_indrct.10.oracle,
  SE_indrct.11.oracle=SE_indrct.11.oracle,
  cp_indrct.10.oracle=cp_indrct.10.oracle,
  cp_indrct.11.oracle=cp_indrct.11.oracle,
  reject_indrct.10.oracle=reject_indrct.10.oracle,
  reject_indrct.11.oracle=reject_indrct.11.oracle,
  est_drct.runze=est_drct.runze,
  SE_drct.runze=SE_drct.runze,
  cp_drct10.runze=cp_drct10.runze,
  cp_drct11.runze=cp_drct11.runze,
  reject_drct.runze=reject_drct.runze,
  est_indrct.runze=est_indrct.runze,
  SE_indrct.runze=SE_indrct.runze,
  cp_indrct10.runze=cp_indrct10.runze,
  cp_indrct11.runze=cp_indrct11.runze,
  reject_indrct.runze=reject_indrct.runze,
  bias_theta.lasso=bias_theta.lasso,
  sr_theta.lasso=sr_theta.lasso,
  est_drct.10.lasso=est_drct.10.lasso,
  est_drct.11.lasso=est_drct.11.lasso,
  bias_drct.10.lasso=bias_drct.10.lasso,
  bias_drct.11.lasso=bias_drct.11.lasso,
  est_indrct.10.lasso=est_indrct.10.lasso,
  est_indrct.11.lasso=est_indrct.11.lasso,
  bias_indrct.10.lasso=bias_indrct.10.lasso,
  bias_indrct.11.lasso=bias_indrct.11.lasso,
  est_drct.dave.cY=est_drct.dave.cY,
  SE_drct.dave.cY=SE_drct.dave.cY,
  cp_drct10.dave.cY=cp_drct10.dave.cY,
  cp_drct11.dave.cY=cp_drct11.dave.cY,
  reject_drct.dave.cY=reject_drct.dave.cY,
  est_indrct.dave.cY=est_indrct.dave.cY,
  SE_indrct.dave.cY=SE_indrct.dave.cY,
  cp_indrct10.dave.cY=cp_indrct10.dave.cY,
  cp_indrct11.dave.cY=cp_indrct11.dave.cY,
  reject_indrct.dave.cY=reject_indrct.dave.cY,
  est_indrct.10.Huang= est_indrct.10.Huang,
  est_indrct.11.Huang= est_indrct.11.Huang,
  SE_indrct.10.Huang= SE_indrct.10.Huang,
  SE_indrct.11.Huang= SE_indrct.11.Huang,
  cp_indrct.10.Huang= cp_indrct.10.Huang,
  cp_indrct.11.Huang= cp_indrct.11.Huang,
  reject_indrct.10.Huang= reject_indrct.10.Huang,
  reject_indrct.11.Huang= reject_indrct.11.Huang,
  SE_indrct.10.Huang.boot= SE_indrct.10.Huang.boot,
  SE_indrct.11.Huang.boot= SE_indrct.11.Huang.boot,
  cp_indrct.10.Huang.boot= cp_indrct.10.Huang.boot,
  cp_indrct.11.Huang.boot= cp_indrct.11.Huang.boot,
  reject_indrct.10.Huang.boot= reject_indrct.10.Huang.boot,
  reject_indrct.11.Huang.boot= reject_indrct.11.Huang.boot,
  PCAbias_indrct.10.Huang = PCAbias_indrct.10.Huang,
  PCAbias_indrct.11.Huang = PCAbias_indrct.11.Huang
)


# tracemem(result)

library(gdata)

mv(from = "result", to = paste0("result", simnum))

save(list = c( paste0("result", simnum) ),
     file = paste0("R_out/cD_p500_1_", simnum, ".Rdata"))


