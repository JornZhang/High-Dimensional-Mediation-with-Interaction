



########################################## Glasso for Sigma_G #####################################

############# This part adjust the parameters from the original real data with p = 1008
rm(list = ls())
library(ncvreg)
library(glasso)
load("parameters from real data.Rdata")




## since we set p=1000, we only use the first 1000 eps_Gj to generate the covariance matrix of eps_G
#################### use glasso to choose sigma_G
eps_G1 = eps_G[, c(1:1000) ]
sigma_G1 = t(eps_G1)%*%eps_G1 / nrow(eps_G1)

tuning_glasso = exp(seq(-2.5,0,length.out=200)) # the penalty rho least in glasso

glasso_path = glassopath(sigma_G1, rholist = tuning_glasso  ) 

prop0_G = c() # the proportion of nonzero entries in the covariance matrix
sum_abs_G=c() # the sum of absolute values of the entries in the covariance matrix
# minmax_eigen_G = matrix(0, nrow = length(tuning_glasso), ncol = 2)
# colnames(minmax_eigen_G) = c("min_eigen_value", "max_eigen_value")
for (i in 1:length(tuning_glasso)) {
  prop0_G[i] = sum(glasso_path$w[,,i]!=0) / ( nrow(glasso_path$w[,,i])^2 )
  # tmp_eigenG = eigen(glasso_path$w[,,i])$values
  # minmax_eigen_G[i,"min_eigen_value"] = tmp_eigenG[length(tmp_eigenG)]
  # minmax_eigen_G[i,"max_eigen_value"] = tmp_eigenG[1]
  sum_abs_G[i] = sum(abs(glasso_path$w[,,i]))
  print(i)
}
# rm(tmp_eigenG)
prop0_G

boxplot(diag(glasso_path$w[,,1]), 
        diag(glasso_path$w[,,2]), 
        diag(glasso_path$w[,,3]),
        diag(glasso_path$w[,,4]),
        diag(glasso_path$w[,,5]),
        diag(glasso_path$w[,,6]),
        diag(glasso_path$w[,,7]),
        diag(glasso_path$w[,,8]),
        diag(glasso_path$w[,,9]),
        diag(glasso_path$w[,,10]),
        diag(glasso_path$w[,,11]),
        diag(glasso_path$w[,,20]))

plot(sum_abs_G)
rm(prop0_G, sum_abs_G)

################# We will use the correlation matrix of glasso_path$w[,, 1] ######################
Sigma_G_gl = cov2cor(glasso_path$w[,, 1])



boxplot(cbind(beta0, beta1, beta2, diag(Sigma_G_gl), diag( t(eps_G)%*%eps_G / nrow(eps_G) ) ), 
        names = c( expression(beta[0]), expression(beta[1]),
                   expression(beta[2]), expression(diag~(Sigma[G]^glasso )),
                   expression(diag~(Sigma[G]^real))   ),
        col= 4, main = "boxplots of the coefficients in the mediation model"  )

rm(eps_G, eps_G1, sigma_G1, tuning_glasso, glasso_path)


for (i in 1:(nrow(Sigma_G_gl)-1)) {
  for (j in (i+1):(ncol(Sigma_G_gl)-1)) {
    if(Sigma_G_gl[i,j]!=0){
            print(paste0(  "location:",i,", ", j,", value:", Sigma_G_gl[i,j]))
    }
  }
}

#######################################################################################

## Since p=1000 or 2000 rather than 1008, we need to generate other quantities by our selves.
save(list = c("Sigma_G_gl"), file = "Glasso_Sigma_G_p=1000.Rdata")










############################################## fisrt try p=1000 ################################
rm(list = ls())

p=1000
source("MyMethod_functions.R")


### glasso Sigma_G ###
load("Glasso_Sigma_G_p=1000.Rdata") # a 1000 by 1000 matrix (need adjustment when p=2000)


###### AR1 Sigma_G (do not run this if you use glasso Sigma_G) #####
rho_G=0.5 # if rho_G=0, Sigma_G=I_p * c_G
c_sigmaG=1 # c_G is the magnitude of the noise
# AR1
Gamma_G=matrix(0, p, p)
for (i in 1:p) {
  for (j in 1:p) {
    if(i == j ){
      Gamma_G[i,j] = 1
    } else {
      Gamma_G[i,j] = rho_G^(abs(i-j))
    }
  }
}
Sigma_G_AR1=c_sigmaG*Gamma_G
rm(Gamma_G, c_sigmaG)
############################################################




# Sigma_G_gl = 2 * Sigma_G_gl # change the magnitude of Sigma_G


#### find glasso Sigma)_G^{1/2}
Q_G=eigen(Sigma_G_gl)$vectors
Lambda0.5_G=diag(c(sqrt(eigen(Sigma_G_gl)$values)))
Sigma0.5_G_gl=Q_G%*%Lambda0.5_G%*%t(Q_G)
# check whether it is correct
print(max(abs(Sigma0.5_G_gl%*%Sigma0.5_G_gl-Sigma_G_gl)))
rm(Q_G, Lambda0.5_G)



#run = 500
run=100

n.seq=c(120)


# define the paramters
r=2
q=2


theta0=-2


theta.T=c(1.5, -0.5)

A_G = c(20, 27, 106,  146,  149,  268 , 497,  589,  807,  879)
A_D1 = c(106, 879)
A_D2 = c(807)


theta.G = rep(0, times=p)
# theta.G[A_G] = c(-3, 2, 1,-1, -2, 4, -2, 4, 1, -1)
theta.G[A_G] = c(-3.3, 2.3, 0.2,-1.3, -2.5, 4.3, -2.2, 3.6,-0.6, 0.3)

vtheta.D = rep(0, times=(2*p))

#vtheta.D[A_D1] = c(5, 4)
vtheta.D[A_D1] = c(5.3, 2.9)

#vtheta.D[p+A_D2] = c(5)
vtheta.D[p+A_D2] = c(3.3)

which(vtheta.D!=0)


theta.X = c(3, -3)

theta = c(theta0, theta.T, theta.X, theta.G, vtheta.D)
theta.D1 = vtheta.D[1:p]
theta.D2 = vtheta.D[(p+1):(p+p)]

sigma_Y = 1

set.seed(1234)

# beta0 = matrix( rnorm(p, sd = 0.1) , ncol=1)
beta0 = matrix( rnorm(p, sd = 2) , ncol=1)
beta0[A_G,1] = seq(0.1, 1, by=0.1)
plot(beta0)

# beta1 = matrix( rnorm(p, sd = 0.1) , ncol=1)
beta1 = matrix( rnorm(p, sd = 2) , ncol=1)

#beta1[A_G,1] = (-1)*seq(0.05, 0.5, by=0.05)
beta1[A_G,1] = (-1)*seq(0.1, 1, by=0.1)

# beta2 = matrix( rnorm(p, sd = 0.1) , ncol=1)
beta2 = matrix( rnorm(p, sd = 2) , ncol=1)

# beta2[A_G,1] = (-1)*seq(0.05, 0.5, by=0.05)
beta2[A_G,1] = (-1)*seq(0.1, 1, by=0.1)


B = cbind(beta1, beta2)

mu.X = matrix(c(0.5,0), nrow=2)
sigma.X = matrix(c(0.25,0,0,1), 2,2)
A = matrix(rnorm((p*r), sd=0.1), nrow = p, ncol = r)
#A = matrix(rnorm((p*r), sd=2), nrow = p, ncol = r)


p.fac<- rep(1, (q+r+p+q*p))
p.fac[1:(q+r)] <- 0

p.fac.runze<- rep(1, (q+r+p))
p.fac.runze[1:(q+r)] <- 0




######## proposed cv
bias_theta.scad.cv=matrix(0, nrow = run, ncol = length(n.seq))
sr_theta.scad.cv=matrix(0, nrow = run, ncol = length(n.seq))

est_indrct.10.scad.cv=matrix(0, nrow = run, ncol = length(n.seq))
est_indrct.11.scad.cv=matrix(0, nrow = run, ncol = length(n.seq))
est_indrct.20.scad.cv=matrix(0, nrow = run, ncol = length(n.seq))
est_indrct.21.scad.cv=matrix(0, nrow = run, ncol = length(n.seq))

est_drct.10.scad.cv=matrix(0, nrow = run, ncol = length(n.seq))
est_drct.11.scad.cv=matrix(0, nrow = run, ncol = length(n.seq))
est_drct.20.scad.cv=matrix(0, nrow = run, ncol = length(n.seq))
est_drct.21.scad.cv=matrix(0, nrow = run, ncol = length(n.seq))

se_indrct.10.scad.cv=matrix(0, nrow = run, ncol = length(n.seq))
se_indrct.11.scad.cv=matrix(0, nrow = run, ncol = length(n.seq))
se_indrct.20.scad.cv=matrix(0, nrow = run, ncol = length(n.seq))
se_indrct.21.scad.cv=matrix(0, nrow = run, ncol = length(n.seq))

se_drct.10.scad.cv=matrix(0, nrow = run, ncol = length(n.seq))
se_drct.11.scad.cv=matrix(0, nrow = run, ncol = length(n.seq))
se_drct.20.scad.cv=matrix(0, nrow = run, ncol = length(n.seq))
se_drct.21.scad.cv=matrix(0, nrow = run, ncol = length(n.seq))

cp_indrct.10.scad.cv=matrix(0, nrow = run, ncol = length(n.seq))
cp_indrct.11.scad.cv=matrix(0, nrow = run, ncol = length(n.seq))
cp_indrct.20.scad.cv=matrix(0, nrow = run, ncol = length(n.seq))
cp_indrct.21.scad.cv=matrix(0, nrow = run, ncol = length(n.seq))

cp_drct.10.scad.cv=matrix(0, nrow = run, ncol = length(n.seq))
cp_drct.11.scad.cv=matrix(0, nrow = run, ncol = length(n.seq))
cp_drct.20.scad.cv=matrix(0, nrow = run, ncol = length(n.seq))
cp_drct.21.scad.cv=matrix(0, nrow = run, ncol = length(n.seq))

reject_indrct.10.scad.cv=matrix(0, nrow = run, ncol = length(n.seq))
reject_indrct.11.scad.cv=matrix(0, nrow = run, ncol = length(n.seq))
reject_indrct.20.scad.cv=matrix(0, nrow = run, ncol = length(n.seq))
reject_indrct.21.scad.cv=matrix(0, nrow = run, ncol = length(n.seq))

reject_drct.10.scad.cv=matrix(0, nrow = run, ncol = length(n.seq))
reject_drct.11.scad.cv=matrix(0, nrow = run, ncol = length(n.seq))
reject_drct.20.scad.cv=matrix(0, nrow = run, ncol = length(n.seq))
reject_drct.21.scad.cv=matrix(0, nrow = run, ncol = length(n.seq))

lambda_scad.cv=matrix(0, nrow = run, ncol = length(n.seq))






######## SCAD HBIC
bias_theta.scad.hbic=matrix(0, nrow = run, ncol = length(n.seq))
sr_theta.scad.hbic=matrix(0, nrow = run, ncol = length(n.seq))

est_indrct.10.scad.hbic=matrix(0, nrow = run, ncol = length(n.seq))
est_indrct.11.scad.hbic=matrix(0, nrow = run, ncol = length(n.seq))
est_indrct.20.scad.hbic=matrix(0, nrow = run, ncol = length(n.seq))
est_indrct.21.scad.hbic=matrix(0, nrow = run, ncol = length(n.seq))

est_drct.10.scad.hbic=matrix(0, nrow = run, ncol = length(n.seq))
est_drct.11.scad.hbic=matrix(0, nrow = run, ncol = length(n.seq))
est_drct.20.scad.hbic=matrix(0, nrow = run, ncol = length(n.seq))
est_drct.21.scad.hbic=matrix(0, nrow = run, ncol = length(n.seq))


se_indrct.10.scad.hbic=matrix(0, nrow = run, ncol = length(n.seq))
se_indrct.11.scad.hbic=matrix(0, nrow = run, ncol = length(n.seq))
se_indrct.20.scad.hbic=matrix(0, nrow = run, ncol = length(n.seq))
se_indrct.21.scad.hbic=matrix(0, nrow = run, ncol = length(n.seq))

se_drct.10.scad.hbic=matrix(0, nrow = run, ncol = length(n.seq))
se_drct.11.scad.hbic=matrix(0, nrow = run, ncol = length(n.seq))
se_drct.20.scad.hbic=matrix(0, nrow = run, ncol = length(n.seq))
se_drct.21.scad.hbic=matrix(0, nrow = run, ncol = length(n.seq))

cp_indrct.10.scad.hbic=matrix(0, nrow = run, ncol = length(n.seq))
cp_indrct.11.scad.hbic=matrix(0, nrow = run, ncol = length(n.seq))
cp_indrct.20.scad.hbic=matrix(0, nrow = run, ncol = length(n.seq))
cp_indrct.21.scad.hbic=matrix(0, nrow = run, ncol = length(n.seq))

cp_drct.10.scad.hbic=matrix(0, nrow = run, ncol = length(n.seq))
cp_drct.11.scad.hbic=matrix(0, nrow = run, ncol = length(n.seq))
cp_drct.20.scad.hbic=matrix(0, nrow = run, ncol = length(n.seq))
cp_drct.21.scad.hbic=matrix(0, nrow = run, ncol = length(n.seq))

reject_indrct.10.scad.hbic=matrix(0, nrow = run, ncol = length(n.seq))
reject_indrct.11.scad.hbic=matrix(0, nrow = run, ncol = length(n.seq))
reject_indrct.20.scad.hbic=matrix(0, nrow = run, ncol = length(n.seq))
reject_indrct.21.scad.hbic=matrix(0, nrow = run, ncol = length(n.seq))

reject_drct.10.scad.hbic=matrix(0, nrow = run, ncol = length(n.seq))
reject_drct.11.scad.hbic=matrix(0, nrow = run, ncol = length(n.seq))
reject_drct.20.scad.hbic=matrix(0, nrow = run, ncol = length(n.seq))
reject_drct.21.scad.hbic=matrix(0, nrow = run, ncol = length(n.seq))

lambda_scad.hbic=matrix(0, nrow = run, ncol = length(n.seq))





## Runze cv
est_indrct.1.runze.cv=matrix(0, nrow = run, ncol = length(n.seq))
est_indrct.2.runze.cv=matrix(0, nrow = run, ncol = length(n.seq))
se_indrct.1.runze.cv=matrix(0, nrow = run, ncol = length(n.seq))
se_indrct.2.runze.cv=matrix(0, nrow = run, ncol = length(n.seq))

cp_indrct.10.runze.cv=matrix(0, nrow = run, ncol = length(n.seq))
cp_indrct.11.runze.cv=matrix(0, nrow = run, ncol = length(n.seq))
cp_indrct.20.runze.cv=matrix(0, nrow = run, ncol = length(n.seq))
cp_indrct.21.runze.cv=matrix(0, nrow = run, ncol = length(n.seq))

reject_indrct.10.runze.cv=matrix(0, nrow = run, ncol = length(n.seq))
reject_indrct.11.runze.cv=matrix(0, nrow = run, ncol = length(n.seq))
reject_indrct.20.runze.cv=matrix(0, nrow = run, ncol = length(n.seq))
reject_indrct.21.runze.cv=matrix(0, nrow = run, ncol = length(n.seq))

est_drct.1.runze.cv=matrix(0, nrow = run, ncol = length(n.seq))
est_drct.2.runze.cv=matrix(0, nrow = run, ncol = length(n.seq))
se_drct.1.runze.cv=matrix(0, nrow = run, ncol = length(n.seq))
se_drct.2.runze.cv=matrix(0, nrow = run, ncol = length(n.seq))

cp_drct.10.runze.cv=matrix(0, nrow = run, ncol = length(n.seq))
cp_drct.11.runze.cv=matrix(0, nrow = run, ncol = length(n.seq))
cp_drct.20.runze.cv=matrix(0, nrow = run, ncol = length(n.seq))
cp_drct.21.runze.cv=matrix(0, nrow = run, ncol = length(n.seq))

reject_drct.10.runze.cv=matrix(0, nrow = run, ncol = length(n.seq))
reject_drct.11.runze.cv=matrix(0, nrow = run, ncol = length(n.seq))
reject_drct.20.runze.cv=matrix(0, nrow = run, ncol = length(n.seq))
reject_drct.21.runze.cv=matrix(0, nrow = run, ncol = length(n.seq))


## Runze hbic
est_indrct.1.runze.hbic=matrix(0, nrow = run, ncol = length(n.seq))
est_indrct.2.runze.hbic=matrix(0, nrow = run, ncol = length(n.seq))
se_indrct.1.runze.hbic=matrix(0, nrow = run, ncol = length(n.seq))
se_indrct.2.runze.hbic=matrix(0, nrow = run, ncol = length(n.seq))

cp_indrct.10.runze.hbic=matrix(0, nrow = run, ncol = length(n.seq))
cp_indrct.11.runze.hbic=matrix(0, nrow = run, ncol = length(n.seq))
cp_indrct.20.runze.hbic=matrix(0, nrow = run, ncol = length(n.seq))
cp_indrct.21.runze.hbic=matrix(0, nrow = run, ncol = length(n.seq))

reject_indrct.10.runze.hbic=matrix(0, nrow = run, ncol = length(n.seq))
reject_indrct.11.runze.hbic=matrix(0, nrow = run, ncol = length(n.seq))
reject_indrct.20.runze.hbic=matrix(0, nrow = run, ncol = length(n.seq))
reject_indrct.21.runze.hbic=matrix(0, nrow = run, ncol = length(n.seq))

est_drct.1.runze.hbic=matrix(0, nrow = run, ncol = length(n.seq))
est_drct.2.runze.hbic=matrix(0, nrow = run, ncol = length(n.seq))
se_drct.1.runze.hbic=matrix(0, nrow = run, ncol = length(n.seq))
se_drct.2.runze.hbic=matrix(0, nrow = run, ncol = length(n.seq))

cp_drct.10.runze.hbic=matrix(0, nrow = run, ncol = length(n.seq))
cp_drct.11.runze.hbic=matrix(0, nrow = run, ncol = length(n.seq))
cp_drct.20.runze.hbic=matrix(0, nrow = run, ncol = length(n.seq))
cp_drct.21.runze.hbic=matrix(0, nrow = run, ncol = length(n.seq))

reject_drct.10.runze.hbic=matrix(0, nrow = run, ncol = length(n.seq))
reject_drct.11.runze.hbic=matrix(0, nrow = run, ncol = length(n.seq))
reject_drct.20.runze.hbic=matrix(0, nrow = run, ncol = length(n.seq))
reject_drct.21.runze.hbic=matrix(0, nrow = run, ncol = length(n.seq))


NIE10 = as.numeric(t(beta1)%*%theta.G)
NIE11 = as.numeric(t(beta1)%*%theta.G) + as.numeric(t(beta1)%*%theta.D1)
NIE20 = as.numeric(t(beta2)%*%theta.G)
NIE21 = as.numeric(t(beta2)%*%theta.G) + as.numeric(t(beta2)%*%theta.D2)

NDE10 = theta.T[1] + t(beta0)%*%theta.D1 + t(mu.X)%*%t(A)%*%theta.D1
NDE11 = theta.T[1] + t(beta0)%*%theta.D1 + t(mu.X)%*%t(A)%*%theta.D1 + as.numeric(t(beta1)%*%theta.D1) 
NDE20 = theta.T[2] + t(beta0)%*%theta.D2 + t(mu.X)%*%t(A)%*%theta.D2
NDE21 = theta.T[2] + t(beta0)%*%theta.D2 + t(mu.X)%*%t(A)%*%theta.D2 + as.numeric(t(beta2)%*%theta.D2) 



for (loop.n in 1:length(n.seq)) {
  
  n=n.seq[loop.n]
  
  one.n=as.matrix(rep(1, times=n))
  
  TRT=matrix(0, nrow = n, ncol = q)
  TRT[,1]=c( rep(1, times=(n/3)), rep(0, times=(2*n/3))  )
  TRT[,2]=c( rep(0, times=(n/3)), rep(1, times=(n/3)), rep(0, times=(n/3))  )
  
  colnames(TRT) = c("T1", "T2")
  
  for (loop.run in 1:run) {
    
    # generate mediation and response
    
    White_G=matrix(rnorm(n*p), nrow = n, ncol = p)
    #White_G=matrix(runif(n*p, min = -sqrt(3), max = sqrt(3)), n,p)
    #White_G=matrix((rpois(n*p, lambda = 1)-1), n, p)
    #White_G=matrix((2*rbinom(n*p,1,0.5)-1), n,p)
    #White_G=matrix( (sqrt(3/5)*rt(n*p, df=5)) , nrow = n, ncol = p)
    #White_G=matrix( (sqrt(6/8)*rt(n*p, df=8)) , nrow = n, ncol = p)
    
    eps_G=White_G %*%Sigma0.5_G_gl
    
    X = matrix( rnorm(n*r, mean = 0, sd=1), n, r)
    X[,1] = rbinom(n, 1, 0.5 ) # the first column of X is from sex
    colnames(X) = paste0("X", 1:ncol(X))
    
    G=one.n%*%t(beta0)+TRT%*%t(B) + X%*%t(A) +eps_G # G is n by p matrix whose ith row is the ith obs of the mediation
    colnames(G) = paste0("G", 1:ncol(G))
    
    eps_Y=as.matrix(rnorm(n, mean = 0, sd=sigma_Y)) # normally distributed
    #eps_Y= sigma_Y * as.matrix(runif(n, min = -sqrt(3), max = sqrt(3))) # uniformly distributed
    #eps_Y = sigma_Y * as.matrix(rpois(n, lambda = 1)-1) # poisson centered
    #eps_Y=sigma_Y * as.matrix(2*rbinom(n,1,0.5)-1) # symmetric bernoulli
    #eps_Y=sigma_Y * as.matrix( sqrt(3/5)*rt(n, df=5) )
    #eps_Y=sigma_Y * as.matrix( sqrt(6/8)*rt(n, df=8) )
    
    
    G.T=matrix(0, nrow = n, ncol = q*p)
    G.T[1:(n/3), 1:p] = G[1:(n/3), ]
    G.T[((n/3)+1): (2*n/3), (p+1):(2*p)]=G[((n/3)+1): (2*n/3), ]
    colnames(G.T) = c( paste0("T1.G", 1:ncol(G)), paste0("T2.G", 1:ncol(G)) )
    
    Y=theta0*one.n + TRT%*%theta.T + X%*%theta.X + G%*%theta.G + G.T%*%as.matrix(vtheta.D) + eps_Y  # B is n by 1 outcome
    colnames(Y) = "Y"
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    ############################################### proposed ##########################################################################
    W_outcome=cbind(TRT, X, G, G.T)
    
    # max( abs( Y - cbind(one.n, W_outcome)%*%theta - eps_Y ) )
    
    fitp.scad=cv.ncvreg(W_outcome,Y, penalty = "SCAD", penalty.factor = p.fac ) # use default lambda
    
    ############################# the following annotations are about the output of cv.ncvreg() and ncvreg() functions #############################
    # lambda.seq = exp(seq(8,-20, length.out=100))
    
    # fitp.scad=cv.ncvreg(W_outcome,Y, penalty = "SCAD", penalty.factor = p.fac, lambda = lambda.seq )
    
    # fitp.scad$lambda - lambda.seq  # the outputs are all zeros, so the lambdas used in cv.ncvreg is exactly the argument "lambda = ".
    
    # fitp.scad$lambda.min - fitp.scad$lambda[fitp.scad$min] # fitp.scad$min is the number location (index) of fitp.scad$lambda which achives the smallest cve.
    
    # fitp.scad$cve # the red points (mean of cross-validation error under several (10) folds) in plot(fitp.scad)
    
    ## ncvreg(...) is the same as cv.ncvreg()$fit. To see this:
    # fit.whole = ncvreg(W_outcome,Y, penalty = "SCAD", penalty.factor = p.fac, lambda = lambda.seq ) # fit the whole data directly without cross-validation
    
    ## lambda are the same:
    # fitp.scad$fit$lambda - fit.whole$lambda # the results are all zero (they are using the same set of lambda)
    
    ## the coefficients are the same:
    # dim(coef(fit.whole)) # 3035  100  which is (1+q+r+p+p*q) (including the intercept) and length(lambda.seq)
    # coef(fitp.scad) # it produces the estimated theta's corresponding to lambda.min
    # dim(fitp.scad$fit$beta) # 3035  100 which is (1+q+r+p+p*q) (including the intercept) and length(lambda.seq)
    
    # dim( fitp.scad$fit$beta - coef(fit.whole) ) # 3035  100 which is (1+q+r+p+p*q) and length(lambda.seq)
    # min(abs( fitp.scad$fit$beta - coef(fit.whole) )) # the output is 0
    ## So we do not need to fit "fit.whole = ncvreg(...)" out of the whole data because fitp.scad = cv.ncvreg(...)$fit$beta has already produce the same coefficients from the whole data set
    
    
    ## To verify the statement above, you can also try the following codes to see the coefficients corresponding the lambda.min
    # max( abs( coef(fit.whole)[, fitp.scad$min] - coef(fitp.scad) ) ) # the output is 0
    # max( abs( coef(fit.whole)[, fitp.scad$min] - coef(fitp.scad, s = "lambda.min") ) ) # the output is 0
    # max( abs( coef(fit.whole)[, fitp.scad$min] - fitp.scad$fit$beta[, fitp.scad$min] ) ) # the output is 0, which means that coef(fitp.scad, s = "lambda.min") is calculated from the whole data set.
    
    ## One Caveat:
    # summary(fitp.scad)$sigma # This is not the sum of squared residuals/n, this is sqrt(fitp.scad$cve). See the next row:
    # summary(fitp.scad)$sigma - sqrt(fitp.scad$cve) # the outputs are all zeros
    
    ## So there is no direct result for the ( sum of squared residuals / n ) from each lambda
    ########################################################################################################
    
    
    
  
    ###### estimation by cv
    
    # step 1 outcome model
    
    lambda_scad.cv[loop.run, loop.n] = fitp.scad$lambda.min
    
    theta0.hat.cv=coef(fitp.scad, s = "lambda.min")[1]
    
    theta.T.hat.cv=coef(fitp.scad, s = "lambda.min")[(1+1): (1+q)]
    
    theta.X.hat.cv=coef(fitp.scad, s = "lambda.min")[(1+q+1) : (1+q+r)]
    
    theta.G.hat.cv=coef(fitp.scad, s = "lambda.min")[(1+q+r+1):(1+q+r+p)]
    
    theta.D1.hat.cv=coef(fitp.scad, s = "lambda.min")[(1+q+r+p+1):(1+q+r+p+p)]
    
    theta.D2.hat.cv=coef(fitp.scad, s = "lambda.min")[(1+q+r+p+p+1):(1+q+r+p+p+p)]
    
    vtheta.D_hat.cv=c(theta.D1.hat.cv, theta.D2.hat.cv)
    
    names(theta0.hat.cv)=NULL; names(theta.T.hat.cv)=NULL; names(theta.X.hat.cv)=NULL 
    names(theta.G.hat.cv)=NULL; names(vtheta.D_hat.cv)=NULL; names(theta.D1.hat.cv)=NULL; names(theta.D2.hat.cv)=NULL;
    
    bias_theta.scad.cv[loop.run,loop.n]=sum(abs(coef(fitp.scad, s = "lambda.min") - theta))
    sr_theta.scad.cv[loop.run,loop.n]=sum(sign(coef(fitp.scad, s = "lambda.min"))!=sign(theta))
    
    # check the error by different parts
    # sum(sign(coef(fitp.scad, s = "lambda.min"))!=sign(theta)) # support recovery
    # sum(abs(theta.T - theta.T.hat.cv))
    # sum(abs(theta.G.hat.cv - theta.G))
    # sum(abs(theta.X.hat.cv - theta.X))
    # sum(abs(theta.D1.hat.cv - theta.D1))
    # sum(abs(theta.D2.hat.cv - theta.D2))
    
    # step 2 mediation model
    
    W_mediator=cbind(one.n, TRT, X)
    
    # betahat = solve(t(W_mediator)%*%W_mediator) %*% t(W_mediator) %*% G
    
    # sum(abs(betahat - t(cbind(beta0,beta1, beta2, A) ) ))
    
    # beta.theta.G.hat.cv = betahat %*% theta.G.hat.cv
    
    # beta.theta.D1.hat.cv = betahat %*% theta.D1.hat.cv
    
    # beta.theta.D2.hat.cv = betahat %*% theta.D2.hat.cv
    
    # calculate the estimation error
    
    SE_smple.cv = OraProjInt_SE_sample_X(Y = Y, Trt = TRT, Mediator = G, Covariate = X, 
                                         Trt.Medator = G.T,
                                         theta0_hat = theta0.hat.cv,
                                         theta.T_hat = theta.T.hat.cv,
                                         theta.X_hat = theta.X.hat.cv,
                                         theta.G_hat = theta.G.hat.cv,
                                         vtheta.D_hat = vtheta.D_hat.cv)
   
    
    est_indrct.10.scad.cv[loop.run,loop.n]= SE_smple.cv$NIE_inference["NIE(T 1 = 0)", "Estimate"] 
    est_indrct.20.scad.cv[loop.run,loop.n]= SE_smple.cv$NIE_inference["NIE(T 2 = 0)", "Estimate"] 
    
    est_indrct.11.scad.cv[loop.run,loop.n]= SE_smple.cv$NIE_inference["NIE(T 1 = 1)", "Estimate"] 
    est_indrct.21.scad.cv[loop.run,loop.n]= SE_smple.cv$NIE_inference["NIE(T 2 = 1)", "Estimate"] 
    
    est_drct.10.scad.cv[loop.run,loop.n]= SE_smple.cv$NDE_inference["NDE(T 1 = 0)", "Estimate"] 
    est_drct.20.scad.cv[loop.run,loop.n]= SE_smple.cv$NDE_inference["NDE(T 2 = 0)", "Estimate"] 
    
    est_drct.11.scad.cv[loop.run,loop.n]= SE_smple.cv$NDE_inference["NDE(T 1 = 1)", "Estimate"] 
    est_drct.21.scad.cv[loop.run,loop.n]= SE_smple.cv$NDE_inference["NDE(T 2 = 1)", "Estimate"] 
    
    
    
    se_indrct.10.scad.cv[loop.run,loop.n]= SE_smple.cv$NIE_inference["NIE(T 1 = 0)", "SE"] 
    se_indrct.20.scad.cv[loop.run,loop.n]= SE_smple.cv$NIE_inference["NIE(T 2 = 0)", "SE"] 
    
    se_indrct.11.scad.cv[loop.run,loop.n]= SE_smple.cv$NIE_inference["NIE(T 1 = 1)", "SE"] 
    se_indrct.21.scad.cv[loop.run,loop.n]= SE_smple.cv$NIE_inference["NIE(T 2 = 1)", "SE"] 
    
    se_drct.10.scad.cv[loop.run,loop.n]= SE_smple.cv$NDE_inference["NDE(T 1 = 0)", "SE"] 
    se_drct.20.scad.cv[loop.run,loop.n]= SE_smple.cv$NDE_inference["NDE(T 2 = 0)", "SE"] 
    
    se_drct.11.scad.cv[loop.run,loop.n]= SE_smple.cv$NDE_inference["NDE(T 1 = 1)", "SE"] 
    se_drct.21.scad.cv[loop.run,loop.n]= SE_smple.cv$NDE_inference["NDE(T 2 = 1)", "SE"] 
    
    # cp
    if ( ( ( NIE10 ) <= ( est_indrct.10.scad.cv[loop.run,loop.n] + qnorm(0.975)*se_indrct.10.scad.cv[loop.run, loop.n] ) ) & ( ( NIE10 ) >= ( est_indrct.10.scad.cv[loop.run,loop.n] - qnorm(0.975)*se_indrct.10.scad.cv[loop.run, loop.n] ) )  ){
      cp_indrct.10.scad.cv[loop.run, loop.n] = 1
    }
    
    if ( ( ( NIE20) <= ( est_indrct.20.scad.cv[loop.run,loop.n] + qnorm(0.975)*se_indrct.20.scad.cv[loop.run, loop.n] ) ) & ( ( NIE20 ) >= ( est_indrct.20.scad.cv[loop.run,loop.n] - qnorm(0.975)*se_indrct.20.scad.cv[loop.run, loop.n] ) )  ){
      cp_indrct.20.scad.cv[loop.run, loop.n] = 1
    }
    
    
    if ( ( ( NIE11 ) <= ( est_indrct.11.scad.cv[loop.run,loop.n] + qnorm(0.975)*se_indrct.11.scad.cv[loop.run, loop.n] ) ) & ( ( NIE11 ) >= ( est_indrct.11.scad.cv[loop.run,loop.n] - qnorm(0.975)*se_indrct.11.scad.cv[loop.run, loop.n] ) )  ){
      cp_indrct.11.scad.cv[loop.run, loop.n] = 1
    }
    
    if ( ( ( NIE21 ) <= ( est_indrct.21.scad.cv[loop.run,loop.n] + qnorm(0.975)*se_indrct.21.scad.cv[loop.run, loop.n] ) ) & ( ( NIE21 ) >= ( est_indrct.21.scad.cv[loop.run,loop.n] - qnorm(0.975)*se_indrct.21.scad.cv[loop.run, loop.n] ) )  ){
      cp_indrct.21.scad.cv[loop.run, loop.n] = 1
    }
    
    #
    if ( ( ( NDE10 ) <= ( est_drct.10.scad.cv[loop.run,loop.n] + qnorm(0.975)*se_drct.10.scad.cv[loop.run, loop.n] ) ) & ( ( NDE10 ) >= ( est_drct.10.scad.cv[loop.run,loop.n] - qnorm(0.975)*se_drct.10.scad.cv[loop.run, loop.n] ) )  ){
      cp_drct.10.scad.cv[loop.run, loop.n] = 1
    }
    
    if ( ( ( NDE20 ) <= ( est_drct.20.scad.cv[loop.run,loop.n] + qnorm(0.975)*se_drct.20.scad.cv[loop.run, loop.n] ) ) & ( ( NDE20 ) >= ( est_drct.20.scad.cv[loop.run,loop.n] - qnorm(0.975)*se_drct.20.scad.cv[loop.run, loop.n] ) )  ){
      cp_drct.20.scad.cv[loop.run, loop.n] = 1
    }
    
    if ( ( ( NDE11 ) <= ( est_drct.11.scad.cv[loop.run,loop.n] + qnorm(0.975)*se_drct.11.scad.cv[loop.run, loop.n] ) ) & ( ( NDE11 ) >= ( est_drct.11.scad.cv[loop.run,loop.n] - qnorm(0.975)*se_drct.11.scad.cv[loop.run, loop.n] ) )  ){
      cp_drct.11.scad.cv[loop.run, loop.n] = 1
    }
    
    if ( ( ( NDE21 ) <= ( est_drct.21.scad.cv[loop.run,loop.n] + qnorm(0.975)*se_drct.21.scad.cv[loop.run, loop.n] ) ) & ( ( NDE21 ) >= ( est_drct.21.scad.cv[loop.run,loop.n] - qnorm(0.975)*se_drct.21.scad.cv[loop.run, loop.n] ) )  ){
      cp_drct.21.scad.cv[loop.run, loop.n] = 1
    }
    
    # power
    if ( ( 0 > ( est_indrct.10.scad.cv[loop.run,loop.n] + qnorm(0.975)*se_indrct.10.scad.cv[loop.run, loop.n] ) ) | ( 0 < ( est_indrct.10.scad.cv[loop.run,loop.n] - qnorm(0.975)*se_indrct.10.scad.cv[loop.run, loop.n] ) )  ){
      reject_indrct.10.scad.cv[loop.run, loop.n] = 1
    }
    
    if ( ( 0 > ( est_indrct.20.scad.cv[loop.run,loop.n] + qnorm(0.975)*se_indrct.20.scad.cv[loop.run, loop.n] ) ) | ( 0 < ( est_indrct.20.scad.cv[loop.run,loop.n] - qnorm(0.975)*se_indrct.20.scad.cv[loop.run, loop.n] ) )  ){
      reject_indrct.20.scad.cv[loop.run, loop.n] = 1
    }
    
    
    if ( ( 0 > ( est_indrct.11.scad.cv[loop.run,loop.n] + qnorm(0.975)*se_indrct.11.scad.cv[loop.run, loop.n] ) ) | ( 0 < ( est_indrct.11.scad.cv[loop.run,loop.n] - qnorm(0.975)*se_indrct.11.scad.cv[loop.run, loop.n] ) )  ){
      reject_indrct.11.scad.cv[loop.run, loop.n] = 1
    }
    
    if ( ( 0 > ( est_indrct.21.scad.cv[loop.run,loop.n] + qnorm(0.975)*se_indrct.21.scad.cv[loop.run, loop.n] ) ) | ( 0 < ( est_indrct.21.scad.cv[loop.run,loop.n] - qnorm(0.975)*se_indrct.21.scad.cv[loop.run, loop.n] ) )  ){
      reject_indrct.21.scad.cv[loop.run, loop.n] = 1
    }
    
    #
    if ( ( 0 > ( est_drct.10.scad.cv[loop.run,loop.n] + qnorm(0.975)*se_drct.10.scad.cv[loop.run, loop.n] ) ) | ( 0 < ( est_drct.10.scad.cv[loop.run,loop.n] - qnorm(0.975)*se_drct.10.scad.cv[loop.run, loop.n] ) )  ){
      reject_drct.10.scad.cv[loop.run, loop.n] = 1
    }
    
    if ( ( 0 > ( est_drct.20.scad.cv[loop.run,loop.n] + qnorm(0.975)*se_drct.20.scad.cv[loop.run, loop.n] ) ) | ( 0 < ( est_drct.20.scad.cv[loop.run,loop.n] - qnorm(0.975)*se_drct.20.scad.cv[loop.run, loop.n] ) )  ){
      reject_drct.20.scad.cv[loop.run, loop.n] = 1
    }
    
    if ( ( 0 > ( est_drct.11.scad.cv[loop.run,loop.n] + qnorm(0.975)*se_drct.11.scad.cv[loop.run, loop.n] ) ) | ( 0 < ( est_drct.11.scad.cv[loop.run,loop.n] - qnorm(0.975)*se_drct.11.scad.cv[loop.run, loop.n] ) )  ){
      reject_drct.11.scad.cv[loop.run, loop.n] = 1
    }
    
    if ( ( 0 > ( est_drct.21.scad.cv[loop.run,loop.n] + qnorm(0.975)*se_drct.21.scad.cv[loop.run, loop.n] ) ) | ( 0 < ( est_drct.21.scad.cv[loop.run,loop.n] - qnorm(0.975)*se_drct.21.scad.cv[loop.run, loop.n] ) )  ){
      reject_drct.21.scad.cv[loop.run, loop.n] = 1
    }
    
    
    ### HBIC ###
    
    # First calculate hbic corresponding to each lambda
    lambda.seq = fitp.scad$lambda
    
    hbic = c()
    for (loop.lambda in 1:length(lambda.seq)) {
      sigma2_Yhat.tmp =  sum( ( Y - cbind(one.n, W_outcome) %*% fitp.scad$fit$beta[,loop.lambda] )^2 ) / n 
      theta.G.tmp = (fitp.scad$fit$beta[,loop.lambda])[(1+q+r+1):(1+q+r+p)]
      vtheta.D.tmp = (fitp.scad$fit$beta[,loop.lambda])[(1+q+r+p+1):(1+q+r+p+p*q)]
      df = 1+q+r+length(which(theta.G.tmp!= 0))+length(which(vtheta.D.tmp!= 0))
      hbic[loop.lambda] = log( sigma2_Yhat.tmp ) +  df*log(log(n))*log(1+q+r+p+p*q)/n
    }
    
    ## aa is the log sum of squared residuals / n, which is decreasing linearly after it hits the lambda.min from cross-validation
    # aa = c()
    # for (loop.lambda in 1:length(lambda.seq)) {
    #   sigma2_Yhat.tmp =  sum( ( Y - cbind(one.n, W_outcome) %*% fitp.scad$fit$beta[,loop.lambda] )^2 ) / n
    #   aa[loop.lambda] = log( sigma2_Yhat.tmp )
    # }
    ## bb is the |M_lambda|_0 * Cn * log(1+q+r+p+p*q) / n in eqn(2.8) of the HBIC paper "CALIBRATING NONCONVEX PENALIZED REGRESSION IN ULTRA-HIGH DIMENSION"
    ## bb is increasing slower and slower
    # bb = c()
    # for (loop.lambda in 1:length(lambda.seq)) {
    #  theta.G.tmp = (fitp.scad$fit$beta[,loop.lambda])[(1+q+r+1):(1+q+r+p)]
    #  vtheta.D.tmp = (fitp.scad$fit$beta[,loop.lambda])[(1+q+r+p+1):(1+q+r+p+p*q)]
    #  df = 1+q+r+length(which(theta.G.tmp!= 0))+length(which(vtheta.D.tmp!= 0))
    #  bb[loop.lambda] = df*log(log(n))*log(1+q+r+p+p*q)/n
    # }
    ## since aa+bb=hbic, as lambda decreases (the entries in lambda.seq is decreasing), aa will dominate, so it tends to select a smaller lambda, if the range of lambda.seq is too large. So we use the lambda sequence generated by cross-validation.
    
    hbic.min = which.min(hbic)
    
    lambda_scad.hbic[loop.run, loop.n] = lambda.seq[hbic.min]
    
    theta0.hat.hbic=fitp.scad$fit$beta[,hbic.min][1]
    
    theta.T.hat.hbic=fitp.scad$fit$beta[,hbic.min][(1+1): (1+q)]
    
    theta.X.hat.hbic=fitp.scad$fit$beta[,hbic.min][(1+q+1) : (1+q+r)]
    
    theta.G.hat.hbic=fitp.scad$fit$beta[,hbic.min][(1+q+r+1):(1+q+r+p)]
    
    theta.D1.hat.hbic=fitp.scad$fit$beta[,hbic.min][(1+q+r+p+1):(1+q+r+p+p)]
    
    theta.D2.hat.hbic=fitp.scad$fit$beta[,hbic.min][(1+q+r+p+p+1):(1+q+r+p+p+p)]
    
    vtheta.D_hat.hbic=c(theta.D1.hat.hbic, theta.D2.hat.hbic)
    
    names(theta0.hat.hbic)=NULL; names(theta.T.hat.hbic)=NULL; names(theta.X.hat.hbic)=NULL 
    names(theta.G.hat.hbic)=NULL; names(vtheta.D_hat.hbic)=NULL; names(theta.D1.hat.hbic)=NULL; names(theta.D2.hat.hbic)=NULL;
    
    
    bias_theta.scad.hbic[loop.run,loop.n]=sum(abs(fitp.scad$fit$beta[,hbic.min] - theta))
    sr_theta.scad.hbic[loop.run,loop.n]=sum(sign(fitp.scad$fit$beta[,hbic.min])!=sign(theta))
    
    # step 2 mediation model
    
    W_mediator=cbind(one.n, TRT, X)
    
    # betahat = solve(t(W_mediator)%*%W_mediator) %*% t(W_mediator) %*% G
    
    # beta.theta.G.hat.hbic = betahat %*% theta.G.hat.hbic
    
    # beta.theta.D1.hat.hbic = betahat %*% theta.D1.hat.hbic
    
    # beta.theta.D2.hat.hbic = betahat %*% theta.D2.hat.hbic
    
    # calculate the estimation error
    
    SE_smple.hbic = OraProjInt_SE_sample_X(Y = Y, Trt = TRT, Mediator = G, Covariate = X, 
                                           Trt.Medator = G.T,
                                           theta0_hat = theta0.hat.hbic,
                                           theta.T_hat = theta.T.hat.hbic,
                                           theta.X_hat = theta.X.hat.hbic,
                                           theta.G_hat = theta.G.hat.hbic,
                                           vtheta.D_hat = vtheta.D_hat.hbic)
    
    
    est_indrct.10.scad.hbic[loop.run,loop.n]= SE_smple.hbic$NIE_inference["NIE(T 1 = 0)", "Estimate"] 
    est_indrct.20.scad.hbic[loop.run,loop.n]= SE_smple.hbic$NIE_inference["NIE(T 2 = 0)", "Estimate"] 
    
    est_indrct.11.scad.hbic[loop.run,loop.n]= SE_smple.hbic$NIE_inference["NIE(T 1 = 1)", "Estimate"] 
    est_indrct.21.scad.hbic[loop.run,loop.n]= SE_smple.hbic$NIE_inference["NIE(T 2 = 1)", "Estimate"] 
    
    est_drct.10.scad.hbic[loop.run,loop.n]= SE_smple.hbic$NDE_inference["NDE(T 1 = 0)", "Estimate"] 
    est_drct.20.scad.hbic[loop.run,loop.n]= SE_smple.hbic$NDE_inference["NDE(T 2 = 0)", "Estimate"] 
    
    est_drct.11.scad.hbic[loop.run,loop.n]= SE_smple.hbic$NDE_inference["NDE(T 1 = 1)", "Estimate"] 
    est_drct.21.scad.hbic[loop.run,loop.n]= SE_smple.hbic$NDE_inference["NDE(T 2 = 1)", "Estimate"] 
    
    
    
    se_indrct.10.scad.hbic[loop.run,loop.n]= SE_smple.hbic$NIE_inference["NIE(T 1 = 0)", "SE"] 
    se_indrct.20.scad.hbic[loop.run,loop.n]= SE_smple.hbic$NIE_inference["NIE(T 2 = 0)", "SE"] 
    
    se_indrct.11.scad.hbic[loop.run,loop.n]= SE_smple.hbic$NIE_inference["NIE(T 1 = 1)", "SE"] 
    se_indrct.21.scad.hbic[loop.run,loop.n]= SE_smple.hbic$NIE_inference["NIE(T 2 = 1)", "SE"] 
    
    se_drct.10.scad.hbic[loop.run,loop.n]= SE_smple.hbic$NDE_inference["NDE(T 1 = 0)", "SE"] 
    se_drct.20.scad.hbic[loop.run,loop.n]= SE_smple.hbic$NDE_inference["NDE(T 2 = 0)", "SE"] 
    
    se_drct.11.scad.hbic[loop.run,loop.n]= SE_smple.hbic$NDE_inference["NDE(T 1 = 1)", "SE"] 
    se_drct.21.scad.hbic[loop.run,loop.n]= SE_smple.hbic$NDE_inference["NDE(T 2 = 1)", "SE"] 
    
    # cp
    if ( ( ( NIE10 ) <= ( est_indrct.10.scad.hbic[loop.run,loop.n] + qnorm(0.975)*se_indrct.10.scad.hbic[loop.run, loop.n] ) ) & ( ( NIE10 ) >= ( est_indrct.10.scad.hbic[loop.run,loop.n] - qnorm(0.975)*se_indrct.10.scad.hbic[loop.run, loop.n] ) )  ){
      cp_indrct.10.scad.hbic[loop.run, loop.n] = 1
    }
    
    if ( ( ( NIE20) <= ( est_indrct.20.scad.hbic[loop.run,loop.n] + qnorm(0.975)*se_indrct.20.scad.hbic[loop.run, loop.n] ) ) & ( ( NIE20 ) >= ( est_indrct.20.scad.hbic[loop.run,loop.n] - qnorm(0.975)*se_indrct.20.scad.hbic[loop.run, loop.n] ) )  ){
      cp_indrct.20.scad.hbic[loop.run, loop.n] = 1
    }
    
    
    if ( ( ( NIE11 ) <= ( est_indrct.11.scad.hbic[loop.run,loop.n] + qnorm(0.975)*se_indrct.11.scad.hbic[loop.run, loop.n] ) ) & ( ( NIE11 ) >= ( est_indrct.11.scad.hbic[loop.run,loop.n] - qnorm(0.975)*se_indrct.11.scad.hbic[loop.run, loop.n] ) )  ){
      cp_indrct.11.scad.hbic[loop.run, loop.n] = 1
    }
    
    if ( ( ( NIE21 ) <= ( est_indrct.21.scad.hbic[loop.run,loop.n] + qnorm(0.975)*se_indrct.21.scad.hbic[loop.run, loop.n] ) ) & ( ( NIE21 ) >= ( est_indrct.21.scad.hbic[loop.run,loop.n] - qnorm(0.975)*se_indrct.21.scad.hbic[loop.run, loop.n] ) )  ){
      cp_indrct.21.scad.hbic[loop.run, loop.n] = 1
    }
    
    #
    if ( ( ( NDE10 ) <= ( est_drct.10.scad.hbic[loop.run,loop.n] + qnorm(0.975)*se_drct.10.scad.hbic[loop.run, loop.n] ) ) & ( ( NDE10 ) >= ( est_drct.10.scad.hbic[loop.run,loop.n] - qnorm(0.975)*se_drct.10.scad.hbic[loop.run, loop.n] ) )  ){
      cp_drct.10.scad.hbic[loop.run, loop.n] = 1
    }
    
    if ( ( ( NDE20 ) <= ( est_drct.20.scad.hbic[loop.run,loop.n] + qnorm(0.975)*se_drct.20.scad.hbic[loop.run, loop.n] ) ) & ( ( NDE20 ) >= ( est_drct.20.scad.hbic[loop.run,loop.n] - qnorm(0.975)*se_drct.20.scad.hbic[loop.run, loop.n] ) )  ){
      cp_drct.20.scad.hbic[loop.run, loop.n] = 1
    }
    
    if ( ( ( NDE11 ) <= ( est_drct.11.scad.hbic[loop.run,loop.n] + qnorm(0.975)*se_drct.11.scad.hbic[loop.run, loop.n] ) ) & ( ( NDE11 ) >= ( est_drct.11.scad.hbic[loop.run,loop.n] - qnorm(0.975)*se_drct.11.scad.hbic[loop.run, loop.n] ) )  ){
      cp_drct.11.scad.hbic[loop.run, loop.n] = 1
    }
    
    if ( ( ( NDE21 ) <= ( est_drct.21.scad.hbic[loop.run,loop.n] + qnorm(0.975)*se_drct.21.scad.hbic[loop.run, loop.n] ) ) & ( ( NDE21 ) >= ( est_drct.21.scad.hbic[loop.run,loop.n] - qnorm(0.975)*se_drct.21.scad.hbic[loop.run, loop.n] ) )  ){
      cp_drct.21.scad.hbic[loop.run, loop.n] = 1
    }
    
    # power
    if ( ( 0 > ( est_indrct.10.scad.hbic[loop.run,loop.n] + qnorm(0.975)*se_indrct.10.scad.hbic[loop.run, loop.n] ) ) | ( 0 < ( est_indrct.10.scad.hbic[loop.run,loop.n] - qnorm(0.975)*se_indrct.10.scad.hbic[loop.run, loop.n] ) )  ){
      reject_indrct.10.scad.hbic[loop.run, loop.n] = 1
    }
    
    if ( ( 0 > ( est_indrct.20.scad.hbic[loop.run,loop.n] + qnorm(0.975)*se_indrct.20.scad.hbic[loop.run, loop.n] ) ) | ( 0 < ( est_indrct.20.scad.hbic[loop.run,loop.n] - qnorm(0.975)*se_indrct.20.scad.hbic[loop.run, loop.n] ) )  ){
      reject_indrct.20.scad.hbic[loop.run, loop.n] = 1
    }
    
    
    if ( ( 0 > ( est_indrct.11.scad.hbic[loop.run,loop.n] + qnorm(0.975)*se_indrct.11.scad.hbic[loop.run, loop.n] ) ) | ( 0 < ( est_indrct.11.scad.hbic[loop.run,loop.n] - qnorm(0.975)*se_indrct.11.scad.hbic[loop.run, loop.n] ) )  ){
      reject_indrct.11.scad.hbic[loop.run, loop.n] = 1
    }
    
    if ( ( 0 > ( est_indrct.21.scad.hbic[loop.run,loop.n] + qnorm(0.975)*se_indrct.21.scad.hbic[loop.run, loop.n] ) ) | ( 0 < ( est_indrct.21.scad.hbic[loop.run,loop.n] - qnorm(0.975)*se_indrct.21.scad.hbic[loop.run, loop.n] ) )  ){
      reject_indrct.21.scad.hbic[loop.run, loop.n] = 1
    }
    
    #
    if ( ( 0 > ( est_drct.10.scad.hbic[loop.run,loop.n] + qnorm(0.975)*se_drct.10.scad.hbic[loop.run, loop.n] ) ) | ( 0 < ( est_drct.10.scad.hbic[loop.run,loop.n] - qnorm(0.975)*se_drct.10.scad.hbic[loop.run, loop.n] ) )  ){
      reject_drct.10.scad.hbic[loop.run, loop.n] = 1
    }
    
    if ( ( 0 > ( est_drct.20.scad.hbic[loop.run,loop.n] + qnorm(0.975)*se_drct.20.scad.hbic[loop.run, loop.n] ) ) | ( 0 < ( est_drct.20.scad.hbic[loop.run,loop.n] - qnorm(0.975)*se_drct.20.scad.hbic[loop.run, loop.n] ) )  ){
      reject_drct.20.scad.hbic[loop.run, loop.n] = 1
    }
    
    if ( ( 0 > ( est_drct.11.scad.hbic[loop.run,loop.n] + qnorm(0.975)*se_drct.11.scad.hbic[loop.run, loop.n] ) ) | ( 0 < ( est_drct.11.scad.hbic[loop.run,loop.n] - qnorm(0.975)*se_drct.11.scad.hbic[loop.run, loop.n] ) )  ){
      reject_drct.11.scad.hbic[loop.run, loop.n] = 1
    }
    
    if ( ( 0 > ( est_drct.21.scad.hbic[loop.run,loop.n] + qnorm(0.975)*se_drct.21.scad.hbic[loop.run, loop.n] ) ) | ( 0 < ( est_drct.21.scad.hbic[loop.run,loop.n] - qnorm(0.975)*se_drct.21.scad.hbic[loop.run, loop.n] ) )  ){
      reject_drct.21.scad.hbic[loop.run, loop.n] = 1
    }
    #########################################################################################################################
    
    
    
    
    
    
    
    
    
    
    
    
    ############### runze ########################################################################
    W_outcome_runze=cbind(TRT, X, G)
    
    # max( abs( Y - cbind(one.n, W_outcome)%*%theta - eps_Y ) )
    
    fitp.runze=cv.ncvreg(W_outcome_runze, Y, penalty = "SCAD", penalty.factor = p.fac.runze ) # use default lambda
    # plot(fitp.scad)

    
    ###### estimation by cv
    
    # step 1 outcome model
    
    theta0.hat.runze.cv=coef(fitp.runze, s = "lambda.min")[1]
    
    theta.T.hat.runze.cv=coef(fitp.runze, s = "lambda.min")[(1+1): (1+q)]
    
    theta.X.hat.runze.cv=coef(fitp.runze, s = "lambda.min")[(1+q+1) : (1+q+r)]
    
    theta.G.hat.runze.cv=coef(fitp.runze, s = "lambda.min")[(1+q+r+1):(1+q+r+p)]
    
    names(theta0.hat.runze.cv)=NULL; names(theta.T.hat.runze.cv)=NULL; names(theta.X.hat.runze.cv)=NULL 
    names(theta.G.hat.runze.cv)=NULL
    
    # step 2 mediation model
    
    SE_smple.runze.cv = Runze_SE_sample_noX(Y = Y, Trt = cbind(TRT,X), Mediator = G, 
                                      theta0_hat = theta0.hat.runze.cv,
                                      theta.T_hat = c(theta.T.hat.runze.cv, theta.X.hat.runze.cv),
                                      theta.G_hat = theta.G.hat.runze.cv)
    
    
    est_indrct.1.runze.cv[loop.run,loop.n] = SE_smple.runze.cv$NIE_est["T1", 1] 
    est_indrct.2.runze.cv[loop.run,loop.n] = SE_smple.runze.cv$NIE_est["T2", 1]
    
    se_indrct.1.runze.cv[loop.run,loop.n]= SE_smple.runze.cv$NIE_SE["T1", 1] 
    se_indrct.2.runze.cv[loop.run,loop.n]= SE_smple.runze.cv$NIE_SE["T2", 1] 
    
    est_drct.1.runze.cv[loop.run, loop.n] = SE_smple.runze.cv$NDE_est["T1", 1]
    est_drct.2.runze.cv[loop.run, loop.n] = SE_smple.runze.cv$NDE_est["T2", 1]
    
    se_drct.1.runze.cv[loop.run, loop.n] = SE_smple.runze.cv$NDE_SE["T1", 1]
    se_drct.1.runze.cv[loop.run, loop.n] = SE_smple.runze.cv$NDE_SE["T2", 1]
 
    # cp
    if ( ( ( NIE10 ) <= ( est_indrct.1.runze.cv[loop.run,loop.n] + qnorm(0.975)*se_indrct.1.runze.cv[loop.run, loop.n] ) ) & ( ( NIE10 ) >= ( est_indrct.1.runze.cv[loop.run,loop.n] - qnorm(0.975)*se_indrct.1.runze.cv[loop.run, loop.n] ) )  ){
      cp_indrct.10.runze.cv[loop.run, loop.n] = 1
    }
    
    if ( ( ( NIE20) <= ( est_indrct.2.runze.cv[loop.run,loop.n] + qnorm(0.975)*se_indrct.2.runze.cv[loop.run, loop.n] ) ) & ( ( NIE20 ) >= ( est_indrct.2.runze.cv[loop.run,loop.n] - qnorm(0.975)*se_indrct.2.runze.cv[loop.run, loop.n] ) )  ){
      cp_indrct.20.runze.cv[loop.run, loop.n] = 1
    }
    
    
    if ( ( ( NIE11 ) <= ( est_indrct.1.runze.cv[loop.run,loop.n] + qnorm(0.975)*se_indrct.1.runze.cv[loop.run, loop.n] ) ) & ( ( NIE11 ) >= ( est_indrct.1.runze.cv[loop.run,loop.n] - qnorm(0.975)*se_indrct.1.runze.cv[loop.run, loop.n] ) )  ){
      cp_indrct.11.runze.cv[loop.run, loop.n] = 1
    }
    
    if ( ( ( NIE21 ) <= ( est_indrct.2.runze.cv[loop.run,loop.n] + qnorm(0.975)*se_indrct.2.runze.cv[loop.run, loop.n] ) ) & ( ( NIE21 ) >= ( est_indrct.2.runze.cv[loop.run,loop.n] - qnorm(0.975)*se_indrct.2.runze.cv[loop.run, loop.n] ) )  ){
      cp_indrct.21.runze.cv[loop.run, loop.n] = 1
    }
    
    #
    if ( ( ( NDE10 ) <= ( est_drct.1.runze.cv[loop.run,loop.n] + qnorm(0.975)*se_drct.1.runze.cv[loop.run, loop.n] ) ) & ( ( NDE10 ) >= ( est_drct.1.runze.cv[loop.run,loop.n] - qnorm(0.975)*se_drct.1.runze.cv[loop.run, loop.n] ) )  ){
      cp_drct.10.runze.cv[loop.run, loop.n] = 1
    }
    
    if ( ( ( NDE20 ) <= ( est_drct.2.runze.cv[loop.run,loop.n] + qnorm(0.975)*se_drct.2.runze.cv[loop.run, loop.n] ) ) & ( ( NDE20 ) >= ( est_drct.2.runze.cv[loop.run,loop.n] - qnorm(0.975)*se_drct.2.runze.cv[loop.run, loop.n] ) )  ){
      cp_drct.20.runze.cv[loop.run, loop.n] = 1
    }
    
    if ( ( ( NDE11 ) <= ( est_drct.1.runze.cv[loop.run,loop.n] + qnorm(0.975)*se_drct.1.runze.cv[loop.run, loop.n] ) ) & ( ( NDE11 ) >= ( est_drct.1.runze.cv[loop.run,loop.n] - qnorm(0.975)*se_drct.1.runze.cv[loop.run, loop.n] ) )  ){
      cp_drct.11.runze.cv[loop.run, loop.n] = 1
    }
    
    if ( ( ( NDE21 ) <= ( est_drct.2.runze.cv[loop.run,loop.n] + qnorm(0.975)*se_drct.2.runze.cv[loop.run, loop.n] ) ) & ( ( NDE21 ) >= ( est_drct.2.runze.cv[loop.run,loop.n] - qnorm(0.975)*se_drct.2.runze.cv[loop.run, loop.n] ) )  ){
      cp_drct.21.runze.cv[loop.run, loop.n] = 1
    }
    
    # power
    if ( ( 0 > ( est_indrct.1.runze.cv[loop.run,loop.n] + qnorm(0.975)*se_indrct.1.runze.cv[loop.run, loop.n] ) ) | ( 0 < ( est_indrct.1.runze.cv[loop.run,loop.n] - qnorm(0.975)*se_indrct.1.runze.cv[loop.run, loop.n] ) )  ){
      reject_indrct.10.runze.cv[loop.run, loop.n] = 1
      reject_indrct.11.runze.cv[loop.run, loop.n] = 1
    }
    
    if ( ( 0 > ( est_indrct.2.runze.cv[loop.run,loop.n] + qnorm(0.975)*se_indrct.2.runze.cv[loop.run, loop.n] ) ) | ( 0 < ( est_indrct.2.runze.cv[loop.run,loop.n] - qnorm(0.975)*se_indrct.2.runze.cv[loop.run, loop.n] ) )  ){
      reject_indrct.20.runze.cv[loop.run, loop.n] = 1
      reject_indrct.21.runze.cv[loop.run, loop.n] = 1
    }
    
    #
    if ( ( 0 > ( est_drct.1.runze.cv[loop.run,loop.n] + qnorm(0.975)*se_drct.1.runze.cv[loop.run, loop.n] ) ) | ( 0 < ( est_drct.1.runze.cv[loop.run,loop.n] - qnorm(0.975)*se_drct.1.runze.cv[loop.run, loop.n] ) )  ){
      reject_drct.10.runze.cv[loop.run, loop.n] = 1
      reject_drct.11.runze.cv[loop.run, loop.n] = 1
    }
    
    if ( ( 0 > ( est_drct.2.runze.cv[loop.run,loop.n] + qnorm(0.975)*se_drct.2.runze.cv[loop.run, loop.n] ) ) | ( 0 < ( est_drct.2.runze.cv[loop.run,loop.n] - qnorm(0.975)*se_drct.2.runze.cv[loop.run, loop.n] ) )  ){
      reject_drct.20.runze.cv[loop.run, loop.n] = 1
      reject_drct.21.runze.cv[loop.run, loop.n] = 1
    }
    
    
    
    
    ### HBIC ###
    
    # First calculate hbic corresponding to each lambda
    lambda.seq.runze = fitp.runze$lambda
    
    hbic.runze = c()
    for (loop.lambda in 1:length(lambda.seq.runze)) {
      sigma2_Yhat.tmp =  sum( ( Y - cbind(one.n, W_outcome_runze) %*% fitp.runze$fit$beta[,loop.lambda] )^2 ) / n 
      theta.G.tmp = (fitp.runze$fit$beta[,loop.lambda])[(1+q+r+1):(1+q+r+p)]
      df = 1+q+r+length(which(theta.G.tmp!= 0))
      hbic.runze[loop.lambda] = log( sigma2_Yhat.tmp ) +  df*log(log(n))*log(1+q+r+p)/n
    }
    
    
    hbic.min.runze = which.min(hbic.runze)
    
    theta0.hat.runze.hbic=fitp.runze$fit$beta[,hbic.min.runze][1]
    
    theta.T.hat.runze.hbic=fitp.runze$fit$beta[,hbic.min.runze][(1+1): (1+q)]
    
    theta.X.hat.runze.hbic=fitp.runze$fit$beta[,hbic.min.runze][(1+q+1) : (1+q+r)]
    
    theta.G.hat.runze.hbic=fitp.runze$fit$beta[,hbic.min.runze][(1+q+r+1):(1+q+r+p)]
    
    names(theta0.hat.runze.hbic)=NULL; names(theta.T.hat.runze.hbic)=NULL; names(theta.X.hat.runze.hbic)=NULL 
    names(theta.G.hat.runze.hbic)=NULL
    
    # step 2 mediation model
    
    SE_smple.runze.hbic = Runze_SE_sample_noX(Y = Y, Trt = cbind(TRT,X), Mediator = G, 
                                            theta0_hat = theta0.hat.runze.hbic,
                                            theta.T_hat = c(theta.T.hat.runze.hbic, theta.X.hat.runze.hbic),
                                            theta.G_hat = theta.G.hat.runze.hbic)
    
    
    est_indrct.1.runze.hbic[loop.run,loop.n] = SE_smple.runze.hbic$NIE_est["T1", 1] 
    est_indrct.2.runze.hbic[loop.run,loop.n] = SE_smple.runze.hbic$NIE_est["T2", 1]
    
    se_indrct.1.runze.hbic[loop.run,loop.n]= SE_smple.runze.hbic$NIE_SE["T1", 1] 
    se_indrct.2.runze.hbic[loop.run,loop.n]= SE_smple.runze.hbic$NIE_SE["T2", 1] 
    
    est_drct.1.runze.hbic[loop.run, loop.n] = SE_smple.runze.hbic$NDE_est["T1", 1]
    est_drct.2.runze.hbic[loop.run, loop.n] = SE_smple.runze.hbic$NDE_est["T2", 1]
    
    se_drct.1.runze.hbic[loop.run, loop.n] = SE_smple.runze.hbic$NDE_SE["T1", 1]
    se_drct.1.runze.hbic[loop.run, loop.n] = SE_smple.runze.hbic$NDE_SE["T2", 1]
    
    # cp
    if ( ( ( NIE10 ) <= ( est_indrct.1.runze.hbic[loop.run,loop.n] + qnorm(0.975)*se_indrct.1.runze.hbic[loop.run, loop.n] ) ) & ( ( NIE10 ) >= ( est_indrct.1.runze.hbic[loop.run,loop.n] - qnorm(0.975)*se_indrct.1.runze.hbic[loop.run, loop.n] ) )  ){
      cp_indrct.10.runze.hbic[loop.run, loop.n] = 1
    }
    
    if ( ( ( NIE20) <= ( est_indrct.2.runze.hbic[loop.run,loop.n] + qnorm(0.975)*se_indrct.2.runze.hbic[loop.run, loop.n] ) ) & ( ( NIE20 ) >= ( est_indrct.2.runze.hbic[loop.run,loop.n] - qnorm(0.975)*se_indrct.2.runze.hbic[loop.run, loop.n] ) )  ){
      cp_indrct.20.runze.hbic[loop.run, loop.n] = 1
    }
    
    
    if ( ( ( NIE11 ) <= ( est_indrct.1.runze.hbic[loop.run,loop.n] + qnorm(0.975)*se_indrct.1.runze.hbic[loop.run, loop.n] ) ) & ( ( NIE11 ) >= ( est_indrct.1.runze.hbic[loop.run,loop.n] - qnorm(0.975)*se_indrct.1.runze.hbic[loop.run, loop.n] ) )  ){
      cp_indrct.11.runze.hbic[loop.run, loop.n] = 1
    }
    
    if ( ( ( NIE21 ) <= ( est_indrct.2.runze.hbic[loop.run,loop.n] + qnorm(0.975)*se_indrct.2.runze.hbic[loop.run, loop.n] ) ) & ( ( NIE21 ) >= ( est_indrct.2.runze.hbic[loop.run,loop.n] - qnorm(0.975)*se_indrct.2.runze.hbic[loop.run, loop.n] ) )  ){
      cp_indrct.21.runze.hbic[loop.run, loop.n] = 1
    }
    
    #
    if ( ( ( NDE10 ) <= ( est_drct.1.runze.hbic[loop.run,loop.n] + qnorm(0.975)*se_drct.1.runze.hbic[loop.run, loop.n] ) ) & ( ( NDE10 ) >= ( est_drct.1.runze.hbic[loop.run,loop.n] - qnorm(0.975)*se_drct.1.runze.hbic[loop.run, loop.n] ) )  ){
      cp_drct.10.runze.hbic[loop.run, loop.n] = 1
    }
    
    if ( ( ( NDE20 ) <= ( est_drct.2.runze.hbic[loop.run,loop.n] + qnorm(0.975)*se_drct.2.runze.hbic[loop.run, loop.n] ) ) & ( ( NDE20 ) >= ( est_drct.2.runze.hbic[loop.run,loop.n] - qnorm(0.975)*se_drct.2.runze.hbic[loop.run, loop.n] ) )  ){
      cp_drct.20.runze.hbic[loop.run, loop.n] = 1
    }
    
    if ( ( ( NDE11 ) <= ( est_drct.1.runze.hbic[loop.run,loop.n] + qnorm(0.975)*se_drct.1.runze.hbic[loop.run, loop.n] ) ) & ( ( NDE11 ) >= ( est_drct.1.runze.hbic[loop.run,loop.n] - qnorm(0.975)*se_drct.1.runze.hbic[loop.run, loop.n] ) )  ){
      cp_drct.11.runze.hbic[loop.run, loop.n] = 1
    }
    
    if ( ( ( NDE21 ) <= ( est_drct.2.runze.hbic[loop.run,loop.n] + qnorm(0.975)*se_drct.2.runze.hbic[loop.run, loop.n] ) ) & ( ( NDE21 ) >= ( est_drct.2.runze.hbic[loop.run,loop.n] - qnorm(0.975)*se_drct.2.runze.hbic[loop.run, loop.n] ) )  ){
      cp_drct.21.runze.hbic[loop.run, loop.n] = 1
    }
    
    # power
    if ( ( 0 > ( est_indrct.1.runze.hbic[loop.run,loop.n] + qnorm(0.975)*se_indrct.1.runze.hbic[loop.run, loop.n] ) ) | ( 0 < ( est_indrct.1.runze.hbic[loop.run,loop.n] - qnorm(0.975)*se_indrct.1.runze.hbic[loop.run, loop.n] ) )  ){
      reject_indrct.10.runze.hbic[loop.run, loop.n] = 1
      reject_indrct.11.runze.hbic[loop.run, loop.n] = 1
    }
    
    if ( ( 0 > ( est_indrct.2.runze.hbic[loop.run,loop.n] + qnorm(0.975)*se_indrct.2.runze.hbic[loop.run, loop.n] ) ) | ( 0 < ( est_indrct.2.runze.hbic[loop.run,loop.n] - qnorm(0.975)*se_indrct.2.runze.hbic[loop.run, loop.n] ) )  ){
      reject_indrct.20.runze.hbic[loop.run, loop.n] = 1
      reject_indrct.21.runze.hbic[loop.run, loop.n] = 1
    }
    
    #
    if ( ( 0 > ( est_drct.1.runze.hbic[loop.run,loop.n] + qnorm(0.975)*se_drct.1.runze.hbic[loop.run, loop.n] ) ) | ( 0 < ( est_drct.1.runze.hbic[loop.run,loop.n] - qnorm(0.975)*se_drct.1.runze.hbic[loop.run, loop.n] ) )  ){
      reject_drct.10.runze.hbic[loop.run, loop.n] = 1
      reject_drct.11.runze.hbic[loop.run, loop.n] = 1
    }
    
    if ( ( 0 > ( est_drct.2.runze.hbic[loop.run,loop.n] + qnorm(0.975)*se_drct.2.runze.hbic[loop.run, loop.n] ) ) | ( 0 < ( est_drct.2.runze.hbic[loop.run,loop.n] - qnorm(0.975)*se_drct.2.runze.hbic[loop.run, loop.n] ) )  ){
      reject_drct.20.runze.hbic[loop.run, loop.n] = 1
      reject_drct.21.runze.hbic[loop.run, loop.n] = 1
    }
    
    
    
    print(loop.run)
    
  }
  print(loop.n)
}


plot(lambda_scad.cv[,1], lambda_scad.hbic[,1])
abline(a=0,b=1,col=2)


cp_pwr_cv = matrix(NA, 8, 4)
colnames(cp_pwr_cv) = c("cp proposed", "cp Runze", "power proposed", "power Runze")
rownames(cp_pwr_cv) = c("NIE10", "NIE11", "NIE20", "NIE21", "NDE10", "NDE11", "NDE20", "NDE21")
cp_pwr_cv["NIE10", ] = c( mean(cp_indrct.10.scad.cv), mean(cp_indrct.10.runze.cv), 
                       mean(reject_indrct.10.scad.cv), mean(reject_indrct.10.runze.cv)  )
cp_pwr_cv["NIE11", ] = c( mean(cp_indrct.11.scad.cv), mean(cp_indrct.11.runze.cv), 
                       mean(reject_indrct.11.scad.cv), mean(reject_indrct.11.runze.cv)  )
cp_pwr_cv["NIE20", ] = c( mean(cp_indrct.20.scad.cv), mean(cp_indrct.20.runze.cv), 
                       mean(reject_indrct.20.scad.cv), mean(reject_indrct.20.runze.cv)  )
cp_pwr_cv["NIE21", ] = c( mean(cp_indrct.21.scad.cv), mean(cp_indrct.21.runze.cv), 
                       mean(reject_indrct.21.scad.cv), mean(reject_indrct.21.runze.cv)  )
cp_pwr_cv["NDE10", ] = c( mean(cp_drct.10.scad.cv), mean(cp_drct.10.runze.cv), 
                       mean(reject_drct.10.scad.cv), mean(reject_drct.10.runze.cv)  )
cp_pwr_cv["NDE11", ] = c( mean(cp_drct.11.scad.cv), mean(cp_drct.11.runze.cv), 
                       mean(reject_drct.11.scad.cv), mean(reject_drct.11.runze.cv)  )
cp_pwr_cv["NDE20", ] = c( mean(cp_drct.20.scad.cv), mean(cp_drct.20.runze.cv), 
                       mean(reject_drct.20.scad.cv), mean(reject_drct.20.runze.cv)  )
cp_pwr_cv["NDE21", ] = c( mean(cp_drct.21.scad.cv), mean(cp_drct.21.runze.cv), 
                       mean(reject_drct.21.scad.cv), mean(reject_drct.21.runze.cv)  )


cp_pwr_cv



cp_pwr_hbic = matrix(NA, 8, 4)
colnames(cp_pwr_hbic) = c("cp proposed", "cp Runze", "power proposed", "power Runze")
rownames(cp_pwr_hbic) = c("NIE10", "NIE11", "NIE20", "NIE21", "NDE10", "NDE11", "NDE20", "NDE21")
cp_pwr_hbic["NIE10", ] = c( mean(cp_indrct.10.scad.hbic), mean(cp_indrct.10.runze.hbic), 
                          mean(reject_indrct.10.scad.hbic), mean(reject_indrct.10.runze.hbic)  )
cp_pwr_hbic["NIE11", ] = c( mean(cp_indrct.11.scad.hbic), mean(cp_indrct.11.runze.hbic), 
                          mean(reject_indrct.11.scad.hbic), mean(reject_indrct.11.runze.hbic)  )
cp_pwr_hbic["NIE20", ] = c( mean(cp_indrct.20.scad.hbic), mean(cp_indrct.20.runze.hbic), 
                          mean(reject_indrct.20.scad.hbic), mean(reject_indrct.20.runze.hbic)  )
cp_pwr_hbic["NIE21", ] = c( mean(cp_indrct.21.scad.hbic), mean(cp_indrct.21.runze.hbic), 
                          mean(reject_indrct.21.scad.hbic), mean(reject_indrct.21.runze.hbic)  )
cp_pwr_hbic["NDE10", ] = c( mean(cp_drct.10.scad.hbic), mean(cp_drct.10.runze.hbic), 
                          mean(reject_drct.10.scad.hbic), mean(reject_drct.10.runze.hbic)  )
cp_pwr_hbic["NDE11", ] = c( mean(cp_drct.11.scad.hbic), mean(cp_drct.11.runze.hbic), 
                          mean(reject_drct.11.scad.hbic), mean(reject_drct.11.runze.hbic)  )
cp_pwr_hbic["NDE20", ] = c( mean(cp_drct.20.scad.hbic), mean(cp_drct.20.runze.hbic), 
                          mean(reject_drct.20.scad.hbic), mean(reject_drct.20.runze.hbic)  )
cp_pwr_hbic["NDE21", ] = c( mean(cp_drct.21.scad.hbic), mean(cp_drct.21.runze.hbic), 
                          mean(reject_drct.21.scad.hbic), mean(reject_drct.21.runze.hbic)  )


cp_pwr_hbic

round(cp_pwr_hbic, digits = 3)




save(list = c("A", "B", "beta0", "beta1", "beta2", "cp_pwr_cv", "cp_pwr_hbic",
              "lambda_scad.cv", "lambda_scad.hbic", "mu.X", "NDE10", "NDE11", 
              "NDE20", "NDE21", "Sigma_G_gl", "sigma.X", "A_D1", "A_D2", "A_G", "n.seq", 
              "NIE10", "NIE11", "NIE20", "NIE21", "p", "p.fac", "p.fac.runze", "q", "r", "run", "sigma_Y",
              "theta", "theta.D1", "theta.D2", "theta.G", "theta.T", "theta.X", "theta0", "vtheta.D"),
     file = "p=1000_n=120_glasso_NormalepsG_NormalepsY.Rdata")





















############################################## Then try p=2000 ################################
rm(list = ls())

p=2000
source("MyMethod_functions.R")


### glasso Sigma_G ###
load("Glasso_Sigma_G_p=1000.Rdata") # a 1000 by 1000 matrix (need adjustment when p=2000)
Sigma_G_gl = bdiag(Sigma_G_gl, Sigma_G_gl)


###### AR1 Sigma_G (do not run this if you use glasso Sigma_G) #####
rho_G=0.5 # if rho_G=0, Sigma_G=I_p * c_G
c_sigmaG=1 # c_G is the magnitude of the noise
# AR1
Gamma_G=matrix(0, p, p)
for (i in 1:p) {
  for (j in 1:p) {
    if(i == j ){
      Gamma_G[i,j] = 1
    } else {
      Gamma_G[i,j] = rho_G^(abs(i-j))
    }
  }
}
Sigma_G_AR1=c_sigmaG*Gamma_G
rm(Gamma_G, c_sigmaG)
############################################################






#### find glasso Sigma)_G^{1/2}
Q_G=eigen(Sigma_G_gl)$vectors
Lambda0.5_G=diag(c(sqrt(eigen(Sigma_G_gl)$values)))
Sigma0.5_G_gl=Q_G%*%Lambda0.5_G%*%t(Q_G)
# check whether it is correct
print(max(abs(Sigma0.5_G_gl%*%Sigma0.5_G_gl-Sigma_G_gl)))
rm(Q_G, Lambda0.5_G, Sigma_G_gl)



run = 500
n.seq=c(120)


# define the paramters
r=2
q=2

theta0=-2


theta.T=c(1.5, -0.5)

A_G = c(20, 27, 106,  146,  149,  268 , 497,  589,  807,  879)
A_D1 = c(106, 879)
A_D2 = c(807)

theta.G = rep(0, times=p)
theta.G = rep(0, times=p)
theta.G[A_G] = c(-3, 2, 1,-1, -2, 4, -2, 4, 1, -1)

vtheta.D = rep(0, times=(2*p))
vtheta.D[A_D1] = c(5, 4)
vtheta.D[p+A_D2] = c(5)
which(vtheta.D!=0)


theta.X = c(3, -3)

theta = c(theta0, theta.T, theta.X, theta.G, vtheta.D)
theta.D1 = vtheta.D[1:p]
theta.D2 = vtheta.D[(p+1):(p+p)]

sigma_Y = 1

set.seed(1234)

beta0 = matrix( rnorm(p, sd = 0.1) , ncol=1)
beta0[A_G,1] = seq(0.1, 1, by=0.1)
plot(beta0)

beta1 = matrix( rnorm(p, sd = 0.1) , ncol=1)
beta1[A_G,1] = (-1)*seq(0.05, 0.5, by=0.05)

beta2 = matrix( rnorm(p, sd = 0.1) , ncol=1)
beta2[A_G,1] = (-1)*seq(0.05, 0.5, by=0.05)

B = cbind(beta1, beta2)

mu.X = matrix(c(0.5,0), nrow=2)
sigma.X = matrix(c(0.25,0,0,1), 2,2)
A = matrix(rnorm((p*r), sd=0.1), nrow = p, ncol = r)


p.fac<- rep(1, (q+r+p+q*p))
p.fac[1:(q+r)] <- 0

p.fac.runze<- rep(1, (q+r+p))
p.fac.runze[1:(q+r)] <- 0




######## proposed cv
bias_theta.scad.cv=matrix(0, nrow = run, ncol = length(n.seq))
sr_theta.scad.cv=matrix(0, nrow = run, ncol = length(n.seq))

est_indrct.10.scad.cv=matrix(0, nrow = run, ncol = length(n.seq))
est_indrct.11.scad.cv=matrix(0, nrow = run, ncol = length(n.seq))
est_indrct.20.scad.cv=matrix(0, nrow = run, ncol = length(n.seq))
est_indrct.21.scad.cv=matrix(0, nrow = run, ncol = length(n.seq))

est_drct.10.scad.cv=matrix(0, nrow = run, ncol = length(n.seq))
est_drct.11.scad.cv=matrix(0, nrow = run, ncol = length(n.seq))
est_drct.20.scad.cv=matrix(0, nrow = run, ncol = length(n.seq))
est_drct.21.scad.cv=matrix(0, nrow = run, ncol = length(n.seq))

se_indrct.10.scad.cv=matrix(0, nrow = run, ncol = length(n.seq))
se_indrct.11.scad.cv=matrix(0, nrow = run, ncol = length(n.seq))
se_indrct.20.scad.cv=matrix(0, nrow = run, ncol = length(n.seq))
se_indrct.21.scad.cv=matrix(0, nrow = run, ncol = length(n.seq))

se_drct.10.scad.cv=matrix(0, nrow = run, ncol = length(n.seq))
se_drct.11.scad.cv=matrix(0, nrow = run, ncol = length(n.seq))
se_drct.20.scad.cv=matrix(0, nrow = run, ncol = length(n.seq))
se_drct.21.scad.cv=matrix(0, nrow = run, ncol = length(n.seq))

cp_indrct.10.scad.cv=matrix(0, nrow = run, ncol = length(n.seq))
cp_indrct.11.scad.cv=matrix(0, nrow = run, ncol = length(n.seq))
cp_indrct.20.scad.cv=matrix(0, nrow = run, ncol = length(n.seq))
cp_indrct.21.scad.cv=matrix(0, nrow = run, ncol = length(n.seq))

cp_drct.10.scad.cv=matrix(0, nrow = run, ncol = length(n.seq))
cp_drct.11.scad.cv=matrix(0, nrow = run, ncol = length(n.seq))
cp_drct.20.scad.cv=matrix(0, nrow = run, ncol = length(n.seq))
cp_drct.21.scad.cv=matrix(0, nrow = run, ncol = length(n.seq))

reject_indrct.10.scad.cv=matrix(0, nrow = run, ncol = length(n.seq))
reject_indrct.11.scad.cv=matrix(0, nrow = run, ncol = length(n.seq))
reject_indrct.20.scad.cv=matrix(0, nrow = run, ncol = length(n.seq))
reject_indrct.21.scad.cv=matrix(0, nrow = run, ncol = length(n.seq))

reject_drct.10.scad.cv=matrix(0, nrow = run, ncol = length(n.seq))
reject_drct.11.scad.cv=matrix(0, nrow = run, ncol = length(n.seq))
reject_drct.20.scad.cv=matrix(0, nrow = run, ncol = length(n.seq))
reject_drct.21.scad.cv=matrix(0, nrow = run, ncol = length(n.seq))

lambda_scad.cv=matrix(0, nrow = run, ncol = length(n.seq))






######## SCAD HBIC
bias_theta.scad.hbic=matrix(0, nrow = run, ncol = length(n.seq))
sr_theta.scad.hbic=matrix(0, nrow = run, ncol = length(n.seq))

est_indrct.10.scad.hbic=matrix(0, nrow = run, ncol = length(n.seq))
est_indrct.11.scad.hbic=matrix(0, nrow = run, ncol = length(n.seq))
est_indrct.20.scad.hbic=matrix(0, nrow = run, ncol = length(n.seq))
est_indrct.21.scad.hbic=matrix(0, nrow = run, ncol = length(n.seq))

est_drct.10.scad.hbic=matrix(0, nrow = run, ncol = length(n.seq))
est_drct.11.scad.hbic=matrix(0, nrow = run, ncol = length(n.seq))
est_drct.20.scad.hbic=matrix(0, nrow = run, ncol = length(n.seq))
est_drct.21.scad.hbic=matrix(0, nrow = run, ncol = length(n.seq))


se_indrct.10.scad.hbic=matrix(0, nrow = run, ncol = length(n.seq))
se_indrct.11.scad.hbic=matrix(0, nrow = run, ncol = length(n.seq))
se_indrct.20.scad.hbic=matrix(0, nrow = run, ncol = length(n.seq))
se_indrct.21.scad.hbic=matrix(0, nrow = run, ncol = length(n.seq))

se_drct.10.scad.hbic=matrix(0, nrow = run, ncol = length(n.seq))
se_drct.11.scad.hbic=matrix(0, nrow = run, ncol = length(n.seq))
se_drct.20.scad.hbic=matrix(0, nrow = run, ncol = length(n.seq))
se_drct.21.scad.hbic=matrix(0, nrow = run, ncol = length(n.seq))

cp_indrct.10.scad.hbic=matrix(0, nrow = run, ncol = length(n.seq))
cp_indrct.11.scad.hbic=matrix(0, nrow = run, ncol = length(n.seq))
cp_indrct.20.scad.hbic=matrix(0, nrow = run, ncol = length(n.seq))
cp_indrct.21.scad.hbic=matrix(0, nrow = run, ncol = length(n.seq))

cp_drct.10.scad.hbic=matrix(0, nrow = run, ncol = length(n.seq))
cp_drct.11.scad.hbic=matrix(0, nrow = run, ncol = length(n.seq))
cp_drct.20.scad.hbic=matrix(0, nrow = run, ncol = length(n.seq))
cp_drct.21.scad.hbic=matrix(0, nrow = run, ncol = length(n.seq))

reject_indrct.10.scad.hbic=matrix(0, nrow = run, ncol = length(n.seq))
reject_indrct.11.scad.hbic=matrix(0, nrow = run, ncol = length(n.seq))
reject_indrct.20.scad.hbic=matrix(0, nrow = run, ncol = length(n.seq))
reject_indrct.21.scad.hbic=matrix(0, nrow = run, ncol = length(n.seq))

reject_drct.10.scad.hbic=matrix(0, nrow = run, ncol = length(n.seq))
reject_drct.11.scad.hbic=matrix(0, nrow = run, ncol = length(n.seq))
reject_drct.20.scad.hbic=matrix(0, nrow = run, ncol = length(n.seq))
reject_drct.21.scad.hbic=matrix(0, nrow = run, ncol = length(n.seq))

lambda_scad.hbic=matrix(0, nrow = run, ncol = length(n.seq))





## Runze cv
est_indrct.1.runze.cv=matrix(0, nrow = run, ncol = length(n.seq))
est_indrct.2.runze.cv=matrix(0, nrow = run, ncol = length(n.seq))
se_indrct.1.runze.cv=matrix(0, nrow = run, ncol = length(n.seq))
se_indrct.2.runze.cv=matrix(0, nrow = run, ncol = length(n.seq))

cp_indrct.10.runze.cv=matrix(0, nrow = run, ncol = length(n.seq))
cp_indrct.11.runze.cv=matrix(0, nrow = run, ncol = length(n.seq))
cp_indrct.20.runze.cv=matrix(0, nrow = run, ncol = length(n.seq))
cp_indrct.21.runze.cv=matrix(0, nrow = run, ncol = length(n.seq))

reject_indrct.10.runze.cv=matrix(0, nrow = run, ncol = length(n.seq))
reject_indrct.11.runze.cv=matrix(0, nrow = run, ncol = length(n.seq))
reject_indrct.20.runze.cv=matrix(0, nrow = run, ncol = length(n.seq))
reject_indrct.21.runze.cv=matrix(0, nrow = run, ncol = length(n.seq))

est_drct.1.runze.cv=matrix(0, nrow = run, ncol = length(n.seq))
est_drct.2.runze.cv=matrix(0, nrow = run, ncol = length(n.seq))
se_drct.1.runze.cv=matrix(0, nrow = run, ncol = length(n.seq))
se_drct.2.runze.cv=matrix(0, nrow = run, ncol = length(n.seq))

cp_drct.10.runze.cv=matrix(0, nrow = run, ncol = length(n.seq))
cp_drct.11.runze.cv=matrix(0, nrow = run, ncol = length(n.seq))
cp_drct.20.runze.cv=matrix(0, nrow = run, ncol = length(n.seq))
cp_drct.21.runze.cv=matrix(0, nrow = run, ncol = length(n.seq))

reject_drct.10.runze.cv=matrix(0, nrow = run, ncol = length(n.seq))
reject_drct.11.runze.cv=matrix(0, nrow = run, ncol = length(n.seq))
reject_drct.20.runze.cv=matrix(0, nrow = run, ncol = length(n.seq))
reject_drct.21.runze.cv=matrix(0, nrow = run, ncol = length(n.seq))


## Runze hbic
est_indrct.1.runze.hbic=matrix(0, nrow = run, ncol = length(n.seq))
est_indrct.2.runze.hbic=matrix(0, nrow = run, ncol = length(n.seq))
se_indrct.1.runze.hbic=matrix(0, nrow = run, ncol = length(n.seq))
se_indrct.2.runze.hbic=matrix(0, nrow = run, ncol = length(n.seq))

cp_indrct.10.runze.hbic=matrix(0, nrow = run, ncol = length(n.seq))
cp_indrct.11.runze.hbic=matrix(0, nrow = run, ncol = length(n.seq))
cp_indrct.20.runze.hbic=matrix(0, nrow = run, ncol = length(n.seq))
cp_indrct.21.runze.hbic=matrix(0, nrow = run, ncol = length(n.seq))

reject_indrct.10.runze.hbic=matrix(0, nrow = run, ncol = length(n.seq))
reject_indrct.11.runze.hbic=matrix(0, nrow = run, ncol = length(n.seq))
reject_indrct.20.runze.hbic=matrix(0, nrow = run, ncol = length(n.seq))
reject_indrct.21.runze.hbic=matrix(0, nrow = run, ncol = length(n.seq))

est_drct.1.runze.hbic=matrix(0, nrow = run, ncol = length(n.seq))
est_drct.2.runze.hbic=matrix(0, nrow = run, ncol = length(n.seq))
se_drct.1.runze.hbic=matrix(0, nrow = run, ncol = length(n.seq))
se_drct.2.runze.hbic=matrix(0, nrow = run, ncol = length(n.seq))

cp_drct.10.runze.hbic=matrix(0, nrow = run, ncol = length(n.seq))
cp_drct.11.runze.hbic=matrix(0, nrow = run, ncol = length(n.seq))
cp_drct.20.runze.hbic=matrix(0, nrow = run, ncol = length(n.seq))
cp_drct.21.runze.hbic=matrix(0, nrow = run, ncol = length(n.seq))

reject_drct.10.runze.hbic=matrix(0, nrow = run, ncol = length(n.seq))
reject_drct.11.runze.hbic=matrix(0, nrow = run, ncol = length(n.seq))
reject_drct.20.runze.hbic=matrix(0, nrow = run, ncol = length(n.seq))
reject_drct.21.runze.hbic=matrix(0, nrow = run, ncol = length(n.seq))


NIE10 = as.numeric(t(beta1)%*%theta.G)
NIE11 = as.numeric(t(beta1)%*%theta.G) + as.numeric(t(beta1)%*%theta.D1)
NIE20 = as.numeric(t(beta2)%*%theta.G)
NIE21 = as.numeric(t(beta2)%*%theta.G) + as.numeric(t(beta2)%*%theta.D2)

NDE10 = theta.T[1] + t(beta0)%*%theta.D1 + t(mu.X)%*%t(A)%*%theta.D1
NDE11 = theta.T[1] + t(beta0)%*%theta.D1 + t(mu.X)%*%t(A)%*%theta.D1 + as.numeric(t(beta1)%*%theta.D1) 
NDE20 = theta.T[2] + t(beta0)%*%theta.D2 + t(mu.X)%*%t(A)%*%theta.D2
NDE21 = theta.T[2] + t(beta0)%*%theta.D2 + t(mu.X)%*%t(A)%*%theta.D2 + as.numeric(t(beta2)%*%theta.D2) 



for (loop.n in 1:length(n.seq)) {
  
  n=n.seq[loop.n]
  
  one.n=as.matrix(rep(1, times=n))
  
  TRT=matrix(0, nrow = n, ncol = q)
  TRT[,1]=c( rep(1, times=(n/3)), rep(0, times=(2*n/3))  )
  TRT[,2]=c( rep(0, times=(n/3)), rep(1, times=(n/3)), rep(0, times=(n/3))  )
  
  colnames(TRT) = c("T1", "T2")
  
  for (loop.run in 1:run) {
    
    # generate mediation and response
    
    #White_G=matrix(rnorm(n*p), nrow = n, ncol = p)
    #White_G=matrix(runif(n*p, min = -sqrt(3), max = sqrt(3)), n,p)
    #White_G=matrix((rpois(n*p, lambda = 1)-1), n, p)
    #White_G=matrix((2*rbinom(n*p,1,0.5)-1), n,p)
    #White_G=matrix( (sqrt(3/5)*rt(n*p, df=5)) , nrow = n, ncol = p)
    White_G=matrix( (sqrt(6/8)*rt(n*p, df=8)) , nrow = n, ncol = p)
    
    eps_G=White_G %*%Sigma0.5_G_gl
    
    X = matrix( rnorm(n*r, mean = 0, sd=1), n, r)
    X[,1] = rbinom(n, 1, 0.5 ) # the first column of X is from sex
    colnames(X) = paste0("X", 1:ncol(X))
    
    G=one.n%*%t(beta0)+TRT%*%t(B) + X%*%t(A) +eps_G # G is n by p matrix whose ith row is the ith obs of the mediation
    colnames(G) = paste0("G", 1:ncol(G))
    
    eps_Y=as.matrix(rnorm(n, mean = 0, sd=sigma_Y)) # normally distributed
    #eps_Y= sigma_Y * as.matrix(runif(n, min = -sqrt(3), max = sqrt(3))) # uniformly distributed
    #eps_Y = sigma_Y * as.matrix(rpois(n, lambda = 1)-1) # poisson centered
    #eps_Y=sigma_Y * as.matrix(2*rbinom(n,1,0.5)-1) # symmetric bernoulli
    #eps_Y=sigma_Y * as.matrix( sqrt(3/5)*rt(n, df=5) )
    #eps_Y=sigma_Y * as.matrix( sqrt(6/8)*rt(n, df=8) )
    
    
    G.T=matrix(0, nrow = n, ncol = q*p)
    G.T[1:(n/3), 1:p] = G[1:(n/3), ]
    G.T[((n/3)+1): (2*n/3), (p+1):(2*p)]=G[((n/3)+1): (2*n/3), ]
    colnames(G.T) = c( paste0("T1.G", 1:ncol(G)), paste0("T2.G", 1:ncol(G)) )
    
    Y=theta0*one.n + TRT%*%theta.T + X%*%theta.X + G%*%theta.G + G.T%*%as.matrix(vtheta.D) + eps_Y  # B is n by 1 outcome
    colnames(Y) = "Y"
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    ############################################### proposed ##########################################################################
    W_outcome=cbind(TRT, X, G, G.T)
    
    # max( abs( Y - cbind(one.n, W_outcome)%*%theta - eps_Y ) )
    
    fitp.scad=cv.ncvreg(W_outcome,Y, penalty = "SCAD", penalty.factor = p.fac ) # use default lambda
    
    ############################# the following annotations are about the output of cv.ncvreg() and ncvreg() functions #############################
    # lambda.seq = exp(seq(8,-20, length.out=100))
    
    # fitp.scad=cv.ncvreg(W_outcome,Y, penalty = "SCAD", penalty.factor = p.fac, lambda = lambda.seq )
    
    # fitp.scad$lambda - lambda.seq  # the outputs are all zeros, so the lambdas used in cv.ncvreg is exactly the argument "lambda = ".
    
    # fitp.scad$lambda.min - fitp.scad$lambda[fitp.scad$min] # fitp.scad$min is the number location (index) of fitp.scad$lambda which achives the smallest cve.
    
    # fitp.scad$cve # the red points (mean of cross-validation error under several (10) folds) in plot(fitp.scad)
    
    ## ncvreg(...) is the same as cv.ncvreg()$fit. To see this:
    # fit.whole = ncvreg(W_outcome,Y, penalty = "SCAD", penalty.factor = p.fac, lambda = lambda.seq ) # fit the whole data directly without cross-validation
    
    ## lambda are the same:
    # fitp.scad$fit$lambda - fit.whole$lambda # the results are all zero (they are using the same set of lambda)
    
    ## the coefficients are the same:
    # dim(coef(fit.whole)) # 3035  100  which is (1+q+r+p+p*q) (including the intercept) and length(lambda.seq)
    # coef(fitp.scad) # it produces the estimated theta's corresponding to lambda.min
    # dim(fitp.scad$fit$beta) # 3035  100 which is (1+q+r+p+p*q) (including the intercept) and length(lambda.seq)
    
    # dim( fitp.scad$fit$beta - coef(fit.whole) ) # 3035  100 which is (1+q+r+p+p*q) and length(lambda.seq)
    # min(abs( fitp.scad$fit$beta - coef(fit.whole) )) # the output is 0
    ## So we do not need to fit "fit.whole = ncvreg(...)" out of the whole data because fitp.scad = cv.ncvreg(...)$fit$beta has already produce the same coefficients from the whole data set
    
    
    ## To verify the statement above, you can also try the following codes to see the coefficients corresponding the lambda.min
    # max( abs( coef(fit.whole)[, fitp.scad$min] - coef(fitp.scad) ) ) # the output is 0
    # max( abs( coef(fit.whole)[, fitp.scad$min] - coef(fitp.scad, s = "lambda.min") ) ) # the output is 0
    # max( abs( coef(fit.whole)[, fitp.scad$min] - fitp.scad$fit$beta[, fitp.scad$min] ) ) # the output is 0, which means that coef(fitp.scad, s = "lambda.min") is calculated from the whole data set.
    
    ## One Caveat:
    # summary(fitp.scad)$sigma # This is not the sum of squared residuals/n, this is sqrt(fitp.scad$cve). See the next row:
    # summary(fitp.scad)$sigma - sqrt(fitp.scad$cve) # the outputs are all zeros
    
    ## So there is no direct result for the ( sum of squared residuals / n ) from each lambda
    ########################################################################################################
    
    
    
    
    ###### estimation by cv
    
    # step 1 outcome model
    
    lambda_scad.cv[loop.run, loop.n] = fitp.scad$lambda.min
    
    theta0.hat.cv=coef(fitp.scad, s = "lambda.min")[1]
    
    theta.T.hat.cv=coef(fitp.scad, s = "lambda.min")[(1+1): (1+q)]
    
    theta.X.hat.cv=coef(fitp.scad, s = "lambda.min")[(1+q+1) : (1+q+r)]
    
    theta.G.hat.cv=coef(fitp.scad, s = "lambda.min")[(1+q+r+1):(1+q+r+p)]
    
    theta.D1.hat.cv=coef(fitp.scad, s = "lambda.min")[(1+q+r+p+1):(1+q+r+p+p)]
    
    theta.D2.hat.cv=coef(fitp.scad, s = "lambda.min")[(1+q+r+p+p+1):(1+q+r+p+p+p)]
    
    vtheta.D_hat.cv=c(theta.D1.hat.cv, theta.D2.hat.cv)
    
    names(theta0.hat.cv)=NULL; names(theta.T.hat.cv)=NULL; names(theta.X.hat.cv)=NULL 
    names(theta.G.hat.cv)=NULL; names(vtheta.D_hat.cv)=NULL; names(theta.D1.hat.cv)=NULL; names(theta.D2.hat.cv)=NULL;
    
    bias_theta.scad.cv[loop.run,loop.n]=sum(abs(coef(fitp.scad, s = "lambda.min") - theta))
    sr_theta.scad.cv[loop.run,loop.n]=sum(sign(coef(fitp.scad, s = "lambda.min"))!=sign(theta))
    
    # check the error by different parts
    # sum(sign(coef(fitp.scad, s = "lambda.min"))!=sign(theta)) # support recovery
    # sum(abs(theta.T - theta.T.hat.cv))
    # sum(abs(theta.G.hat.cv - theta.G))
    # sum(abs(theta.X.hat.cv - theta.X))
    # sum(abs(theta.D1.hat.cv - theta.D1))
    # sum(abs(theta.D2.hat.cv - theta.D2))
    
    # step 2 mediation model
    
    W_mediator=cbind(one.n, TRT, X)
    
    # betahat = solve(t(W_mediator)%*%W_mediator) %*% t(W_mediator) %*% G
    
    # sum(abs(betahat - t(cbind(beta0,beta1, beta2, A) ) ))
    
    # beta.theta.G.hat.cv = betahat %*% theta.G.hat.cv
    
    # beta.theta.D1.hat.cv = betahat %*% theta.D1.hat.cv
    
    # beta.theta.D2.hat.cv = betahat %*% theta.D2.hat.cv
    
    # calculate the estimation error
    
    SE_smple.cv = OraProjInt_SE_sample_X(Y = Y, Trt = TRT, Mediator = G, Covariate = X, 
                                         Trt.Medator = G.T,
                                         theta0_hat = theta0.hat.cv,
                                         theta.T_hat = theta.T.hat.cv,
                                         theta.X_hat = theta.X.hat.cv,
                                         theta.G_hat = theta.G.hat.cv,
                                         vtheta.D_hat = vtheta.D_hat.cv)
    
    
    est_indrct.10.scad.cv[loop.run,loop.n]= SE_smple.cv$NIE_inference["NIE(T 1 = 0)", "Estimate"] 
    est_indrct.20.scad.cv[loop.run,loop.n]= SE_smple.cv$NIE_inference["NIE(T 2 = 0)", "Estimate"] 
    
    est_indrct.11.scad.cv[loop.run,loop.n]= SE_smple.cv$NIE_inference["NIE(T 1 = 1)", "Estimate"] 
    est_indrct.21.scad.cv[loop.run,loop.n]= SE_smple.cv$NIE_inference["NIE(T 2 = 1)", "Estimate"] 
    
    est_drct.10.scad.cv[loop.run,loop.n]= SE_smple.cv$NDE_inference["NDE(T 1 = 0)", "Estimate"] 
    est_drct.20.scad.cv[loop.run,loop.n]= SE_smple.cv$NDE_inference["NDE(T 2 = 0)", "Estimate"] 
    
    est_drct.11.scad.cv[loop.run,loop.n]= SE_smple.cv$NDE_inference["NDE(T 1 = 1)", "Estimate"] 
    est_drct.21.scad.cv[loop.run,loop.n]= SE_smple.cv$NDE_inference["NDE(T 2 = 1)", "Estimate"] 
    
    
    
    se_indrct.10.scad.cv[loop.run,loop.n]= SE_smple.cv$NIE_inference["NIE(T 1 = 0)", "SE"] 
    se_indrct.20.scad.cv[loop.run,loop.n]= SE_smple.cv$NIE_inference["NIE(T 2 = 0)", "SE"] 
    
    se_indrct.11.scad.cv[loop.run,loop.n]= SE_smple.cv$NIE_inference["NIE(T 1 = 1)", "SE"] 
    se_indrct.21.scad.cv[loop.run,loop.n]= SE_smple.cv$NIE_inference["NIE(T 2 = 1)", "SE"] 
    
    se_drct.10.scad.cv[loop.run,loop.n]= SE_smple.cv$NDE_inference["NDE(T 1 = 0)", "SE"] 
    se_drct.20.scad.cv[loop.run,loop.n]= SE_smple.cv$NDE_inference["NDE(T 2 = 0)", "SE"] 
    
    se_drct.11.scad.cv[loop.run,loop.n]= SE_smple.cv$NDE_inference["NDE(T 1 = 1)", "SE"] 
    se_drct.21.scad.cv[loop.run,loop.n]= SE_smple.cv$NDE_inference["NDE(T 2 = 1)", "SE"] 
    
    # cp
    if ( ( ( NIE10 ) <= ( est_indrct.10.scad.cv[loop.run,loop.n] + qnorm(0.975)*se_indrct.10.scad.cv[loop.run, loop.n] ) ) & ( ( NIE10 ) >= ( est_indrct.10.scad.cv[loop.run,loop.n] - qnorm(0.975)*se_indrct.10.scad.cv[loop.run, loop.n] ) )  ){
      cp_indrct.10.scad.cv[loop.run, loop.n] = 1
    }
    
    if ( ( ( NIE20) <= ( est_indrct.20.scad.cv[loop.run,loop.n] + qnorm(0.975)*se_indrct.20.scad.cv[loop.run, loop.n] ) ) & ( ( NIE20 ) >= ( est_indrct.20.scad.cv[loop.run,loop.n] - qnorm(0.975)*se_indrct.20.scad.cv[loop.run, loop.n] ) )  ){
      cp_indrct.20.scad.cv[loop.run, loop.n] = 1
    }
    
    
    if ( ( ( NIE11 ) <= ( est_indrct.11.scad.cv[loop.run,loop.n] + qnorm(0.975)*se_indrct.11.scad.cv[loop.run, loop.n] ) ) & ( ( NIE11 ) >= ( est_indrct.11.scad.cv[loop.run,loop.n] - qnorm(0.975)*se_indrct.11.scad.cv[loop.run, loop.n] ) )  ){
      cp_indrct.11.scad.cv[loop.run, loop.n] = 1
    }
    
    if ( ( ( NIE21 ) <= ( est_indrct.21.scad.cv[loop.run,loop.n] + qnorm(0.975)*se_indrct.21.scad.cv[loop.run, loop.n] ) ) & ( ( NIE21 ) >= ( est_indrct.21.scad.cv[loop.run,loop.n] - qnorm(0.975)*se_indrct.21.scad.cv[loop.run, loop.n] ) )  ){
      cp_indrct.21.scad.cv[loop.run, loop.n] = 1
    }
    
    #
    if ( ( ( NDE10 ) <= ( est_drct.10.scad.cv[loop.run,loop.n] + qnorm(0.975)*se_drct.10.scad.cv[loop.run, loop.n] ) ) & ( ( NDE10 ) >= ( est_drct.10.scad.cv[loop.run,loop.n] - qnorm(0.975)*se_drct.10.scad.cv[loop.run, loop.n] ) )  ){
      cp_drct.10.scad.cv[loop.run, loop.n] = 1
    }
    
    if ( ( ( NDE20 ) <= ( est_drct.20.scad.cv[loop.run,loop.n] + qnorm(0.975)*se_drct.20.scad.cv[loop.run, loop.n] ) ) & ( ( NDE20 ) >= ( est_drct.20.scad.cv[loop.run,loop.n] - qnorm(0.975)*se_drct.20.scad.cv[loop.run, loop.n] ) )  ){
      cp_drct.20.scad.cv[loop.run, loop.n] = 1
    }
    
    if ( ( ( NDE11 ) <= ( est_drct.11.scad.cv[loop.run,loop.n] + qnorm(0.975)*se_drct.11.scad.cv[loop.run, loop.n] ) ) & ( ( NDE11 ) >= ( est_drct.11.scad.cv[loop.run,loop.n] - qnorm(0.975)*se_drct.11.scad.cv[loop.run, loop.n] ) )  ){
      cp_drct.11.scad.cv[loop.run, loop.n] = 1
    }
    
    if ( ( ( NDE21 ) <= ( est_drct.21.scad.cv[loop.run,loop.n] + qnorm(0.975)*se_drct.21.scad.cv[loop.run, loop.n] ) ) & ( ( NDE21 ) >= ( est_drct.21.scad.cv[loop.run,loop.n] - qnorm(0.975)*se_drct.21.scad.cv[loop.run, loop.n] ) )  ){
      cp_drct.21.scad.cv[loop.run, loop.n] = 1
    }
    
    # power
    if ( ( 0 > ( est_indrct.10.scad.cv[loop.run,loop.n] + qnorm(0.975)*se_indrct.10.scad.cv[loop.run, loop.n] ) ) | ( 0 < ( est_indrct.10.scad.cv[loop.run,loop.n] - qnorm(0.975)*se_indrct.10.scad.cv[loop.run, loop.n] ) )  ){
      reject_indrct.10.scad.cv[loop.run, loop.n] = 1
    }
    
    if ( ( 0 > ( est_indrct.20.scad.cv[loop.run,loop.n] + qnorm(0.975)*se_indrct.20.scad.cv[loop.run, loop.n] ) ) | ( 0 < ( est_indrct.20.scad.cv[loop.run,loop.n] - qnorm(0.975)*se_indrct.20.scad.cv[loop.run, loop.n] ) )  ){
      reject_indrct.20.scad.cv[loop.run, loop.n] = 1
    }
    
    
    if ( ( 0 > ( est_indrct.11.scad.cv[loop.run,loop.n] + qnorm(0.975)*se_indrct.11.scad.cv[loop.run, loop.n] ) ) | ( 0 < ( est_indrct.11.scad.cv[loop.run,loop.n] - qnorm(0.975)*se_indrct.11.scad.cv[loop.run, loop.n] ) )  ){
      reject_indrct.11.scad.cv[loop.run, loop.n] = 1
    }
    
    if ( ( 0 > ( est_indrct.21.scad.cv[loop.run,loop.n] + qnorm(0.975)*se_indrct.21.scad.cv[loop.run, loop.n] ) ) | ( 0 < ( est_indrct.21.scad.cv[loop.run,loop.n] - qnorm(0.975)*se_indrct.21.scad.cv[loop.run, loop.n] ) )  ){
      reject_indrct.21.scad.cv[loop.run, loop.n] = 1
    }
    
    #
    if ( ( 0 > ( est_drct.10.scad.cv[loop.run,loop.n] + qnorm(0.975)*se_drct.10.scad.cv[loop.run, loop.n] ) ) | ( 0 < ( est_drct.10.scad.cv[loop.run,loop.n] - qnorm(0.975)*se_drct.10.scad.cv[loop.run, loop.n] ) )  ){
      reject_drct.10.scad.cv[loop.run, loop.n] = 1
    }
    
    if ( ( 0 > ( est_drct.20.scad.cv[loop.run,loop.n] + qnorm(0.975)*se_drct.20.scad.cv[loop.run, loop.n] ) ) | ( 0 < ( est_drct.20.scad.cv[loop.run,loop.n] - qnorm(0.975)*se_drct.20.scad.cv[loop.run, loop.n] ) )  ){
      reject_drct.20.scad.cv[loop.run, loop.n] = 1
    }
    
    if ( ( 0 > ( est_drct.11.scad.cv[loop.run,loop.n] + qnorm(0.975)*se_drct.11.scad.cv[loop.run, loop.n] ) ) | ( 0 < ( est_drct.11.scad.cv[loop.run,loop.n] - qnorm(0.975)*se_drct.11.scad.cv[loop.run, loop.n] ) )  ){
      reject_drct.11.scad.cv[loop.run, loop.n] = 1
    }
    
    if ( ( 0 > ( est_drct.21.scad.cv[loop.run,loop.n] + qnorm(0.975)*se_drct.21.scad.cv[loop.run, loop.n] ) ) | ( 0 < ( est_drct.21.scad.cv[loop.run,loop.n] - qnorm(0.975)*se_drct.21.scad.cv[loop.run, loop.n] ) )  ){
      reject_drct.21.scad.cv[loop.run, loop.n] = 1
    }
    
    
    ### HBIC ###
    
    # First calculate hbic corresponding to each lambda
    lambda.seq = fitp.scad$lambda
    
    hbic = c()
    for (loop.lambda in 1:length(lambda.seq)) {
      sigma2_Yhat.tmp =  sum( ( Y - cbind(one.n, W_outcome) %*% fitp.scad$fit$beta[,loop.lambda] )^2 ) / n 
      theta.G.tmp = (fitp.scad$fit$beta[,loop.lambda])[(1+q+r+1):(1+q+r+p)]
      vtheta.D.tmp = (fitp.scad$fit$beta[,loop.lambda])[(1+q+r+p+1):(1+q+r+p+p*q)]
      df = 1+q+r+length(which(theta.G.tmp!= 0))+length(which(vtheta.D.tmp!= 0))
      hbic[loop.lambda] = log( sigma2_Yhat.tmp ) +  df*log(log(n))*log(1+q+r+p+p*q)/n
    }
    
    ## aa is the log sum of squared residuals / n, which is decreasing linearly after it hits the lambda.min from cross-validation
    # aa = c()
    # for (loop.lambda in 1:length(lambda.seq)) {
    #   sigma2_Yhat.tmp =  sum( ( Y - cbind(one.n, W_outcome) %*% fitp.scad$fit$beta[,loop.lambda] )^2 ) / n
    #   aa[loop.lambda] = log( sigma2_Yhat.tmp )
    # }
    ## bb is the |M_lambda|_0 * Cn * log(1+q+r+p+p*q) / n in eqn(2.8) of the HBIC paper "CALIBRATING NONCONVEX PENALIZED REGRESSION IN ULTRA-HIGH DIMENSION"
    ## bb is increasing slower and slower
    # bb = c()
    # for (loop.lambda in 1:length(lambda.seq)) {
    #  theta.G.tmp = (fitp.scad$fit$beta[,loop.lambda])[(1+q+r+1):(1+q+r+p)]
    #  vtheta.D.tmp = (fitp.scad$fit$beta[,loop.lambda])[(1+q+r+p+1):(1+q+r+p+p*q)]
    #  df = 1+q+r+length(which(theta.G.tmp!= 0))+length(which(vtheta.D.tmp!= 0))
    #  bb[loop.lambda] = df*log(log(n))*log(1+q+r+p+p*q)/n
    # }
    ## since aa+bb=hbic, as lambda decreases (the entries in lambda.seq is decreasing), aa will dominate, so it tends to select a smaller lambda, if the range of lambda.seq is too large. So we use the lambda sequence generated by cross-validation.
    
    hbic.min = which.min(hbic)
    
    lambda_scad.hbic[loop.run, loop.n] = lambda.seq[hbic.min]
    
    theta0.hat.hbic=fitp.scad$fit$beta[,hbic.min][1]
    
    theta.T.hat.hbic=fitp.scad$fit$beta[,hbic.min][(1+1): (1+q)]
    
    theta.X.hat.hbic=fitp.scad$fit$beta[,hbic.min][(1+q+1) : (1+q+r)]
    
    theta.G.hat.hbic=fitp.scad$fit$beta[,hbic.min][(1+q+r+1):(1+q+r+p)]
    
    theta.D1.hat.hbic=fitp.scad$fit$beta[,hbic.min][(1+q+r+p+1):(1+q+r+p+p)]
    
    theta.D2.hat.hbic=fitp.scad$fit$beta[,hbic.min][(1+q+r+p+p+1):(1+q+r+p+p+p)]
    
    vtheta.D_hat.hbic=c(theta.D1.hat.hbic, theta.D2.hat.hbic)
    
    names(theta0.hat.hbic)=NULL; names(theta.T.hat.hbic)=NULL; names(theta.X.hat.hbic)=NULL 
    names(theta.G.hat.hbic)=NULL; names(vtheta.D_hat.hbic)=NULL; names(theta.D1.hat.hbic)=NULL; names(theta.D2.hat.hbic)=NULL;
    
    
    bias_theta.scad.hbic[loop.run,loop.n]=sum(abs(fitp.scad$fit$beta[,hbic.min] - theta))
    sr_theta.scad.hbic[loop.run,loop.n]=sum(sign(fitp.scad$fit$beta[,hbic.min])!=sign(theta))
    
    # step 2 mediation model
    
    W_mediator=cbind(one.n, TRT, X)
    
    # betahat = solve(t(W_mediator)%*%W_mediator) %*% t(W_mediator) %*% G
    
    # beta.theta.G.hat.hbic = betahat %*% theta.G.hat.hbic
    
    # beta.theta.D1.hat.hbic = betahat %*% theta.D1.hat.hbic
    
    # beta.theta.D2.hat.hbic = betahat %*% theta.D2.hat.hbic
    
    # calculate the estimation error
    
    SE_smple.hbic = OraProjInt_SE_sample_X(Y = Y, Trt = TRT, Mediator = G, Covariate = X, 
                                           Trt.Medator = G.T,
                                           theta0_hat = theta0.hat.hbic,
                                           theta.T_hat = theta.T.hat.hbic,
                                           theta.X_hat = theta.X.hat.hbic,
                                           theta.G_hat = theta.G.hat.hbic,
                                           vtheta.D_hat = vtheta.D_hat.hbic)
    
    
    est_indrct.10.scad.hbic[loop.run,loop.n]= SE_smple.hbic$NIE_inference["NIE(T 1 = 0)", "Estimate"] 
    est_indrct.20.scad.hbic[loop.run,loop.n]= SE_smple.hbic$NIE_inference["NIE(T 2 = 0)", "Estimate"] 
    
    est_indrct.11.scad.hbic[loop.run,loop.n]= SE_smple.hbic$NIE_inference["NIE(T 1 = 1)", "Estimate"] 
    est_indrct.21.scad.hbic[loop.run,loop.n]= SE_smple.hbic$NIE_inference["NIE(T 2 = 1)", "Estimate"] 
    
    est_drct.10.scad.hbic[loop.run,loop.n]= SE_smple.hbic$NDE_inference["NDE(T 1 = 0)", "Estimate"] 
    est_drct.20.scad.hbic[loop.run,loop.n]= SE_smple.hbic$NDE_inference["NDE(T 2 = 0)", "Estimate"] 
    
    est_drct.11.scad.hbic[loop.run,loop.n]= SE_smple.hbic$NDE_inference["NDE(T 1 = 1)", "Estimate"] 
    est_drct.21.scad.hbic[loop.run,loop.n]= SE_smple.hbic$NDE_inference["NDE(T 2 = 1)", "Estimate"] 
    
    
    
    se_indrct.10.scad.hbic[loop.run,loop.n]= SE_smple.hbic$NIE_inference["NIE(T 1 = 0)", "SE"] 
    se_indrct.20.scad.hbic[loop.run,loop.n]= SE_smple.hbic$NIE_inference["NIE(T 2 = 0)", "SE"] 
    
    se_indrct.11.scad.hbic[loop.run,loop.n]= SE_smple.hbic$NIE_inference["NIE(T 1 = 1)", "SE"] 
    se_indrct.21.scad.hbic[loop.run,loop.n]= SE_smple.hbic$NIE_inference["NIE(T 2 = 1)", "SE"] 
    
    se_drct.10.scad.hbic[loop.run,loop.n]= SE_smple.hbic$NDE_inference["NDE(T 1 = 0)", "SE"] 
    se_drct.20.scad.hbic[loop.run,loop.n]= SE_smple.hbic$NDE_inference["NDE(T 2 = 0)", "SE"] 
    
    se_drct.11.scad.hbic[loop.run,loop.n]= SE_smple.hbic$NDE_inference["NDE(T 1 = 1)", "SE"] 
    se_drct.21.scad.hbic[loop.run,loop.n]= SE_smple.hbic$NDE_inference["NDE(T 2 = 1)", "SE"] 
    
    # cp
    if ( ( ( NIE10 ) <= ( est_indrct.10.scad.hbic[loop.run,loop.n] + qnorm(0.975)*se_indrct.10.scad.hbic[loop.run, loop.n] ) ) & ( ( NIE10 ) >= ( est_indrct.10.scad.hbic[loop.run,loop.n] - qnorm(0.975)*se_indrct.10.scad.hbic[loop.run, loop.n] ) )  ){
      cp_indrct.10.scad.hbic[loop.run, loop.n] = 1
    }
    
    if ( ( ( NIE20) <= ( est_indrct.20.scad.hbic[loop.run,loop.n] + qnorm(0.975)*se_indrct.20.scad.hbic[loop.run, loop.n] ) ) & ( ( NIE20 ) >= ( est_indrct.20.scad.hbic[loop.run,loop.n] - qnorm(0.975)*se_indrct.20.scad.hbic[loop.run, loop.n] ) )  ){
      cp_indrct.20.scad.hbic[loop.run, loop.n] = 1
    }
    
    
    if ( ( ( NIE11 ) <= ( est_indrct.11.scad.hbic[loop.run,loop.n] + qnorm(0.975)*se_indrct.11.scad.hbic[loop.run, loop.n] ) ) & ( ( NIE11 ) >= ( est_indrct.11.scad.hbic[loop.run,loop.n] - qnorm(0.975)*se_indrct.11.scad.hbic[loop.run, loop.n] ) )  ){
      cp_indrct.11.scad.hbic[loop.run, loop.n] = 1
    }
    
    if ( ( ( NIE21 ) <= ( est_indrct.21.scad.hbic[loop.run,loop.n] + qnorm(0.975)*se_indrct.21.scad.hbic[loop.run, loop.n] ) ) & ( ( NIE21 ) >= ( est_indrct.21.scad.hbic[loop.run,loop.n] - qnorm(0.975)*se_indrct.21.scad.hbic[loop.run, loop.n] ) )  ){
      cp_indrct.21.scad.hbic[loop.run, loop.n] = 1
    }
    
    #
    if ( ( ( NDE10 ) <= ( est_drct.10.scad.hbic[loop.run,loop.n] + qnorm(0.975)*se_drct.10.scad.hbic[loop.run, loop.n] ) ) & ( ( NDE10 ) >= ( est_drct.10.scad.hbic[loop.run,loop.n] - qnorm(0.975)*se_drct.10.scad.hbic[loop.run, loop.n] ) )  ){
      cp_drct.10.scad.hbic[loop.run, loop.n] = 1
    }
    
    if ( ( ( NDE20 ) <= ( est_drct.20.scad.hbic[loop.run,loop.n] + qnorm(0.975)*se_drct.20.scad.hbic[loop.run, loop.n] ) ) & ( ( NDE20 ) >= ( est_drct.20.scad.hbic[loop.run,loop.n] - qnorm(0.975)*se_drct.20.scad.hbic[loop.run, loop.n] ) )  ){
      cp_drct.20.scad.hbic[loop.run, loop.n] = 1
    }
    
    if ( ( ( NDE11 ) <= ( est_drct.11.scad.hbic[loop.run,loop.n] + qnorm(0.975)*se_drct.11.scad.hbic[loop.run, loop.n] ) ) & ( ( NDE11 ) >= ( est_drct.11.scad.hbic[loop.run,loop.n] - qnorm(0.975)*se_drct.11.scad.hbic[loop.run, loop.n] ) )  ){
      cp_drct.11.scad.hbic[loop.run, loop.n] = 1
    }
    
    if ( ( ( NDE21 ) <= ( est_drct.21.scad.hbic[loop.run,loop.n] + qnorm(0.975)*se_drct.21.scad.hbic[loop.run, loop.n] ) ) & ( ( NDE21 ) >= ( est_drct.21.scad.hbic[loop.run,loop.n] - qnorm(0.975)*se_drct.21.scad.hbic[loop.run, loop.n] ) )  ){
      cp_drct.21.scad.hbic[loop.run, loop.n] = 1
    }
    
    # power
    if ( ( 0 > ( est_indrct.10.scad.hbic[loop.run,loop.n] + qnorm(0.975)*se_indrct.10.scad.hbic[loop.run, loop.n] ) ) | ( 0 < ( est_indrct.10.scad.hbic[loop.run,loop.n] - qnorm(0.975)*se_indrct.10.scad.hbic[loop.run, loop.n] ) )  ){
      reject_indrct.10.scad.hbic[loop.run, loop.n] = 1
    }
    
    if ( ( 0 > ( est_indrct.20.scad.hbic[loop.run,loop.n] + qnorm(0.975)*se_indrct.20.scad.hbic[loop.run, loop.n] ) ) | ( 0 < ( est_indrct.20.scad.hbic[loop.run,loop.n] - qnorm(0.975)*se_indrct.20.scad.hbic[loop.run, loop.n] ) )  ){
      reject_indrct.20.scad.hbic[loop.run, loop.n] = 1
    }
    
    
    if ( ( 0 > ( est_indrct.11.scad.hbic[loop.run,loop.n] + qnorm(0.975)*se_indrct.11.scad.hbic[loop.run, loop.n] ) ) | ( 0 < ( est_indrct.11.scad.hbic[loop.run,loop.n] - qnorm(0.975)*se_indrct.11.scad.hbic[loop.run, loop.n] ) )  ){
      reject_indrct.11.scad.hbic[loop.run, loop.n] = 1
    }
    
    if ( ( 0 > ( est_indrct.21.scad.hbic[loop.run,loop.n] + qnorm(0.975)*se_indrct.21.scad.hbic[loop.run, loop.n] ) ) | ( 0 < ( est_indrct.21.scad.hbic[loop.run,loop.n] - qnorm(0.975)*se_indrct.21.scad.hbic[loop.run, loop.n] ) )  ){
      reject_indrct.21.scad.hbic[loop.run, loop.n] = 1
    }
    
    #
    if ( ( 0 > ( est_drct.10.scad.hbic[loop.run,loop.n] + qnorm(0.975)*se_drct.10.scad.hbic[loop.run, loop.n] ) ) | ( 0 < ( est_drct.10.scad.hbic[loop.run,loop.n] - qnorm(0.975)*se_drct.10.scad.hbic[loop.run, loop.n] ) )  ){
      reject_drct.10.scad.hbic[loop.run, loop.n] = 1
    }
    
    if ( ( 0 > ( est_drct.20.scad.hbic[loop.run,loop.n] + qnorm(0.975)*se_drct.20.scad.hbic[loop.run, loop.n] ) ) | ( 0 < ( est_drct.20.scad.hbic[loop.run,loop.n] - qnorm(0.975)*se_drct.20.scad.hbic[loop.run, loop.n] ) )  ){
      reject_drct.20.scad.hbic[loop.run, loop.n] = 1
    }
    
    if ( ( 0 > ( est_drct.11.scad.hbic[loop.run,loop.n] + qnorm(0.975)*se_drct.11.scad.hbic[loop.run, loop.n] ) ) | ( 0 < ( est_drct.11.scad.hbic[loop.run,loop.n] - qnorm(0.975)*se_drct.11.scad.hbic[loop.run, loop.n] ) )  ){
      reject_drct.11.scad.hbic[loop.run, loop.n] = 1
    }
    
    if ( ( 0 > ( est_drct.21.scad.hbic[loop.run,loop.n] + qnorm(0.975)*se_drct.21.scad.hbic[loop.run, loop.n] ) ) | ( 0 < ( est_drct.21.scad.hbic[loop.run,loop.n] - qnorm(0.975)*se_drct.21.scad.hbic[loop.run, loop.n] ) )  ){
      reject_drct.21.scad.hbic[loop.run, loop.n] = 1
    }
    #########################################################################################################################
    
    
    
    
    
    
    
    
    
    
    
    
    ############### runze ########################################################################
    W_outcome_runze=cbind(TRT, X, G)
    
    # max( abs( Y - cbind(one.n, W_outcome)%*%theta - eps_Y ) )
    
    fitp.runze=cv.ncvreg(W_outcome_runze, Y, penalty = "SCAD", penalty.factor = p.fac.runze ) # use default lambda
    # plot(fitp.scad)
    
    
    ###### estimation by cv
    
    # step 1 outcome model
    
    theta0.hat.runze.cv=coef(fitp.runze, s = "lambda.min")[1]
    
    theta.T.hat.runze.cv=coef(fitp.runze, s = "lambda.min")[(1+1): (1+q)]
    
    theta.X.hat.runze.cv=coef(fitp.runze, s = "lambda.min")[(1+q+1) : (1+q+r)]
    
    theta.G.hat.runze.cv=coef(fitp.runze, s = "lambda.min")[(1+q+r+1):(1+q+r+p)]
    
    names(theta0.hat.runze.cv)=NULL; names(theta.T.hat.runze.cv)=NULL; names(theta.X.hat.runze.cv)=NULL 
    names(theta.G.hat.runze.cv)=NULL
    
    # step 2 mediation model
    
    SE_smple.runze.cv = Runze_SE_sample_noX(Y = Y, Trt = cbind(TRT,X), Mediator = G, 
                                            theta0_hat = theta0.hat.runze.cv,
                                            theta.T_hat = c(theta.T.hat.runze.cv, theta.X.hat.runze.cv),
                                            theta.G_hat = theta.G.hat.runze.cv)
    
    
    est_indrct.1.runze.cv[loop.run,loop.n] = SE_smple.runze.cv$NIE_est["T1", 1] 
    est_indrct.2.runze.cv[loop.run,loop.n] = SE_smple.runze.cv$NIE_est["T2", 1]
    
    se_indrct.1.runze.cv[loop.run,loop.n]= SE_smple.runze.cv$NIE_SE["T1", 1] 
    se_indrct.2.runze.cv[loop.run,loop.n]= SE_smple.runze.cv$NIE_SE["T2", 1] 
    
    est_drct.1.runze.cv[loop.run, loop.n] = SE_smple.runze.cv$NDE_est["T1", 1]
    est_drct.2.runze.cv[loop.run, loop.n] = SE_smple.runze.cv$NDE_est["T2", 1]
    
    se_drct.1.runze.cv[loop.run, loop.n] = SE_smple.runze.cv$NDE_SE["T1", 1]
    se_drct.1.runze.cv[loop.run, loop.n] = SE_smple.runze.cv$NDE_SE["T2", 1]
    
    # cp
    if ( ( ( NIE10 ) <= ( est_indrct.1.runze.cv[loop.run,loop.n] + qnorm(0.975)*se_indrct.1.runze.cv[loop.run, loop.n] ) ) & ( ( NIE10 ) >= ( est_indrct.1.runze.cv[loop.run,loop.n] - qnorm(0.975)*se_indrct.1.runze.cv[loop.run, loop.n] ) )  ){
      cp_indrct.10.runze.cv[loop.run, loop.n] = 1
    }
    
    if ( ( ( NIE20) <= ( est_indrct.2.runze.cv[loop.run,loop.n] + qnorm(0.975)*se_indrct.2.runze.cv[loop.run, loop.n] ) ) & ( ( NIE20 ) >= ( est_indrct.2.runze.cv[loop.run,loop.n] - qnorm(0.975)*se_indrct.2.runze.cv[loop.run, loop.n] ) )  ){
      cp_indrct.20.runze.cv[loop.run, loop.n] = 1
    }
    
    
    if ( ( ( NIE11 ) <= ( est_indrct.1.runze.cv[loop.run,loop.n] + qnorm(0.975)*se_indrct.1.runze.cv[loop.run, loop.n] ) ) & ( ( NIE11 ) >= ( est_indrct.1.runze.cv[loop.run,loop.n] - qnorm(0.975)*se_indrct.1.runze.cv[loop.run, loop.n] ) )  ){
      cp_indrct.11.runze.cv[loop.run, loop.n] = 1
    }
    
    if ( ( ( NIE21 ) <= ( est_indrct.2.runze.cv[loop.run,loop.n] + qnorm(0.975)*se_indrct.2.runze.cv[loop.run, loop.n] ) ) & ( ( NIE21 ) >= ( est_indrct.2.runze.cv[loop.run,loop.n] - qnorm(0.975)*se_indrct.2.runze.cv[loop.run, loop.n] ) )  ){
      cp_indrct.21.runze.cv[loop.run, loop.n] = 1
    }
    
    #
    if ( ( ( NDE10 ) <= ( est_drct.1.runze.cv[loop.run,loop.n] + qnorm(0.975)*se_drct.1.runze.cv[loop.run, loop.n] ) ) & ( ( NDE10 ) >= ( est_drct.1.runze.cv[loop.run,loop.n] - qnorm(0.975)*se_drct.1.runze.cv[loop.run, loop.n] ) )  ){
      cp_drct.10.runze.cv[loop.run, loop.n] = 1
    }
    
    if ( ( ( NDE20 ) <= ( est_drct.2.runze.cv[loop.run,loop.n] + qnorm(0.975)*se_drct.2.runze.cv[loop.run, loop.n] ) ) & ( ( NDE20 ) >= ( est_drct.2.runze.cv[loop.run,loop.n] - qnorm(0.975)*se_drct.2.runze.cv[loop.run, loop.n] ) )  ){
      cp_drct.20.runze.cv[loop.run, loop.n] = 1
    }
    
    if ( ( ( NDE11 ) <= ( est_drct.1.runze.cv[loop.run,loop.n] + qnorm(0.975)*se_drct.1.runze.cv[loop.run, loop.n] ) ) & ( ( NDE11 ) >= ( est_drct.1.runze.cv[loop.run,loop.n] - qnorm(0.975)*se_drct.1.runze.cv[loop.run, loop.n] ) )  ){
      cp_drct.11.runze.cv[loop.run, loop.n] = 1
    }
    
    if ( ( ( NDE21 ) <= ( est_drct.2.runze.cv[loop.run,loop.n] + qnorm(0.975)*se_drct.2.runze.cv[loop.run, loop.n] ) ) & ( ( NDE21 ) >= ( est_drct.2.runze.cv[loop.run,loop.n] - qnorm(0.975)*se_drct.2.runze.cv[loop.run, loop.n] ) )  ){
      cp_drct.21.runze.cv[loop.run, loop.n] = 1
    }
    
    # power
    if ( ( 0 > ( est_indrct.1.runze.cv[loop.run,loop.n] + qnorm(0.975)*se_indrct.1.runze.cv[loop.run, loop.n] ) ) | ( 0 < ( est_indrct.1.runze.cv[loop.run,loop.n] - qnorm(0.975)*se_indrct.1.runze.cv[loop.run, loop.n] ) )  ){
      reject_indrct.10.runze.cv[loop.run, loop.n] = 1
      reject_indrct.11.runze.cv[loop.run, loop.n] = 1
    }
    
    if ( ( 0 > ( est_indrct.2.runze.cv[loop.run,loop.n] + qnorm(0.975)*se_indrct.2.runze.cv[loop.run, loop.n] ) ) | ( 0 < ( est_indrct.2.runze.cv[loop.run,loop.n] - qnorm(0.975)*se_indrct.2.runze.cv[loop.run, loop.n] ) )  ){
      reject_indrct.20.runze.cv[loop.run, loop.n] = 1
      reject_indrct.21.runze.cv[loop.run, loop.n] = 1
    }
    
    #
    if ( ( 0 > ( est_drct.1.runze.cv[loop.run,loop.n] + qnorm(0.975)*se_drct.1.runze.cv[loop.run, loop.n] ) ) | ( 0 < ( est_drct.1.runze.cv[loop.run,loop.n] - qnorm(0.975)*se_drct.1.runze.cv[loop.run, loop.n] ) )  ){
      reject_drct.10.runze.cv[loop.run, loop.n] = 1
      reject_drct.11.runze.cv[loop.run, loop.n] = 1
    }
    
    if ( ( 0 > ( est_drct.2.runze.cv[loop.run,loop.n] + qnorm(0.975)*se_drct.2.runze.cv[loop.run, loop.n] ) ) | ( 0 < ( est_drct.2.runze.cv[loop.run,loop.n] - qnorm(0.975)*se_drct.2.runze.cv[loop.run, loop.n] ) )  ){
      reject_drct.20.runze.cv[loop.run, loop.n] = 1
      reject_drct.21.runze.cv[loop.run, loop.n] = 1
    }
    
    
    
    
    ### HBIC ###
    
    # First calculate hbic corresponding to each lambda
    lambda.seq.runze = fitp.runze$lambda
    
    hbic.runze = c()
    for (loop.lambda in 1:length(lambda.seq.runze)) {
      sigma2_Yhat.tmp =  sum( ( Y - cbind(one.n, W_outcome_runze) %*% fitp.runze$fit$beta[,loop.lambda] )^2 ) / n 
      theta.G.tmp = (fitp.runze$fit$beta[,loop.lambda])[(1+q+r+1):(1+q+r+p)]
      df = 1+q+r+length(which(theta.G.tmp!= 0))
      hbic.runze[loop.lambda] = log( sigma2_Yhat.tmp ) +  df*log(log(n))*log(1+q+r+p)/n
    }
    
    
    hbic.min.runze = which.min(hbic.runze)
    
    theta0.hat.runze.hbic=fitp.runze$fit$beta[,hbic.min.runze][1]
    
    theta.T.hat.runze.hbic=fitp.runze$fit$beta[,hbic.min.runze][(1+1): (1+q)]
    
    theta.X.hat.runze.hbic=fitp.runze$fit$beta[,hbic.min.runze][(1+q+1) : (1+q+r)]
    
    theta.G.hat.runze.hbic=fitp.runze$fit$beta[,hbic.min.runze][(1+q+r+1):(1+q+r+p)]
    
    names(theta0.hat.runze.hbic)=NULL; names(theta.T.hat.runze.hbic)=NULL; names(theta.X.hat.runze.hbic)=NULL 
    names(theta.G.hat.runze.hbic)=NULL
    
    # step 2 mediation model
    
    SE_smple.runze.hbic = Runze_SE_sample_noX(Y = Y, Trt = cbind(TRT,X), Mediator = G, 
                                              theta0_hat = theta0.hat.runze.hbic,
                                              theta.T_hat = c(theta.T.hat.runze.hbic, theta.X.hat.runze.hbic),
                                              theta.G_hat = theta.G.hat.runze.hbic)
    
    
    est_indrct.1.runze.hbic[loop.run,loop.n] = SE_smple.runze.hbic$NIE_est["T1", 1] 
    est_indrct.2.runze.hbic[loop.run,loop.n] = SE_smple.runze.hbic$NIE_est["T2", 1]
    
    se_indrct.1.runze.hbic[loop.run,loop.n]= SE_smple.runze.hbic$NIE_SE["T1", 1] 
    se_indrct.2.runze.hbic[loop.run,loop.n]= SE_smple.runze.hbic$NIE_SE["T2", 1] 
    
    est_drct.1.runze.hbic[loop.run, loop.n] = SE_smple.runze.hbic$NDE_est["T1", 1]
    est_drct.2.runze.hbic[loop.run, loop.n] = SE_smple.runze.hbic$NDE_est["T2", 1]
    
    se_drct.1.runze.hbic[loop.run, loop.n] = SE_smple.runze.hbic$NDE_SE["T1", 1]
    se_drct.1.runze.hbic[loop.run, loop.n] = SE_smple.runze.hbic$NDE_SE["T2", 1]
    
    # cp
    if ( ( ( NIE10 ) <= ( est_indrct.1.runze.hbic[loop.run,loop.n] + qnorm(0.975)*se_indrct.1.runze.hbic[loop.run, loop.n] ) ) & ( ( NIE10 ) >= ( est_indrct.1.runze.hbic[loop.run,loop.n] - qnorm(0.975)*se_indrct.1.runze.hbic[loop.run, loop.n] ) )  ){
      cp_indrct.10.runze.hbic[loop.run, loop.n] = 1
    }
    
    if ( ( ( NIE20) <= ( est_indrct.2.runze.hbic[loop.run,loop.n] + qnorm(0.975)*se_indrct.2.runze.hbic[loop.run, loop.n] ) ) & ( ( NIE20 ) >= ( est_indrct.2.runze.hbic[loop.run,loop.n] - qnorm(0.975)*se_indrct.2.runze.hbic[loop.run, loop.n] ) )  ){
      cp_indrct.20.runze.hbic[loop.run, loop.n] = 1
    }
    
    
    if ( ( ( NIE11 ) <= ( est_indrct.1.runze.hbic[loop.run,loop.n] + qnorm(0.975)*se_indrct.1.runze.hbic[loop.run, loop.n] ) ) & ( ( NIE11 ) >= ( est_indrct.1.runze.hbic[loop.run,loop.n] - qnorm(0.975)*se_indrct.1.runze.hbic[loop.run, loop.n] ) )  ){
      cp_indrct.11.runze.hbic[loop.run, loop.n] = 1
    }
    
    if ( ( ( NIE21 ) <= ( est_indrct.2.runze.hbic[loop.run,loop.n] + qnorm(0.975)*se_indrct.2.runze.hbic[loop.run, loop.n] ) ) & ( ( NIE21 ) >= ( est_indrct.2.runze.hbic[loop.run,loop.n] - qnorm(0.975)*se_indrct.2.runze.hbic[loop.run, loop.n] ) )  ){
      cp_indrct.21.runze.hbic[loop.run, loop.n] = 1
    }
    
    #
    if ( ( ( NDE10 ) <= ( est_drct.1.runze.hbic[loop.run,loop.n] + qnorm(0.975)*se_drct.1.runze.hbic[loop.run, loop.n] ) ) & ( ( NDE10 ) >= ( est_drct.1.runze.hbic[loop.run,loop.n] - qnorm(0.975)*se_drct.1.runze.hbic[loop.run, loop.n] ) )  ){
      cp_drct.10.runze.hbic[loop.run, loop.n] = 1
    }
    
    if ( ( ( NDE20 ) <= ( est_drct.2.runze.hbic[loop.run,loop.n] + qnorm(0.975)*se_drct.2.runze.hbic[loop.run, loop.n] ) ) & ( ( NDE20 ) >= ( est_drct.2.runze.hbic[loop.run,loop.n] - qnorm(0.975)*se_drct.2.runze.hbic[loop.run, loop.n] ) )  ){
      cp_drct.20.runze.hbic[loop.run, loop.n] = 1
    }
    
    if ( ( ( NDE11 ) <= ( est_drct.1.runze.hbic[loop.run,loop.n] + qnorm(0.975)*se_drct.1.runze.hbic[loop.run, loop.n] ) ) & ( ( NDE11 ) >= ( est_drct.1.runze.hbic[loop.run,loop.n] - qnorm(0.975)*se_drct.1.runze.hbic[loop.run, loop.n] ) )  ){
      cp_drct.11.runze.hbic[loop.run, loop.n] = 1
    }
    
    if ( ( ( NDE21 ) <= ( est_drct.2.runze.hbic[loop.run,loop.n] + qnorm(0.975)*se_drct.2.runze.hbic[loop.run, loop.n] ) ) & ( ( NDE21 ) >= ( est_drct.2.runze.hbic[loop.run,loop.n] - qnorm(0.975)*se_drct.2.runze.hbic[loop.run, loop.n] ) )  ){
      cp_drct.21.runze.hbic[loop.run, loop.n] = 1
    }
    
    # power
    if ( ( 0 > ( est_indrct.1.runze.hbic[loop.run,loop.n] + qnorm(0.975)*se_indrct.1.runze.hbic[loop.run, loop.n] ) ) | ( 0 < ( est_indrct.1.runze.hbic[loop.run,loop.n] - qnorm(0.975)*se_indrct.1.runze.hbic[loop.run, loop.n] ) )  ){
      reject_indrct.10.runze.hbic[loop.run, loop.n] = 1
      reject_indrct.11.runze.hbic[loop.run, loop.n] = 1
    }
    
    if ( ( 0 > ( est_indrct.2.runze.hbic[loop.run,loop.n] + qnorm(0.975)*se_indrct.2.runze.hbic[loop.run, loop.n] ) ) | ( 0 < ( est_indrct.2.runze.hbic[loop.run,loop.n] - qnorm(0.975)*se_indrct.2.runze.hbic[loop.run, loop.n] ) )  ){
      reject_indrct.20.runze.hbic[loop.run, loop.n] = 1
      reject_indrct.21.runze.hbic[loop.run, loop.n] = 1
    }
    
    #
    if ( ( 0 > ( est_drct.1.runze.hbic[loop.run,loop.n] + qnorm(0.975)*se_drct.1.runze.hbic[loop.run, loop.n] ) ) | ( 0 < ( est_drct.1.runze.hbic[loop.run,loop.n] - qnorm(0.975)*se_drct.1.runze.hbic[loop.run, loop.n] ) )  ){
      reject_drct.10.runze.hbic[loop.run, loop.n] = 1
      reject_drct.11.runze.hbic[loop.run, loop.n] = 1
    }
    
    if ( ( 0 > ( est_drct.2.runze.hbic[loop.run,loop.n] + qnorm(0.975)*se_drct.2.runze.hbic[loop.run, loop.n] ) ) | ( 0 < ( est_drct.2.runze.hbic[loop.run,loop.n] - qnorm(0.975)*se_drct.2.runze.hbic[loop.run, loop.n] ) )  ){
      reject_drct.20.runze.hbic[loop.run, loop.n] = 1
      reject_drct.21.runze.hbic[loop.run, loop.n] = 1
    }
    
    
    
    print(loop.run)
    
  }
  print(loop.n)
}




plot(lambda_scad.cv[,1], lambda_scad.hbic[,1])
abline(a=0,b=1,col=2)


cp_pwr_cv = matrix(NA, 8, 4)
colnames(cp_pwr_cv) = c("cp proposed", "cp Runze", "power proposed", "power Runze")
rownames(cp_pwr_cv) = c("NIE10", "NIE11", "NIE20", "NIE21", "NDE10", "NDE11", "NDE20", "NDE21")
cp_pwr_cv["NIE10", ] = c( mean(cp_indrct.10.scad.cv), mean(cp_indrct.10.runze.cv), 
                          mean(reject_indrct.10.scad.cv), mean(reject_indrct.10.runze.cv)  )
cp_pwr_cv["NIE11", ] = c( mean(cp_indrct.11.scad.cv), mean(cp_indrct.11.runze.cv), 
                          mean(reject_indrct.11.scad.cv), mean(reject_indrct.11.runze.cv)  )
cp_pwr_cv["NIE20", ] = c( mean(cp_indrct.20.scad.cv), mean(cp_indrct.20.runze.cv), 
                          mean(reject_indrct.20.scad.cv), mean(reject_indrct.20.runze.cv)  )
cp_pwr_cv["NIE21", ] = c( mean(cp_indrct.21.scad.cv), mean(cp_indrct.21.runze.cv), 
                          mean(reject_indrct.21.scad.cv), mean(reject_indrct.21.runze.cv)  )
cp_pwr_cv["NDE10", ] = c( mean(cp_drct.10.scad.cv), mean(cp_drct.10.runze.cv), 
                          mean(reject_drct.10.scad.cv), mean(reject_drct.10.runze.cv)  )
cp_pwr_cv["NDE11", ] = c( mean(cp_drct.11.scad.cv), mean(cp_drct.11.runze.cv), 
                          mean(reject_drct.11.scad.cv), mean(reject_drct.11.runze.cv)  )
cp_pwr_cv["NDE20", ] = c( mean(cp_drct.20.scad.cv), mean(cp_drct.20.runze.cv), 
                          mean(reject_drct.20.scad.cv), mean(reject_drct.20.runze.cv)  )
cp_pwr_cv["NDE21", ] = c( mean(cp_drct.21.scad.cv), mean(cp_drct.21.runze.cv), 
                          mean(reject_drct.21.scad.cv), mean(reject_drct.21.runze.cv)  )


cp_pwr_cv



cp_pwr_hbic = matrix(NA, 8, 4)
colnames(cp_pwr_hbic) = c("cp proposed", "cp Runze", "power proposed", "power Runze")
rownames(cp_pwr_hbic) = c("NIE10", "NIE11", "NIE20", "NIE21", "NDE10", "NDE11", "NDE20", "NDE21")
cp_pwr_hbic["NIE10", ] = c( mean(cp_indrct.10.scad.hbic), mean(cp_indrct.10.runze.hbic), 
                            mean(reject_indrct.10.scad.hbic), mean(reject_indrct.10.runze.hbic)  )
cp_pwr_hbic["NIE11", ] = c( mean(cp_indrct.11.scad.hbic), mean(cp_indrct.11.runze.hbic), 
                            mean(reject_indrct.11.scad.hbic), mean(reject_indrct.11.runze.hbic)  )
cp_pwr_hbic["NIE20", ] = c( mean(cp_indrct.20.scad.hbic), mean(cp_indrct.20.runze.hbic), 
                            mean(reject_indrct.20.scad.hbic), mean(reject_indrct.20.runze.hbic)  )
cp_pwr_hbic["NIE21", ] = c( mean(cp_indrct.21.scad.hbic), mean(cp_indrct.21.runze.hbic), 
                            mean(reject_indrct.21.scad.hbic), mean(reject_indrct.21.runze.hbic)  )
cp_pwr_hbic["NDE10", ] = c( mean(cp_drct.10.scad.hbic), mean(cp_drct.10.runze.hbic), 
                            mean(reject_drct.10.scad.hbic), mean(reject_drct.10.runze.hbic)  )
cp_pwr_hbic["NDE11", ] = c( mean(cp_drct.11.scad.hbic), mean(cp_drct.11.runze.hbic), 
                            mean(reject_drct.11.scad.hbic), mean(reject_drct.11.runze.hbic)  )
cp_pwr_hbic["NDE20", ] = c( mean(cp_drct.20.scad.hbic), mean(cp_drct.20.runze.hbic), 
                            mean(reject_drct.20.scad.hbic), mean(reject_drct.20.runze.hbic)  )
cp_pwr_hbic["NDE21", ] = c( mean(cp_drct.21.scad.hbic), mean(cp_drct.21.runze.hbic), 
                            mean(reject_drct.21.scad.hbic), mean(reject_drct.21.runze.hbic)  )


cp_pwr_hbic




save(list = c("A", "B", "beta0", "beta1", "beta2", "cp_pwr_cv", "cp_pwr_hbic",
              "lambda_scad.cv", "lambda_scad.hbic", "mu.X", "NDE10", "NDE11", 
              "NDE20", "NDE21", "Sigma0.5_G_gl", "sigma.X", "A_D1", "A_D2", "A_G", "n.seq", 
              "NIE10", "NIE11", "NIE20", "NIE21", "p", "p.fac", "p.fac.runze", "q", "r", "run", "sigma_Y",
              "theta", "theta.D1", "theta.D2", "theta.G", "theta.T", "theta.X", "theta0", "vtheta.D"),
     file = "p=2000 n=120 glasso t8epsG NormalepsY.Rdata")








