rm(list = ls())

for(simnum in 1:50) {
  load(paste0("cD_p500_1_", simnum, ".Rdata"))
}



library(gdata)

max(abs(result1$est_drct.runze -result5$est_drct.runze))
max(abs(result1$B -result5$B))


run=500


# copy directly from the .R file in server
set.seed(12301)


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







max(abs(B - result1$B))
max(abs(beta0 - result10$beta0))


each.run=10

for (simnum in 1:50) {
  
  mv(from = paste0("result", simnum), to = "tmp")
  
  #### scad
  bias_theta.scad[ ((simnum-1)*each.run+1):(simnum*each.run) ,]=tmp$bias_theta.scad
  sr_theta.scad[ ((simnum-1)*each.run+1):(simnum*each.run) ,]=tmp$sr_theta.scad
  
  est_drct.10.scad[ ((simnum-1)*each.run+1):(simnum*each.run) ,]=tmp$est_drct.10.scad
  est_drct.11.scad[ ((simnum-1)*each.run+1):(simnum*each.run) ,]=tmp$est_drct.11.scad
  SE_drct.10.scad[ ((simnum-1)*each.run+1):(simnum*each.run) ,]=tmp$SE_drct.10.scad
  SE_drct.11.scad[ ((simnum-1)*each.run+1):(simnum*each.run) ,]=tmp$SE_drct.11.scad
  cp_drct.10.scad[ ((simnum-1)*each.run+1):(simnum*each.run) ,]=tmp$cp_drct.10.scad
  cp_drct.11.scad[ ((simnum-1)*each.run+1):(simnum*each.run) ,]=tmp$cp_drct.11.scad
  reject_drct.10.scad[ ((simnum-1)*each.run+1):(simnum*each.run) ,]=tmp$reject_drct.10.scad
  reject_drct.11.scad[ ((simnum-1)*each.run+1):(simnum*each.run) ,]=tmp$reject_drct.11.scad
  
  est_indrct.10.scad[ ((simnum-1)*each.run+1):(simnum*each.run) ,] = tmp$est_indrct.10.scad
  est_indrct.11.scad[ ((simnum-1)*each.run+1):(simnum*each.run) ,] = tmp$est_indrct.11.scad
  SE_indrct.10.scad[ ((simnum-1)*each.run+1):(simnum*each.run) ,] = tmp$SE_indrct.10.scad
  SE_indrct.11.scad[ ((simnum-1)*each.run+1):(simnum*each.run) ,] = tmp$SE_indrct.11.scad
  cp_indrct.10.scad[ ((simnum-1)*each.run+1):(simnum*each.run) ,] = tmp$cp_indrct.10.scad
  cp_indrct.11.scad[ ((simnum-1)*each.run+1):(simnum*each.run) ,] = tmp$cp_indrct.11.scad
  reject_indrct.10.scad[ ((simnum-1)*each.run+1):(simnum*each.run) ,] = tmp$reject_indrct.10.scad
  reject_indrct.11.scad[ ((simnum-1)*each.run+1):(simnum*each.run) ,] = tmp$reject_indrct.11.scad
  
  ##### Oracle
  est_drct.10.oracle[ ((simnum-1)*each.run+1):(simnum*each.run) ,] = tmp$est_drct.10.oracle
  est_drct.11.oracle[ ((simnum-1)*each.run+1):(simnum*each.run) ,] = tmp$est_drct.11.oracle
  SE_drct.10.oracle[ ((simnum-1)*each.run+1):(simnum*each.run) ,] = tmp$SE_drct.10.oracle
  SE_drct.11.oracle[ ((simnum-1)*each.run+1):(simnum*each.run) ,] = tmp$SE_drct.11.oracle
  cp_drct.10.oracle[ ((simnum-1)*each.run+1):(simnum*each.run) ,] = tmp$cp_drct.10.oracle
  cp_drct.11.oracle[ ((simnum-1)*each.run+1):(simnum*each.run) ,] = tmp$cp_drct.11.oracle
  reject_drct.10.oracle[ ((simnum-1)*each.run+1):(simnum*each.run) ,] = tmp$reject_drct.10.oracle
  reject_drct.11.oracle[ ((simnum-1)*each.run+1):(simnum*each.run) ,] = tmp$reject_drct.11.oracle
  
  est_indrct.10.oracle[ ((simnum-1)*each.run+1):(simnum*each.run) ,] = tmp$est_indrct.10.oracle
  est_indrct.11.oracle[ ((simnum-1)*each.run+1):(simnum*each.run) ,] = tmp$est_indrct.11.oracle
  SE_indrct.10.oracle[ ((simnum-1)*each.run+1):(simnum*each.run) ,] = tmp$SE_indrct.10.oracle
  SE_indrct.11.oracle[ ((simnum-1)*each.run+1):(simnum*each.run) ,] = tmp$SE_indrct.11.oracle
  cp_indrct.10.oracle[ ((simnum-1)*each.run+1):(simnum*each.run) ,] = tmp$cp_indrct.10.oracle
  cp_indrct.11.oracle[ ((simnum-1)*each.run+1):(simnum*each.run) ,] = tmp$cp_indrct.11.oracle
  reject_indrct.10.oracle[ ((simnum-1)*each.run+1):(simnum*each.run) ,] = tmp$reject_indrct.10.oracle
  reject_indrct.11.oracle[ ((simnum-1)*each.run+1):(simnum*each.run) ,] = tmp$reject_indrct.11.oracle
  
  
  ##### Runze
  est_drct.runze[ ((simnum-1)*each.run+1):(simnum*each.run) ,] = tmp$est_drct.runze
  SE_drct.runze[ ((simnum-1)*each.run+1):(simnum*each.run) ,] = tmp$SE_drct.runze
  cp_drct10.runze[ ((simnum-1)*each.run+1):(simnum*each.run) ,] = tmp$cp_drct10.runze
  cp_drct11.runze[ ((simnum-1)*each.run+1):(simnum*each.run) ,] = tmp$cp_drct11.runze
  reject_drct.runze[ ((simnum-1)*each.run+1):(simnum*each.run) ,] = tmp$reject_drct.runze
  
  est_indrct.runze[ ((simnum-1)*each.run+1):(simnum*each.run) ,] = tmp$est_indrct.runze
  SE_indrct.runze[ ((simnum-1)*each.run+1):(simnum*each.run) ,] = tmp$SE_indrct.runze
  cp_indrct10.runze[ ((simnum-1)*each.run+1):(simnum*each.run) ,] = tmp$cp_indrct10.runze
  cp_indrct11.runze[ ((simnum-1)*each.run+1):(simnum*each.run) ,] = tmp$cp_indrct11.runze
  reject_indrct.runze[ ((simnum-1)*each.run+1):(simnum*each.run) ,] = tmp$reject_indrct.runze
  
  
  
  
  ##### Lasso
  bias_theta.lasso[ ((simnum-1)*each.run+1):(simnum*each.run) ,] = tmp$bias_theta.lasso
  sr_theta.lasso[ ((simnum-1)*each.run+1):(simnum*each.run) ,] = tmp$sr_theta.lasso
  
  est_drct.10.lasso[ ((simnum-1)*each.run+1):(simnum*each.run) ,] = tmp$est_drct.10.lasso
  est_drct.11.lasso[ ((simnum-1)*each.run+1):(simnum*each.run) ,] = tmp$est_drct.11.lasso
  bias_drct.10.lasso[ ((simnum-1)*each.run+1):(simnum*each.run) ,] = tmp$bias_drct.10.lasso
  bias_drct.11.lasso[ ((simnum-1)*each.run+1):(simnum*each.run) ,] = tmp$bias_drct.11.lasso
  
  est_indrct.10.lasso[ ((simnum-1)*each.run+1):(simnum*each.run) ,] = tmp$est_indrct.10.lasso
  est_indrct.11.lasso[ ((simnum-1)*each.run+1):(simnum*each.run) ,] = tmp$est_indrct.11.lasso
  bias_indrct.10.lasso[ ((simnum-1)*each.run+1):(simnum*each.run) ,] = tmp$bias_indrct.10.lasso
  bias_indrct.11.lasso[ ((simnum-1)*each.run+1):(simnum*each.run) ,] = tmp$bias_indrct.11.lasso
  
  
  
  
  ##### Dave
  est_drct.dave.cY[ ((simnum-1)*each.run+1):(simnum*each.run) ,] = tmp$est_drct.dave.cY
  SE_drct.dave.cY[ ((simnum-1)*each.run+1):(simnum*each.run) ,] = tmp$SE_drct.dave.cY
  cp_drct10.dave.cY[ ((simnum-1)*each.run+1):(simnum*each.run) ,] = tmp$cp_drct10.dave.cY
  cp_drct11.dave.cY[ ((simnum-1)*each.run+1):(simnum*each.run) ,] = tmp$cp_drct11.dave.cY
  reject_drct.dave.cY[ ((simnum-1)*each.run+1):(simnum*each.run) ,] = tmp$reject_drct.dave.cY
  
  est_indrct.dave.cY[ ((simnum-1)*each.run+1):(simnum*each.run) ,] = tmp$est_indrct.dave.cY
  SE_indrct.dave.cY[ ((simnum-1)*each.run+1):(simnum*each.run) ,] = tmp$SE_indrct.dave.cY
  cp_indrct10.dave.cY[ ((simnum-1)*each.run+1):(simnum*each.run) ,] = tmp$cp_indrct10.dave.cY
  cp_indrct11.dave.cY[ ((simnum-1)*each.run+1):(simnum*each.run) ,] = tmp$cp_indrct11.dave.cY
  reject_indrct.dave.cY[ ((simnum-1)*each.run+1):(simnum*each.run) ,] = tmp$reject_indrct.dave.cY
  
  ### Huang
  est_indrct.10.Huang[ ((simnum-1)*each.run+1):(simnum*each.run) ,] = tmp$est_indrct.10.Huang
  est_indrct.11.Huang[ ((simnum-1)*each.run+1):(simnum*each.run) ,] = tmp$est_indrct.11.Huang
  
  SE_indrct.10.Huang[ ((simnum-1)*each.run+1):(simnum*each.run) ,] = tmp$SE_indrct.10.Huang
  SE_indrct.11.Huang[ ((simnum-1)*each.run+1):(simnum*each.run) ,] = tmp$SE_indrct.11.Huang
  cp_indrct.10.Huang[ ((simnum-1)*each.run+1):(simnum*each.run) ,] = tmp$cp_indrct.10.Huang
  cp_indrct.11.Huang[ ((simnum-1)*each.run+1):(simnum*each.run) ,] = tmp$cp_indrct.11.Huang
  reject_indrct.10.Huang[ ((simnum-1)*each.run+1):(simnum*each.run) ,] = tmp$reject_indrct.10.Huang
  reject_indrct.11.Huang[ ((simnum-1)*each.run+1):(simnum*each.run) ,] = tmp$reject_indrct.11.Huang
  
  SE_indrct.10.Huang.boot[ ((simnum-1)*each.run+1):(simnum*each.run) ,] = tmp$SE_indrct.10.Huang.boot
  SE_indrct.11.Huang.boot[ ((simnum-1)*each.run+1):(simnum*each.run) ,] = tmp$SE_indrct.11.Huang.boot
  cp_indrct.10.Huang.boot[ ((simnum-1)*each.run+1):(simnum*each.run) ,] = tmp$cp_indrct.10.Huang.boot
  cp_indrct.11.Huang.boot[ ((simnum-1)*each.run+1):(simnum*each.run) ,] = tmp$cp_indrct.11.Huang.boot
  reject_indrct.10.Huang.boot[ ((simnum-1)*each.run+1):(simnum*each.run) ,] = tmp$reject_indrct.10.Huang.boot
  reject_indrct.11.Huang.boot[ ((simnum-1)*each.run+1):(simnum*each.run) ,] = tmp$reject_indrct.11.Huang.boot
  
  PCAbias_indrct.10.Huang[ ((simnum-1)*each.run+1):(simnum*each.run) ,] = tmp$PCAbias_indrct.10.Huang
  PCAbias_indrct.11.Huang[ ((simnum-1)*each.run+1):(simnum*each.run) ,] = tmp$PCAbias_indrct.11.Huang
  
  
  
  
  
  rm(tmp)
  
}

save(list = c( "run",
               "n",
               "p",
               "q",
               "one.n",
               "TRT",
               "theta0",
               "c.T",
               "theta.T",
               "c.G",
               "s.G",
               "theta.G",
               "cstar.seq",
               "s.D",
               "sigma_Y",
               "c_beta0",
               "beta0",
               "c_B",
               "B",
               "rho_G",
               "c_sigmaG",
               "Sigma_G",
               "p.fac",
               "p.fac_runze",
               "p.fac_lasso" ,
               "bias_theta.scad",
               "sr_theta.scad",
               "est_drct.10.scad",
               "est_drct.11.scad",
               "SE_drct.10.scad",
               "SE_drct.11.scad",
               "cp_drct.10.scad",
               "cp_drct.11.scad",
               "reject_drct.10.scad",
               "reject_drct.11.scad",
               "est_indrct.10.scad",
               "est_indrct.11.scad",
               "SE_indrct.10.scad",
               "SE_indrct.11.scad",
               "cp_indrct.10.scad",
               "cp_indrct.11.scad",
               "reject_indrct.10.scad",
               "reject_indrct.11.scad",
               "est_drct.10.oracle",
               "est_drct.11.oracle",
               "SE_drct.10.oracle",
               "SE_drct.11.oracle",
               "cp_drct.10.oracle",
               "cp_drct.11.oracle",
               "reject_drct.10.oracle",
               "reject_drct.11.oracle",
               "est_indrct.10.oracle",
               "est_indrct.11.oracle",
               "SE_indrct.10.oracle",
               "SE_indrct.11.oracle",
               "cp_indrct.10.oracle",
               "cp_indrct.11.oracle",
               "reject_indrct.10.oracle",
               "reject_indrct.11.oracle",
               "est_drct.runze",
               "SE_drct.runze",
               "cp_drct10.runze",
               "cp_drct11.runze",
               "reject_drct.runze",
               "est_indrct.runze",
               "SE_indrct.runze",
               "cp_indrct10.runze",
               "cp_indrct11.runze",
               "reject_indrct.runze",
               "bias_theta.lasso",
               "sr_theta.lasso",
               "est_drct.10.lasso",
               "est_drct.11.lasso",
               "bias_drct.10.lasso",
               "bias_drct.11.lasso",
               "est_indrct.10.lasso",
               "est_indrct.11.lasso",
               "bias_indrct.10.lasso",
               "bias_indrct.11.lasso",
               "est_drct.dave.cY",
               "SE_drct.dave.cY",
               "cp_drct10.dave.cY",
               "cp_drct11.dave.cY",
               "reject_drct.dave.cY",
               "est_indrct.dave.cY",
               "SE_indrct.dave.cY",
               "cp_indrct10.dave.cY",
               "cp_indrct11.dave.cY",
               "reject_indrct.dave.cY",
               "est_indrct.10.Huang",
               "est_indrct.11.Huang",
               "SE_indrct.10.Huang",
               "SE_indrct.11.Huang",
               "cp_indrct.10.Huang",
               "cp_indrct.11.Huang",
               "reject_indrct.10.Huang",
               "reject_indrct.11.Huang",
               "SE_indrct.10.Huang.boot",
               "SE_indrct.11.Huang.boot",
               "cp_indrct.10.Huang.boot",
               "cp_indrct.11.Huang.boot",
               "reject_indrct.10.Huang.boot",
               "reject_indrct.11.Huang.boot",
               "PCAbias_indrct.10.Huang",
               "PCAbias_indrct.11.Huang"
                 ),
       file = paste0("power cp q=1 1.1.1.Rdata") )
