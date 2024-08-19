library(ggplot2)

cstar.seq




############## NDE part ##################################################################################


boxplot(sqrt(n)*bias_theta.scad, main="SCAD sqrt(n)*bias l1 norm",
        names = cstar.seq)
boxplot(sr_theta.scad, main = "SCAD support recovery error",
        names = cstar.seq)







plot(x=cstar.seq, y = apply(cp_drct.10.scad, 2, mean), type = "b", col=1, pch=18, lwd=3,
     xlab = expression(c[D]),
     ylab = "Empirical Coverage Probability", ylim = c(0,1), yaxt="n",
     main = expression(zeta(T[1]==0) == 3.15 - 4.2~ c[D]) , cex=0.6)
lines(x=cstar.seq, y = apply(cp_drct.10.oracle, 2, mean), type = "b", col=2, pch=5, lwd=2, lty=2, cex = 0.6)
lines(x=cstar.seq, y= apply(cp_drct10.runze, 2, mean), type = "b", col=6, pch =5, lwd = 2, cex=0.6)
lines(x=cstar.seq, y= apply(cp_drct10.dave.cY, 2, mean), type = "b", col=4, pch =5, lwd = 2, cex=0.6)
#lines(x=cstar.seq, y= apply(cp_drct10.dave, 2, mean), type = "b", col="gray", pch =5, lwd = 2, cex=0.6)
#lines(x=cstar.seq, y= apply(cp_drct.10.Huang, 2, mean), type = "b", col=7, pch =5, lwd = 2, cex=0.6)
abline(h=0.95, lty=3)
abline(h=0, lty=3)
axis(side = 2, at = c(0,  0.5, 0.95))
legend("right", legend = c("Proposed", 
                           "Oracle", "Runze", "Dave"),
       pch=c(18, 5, 5,5,5), col=c(1, 2, 6, 4), 
       lty=c(1,2,1, 1), lwd = c(3,2,2,2),
       cex = 0.8)

plot(x=cstar.seq[c(1, 2, 7, 12, 13, 14, 24, 29, 30)], 
     y = apply(cp_drct.10.scad[,c(1, 2, 7, 12, 13, 14, 24, 29, 30) ], 2, mean), 
     type = "b", col=1, pch=18, lwd=3,
     xlab = expression(c[D]),
     ylab = "Empirical Coverage Probability", ylim = c(0,1), yaxt="n",
     main = expression(zeta(T[1]==0) == 3.15 - 4.2~ c[D]) , cex=0.6)
lines(x=cstar.seq[c(1, 2, 7, 12, 13, 14, 24, 29, 30)], 
      y = apply(cp_drct.10.oracle[,c(1, 2, 7, 12, 13, 14, 24, 29, 30) ], 2, mean), 
      type = "b", col=2, pch=5, lwd=2, lty=2, cex = 0.6)
lines(x=cstar.seq[c(1, 2, 7, 12, 13, 14, 24, 29, 30)], 
      y= apply(cp_drct10.runze[,c(1, 2, 7, 12, 13, 14, 24, 29, 30) ], 2, mean), 
      type = "b", col=6, pch =5, lwd = 2, cex=0.6)
lines(x=cstar.seq[c(1, 2, 7, 12, 13, 14, 24, 29, 30)], 
      y= apply(cp_drct10.dave.cY[,c(1, 2, 7, 12, 13, 14, 24, 29, 30) ], 2, mean), 
      type = "b", col=4, pch =5, lwd = 2, cex=0.6)
abline(h=0.95, lty=3)
abline(h=0, lty=3)
axis(side = 2, at = c(0,  0.5, 0.95))
legend("right", legend = c("Proposed", 
                           "Oracle", "Runze", "Dave"),
       pch=c(18, 5, 5,5,5), col=c(1, 2, 6, 4), 
       lty=c(1,2,1, 1), lwd = c(3,2,2,2),
       cex = 0.8)



plot(x=cstar.seq, y = apply(cp_drct.11.scad, 2, mean), type = "b", col=1, pch=18, lwd=3,
     xlab = expression(c[D]),
     ylab = "Empirical Coverage Probability", ylim = c(0,1), yaxt="n",
     main = expression( zeta(T[1]==1) == 3.15 ) , cex=0.6)
lines(x=cstar.seq, y = apply(cp_drct.11.oracle, 2, mean), type = "b", col=2, pch=5, lwd=2, lty=2, cex = 0.6)
lines(x=cstar.seq, y= apply(cp_drct11.runze, 2, mean), type = "b", col=6, pch =5, lwd = 2, cex=0.6)
lines(x=cstar.seq, y= apply(cp_drct11.dave.cY, 2, mean), type = "b", col=4, pch =5, lwd = 2, cex=0.6)
#lines(x=cstar.seq, y= apply(cp_drct11.dave, 2, mean), type = "b", col="gray", pch =5, lwd = 2, cex=0.6)
#lines(x=cstar.seq, y= apply(cp_drct.11.Huang, 2, mean), type = "b", col=7, pch =5, lwd = 2, cex=0.6)
abline(h=0.95, lty=3)
abline(h=0, lty=3)
axis(side = 2, at = c(0,  0.5, 0.95))
legend("left", legend = c("Proposed", 
                          "Oracle", "Runze", "Dave"),
       pch=c(18, 5, 5,5,5), col=c(1, 2, 6, 4), 
       lty=c(1,2,1, 1), lwd = c(3,2,2,2),
       cex = 0.8)



plot(x=cstar.seq[c(1, 2, 7, 12, 13, 14, 24, 29, 30)], 
     y = apply(cp_drct.11.scad[,c(1, 2, 7, 12, 13, 14, 24, 29, 30)], 2, mean), 
     type = "b", col=1, pch=18, lwd=3,
     xlab = expression(c[D]),
     ylab = "Empirical Coverage Probability", ylim = c(0,1), yaxt="n",
     main = expression( zeta(T[1]==1) == 0.1 ) , cex=0.6)
lines(x=cstar.seq[c(1, 2, 7, 12, 13, 14, 24, 29, 30)], 
      y = apply(cp_drct.11.oracle[,c(1, 2, 7, 12, 13, 14, 24, 29, 30)], 2, mean), 
      type = "b", col=2, pch=5, lwd=2, lty=2, cex = 0.6)
lines(x=cstar.seq[c(1, 2, 7, 12, 13, 14, 24, 29, 30)], 
      y= apply(cp_drct11.runze[,c(1, 2, 7, 12, 13, 14, 24, 29, 30)], 2, mean), 
      type = "b", col=6, pch =5, lwd = 2, cex=0.6)
lines(x=cstar.seq[c(1, 2, 7, 12, 13, 14, 24, 29, 30)], 
      y= apply(cp_drct11.dave.cY[,c(1, 2, 7, 12, 13, 14, 24, 29, 30)], 2, mean), 
      type = "b", col=4, pch =5, lwd = 2, cex=0.6)
abline(h=0.95, lty=3)
abline(h=0, lty=3)
axis(side = 2, at = c(0,  0.5, 0.95))
legend("left", legend = c("Proposed", 
                          "Oracle", "Runze", "Dave"),
       pch=c(18, 5, 5,5,5), col=c(1, 2, 6, 4), 
       lty=c(1,2,1, 1), lwd = c(3,2,2,2),
       cex = 0.8)







plot(x=cstar.seq, y = apply(reject_drct.10.scad, 2, mean), type = "l", col=1, pch=18, lwd=3,
     xlab = expression(c[D]),
     ylab = "Empirical Power", ylim = c(0,1), yaxt="n",
     main = expression(zeta(T[1]==0) == 3.15 - 4.2~ c[D] ) , cex=0.6)
lines(x=cstar.seq, y = apply(reject_drct.10.oracle, 2, mean), type = "b", col=2, pch=5, lwd=2, lty=2, cex = 0.6)
lines(x=cstar.seq, y= apply(reject_drct.runze, 2, mean), type = "b", col=6, pch =5, lwd = 2, cex=0.6)
lines(x=cstar.seq, y= apply(reject_drct.dave.cY, 2, mean), type = "b", col=4, pch =5, lwd = 2, cex=0.6)
#lines(x=cstar.seq, y= apply(reject_drct.dave, 2, mean), type = "b", col="gray", pch =5, lwd = 2, cex=0.6)
#lines(x=cstar.seq, y= apply(reject_drct.10.Huang, 2, mean), type = "b", col=7, pch =5, lwd = 2, cex=0.6)
abline(h=0.95, lty=3)
abline(h=0, lty=3)
abline(h=0.05, lty=3)
axis(side = 2, at = c(0,  0.5, 0.95))
legend("left", legend = c("Proposed", 
                          "Oracle", "Runze", "Dave"),
       pch=c(18, 5, 5,5,5), col=c(1, 2, 6, 4), 
       lty=c(1,2,1, 1), lwd = c(3,2,2,2),
       cex = 0.8)

plot(x=cstar.seq[c(14:28)], 
     y = apply(reject_drct.10.scad[,c(14:28)], 2, mean), 
     type = "b", col=1, pch=18, lwd=3,
     xlab = expression(c[D]),
     ylab = "Empirical Power", ylim = c(0,1), yaxt="n",
     main = expression(zeta(T[1]==0) == 3.15 - 4.2~ c[D] ) , cex=0.6)
lines(x=cstar.seq[c(14:28)], 
      y = apply(reject_drct.10.oracle[,c(14:28)], 2, mean), 
      type = "b", col=2, pch=5, lwd=2, lty=2, cex = 0.6)
lines(x=cstar.seq[c(14:28)], 
      y= apply(reject_drct.runze[,c(14:28)], 2, mean), 
      type = "b", col=6, pch =5, lwd = 2, cex=0.6)
lines(x=cstar.seq[c(14:28)], 
      y= apply(reject_drct.dave.cY[,c(14:28)], 2, mean), 
      type = "b", col=4, pch =5, lwd = 2, cex=0.6)
abline(h=0.95, lty=3)
abline(h=0.05, lty=3)
abline(v=0.75, lty=3)
axis(side = 2, at = c(0.05,  0.5, 0.95))
axis(side = 1, at = c(round(0.75, digits = 3)))
legend("bottomleft", legend = c("Proposed", 
                          "Oracle", "Runze", "Dave"),
       pch=c(18, 5, 5,5,5), col=c(1, 2, 6, 4), 
       lty=c(1,2,1, 1), lwd = c(3,2,2,2),
       cex = 0.8)





plot(x=cstar.seq, y = apply(reject_drct.11.scad, 2, mean), type = "b", col=1, pch=18, lwd=3,
     xlab = expression(c[D]),
     ylab = "Empirical Power", ylim = c(0,1), yaxt="n",
     main = expression( zeta(T[1]==1) == 3.15 ) , cex=0.6)
lines(x=cstar.seq, y = apply(reject_drct.11.oracle, 2, mean), type = "b", col=2, pch=5, lwd=2, lty=2, cex = 0.6)
lines(x=cstar.seq, y= apply(reject_drct.runze, 2, mean), type = "b", col=6, pch =5, lwd = 2, cex=0.6)
lines(x=cstar.seq, y= apply(reject_drct.dave.cY, 2, mean), type = "b", col=4, pch =5, lwd = 2, cex=0.6)
#lines(x=cstar.seq, y= apply(reject_drct.dave, 2, mean), type = "b", col="gray", pch =5, lwd = 2, cex=0.6)
#lines(x=cstar.seq, y= apply(reject_drct.11.Huang, 2, mean), type = "b", col=7, pch =5, lwd = 2, cex=0.6)
abline(h=0.05, lty=3)
abline(h=0.95, lty=3)
axis(side = 2, at = c(0.05,  0.5, 0.95))
legend("left", legend = c("Proposed", 
                          "Oracle", "Runze", "Dave"),
       pch=c(18, 5, 5,5,5), col=c(1, 2, 6, 4), 
       lty=c(1,2,1, 1), lwd = c(3,2,2,2),
       cex = 0.8)

















### lasso ###
# compare bias from lasso and our method
# first calculate bias of our method
bias_drct.10.scad=matrix(0, nrow = run, ncol = length(cstar.seq))
bias_drct.11.scad=matrix(0, nrow = run, ncol = length(cstar.seq))
for (loop.c in 1:length(cstar.seq)) {
  
  c.D = cstar.seq[loop.c]
  Theta.D=c.D*matrix( c( c(1, 0.9, 0.8, 0.7, 0.6), 
                         c(1, 0.9, 0.8, 0.7, 0.6),
                         rep(0, p - s.D) ), 
                      ncol = q )
  NDE0 = as.numeric( theta.T + t(beta0)%*%Theta.D )
  NDE1 = as.numeric( theta.T + t(beta0 + B)%*%Theta.D  )
  
  print(paste0("c_D = ", c.D, ", NDE0=", NDE0, ", NDE1=", NDE1))
  
  
  for (loop.run in 1:run) {
    
    bias_drct.10.scad[loop.run, loop.c] = est_drct.10.scad[loop.run, loop.c] - NDE0
    bias_drct.11.scad[loop.run, loop.c] = est_drct.11.scad[loop.run, loop.c] - NDE1
    
  }
  
}


# draw the boxplot comparing the bias of lasso and scad
value=c()

# NDE(T=0)
for (i in 1:length(cstar.seq)) {
  value=c(value, (bias_drct.10.lasso[,i]), (bias_drct.10.scad[,i]))
}

cD_nor=rep(cstar.seq,each=2*run)
cD_nor=factor(cD_nor, 
              levels = as.character(cstar.seq) )

method=rep( rep(c("Lasso","SCAD"), each=run), times=length(cstar.seq) )

total=data.frame(value, cD_nor, method)


# NDE(T=0)
ggplot(total, aes(x=cD_nor, y=(value), fill=method)) +
  geom_boxplot() + scale_fill_manual(values=c("turquoise", "tomato"))  +
  labs(x=expression(c[D]),y=NULL,title = NULL)+ 
  theme(legend.title = element_blank()) + 
  ggtitle(expression( Bias~of~zeta(T[1]==0)==0.1 - 4.2~ c[D]  ))+
  geom_hline(yintercept=0,  color = "red")







value=c()

# NDE(T=1)
for (i in 1:length(cstar.seq)) {
  value=c(value, (bias_drct.11.lasso[,i]), (bias_drct.11.scad[,i]))
}



cD_nor=rep(cstar.seq,each=2*run)
cD_nor=factor(cD_nor, 
              levels = as.character(cstar.seq) )

method=rep( rep(c("Lasso","SCAD"), each=run), times=length(cstar.seq) )

total=data.frame(value, cD_nor, method)

# NDE(T=1)
ggplot(total, aes(x=cD_nor, y=(value), fill=method)) +
  geom_boxplot() + scale_fill_manual(values=c("turquoise", "tomato"))  +
  labs(x=expression(c[D]),y=NULL,title = NULL)+ 
  theme(legend.title = element_blank()) + 
  ggtitle(expression( Bias~of~zeta(T[1]==1)==0.1  ))+
  geom_hline(yintercept=0,  color = "red")










### compare Runze and Dave estimation bias
boxplot(est_drct.runze)
boxplot(est_drct.dave.cY)


value=c()

for (i in 1:length(cstar.seq)) {
  value=c(value, (est_drct.runze[,i]), (est_drct.dave.cY[,i]) )
}

cD_nor=rep(cstar.seq,each=( nrow(est_drct.runze) + nrow(est_drct.dave.cY) ) )
cD_nor=factor(cD_nor, 
              levels = as.character(cstar.seq) )

method=rep( c( rep("Runze", nrow(est_drct.runze)), 
               rep("Dave", nrow(est_drct.dave.cY)) ) , 
            times=length(cstar.seq) )

total=data.frame(value, cD_nor, method)



ggplot(total, aes(x=cD_nor, y=(value), fill=method)) +
  geom_boxplot() + scale_fill_manual(values=c("cyan", "pink"))  +
  labs(x=expression(c[D]),y=NULL,title = NULL)+ 
  theme(legend.title = element_blank()) + 
  ggtitle(expression(Point~Estimates~of~NDE~by~Runze~and~Dave))


boxplot(est_drct.dave.cY[,13], est_drct.runze[,13])
abline(h = theta.T, col=2)



















###### NIE part #####################################################################

cstar.seq



boxplot(sqrt(n)*bias_theta.scad, main="SCAD sqrt(n)*bias l1 norm")
boxplot(sr_theta.scad, main = "SCAD support recovery error")



plot(x=cstar.seq, y = apply(cp_indrct.10.scad, 2, mean), type = "b", col=1, pch=18, lwd=3,
     xlab = expression(c[D]),
     ylab = "Empirical Coverage Probability", ylim = c(0,1), yaxt="n",
     main = expression(delta(T[1]==0) == 4.2) , cex=0.6)
lines(x=cstar.seq, y = apply(cp_indrct.10.oracle, 2, mean), type = "b", col=2, pch=5, lwd=2, lty=2, cex = 0.6)
lines(x=cstar.seq, y= apply(cp_indrct10.runze, 2, mean), type = "b", col=6, pch =5, lwd = 2, cex=0.6)
lines(x=cstar.seq, y= apply(cp_indrct10.dave.cY, 2, mean), type = "b", col=4, pch =5, lwd = 2, cex=0.6)
#lines(x=cstar.seq, y= apply(cp_indrct10.dave, 2, mean), type = "b", col="gray", pch =5, lwd = 2, cex=0.6)
lines(x=cstar.seq, y= apply(cp_indrct.10.Huang.boot, 2, mean), type = "b", col=7, pch =5, lwd = 2, cex=0.6)
abline(h=0.95, lty=3)
abline(h=0, lty=3)
axis(side = 2, at = c(0,  0.5, 0.95))
legend("left", legend = c("Proposed", 
                           "Oracle", "Runze", "Dave", "Huang"),
       pch=c(18, 5, 5,5,5,5), col=c(1, 2, 6, 4, 7), 
       lty=c(1,2,1, 1, 1), lwd = c(3,2,2,2,2),
       cex = 0.8)


plot(x=cstar.seq[c(1, 2, 7, 12, 18, 26, 27, 28, 29)], 
     y = apply(cp_indrct.10.scad[,c(1, 2, 7, 12, 18, 26, 27, 28, 29)], 2, mean), 
     type = "b", col=1, pch=18, lwd=3,
     xlab = expression(c[D]),
     ylab = "Empirical Coverage Probability", ylim = c(0,1), yaxt="n",
     main = expression(delta(T[1]==0) == 4.2) , cex=0.6)
lines(x=cstar.seq[c(1, 2, 7, 12, 18, 26, 27, 28, 29)], 
      y = apply(cp_indrct.10.oracle[,c(1, 2, 7, 12, 18, 26, 27, 28, 29)], 2, mean), 
      type = "b", col=2, pch=5, lwd=2, lty=2, cex = 0.6)
lines(x=cstar.seq[c(1, 2, 7, 12, 18, 26, 27, 28, 29)], 
      y= apply(cp_indrct10.runze[,c(1, 2, 7, 12, 18, 26, 27, 28, 29)], 2, mean), 
      type = "b", col=6, pch =5, lwd = 2, cex=0.6)
lines(x=cstar.seq[c(1, 2, 7, 12, 18, 26, 27, 28, 29)], 
      y= apply(cp_indrct10.dave.cY[,c(1, 2, 7, 12, 18, 26, 27, 28, 29)], 2, mean), 
      type = "b", col=4, pch =5, lwd = 2, cex=0.6)
lines(x=cstar.seq[c(1, 2, 7, 12, 18, 26, 27, 28, 29)], 
      y= apply(cp_indrct.10.Huang.boot[,c(1, 2, 7, 12, 18, 26, 27, 28, 29)], 2, mean), 
      type = "b", col=7, pch =5, lwd = 2, cex=0.6)
abline(h=0.95, lty=3)
abline(h=0, lty=3)
axis(side = 2, at = c(0,  0.5, 0.95))
legend("left", legend = c("Proposed", 
                           "Oracle", "Runze", "Dave", "Huang"),
       pch=c(18, 5, 5,5,5,5), col=c(1, 2, 6, 4, 7), 
       lty=c(1,2,1, 1, 1), lwd = c(3,2,2,2,2),
       cex = 0.8)




plot(x=cstar.seq, y = apply(cp_indrct.11.scad, 2, mean), type = "b", col=1, pch=18, lwd=3,
     xlab = expression(c[D]),
     ylab = "Empirical Coverage Probability", ylim = c(0,1), yaxt="n",
     main = expression( delta(T[1]==1) == 4.2 (1 + c[D]) ) , cex=0.6)
lines(x=cstar.seq, y = apply(cp_indrct.11.oracle, 2, mean), type = "b", col=2, pch=5, lwd=2, lty=2, cex = 0.6)
lines(x=cstar.seq, y= apply(cp_indrct11.runze, 2, mean), type = "b", col=6, pch =5, lwd = 2, cex=0.6)
lines(x=cstar.seq, y= apply(cp_indrct11.dave.cY, 2, mean), type = "b", col=4, pch =5, lwd = 2, cex=0.6)
lines(x=cstar.seq, y= apply(cp_indrct.11.Huang.boot, 2, mean), type = "b", col=7, pch =5, lwd = 2, cex=0.6)
abline(h=0.95, lty=3)
abline(h=0, lty=3)
axis(side = 2, at = c(0,  0.5, 0.95))
legend("right", legend = c("Proposed", 
                          "Oracle", "Runze", "Dave", "Huang"),
       pch=c(18, 5, 5,5,5,5), col=c(1, 2, 6, 4, 7), 
       lty=c(1,2,1, 1, 1), lwd = c(3,2,2,2,2),
       cex = 0.8)




plot(x=cstar.seq[c(1, 2, 7, 12, 18, 26, 27, 28, 29)], 
     y = apply(cp_indrct.11.scad[,c(1, 2, 7, 12, 18, 26, 27, 28, 29)], 2, mean), 
     type = "b", col=1, pch=18, lwd=3,
     xlab = expression(c[D]),
     ylab = "Empirical Coverage Probability", ylim = c(0,1), yaxt="n",
     main = expression( delta(T[1]==1) == 4.2 (1 + c[D]) ) , cex=0.6)
lines(x=cstar.seq[c(1, 2, 7, 12, 18, 26, 27, 28, 29)], 
      y = apply(cp_indrct.11.oracle[,c(1, 2, 7, 12, 18, 26, 27, 28, 29)], 2, mean), 
      type = "b", col=2, pch=5, lwd=2, lty=2, cex = 0.6)
lines(x=cstar.seq[c(1, 2, 7, 12, 18, 26, 27, 28, 29)], 
      y= apply(cp_indrct11.runze[,c(1, 2, 7, 12, 18, 26, 27, 28, 29)], 2, mean), 
      type = "b", col=6, pch =5, lwd = 2, cex=0.6)
lines(x=cstar.seq[c(1, 2, 7, 12, 18, 26, 27, 28, 29)], 
      y= apply(cp_indrct11.dave.cY[,c(1, 2, 7, 12, 18, 26, 27, 28, 29)], 2, mean), 
      type = "b", col=4, pch =5, lwd = 2, cex=0.6)
lines(x=cstar.seq[c(1, 2, 7, 12, 18, 26, 27, 28, 29)], 
      y= apply(cp_indrct.11.Huang.boot[,c(1, 2, 7, 12, 18, 26, 27, 28, 29)], 2, mean), 
      type = "b", col=7, pch =5, lwd = 2, cex=0.6)
abline(h=0.95, lty=3)
abline(h=0, lty=3)
axis(side = 2, at = c(0,  0.5, 0.95))
legend("right", legend = c("Proposed", 
                           "Oracle", "Runze", "Dave", "Huang"),
       pch=c(18, 5, 5,5,5,5), col=c(1, 2, 6, 4, 7), 
       lty=c(1,2,1, 1, 1), lwd = c(3,2,2,2,2),
       cex = 0.8)




plot(x=cstar.seq, y = apply(reject_indrct.10.scad, 2, mean), type = "b", col=1, pch=18, lwd=3,
     xlab = expression(c[D]),
     ylab = "Empirical Power", ylim = c(0,1), yaxt="n",
     main = expression(delta(T[1]==0) == 4.2) , cex=0.6)
lines(x=cstar.seq, y = apply(reject_indrct.10.oracle, 2, mean), type = "b", col=2, pch=5, lwd=2, lty=2, cex = 0.6)
lines(x=cstar.seq, y= apply(reject_indrct.runze, 2, mean), type = "b", col=6, pch =5, lwd = 2, cex=0.6)
lines(x=cstar.seq, y= apply(reject_indrct.dave.cY, 2, mean), type = "b", col=4, pch =5, lwd = 2, cex=0.6)
lines(x=cstar.seq, y= apply(reject_indrct.10.Huang.boot, 2, mean), type = "b", col=7, pch =5, lwd = 2, cex=0.6)
abline(h=0.95, lty=3)
abline(h=0, lty=3)
axis(side = 2, at = c(0,  0.5, 0.95))
legend("right", legend = c("Proposed", 
                           "Oracle", "Runze", "Dave", "Huang"),
       pch=c(18, 5, 5,5,5,5), col=c(1, 2, 6, 4, 7), 
       lty=c(1,2,1, 1, 1), lwd = c(3,2,2,2,2),
       cex = 0.8)

plot(x=cstar.seq[c(1, 2, 7, 12, 18, 26, 27, 28, 29)], 
     y = apply(reject_indrct.10.scad[,c(1, 2, 7, 12, 18, 26, 27, 28, 29)], 2, mean), 
     type = "b", col=1, pch=18, lwd=3,
     xlab = expression(c[D]),
     ylab = "Empirical Power", ylim = c(0,1), yaxt="n",
     main = expression(delta(T[1]==0) == 4.2) , cex=0.6)
lines(x=cstar.seq[c(1, 2, 7, 12, 18, 26, 27, 28, 29)], 
      y = apply(reject_indrct.10.oracle[,c(1, 2, 7, 12, 18, 26, 27, 28, 29)], 2, mean), 
      type = "b", col=2, pch=5, lwd=2, lty=2, cex = 0.6)
lines(x=cstar.seq[c(1, 2, 7, 12, 18, 26, 27, 28, 29)], 
      y= apply(reject_indrct.runze[,c(1, 2, 7, 12, 18, 26, 27, 28, 29)], 2, mean), 
      type = "b", col=6, pch =5, lwd = 2, cex=0.6)
lines(x=cstar.seq[c(1, 2, 7, 12, 18, 26, 27, 28, 29)], 
      y= apply(reject_indrct.dave.cY[,c(1, 2, 7, 12, 18, 26, 27, 28, 29)], 2, mean), 
      type = "b", col=4, pch =5, lwd = 2, cex=0.6)
lines(x=cstar.seq[c(1, 2, 7, 12, 18, 26, 27, 28, 29)], 
      y= apply(reject_indrct.10.Huang.boot[,c(1, 2, 7, 12, 18, 26, 27, 28, 29)], 2, mean), 
      type = "b", col=7, pch =5, lwd = 2, cex=0.6)
abline(h=0.95, lty=3)
abline(h=0, lty=3)
axis(side = 2, at = c(0,  0.5, 0.95))
legend("right", legend = c("Proposed", 
                           "Oracle", "Runze", "Dave", "Huang"),
       pch=c(18, 5, 5,5,5,5), col=c(1, 2, 6, 4, 7), 
       lty=c(1,2,1, 1, 1), lwd = c(3,2,2,2,2),
       cex = 0.8)






plot(x=cstar.seq, y = apply(reject_indrct.11.scad, 2, mean), type = "b", col=1, pch=18, lwd=3,
     xlab = expression(c[D]),
     ylab = "Empirical Power", ylim = c(0,1), yaxt="n",
     main = expression( delta(T[1]==1) == 4.2 (1 + c[D]) ) , cex=0.6)
lines(x=cstar.seq, y = apply(reject_indrct.11.oracle, 2, mean), type = "b", col=2, pch=5, lwd=2, lty=2, cex = 0.6)
lines(x=cstar.seq, y= apply(reject_indrct.runze, 2, mean), type = "b", col=6, pch =5, lwd = 2, cex=0.6)
lines(x=cstar.seq, y= apply(reject_indrct.dave.cY, 2, mean), type = "b", col=4, pch =5, lwd = 2, cex=0.6)
#lines(x=cstar.seq, y= apply(reject_indrct.dave, 2, mean), type = "b", col="gray", pch =5, lwd = 2, cex=0.6)
lines(x=cstar.seq, y= apply(reject_indrct.11.Huang.boot, 2, mean), type = "b", col=7, pch =5, lwd = 2, cex=0.6)
abline(h=0.05, lty=3)
abline(h=0.95, lty=3)
axis(side = 2, at = c(0.05,  0.5, 0.95))
legend("right", legend = c("Proposed", 
                          "Oracle", "Runze", "Dave", "Huang"),
       pch=c(18, 5, 5,5,5,5), col=c(1, 2, 6, 4, 7), 
       lty=c(1,2,1, 1, 1), lwd = c(3,2,2,2,2),
       cex = 0.8)



plot(x=cstar.seq[3:11], 
     y = apply(reject_indrct.11.scad[,3:11], 2, mean), 
     type = "b", col=1, pch=18, lwd=3,
     xlab = expression(c[D]),
     ylab = "Empirical Power", ylim = c(0,1), yaxt="n",
     main = expression( delta(T[1]==1) == 4.2 (1 + c[D]) ) , cex=0.6)
lines(x=cstar.seq[3:11], y = apply(reject_indrct.11.oracle[,3:11], 2, mean), type = "b", col=2, pch=5, lwd=2, lty=2, cex = 0.6)
lines(x=cstar.seq[3:11], y= apply(reject_indrct.runze[,3:11], 2, mean), type = "b", col=6, pch =5, lwd = 2, cex=0.6)
lines(x=cstar.seq[3:11], y= apply(reject_indrct.dave.cY[,3:11], 2, mean), type = "b", col=4, pch =5, lwd = 2, cex=0.6)
lines(x=cstar.seq[3:11], y= apply(reject_indrct.11.Huang.boot[,3:11], 2, mean), type = "b", col=7, pch =5, lwd = 2, cex=0.6)
abline(h=0.05, lty=3)
abline(h=0.95, lty=3)
axis(side = 2, at = c(0.05,  0.5, 0.95))
legend("left", legend = c("Proposed", 
                          "Oracle", "Runze", "Dave", "Huang"),
       pch=c(18, 5, 5,5,5,5), col=c(1, 2, 6, 4, 7), 
       lty=c(1,2,1, 1, 1), lwd = c(3,2,2,2,2),
       cex = 0.8)


boxplot(cbind(SE_indrct.10.Huang, SE_indrct.10.Huang.boot))





### lasso ###
# compare bias from lasso and our method
# first calculate bias of our method
bias_indrct.10.scad=matrix(0, nrow = run, ncol = length(cstar.seq))
bias_indrct.11.scad=matrix(0, nrow = run, ncol = length(cstar.seq))
for (loop.c in 1:length(cstar.seq)) {
  
  c.D = cstar.seq[loop.c]
  Theta.D=c.D*matrix( c( c(1, 0.9, 0.8, 0.7, 0.6), 
                         c(1, 0.9, 0.8, 0.7, 0.6),
                         rep(0, p - s.D) ), 
                      ncol = q )
  print(paste0("NIE0=", (t(B)%*%theta.G) ))
  print(paste0("NIE1=",(t(B)%*%theta.G + t(Theta.D)%*%B) ))
  
  for (loop.run in 1:run) {
    
    bias_indrct.10.scad[loop.run, loop.c] = est_indrct.10.scad[loop.run, loop.c] - t(B)%*%theta.G
    bias_indrct.11.scad[loop.run, loop.c] = est_indrct.11.scad[loop.run, loop.c] - (t(B)%*%theta.G + t(Theta.D)%*%B)
    
  }
  
}


##ggplot
value=c()

# NIE(T=0)
for (i in 1:length(cstar.seq)) {
  value=c(value, (bias_indrct.10.lasso[,i]), (bias_indrct.10.scad[,i]))
}

cD_nor=rep(cstar.seq,each=2*run)
cD_nor=factor(cD_nor, 
              levels = as.character(cstar.seq) )

method=rep( rep(c("Lasso","SCAD"), each=run), times=length(cstar.seq) )

total=data.frame(value, cD_nor, method)


# NIE(T=0)
ggplot(total, aes(x=cD_nor, y=(value), fill=method)) +
  geom_boxplot() + scale_fill_manual(values=c("turquoise", "tomato"))  +
  labs(x=expression(c[D]),y=NULL,title = NULL)+ 
  theme(legend.title = element_blank()) + 
  ggtitle(expression( Bias~of~delta(T[1]==0)==4.2  ))+
  geom_hline(yintercept=0,  color = "red")






##ggplot
value=c()

# NIE(T=1)
for (i in 1:length(cstar.seq)) {
  value=c(value, (bias_indrct.11.lasso[,i]), (bias_indrct.11.scad[,i]))
}

cD_nor=rep(cstar.seq,each=2*run)
cD_nor=factor(cD_nor, 
              levels = as.character(cstar.seq) )

method=rep( rep(c("Lasso","SCAD"), each=run), times=length(cstar.seq) )

total=data.frame(value, cD_nor, method)

# NIE(T=1)
ggplot(total, aes(x=cD_nor, y=(value), fill=method)) +
  geom_boxplot() + scale_fill_manual(values=c("turquoise", "tomato"))  +
  labs(x=expression(c[D]),y=NULL,title = NULL)+ 
  theme(legend.title = element_blank()) + 
  ggtitle(expression( Bias~of~delta(T[1]==1)==4.2~(1+c[D])  ))+
  geom_hline(yintercept=0,  color = "red")




### compare Runze and Dave estimation bias
value=c()

for (i in 1:length(cstar.seq)) {
  value=c(value, (est_indrct.runze[,i]), (est_indrct.dave.cY[,i]) )
}

cD_nor=rep(cstar.seq,each=( nrow(est_indrct.runze) + nrow(est_indrct.dave.cY) ) )
cD_nor=factor(cD_nor, 
              levels = as.character(cstar.seq) )

method=rep( c( rep("Runze", nrow(est_indrct.runze)), 
               rep("Dave", nrow(est_indrct.dave.cY)) ) , 
            times=length(cstar.seq) )

total=data.frame(value, cD_nor, method)



ggplot(total, aes(x=cD_nor, y=(value), fill=method)) +
  geom_boxplot() + scale_fill_manual(values=c("cyan", "pink", "orange"))  +
  labs(x=expression(c[D]),y=NULL,title = NULL)+ 
  theme(legend.title = element_blank()) + 
  ggtitle(expression(Point~Estimates~of~NIE~by~Runze~and~Dave))




### Huang 2016 bias
library(ggplot2)
value=c()

# NIE(T=0)
for (i in 1:length(cstar.seq)) {
  value=c(value, est_indrct.10.Huang[,i] )
}

cD_nor=rep(cstar.seq,each=(run))
cD_nor=factor(cD_nor, 
              levels = as.character(cstar.seq) )


total=data.frame(value, cD_nor)

# NIE(T=0)
truevalue = data.frame(location=c(1:length(cstar.seq)), true.value=rep(4.2, length(cstar.seq)))

# NIE(T=0)
ggplot(total, aes(x=cD_nor, y=(value))) +
  geom_boxplot(fill="cyan")   +
  labs(x=expression(c[D]),y=NULL,title = NULL)+ 
  theme(legend.title = element_blank()) + 
  ggtitle(expression(Point~Estimates~of~delta(T[1]==0) == 4.2 ~by~Huang) ) +
  geom_point(data=truevalue, aes(x=location,y=true.value),colour="red", size=2, shape = 24) # shape=24 is NIE0








### Huang 2016 bias
library(ggplot2)
value=c()

# NIE(T=1)
for (i in 1:length(cstar.seq)) {
  value=c(value, est_indrct.11.Huang[,i] )
}

cD_nor=rep(cstar.seq,each=(run))
cD_nor=factor(cD_nor, 
              levels = as.character(cstar.seq) )


total=data.frame(value, cD_nor)

# NIE(T=1)
truevalue = data.frame(location=c(1:length(cstar.seq)), true.value=4.2*(1+cstar.seq))

# NIE(T=1)
ggplot(total, aes(x=cD_nor, y=(value))) +
  geom_boxplot(fill="cyan")   +
  labs(x=expression(c[D]),y=NULL,title = NULL)+ 
  theme(legend.title = element_blank()) + 
  ggtitle(expression(Point~Estimates~of~delta(T[1]==1) == 4.2 (1 + c[D]) ~by~Huang) ) +
  geom_point(data=truevalue, aes(x=location,y=true.value),colour="red", size=2, shape = 17) # shape=17 is NIE1


boxplot(SE_indrct.10.Huang, names = cstar.seq)
plot(x=cstar.seq, y = apply(est_indrct.10.Huang, 2, sd), type = "b", pch=19)
lines(x = cstar.seq, y = apply(SE_indrct.10.Huang, 2, mean), type="b", pch=1, col=2)
lines(x = cstar.seq, y = apply(SE_indrct.10.Huang.boot, 2, mean), type="b", pch=1, col=4)

boxplot(SE_indrct.11.Huang, names = cstar.seq)
plot(x=cstar.seq, y = apply(est_indrct.11.Huang, 2, sd), type = "b", pch=19)
lines(x = cstar.seq, y = apply(SE_indrct.11.Huang, 2, mean), type="b", pch=1, col=2)
lines(x = cstar.seq, y = apply(SE_indrct.11.Huang.boot, 2, mean), type="b", pch=1, col=4)


############################### draw the graph showing cp and powers from all methods #################
library(ggplot2)
library("ggpubr")


# NDE
value = c(apply(cp_drct.10.scad[,c(1, 2, 7, 12, 13, 14, 24, 29, 30) ], 2, mean), 
          apply(cp_drct.10.oracle[,c(1, 2, 7, 12, 13, 14, 24, 29, 30) ], 2, mean), 
          apply(cp_drct10.runze[,c(1, 2, 12, 13, 14, 24, 29, 30) ], 2, mean), 
          apply(cp_drct10.dave.cY[,c(1, 2, 7, 12, 13, 14, 24, 29, 30) ], 2, mean),
          apply(cp_drct.11.scad[,c(1, 2, 7, 12, 13, 14, 24, 29, 30) ], 2, mean),
          apply(cp_drct.11.oracle[,c(1, 2, 7, 12, 13, 14, 24, 29, 30) ], 2, mean), 
          apply(cp_drct11.runze[,c(1, 2, 7, 12, 13, 14, 24, 29, 30) ], 2, mean), 
          apply(cp_drct11.dave.cY[,c(1, 2, 7, 12, 13, 14, 24, 29, 30) ], 2, mean),
          apply(reject_drct.10.scad[,c(14:29) ], 2, mean),
          apply(reject_drct.10.oracle[,c(14:29) ], 2, mean),
          apply(reject_drct.runze[,c(14:29)], 2, mean),
          apply(reject_drct.dave.cY[,c(14:29) ], 2, mean) )
c_D = c(
  cstar.seq[c(1, 2, 7, 12, 13, 14, 24, 29, 30)],
  cstar.seq[c(1, 2, 7, 12, 13, 14, 24, 29, 30)],
  cstar.seq[c(1, 2, 12, 13, 14, 24, 29, 30)],
  cstar.seq[c(1, 2, 7, 12, 13, 14, 24, 29, 30)],
  cstar.seq[c(1, 2, 7, 12, 13, 14, 24, 29, 30)],
  cstar.seq[c(1, 2, 7, 12, 13, 14, 24, 29, 30)],
  cstar.seq[c(1, 2, 7, 12, 13, 14, 24, 29, 30)],
  cstar.seq[c(1, 2, 7, 12, 13, 14, 24, 29, 30)],
  cstar.seq[c(14:29)],
  cstar.seq[c(14:29)],
  cstar.seq[c(14:29)],
  cstar.seq[c(14:29)]
)


Method = c(    rep("M1", times=length(c(1, 2, 7, 12, 13, 14, 24, 29, 30))),
               rep("M0", times=length(c(1, 2, 7, 12, 13, 14, 24, 29, 30))),
               rep("M2", times=length(c(1, 2, 12, 13, 14, 24, 29, 30))),
               rep("M3", times=length(c(1, 2, 7, 12, 13, 14, 24, 29, 30))),
               rep("M1", times=length(c(1, 2, 7, 12, 13, 14, 24, 29, 30))),
               rep("M0", times=length(c(1, 2, 7, 12, 13, 14, 24, 29, 30))),
               rep("M2", times=length(c(1, 2, 7, 12, 13, 14, 24, 29, 30))),
               rep("M3", times=length(c(1, 2, 7, 12, 13, 14, 24, 29, 30))),
               rep("M1", times=length(c(14:29))),
               rep("M0", times=length(c(14:29))),
               rep("M2", times=length(c(14:29))),
               rep("M3", times=length(c(14:29)))
           )
             

Measure = c( rep("C0", times=length( c(
                           cstar.seq[c(1, 2, 7, 12, 13, 14, 24, 29, 30)],
                           cstar.seq[c(1, 2, 7, 12, 13, 14, 24, 29, 30)],
                           cstar.seq[c(1, 2, 12, 13, 14, 24, 29, 30)],
                           cstar.seq[c(1, 2, 7, 12, 13, 14, 24, 29, 30)] )
                           ) ),
             rep("C1", times=length( c(
                           cstar.seq[c(1, 2, 7, 12, 13, 14, 24, 29, 30)],
                           cstar.seq[c(1, 2, 7, 12, 13, 14, 24, 29, 30)],
                           cstar.seq[c(1, 2, 7, 12, 13, 14, 24, 29, 30)],
                           cstar.seq[c(1, 2, 7, 12, 13, 14, 24, 29, 30)] )   
                           ) ),
             rep("P0", times=length( c(
                           cstar.seq[c(14:29)],
                           cstar.seq[c(14:29)],
                           cstar.seq[c(14:29)],
                           cstar.seq[c(14:29)] )
               ) )
             )


Effect = rep("NDE", times = length(value) )

NDE.cppwr = data.frame(value=value, c_D = c_D, Method = Method, Measure=Measure, Effect=Effect)


library(ggh4x) # this package make the line not pass through the points by using geom_pointpath() rather than using: geom_line(aes(color = Method)) + geom_point(aes(color = Method, shape = Method))

ggplot(NDE.cppwr, aes(x = c_D, y = value))+ 
  geom_pointpath(aes(color = Method, shape = Method), size=2, linewidth=0.7, stroke=1.4)+
  facet_grid(cols = vars(Measure), scales = "free_x") +
  scale_color_manual(values = c("red", "black", "magenta", "blue")) +
  scale_shape_manual(values = c(1, 19, 2, 3)) + 
  geom_hline(yintercept=0.95,color = "black",linetype=2) +
  geom_hline(yintercept=0,color = "black",linetype=2) +
  geom_hline(yintercept=1,color = "black",linetype=2) +
  geom_hline(yintercept=0.05,color = "black",linetype=2) +
  scale_y_continuous(breaks=c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1)) + 
  labs(x=expression(c[D]),y=" ") 






# NIE 
value = c(apply(cp_indrct.10.scad[,c(1, 2, 7, 12, 13, 14, 24, 29, 30)], 2, mean),
          apply(cp_indrct.10.oracle[,c(1, 2, 7, 12, 13, 14, 24, 29, 30)], 2, mean),
          apply(cp_indrct10.runze[,c(1, 2, 7, 12, 13, 14, 24, 29, 30)], 2, mean),
          apply(cp_indrct10.dave.cY[,c(1, 2, 7, 12, 13, 14, 24, 29, 30)], 2, mean),
          apply(cp_indrct.10.Huang.boot[,c(1, 2, 7, 12, 13, 14, 24, 29, 30)], 2, mean),
          apply(cp_indrct.11.scad[,c(1, 2, 7, 12, 13, 14, 24, 29, 30)], 2, mean),
          apply(cp_indrct.11.oracle[,c(1, 2, 7, 12, 13, 14, 24, 29, 30)], 2, mean),
          apply(cp_indrct11.runze[,c(1, 2, 12, 13, 14, 24, 29, 30)], 2, mean),
          apply(cp_indrct11.dave.cY[,c(1, 2, 7, 12, 13, 14, 24, 29, 30)], 2, mean),
          apply(cp_indrct.11.Huang.boot[,c(1, 2, 7, 12, 13, 14, 24, 29, 30)], 2, mean),
          apply(reject_indrct.11.scad[,3:11], 2, mean),
          apply(reject_indrct.11.oracle[,3:11], 2, mean),
          apply(reject_indrct.runze[,3:11], 2, mean),
          apply(reject_indrct.dave.cY[,3:11], 2, mean),
          apply(reject_indrct.11.Huang.boot[,3:11], 2, mean))



c_D = c(cstar.seq[c(1, 2, 7, 12, 13, 14, 24, 29, 30)],
        cstar.seq[c(1, 2, 7, 12, 13, 14, 24, 29, 30)],
        cstar.seq[c(1, 2, 7, 12, 13, 14, 24, 29, 30)],
        cstar.seq[c(1, 2, 7, 12, 13, 14, 24, 29, 30)],
        cstar.seq[c(1, 2, 7, 12, 13, 14, 24, 29, 30)],
        cstar.seq[c(1, 2, 7, 12, 13, 14, 24, 29, 30)],
        cstar.seq[c(1, 2, 7, 12, 13, 14, 24, 29, 30)],
        cstar.seq[c(1, 2, 12, 13, 14, 24, 29, 30)],
        cstar.seq[c(1, 2, 7, 12, 13, 14, 24, 29, 30)],
        cstar.seq[c(1, 2, 7, 12, 13, 14, 24, 29, 30)],
        cstar.seq[c(3:11)],
        cstar.seq[c(3:11)],
        cstar.seq[c(3:11)],
        cstar.seq[c(3:11)],
        cstar.seq[c(3:11)])

Method = c(    rep("M1", times=length(c(1, 2, 7, 12, 13, 14, 24, 29, 30) ) ),
               rep("M0", times=length(c(1, 2, 7, 12, 13, 14, 24, 29, 30) ) ),
               rep("M2", times=length(c(1, 2, 7, 12, 13, 14, 24, 29, 30) ) ),
               rep("M3", times=length(c(1, 2, 7, 12, 13, 14, 24, 29, 30) ) ),
               rep("M4", times=length(c(1, 2, 7, 12, 13, 14, 24, 29, 30) ) ),
               rep("M1", times=length(c(1, 2, 7, 12, 13, 14, 24, 29, 30) ) ),
               rep("M0", times=length(c(1, 2, 7, 12, 13, 14, 24, 29, 30) ) ),
               rep("M2", times=length(c(1, 2, 12, 13, 14, 24, 29, 30) ) ),
               rep("M3", times=length(c(1, 2, 7, 12, 13, 14, 24, 29, 30) ) ),
               rep("M4", times=length(c(1, 2, 7, 12, 13, 14, 24, 29, 30) ) ),
               rep("M1", times=length(c(3:11) ) ),
               rep("M0", times=length(c(3:11) ) ),
               rep("M2", times=length(c(3:11) ) ),
               rep("M3", times=length(c(3:11) ) ),
               rep("M4", times=length(c(3:11) ) )
)

Measure = c( rep("C0", times=length( c(
  cstar.seq[c(1, 2, 7, 12, 13, 14, 24, 29, 30)],
  cstar.seq[c(1, 2, 7, 12, 13, 14, 24, 29, 30)],
  cstar.seq[c(1, 2, 7, 12, 13, 14, 24, 29, 30)],
  cstar.seq[c(1, 2, 7, 12, 13, 14, 24, 29, 30)],
  cstar.seq[c(1, 2, 7, 12, 13, 14, 24, 29, 30)])
) ),
rep("C1", times=length( c(
  cstar.seq[c(1, 2, 7, 12, 13, 14, 24, 29, 30)],
  cstar.seq[c(1, 2, 7, 12, 13, 14, 24, 29, 30)],
  cstar.seq[c(1, 2, 12, 13, 14, 24, 29, 30)],
  cstar.seq[c(1, 2, 7, 12, 13, 14, 24, 29, 30)],
  cstar.seq[c(1, 2, 7, 12, 13, 14, 24, 29, 30)])   
) ),
rep("P1", times=length( c(
  cstar.seq[c(3:11)],
  cstar.seq[c(3:11)],
  cstar.seq[c(3:11)],
  cstar.seq[c(3:11)],
  cstar.seq[c(3:11)])
) )
)

Effect = rep("NIE", times = length(value) )


NIE.cppwr = data.frame(value=value, c_D = c_D, Method = Method, Measure=Measure, Effect=Effect)


library(ggh4x)




ggplot(NIE.cppwr, aes(x = c_D, y = value))  + 
  facet_grid(cols = vars(Measure), scales = "free_x") +
  geom_pointpath(aes(color = Method, shape = Method)) + 
  scale_color_manual(values = c("red", "black",  "magenta", "blue", "orange")) +
  scale_shape_manual(values = c(1, 19, 2, 3, 4)) + 
  geom_hline(yintercept=0.95,color = "black",linetype=2) +
  geom_hline(yintercept=0,color = "black",linetype=2) +
  geom_hline(yintercept=1,color = "black",linetype=2) +
  geom_hline(yintercept=0.05,color = "black",linetype=2) +
  scale_y_continuous(breaks=c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1)) + 
  labs(x=expression(c[D]),y=" ") +
  facetted_pos_scales(x = list(scale_x_continuous(limits = c(-2, 2)), # this sentence from ggh4x allows us to define different limits of x axis in facet
                               scale_x_continuous(limits = c(-2, 2)),
                               scale_x_continuous(limits = c(-1.08, -0.92))) )




cppwr = rbind(NIE.cppwr, NDE.cppwr)

# cppwr[which(cppwr$Measure == "P0"), "Measure"] = "P"
# cppwr[which(cppwr$Measure == "P1"), "Measure"] = "P"




ggplot(cppwr, aes(x = c_D, y = value))  + 
  facet_grid(cols = vars(Measure), rows = vars(Effect),  scales = "free_x") +
  geom_pointpath(aes(color = Method, shape = Method)) + 
  scale_color_manual(values = c("red", "black", "magenta", "blue", "orange")) +
  scale_shape_manual(values = c(1, 19, 2, 3, 4)) + 
  geom_hline(yintercept=0.95,color = "black",linetype=2) +
  geom_hline(yintercept=0,color = "black",linetype=2) +
  geom_hline(yintercept=1,color = "black",linetype=2) +
  scale_y_continuous(breaks=c(0, 0.25, 0.5, 0.75, 0.95, 1)) + 
  labs(x=expression(c[D]),y=" ") +
  facetted_pos_scales(x = list(
    scale_x_continuous(limits = c(-2, 2)), 
    scale_x_continuous(limits = c(-2, 2)),
    scale_x_continuous(limits = c(-1.08, 0.2))) # the x-axis scale of the third column must be parallel, I haven't find a way to individually define the upper and lower x-axis
                      )




# I have to define a new column, consist of 6 categories "NDE C0", "NDE C1" "NDE P0" "NIE C0", "NIE C1", "NIE P1"
# Then I can define the x-axis of these plots individually, also, ncol = 3 in the facet_grid()
cppwr$Type = rep(NA, times = nrow(cppwr))
for (i in 1:nrow(cppwr)) {
  cppwr[i, "Type"] = paste(cppwr[i,"Effect"], cppwr[i,"Measure"])
}

cppwr$Type=factor(cppwr$Type, 
              levels = c("NIE C0", "NIE C1", "NIE P1",
                         "NDE C0", "NDE C1", "NDE P0") )


ggplot(cppwr, aes(x = c_D, y = value))  + 
  facet_wrap(~Type, ncol = 3,  scales = "free_x") +
  geom_pointpath(aes(color = Method, shape = Method), size=2, linewidth=0.8) + 
  scale_color_manual(values = c("red", "black", "magenta", "blue", "orange")) +
  scale_shape_manual(values = c(1, 19, 2, 3, 4)) + 
  geom_hline(yintercept=0.95,color = "black",linetype=2) +
  geom_hline(yintercept=0.05,color = "black",linetype=2) +
  geom_hline(yintercept=1,color = "black",linetype=2) +
  scale_y_continuous(breaks=c(0, 0.05,0.25, 0.5, 0.75, 0.95, 1)) + 
  labs(x=expression(c[D]),y=" ") +
  facetted_pos_scales(x = list(
    scale_x_continuous(limits = c(-2, 2)), 
    scale_x_continuous(limits = c(-2, 2)),
    scale_x_continuous(limits = c(-1.08, -0.92)), 
    scale_x_continuous(limits = c(-2, 2)),
    scale_x_continuous(limits = c(-2, 2)), 
    scale_x_continuous(limits = c(-1, 1))) # the x-axis scale of the third column must be parallel, I haven't find a way to individually define the upper and lower x-axis
  )

ggplot(cppwr, aes(x = c_D, y = value))  + 
  facet_wrap(~Type, ncol = 3,  scales = "free_x") +
  geom_pointpath(aes(color = Method, shape = Method), size=2, linewidth=0.7, stroke=1.4) + # stroke change the point border thickness
  scale_color_manual(values = c("red", "black", "magenta", "blue", "orange")) +
  scale_shape_manual(values = c(1, 19, 2, 3, 4)) + 
  geom_hline(yintercept=0.95,color = "black",linetype=2) +
  geom_hline(yintercept=0.05,color = "black",linetype=2) +
  geom_hline(yintercept=1,color = "black",linetype=2) +
  scale_y_continuous(breaks=c(0, 0.05,0.25, 0.5, 0.75, 0.95, 1)) + 
  labs(x=expression(c[D]),y=" ") +
  theme(strip.text = element_text(size = 15, color = "black"), # this changes the size of labellers (titles) of each plot
        legend.title = element_text(size=16), # the size of legend title in the legend
        legend.key.size = unit(1.1, 'cm'), # the size of red circle, black dots ... in legend
        legend.text = element_text(size=10), # the size of texts "M0", "M1" ... in the legend
        axis.text.y=element_text(size=11), # the size of tick mark labels on y-axis
        axis.title.x = element_text(size = 15),
        axis.text.x=element_text(size=12) )





################### lasso bias line plot ##########################################
bias_drct.10.scad=matrix(0, nrow = run, ncol = length(cstar.seq))
bias_drct.11.scad=matrix(0, nrow = run, ncol = length(cstar.seq))
for (loop.c in 1:length(cstar.seq)) {
  
  c.D = cstar.seq[loop.c]
  Theta.D=c.D*matrix( c( c(1, 0.9, 0.8, 0.7, 0.6), 
                         c(1, 0.9, 0.8, 0.7, 0.6),
                         rep(0, p - s.D) ), 
                      ncol = q )
  NDE0 = as.numeric( theta.T + t(beta0)%*%Theta.D )
  NDE1 = as.numeric( theta.T + t(beta0 + B)%*%Theta.D  )
  
  print(paste0( "c_D=", c.D, ", NDE0=", NDE0, ", NDE1=", NDE1))
  
  for (loop.run in 1:run) {
    
    bias_drct.10.scad[loop.run, loop.c] = est_drct.10.scad[loop.run, loop.c] - NDE0
    bias_drct.11.scad[loop.run, loop.c] = est_drct.11.scad[loop.run, loop.c] - NDE1
    
  }
  
}


bias_indrct.10.scad=matrix(0, nrow = run, ncol = length(cstar.seq))
bias_indrct.11.scad=matrix(0, nrow = run, ncol = length(cstar.seq))
for (loop.c in 1:length(cstar.seq)) {
  
  c.D = cstar.seq[loop.c]
  Theta.D=c.D*matrix( c( c(1, 0.9, 0.8, 0.7, 0.6), 
                         c(1, 0.9, 0.8, 0.7, 0.6),
                         rep(0, p - s.D) ), 
                      ncol = q )
  print(paste0( "c_D = ", c.D,  ", NIE0=", (t(B)%*%theta.G), ", NIE1=", (t(B)%*%theta.G + t(Theta.D)%*%B) ))
  
  for (loop.run in 1:run) {
    
    bias_indrct.10.scad[loop.run, loop.c] = est_indrct.10.scad[loop.run, loop.c] - t(B)%*%theta.G
    bias_indrct.11.scad[loop.run, loop.c] = est_indrct.11.scad[loop.run, loop.c] - (t(B)%*%theta.G + t(Theta.D)%*%B)
    
  }
  
}

library(ggplot2)
library("ggpubr")
library(ggh4x)


# NDE
value = c(
          apply(bias_drct.10.scad[,c(1, 2, 7, 12, 13, 14, 24, 29, 30) ], 2, mean), 
          apply(bias_drct.10.lasso[,c(1, 2, 7, 12, 13, 14, 24, 29, 30) ], 2, mean), 
          apply(bias_drct.11.scad[,c(1, 2, 7, 12, 13, 14, 24, 29, 30) ], 2, mean), 
          apply(bias_drct.11.lasso[,c(1, 2, 7, 12, 13, 14, 24, 29, 30) ], 2, mean)
          )
c_D = c(
  cstar.seq[c(1, 2, 7, 12, 13, 14, 24, 29, 30)],
  cstar.seq[c(1, 2, 7, 12, 13, 14, 24, 29, 30)],
  cstar.seq[c(1, 2, 7, 12, 13, 14, 24, 29, 30)],
  cstar.seq[c(1, 2, 7, 12, 13, 14, 24, 29, 30)]
)


Method = c(    rep("M1", times=length(c(1, 2, 7, 12, 13, 14, 24, 29, 30))),
               rep("M5", times=length(c(1, 2, 7, 12, 13, 14, 24, 29, 30))),
               rep("M1", times=length(c(1, 2, 7, 12, 13, 14, 24, 29, 30))),
               rep("M5", times=length(c(1, 2, 7, 12, 13, 14, 24, 29, 30)))
)

Type = c ( 
  rep("NDE T=0", times = length( c( 
    c(1, 2, 7, 12, 13, 14, 24, 29, 30),
    c(1, 2, 7, 12, 13, 14, 24, 29, 30) 
                                           ) ) ), 
  rep("NDE T=1", times = length( c( 
    c(1, 2, 7, 12, 13, 14, 24, 29, 30),
    c(1, 2, 7, 12, 13, 14, 24, 29, 30) 
                                           ) ) )
)

NDE.bias = data.frame(value=value, c_D = c_D, Method = Method, Type = Type)


# NIE
value = c(
  apply(bias_indrct.10.scad[,c(1, 2, 7, 12, 13, 14, 24, 29, 30) ], 2, mean), 
  apply(bias_indrct.10.lasso[,c(1, 2, 7, 12, 13, 14, 24, 29, 30) ], 2, mean), 
  apply(bias_indrct.11.scad[,c(1, 2, 7, 12, 13, 14, 24, 29, 30) ], 2, mean), 
  apply(bias_indrct.11.lasso[,c(1, 2, 7, 12, 13, 14, 24, 29, 30) ], 2, mean)
)
c_D = c(
  cstar.seq[c(1, 2, 7, 12, 13, 14, 24, 29, 30)],
  cstar.seq[c(1, 2, 7, 12, 13, 14, 24, 29, 30)],
  cstar.seq[c(1, 2, 7, 12, 13, 14, 24, 29, 30)],
  cstar.seq[c(1, 2, 7, 12, 13, 14, 24, 29, 30)]
)


Method = c(    rep("M1", times=length(c(1, 2, 7, 12, 13, 14, 24, 29, 30))),
               rep("M5", times=length(c(1, 2, 7, 12, 13, 14, 24, 29, 30))),
               rep("M1", times=length(c(1, 2, 7, 12, 13, 14, 24, 29, 30))),
               rep("M5", times=length(c(1, 2, 7, 12, 13, 14, 24, 29, 30)))
)

Type = c ( 
  rep("NIE T=0", times = length( c( 
    c(1, 2, 7, 12, 13, 14, 24, 29, 30),
    c(1, 2, 7, 12, 13, 14, 24, 29, 30) 
  ) ) ), 
  rep("NIE T=1", times = length( c( 
    c(1, 2, 7, 12, 13, 14, 24, 29, 30),
    c(1, 2, 7, 12, 13, 14, 24, 29, 30) 
  ) ) )
)

NIE.bias = data.frame(value=value, c_D = c_D, Method = Method, Type = Type)

bias = rbind(NDE.bias, NIE.bias)


ggplot(bias, aes(x = c_D, y = value))  + 
  facet_wrap(~Type, ncol = 4) +
  geom_pointpath(aes(color = Method, shape = Method), size=2, linewidth=0.8) + 
  scale_color_manual(values = c("black", "brown")) +
  scale_shape_manual(values = c(19, 5)) + 
  geom_hline(yintercept=0, color = "black",linetype=2) + 
  labs(x=expression(c[D]),y=" ") 




## the bias of lasso seems very weird ##

plot(x=cstar.seq, y=apply(bias_drct.11.lasso, 2, mean), 
     col=1, pch=19, type = "b", ylim = c(-1,1))
lines(x=cstar.seq, y=apply(bias_indrct.10.lasso, 2, mean), 
      col=2, pch=19, type = "b")

apply((bias_drct.11.lasso+bias_indrct.10.lasso), 2, mean)
apply((bias_drct.11.scad+bias_indrct.10.scad), 2, mean)

plot(x=cstar.seq, y=apply((bias_drct.11.lasso+bias_indrct.10.lasso), 2, mean), 
     col=2, pch=19, type = "b")
lines(x=cstar.seq, y=apply((bias_drct.11.scad+bias_indrct.10.scad), 2, mean), 
      col=1, pch=19, type = "b")




plot(x=cstar.seq, y=apply(bias_drct.10.lasso, 2, mean), 
     col=1, pch=19, type = "b", ylim = c(-1,1))
lines(x=cstar.seq, y=apply(bias_indrct.11.lasso, 2, mean), 
      col=2, pch=19, type = "b")

apply((bias_drct.10.lasso+bias_indrct.11.lasso), 2, mean)
apply((bias_drct.10.scad+bias_indrct.11.scad), 2, mean)

plot(x=cstar.seq, y=apply((bias_drct.10.lasso+bias_indrct.11.lasso), 2, mean), 
     col=2, pch=19, type = "b")
lines(x=cstar.seq, y=apply((bias_drct.10.scad+bias_indrct.11.scad), 2, mean), 
      col=1, pch=19, type = "b")
