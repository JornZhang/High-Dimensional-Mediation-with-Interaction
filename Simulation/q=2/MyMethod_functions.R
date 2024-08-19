
### Our method is our Method: oracal penalty + projection with interaction: OraProjInt



### The function to calculate NIE, NDE, and individual theta, SE's (sample version) ### 
## this function contains the covariate X ##
OraProjInt_SE_sample_X = function(Y, Trt, Mediator, Covariate, Trt.Medator,
                          theta0_hat, theta.T_hat, theta.X_hat, theta.G_hat, vtheta.D_hat){
  
  #' @param Y vector, n-dimensional outcome vector (with column name such as "Y").
  #' @param Trt matrix or dataframe, n by q treatment matrix without intercept (with column name such as "T1", "T2",.., "Tq").
  #' @param Covariate matrix or dataframe, n by r covariate matrix without intercept (with column name such as "X1", ...).
  #' @param Mediator matrix or dataframe, n by p mediator matrix (with column name, "G1",...).
  #' @param Trt.Medator matrix or dataframe, n by qp interaction matrix (with column name such as "T1.G1", "T1.G2",...).
  #' @param theta0_hat a vector with length 1 from SCAD estimation
  #' @param theta.T_hat a vector with length q from SCAD estimation
  #' @param theta.X_hat a vector with length r from SCAD estimation
  #' @param theta.G_hat a vector of length p from SCAD estimation
  #' @param vtheta.D_hat a vector of length pq, which is vec{Theta_D}, from SCAD estimation
  #' Note that Trt, Mediator, Covariate, Trt.Medator must have column names in words!
  #' theta0_hat, theta.T_hat, theta.X_hat, theta.G_hat, vtheta.D_hat are vectors with no names! Just pure vector
  #' If you use ncvreg packages, it will produce names of theta0_hat, theta.T_hat...
  #' To remove those names of theta, use the code before plugging them in this function such as: names(theta0.hat)=NULL; names(theta.T.hat)=NULL; names(theta.X.hat)=NULL; names(theta.G.hat)=NULL; names(vtheta.D.hat)=NULL
  #' @return outcome.theta_inference is a matrix of Estimate, SE, test-stat, p value of the SCAD outcome model inference for nonzero theta's
  #' @return NIE_inference, NDE_inference are the matrix of Estimate, SE, test-stat, p value of the projected estimators
  #' @return W_A = marix(one.n, T, X, G[,A_G], T.G[,A_D1,..., A_Dq])
  #' @return A_G: a vector of length \hat{s}_G, nonzero locations (number index in 1:p) in theta_G
  #' @return A_D: a list of vectors. Each vector is the nonzero location (number index in 1:p) in theta_{D_j} with length \hat{s}_{D_j}
  
  
  p = ncol(Mediator)
  n = nrow(Mediator)
  q = ncol(Trt)
  r = ncol(Covariate)
  
  A_G = which(theta.G_hat!=0) # A_G contains number location in 1:p with length \hat{s_G}
  A_G_cpg = colnames(Mediator)[A_G] # A_G contains word names of those support mediators with length \hat{s.G}
  
  A_D = list() # A_D has q components. Each components contains number location in 1:p with length \hat{s_{D_j}}
  A_D_cpg = list() # This contains the word names of those support mediators with length \hat{s_{D_j}}. These names are in the form of Trt*Medator name
  A_D_cpg_forunion=list() # This contains the word names of those support mediators with length \hat{s_{D_j}}. These names are still the mediator name rather than Trt*Medator name
  for (j in 1:q) {
    A_D[[j]] = which(vtheta.D_hat[ ((j-1)*p+1) : (j*p) ]!=0)
    tmp = Trt.Medator[, ((j-1)*p+1) : (j*p) ]
    A_D_cpg[[j]] = colnames(tmp)[A_D[[j]]]
    A_D_cpg_forunion[[j]] = colnames(Mediator)[A_D[[j]]] # if there is no support, length(A_D[[j]]) = 0, length(A_D_cpg_forunion[[j]])=0
  }
  rm(tmp)
  
  one.n = matrix(1, nrow=nrow(Mediator), ncol=1)
  
  W_A = cbind(one.n, Trt, Covariate, Mediator[,A_G])
  for (j in 1:q) {
    tmp = Trt.Medator[, ((j-1)*p+1) : (j*p) ]
    W_A = cbind(W_A, tmp[, A_D[[j]] ])
  }
  
  rm(tmp)
  
  tmpname = c()
  tmpname[1] = "one"
  tmpname[(1+1):(1+q)] = colnames(Trt)
  tmpname[(1+q+1):(1+q+r)] = colnames(Covariate)
  if(length(A_G) > 0){
    tmpname[(1+q+r+1):(1+q+r+length(A_G))] = A_G_cpg
    for (j in 1:q) {
      tmpname = c(tmpname, A_D_cpg[[j]]) # if length(A_D_cpg[[j]])=0, this step will not change the previous tmpname
    }
  } else {
    for (j in 1:q) {
      tmpname = c(tmpname, A_D_cpg[[j]])
    }
  }
  colnames(W_A) = tmpname
  rm(tmpname)
  
  
  Sigmahat = (t(W_A)%*%W_A) / nrow(Mediator)
  Sigmahat.inv = solve(Sigmahat)
  
  eps_Y = Y - one.n%*%theta0_hat - Trt%*%theta.T_hat - Covariate%*%theta.X_hat - Mediator%*%theta.G_hat - Trt.Medator%*%vtheta.D_hat
  sigma2_Y = as.numeric( t(eps_Y)%*%eps_Y/n ) # the estimate of sigma_Y^2
  rm(eps_Y)
  
  theta.G_hat_AG = theta.G_hat[A_G]
  Theta.D_hat_AD = list()
  for (j in 1:q) {
    tmp = vtheta.D_hat[ ((j-1)*p+1) : (j*p) ]
    Theta.D_hat_AD[[j]] = tmp[A_D[[j]]] # if A_D[[j]] is empty, length(Theta.D_hat_AD[[j]])=0
  }
  rm(tmp)
  
  # thetahat only store those nonzero theta's
  thetahat = c(theta0_hat, theta.T_hat, theta.X_hat, theta.G_hat_AG)
  for (j in 1:q) {
    thetahat = c(thetahat, Theta.D_hat_AD[[j]])
  }
  
  asymp.var_SCAD = sigma2_Y * Sigmahat.inv
  
  ### Outcome.theta_inference is the individual theta's inference from the SCAD step in the outcome model
  Outcome.theta_inference = cbind( thetahat, 
                                   sqrt(diag(asymp.var_SCAD))/sqrt(n), 
                                   sqrt(n)*thetahat /sqrt(diag(asymp.var_SCAD)),  
                                   2*pnorm(-abs( sqrt(n)*thetahat /sqrt(diag(asymp.var_SCAD)) )) )
  
  colnames(Outcome.theta_inference) = c("Estimate", "SE", "test-stat", "p value")
  
  
  ###
  H = cbind(one.n, Trt, Covariate)
  
  betahat = solve(t(H)%*%H)%*%t(H)%*%Mediator
  
  eps_G_hat = Mediator - H%*%betahat
  Sigma_G_hat = t(eps_G_hat)%*%eps_G_hat / n
  
  J = Sigmahat[1:(1+q+r),1:(1+q+r)]
  
  
  beta.theta.G = betahat %*% theta.G_hat # (1+q+r) by 1
  beta.Theta.D = list() # Each component is (1+q+r) by 1
  for (j in 1:q) {
    beta.Theta.D[[j]] = betahat %*% vtheta.D_hat[ ((j-1)*p+1) : (j*p) ] # if there is no active support, this vector could be zero
  }
  
  NIE_inference = matrix(NA, nrow = (q*2), ncol = 4)
  colnames(NIE_inference) = c("Estimate", "SE", "test-stat", "p value")
  rownames(NIE_inference) = rep("NIE", times = (q*2))
  
  NDE_inference = matrix(NA, nrow = (q*2), ncol = 4)
  colnames(NDE_inference) = c("Estimate", "SE", "test-stat", "p value")
  rownames(NDE_inference) = rep("NDE", times = (q*2))
  
  for (j in 1:q) {
    rownames(NIE_inference)[(2*j-1)] = paste("NIE(T", j, "= 0)")
    rownames(NIE_inference)[(2*j)] = paste("NIE(T", j, "= 1)")
    
    
    if(length(A_G)==0){
      NIE_inference[(2*j-1),"Estimate"] = 0
      NIE_inference[(2*j-1),"SE"] = 0
      NIE_inference[(2*j-1),"p value"] = 1
      
      if(length(A_D[[j]])==0){
        NIE_inference[(2*j),] = NIE_inference[(2*j-1),]
      } else {
        
        E_A_Dj = matrix(0, nrow = length(A_D[[j]]), ncol = nrow(Sigmahat))
        colnames(E_A_Dj) = colnames(Sigmahat)
        for (h in 1:nrow(E_A_Dj)) {
          E_A_Dj[h, A_D_cpg[[j]][h] ] = 1
        }
        
        theta.GDj_A = c(theta.G_hat_AG, Theta.D_hat_AD[[j]]) 
        Sigma_G_hat_GDj.A = Sigma_G_hat[c(A_G_cpg, A_D_cpg_forunion[[j]]), c(A_G_cpg, A_D_cpg_forunion[[j]])]
        
        var_indrct_j1 = as.numeric( sigma2_Y * (t(betahat[(j+1),A_D[[j]]])%*%E_A_Dj ) %*% (Sigmahat.inv) %*% ( t(E_A_Dj)%*%betahat[(j+1),A_D[[j]] ] ) ) +
          as.numeric((solve(J)[(j+1),(j+1)])* t(theta.GDj_A) %*% (Sigma_G_hat_GDj.A)%*% (theta.GDj_A))
        
        NIE_inference[(2*j),"Estimate"] = (beta.theta.G[(j+1),1] + beta.Theta.D[[j]][(j+1),1] )
        NIE_inference[(2*j),"SE"] = sqrt(var_indrct_j1)/sqrt(n)
        NIE_inference[(2*j),"test-stat"] = NIE_inference[(2*j),"Estimate"] / NIE_inference[(2*j),"SE"]
        NIE_inference[(2*j),"p value"] = 2*pnorm( -abs( NIE_inference[(2*j),"test-stat"] )  )
        
      }
      
      
    } else {
      
      var_indrct_j0 = as.numeric( (sigma2_Y) * t(betahat[(j+1),A_G]) %*% (Sigmahat.inv[A_G_cpg, A_G_cpg]) %*% (betahat[(j+1),A_G]) ) +
        as.numeric( (solve(J)[(j+1),(j+1)]) * t(theta.G_hat[A_G]) %*% (Sigma_G_hat[A_G, A_G])%*% (theta.G_hat[A_G])) 
      NIE_inference[(2*j-1),"Estimate"] = beta.theta.G[(j+1),1]
      NIE_inference[(2*j-1),"SE"] = sqrt(var_indrct_j0)/sqrt(n)
      NIE_inference[(2*j-1),"test-stat"] = NIE_inference[(2*j-1),"Estimate"] / NIE_inference[(2*j-1),"SE"]
      NIE_inference[(2*j-1),"p value"] = 2*pnorm( -abs( NIE_inference[(2*j-1),"test-stat"] )  )
      
      if(length(A_D[[j]])==0){
        NIE_inference[(2*j),] = NIE_inference[(2*j-1),]
      } else {
        
        E_A_G = matrix(0, nrow = length(A_G), ncol = nrow(Sigmahat))
        colnames(E_A_G) = colnames(Sigmahat)
        for (h in 1:nrow(E_A_G)) {
          E_A_G[h, A_G_cpg[h]]=1
        }
        
        
        E_A_Dj = matrix(0, nrow = length(A_D[[j]]), ncol = nrow(Sigmahat))
        colnames(E_A_Dj) = colnames(Sigmahat)
        for (h in 1:nrow(E_A_Dj)) {
          E_A_Dj[h, A_D_cpg[[j]][h] ] = 1
        }
        
        theta.GDj_A = c(theta.G_hat_AG, Theta.D_hat_AD[[j]]) 
        Sigma_G_hat_GDj.A = Sigma_G_hat[c(A_G_cpg, A_D_cpg_forunion[[j]]), c(A_G_cpg, A_D_cpg_forunion[[j]])]
        
        var_indrct_j1 = as.numeric( sigma2_Y * (t(betahat[(j+1),A_G])%*%E_A_G + t(betahat[(j+1),A_D[[j]]])%*%E_A_Dj ) %*% (Sigmahat.inv) %*% ( t(E_A_G)%*%betahat[(j+1),A_G] + t(E_A_Dj)%*%betahat[(j+1),A_D[[j]] ] ) ) +
          as.numeric((solve(J)[(j+1),(j+1)])* t(theta.GDj_A) %*% (Sigma_G_hat_GDj.A)%*% (theta.GDj_A))
        
        NIE_inference[(2*j),"Estimate"] = (beta.theta.G[(j+1),1] + beta.Theta.D[[j]][(j+1),1] )
        NIE_inference[(2*j),"SE"] = sqrt(var_indrct_j1)/sqrt(n)
        NIE_inference[(2*j),"test-stat"] = NIE_inference[(2*j),"Estimate"] / NIE_inference[(2*j),"SE"]
        NIE_inference[(2*j),"p value"] = 2*pnorm( -abs( NIE_inference[(2*j),"test-stat"] )  )
        
      }
      
    }
    
    
    
    
  }
  
  
  for (j in 1:q) {
    rownames(NDE_inference)[(2*j-1)] = paste("NDE(T", j, "= 0)")
    rownames(NDE_inference)[(2*j)] = paste("NDE(T", j, "= 1)")
    
    if(length(A_D[[j]])==0){
      
      var_drct_j = as.numeric(sigma2_Y *( Sigmahat.inv[ (j+1),(j+1) ] ) )
      NDE_inference[(2*j-1), "Estimate"] = theta.T_hat[(j)]
      NDE_inference[(2*j-1), "SE"] = sqrt(var_drct_j)/sqrt(n)
      NDE_inference[(2*j-1), "test-stat"] = NDE_inference[(2*j-1), "Estimate"] / NDE_inference[(2*j-1), "SE"]
      NDE_inference[(2*j-1), "p value"] = 2*pnorm(-abs( NDE_inference[(2*j-1), "test-stat"] ))
      NDE_inference[(2*j),] = NDE_inference[(2*j-1),]
      
    } else {
      
      Xtilde = apply(Covariate, 2, mean)
      
      e_tilde.jp1 = matrix(0, nrow=nrow(Sigmahat), ncol=1)
      e_tilde.jp1[(j+1),1]=1
      
      E_A_Dj = matrix(0, nrow = length(A_D[[j]]), ncol = nrow(Sigmahat))
      colnames(E_A_Dj) = colnames(Sigmahat)
      for (h in 1:nrow(E_A_Dj)) {
        E_A_Dj[h, A_D_cpg[[j]][h]] = 1
      }
      
      alpha.hat.j = e_tilde.jp1 + t(E_A_Dj)%*%( betahat[1,A_D_cpg_forunion[[j]] ]) + t(E_A_Dj)%*%t(betahat[-c(1:(1+q)), A_D_cpg_forunion[[j]] ])%*% Xtilde
      
      
      e1 = matrix(0, nrow=nrow(J), ncol=1)
      e1[1,1]=1
      
      E_X = matrix(0, nrow= ncol(Covariate), ncol=ncol(J))
      for (h in 1:nrow(E_X)) {
        E_X[h, (h+1+q)] = 1
      }
      
      var_drct_j0 = as.numeric( sigma2_Y *( t(alpha.hat.j)%*%Sigmahat.inv %*%alpha.hat.j)) + 
        as.numeric(t(e1 + t(E_X)%*%Xtilde)%*%solve(J)%*%(e1 + t(E_X)%*%Xtilde))* as.numeric(t( Theta.D_hat_AD[[j]]  ) %*% (Sigma_G_hat[A_D_cpg_forunion[[j]], A_D_cpg_forunion[[j]] ])%*% ( Theta.D_hat_AD[[j]] ))
      
      NDE_inference[(2*j-1), "Estimate"] = theta.T_hat[j] + beta.Theta.D[[j]][1,1] + t(Xtilde)%*%( beta.Theta.D[[j]][-c(1:(1+q)),1] )
      NDE_inference[(2*j-1), "SE"] = sqrt(var_drct_j0)/sqrt(n)
      NDE_inference[(2*j-1), "test-stat"] = NDE_inference[(2*j-1), "Estimate"] / NDE_inference[(2*j-1), "SE"]
      NDE_inference[(2*j-1), "p value"] = 2*pnorm(-abs( NDE_inference[(2*j-1), "test-stat"] ))
      
      
      e.jp1 = matrix(0, nrow=nrow(J), ncol=1)
      e.jp1[(j+1),1]=1
      
      var_drct_j1 = as.numeric(sigma2_Y *( t( alpha.hat.j + t(E_A_Dj)%*%betahat[(j+1),A_D[[j]] ] )%*%Sigmahat.inv %*%( alpha.hat.j + t(E_A_Dj)%*%betahat[(j+1),A_D[[j]] ] ) )) + 
        as.numeric(t(e1 + t(E_X)%*%Xtilde + e.jp1)%*%solve(J)%*%(e1 + t(E_X)%*%Xtilde + e.jp1)) * as.numeric(t( Theta.D_hat_AD[[j]]  ) %*% (Sigma_G_hat[A_D_cpg_forunion[[j]], A_D_cpg_forunion[[j]] ])%*% ( Theta.D_hat_AD[[j]] ))
      
      NDE_inference[(2*j), "Estimate"] = theta.T_hat[j] + beta.Theta.D[[j]][1,1] + t(Xtilde)%*%( beta.Theta.D[[j]][-c(1:(1+q)),1] ) + beta.Theta.D[[j]][(1+j),1]
      NDE_inference[(2*j), "SE"] = sqrt(var_drct_j1)/sqrt(n)
      NDE_inference[(2*j), "test-stat"] = NDE_inference[(2*j), "Estimate"] / NDE_inference[(2*j), "SE"]
      NDE_inference[(2*j), "p value"] = 2*pnorm(-abs( NDE_inference[(2*j), "test-stat"] ))
      
      
      
    }
    
    
    
    
    
    
  }
  
  return( list(outcome.theta_inference = Outcome.theta_inference,
               NIE_inference = NIE_inference,
               NDE_inference = NDE_inference,
               W_A = W_A, A_G = A_G, A_D = A_D)
  )
  
}




















##### The function to calculate NIE, NDE SE's (population version) with known beta, theta, and support ###
### This function is only used in theoretical simulation ###
### Note that the oracle method has the same NIE, NDE, individual theta, SE but the point estimation is calculated by two low-dim OLSs rather than SCAD + post selected OLS ###
### Note that this SE is the SE for \sqrt{n} (effecthat - effect)
### This function contains the covariate X
OraProjInt_SE_popu_X = function(theta0, theta.T, theta.G, Theta.D, theta.X,
                              beta0, B, A, Sigma_G, mu.X, sigma.X, sigma_Y,
                              Xtilde,
                              prop){
  
  #' @param theta0 a 1 by 1 matrix 
  #' @param theta.T a matrix with nrow = q, ncol = 1
  #' @param theta.X a matrix with nrow = r, ncol = 1
  #' @param theta.G a matrix with nrow = p, ncol =  1
  #' @param Theta.D a matrix with nrow = p, ncol = q
  #' @param beta0 a matrix with nrow = p, ncol = 1
  #' @param B a matrix with nrow = p, ncol = q
  #' @param A a matrix with nrow = p, ncol= r
  #' @param Sigma_G a p by p matrix, covariance matrix of bepsilon_G
  #' @param mu.X a matrix with nrow = r, ncol = 1
  #' @param sigma.X a matrix with nrow = r, ncol = r
  #' @param sigma_Y, a scalar (vector with length 1) which is the standard deviation of epsilon_Yi
  #' @param Xtilde, a matrix with nrow= r, ncol =1 used in the expression of direct effect, usually this value is mu.X
  #' @param prop is a vector with length q, whose entries are the p_j = lim( n_j/n ), j=1,...,q, the proportion of the j-th treatment. Note that sum(prop) < 1 since we have p_0 not included.
  #' The entries in theta0, theta.T, theta.X, theta.G, Theta.D, beta0, B, A, Sigma_G, mu.X, sigma.X have no rownames or colnames
  #' The theta's can have different nonzero entry locations, not limited to heredity
  #' @return NIE_SE, NDE_SE are the matrix of SE of the NIE and NDE for different j

  
  q = ncol(Theta.D)
  r = ncol(A)
  
  # First find the support of theta.G and Theta.D
  # which function will first transform the argument theta.G or Theta.D into a vector connected by columns
  # Then give an vector output with nonzero locations in numbers
  A_G = which(theta.G != 0) # A_G is a vector, containing number location in 1:p with length s.G
  A_D = list() # A_D has q components. Each components contains number location in 1:p with length s.Dj
  for (j in 1:q) {
    A_D[[j]] = which(Theta.D[,j]!=0)
  }
  
  s.G = length(A_G)
  
  s.D = c()
  for (j in 1:q) {
    s.D[j] = length(A_D[[j]])
  }
  
  Sigma11 = matrix(0, nrow = (1+q), ncol = (1+q))
  Sigma11[1,] = c(1, prop)
  for (j in 1:q) {
    Sigma11[(j+1), 1] = prop[j]
    Sigma11[(j+1), (j+1)] = prop[j]
  }
  
  prop0 = 1 -  sum(prop)
  
  Sigma12 = matrix(0, nrow = (1+q), ncol = r)
  Sigma12[1,] = t(mu.X)
  for (j in 1:q) {
    Sigma12[(j+1), ] = prop[j]*t(mu.X)
  }
  
  Sigma13 = matrix(0, nrow = (1+q), ncol = s.G)
  Sigma13[1, ] = prop0 * t( beta0[A_G,1] + A[A_G, ]%*%mu.X )
  for (j in 1:q) {
    Sigma13[1, ] = Sigma13[1, ] + prop[j] * t( beta0[A_G,1] + B[A_G,j] + A[A_G, ]%*%mu.X  )
  }
  for (j in 1:q) {
    Sigma13[(j+1), ] = prop[j] * t( beta0[A_G,1] + B[A_G,j] + A[A_G, ]%*%mu.X  )
  }
  
  Sigma14 = matrix(0, nrow = (1+q), ncol = sum(s.D))
  if(q == 1){
    Sigma14[1, 1:s.D[1]] = prop[1] * t( beta0[A_D[[1]], 1] + B[ A_D[[1]] , 1] + A[ A_D[[1]] ,  ]%*%mu.X )
    Sigma14[2, 1:s.D[1]] = prop[1] * t( beta0[A_D[[1]], 1] + B[ A_D[[1]] , 1] + A[ A_D[[1]] ,  ]%*%mu.X )
  } else {
    Sigma14[1, 1:s.D[1]] = prop[1] * t( beta0[A_D[[1]], 1] + B[ A_D[[1]] , 1] + A[ A_D[[1]] ,  ]%*%mu.X )
    Sigma14[2, 1:s.D[1]] = prop[1] * t( beta0[A_D[[1]], 1] + B[ A_D[[1]] , 1] + A[ A_D[[1]] ,  ]%*%mu.X )
    for (j in 2:q) {
      Sigma14[1, (sum(s.D[1:(j-1)])+1) : sum(s.D[1:j]) ] = prop[j] * t( beta0[A_D[[j]], 1] + B[ A_D[[j]] , j] + A[ A_D[[j]] ,  ]%*%mu.X ) # just because it has 1:(j-1), if j=1, 1:(j-1) is meaningless, so I have to separate j=1 and j>=2. Actually, I do not need such procedure. I can first define j=1 part. Then add if (q>1){ loop the j>=2 part  } like Sigma34. But I have already done it here so I won't change.
  Sigma14[(j+1), (sum(s.D[1:(j-1)])+1) : sum(s.D[1:j]) ] = prop[j] * t( beta0[A_D[[j]], 1] + B[ A_D[[j]] , j] + A[ A_D[[j]] ,  ]%*%mu.X )
    }
  }
  
  
  Sigma21 = t(Sigma12)
  
  Sigma22 = matrix(0, r, r)
  Sigma22 = sigma.X + mu.X%*%t(mu.X)
  
  Sigma23 = prop0 * ( mu.X%*%t( beta0[ A_G , 1] + A[ A_G , ]%*%mu.X  ) + sigma.X%*%t(A[ A_G , ]) )
  for (j in 1:q) {
    Sigma23 = Sigma23 + prop[j] * ( mu.X%*%t( beta0[ A_G , 1] + B[A_G,j] + A[ A_G , ]%*%mu.X  ) + sigma.X%*%t(A[ A_G , ]) )
  }
  
  
  Sigma24 = matrix(0, r, sum(s.D))
  if(q == 1){
    Sigma24[ , 1:s.D[1]] = prop[1] * ( mu.X%*%t( beta0[ A_D[[1]], 1] + B[A_D[[1]],1] + A[ A_D[[1]], ]%*%mu.X  ) + sigma.X%*%t(A[A_D[[1]],])  )
  } else {
    Sigma24[ , 1:s.D[1]] = prop[1] * ( mu.X%*%t( beta0[ A_D[[1]], 1] + B[A_D[[1]],1] + A[ A_D[[1]], ]%*%mu.X  ) + sigma.X%*%t(A[A_D[[1]],])  ) 
    for (j in 2:q) {
      if(length(A_D[[j]]) == 1){
        Sigma24[, (sum(s.D[1:(j-1)])+1) : sum(s.D[1:j])] = prop[j] * ( mu.X%*%t( beta0[ A_D[[j]], 1] + B[A_D[[j]],j] + A[ A_D[[j]], ]%*%mu.X  ) + sigma.X%*%A[A_D[[j]],]  ) 
      } else {
        Sigma24[, (sum(s.D[1:(j-1)])+1) : sum(s.D[1:j])] = prop[j] * ( mu.X%*%t( beta0[ A_D[[j]], 1] + B[A_D[[j]],j] + A[ A_D[[j]], ]%*%mu.X  ) + sigma.X%*%t(A[A_D[[j]],])  ) 
      }
    }
  }
  
  
  Sigma31 = t(Sigma13)
  
  Sigma32 = t(Sigma23)
  
  Sigma33 = prop0 * ( ( beta0[A_G,1] + A[A_G,]%*%mu.X )%*%t( beta0[A_G,1] + A[A_G,]%*%mu.X ) + A[A_G, ]%*%sigma.X%*%t(A[A_G,]) + Sigma_G[A_G, A_G]  )
  for (j in 1:q) {
    Sigma33 =  Sigma33 + prop[j] * ( ( beta0[A_G,1] + B[A_G,j] + A[A_G,]%*%mu.X )%*%t( beta0[A_G,1] + B[A_G,j] + A[A_G,]%*%mu.X ) + A[A_G, ]%*%sigma.X%*%t(A[A_G,]) + Sigma_G[A_G, A_G]  )
  }
  
  
  Sigma34 = matrix(0, nrow = s.G, ncol=sum(s.D) )
  Sigma34[ , 1:s.D[1]] = prop[1] *  ( ( beta0[ A_G, 1] + B[A_G, 1] + A[A_G, ]%*%mu.X  )%*%t( beta0[ A_D[[1]], 1] + B[A_D[[1]],1] + A[ A_D[[1]], ]%*%mu.X  ) + A[A_G, ]%*%sigma.X%*%t(A[A_D[[1]],]) + Sigma_G[ A_G, A_D[[1]] ]  ) 
  if(q > 1){
    for (j in 2:q) {
      if(length(A_D[[j]])==1){
        Sigma34[, (sum(s.D[1:(j-1)])+1) : sum(s.D[1:j])] = prop[j] *  ( ( beta0[ A_G, 1] + B[A_G, j] + A[A_G, ]%*%mu.X  )%*%t( beta0[ A_D[[j]], 1] + B[A_D[[j]],j] + A[ A_D[[j]], ]%*%mu.X  ) + A[A_G, ]%*%sigma.X%*%A[A_D[[j]],] + Sigma_G[ A_G, A_D[[j]] ]  ) 
        
      } else {
        Sigma34[, (sum(s.D[1:(j-1)])+1) : sum(s.D[1:j])] = prop[j] *  ( ( beta0[ A_G, 1] + B[A_G, j] + A[A_G, ]%*%mu.X  )%*%t( beta0[ A_D[[j]], 1] + B[A_D[[j]],j] + A[ A_D[[j]], ]%*%mu.X  ) + A[A_G, ]%*%sigma.X%*%t(A[A_D[[j]],]) + Sigma_G[ A_G, A_D[[j]] ]  ) 
        
      }
        
    }
  }
  
  
  Sigma41 = t(Sigma14)
  Sigma42 = t(Sigma24)
  Sigma43 = t(Sigma34)
  
  Sigma44 = matrix(0, nrow = sum(s.D), ncol = sum(s.D)  )
  Sigma44[ 1:s.D[1], 1:s.D[1] ] = prop[1] *  ( ( beta0[ A_D[[1]], 1] + B[A_D[[1]], 1] + A[A_D[[1]], ]%*%mu.X  )%*%t( beta0[ A_D[[1]], 1] + B[A_D[[1]],1] + A[ A_D[[1]], ]%*%mu.X  ) + A[A_D[[1]], ]%*%sigma.X%*%t(A[A_D[[1]],]) + Sigma_G[ A_D[[1]], A_D[[1]] ]  ) 
  if(q > 1){
    for (j in 2:q) {
      if(length(A_D[[j]]) == 1){
        Sigma44[ (sum(s.D[1:(j-1)])+1) : sum(s.D[1:j]) , (sum(s.D[1:(j-1)])+1) : sum(s.D[1:j]) ] = prop[j] *  ( ( beta0[ A_D[[j]], 1] + B[A_D[[j]], j] + A[A_D[[j]], ]%*%mu.X  )%*%t( beta0[ A_D[[j]], 1] + B[A_D[[j]],j] + A[ A_D[[j]], ]%*%mu.X  ) + t(A[A_D[[j]], ])%*%sigma.X%*%A[A_D[[j]],] + Sigma_G[ A_D[[j]], A_D[[j]] ]  ) 
      } else {
        Sigma44[ (sum(s.D[1:(j-1)])+1) : sum(s.D[1:j]) , (sum(s.D[1:(j-1)])+1) : sum(s.D[1:j]) ] = prop[j] *  ( ( beta0[ A_D[[j]], 1] + B[A_D[[j]], j] + A[A_D[[j]], ]%*%mu.X  )%*%t( beta0[ A_D[[j]], 1] + B[A_D[[j]],j] + A[ A_D[[j]], ]%*%mu.X  ) + A[A_D[[j]], ]%*%sigma.X%*%t(A[A_D[[j]],]) + Sigma_G[ A_D[[j]], A_D[[j]] ]  ) 
      }
      
    }
  }
  
  Sigma = rbind(cbind(Sigma11, Sigma12, Sigma13, Sigma14),
                cbind(Sigma21, Sigma22, Sigma23, Sigma24),
                cbind(Sigma31, Sigma32, Sigma33, Sigma34),
                cbind(Sigma41, Sigma42, Sigma43, Sigma44))
  
  Sigma.inv=solve(Sigma)
  J=rbind(cbind(Sigma11, Sigma12),
          cbind(Sigma21, Sigma22))
  J.inv=solve(J)
  
  
  NIE_SE = matrix(NA, nrow = q, ncol = 2)
  colnames(NIE_SE) = c("Tj=0", "Tj=1")
  rownames(NIE_SE) = paste0("j=", 1:q)
  
  NDE_SE = matrix(NA, nrow = q, ncol = 2)
  colnames(NDE_SE) = c("Tj=0", "Tj=1")
  rownames(NDE_SE) = paste0("j=", 1:q)
  
  E_AG = matrix(0, nrow = s.G, ncol = ncol(Sigma))
  for (h in 1:s.G) {
    E_AG[h, 1+q+r+h]=1
  }
  
  E_AD1 = matrix(0, nrow = s.D[1], ncol = ncol(Sigma))
  for (h in 1:s.D[1]) {
    E_AD1[h, 1+q+r+s.G+h] = 1
  }
  
  var.NIE0 = (sigma_Y^2)* ( t(B[A_G,1])%*% E_AG %*% (Sigma.inv) %*% t(E_AG)%*%B[A_G,1]) + (J.inv[2,2]) * (t(theta.G[A_G,1])%*%(Sigma_G[A_G,A_G])%*%theta.G[A_G,1])
  
  NIE_SE[1,"Tj=0"] = sqrt(var.NIE0)
  rm(var.NIE0)
  
  var.NIE1 = (sigma_Y^2)* ( ( t(B[A_G,1])%*%E_AG + t(B[A_D[[1]],1])%*%E_AD1 )%*%Sigma.inv%*%( t(E_AG)%*%B[A_G,1] + t(E_AD1)%*%B[A_D[[1]],1] ) ) + 
    (J.inv[2,2]) *  ( t(c(theta.G[A_G,1], Theta.D[A_D[[1]],1]))%*% ( rbind( cbind(Sigma_G[A_G,A_G], Sigma_G[A_G,A_D[[1]]]),
                                                                     cbind(Sigma_G[A_D[[1]],A_G], Sigma_G[A_D[[1]],A_D[[1]]])) )%*%c(theta.G[A_G,1], Theta.D[A_D[[1]],1]) )
  
  NIE_SE[1, "Tj=1"] = sqrt(var.NIE1)
  rm(var.NIE1)
  
  etilde.2 = matrix(0, nrow = nrow(Sigma), ncol=1)
  etilde.2[2,1]=1
  
  e1 = matrix(0, nrow=(1+q+r), ncol=1)
  e1[1,1] = 1
  
  e2 = matrix(0, nrow=(1+q+r), ncol=1)
  e2[2,1] = 1
  
  alpha = etilde.2 + t(E_AD1)%*%beta0[A_D[[1]],1] + t(E_AD1)%*%A[A_D[[1]],]%*%Xtilde
  
  E_X = matrix(0, r, (1+q+r))
  for (m in 1:r) {
    E_X[m, 1+q+m] = 1
  }
  
  var.NDE0 = (sigma_Y^2) * (t(alpha)%*%Sigma.inv%*%alpha) + ( t(e1 + t(E_X)%*%Xtilde)%*%J.inv%*%(e1 + t(E_X)%*%Xtilde) ) * ( t(Theta.D[A_D[[1]],1])%*%Sigma_G[ A_D[[1]],A_D[[1]] ]%*%Theta.D[A_D[[1]],1] )
  NDE_SE[1, "Tj=0"] = sqrt(var.NDE0)
  rm(var.NDE0)
  
  var.NDE1 = (sigma_Y^2) * ( t(alpha + t(E_AD1)%*%B[A_D[[1]],1])%*%Sigma.inv%*%(alpha + t(E_AD1)%*%B[A_D[[1]],1]) ) + 
    ( t(e1 + t(E_X)%*%Xtilde + e2)%*%J.inv%*%(e1 + t(E_X)%*%Xtilde + e2) ) * ( t(Theta.D[A_D[[1]],1])%*%Sigma_G[A_D[[1]],A_D[[1]]]%*%Theta.D[A_D[[1]],1] )
  NDE_SE[1, "Tj=1"] = sqrt(var.NDE1)
  rm(var.NDE1)
  
  rm(alpha)
  
  if(q > 1){
    for (j in 2:q) {
      
      E_ADj = matrix(0, nrow = s.D[j], ncol = ncol(Sigma))
      for (h in 1:s.D[j]) {
        E_ADj[h, ( 1+q+r+s.G + sum(s.D[1:(j-1)]) + h ) ] = 1
      }
      
      var.NIE0 = (sigma_Y^2)* ( t(B[A_G,j])%*% E_AG %*% (Sigma.inv) %*% t(E_AG)%*%B[A_G,j]) + (J.inv[(j+1),(j+1)]) * (t(theta.G[A_G,1])%*%(Sigma_G[A_G,A_G])%*%theta.G[A_G,1])
      
      NIE_SE[j,"Tj=0"] = sqrt(var.NIE0)
      rm(var.NIE0)
      
      if(length(A_D[[j]]) == 1){
        var.NIE1 = (sigma_Y^2)* ( ( t(B[A_G,j])%*%E_AG + t(B[A_D[[j]],j])%*%E_ADj )%*%Sigma.inv%*%( t(E_AG)%*%B[A_G,j] + t(E_ADj)%*%B[A_D[[j]],j] ) ) + 
          (J.inv[(j+1),(j+1)]) *  ( t(c(theta.G[A_G,1], Theta.D[A_D[[j]],j]))%*% ( rbind( cbind(Sigma_G[A_G,A_G], Sigma_G[A_G,A_D[[j]]]),
                                                                                          cbind( t(as.matrix(Sigma_G[A_D[[j]],A_G])), Sigma_G[A_D[[j]],A_D[[j]]]) ) )%*%c(theta.G[A_G,1], Theta.D[A_D[[j]],j]) )
        
      } else {
        var.NIE1 = (sigma_Y^2)* ( ( t(B[A_G,j])%*%E_AG + t(B[A_D[[j]],j])%*%E_ADj )%*%Sigma.inv%*%( t(E_AG)%*%B[A_G,j] + t(E_ADj)%*%B[A_D[[j]],j] ) ) + 
          (J.inv[(j+1),(j+1)]) *  ( t(c(theta.G[A_G,1], Theta.D[A_D[[j]],j]))%*% ( rbind( cbind(Sigma_G[A_G,A_G], Sigma_G[A_G,A_D[[j]]]),
                                                                                          cbind(Sigma_G[A_D[[j]],A_G], Sigma_G[A_D[[j]],A_D[[j]]]) ) )%*%c(theta.G[A_G,1], Theta.D[A_D[[j]],j]) )
        
      }
      
      NIE_SE[j, "Tj=1"] = sqrt(var.NIE1)
      rm(var.NIE1)
      
      etilde.jp1 = matrix(0, nrow = nrow(Sigma), ncol=1)
      etilde.jp1[(j+1),1]=1
      
      ejp1 = matrix(0, nrow=(1+q+r), ncol=1)
      ejp1[(j+1),1] = 1
      
      alpha = etilde.jp1 + t(E_ADj)%*%beta0[A_D[[j]],1] + t(E_ADj)%*%A[A_D[[j]],]%*%Xtilde
      
      var.NDE0 = (sigma_Y^2) * (t(alpha)%*%Sigma.inv%*%alpha) + ( t(e1 + t(E_X)%*%Xtilde)%*%J.inv%*%(e1 + t(E_X)%*%Xtilde) ) * ( t(Theta.D[A_D[[j]],j])%*%Sigma_G[ A_D[[j]],A_D[[j]] ]%*%Theta.D[A_D[[j]],j] )
      NDE_SE[j, "Tj=0"] = sqrt(var.NDE0)
      rm(var.NDE0)
      
      var.NDE1 = (sigma_Y^2) * ( t(alpha + t(E_ADj)%*%B[A_D[[j]],j])%*%Sigma.inv%*%(alpha + t(E_ADj)%*%B[A_D[[j]],j]) ) + 
        ( t(e1 + t(E_X)%*%Xtilde + ejp1)%*%J.inv%*%(e1 + t(E_X)%*%Xtilde + ejp1) ) * ( t(Theta.D[A_D[[j]],j])%*%Sigma_G[A_D[[j]],A_D[[j]]]%*%Theta.D[A_D[[j]],j] )
      NDE_SE[j, "Tj=1"] = sqrt(var.NDE1)
      rm(var.NDE1)
      
      rm(alpha)
      
    }
  }
  
  return( list(NIE_SE = NIE_SE, NDE_SE = NDE_SE)
  )
  
}
























### The function to calculate NIE, NDE, and individual theta, SE's (sample version) ### 
## this function does not the covariate X ##
OraProjInt_SE_sample_noX = function(Y, Trt, Mediator, Trt.Medator,
                                  theta0_hat, theta.T_hat, theta.G_hat, vtheta.D_hat){
  
  #' @param Y vector, n-dimensional outcome vector (with column name such as "Y").
  #' @param Trt matrix or dataframe, n by q treatment matrix without intercept (with column name such as "T1", "T2",.., "Tq").
  #' @param Mediator matrix or dataframe, n by p mediator matrix (with column name, "G1",...).
  #' @param Trt.Medator matrix or dataframe, n by qp interaction matrix (with column name such as "T1.G1", "T1.G2",...).
  #' @param theta0_hat a vector with length 1 from SCAD estimation
  #' @param theta.T_hat a vector with length q from SCAD estimation
  #' @param theta.G_hat a vector of length p from SCAD estimation
  #' @param vtheta.D_hat a vector of length pq, which is vec{Theta_D}, from SCAD estimation
  #' Note that Trt, Mediator, Trt.Medator must have column names in words!
  #' theta0_hat, theta.T_hat, theta.G_hat, vtheta.D_hat are vectors with no names! Just pure vector
  #' If you use ncvreg packages, it will produce names of theta0_hat, theta.T_hat...
  #' To remove those names of theta, use the code before plugging them in this function such as: names(theta0.hat)=NULL; names(theta.T.hat)=NULL; names(theta.X.hat)=NULL; names(theta.G.hat)=NULL; names(vtheta.D.hat)=NULL
  #' @return outcome.theta_inference is a matrix of Estimate, SE, test-stat, p value of the SCAD outcome model inference for nonzero theta's
  #' @return NIE_inference, the matrix of Estimate, SE, test-stat, p value of the NIE estimators
  #' @return NDE_inference, the matrix of Estimate, SE, tests-stat, p value of NDE estimators
  #' @return W_A = marix(one.n, T, G[,A_G], T.G[,A_D1,..., A_Dq])
  #' @return A_G: a vector of length \hat{s}_G, nonzero locations (number index in 1:p) in theta_G
  #' @return A_D: a list of vectors. Each vector is the nonzero location (number index in 1:p) in theta_{D_j} with length \hat{s}_{D_j}
  
  
  p = ncol(Mediator)
  n = nrow(Mediator)
  q = ncol(Trt)
 
  
  A_G = which(theta.G_hat!=0) # A_G contains number location in 1:p with length \hat{s_G}
  A_G_cpg = colnames(Mediator)[A_G] # A_G contains word names of those support mediators with length \hat{s.G}
  
  A_D = list() # A_D has q components. Each components contains number location in 1:p with length \hat{s_{D_j}}
  A_D_cpg = list() # This contains the word names of those support mediators with length \hat{s_{D_j}}. These names are in the form of Trt*Medator name
  A_D_cpg_forunion=list() # This contains the word names of those support mediators with length \hat{s_{D_j}}. These names are still the mediator name rather than Trt*Medator name
  for (j in 1:q) {
    A_D[[j]] = which(vtheta.D_hat[ ((j-1)*p+1) : (j*p) ]!=0)
    tmp = Trt.Medator[, ((j-1)*p+1) : (j*p) ]
    A_D_cpg[[j]] = colnames(tmp)[A_D[[j]]]
    A_D_cpg_forunion[[j]] = colnames(Mediator)[A_D[[j]]] # if there is no support, length(A_D[[j]]) = 0, length(A_D_cpg_forunion[[j]])=0
  }
  rm(tmp)
  
  one.n = matrix(1, nrow=nrow(Mediator), ncol=1)
  
  W_A = cbind(one.n, Trt, Mediator[,A_G])
  for (j in 1:q) {
    tmp = Trt.Medator[, ((j-1)*p+1) : (j*p) ]
    W_A = cbind(W_A, tmp[, A_D[[j]] ])
  }
  
  rm(tmp)
  
  tmpname = c()
  tmpname[1] = "one"
  tmpname[(1+1):(1+q)] = colnames(Trt)
  if(length(A_G) > 0){
    tmpname[(1+q+1):(1+q+length(A_G))] = A_G_cpg
    for (j in 1:q) {
      tmpname = c(tmpname, A_D_cpg[[j]]) # if length(A_D_cpg[[j]])=0, this step will not change the previous tmpname
    }
  } else {
    for (j in 1:q) {
      tmpname = c(tmpname, A_D_cpg[[j]])
    }
  }
  colnames(W_A) = tmpname
  rm(tmpname)
  
  
  Sigmahat = (t(W_A)%*%W_A) / nrow(Mediator)
  Sigmahat.inv = solve(Sigmahat)
  
  eps_Y = Y - one.n%*%theta0_hat - Trt%*%theta.T_hat - Mediator%*%theta.G_hat - Trt.Medator%*%vtheta.D_hat
  sigma2_Y = as.numeric( t(eps_Y)%*%eps_Y/n ) # the estimate of sigma_Y^2
  rm(eps_Y)
  
  theta.G_hat_AG = theta.G_hat[A_G]
  Theta.D_hat_AD = list()
  for (j in 1:q) {
    tmp = vtheta.D_hat[ ((j-1)*p+1) : (j*p) ]
    Theta.D_hat_AD[[j]] = tmp[A_D[[j]]] # if A_D[[j]] is empty, length(Theta.D_hat_AD[[j]])=0
  }
  rm(tmp)
  
  # thetahat only store those nonzero theta's
  thetahat = c(theta0_hat, theta.T_hat,  theta.G_hat_AG)
  for (j in 1:q) {
    thetahat = c(thetahat, Theta.D_hat_AD[[j]])
  }
  
  asymp.var_SCAD = sigma2_Y * Sigmahat.inv
  
  ### Outcome.theta_inference is the individual theta's inference from the SCAD step in the outcome model
  Outcome.theta_inference = cbind( thetahat, 
                                   sqrt(diag(asymp.var_SCAD))/sqrt(n), 
                                   sqrt(n)*thetahat /sqrt(diag(asymp.var_SCAD)),  
                                   2*pnorm(-abs( sqrt(n)*thetahat /sqrt(diag(asymp.var_SCAD)) )) )
  
  colnames(Outcome.theta_inference) = c("Estimate", "SE", "test-stat", "p value")
  
  
  ###
  H = cbind(one.n, Trt)
  
  betahat = solve(t(H)%*%H)%*%t(H)%*%Mediator
  
  eps_G_hat = Mediator - H%*%betahat
  Sigma_G_hat = t(eps_G_hat)%*%eps_G_hat / n
  
  J = Sigmahat[1:(1+q),1:(1+q)]
  
  
  beta.theta.G = betahat %*% theta.G_hat # (1+q) by 1
  beta.Theta.D = list() # Each component is (1+q) by 1
  for (j in 1:q) {
    beta.Theta.D[[j]] = betahat %*% vtheta.D_hat[ ((j-1)*p+1) : (j*p) ] # if there is no active support, this vector could be zero
  }
  
  NIE_inference = matrix(NA, nrow = (q*2), ncol = 4)
  colnames(NIE_inference) = c("Estimate", "SE", "test-stat", "p value")
  rownames(NIE_inference) = rep("NIE", times = (q*2))
  
  for (j in 1:q) {
    rownames(NIE_inference)[(2*j-1)] = paste("NIE(T", j, "= 0)")
    rownames(NIE_inference)[(2*j)] = paste("NIE(T", j, "= 1)")
    
    if(length(A_G)==0){
      NIE_inference[(2*j-1),"Estimate"] = 0
      NIE_inference[(2*j-1),"SE"] = 0
      NIE_inference[(2*j-1),"p value"] = 1
      
      if(length(A_D[[j]])==0){
        NIE_inference[(2*j),] = NIE_inference[(2*j-1),]
      } else {
        
        E_A_Dj = matrix(0, nrow = length(A_D[[j]]), ncol = nrow(Sigmahat))
        colnames(E_A_Dj) = colnames(Sigmahat)
        for (h in 1:nrow(E_A_Dj)) {
          E_A_Dj[h, A_D_cpg[[j]][h] ] = 1
        }
        
        theta.GDj_A = c(theta.G_hat_AG, Theta.D_hat_AD[[j]]) 
        Sigma_G_hat_GDj.A = Sigma_G_hat[c(A_G_cpg, A_D_cpg_forunion[[j]]), c(A_G_cpg, A_D_cpg_forunion[[j]])]
        
        var_indrct_j1 = as.numeric( sigma2_Y * (t(betahat[(j+1),A_D[[j]]])%*%E_A_Dj ) %*% (Sigmahat.inv) %*% ( t(E_A_Dj)%*%betahat[(j+1),A_D[[j]] ] ) ) +
          as.numeric((solve(J)[(j+1),(j+1)])* t(theta.GDj_A) %*% (Sigma_G_hat_GDj.A)%*% (theta.GDj_A))
        
        NIE_inference[(2*j),"Estimate"] = (beta.theta.G[(j+1),1] + beta.Theta.D[[j]][(j+1),1] )
        NIE_inference[(2*j),"SE"] = sqrt(var_indrct_j1)/sqrt(n)
        NIE_inference[(2*j),"test-stat"] = NIE_inference[(2*j),"Estimate"] / NIE_inference[(2*j),"SE"]
        NIE_inference[(2*j),"p value"] = 2*pnorm( -abs( NIE_inference[(2*j),"test-stat"] )  )
        
      }
      
      
    } else {
      
      var_indrct_j0 = as.numeric( (sigma2_Y) * t(betahat[(j+1),A_G]) %*% (Sigmahat.inv[A_G_cpg, A_G_cpg]) %*% (betahat[(j+1),A_G]) ) +
        as.numeric( (solve(J)[(j+1),(j+1)]) * t(theta.G_hat[A_G]) %*% (Sigma_G_hat[A_G, A_G])%*% (theta.G_hat[A_G])) 
      NIE_inference[(2*j-1),"Estimate"] = beta.theta.G[(j+1),1]
      NIE_inference[(2*j-1),"SE"] = sqrt(var_indrct_j0)/sqrt(n)
      NIE_inference[(2*j-1),"test-stat"] = NIE_inference[(2*j-1),"Estimate"] / NIE_inference[(2*j-1),"SE"]
      NIE_inference[(2*j-1),"p value"] = 2*pnorm( -abs( NIE_inference[(2*j-1),"test-stat"] )  )
      
      if(length(A_D[[j]])==0){
        NIE_inference[(2*j),] = NIE_inference[(2*j-1),]
      } else {
        
        E_A_G = matrix(0, nrow = length(A_G), ncol = nrow(Sigmahat))
        colnames(E_A_G) = colnames(Sigmahat)
        for (h in 1:nrow(E_A_G)) {
          E_A_G[h, A_G_cpg[h]]=1
        }
        
        
        E_A_Dj = matrix(0, nrow = length(A_D[[j]]), ncol = nrow(Sigmahat))
        colnames(E_A_Dj) = colnames(Sigmahat)
        for (h in 1:nrow(E_A_Dj)) {
          E_A_Dj[h, A_D_cpg[[j]][h] ] = 1
        }
        
        theta.GDj_A = c(theta.G_hat_AG, Theta.D_hat_AD[[j]]) 
        Sigma_G_hat_GDj.A = Sigma_G_hat[c(A_G_cpg, A_D_cpg_forunion[[j]]), c(A_G_cpg, A_D_cpg_forunion[[j]])]
        
        var_indrct_j1 = as.numeric( sigma2_Y * (t(betahat[(j+1),A_G])%*%E_A_G + t(betahat[(j+1),A_D[[j]]])%*%E_A_Dj ) %*% (Sigmahat.inv) %*% ( t(E_A_G)%*%betahat[(j+1),A_G] + t(E_A_Dj)%*%betahat[(j+1),A_D[[j]] ] ) ) +
          as.numeric((solve(J)[(j+1),(j+1)])* t(theta.GDj_A) %*% (Sigma_G_hat_GDj.A)%*% (theta.GDj_A))
        
        NIE_inference[(2*j),"Estimate"] = (beta.theta.G[(j+1),1] + beta.Theta.D[[j]][(j+1),1] )
        NIE_inference[(2*j),"SE"] = sqrt(var_indrct_j1)/sqrt(n)
        NIE_inference[(2*j),"test-stat"] = NIE_inference[(2*j),"Estimate"] / NIE_inference[(2*j),"SE"]
        NIE_inference[(2*j),"p value"] = 2*pnorm( -abs( NIE_inference[(2*j),"test-stat"] )  )
        
      }
      
    }
    

  }
  
  
  NDE_inference = matrix(NA, nrow = (q*2), ncol = 4)
  colnames(NDE_inference) = c("Estimate", "SE", "test-stat", "p value")
  rownames(NDE_inference) = rep("NDE", times = (q*2))
  
  
  for (j in 1:q) {
    rownames(NDE_inference)[(2*j-1)] = paste("NDE(T", j, "= 0)")
    rownames(NDE_inference)[(2*j)] = paste("NDE(T", j, "= 1)")
    
    if(length(A_D[[j]])==0){
      
      var_drct_j = as.numeric(sigma2_Y *( Sigmahat.inv[ (j+1),(j+1) ] ) )
      NDE_inference[(2*j-1), "Estimate"] = theta.T_hat[(j)]
      NDE_inference[(2*j-1), "SE"] = sqrt(var_drct_j)/sqrt(n)
      NDE_inference[(2*j-1), "test-stat"] = NDE_inference[(2*j-1), "Estimate"] / NDE_inference[(2*j-1), "SE"]
      NDE_inference[(2*j-1), "p value"] = 2*pnorm(-abs( NDE_inference[(2*j-1), "test-stat"] ))
      NDE_inference[(2*j),] = NDE_inference[(2*j-1),]
      
    } else {
      
      
      
      e_tilde.jp1 = matrix(0, nrow=nrow(Sigmahat), ncol=1)
      e_tilde.jp1[(j+1),1]=1
      
      E_A_Dj = matrix(0, nrow = length(A_D[[j]]), ncol = nrow(Sigmahat))
      colnames(E_A_Dj) = colnames(Sigmahat)
      for (h in 1:nrow(E_A_Dj)) {
        E_A_Dj[h, A_D_cpg[[j]][h]] = 1
      }
      
      alpha.hat.j = e_tilde.jp1 + t(E_A_Dj)%*%( betahat[1,A_D_cpg_forunion[[j]] ])
      
      
      e1 = matrix(0, nrow=nrow(J), ncol=1)
      e1[1,1]=1
      
      var_drct_j0 = as.numeric( sigma2_Y *( t(alpha.hat.j)%*%Sigmahat.inv %*%alpha.hat.j)) + 
        as.numeric(t(e1)%*%solve(J)%*%(e1))* as.numeric(t( Theta.D_hat_AD[[j]]  ) %*% (Sigma_G_hat[A_D_cpg_forunion[[j]], A_D_cpg_forunion[[j]] ])%*% ( Theta.D_hat_AD[[j]] ))
      
      NDE_inference[(2*j-1), "Estimate"] = theta.T_hat[j] + beta.Theta.D[[j]][1,1]
      NDE_inference[(2*j-1), "SE"] = sqrt(var_drct_j0)/sqrt(n)
      NDE_inference[(2*j-1), "test-stat"] = NDE_inference[(2*j-1), "Estimate"] / NDE_inference[(2*j-1), "SE"]
      NDE_inference[(2*j-1), "p value"] = 2*pnorm(-abs( NDE_inference[(2*j-1), "test-stat"] ))
      
      
      e.jp1 = matrix(0, nrow=nrow(J), ncol=1)
      e.jp1[(j+1),1]=1
      
      var_drct_j1 = as.numeric(sigma2_Y *( t( alpha.hat.j + t(E_A_Dj)%*%betahat[(j+1),A_D[[j]] ] )%*%Sigmahat.inv %*%( alpha.hat.j + t(E_A_Dj)%*%betahat[(j+1),A_D[[j]] ] ) )) + 
        as.numeric(t(e1 + e.jp1)%*%solve(J)%*%(e1 + e.jp1)) * as.numeric(t( Theta.D_hat_AD[[j]]  ) %*% (Sigma_G_hat[A_D_cpg_forunion[[j]], A_D_cpg_forunion[[j]] ])%*% ( Theta.D_hat_AD[[j]] ))
      
      NDE_inference[(2*j), "Estimate"] = theta.T_hat[j] + beta.Theta.D[[j]][1,1] + beta.Theta.D[[j]][(1+j),1]
      NDE_inference[(2*j), "SE"] = sqrt(var_drct_j1)/sqrt(n)
      NDE_inference[(2*j), "test-stat"] = NDE_inference[(2*j), "Estimate"] / NDE_inference[(2*j), "SE"]
      NDE_inference[(2*j), "p value"] = 2*pnorm(-abs( NDE_inference[(2*j), "test-stat"] ))
      
    }
    
  }
  

  
  return( list(outcome.theta_inference = Outcome.theta_inference,
               NIE_inference = NIE_inference,
               NDE_inference = NDE_inference,
               W_A = W_A, A_G = A_G, A_D = A_D)
  )
  
}
























































##### The function to calculate NIE SE's (population version) with known beta, theta, and support ###
### This function is only used in theoretical simulation ###
### Note that the oracle method has the same NIE, individual theta, SE but the point estimation is calculated by two low-dim OLSs rather than SCAD + post selected OLS ###
### Note that this SE is the SE for \sqrt{n} (effecthat - effect)
### This function does NOT the covariate X
OraProjInt_SE_popu_noX = function(theta0, theta.T, theta.G, Theta.D,
                                beta0, B, Sigma_G, sigma_Y,
                                prop){
  
  #' @param theta0 a 1 by 1 matrix 
  #' @param theta.T a matrix with nrow = q, ncol = 1
  #' @param theta.G a matrix with nrow = p, ncol =  1
  #' @param Theta.D a matrix with nrow = p, ncol = q
  #' @param beta0 a matrix with nrow = p, ncol = 1
  #' @param B a matrix with nrow = p, ncol = q
  #' @param Sigma_G a p by p matrix, covariance matrix of bepsilon_G
  #' @param sigma_Y, a scalar (vector with length 1) which is the standard deviation of epsilon_Yi
  #' @param prop is a vector with length q, whose entries are the p_j = lim( n_j/n ), j=1,...,q, the proportion of the j-th treatment. Note that sum(prop) < 1 since we have p_0 not included.
  #' The entries in theta0, theta.T, theta.G, Theta.D, beta0, B, Sigma_G have no rownames or colnames
  #' The theta's can have different nonzero entry locations, not limited to heredity
  #' @return NIE_SE is the matrix of SE of the NIE and NDE for different j
  
  
  q = ncol(Theta.D)
  
  # First find the support of theta.G and Theta.D
  # which function will first transform the argument theta.G or Theta.D into a vector connected by columns
  # Then give an vector output with nonzero locations in numbers
  A_G = which(theta.G != 0) # A_G is a vector, containing number location in 1:p with length s.G
  A_D = list() # A_D has q components. Each components contains number location in 1:p with length s.Dj
  for (j in 1:q) {
    A_D[[j]] = which(Theta.D[,j]!=0)
  }
  
  s.G = length(A_G)
  
  s.D = c()
  for (j in 1:q) {
    s.D[j] = length(A_D[[j]])
  }
  
  Sigma11 = matrix(0, nrow = (1+q), ncol = (1+q))
  Sigma11[1,] = c(1, prop)
  for (j in 1:q) {
    Sigma11[(j+1), 1] = prop[j]
    Sigma11[(j+1), (j+1)] = prop[j]
  }
  
  prop0 = 1 -  sum(prop)
  
 
  Sigma13 = matrix(0, nrow = (1+q), ncol = s.G)
  Sigma13[1, ] = prop0 * t( beta0[A_G,1] )
  for (j in 1:q) {
    Sigma13[1, ] = Sigma13[1, ] + prop[j] * t( beta0[A_G,1] + B[A_G,j]  )
  }
  for (j in 1:q) {
    Sigma13[(j+1), ] = prop[j] * t( beta0[A_G,1] + B[A_G,j]  )
  }
  
  Sigma14 = matrix(0, nrow = (1+q), ncol = sum(s.D))
  if(q == 1){
    Sigma14[1, 1:s.D[1]] = prop[1] * t( beta0[A_D[[1]], 1] + B[ A_D[[1]] , 1] )
    Sigma14[2, 1:s.D[1]] = prop[1] * t( beta0[A_D[[1]], 1] + B[ A_D[[1]] , 1] )
  } else {
    Sigma14[1, 1:s.D[1]] = prop[1] * t( beta0[A_D[[1]], 1] + B[ A_D[[1]] , 1] )
    Sigma14[2, 1:s.D[1]] = prop[1] * t( beta0[A_D[[1]], 1] + B[ A_D[[1]] , 1] )
    for (j in 2:q) {
      Sigma14[1, (sum(s.D[1:(j-1)])+1) : sum(s.D[1:j]) ] = prop[j] * t( beta0[A_D[[j]], 1] + B[ A_D[[j]] , j] ) # just because it has 1:(j-1), if j=1, 1:(j-1) is meaningless, so I have to separate j=1 and j>=2. Actually, I do not need such procedure. I can first define j=1 part. Then add if (q>1){ loop the j>=2 part  } like Sigma34. But I have already done it here so I won't change.
      Sigma14[(j+1), (sum(s.D[1:(j-1)])+1) : sum(s.D[1:j]) ] = prop[j] * t( beta0[A_D[[j]], 1] + B[ A_D[[j]] , j] )
    }
  }
  
  
  Sigma31 = t(Sigma13)
  
  Sigma33 = prop0 * ( ( beta0[A_G,1] )%*%t( beta0[A_G,1] ) + Sigma_G[A_G, A_G]  )
  for (j in 1:q) {
    Sigma33 =  Sigma33 + prop[j] * ( ( beta0[A_G,1] + B[A_G,j] )%*%t( beta0[A_G,1] + B[A_G,j] ) + Sigma_G[A_G, A_G] )
  }
  
  
  Sigma34 = matrix(0, nrow = s.G, ncol=sum(s.D) )
  Sigma34[ , 1:s.D[1]] = prop[1] *  ( ( beta0[ A_G, 1] + B[A_G, 1] )%*%t( beta0[ A_D[[1]], 1] + B[A_D[[1]],1] ) + Sigma_G[ A_G, A_D[[1]] ]  ) 
  if(q > 1){
    for (j in 2:q) {
      Sigma34[, (sum(s.D[1:(j-1)])+1) : sum(s.D[1:j])] = prop[j] *  ( ( beta0[ A_G, 1] + B[A_G, j] )%*%t( beta0[ A_D[[j]], 1] + B[A_D[[j]],j] ) + Sigma_G[ A_G, A_D[[j]] ]  ) 
    }
  }
  
  
  Sigma41 = t(Sigma14)
  
  Sigma43 = t(Sigma34)
  
  Sigma44 = matrix(0, nrow = sum(s.D), ncol = sum(s.D)  )
  Sigma44[ 1:s.D[1], 1:s.D[1] ] = prop[1] *  ( ( beta0[ A_D[[1]], 1] + B[A_D[[1]], 1] )%*%t( beta0[ A_D[[1]], 1] + B[A_D[[1]],1] ) + Sigma_G[ A_D[[1]], A_D[[1]] ]  ) 
  if(q > 1){
    for (j in 2:q) {
      Sigma44[ (sum(s.D[1:(j-1)])+1) : sum(s.D[1:j]) , (sum(s.D[1:(j-1)])+1) : sum(s.D[1:j]) ] = prop[j] *  ( ( beta0[ A_D[[j]], 1] + B[A_D[[j]], j] )%*%t( beta0[ A_D[[j]], 1] + B[A_D[[j]],j] ) + Sigma_G[ A_D[[j]], A_D[[j]] ]  ) 
    }
  }
  
  Sigma = rbind(cbind(Sigma11, Sigma13, Sigma14),
                cbind(Sigma31, Sigma33, Sigma34),
                cbind(Sigma41, Sigma43, Sigma44))
  
  Sigma.inv=solve(Sigma)
  
  Sigma11.inv = solve(Sigma11)
  
  
  NIE_SE = matrix(NA, nrow = q, ncol = 2)
  colnames(NIE_SE) = c("Tj=0", "Tj=1")
  rownames(NIE_SE) = paste0("j=", 1:q)
  
  E_AG = matrix(0, nrow = s.G, ncol = ncol(Sigma))
  for (h in 1:s.G) {
    E_AG[h, 1+q+h]=1
  }
  
  E_AD1 = matrix(0, nrow = s.D[1], ncol = ncol(Sigma))
  for (h in 1:s.D[1]) {
    E_AD1[h, 1+q+s.G+h] = 1
  }
  
  var.NIE0 = (sigma_Y^2)* ( t(B[A_G,1])%*% E_AG %*% (Sigma.inv) %*% t(E_AG)%*%B[A_G,1] ) + (Sigma11.inv[2,2]) * (t(theta.G[A_G,1])%*%(Sigma_G[A_G,A_G])%*%theta.G[A_G,1])
  
  NIE_SE[1,"Tj=0"] = sqrt(var.NIE0)
  rm(var.NIE0)
  
  var.NIE1 = (sigma_Y^2)* ( ( t(B[A_G,1])%*%E_AG + t(B[A_D[[1]],1])%*%E_AD1 )%*%Sigma.inv%*%( t(E_AG)%*%B[A_G,1] + t(E_AD1)%*%B[A_D[[1]],1] ) ) + 
    (Sigma11.inv[2,2]) *  ( t(c(theta.G[A_G,1], Theta.D[A_D[[1]],1]))%*% ( rbind( cbind(Sigma_G[A_G,A_G], Sigma_G[A_G,A_D[[1]]]),
                                                                            cbind(Sigma_G[A_D[[1]],A_G], Sigma_G[A_D[[1]],A_D[[1]]])) )%*%c(theta.G[A_G,1], Theta.D[A_D[[1]],1]) )
  
  NIE_SE[1, "Tj=1"] = sqrt(var.NIE1)
  rm(var.NIE1)
  
  if(q > 1){
    for (j in 2:q) {
      
      E_ADj = matrix(0, nrow = s.D[j], ncol = ncol(Sigma))
      for (h in 1:s.D[j]) {
        E_ADj[h, ( 1+q+s.G + sum(s.D[1:(j-1)]) + h ) ] = 1
      }
      
      var.NIE0 = (sigma_Y^2)* ( t(B[A_G,j])%*% E_AG %*% (Sigma.inv) %*% t(E_AG)%*%B[A_G,j]) + (Sigma11.inv[(j+1),(j+1)]) * (t(theta.G[A_G,1])%*%(Sigma_G[A_G,A_G])%*%theta.G[A_G,1])
      
      NIE_SE[j,"Tj=0"] = sqrt(var.NIE0)
      rm(var.NIE0)
      
      var.NIE1 = (sigma_Y^2)* ( ( t(B[A_G,j])%*%E_AG + t(B[A_D[[j]],j])%*%E_ADj )%*%Sigma.inv%*%( t(E_AG)%*%B[A_G,j] + t(E_ADj)%*%B[A_D[[j]],j] ) ) + 
        (Sigma11.inv[(j+1),(j+1)]) *  ( t(c(theta.G[A_G,1], Theta.D[A_D[[j]],j]))%*% ( rbind( cbind(Sigma_G[A_G,A_G], Sigma_G[A_G,A_D[[j]]]),
                                                                                        cbind(Sigma_G[A_D[[j]],A_G], Sigma_G[A_D[[j]],A_D[[j]]])) )%*%c(theta.G[A_G,1], Theta.D[A_D[[j]],j]) )
      
      NIE_SE[j, "Tj=1"] = sqrt(var.NIE1)
      rm(var.NIE1)
      
    }
  }
  
  return( list(NIE_SE = NIE_SE)
  )
  
}















### The function to calculate Runze NIE, NDE point estimates and SE's (sample version) ### 
## this function does not have the covariate X ##
Runze_SE_sample_noX = function(Y, Trt, Mediator, 
                              theta0_hat, theta.T_hat, theta.G_hat){
  
  #' @param Y vector, n-dimensional outcome vector (with column name such as "Y").
  #' @param Trt matrix or dataframe, n by q treatment matrix without intercept (with column name such as "T1", "T2",.., "Tq").
  #' @param Mediator matrix or dataframe, n by p mediator matrix (with column name, "G1",...).
  #' @param theta0_hat a vector with length 1 from SCAD estimation
  #' @param theta.T_hat a vector with length q from SCAD estimation
  #' @param theta.G_hat a vector of length p from SCAD estimation
  #' Note that Trt, Mediator, Trt.Medator must have column names in words!
  #' theta0_hat, theta.T_hat, theta.G_hat are vectors with no names! Just pure vector
  #' If you use ncvreg packages, it will produce names of theta0_hat, theta.T_hat...
  #' To remove those names of theta, use the code before plugging them in this function such as: names(theta0.hat)=NULL; names(theta.T.hat)=NULL; names(theta.X.hat)=NULL; names(theta.G.hat)=NULL; names(vtheta.D.hat)=NULL
  #' @return NIE_est, the point estimator of beta0^T theta_G, beta1^T theta_G,..., beta_q^T theta_G. Its (j+1)th entry is for betaj^T theta_G, j=1,...,q.
  #' @return NIE_SE, see below:
  #' the SE of NIE_est, it has already been divided by sqrt(n)!
  #' Its (j+1) -th entry is the SE for beta_j^T theta_G (because we see intercept one.n as a column of exposure in Runze's paper)
  #' @return NDE_est, the point estimator of theta.T1,...,theta.Tq. Its (j+1)th entry is for theta.Tj, j=1...q.
  #' @return NDE_SE, the SE, it has already been divided by sqrt(n). Its (j+1)th entry is for theta.Tj, j=1...q.
 
  
  
  p = ncol(Mediator)
  n = nrow(Mediator)
  q = ncol(Trt)
  
  
  A_G = which(theta.G_hat!=0) # A_G contains number location in 1:p with length \hat{s_G}
  A_G_cpg = colnames(Mediator)[A_G] # A_G contains word names of those support mediators with length \hat{s.G}
  
  one.n = matrix(1, nrow=nrow(Mediator), ncol=1)
  
  W_A = cbind(one.n, Trt, Mediator[,A_G])
  
  
  tmpname = c()
  tmpname[1] = "one"
  tmpname[(1+1):(1+q)] = colnames(Trt)
  if(length(A_G) > 0){
    tmpname[(1+q+1):(1+q+length(A_G))] = A_G_cpg
  } 
  colnames(W_A) = tmpname
  rm(tmpname)
  
  
  Sigmahat = (t(W_A)%*%W_A) / nrow(Mediator)
  Sigmahat.inv = solve(Sigmahat)
  
  eps_Y = Y - one.n%*%theta0_hat - Trt%*%theta.T_hat - Mediator%*%theta.G_hat 
  sigma2_Y = as.numeric( t(eps_Y)%*%eps_Y/n ) # the estimate of sigma_Y^2
  rm(eps_Y)
  
  
  
  ###
  H = cbind(one.n, Trt)
  colnames(H)[1] = "one"
  
  J = (t(H)%*%H)/n
  J.inv = solve(J)
  
  total_runze = solve(t(H)%*%H)%*%t(H)%*%Y
  
  NIE.est = total_runze - as.matrix(c(theta0_hat, theta.T_hat))
  
  eps_3 = Y - H%*%total_runze
  sigma2_H = as.numeric( t(eps_3)%*%eps_3/n )
  
  sigma2_2 = c()
  if(sigma2_H > sigma2_Y){
    sigma2_2 = sigma2_H - sigma2_Y
  } else {
    sigma2_2 = 0
  }
  
  # the covariance matrix of NIE in eqn (2.8) of Runze JoE
  COV.NIE = (1/(n)) * ( sigma2_Y * (Sigmahat.inv[1:(1+q), 1:(1+q)] - J.inv ) + sigma2_2 * J.inv  )
  NIE.se = as.matrix(sqrt(diag(COV.NIE)))
  
  NDE.est = as.matrix(c(theta0_hat, theta.T_hat))
  rownames(NDE.est) = c("one", colnames(Trt)  )
  
  COV.NDE = (1/(n)) * ( sigma2_Y * ( Sigmahat.inv[1:(1+q), 1:(1+q)] )  )
  NDE.se = as.matrix(sqrt(diag(COV.NDE)))
  
 
  return( list(NIE_est = NIE.est, NIE_SE = NIE.se, 
               NDE_est = NDE.est, NDE_SE = NDE.se) )
  
}








































### Dave
library(scalreg)
library(flare)

# This function only output Dave's debiased result
Dave = function (Y, G, S, lam_list = NA, center.SG = TRUE, center.Y = TRUE)  
{ 
  
  #' @param Y n by one vector of outcome
  #' @param G n by p matrix of mediators
  #' @param S n by q matrix of exposures, q could > 1
  #' Note that mediation_setting is always "incomplete" i.e. Treatment or exposure is presented in the outcome model.
  #' @param lam_list is a single number, the tau_n used in Omegahat, usually is just sqrt(log(p)/n)/3
  #' @param center.SG: Whether center S and G or not
  #' @param center.Y: Whether center Y or not
  #' @return infer_out:  inference the matrix for the inference of indirect effect beta_0 and direct effect alpha_1
  
  n = dim(G)[1]
  p = dim(G)[2]
  q = dim(S)[2]
  
  if (center.Y == TRUE) {
    Y = Y - mean(Y)
  }
  
  
  if (center.SG == TRUE) {
    sm = matrix(rep(colMeans(S), n), nrow = n, ncol = q, 
                byrow = T)
    S = S - sm
    gm = matrix(rep(colMeans(G), n), nrow = n, ncol = p, 
                byrow = T)
    G = G - gm
  }
  
  X = cbind(G, S)
  sigma_SS_hat = t(S) %*% S/n
  Sigma_SG_hat = t(S) %*% G/n
  Sigma_GG_hat = t(G) %*% G/n
  Sigma_XX_hat = t(X) %*% X/n
  sigma_SS_hat_inverse = solve(sigma_SS_hat)
  result_scalreg = scalreg(X, Y)
  alpha_hat = result_scalreg$co # it does not contain the intercept, length(alpha_hat) = ncol(X)
  sigma1_hat = result_scalreg$hsigma
  if (sigma1_hat > 10) {
    result_scalreg = scalreg(round(X, 2), round(Y, 2))
    alpha_hat = result_scalreg$co
    sigma1_hat = result_scalreg$hsigma
  }
  
  sigma_hat = summary(lm(Y ~ S))$sigma
  sigma2_hat = sqrt(max(sigma_hat^2 - sigma1_hat^2, 0))
  lambda.k_hat = t(X) %*% (Y - X %*% alpha_hat)/n
  Dhat = matrix(0, 2 * q, p + q)
  Dhat[1:q, 1:p] = Sigma_SG_hat
  Dhat[(q + 1):(2 * q), ((p + 1):(p + q))] = sigma_SS_hat
  Beta.list = list()
  for (qj in 1:(2 * q)) {
    outj = slim(X, t(n * Dhat), Y_rachel_column = qj, 
                method = "dantzig", lambda = lam_list)
    Beta.list[[qj]] = outj$beta
  }
  
  
  Omega_hat = matrix(0, 2 * q, p + q)
  for (qj in 1:(2 * q)) {
    Omega_hat[qj, ] = Beta.list[[qj]][, 1]
  }
  
  
  # eqn (6)
  beta_part1 = kronecker(diag(2), sigma_SS_hat_inverse) %*% 
    Omega_hat %*% lambda.k_hat
  
  # \hat{b} in Theorem 2 (indirect effect estimator)
  beta_hat = beta_part1[1:q, 1] + sigma_SS_hat_inverse %*% 
    Sigma_SG_hat %*% alpha_hat[1:p]
  # \hat{a} in Theorem 2 (direct effect estimator)
  alpha1_hat = beta_part1[(q + 1):(2 * q), 1] + alpha_hat[(p+1):(p+q)] # I changed the original function alpha_hat[p+1] here
  
  
  total_hat = alpha1_hat + beta_hat
  
  # the asymptotic co-variance matrix below Theorem 2
  Cov_part1 = sigma1_hat^2 * kronecker(diag(2), sigma_SS_hat_inverse) %*% 
    Omega_hat %*% Sigma_XX_hat %*% t(Omega_hat) %*% kronecker(diag(2), 
                                                              sigma_SS_hat_inverse)
  Sigma_11 = Cov_part1[1:q, 1:q] + sigma2_hat^2 * sigma_SS_hat_inverse
  Sigma_12 = Cov_part1[(q + 1):(2 * q), 1:q]
  Sigma_22 = Cov_part1[(q + 1):(2 * q), (q + 1):(2 * q)]
  
  sigma_beta_hat = Sigma_11 # covariance matrix for \sqrt{n}(\hat{b} - beta_0) in Theorem 2
  sigma_alpha1_hat = Sigma_22 # covariance matrix for \sqrt{n}(\hat{a} - alpha_1) in Theorem 2
  
  infer_out = list()
  infer_out$tau_n = lam_list
  infer_out$inference = matrix(NA, nrow=2*q, ncol = 4)
  rownames(infer_out$inference) = c(paste0("beta_0,",1:q), paste0("alpha_1,",1:q))
  colnames(infer_out$inference) = c("est", "SE", "test stat", "p value")
  
  for (j in 1:q) {
    infer_out$inference[j, "est"] = beta_hat[j]
    infer_out$inference[j, "SE"] = sqrt( (sigma_beta_hat[j,j])/n )
    if(infer_out$inference[j, "SE"] != 0){
      infer_out$inference[j, "test stat"] = beta_hat[j] / sqrt( (sigma_beta_hat[j,j])/n )
      infer_out$inference[j, "p value"] = 2 * (1 - pnorm(abs( infer_out$inference[j, "test stat"] )))
    }
  }
  
  if (q ==1 ) {
    infer_out$inference[q+1, "est"] = alpha1_hat
    infer_out$inference[q+1, "SE"] = sqrt( (sigma_alpha1_hat)/n )
    if(infer_out$inference[q+1, "SE"]!=0){
      infer_out$inference[q+1, "test stat"] = alpha1_hat / sqrt( (sigma_alpha1_hat)/n )
      infer_out$inference[q+1, "p value"] = 2 * (1 - pnorm(abs( infer_out$inference[q+1, "test stat"] )))
    }
  } else {
    for (j in 1:q) {
      infer_out$inference[q+j, "est"] = alpha1_hat[j]
      infer_out$inference[q+j, "SE"] = sqrt( (sigma_alpha1_hat[j,j])/n )
      if(infer_out$inference[q+j, "SE"] != 0){
        infer_out$inference[q+j, "test stat"] = alpha1_hat[j] / sqrt( (sigma_alpha1_hat[j,j])/n )
        infer_out$inference[q+j, "p value"] = 2 * (1 - pnorm(abs( infer_out$inference[q+j, "test stat"] )))
      }
    }
  }
  return(infer_out)
}

# This function will output both Dave's debiased result and naive non-debiased (scaled lasso) point estimate
Dave1 = function (Y, G, S, lam_list = NA, center.SG = TRUE, center.Y = TRUE)  
{ 
  
  #' @param Y n by one vector of outcome
  #' @param G n by p matrix of mediators
  #' @param S n by q matrix of exposures, q could > 1
  #' Note that mediation_setting is always "incomplete" i.e. Treatment or exposure is presented in the outcome model.
  #' @param lam_list is a single number, the tau_n used in Omegahat, usually is just sqrt(log(p)/n)/3
  #' @param center.SG: Whether center S and G or not
  #' @param center.Y: Whether center Y or not
  #' @return debiased_inference:  inference the matrix for the inference of indirect effect beta_0 and direct effect alpha_1
  #' @return non_debiased_NIE_est: point estimate of NIE without Omega part in eqn(6) of Dave's paper
  #' @return non_debiased_NDE_est: point estimate of NDE without Omega part in eqn(6) of Dave's paper
  
  n = dim(G)[1]
  p = dim(G)[2]
  q = dim(S)[2]
  
  if (center.Y == TRUE) {
    Y = Y - mean(Y)
  }
  
  
  if (center.SG == TRUE) {
    sm = matrix(rep(colMeans(S), n), nrow = n, ncol = q, 
                byrow = T)
    S = S - sm
    gm = matrix(rep(colMeans(G), n), nrow = n, ncol = p, 
                byrow = T)
    G = G - gm
  }
  
  X = cbind(G, S)
  sigma_SS_hat = t(S) %*% S/n
  Sigma_SG_hat = t(S) %*% G/n
  Sigma_GG_hat = t(G) %*% G/n
  Sigma_XX_hat = t(X) %*% X/n
  sigma_SS_hat_inverse = solve(sigma_SS_hat)
  result_scalreg = scalreg(X, Y)
  alpha_hat = result_scalreg$co # it does not contain the intercept, length(alpha_hat) = ncol(X)
  sigma1_hat = result_scalreg$hsigma
  if (sigma1_hat > 10) {
    result_scalreg = scalreg(round(X, 2), round(Y, 2))
    alpha_hat = result_scalreg$co
    sigma1_hat = result_scalreg$hsigma
  }
  
  sigma_hat = summary(lm(Y ~ S))$sigma
  sigma2_hat = sqrt(max(sigma_hat^2 - sigma1_hat^2, 0))
  lambda.k_hat = t(X) %*% (Y - X %*% alpha_hat)/n
  Dhat = matrix(0, 2 * q, p + q)
  Dhat[1:q, 1:p] = Sigma_SG_hat
  Dhat[(q + 1):(2 * q), ((p + 1):(p + q))] = sigma_SS_hat
  Beta.list = list()
  for (qj in 1:(2 * q)) {
    outj = slim(X, t(n * Dhat), Y_rachel_column = qj, 
                method = "dantzig", lambda = lam_list)
    Beta.list[[qj]] = outj$beta
  }
  
  
  Omega_hat = matrix(0, 2 * q, p + q)
  for (qj in 1:(2 * q)) {
    Omega_hat[qj, ] = Beta.list[[qj]][, 1]
  }
  
  
  # eqn (6)
  beta_part1 = kronecker(diag(2), sigma_SS_hat_inverse) %*% 
    Omega_hat %*% lambda.k_hat
  
  # \hat{b} in Theorem 2 (indirect effect estimator)
  beta_hat = beta_part1[1:q, 1] + sigma_SS_hat_inverse %*% 
    Sigma_SG_hat %*% alpha_hat[1:p]
  # \hat{a} in Theorem 2 (direct effect estimator)
  alpha1_hat = beta_part1[(q + 1):(2 * q), 1] + alpha_hat[(p+1):(p+q)] # I changed the original function alpha_hat[p+1] here
  
  
  total_hat = alpha1_hat + beta_hat
  
  # the asymptotic co-variance matrix below Theorem 2
  Cov_part1 = sigma1_hat^2 * kronecker(diag(2), sigma_SS_hat_inverse) %*% 
    Omega_hat %*% Sigma_XX_hat %*% t(Omega_hat) %*% kronecker(diag(2), 
                                                              sigma_SS_hat_inverse)
  Sigma_11 = Cov_part1[1:q, 1:q] + sigma2_hat^2 * sigma_SS_hat_inverse
  Sigma_12 = Cov_part1[(q + 1):(2 * q), 1:q]
  Sigma_22 = Cov_part1[(q + 1):(2 * q), (q + 1):(2 * q)]
  
  sigma_beta_hat = Sigma_11 # covariance matrix for \sqrt{n}(\hat{b} - beta_0) in Theorem 2
  sigma_alpha1_hat = Sigma_22 # covariance matrix for \sqrt{n}(\hat{a} - alpha_1) in Theorem 2
  
  infer_out = list()
  infer_out$tau_n = lam_list
  infer_out$inference = matrix(NA, nrow=2*q, ncol = 4)
  rownames(infer_out$inference) = c(paste0("beta_0,",1:q), paste0("alpha_1,",1:q))
  colnames(infer_out$inference) = c("est", "SE", "test stat", "p value")
  
  for (j in 1:q) {
    infer_out$inference[j, "est"] = beta_hat[j]
    infer_out$inference[j, "SE"] = sqrt( (sigma_beta_hat[j,j])/n )
    if(infer_out$inference[j, "SE"] != 0){
      infer_out$inference[j, "test stat"] = beta_hat[j] / sqrt( (sigma_beta_hat[j,j])/n )
      infer_out$inference[j, "p value"] = 2 * (1 - pnorm(abs( infer_out$inference[j, "test stat"] )))
    }
  }
  
  if (q ==1 ) {
    infer_out$inference[q+1, "est"] = alpha1_hat
    infer_out$inference[q+1, "SE"] = sqrt( (sigma_alpha1_hat)/n )
    if(infer_out$inference[q+1, "SE"]!=0){
      infer_out$inference[q+1, "test stat"] = alpha1_hat / sqrt( (sigma_alpha1_hat)/n )
      infer_out$inference[q+1, "p value"] = 2 * (1 - pnorm(abs( infer_out$inference[q+1, "test stat"] )))
    }
  } else {
    for (j in 1:q) {
      infer_out$inference[q+j, "est"] = alpha1_hat[j]
      infer_out$inference[q+j, "SE"] = sqrt( (sigma_alpha1_hat[j,j])/n )
      if(infer_out$inference[q+j, "SE"] != 0){
        infer_out$inference[q+j, "test stat"] = alpha1_hat[j] / sqrt( (sigma_alpha1_hat[j,j])/n )
        infer_out$inference[q+j, "p value"] = 2 * (1 - pnorm(abs( infer_out$inference[q+j, "test stat"] )))
      }
    }
  }
  
  non.debiased.nie = sigma_SS_hat_inverse %*% Sigma_SG_hat %*% alpha_hat[1:p]
  non.debiased.nde = alpha_hat[(p+1):(p+q)]
  
  return( list(debiased_inference = infer_out,
              non_debiased_NIE_est = non.debiased.nie,
              non_debiased_NDE_est = non.debiased.nde) )
}




### slim function in Dave
#----------------------------------------------------------------------------------#
# Package: flare                                                                   #
# flare.slim(): The user interface for slim()                                      #
# Author: Xingguo Li                                                               #
# Email: <xingguo.leo@gmail.com>                                                   #
# Date: Mar 16th 2014                                                              #
# Version: 1.2.0                                                                   #
#----------------------------------------------------------------------------------#

slim <- function(X,
                 Y_rachel,
                 Y_rachel_column = 1,
                 lambda = NULL,
                 nlambda = NULL,
                 lambda.min.value = NULL,
                 lambda.min.ratio = NULL,
                 rho = 1,
                 method="dantzig",
                 q = 2,
                 res.sd = FALSE,
                 prec = 1e-5,
                 max.ite = 1e5,
                 verbose = TRUE)
{
  if(method!="dantzig" && method!="lq" && method!="lasso"){
    cat("\"method\" must be dantzig, lasso or lq.\n")
    return(NULL)
  }
  if(method=="lq"){
    if(q<1 || q>2){
      cat("q must be in [1, 2] when method = \"lq\".\n")
      return(NULL)
    }
  }
  if(verbose) {
    cat("Sparse Linear Regression with L1 Regularization.\n")
  }
  
  
  
  n = nrow(X)  #*****************************************
  d = ncol(X)  #****************************************
  
  if(n==0 || d==0) {
    cat("No data input.\n")
    return(NULL)
  }
  
  #rachel read below***************************
  maxdf = max(n,d)
  xm=matrix(rep(colMeans(X),n),nrow=n,ncol=d,byrow=T)
  x1=X-xm
  sdxinv=1/sqrt(colSums(x1^2)/(n-1))
  xx=x1*matrix(rep(sdxinv,n),nrow=n,ncol=d,byrow=T)
  ym=NA
  
  #y1=Y-ym
  Y_rachel_ncol = dim(Y_rachel)[2]
  y1 = Y_rachel * matrix(rep(sdxinv,Y_rachel_ncol),nrow=d,ncol=Y_rachel_ncol)
  
  #rachel read above***************************
  if(res.sd == TRUE){
    sdy=sqrt(sum(y1^2)/(n-1))
    yy=y1/sdy
  }else{
    sdy = 1
    yy = y1
  }
  intercept = FALSE
  
  if(intercept){
    xx = cbind(rep(1, nrow(xx)), xx)
    X = cbind(rep(1, nrow(X)), X)
    d = d+1
  }
  
  if(!is.null(lambda)) nlambda = length(lambda)
  if(is.null(lambda)){
    if(is.null(nlambda))
      nlambda = 5
    if(method=="dantzig"){
      if(intercept)
        lambda.max = max(abs(crossprod(xx[,2:d],yy/n)))
      else
        #lambda.max = max(abs(crossprod(xx,yy/n)))
        lambda.max = max(abs(yy/n))
    }
    
    if(method=="dantzig"){
      if(is.null(lambda.min.ratio)){
        lambda.min.ratio = 0.5
      }
      if(is.null(lambda.min.value)){
        lambda.min.value = lambda.min.ratio*lambda.max
      }
    }else{
      if(is.null(lambda.min.value)){
        lambda.min.value = sqrt(log(d)/n)
      }else{
        if(is.null(lambda.min.ratio)){
          lambda.min.ratio = lambda.min.value/lambda.max
        }
      }
    }
    if(lambda.max<lambda.min.value){
      lambda.max = 1
      lambda.min.value = 0.4
    }
    lambda = exp(seq(log(lambda.max), log(lambda.min.value), length = nlambda))
    rm(lambda.max,lambda.min.value,lambda.min.ratio)
    gc()
  }
  if(is.null(rho))
    rho = 1
  begt=Sys.time()
  
  yy_run = as.matrix(yy[,Y_rachel_column])
  
  if(method=="dantzig"){ # dantzig
    if(d>=n)
      
      out = slim.dantzig.ladm.scr(yy_run, xx, lambda, nlambda, n, d, maxdf, rho, max.ite, prec, intercept, verbose)
    else
      out = slim.dantzig.ladm.scr2(yy_run, xx, lambda, nlambda, n, d, maxdf, rho, max.ite, prec, intercept, verbose)
    q = "infty"
  }
  
  runt=Sys.time()-begt
  
  df=rep(0,nlambda)
  if(intercept){
    for(i in 1:nlambda)
      df[i] = sum(out$beta[[i]][2:d]!=0)
  }else{
    for(i in 1:nlambda)
      df[i] = sum(out$beta[[i]]!=0)
  }
  
  est = list()
  intcpt0=matrix(0,nrow=1,ncol=nlambda)
  intcpt=matrix(0,nrow=1,ncol=nlambda)
  if(intercept){
    beta1=matrix(0,nrow=d-1,ncol=nlambda)
    for(k in 1:nlambda){
      tmp.beta = out$beta[[k]][2:d]
      beta1[,k]=sdxinv*tmp.beta*sdy
      intcpt[k] = ym-as.numeric(xm[1,]%*%beta1[,k])+out$beta[[k]][1]*sdy
      intcpt0[k] = intcpt[k]
    }
  }else{
    beta1=matrix(0,nrow=d,ncol=nlambda)
    for(k in 1:nlambda){
      tmp.beta = out$beta[[k]]
      intcpt0[k] = 0
      beta1[,k] = sdxinv*tmp.beta*sdy
      intcpt[k] = ym-as.numeric(xm[1,]%*%beta1[,k])
    }
  }
  
  est$beta0 = out$beta
  est$beta = beta1
  est$intercept = intcpt
  est$intercept0 = intcpt0
  est$Y = NA
  est$X = X
  est$lambda = lambda
  est$nlambda = nlambda
  est$df = df
  est$method = method
  est$q = q
  est$ite =out$ite
  est$verbose = verbose
  est$runtime = runt
  class(est) = "slim"
  #if(verbose) print(est)
  return(est)
}


#----------------------------------------------------------------------------------#
# Package: flare                                                                   #
# slim.dantzig.ladm.scr(): Regression with Dantzig()                               #
# Author: Xingguo Li                                                               #
# Email: <xingguo.leo@gmail.com>                                                   #
# Date: Jul 8th, 2014                                                              #
# Version: 1.4.0                                                                   #
#----------------------------------------------------------------------------------#

slim.dantzig.ladm.scr <- function(Y, X, lambda, nlambda, n, d, maxdf, rho, max.ite, prec, intercept, verbose)
{
  if(verbose==TRUE)
    cat("Dantzig selector with screening.\n")
  #XY = crossprod(X,Y/n)
  XY = Y/n
  XX = crossprod(X,X/n)#+0.05*lambda[nlambda]*diag(d)
  beta = matrix(0,nrow=d,ncol=nlambda)
  ite.int = rep(0,nlambda)
  ite.int1 = rep(0,nlambda)
  ite.int2 = rep(0,nlambda)
  if(intercept) {
    intcep=1
  }else{
    intcep=0
  }
  
  if(n<=3){
    num.scr1 = n
    num.scr2 = n
  }else{
    num.scr1 = ceiling(n/log(n))
    num.scr2 = n-1
  }
  
  order0 = order(abs(XY),decreasing = TRUE)
  idx.scr = order0; num.scr = length(idx.scr)
  idx.scr1 = order0[1:num.scr1]
  idx.scr2 = order0[1:num.scr2]
  X1 = X[,idx.scr]
  XXX = crossprod(X1,crossprod(tcrossprod(X1,X1),X1))/(n^2)
  gamma = max(colSums(abs(XXX)))
  str=.C("slim_dantzig_ladm_scr", as.double(XY), as.double(XX), as.double(XXX),
         as.double(beta), as.integer(n), as.integer(d), as.double(rho),
         as.integer(ite.int), as.integer(ite.int1), as.integer(ite.int2),
         as.integer(num.scr1), as.integer(num.scr2),
         as.integer(idx.scr), as.integer(idx.scr1), as.integer(idx.scr2),
         as.double(gamma), as.double(lambda),
         as.integer(nlambda), as.integer(max.ite), as.double(prec),
         #as.integer(intcep) )
         as.integer(intcep),PACKAGE="freebird")
  beta.list = vector("list", nlambda)
  for(i in 1:nlambda){
    beta.i = unlist(str[4])[((i-1)*d+1):(i*d)]
    beta.list[[i]] = beta.i
  }
  ite.int = unlist(str[8])
  ite.int1 = unlist(str[9])
  ite.int2 = unlist(str[10])
  ite = list()
  ite[[1]] = ite.int1
  ite[[2]] = ite.int2
  ite[[3]] = ite.int
  return(list(beta=beta.list, ite=ite))
}

#----------------------------------------------------------------------------------#
# Package: flare                                                                   #
# slim.dantzig.ladm.scr2(): Regression with Dantzig()                              #
# Author: Xingguo Li                                                               #
# Email: <xingguo.leo@gmail.com>                                                   #
# Date: Jul 8th, 2014                                                              #
# Version: 1.4.0                                                                   #
#----------------------------------------------------------------------------------#

slim.dantzig.ladm.scr2 <- function(Y, X, lambda, nlambda, n, d, maxdf, rho, max.ite, prec, intercept, verbose)
{
  if(verbose==TRUE)
    cat("Dantzig selector with screening.\n")
  #XY = crossprod(X,Y)/n
  XY = Y/n
  XX = crossprod(X)/n#+0.05*lambda[nlambda]*diag(d)
  beta = matrix(0,nrow=d,ncol=nlambda)
  ite.int = rep(0,nlambda)
  ite.int1 = rep(0,nlambda)
  if(intercept) {
    intcep=1
  }else{
    intcep=0
  }
  if(d<=3){
    num.scr1 = d
  }else{
    num.scr1 = ceiling(d/log(d))
  }
  order0 = order(abs(XY),decreasing = TRUE)
  idx.scr = order0; num.scr = length(idx.scr)
  idx.scr1 = order0[1:num.scr1]
  XX1 = XX[idx.scr,idx.scr]
  XXX = crossprod(XX1,XX1)
  gamma = max(colSums(abs(X)))
  begt = Sys.time()
  str=.C("slim_dantzig_ladm_scr2", as.double(XY), as.double(XX), as.double(XXX),
         as.double(beta), as.integer(n), as.integer(d), as.double(rho),
         as.integer(ite.int), as.integer(ite.int1), as.integer(num.scr1),
         as.integer(idx.scr), as.integer(idx.scr1), as.double(gamma), as.double(lambda),
         as.integer(nlambda), as.integer(max.ite), as.double(prec),
         as.integer(intcep),PACKAGE="freebird")
  runt = Sys.time() - begt
  beta.list = vector("list", nlambda)
  for(i in 1:nlambda){
    beta.i = unlist(str[4])[((i-1)*d+1):(i*d)]
    beta.list[[i]] = beta.i
  }
  ite.int = matrix(unlist(str[8]), byrow = FALSE, ncol = nlambda)
  ite.int1 = matrix(unlist(str[9]), byrow = FALSE, ncol = nlambda)
  ite = list()
  ite[[1]] = ite.int1
  ite[[2]] = ite.int
  return(list(beta=beta.list, ite=ite))
}






















library(Matrix)





## apply Huang 2016 with transformed model and delta method
# Sometimes the function Huang2016_noX using svd will give error message when p is large.
# This function use eigen rather than svd
Huang2016_noX1 = function(Y, Trt, Mediator, Trt.Mediator, var_per = 0.8, boot = 0){
  
  #' @param Y: a vector with lenght n or a n by 1 matrix, n-dimensional outcome vector (with column name such as "Y").
  #' @param Trt matrix or dataframe, n by 1 treatment matrix without intercept (with column name such as "T1", "T2",.., "Tq").
  #' @param Mediator matrix or dataframe, n by p mediator matrix (with column name, "G1",...).
  #' Note that we have p>n
  #' @param Trt.Mediator matrix or dataframe, n by p interaction matrix (with column name such as "T1.G1", "T1.G2",...).
  #' @param var_per: a number between 0 to 1, the proportion of principle component in \hat{Var(eps_G)} explained by U
  #' @param boot: the number of bootstrap for estimating the standard error. boot==0 then do NOT use bootstrap. boot>0 then add bootstrap SE estimator with bootstrap sample size boot
  #' @return NIE_inference is the matrix of Estimate, SE, test-stat, p value of the estimators
  #' @return PCA_sigmaG is a list of the matrix of U and the proportion of daig(Sigma_Ghat) explained by U
  #' @return Individual_inference: a list of the point est \hat{\theta} and Cov(\hat{\theta}) in eqn(15)
  
  p = ncol(Mediator)
  n = nrow(Mediator)
  q = ncol(Trt)
  
  fit.g<-lm(Mediator~Trt)
  Sigma_Ghat<-cov(fit.g$residual)
  
  # svds_hat<-svd(Sigma_Ghat) # note that Sigma is NND, so eigen-decomposition should also be the SVD.
  eigen_hat = eigen(Sigma_Ghat)
  # the singular values of Sigma is also eigen values of Sigma by checking the following code
  # max(abs( (eigen_hat$vectors)%*%diag(eigen_hat$values)%*%t(eigen_hat$vectors) - Sigma_Ghat ))
  # max(abs( svds_hat$u %*% diag(svds_hat$d) %*% t(svds_hat$u) - Sigma_Ghat  ))
  # max(abs( svds_hat$v %*% diag(svds_hat$d) %*% t(svds_hat$v) - Sigma_Ghat  ))
  # max(abs( eigen_hat$vectors[,1] - svds_hat$u[,1]  ))
  # max(abs( eigen_hat$vectors[,1] - svds_hat$v[,1]  ))
  # max(abs( eigen_hat$vectors[,2] - (-svds_hat$u[,2])  ))
  # max(abs( eigen_hat$vectors[,2] - (-svds_hat$v[,2])  ))
  # Note that eigen-decomposition is not unique for eigen-vectors corresponding to the same eigenvalue
  # In practice, if Sigma_Ghat is not of full rank, many zero eigenvalues occur
  # So svds_hat$u and svds_hat$v corresponding to those zero singular values (eigenvalues) will be different
  # But svds_hat$u and svds_hat$v corresponding to the same nonzero singular values (eigenvalues) will be the same, and same to those eigen(Sigma)$vectors up to a sign
  # max(abs( Sigma_Ghat - (svds_hat$u) %*% diag(svds_hat$d) %*% t(svds_hat$v) )) # which is almost zero
  
  n.factor<-sum(cumsum(eigen_hat$values)/sum(eigen_hat$values)<var_per)+1
  U<-eigen_hat$vectors[,1:n.factor]
  p.prime<-n.factor
  
  P<-Mediator%*%U # G^\star, independent columns
  
  fit.g2<-lm(P~Trt) # fit P[,j]~cbind(one.n, X, S) individually p times
  
  # (1:p.prime)*2 = 2, 4, ..., 2*p.prime
  s.delta<-diag( vcov(fit.g2) )[(1:p.prime)*2] # the diagonal entries of the upper p by p matrix of \hat{Cov}(\hat{\theta}) in eqn (15)
  
  if (p.prime == 1) {
    cov1.p<-s.delta # the upper p by p matrix of \hat{Cov}(\hat{\theta}) in eqn (15)
  } else {
    cov1.p<-diag(s.delta) # the upper p by p matrix of \hat{Cov}(\hat{\theta}) in eqn (15)
  }
  
  
  thetas.p<-NULL # estimates of theta^*_Gj, theta^*_Dj, j=1,..., p.prime by eqn(10) separately
  mse2<-NULL # sigma^2 for p individual regression of eqn(10)
  
  for (j in 1:p.prime){
    Pj<-P[,j]
    # XS = X*S; S.squared = S^2; SPj = S*Pj; lm(Y~X+S+XS+S.squared+Pj + SPj) # which fitting Y on cbind(one.n, S, XS, S.squared, Pj, SPj) is the same as the following procedure
    fit.y2pj<-lm(Y~Trt+Pj+Trt*Pj) # this is eqn (10) using independence. There should be a Trt^2, but Trt^2 = Trt, so we remove it. 
    thetas.p<-rbind(thetas.p, fit.y2pj$coef[c("Pj", "Trt:Pj")])
    mse2<-c(mse2, summary(fit.y2pj)$sigma^2)
  }
  
  Trt.P = matrix(0, n, p.prime) # the interaction matrix
  for (j in 1:p.prime) {
    Trt.P[,j] = Trt * P[,j] # this is a n by 1 matrix whose coordiates are S[i,1] * G[i,h], i=1,...,n
  }
  
  Z<-cbind(1, Trt, P, Trt.P)
  
  
  info<-( min(mse2)^(-1) ) * (t(Z)%*%Z) # this matrix is singular if 2*p.prime+2>n.
  
  cov2.p<-solve(info)[q+1+1:(2*p.prime), q+1+1:(2*p.prime)] # only pick the sub-matrix corresponding to P and Trt.P
  
  bcov.3p<-bdiag(cov1.p, cov2.p) # the Cov(\hat{\theta}) in eqn (15)
  
  bcov.2p = bcov.3p[1:(2*p.prime), 1:(2*p.prime)] # This is for \sum_{j=1}^p.prime beta_{1j}^\star \theta_{Gj}^\star
  
  theta.sum.p<-apply(thetas.p, 1, sum) # a vector of length p.prime: theta.G_j^star + theta.D_j^star, j=1,...,p.prime
  theta.G.p = thetas.p[, "Pj"] # the individual theta.G_j^star with length j=1:p.prime
  
  if (p.prime == 1) {
    beta.p <- fit.g2$coef["Trt"] # a vector of beta_{1j}^\star, with length j=1:p.prime
  } else {
    beta.p <- fit.g2$coef["Trt",] # a vector of beta_{1j}^\star, with length j=1:p.prime
  }
  
  pam.sum.p<-c(theta.sum.p, rep(beta.p, 2)) # derivative of g(x1,...,xp.prime,y1,...,yp.prim,z1,...,zp.prime) = \sum_{j=1}^p.prime xj(yj + zj) in delta method
  var.nie1<-(t(pam.sum.p)%*%bcov.3p%*%pam.sum.p)[1] # this is the covariance of \sum_{j=1}^p beta_{1j}^\star (theta.Gj^\star + theta.Dj^\star)
  se.nie1 = sqrt(var.nie1)
  pt.nie1 <- sum(beta.p * theta.sum.p)
  test.nie1 = pt.nie1/se.nie1
  pval.nie1 <- 2*pnorm( -abs( test.nie1 )  )
  
  
  pam.G.p<-c(theta.G.p, beta.p) # derivative of g(x1,...,xp.prime,y1,...,yp.prim,z1,...,zp.prime) = \sum_{j=1}^p.prime xj(yj + zj) in delta method
  var.nie0<-(t(pam.G.p)%*%bcov.2p%*%pam.G.p)[1] # this is the covariance of \sum_{j=1}^p beta_{1j}^\star (theta.Gj^\star + theta.Dj^\star)
  se.nie0 = sqrt(var.nie0)
  pt.nie0 <- sum(beta.p * theta.G.p)
  test.nie0 = pt.nie0/se.nie0
  pval.nie0 <- 2*pnorm( -abs( test.nie0 )  )
  
  NIE_inference = matrix(NA, 2, 4)
  rownames(NIE_inference) = c("NIE(T=0)", "NIE(T=1)")
  colnames(NIE_inference) = c("Est", "SE", "t-stat", "p-value")
  NIE_inference["NIE(T=0)",] = c(pt.nie0, se.nie0, test.nie0, pval.nie0)
  NIE_inference["NIE(T=1)",] = c(pt.nie1, se.nie1, test.nie1, pval.nie1)
  
  PCA_sigmaG = list()
  PCA_sigmaG$U = U # the principle components
  PCA_sigmaG$cum_diag = sum(eigen_hat$values[1:p.prime])/sum(eigen_hat$values) # the proportion of diag(Sigma_Ghat) explained by U
  PCA_sigmaG$p.prime = p.prime # the number of selected principle components
  
  Individual_inference = list()
  Individual_inference$beta.thetaG.thetaD = c(beta.p, thetas.p[, "Pj"], thetas.p[, "Trt:Pj"]) # the point estimates of transformed beta_{1j}^\star, j=1:p.prime, theta_{Gj}^\star, j=1:p.prime, theta_{Dj}^\star, j=1:p.prime
  Individual_inference$Cov = bcov.3p # the covariance for the transformed estimates: Individual_inference$beta.thetaG.thetaD as in eqn(15)
  
  
  if(boot == 0){
    return( return( list(NIE_inference = NIE_inference, 
                         PCA_sigmaG = PCA_sigmaG, 
                         Individual_inference = Individual_inference ) ) )
  } else {
    Y =  as.matrix(Y)
    est0.boot = rep(0, times = boot)
    est1.boot = rep(0, times = boot)
    for (b.loop in 1:boot) {
      index = sample(1:n, replace = TRUE)
      
      Y.boot = Y[index, ]; Trt.boot = Trt[index, ] 
      Mediator.boot = Mediator[index, ]; Trt.Mediator.boot = Trt.Mediator[index, ]
      
      fit.g.boot<-lm(Mediator.boot~Trt.boot)
      Sigma_Ghat.boot<-cov(fit.g.boot$residual)
      
      eigen_hat.boot<-eigen(Sigma_Ghat.boot) 
      n.factor.boot<-sum(cumsum(eigen_hat.boot$values)/sum(eigen_hat.boot$values)<var_per)+1
      U.boot<-eigen_hat.boot$vectors[,1:n.factor.boot]
      p.prime.boot<-n.factor.boot
      
      P.boot<-Mediator.boot%*%U.boot # G^\star, independent columns
      
      fit.g2.boot<-lm(P.boot~Trt.boot)
      
      
      thetas.p.boot<-NULL # estimates of theta^*_Gj, theta^*_Dj, j=1,..., p.prime by eqn(10) separately
      
      for (j in 1:p.prime.boot){
        Pj.boot<-P.boot[,j]
        # XS = X*S; S.squared = S^2; SPj = S*Pj; lm(Y~X+S+XS+S.squared+Pj + SPj) # which fitting Y on cbind(one.n, S, XS, S.squared, Pj, SPj) is the same as the following procedure
        fit.y2pj.boot<-lm(Y.boot~Trt.boot+Pj.boot+Trt.boot*Pj.boot) # this is eqn (10) using independence. There should be a Trt^2, but Trt^2 = Trt, so we remove it. 
        thetas.p.boot<-rbind(thetas.p.boot, fit.y2pj.boot$coef[c("Pj.boot", "Trt.boot:Pj.boot")])
      }
      
      theta.sum.p.boot<-apply(thetas.p.boot, 1, sum) # a vector of length p.prime: theta.G_j^star + theta.D_j^star, j=1,...,p.prime
      theta.G.p.boot = thetas.p.boot[, "Pj.boot"] # the individual theta.G_j^star with length j=1:p.prime
      
      if (p.prime.boot == 1) {
        beta.p.boot <- fit.g2.boot$coef["Trt.boot"] # a vector of beta_{1j}^\star, with length j=1:p.prime
      } else {
        beta.p.boot <- fit.g2.boot$coef["Trt.boot",] # a vector of beta_{1j}^\star, with length j=1:p.prime
      }
      
      est1.boot[b.loop] = sum(beta.p.boot * theta.sum.p.boot)
      est0.boot[b.loop] = sum(beta.p.boot * theta.G.p.boot)
      
      #print(b.loop)
    }
    
    return( list(NIE_inference = NIE_inference, 
                 PCA_sigmaG = PCA_sigmaG, 
                 Individual_inference = Individual_inference,
                 bootSE.NIE0 = sd(est0.boot),
                 bootSE.NIE1 = sd(est1.boot) ) ) 
  }
  
  
}
