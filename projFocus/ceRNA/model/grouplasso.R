regfit      = function(X, y, group){
  #Desp:    function to perform group lasso regression  
  #input:   datamatrix with first column as response, other as variables, must be named, group of the variables, must be named
  #output:  coefficients beta, residual, r_square
  #usage:   function internal called by regfitPermu
  require(grpreg)
  if (setequal(colnames(X), names(group))){
    #search for alpha:
    rss   = NULL
    a     = seq(0.01,1,0.1)
    ##tune params
    for (alpha in a){
      grpfit            = grpreg(X,y,group=group,
                                 penalty="grLasso",family="gaussian",
                                 alpha = alpha,eps=0.005,max.iter=10000)
      grpfit_selected   = select(grpfit,criterion="AIC")
      rss               = c(rss,sum(y - grpfit_selected$beta[1] -
                                      X %*% as.matrix(grpfit_selected$beta[-1])))
    }
    alpha = a[which.min(rss)[1]]
    #fit model
    grpfit              = grpreg(X, y, group=group,penalty="grLasso",
                                 alpha = alpha,eps=0.005,max.iter=5000)
    grpfit_selected     = select(grpfit,criterion="AIC")
    
    beta                = grpfit_selected$beta
    residual            = y - beta[1] - X %*% as.matrix(beta[-1])
    residual_sum_squre  = sum((residual)^2)
    reg_sum_square      = sum((X %*% as.matrix(beta[-1]) + beta[1] - mean(y))^2)
    total_sum_square    = reg_sum_square + residual_sum_squre
    r_square            = reg_sum_square / total_sum_square  
    
    return(list(beta    = beta,
                residual          = residual, 
                r_square          = r_square,
                RSS               = residual_sum_squre,
                fit               = grpfit_selected))
  }
  else { print("ERROR! group not match data!")
  }
}