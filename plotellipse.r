plotellipse <- function(mod,newdata,level = 0.95,interval = c("confidence","prediction"))
{
  p <- predict(mod, newdata)

  # center of ellipse
  cent <- c(p[1,1],p[1,2])

  # shape of ellipse
  Z <- model.matrix(mod)
  Y <- mod$model[[1]]
  n <- nrow(Y)         
  m <- ncol(Y)         
  r <- ncol(Z) - 1    
  S <- crossprod(resid(mod))/(n-r-1)

  # ellipse radius
  tt <- terms(mod)
  Terms <- delete.response(tt)
  mf <- model.frame(Terms, newdata, na.action = na.pass,xlev = mod$xlevels)
  z0 <- model.matrix(Terms, mf, contrasts.arg = mod$contrasts)
  
  zz <- z0 %*% solve(t(Z)%*%Z) %*% t(z0)
  if(interval[1]=="prediction") zz <- zz + 1  
  rad <- sqrt((m*(n-r-1)/(n-r-m))*qf(level,m,n-r-m)*zz)
  epoints <- ellipse(center = c(cent), shape = S, radius = c(rad), draw = F)

  # plot  
    resps <- colnames(mod$coefficients)  
    plot(epoints, type = "l",xlab=resps[1],ylab=resps[2],main="",col="red")
    points(x = cent[1],y = cent[2],cex=0.7)
    grid()

}

