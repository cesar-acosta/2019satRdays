# predictmlm.r        # implements Johnson & Wichern formula (7-49)

predictmlm <- function(object,newdata,level=0.95,interval = c("confidence", "prediction"))
{
 form <- as.formula(paste("~",as.character(formula(object))[3]))   # ~cyl + am + carb
 znew <- model.matrix(form, newdata)
 fit <- predict(object, newdata)
 
 Y <- model.frame(object)[,1]   # responses dataframe
 Z <- model.matrix(object)      # matrix with predictors (categorical as binaries)
 n <- nrow(Y)
 m <- ncol(Y)                   # n of responses
 p <- ncol(Z) - 1               # n predictors (counting binaries instead of categoricals)

 # alpha correction
 alpha <- 1 - level
 level <- 1 - m*alpha
 
# residual standard errors  (sigmahat_ii)
 sigmas <- sigma(object)^2 
 fit.var <- diag(znew %*% tcrossprod(solve(crossprod(Z)), znew))
 if(interval[1]=="prediction") fit.var <- fit.var + 1
 constant <- qf(level, df1=m, df2=n-p-m)*m*(n-p-1)/(n-p-m)
 vmat <- (n/(n-p-1)) * outer(fit.var, sigmas)

# boundaries
 lwr <- fit - sqrt(constant) * sqrt(vmat)
 upr <- fit + sqrt(constant) * sqrt(vmat)
 if(nrow(znew)==1)
 {
  ci <- rbind(fit, lwr, upr)
  rownames(ci) <- c("fit", "lwr", "upr")
 }
 else
 {
  ci <- array(0, dim=c(nrow(znew), m, 3))
  dimnames(ci) <- list(1:nrow(znew), colnames(Y), c("fit", "lwr", "upr") )
  ci[,,1] <- fit
  ci[,,2] <- lwr
  ci[,,3] <- upr
 }
ci
}
