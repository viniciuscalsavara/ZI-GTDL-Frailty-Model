
##########
#Gompertz#
##########


###########################
#Fun??o de verossimilhan?a#
###########################

#########################
#Covari?vel em p0 e no b#
#########################

################
#Uma covari?vel#
################

log.like_cov2=function(a,b0,b1,b2,b3){
  
  p0=exp(b0+x*b1)/(1+exp(b0+x*b1))
  b=exp(b2+x*b3)
  
  (-1)*(sum((1-delta0)*log(p0)) +sum(delta0*log(1-p0))+sum(log(b)*delta0*delta) + a*sum(delta0*delta*s)-(1/a)*sum(b*delta0*(exp(a*s)-1))  )
  
}

##################
#Duas covari?veis#
##################

log.like_cov3=function(a,b0,b1,b2,b3,c1,c2){
  
  p0=exp(b0+cov1*b1+cov2*c1)/(1+exp(b0+cov1*b1+cov2*c1))
  b=exp(b2+cov1*b3+cov2*c2)
  
  (-1)*(sum((1-delta0)*log(p0)) +sum(delta0*log(1-p0))+sum(log(b)*delta0*delta) + a*sum(delta0*delta*s)-(1/a)*sum(b*delta0*(exp(a*s)-1))  )
  
}


###########
#Intera??o#
###########


log.like_cov4=function(a,b0,b1,b2,b3,c1,c2,d1,d2){
  
  p0=exp(b0+cov1*b1+cov2*c1+cov3*d1)/(1+exp(b0+cov1*b1+cov2*c1+cov3*d1))
  b=exp(b2+cov1*b3+cov2*c2+cov3*d2)
  
  (-1)*(sum((1-delta0)*log(p0)) +sum(delta0*log(1-p0))+sum(log(b)*delta0*delta) + a*sum(delta0*delta*s)-(1/a)*sum(b*delta0*(exp(a*s)-1))  )
  
}



################



bx=function(b2,b3,cov){
  
  b0=exp(b2+b3*cov[1])
  b1=exp(b2+b3*cov[2])
  aux=cbind(b0,b1)
  return(aux)
}



bx_2=function(b2,b3,c2,cov1,cov2){
  
  b=exp(b2+cov1*b3+cov2*c2)
  return(b)
}




p1x=function(a,p0x,b0x,cov){
  
  
  aux0=(1-p0x[1])*exp(b0x[1]/a)
  aux1=(1-p0x[2])*exp(b0x[2]/a)
  aux=cbind(aux0,aux1)
  return(aux)
}



p1x_2=function(a,p0x,b0x){
  
  aux0=(1-p0x)*exp(b0x/a)
  return(aux0)
}



sg_cov2=function(a,p0x,b0x,t){
  
  sp0=(1-p0x[1])*exp(-(b0x[1]/a)*(exp(a*t)-1))
  sp1=(1-p0x[2])*exp(-(b0x[2]/a)*(exp(a*t)-1))
  aux=cbind(sp0,sp1)
  return(aux)
}



sg_cov2_2=function(a,p0x_0,p0x_1,p0x_2,b0x_0,b0x_1,b0x_2,t){
  
  sp0=(1-p0x_0)*exp(-(b0x_0/a)*(exp(a*t)-1))
  sp1=(1-p0x_1)*exp(-(b0x_1/a)*(exp(a*t)-1))
  sp2=(1-p0x_2)*exp(-(b0x_2/a)*(exp(a*t)-1))
  aux=cbind(sp0,sp1,sp2)
  return(aux)
}


###############
#K covari?veis#
###Gompertz ###
###############

log.like_gomp=function(par){
  
  index=dim(cova)[2]
  index1=length(par)
  a=par[1]
  beta=par[2:(index+1)]
  gama=par[(index+2):index1]
  
  X_p0=tcrossprod(as.matrix(cova),t(beta))
  X_b=tcrossprod(as.matrix(cova),t(gama))
  
  p0=exp(X_p0)/(1+exp(X_p0))
  b=exp(X_b)
  
  (-1)*(sum((1-delta0)*log(p0)) +sum(delta0*log(1-p0))+sum(log(b)*delta0*delta) + a*sum(delta0*delta*s)-(1/a)*sum(b*delta0*(exp(a*s)-1))  )
  
}



##Covariate in a (a(x)=a*alpha)

log.like_gomp_complete <- function(par){
  
  index  <- ncol(cova)
  
  # parâmetros:
  # alpha: efeitos em a(x)   (tamanho = index)
  # beta : efeitos em p0(x)  (tamanho = index)
  # gama : efeitos em b(x)   (tamanho = index)
  alpha <- par[1:index]
  beta  <- par[(index+1):(2*index)]
  gama  <- par[(2*index+1):(3*index)]
  
  # preditores lineares
  X_a  <- as.vector(cova %*% alpha)
  X_p0 <- as.vector(cova %*% beta)
  X_b  <- as.vector(cova %*% gama)
  
  # parâmetros por indivíduo
  a  <- X_a                      # a(x) > 0
  p0 <- exp(X_p0) / (1 + exp(X_p0))   # logistic
  b  <- exp(X_b)                      # b(x) > 0
  
  # contribuições
  # parte de p0 (zero-adjustment / inflação de zero)
  ll_p0 <- sum((1 - delta0) * log(p0)) + sum(delta0 * log(1 - p0))
  
  # parte Gompertz (para os suscetíveis delta0=1)
  # log f(t) = log b + a*t - (b/a)*(exp(a*t) - 1)  (quando delta=1)
  # log S(t) = - (b/a)*(exp(a*t) - 1)             (quando delta=0)
  # unificando com delta:
  ll_gomp <- sum(delta0 * delta * (log(b) + a * s)) -
    sum(delta0 * (b / a) * (exp(a * s) - 1))
  
  # negative log-likelihood
  - (ll_p0 + ll_gomp)
}


###########
##Summary##
###########



summary.optim <- function ( optim.fit ){
  
  Estimate = optim.fit$par
  cov.mat = solve(optim.fit$hessian)
  Std.Error = sqrt( diag(cov.mat)	)
  z.value = Estimate/Std.Error
  p.value = 2*( 1-pnorm(abs(z.value)) )
  
  IC_lower=Estimate-qnorm(0.975)*Std.Error
  IC_upper=Estimate+qnorm(0.975)*Std.Error
  
  ds.out = round(cbind(Estimate, Std.Error,IC_lower,IC_upper, z.value, p.value),4)
  
  colnames(ds.out) <- c('Estimate', 'Std Error', 'IC lower','IC upper','z value', 'p value')
  
  index=length(Estimate)
  index2=(index-1)/2+1
  betas=c("Intercept(p0)",rep("p",(index2-2)))
  gamas=c("Intercept(b)",rep("b",(index2-2)))
  name=c("a",betas,gamas)
  
  rownames(ds.out)<-name
  
  ds.out
  
}






##

summary.optim_complete <- function(optim.fit, covar_names = NULL){
  
  Estimate <- optim.fit$par
  H <- optim.fit$hessian
  
  if (is.null(H) || anyNA(H)) stop("Hessian is missing/NA in optim.fit.")
  if (!is.matrix(H)) H <- as.matrix(H)
  
  # matriz de covariância (com fallback caso Hessiana seja quase singular)
  cov.mat <- tryCatch(
    solve(H),
    error = function(e) {
      message("Hessian is singular/ill-conditioned; using generalized inverse (MASS::ginv).")
      if (!requireNamespace("MASS", quietly = TRUE)) stop("Install MASS to use ginv().")
      MASS::ginv(H)
    }
  )
  
  Std.Error <- sqrt(pmax(diag(cov.mat), 0))  # evita sqrt de negativos numéricos
  z.value   <- Estimate / Std.Error
  p.value   <- 2 * (1 - pnorm(abs(z.value)))
  
  IC_lower <- Estimate - qnorm(0.975) * Std.Error
  IC_upper <- Estimate + qnorm(0.975) * Std.Error
  
  ds.out <- cbind(Estimate, Std.Error, IC_lower, IC_upper, z.value, p.value)
  colnames(ds.out) <- c("Estimate", "Std Error", "IC lower", "IC upper", "z value", "p value")
  
  # ---- nomes dos parâmetros: alpha (a), beta (p0), gamma (b)
  idx <- length(Estimate)
  
  if (idx %% 3 != 0) {
    stop("Expected length(par) multiple of 3: par = c(alpha, beta, gamma).")
  }
  
  p <- idx / 3  # número de covariáveis (inclui intercepto)
  
  if (is.null(covar_names)) {
    covar_names <- paste0("x", seq_len(p))
    covar_names[1] <- "Intercept"
  } else {
    if (length(covar_names) != p) stop("covar_names must have length = ncol(cova).")
  }
  
  rn <- c(paste0("alpha_a(", covar_names, ")"),
          paste0("beta_p0(", covar_names, ")"),
          paste0("gamma_b(", covar_names, ")"))
  
  rownames(ds.out) <- rn
  
  round(ds.out, 4)
}

####Fun??o de sobreviv?ncia dadas as caracter?sticas###


sob_gompertz=function(cova,par,time){
  
  index=length(cova)
  index1=length(par)
  a=par[1]
  beta=par[2:(index+1)]
  gama=par[(index+2):index1]
  
  X_p0=tcrossprod((cova),t(beta))
  X_b=tcrossprod((cova),t(gama))
  
  p0=exp(X_p0)/(1+exp(X_p0))
  b=exp(X_b)
  
  sob=(1-p0)*exp(-(b/a)*(exp(a*time)-1))
  
  
  return(sob)
  
}




###########################
#M?todo delta: Erro padr?o#
###########################

ep_g=function(coef,vcov){
  
  ep=numeric()
  
  ep=sqrt(diag(vcov(v1)))
  
  
  #Erro padr?o para p0(x)=exp(b0+b1*x)/(1+exp(b0+b1*x))#
  
  #Os x s?o os respectivos par?metros da saida v1#
  
  #x1 ? coef(v1)[1] =a
  #x2 ? coef(v1)[2] =beta0
  
  #p0(0)=exp(x1)/(1+exp(x1))
  ep[6]=deltamethod (~ exp(x2) / (1+exp(x2)), coef(v1), vcov(v1)) 
  
  #p0(1)=exp(x1+x2)/(1+exp(x1+x2))
  ep[7]=deltamethod (~ exp(x2+x3) / (1+exp(x2+x3)), coef(v1), vcov(v1)) 
  
  
  
  
  #Cura#
  
  #p1(x)=(1/(1+exp(beta0+beta1*x)))*exp((beta2+beta3*x)/a)
  
  #p1(0)
  ep[8]=deltamethod (~ (1/(1+exp(x2))*exp(exp(x4)/x1)   ), coef(v1), vcov(v1)) 
  
  #p1(1)
  ep[9]=deltamethod (~ (1/(1+exp(x2+x3))*exp(exp(x4+x5)/x1)   ), coef(v1), vcov(v1)) 
  
  
  return(round(ep,3))
  
  
}




#######################################
###erro padr?o para duas covari?veis###
#######################################

ep_g_2=function(coef,vcov){
  
  ep=numeric()
  
  ep=sqrt(diag(vcov(v1)))
  
  
  #Erro padr?o para p0(x)=exp(b0+b1*x)/(1+exp(b0+b1*x))#
  
  #Os x s?o os respectivos par?metros da saida v1#
  
  #x1 ? coef(v1)[1] =a
  #x2 ? coef(v1)[2] =beta0
  
  #p0(0,0)=exp(x1)/(1+exp(x1))
  ep[8]=deltamethod (~ exp(x2) / (1+exp(x2)), coef(v1), vcov(v1)) 
  
  #p0(0,1)=exp(x1+x2)/(1+exp(x1+x2))
  ep[9]=deltamethod (~ exp(x2+x3) / (1+exp(x2+x3)), coef(v1), vcov(v1)) 
  
  
  #p0(1,0)=exp(x1+x6)/(1+exp(x1+x6))
  ep[10]=deltamethod (~ exp(x2+x6) / (1+exp(x2+x6)), coef(v1), vcov(v1)) 
  
  
  
  #Cura#
  
  #p1(x)=(1/(1+exp(beta0+beta1*x+c1*x2)))*exp((beta2+beta3*x+c2*x2)/a)
  
  #p1(0)
  ep[11]=deltamethod (~ (1/(1+exp(x2))*exp(exp(x4)/x1)   ), coef(v1), vcov(v1)) 
  
  #p1(0,1)
  ep[12]=deltamethod (~ (1/(1+exp(x2+x3))*exp(exp(x4+x5)/x1)   ), coef(v1), vcov(v1)) 
  
  #p1(0,1)
  ep[13]=deltamethod (~ (1/(1+exp(x2+x6))*exp(exp(x4+x7)/x1)   ), coef(v1), vcov(v1)) 
  
  
  return(round(ep,3))
  
  
}












#########################################################
###Erro padr?o para p0 e p1 considerando K covari?veis###
#########################################################

ep_gompertz=function(v1,vcov,glicemia1,DMG1,DMG2,tabagismo,idade1,diabetes,macrossomia,hac,IMC1,IMC2){
  
  ep=numeric()
  
  
  #Erro padr?o para p0(x)=exp(beta*X)/(1+exp(beta*X))#
  
  #Os x s?o os respectivos par?metros da saida v1#
  
  #x1 ? coef(v1)[1]=a
  #x2 ? coef(v1)[2]=beta0
  
  ep[1]=deltamethod (~ exp(x2+x3*glicemia1+x4*DMG1+x5*DMG2+x6*tabagismo+x7*idade1+x8*diabetes+x9*macrossomia+x10*hac+x11*IMC1+x12*IMC2) / (1+exp(x2+x3*glicemia1+x4*DMG1+x5*DMG2+x6*tabagismo+x7*idade1+x8*diabetes+x9*macrossomia+x10*hac+x11*IMC1+x12*IMC2)), v1, vcov) 
  
  
  
  #Erro padr?o para a Cura#
  
  #p1(x)=(1/(1+exp(beta*X)))*exp((gama*X)/a)																								
  
  ep[2]=deltamethod (~ (1/(1+exp(x2+x3*glicemia1+x4*DMG1+x5*DMG2+x6*tabagismo+x7*idade1+x8*diabetes+x9*macrossomia+x10*hac+x11*IMC1+x12*IMC2))*exp(exp(x13+x14*glicemia1+x15*DMG1+x16*DMG2+x17*tabagismo+x18*idade1+x19*diabetes+x20*macrossomia+x21*hac+x22*IMC1+x23*IMC2)/x1)   ), v1, vcov) 
  
  erro_padrao=round(cbind(ep[1],ep[2]),3)
  colnames(erro_padrao) <- c('p0', 'p1')
  rownames(erro_padrao) <- c('Standard error')
  
  
  erro_padrao
  
  
}


#########################
#Fun??o de sobreviv?ncia#
#########################


s_g=function(a,b,p0,t){
  
  sp=(1-p0)*exp(-(b/a)*(exp(a*t)-1))
  
  return(sp)
}




s_g_cov_c=function(a,b,p0x,t){
  
  sp=(1-p0x)*exp(-(b/a)*(exp(a*t)-1))
  
  return(sp)
}


s_g_cov=function(a,b,p0x,t){
  
  sp0=(1-p0x[1])*exp(-(b/a)*(exp(a*t)-1))
  sp1=(1-p0x[2])*exp(-(b/a)*(exp(a*t)-1))
  aux=cbind(sp0,sp1)
  return(aux)
}




ic_g=function(est,ep){
  
  ic=matrix(nrow=length(est),ncol=2)
  colnames(ic) <- c("Inferior","Superior")
  
  ic[,1]=round(est-qnorm(0.975)*ep,3)
  ic[,2]=round(est+qnorm(0.975)*ep,3)
  aux=data.frame(ic)
  return(ic)
  
}



######################################
#Estimativa da probabilidade de zeros#
######################################

px=function(b0,b1,cov){
  
  aux=exp(b0+cov*b1)/(1+exp(b0+cov*b1))
  return(aux)
}


px_2=function(b0,b1,c1,cov1,cov2){
  
  aux=exp(b0+cov1*b1+cov2*c1)/(1+exp(b0+cov1*b1+cov2*c1))
  return(aux)
}



prop_zeros=function(est,cova1,v1){
  
  
  index=length(cova1)
  index1=length(est)
  a=est[1]
  beta=est[2:(index+1)]
  gama=est[(index+2):index1]
  
  X_p0=tcrossprod((cova1),t(beta))
  X_b=tcrossprod((cova1),t(gama))
  
  p0=exp(X_p0)/(1+exp(X_p0))
  
  
  return(p0)
}




prop_cura=function(est,cova1,v1){
  
  
  index=length(cova1)
  index1=length(est)
  a=est[1]
  beta=est[2:(index+1)]
  gama=est[(index+2):index1]
  
  X_p0=tcrossprod((cova1),t(beta))
  X_b=tcrossprod((cova1),t(gama))
  
  p0=exp(X_p0)/(1+exp(X_p0))
  b=exp(X_b)
  
  cura=(1-p0)*exp(b/a)
  
  return(cura)
}



