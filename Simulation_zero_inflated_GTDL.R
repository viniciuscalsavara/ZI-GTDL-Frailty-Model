rm(list = ls())
# libraries ---------------------------------------------------------------

library(msm)

# final version -----------------------------------------------------------

# verossimilhança da inflação de zeros
# w e uma matriz da forma (1, x1, ..., x_h)
# tem que entrar com gama = (gama0, gama1, ..., gama_h)
veroZero <- function(w, gamma) {
    
    W0 <- as.matrix(w[w$tempo == 0, 2:ncol(w)])
    W <- as.matrix(w[w$tempo > 0, 2:ncol(w)])
    aux1 <- sum((W0%*%gamma)) - sum(log( 1 + exp(W0%*%gamma))) - sum(log(1 + exp(W%*%gamma)))
    return(-aux1)
}

# verossimilhança do GTDL com fragilidade Gama
# matriz x = (censura, tempo, x1)
# par e um vetor dado por (alpha, log(lambda), log(theta), beta)

veroGTDL <- function(x, par) {
    
    alpha <- par[1]
    lambda <- exp(par[2])
    theta <- exp(par[3])
    beta <- par[4]
    
    cens <- x[x$tempo > 0,1]
    tempo <- x[x$tempo > 0,2]
    X <- as.matrix(x[x$tempo > 0,3])
    
    aux2 <- log(lambda)*sum(cens) +
        sum( cens*(alpha*tempo + X%*%beta)) -
        sum( cens*log( 1 + exp(alpha*tempo + X%*%beta) ) ) - 
        sum( (cens + (1/theta))*log( 1 + ((theta*lambda/alpha)*log( (1 + exp(alpha*tempo + X%*%beta))/(1 + exp(X%*%beta)) )) ) )
    return(-aux2)
}


# função de sobrevivência do modelo com regressão apenas na inflação de zero
sobrevGTDLv1 <- function(t, x, par) {
    
    gamma <- par[1:2]
    alpha <- par[3]
    lambda <- par[4]
    theta <- par[5]
    beta <- par[6]
    
    tempo <- t
    X <- as.matrix(x)
    X0 <- matrix(data = c(rep(1,nrow(X)),X), nrow = nrow(X), ncol = 2)
    
    st <- (1/(1 + exp(X0%*%gamma) )) * (1 + ( (lambda*theta/alpha) * log( (1 + exp(alpha*tempo + X%*%beta) )/(1 + exp(X%*%beta))) ) )^(-1/theta)
    return(st)
}


# função que calcula a proporção de zeros inflacionados dados os níveis da covar
p0f <- function(w, gamma) {
    
    w <- as.matrix(w)
    W <- matrix(data = c(rep(1,nrow(w)),w), nrow = nrow(w), ncol = 2)
    
    p <- exp(W%*%gamma)/(1+exp(W%*%gamma))
    return(p)
}


# função que calcula a fração de cura dados os níveis da covar
p1f <- function(x, par) {
    gamma <- par[1:2]
    alpha <- par[3]
    lambda <- par[4]
    theta <- par[5]
    beta1 <- par[6]
    
    X <- as.matrix(x)
    W <- matrix(data = c(rep(1,nrow(X)),X), nrow = nrow(X), ncol = 2)
    
    p = (1/(1+exp(W%*%gamma))) * (1 + ((lambda*theta)/alpha)*log(1/(1+exp(X%*%beta1))) )^(-1/theta)
    return(p)
}


# função criada para usar no comando uniroot e encontrar a raiz
froot <- function(t,x,par,u) sobrevGTDLv1(t,x,par) - 1 + u


# parâmetros fixados para simulação
set.seed(2024)
gama0 <- -1.0; gama1 <- -0.6
alpha <- -0.1; lambda <- 0.5; theta <- 0.4; beta1 <- -0.9

n <- c(100, 200, 300, 500, 800, 1000, 2000) # tamanho da amostra
repeticao <- 1000 # quantidade de repeticoes

# dimensao 1 (linha) - repreticao 
# dimensao 2 (coluna) - parametro
# dimensao 3 - estimativa pontual ou erro-padrao
# dimensao 4 - varios tamanhos amostrais 

EMV_gtdl <- array(NA, dim = c(repeticao, 4, 2, length(n))) 
EMV_zero <- array(NA, dim = c(repeticao, 2, 2, length(n)))

# salvar o erro padrao para o calculo do IC de p0 
dp_p0 <- array(NA, dim=c(repeticao, 2, length(n)))

# salvar o erro padrao para o calculo do IC de p1 
dp_p1 <- array(NA, dim=c(repeticao, 2, length(n)))

# contagem dos possiveis erros 

erro_optim1 <- erro_optim2 <- erro_optim3 <- numeric()
erro_optim11 <- erro_optim12 <- erro_optim13 <- numeric()


# iteração pra simulação e estimação dos parâmetros
for(k in 1:length(n)){

    cat("n = ", n[k], "\n")
    j <- 1
    while(j <= repeticao) {
    t <- numeric()
    delta <- numeric()
    X <- rbinom(n[k], size = 1, prob = 0.5)
    p0 <- p0f(X, gamma = c(gama0, gama1))
    p1 <- p1f(X, par = c(gama0, gama1, alpha, lambda, theta, beta1))
    
  
      # processo de geração dos tempos
      for (i in 1:n[k]) {
        
        u1 <- runif(1,0,1)
        u2 <- runif(1, min = p0[i], max = (1-p1[i]))
        
        tf <- ifelse(test = u1 <= p0[i],
                     yes = 0,
                     no = ifelse(test = (1 - p1[i]) <= u1,
                                 yes = Inf,
                                 no = uniroot(f = froot, x=X[i], u=u2,
                                              par = c(gama0, gama1, alpha, lambda, theta, beta1),
                                              interval = c(0,1080))$root)
                      )
        
        tf <- as.numeric(tf)
        tc <- runif(1, 0, 65)
        t[i] <- min(tf,tc)
        
        delta[i] <- ifelse(t[i] == tf, 1, 0)
    }
    
    # processo de estimação dos parâmetros
    
    # definicao dos dados 
    dadoSimulado <- data.frame(cens = delta, tempo = t, x1 = X)
    dado0 <- data.frame(tempo = t, w0 = rep(1,n[k]), w1 = X)
    
    # variavel de teste
    TESTE1 <- FALSE
    TESTE2 <- FALSE
    
    # EMV do modelo GTDL com fragilidade gama
    emvg <- try(optim(par=c(-0.1,-0.5,-0.5,1), veroGTDL, x=dadoSimulado, hessian = TRUE, method = "BFGS"),TRUE)
    
    if (inherits(emvg, "try-error")==TRUE){erro_optim1 <- erro_optim1 + 1
    }else{
      v_solver <- try(solve(emvg$hessian), TRUE)
      if (inherits(v_solver, "try-error")==TRUE){erro_optim2 <- erro_optim2 + 1 
      }else{
        min_var <- min(diag(v_solver))
        if(min_var < 0){erro_optim3 <- erro_optim3 + 1
        }else{TESTE1 <- TRUE}
      }
    }

    
    emv0 <- try(optim(par=c(1,1), veroZero, w=dado0, hessian = TRUE, method = "BFGS"),TRUE)
    if (inherits(emv0, "try-error")==TRUE){erro_optim11 <- erro_optim11 + 1
    }else{
      v_solver0 <- try(solve(emv0$hessian), TRUE)
      if (inherits(v_solver0, "try-error")==TRUE){erro_optim12 <- erro_optim12 + 1 
      }else{
        min_var <- min(diag(v_solver0))
        if(min_var < 0){erro_optim13 <- erro_optim13 + 1
        }else{ TESTE2 <- TRUE}
      }
    }
    
    
    if(TESTE1 == TRUE & TESTE2 == TRUE){
        EMV_gtdl[j, , 1, k] <- c(emvg$par[1], exp(emvg$par[2]), exp(emvg$par[3]), emvg$par[4])
        # calculo do erro padrao
        sdg <- sqrt(diag(solve(emvg$hessian)))
        sdg[2] <- deltamethod(g = ~ exp(x1), mean = emvg$par[2], cov = solve(emvg$hessian)[2,2])
        sdg[3] <- deltamethod(g = ~ exp(x1), mean = emvg$par[3], cov = solve(emvg$hessian)[3,3])
        EMV_gtdl[j, , 2, k] <- sdg  
    
        EMV_zero[j, , 1, k] <- c(emv0$par[1], emv0$par[2])
        # calculo do erro padrao
        EMV_zero[j, , 2, k] <- sqrt(diag(v_solver0)) 
        
        
        dp_p0[j, 1, k] <- deltamethod(g = ~ exp(x1+x2) / (1 + exp(x1+x2)), mean = emv0$par, cov = solve(emv0$hessian))
        dp_p0[j, 2, k] <- deltamethod(g = ~ (exp(x1) / (1 + exp(x1))), mean = emv0$par, cov = solve(emv0$hessian))
        
        matriz_aux <- diag(6)
        matriz_aux[1:4, 1:4] <- emvg$hessian
        matriz_aux[5:6, 5:6] <- emv0$hessian
        EMV <- c(emvg$par, emv0$par)
        
        dp_p1[j, 1, k] <- deltamethod(g = ~ (1 - exp(x5+ x6) / (1 + exp(x5 + x6)))*(1 - ((exp(x2)*exp(x3))/x1) * log(1 + exp(x4)))^(-1/exp(x3)), mean = EMV, cov = solve(matriz_aux))
        dp_p1[j, 2, k] <- deltamethod(g = ~ (1 - exp(x5) / (1 + exp(x5)))*(1 - ((exp(x2)*exp(x3))/x1) * log(1 + exp(0)))^(-1/exp(x3)), mean = EMV, cov = solve(matriz_aux))
        
        
        j <- j+1
    }
  }
}  


########## analise dos resultados 

# verificacao dos erros 

# referente ao GTDL com fragilidade gama 
erro_optim1 
erro_optim2 
erro_optim3 
# referente ao zero inflacionado
erro_optim11
erro_optim12
erro_optim13


# salvando o resultado da simulação 
#save.image("teste_com_100_repeticoes.RData")


# analisando os resultados 
parametros <- c(alpha, lambda, theta, beta1, gama0, gama1)


# VIES 
# parâmetros do GTDL com fragilidade 
vies <- matrix(NA, nrow=length(n), ncol=10)
for(j in 1:4){
  for(i in 1:length(n)){
  vies[i,j] <- mean(EMV_gtdl[ , j , 1, i]) - parametros[j]
  }
}

# parâmetros da inflação de zeros 
for(j in 5:6){
  for(i in 1:length(n)){
    vies[i,j] <- mean(EMV_zero[ , (j-4) , 1, i]) - parametros[j]
  }
}

# da proporção de zeros 
p0_verd <- c(p0f(1, c(gama0, gama1)), p0f(0, c(gama0, gama1)))
for(i in 1:length(n)){
    soma0 <- 0
    soma1 <- 0
    for(k in 1:repeticao){
      soma1 <- soma1 + p0f(1, EMV_zero[k , 1: 2 , 1, i])
      soma0 <- soma0 + p0f(0, EMV_zero[k , 1: 2 , 1, i])
    }
    vies[i,7] <- soma1/repeticao - p0_verd[1]
    vies[i,8] <- soma0/repeticao - p0_verd[2]
}


# da fração de cura
p1_verd <- c(p1f(1, c(gama0, gama1, alpha, lambda, theta, beta1)), p1f(0, c(gama0, gama1, alpha, lambda, theta, beta1)))
for(i in 1:length(n)){
  soma0 <- 0
  soma1 <- 0
  for(k in 1:repeticao){
    soma1 <- soma1 + ifelse(EMV_gtdl[k, 1, 1, i]<0, p1f(1, c(EMV_zero[k, 1:2, 1, i], EMV_gtdl[k, 1:4, 1, i])), 0)
    soma0 <- soma0 + ifelse(EMV_gtdl[k, 1, 1, i]<0, p1f(0, c(EMV_zero[k, 1:2, 1, i], EMV_gtdl[k, 1:4, 1, i])), 0)
    }
  vies[i,9] <- soma1/repeticao - p1_verd[1]
  vies[i,10] <- soma0/repeticao - p1_verd[2]
}


library(xtable)
xtable(t(vies), digits = 4)



# EQM
eqm <- matrix(NA, nrow=length(n), ncol=10)
for(j in 1:4){
  for(i in 1:length(n)){
    eqm[i,j] <- mean((EMV_gtdl[ , j , 1, i] - parametros[j])^2)
  }
}

for(j in 5:6){
  for(i in 1:length(n)){
    eqm[i,j] <- mean((EMV_zero[ , (j-4) , 1, i] - parametros[j])^2)
  }
}


# da proporção de zeros 
p0_verd <- c(p0f(1, c(gama0, gama1)), p0f(0, c(gama0, gama1)))
for(i in 1:length(n)){
  soma0 <- 0
  soma1 <- 0
  for(k in 1:repeticao){
    soma1 <- soma1 + (p0f(1, EMV_zero[k , 1: 2 , 1, i]) - p0_verd[1])^2
    soma0 <- soma0 + (p0f(0, EMV_zero[k , 1: 2 , 1, i]) - p0_verd[2])^2
  }
  eqm[i,7] <- soma1/repeticao 
  eqm[i,8] <- soma0/repeticao 
}


# da fração de cura
p1_verd <- c(p1f(1, c(gama0, gama1, alpha, lambda, theta, beta1)), p1f(0, c(gama0, gama1, alpha, lambda, theta, beta1)))
for(i in 1:length(n)){
  soma0 <- 0
  soma1 <- 0
  for(k in 1:repeticao){
    soma1 <- soma1 + (ifelse(EMV_gtdl[k, 1, 1, i]<0, p1f(1, c(EMV_zero[k, 1:2, 1, i], EMV_gtdl[k, 1:4, 1, i])), 0) - p1_verd[1])^2
    soma0 <- soma0 + (ifelse(EMV_gtdl[k, 1, 1, i]<0, p1f(0, c(EMV_zero[k, 1:2, 1, i], EMV_gtdl[k, 1:4, 1, i])), 0) - p1_verd[2])^2
  }
  eqm[i,9] <- soma1/repeticao
  eqm[i,10] <- soma0/repeticao
}



xtable(t(eqm), digits = 4)


# Probabilidade de cobertura 
PC <- matrix(NA, nrow=length(n), ncol=10)
for(j in 1:4){
  for(i in 1:length(n)){
    PC[i,j] <- mean( (EMV_gtdl[ , j , 1, i] - 1.96* EMV_gtdl[ , j , 2, i] < parametros[j]) & 
                       (EMV_gtdl[ , j , 1, i] + 1.96* EMV_gtdl[ , j , 2, i] > parametros[j]) )
  }
}

for(j in 5:6){
  for(i in 1:length(n)){
    PC[i,j] <- mean( (EMV_zero[ , (j-4) , 1, i] - 1.96* EMV_zero[ , (j-4) , 2, i] < parametros[j]) & 
                       (EMV_zero[ , (j-4) , 1, i] + 1.96* EMV_zero[ , (j-4) , 2, i] > parametros[j]) )
  }
}


# da proporção de zeros 
p0_verd 

for(i in 1:length(n)){
  soma1 <- 0
  soma0 <- 0
  LI_1 <- LS_1 <- LI_0 <- LS_0 <- numeric()
    for(k in 1:repeticao){
    gama0 <- EMV_zero[k, 1, 1, i]
    gama1 <- EMV_zero[k, 2, 1, i]
    
      
    LI_1[k] <- exp(gama0 + gama1) / (1 + exp(gama0+gama1)) - 1.96* dp_p0[k,1,i]
    LS_1[k] <- exp(gama0 + gama1) / (1 + exp(gama0+gama1)) + 1.96* dp_p0[k,1,i]
    soma1 <- soma1 + ifelse(LI_1[k] < p0_verd[1] & LS_1[k] > p0_verd[1], 1, 0)  
    
    LI_0[k] <- exp(gama0) / (1 + exp(gama0)) - 1.96* dp_p0[k,2,i]
    LS_0[k] <- exp(gama0) / (1 + exp(gama0)) + 1.96* dp_p0[k,2,i]
    soma0 <- soma0 + ifelse(LI_0[k] < p0_verd[2] & LS_0[k] > p0_verd[2], 1, 0)
    #cat("LI_1=", LI_1, "LS_1=", LS_1, "\n")
    #cat("LI_0=", LI_0, "LS_0=", LS_0, "\n")
  }
  PC[i,7] <- soma1/repeticao 
  PC[i,8] <- soma0/repeticao 
}


head(cbind(LI_0, LS_0))

# da fração de cura
p1_verd 
for(i in 1:length(n)){
  soma1 <- 0
  soma0 <- 0
  for(k in 1:repeticao){
    alpha <- EMV_gtdl[k, 1, 1, i]
    lambda <- EMV_gtdl[k, 2, 1, i]
    theta <- EMV_gtdl[k, 3, 1, i]
    beta1 <- EMV_gtdl[k, 4, 1, i]
    gama0 <- EMV_zero[k, 1, 1, i]
    gama1 <- EMV_zero[k, 2, 1, i]
    if(alpha<0){
      LI_1 <- (1 - exp(gama0 + gama1)/(1+exp(gama0 + gama1)))*(1 - ((lambda*theta)/alpha)*log(1 + exp(beta1)))^(-1/theta) - 1.96* dp_p1[k,1,i]
      LS_1 <- (1 - exp(gama0 + gama1)/(1+exp(gama0 + gama1)))*(1 - ((lambda*theta)/alpha)*log(1 + exp(beta1)))^(-1/theta) + 1.96* dp_p1[k,1,i]
      soma1 <- soma1 + ifelse(LI_1 < p1_verd[1] & LS_1 > p1_verd[1], 1, 0)  
      #cat("k=", k, "soma1=", soma1, "\n")
      LI_1 <- (1 - exp(gama0)/(1+exp(gama0)))*(1 - ((lambda*theta)/alpha)*log(1 + exp(0)))^(-1/theta) - 1.96* dp_p1[k,2,i]
      LS_1 <- (1 - exp(gama0)/(1+exp(gama0)))*(1 - ((lambda*theta)/alpha)*log(1 + exp(0)))^(-1/theta) + 1.96* dp_p1[k,2,i]
      soma0 <- soma0 + ifelse(LI_1 < p1_verd[2] & LS_1 > p1_verd[2], 1, 0)  
      #cat("k=", k, "soma0=", soma0, "\n")
    }
  }
  PC[i,9] <- soma1/repeticao 
  PC[i,10] <- soma0/repeticao 
}


xtable(t(PC), digits = 4)
#round(PC, 4)











