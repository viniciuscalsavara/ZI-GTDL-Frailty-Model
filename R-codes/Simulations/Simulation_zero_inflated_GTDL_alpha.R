#############################################################################
#Simulation study: ZI-GTDL frailty model with covariate into alpha parameter#
#############################################################################


library(foreach)
library(doParallel)
library(doRNG)

script_arg <- grep('^--file=', commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg)) dirname(normalizePath(sub('^--file=', '', script_arg[1]))) else getwd()
source(file.path(script_dir, 'functions_GTDL_alpha.R'))

params_fixed <- list(gamma0=-1.029, gamma1=-0.302,
                     alpha0=-0.134, alpha1=-0.009,
                     lambda=0.493, theta=0.44, beta1=-0.892)
n <- c(100,200,300,500,800,1000,2000)
repeticao <- 1000L
censor_max <- 65
max_attempts_per_n <- 20000L



filter_se <- TRUE
se_threshold <- 20
filter_parameters <- c('alpha0','alpha1','lambda','theta','beta1','gamma0','gamma1')

truth <- c(alpha0=params_fixed$alpha0, alpha1=params_fixed$alpha1,
           lambda=params_fixed$lambda, theta=params_fixed$theta,
           beta1=params_fixed$beta1, gamma0=params_fixed$gamma0,
           gamma1=params_fixed$gamma1,
           p0_x1=p0_total(1,params_fixed), p0_x0=p0_total(0,params_fixed),
           p1_x1=p1_total(1,params_fixed), p1_x0=p1_total(0,params_fixed))
stopifnot(all(filter_parameters %in% names(truth)),
          is.finite(se_threshold), se_threshold > 0)
cat('Filtro de SE: ', if (filter_se) paste0('SE < ',se_threshold) else 'desligado', '\n')

invalid_reason <- function(z, par_names, filter_se,
                           filter_parameters, se_threshold) {
  if (inherits(z,'error')) return(paste0('worker_error: ',conditionMessage(z)))
  if (!is.list(z)) return('invalid_worker_output')
  if (!isTRUE(z$ok)) {
    if (!is.null(z$reason) && length(z$reason)==1L && !is.na(z$reason) &&
        nzchar(z$reason)) return(as.character(z$reason))
    return('replication_not_ok')
  }
  if (is.null(z$est) || is.null(z$se) ||
      !all(par_names %in% names(z$est)) || !all(par_names %in% names(z$se)))
    return('missing_estimates_or_se')
  est <- z$est[par_names]
  se <- z$se[par_names]
  if (any(!is.finite(est)) || any(!is.finite(se)) || any(se < 0))
    return('nonfinite_estimate_or_se')
  if (filter_se) {
    exceeded <- names(se[filter_parameters])[se[filter_parameters] >= se_threshold]
    if (length(exceeded)) {
      return(paste0('se_ge_', se_threshold, ':', paste(exceeded, collapse=';')))
    }
  }
  NA_character_
}

run_simulation <- function(n,repeticao,params_fixed,censor_max,truth,
                           max_attempts_per_n,filter_se,
                           filter_parameters,se_threshold) {
  available <- parallel::detectCores()
  if (is.na(available)) available <- 2L
  n_cores <- max(1L,available-2L)
  cl <- parallel::makeCluster(n_cores)
  on.exit(parallel::stopCluster(cl),add=TRUE)
  doParallel::registerDoParallel(cl)
  doRNG::registerDoRNG(2024)
  results <- vector('list',length(n))
  accounting <- vector('list',length(n))
  reasons <- vector('list',length(n))
  for (k in seq_along(n)) {
    nk <- n[k]
    accepted <- list()
    discarded_reasons <- character(0)
    attempted <- 0L
    batches <- 0L
    cat('\nn =',nk,'\n')
    while(length(accepted)<repeticao && attempted<max_attempts_per_n) {
      remaining <- repeticao-length(accepted)
      batch_size <- min(remaining,max_attempts_per_n-attempted)
      batches <- batches+1L
      batch <- foreach(j=seq_len(batch_size),.errorhandling='pass',
        .export=c('softplus','alpha_params','gtdl_cumhaz','p0_total','p1_total',
                  'surv_total','gtdl_nll','zero_nll','num_grad','delta_se','one_rep_alpha')) %dorng% {
          one_rep_alpha(nk=nk,params_fixed=params_fixed,censor_max=censor_max)
      }
      attempted <- attempted+length(batch)
      reason <- vapply(batch,invalid_reason,character(1),
                       par_names=names(truth),filter_se=filter_se,
                       filter_parameters=filter_parameters,
                       se_threshold=se_threshold)
      good <- is.na(reason)
      if (any(good)) accepted <- c(accepted,batch[good])
      if (any(!good)) discarded_reasons <- c(discarded_reasons,reason[!good])
      cat(sprintf('Lote %d: geradas %d; validas %d/%d; descartadas %d\n',
                  batches,batch_size,length(accepted),repeticao,length(discarded_reasons)))
    }
    accounting[[k]] <- data.frame(n=nk,target_valid=repeticao,
      attempted=attempted,valid=length(accepted),discarded=length(discarded_reasons),
      extra_samples=attempted-repeticao,
      failure_rate=length(discarded_reasons)/attempted,
      failure_percent=100*length(discarded_reasons)/attempted,batches=batches)
    tab <- sort(table(discarded_reasons),decreasing=TRUE)
    if(length(tab)) reasons[[k]] <- data.frame(n=nk,reason=names(tab),count=as.integer(tab))
    if(length(accepted)!=repeticao) {
      write.csv(accounting[[k]],file.path(script_dir,sprintf('GTDL_alpha_INCOMPLETE_SElt20_n_%d.csv',nk)),row.names=FALSE)
      stop(sprintf('n=%d: apenas %d validas apos %d tentativas.',nk,length(accepted),attempted))
    }
    results[[k]] <- accepted
  }
  failure_table <- if(any(!vapply(reasons,is.null,logical(1))))
    do.call(rbind,Filter(Negate(is.null),reasons)) else
    data.frame(n=integer(),reason=character(),count=integer())
  list(res=results,accounting_table=do.call(rbind,accounting),failure_table=failure_table)
}

mc_summary <- function(res,n,truth) {
  rows <- list()
  for(k in seq_along(n)) for(nm in names(truth)) {
    accepted <- res[[k]]
    zhat <- vapply(accepted,function(z) unname(z$est[nm]),numeric(1))
    se <- vapply(accepted,function(z) unname(z$se[nm]),numeric(1))
    rows[[length(rows)+1L]] <- data.frame(n=n[k],parameter=nm,
      true_value=unname(truth[nm]),valid=length(zhat),
      bias=mean(zhat-truth[nm]),rmse=sqrt(mean((zhat-truth[nm])^2)),
      empirical_sd=sd(zhat),mean_se=mean(se),
      coverage95=mean(abs(zhat-truth[nm])<=qnorm(.975)*se))
  }
  do.call(rbind,rows)
}

extract_mc_estimates <- function(res,n,truth,script_dir) {
  par_names <- names(truth)
  all_results <- vector('list',length(n))
  for(k in seq_along(n)) {
    rows <- lapply(seq_along(res[[k]]),function(i) {
      z <- res[[k]][[i]]
      est <- setNames(as.list(unname(z$est[par_names])),paste0('est_',par_names))
      se <- setNames(as.list(unname(z$se[par_names])),paste0('se_',par_names))
      data.frame(n=n[k],replication=i,est,se,check.names=FALSE)
    })
    all_results[[k]] <- do.call(rbind,rows)
    write.csv(all_results[[k]],file.path(script_dir,
      sprintf('GTDL_alpha_MC_estimates_n%d_1000valid_SElt20.csv',n[k])),row.names=FALSE)
  }
  combined <- do.call(rbind,all_results)
  rownames(combined) <- NULL
  write.csv(combined,file.path(script_dir,'GTDL_alpha_MC_estimates_ALL_1000valid.csv'),row.names=FALSE)
  combined
}

sim <- run_simulation(n,repeticao,params_fixed,censor_max,truth,
  max_attempts_per_n,filter_se,filter_parameters,se_threshold)
res <- sim$res
accounting_table <- sim$accounting_table
failure_table <- sim$failure_table
summary_table <- mc_summary(res,n,truth)
mc_estimates <- extract_mc_estimates(res,n,truth,script_dir)

print(accounting_table,row.names=FALSE)
print(failure_table,row.names=FALSE)
print(summary_table,row.names=FALSE,digits=5)
cat('\nTotal de amostras extras:',sum(accounting_table$extra_samples),'\n')

write.csv(accounting_table,file.path(script_dir,'GTDL_alpha_MC_replication_counts.csv'),row.names=FALSE)
write.csv(failure_table,file.path(script_dir,'GTDL_alpha_MC_failure_reasons.csv'),row.names=FALSE)
write.csv(summary_table,file.path(script_dir,'GTDL_alpha_MC_summary_1000valid.csv'),row.names=FALSE)
saveRDS(list(params_fixed=params_fixed,n=n,target_valid=repeticao,
  filter_settings=list(filter_se=filter_se,se_threshold=se_threshold,
    filter_parameters=filter_parameters),
  simulations=res,accounting_table=accounting_table,
  failure_table=failure_table,summary_table=summary_table),
  file.path(script_dir,'GTDL_alpha_MC_results_1000valid.rds'))

