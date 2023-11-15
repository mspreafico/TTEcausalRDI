library(data.table)
library(ipw)
library(survival)

fit.CoxMSM <- function(df){
  iptw_RDIV = ipwpoint(exposure = A, family='multinomial',
                       numerator = ~ V,
                       denominator = ~ trial + age + gender + V + 
                         motox_rule_pre +  motox_gen_pre +  
                         motox_rule_post + motox_gen_post,
                       data = df, trace=F)
  model <- coxph(Surv(efs_time, efs_status) ~ A*V,
                 data=df, weights=iptw_RDIV$ipw.weights, robust=TRUE)
  df_surv <- with(df,
                  data.frame(A = rep(c('0: RDI>=85%','1: 70%<=RDI<85%','2: RDI<70%'),2),
                             V = as.factor(c(0,0,0,1,1,1)) )
  )
  fit <- survfit(model, newdata = df_surv)
  return(fit)
}

CATE.RMST <- function(fit, tau){
  
  time = c(0,fit$time)
  CATE1P = NULL
  CATE1G = NULL
  CATE2P = NULL
  CATE2G = NULL
  for(j in tau){
    nt = length(time[time<=j])
    if(nt>1){
      deltat = c(time[2:nt],j) - time[1:nt]
      # Survival curve up to tau
      s0P = c(1,fit$surv[1:(nt-1),1])
      s0G = c(1,fit$surv[1:(nt-1),4])
      s1P = c(1,fit$surv[1:(nt-1),2])
      s1G = c(1,fit$surv[1:(nt-1),5])
      s2P = c(1,fit$surv[1:(nt-1),3])
      s2G = c(1,fit$surv[1:(nt-1),6])
      
      # Risk: RMST
      R0P = sum(deltat*s0P)
      R0G = sum(deltat*s0G)
      R1P = sum(deltat*s1P)
      R1G = sum(deltat*s1G)
      R2P = sum(deltat*s2P)
      R2G = sum(deltat*s2G)
    }else{
      R0P = R0G = R1P = R1G = R2P = R2G = 1
    }
    
    # CATE
    CATE1P = c(CATE1P, R1P - R0P)
    CATE1G = c(CATE1G, R1G - R0G)
    CATE2P = c(CATE2P, R2P - R0P)
    CATE2G = c(CATE2G, R2G - R0G)
  }
  
  return(list('CATE1P' = CATE1P,
              'CATE1G' = CATE1G,
              'CATE2P' = CATE2P,
              'CATE2G' = CATE2G))
}

df.boot <- function(df, a.var, V.var){
  row = 1:dim(df)[1]
  df = cbind(df,row)
  strategies = unique(df[, get(a.var)])
  effect.mod = unique(df[, get(V.var)])
  sampled.rows = NULL
  for(a in strategies){
    for(v in effect.mod){
      sub.df = df[get(a.var)==a & get(V.var)==v]
      sampled = sample(sub.df$row, length(sub.df$row), 
                       replace = T, prob = sub.df$prob)
      sampled.rows = c(sampled.rows, sampled)
    }
  }
  df.b = df[sampled.rows,]
  return(df.b)
}
