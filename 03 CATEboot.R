#################################################
# Conditional Average Treatment Effects (CATEs) #
#     95% Bootstrap Confidence Intervals        #
#################################################
rm( list = ls() )
graphics.off()
library(data.table)
library(ipw)
library(survminer)
library(survival)

# Set working directory
setwd("~/github/causalRDI")
# Load utils
source('functions/tools_bootstrap.R')
source('functions/tools_reshapes.R')


# Load data & fitted models
load('data/toy_dataset.Rdata')
load('results/fitted_models.Rdata')

## Estimated CATEs
times = seq(1,120)
fit = fit.CoxMSM(df)
resultCATE = CATE.RMST(fit, times)
resultCATE

## Bootstrap CIs

# Step 1-2) Sub-groups
df[A=='0: RDI>=85%' & V=='0', group := 'g=(0,P)']
df[A=='1: 70%<=RDI<85%' & V=='0', group := 'g=(1,P)']
df[A=='2: RDI<70%' & V=='0', group := 'g=(2,P)']
df[A=='0: RDI>=85%' & V=='1', group := 'g=(0,G)']
df[A=='1: 70%<=RDI<85%' & V=='1', group := 'g=(1,G)']
df[A=='2: RDI<70%' & V=='1', group := 'g=(2,G)']

# Step 3) Compute unequal sampling probability based on IPTW
df = cbind(df, SW)
df[group=='g=(0,P)', tot0P := sum(SW)]
df[group=='g=(1,P)', tot1P := sum(SW)]
df[group=='g=(2,P)', tot2P := sum(SW)]
df[group=='g=(0,G)', tot0G := sum(SW)]
df[group=='g=(1,G)', tot1G := sum(SW)]
df[group=='g=(2,G)', tot2G := sum(SW)]

df[group=='g=(0,P)', prob := SW/tot0P]
df[group=='g=(1,P)', prob := SW/tot1P]
df[group=='g=(2,P)', prob := SW/tot2P]
df[group=='g=(0,G)', prob := SW/tot0G]
df[group=='g=(1,G)', prob := SW/tot1G]
df[group=='g=(2,G)', prob := SW/tot2G]


# Step 4) Sampling
B = 100
times = seq(1,120)
boot_CATE1P = matrix(NA, nrow = B, ncol=length(times))
boot_CATE1G = matrix(NA, nrow = B, ncol=length(times))
boot_CATE2P = matrix(NA, nrow = B, ncol=length(times))
boot_CATE2G = matrix(NA, nrow = B, ncol=length(times))
set.seed(3009)
for(b in 1:B){
  df_b = df.boot(df, a.var='A', V.var='V')
  fit_b = fit.CoxMSM(df_b)
  res_b = CATE.RMST(fit_b, times)
  boot_CATE1P[b, ] = res_b$CATE1P
  boot_CATE1G[b, ] = res_b$CATE1G
  boot_CATE2P[b, ] = res_b$CATE2P
  boot_CATE2G[b, ] = res_b$CATE2G
}

# Step 5) define the bounds of the 95% bootstrap CI
Nt = length(times)
CIboot = rbind.data.frame(cbind.data.frame('tau' = times,
                                           'V' = rep(0,Nt),
                                           'Adiff' = rep('1-0',Nt),
                                           'lower' = apply(boot_CATE1P, 2, quantile, probs=0.025, na.rm=T),
                                           'upper' = apply(boot_CATE1P, 2, quantile, probs=0.975, na.rm=T)),
                          cbind.data.frame('tau' = times,
                                           'V' = rep(1,Nt),
                                           'Adiff' = rep('1-0',Nt),
                                           'lower' = apply(boot_CATE1G, 2, quantile, probs=0.025, na.rm=T),
                                           'upper' = apply(boot_CATE1G, 2, quantile, probs=0.975, na.rm=T)),
                          cbind.data.frame('tau' = times,
                                           'V' = rep(0,Nt),
                                           'Adiff' = rep('2-0',Nt),
                                           'lower' = apply(boot_CATE2P, 2, quantile, probs=0.025, na.rm=T),
                                           'upper' = apply(boot_CATE2P, 2, quantile, probs=0.975, na.rm=T)),
                          cbind.data.frame('tau' = times,
                                           'V' = rep(1,Nt),
                                           'Adiff' = rep('2-0',Nt),
                                           'lower' = apply(boot_CATE2G, 2, quantile, probs=0.025, na.rm=T),
                                           'upper' = apply(boot_CATE2G, 2, quantile, probs=0.975, na.rm=T))
)
CIboot = as.data.table(CIboot)
summary(CIboot)

save('CIboot',file='results/CIboot_cate_prob.Rdata')


## Figure 6 ##
##----------##
plotCATE = reshapeCATE(resultCATE, times)
plotCATE

# Poor Responders (PRs)
ggplot(plotCATE[V==0], aes(x=tau, y=RD, group=Adiff, colour=Adiff)) + 
  geom_hline(yintercept=0, color='gray30',  linetype="dashed") +
  geom_ribbon(aes(ymin=CIboot[V==0]$lower, ymax=CIboot[V==0]$upper, group=Adiff, fill=Adiff), 
              linetype=0, alpha=0.2) + geom_line(size=1) +
  scale_colour_manual(values=c('dodgerblue3','darkorange'), name='RDI-strategy') + 
  scale_fill_manual(values=c('dodgerblue3','darkorange'), name='RDI-strategy') +
  labs(y="Estimated months gained (if >0) or lost (if <0)",
       x="Time since end of therapy [months]") +
  ggtitle('Poor Responders: Estimated CATEs over time') + theme_bw() 

# Good Responders (GRs)
ggplot(plotCATE[V==1], aes(x=tau, y=RD, group=Adiff, colour=Adiff)) + 
  geom_hline(yintercept=0, color='gray30',  linetype="dashed") +
  geom_ribbon(aes(ymin=CIboot[V==1]$lower, ymax=CIboot[V==1]$upper, group=Adiff, fill=Adiff), 
              linetype=0, alpha=0.2) + geom_line(size=1) +
  scale_colour_manual(values=c('dodgerblue3','darkorange'), name='RDI-strategy') + 
  scale_fill_manual(values=c('dodgerblue3','darkorange'), name='RDI-strategy') +
  labs(y="Estimated months gained (if >0) or lost (if <0)",
       x="Time since end of therapy [months]") +
  ggtitle('Good Responders: Estimated CATEs over time') + theme_bw() 

