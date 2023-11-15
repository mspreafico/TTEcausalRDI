###########################################################
# Causal inference through marginal structural Cox models #
#                  with effect modifications              #
###########################################################
rm( list = ls() )
graphics.off()
library(data.table)
library(ipw)
library(survival)
library(survminer)

# Set working directory
setwd("~/github/causalRDI")

# Load data
load('data/toy_dataset.Rdata')


# SELECTED IPTW MODEL #
#---------------------#
# IPTW 1
iptw_RDI = ipwpoint(exposure = A, family='multinomial',
                    numerator = ~ V,
                    denominator = ~ trial + age + gender + V + motox_rule_pre + 
                      motox_gen_pre +  motox_rule_post + motox_gen_post,
                    data = df, trace=F)
SW = iptw_RDI$ipw.weights
mean(SW)
sd(SW)
range(SW)


# COX MSM #
#---------#
CoxMSM <- coxph(Surv(efs_time, efs_status) ~ A*V,
                data=df, weights=SW, robust=TRUE)
# Schoenfeld residuals
ggcoxdiagnostics(CoxMSM, type='schoenfeld')


# UNWEIGTHED COX MODEL #
#----------------------#
unw_Cox <- coxph(Surv(efs_time, efs_status) ~ A*V, data=df)

save(SW, unw_Cox, CoxMSM, file='results/fitted_models.Rdata')


## Table 4 ##
##---------##
# Estimated parameters with 95% CIs
# Cox MSM (left Betas)
print(summary(CoxMSM))
# Unweighted Cox model (right Betas_unw)
print(summary(unw_Cox))


## Figure 5 ##
##----------##
# Poor Responders (PRs)
PR_surv <- with(df,
                data.frame(A = c('0: RDI>=85%','1: 70%<=RDI<85%','2: RDI<70%'),
                           V = c('0','0','0') ))
fitPR <- survfit(CoxMSM, newdata = PR_surv)
ggsurvplot(fitPR,  PR_surv, censor=F,
           palette = c('gray50','dodgerblue3','darkorange'),
           legend.labs = c('Standard (a=0)','Reduced (a=1)','Highly-reduced (a=2)'),
           ylab = 'Event-Free Survival probability',
           xlab = 'Time since end of therapy [months]',
           title = 'Poor Responders (PRs)')
# Good Responders (GRs)
GR_surv <- with(df,
                data.frame(A = c('0: RDI>=85%','1: 70%<=RDI<85%','2: RDI<70%'),
                           V = c('1','1','1') ))
fitGR <- survfit(CoxMSM, newdata = GR_surv)
ggsurvplot(fitGR,  GR_surv, censor=F,
           palette = c('gray50','dodgerblue3','darkorange'),
           legend.labs = c('Standard (a=0)','Reduced (a=1)','Highly-reduced (a=2)'),
           ylab = 'Event-Free Survival probability',
           xlab = 'Time since end of therapy [months]',
           title = 'Good Responders (GRs)')
