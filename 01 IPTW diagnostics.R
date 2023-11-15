################################################
#              Model diagnostics for           #
#  Inverse Probability of Treatment Weighting  #
################################################
rm( list = ls() )
library(data.table)
library(ggplot2)
library(ggpubr)
library(ipw)
library(splines)
library(survival)

# Set working directory
setwd("~/github/causalRDI")
# Load utils
source('functions/tools_diagnostics.R')

# Load ddata
load('data/toy_dataset.Rdata')

Npt = length(unique(df$ID))


#----------------
# UNADJUSTED
#----------------
data_cont = df[,.(A, motox_rule_pre, motox_rule_post, motox_gen_pre, motox_gen_post)]
unadj_cont = StdDiff.cont(data_cont, group='A', 
                          covs=c('motox_rule_pre','motox_rule_post','motox_gen_pre','motox_gen_post'),
                          specification = 'Unadjusted')

data_bin = df[,.(A,trial,age,gender,V)]
data_bin$BO06 = ifelse(data_bin$trial=='BO06',1,0)
data_bin$child = ifelse(data_bin$age=='child',1,0)
data_bin$adolescent = ifelse(data_bin$age=='adolescent',1,0)
data_bin$adult = ifelse(data_bin$age=='adult',1,0)
data_bin$male = ifelse(data_bin$gender=="M",1,0)
data_bin$good = ifelse(data_bin$V=="1",1,0)
data_bin = data_bin[,-c(2:5)]

unadj_bin = StdDiff.bin(data_bin, group='A', 
                        covs=c('BO06','child','adolescent','adult','male','good'),
                        specification = 'Unadjusted')

#----------------
# SPECIFICATION 1
#----------------
tmp1 <- ipwpoint(exposure = A, family='multinomial',
                 numerator = ~ V,
                 denominator = ~ trial + age + gender + V +
                   motox_rule_pre + motox_gen_pre +  motox_rule_post + motox_gen_post,
                 data = df, trace=F)

SW_1 = tmp1$ipw.weights

spe1_cont = StdDiff.cont(data_cont, group='A', 
                         sw = tmp1$ipw.weights,
                         covs=c('motox_rule_pre','motox_rule_post','motox_gen_pre','motox_gen_post'),
                         specification = 'IPTW 1')

spe1_bin = StdDiff.bin(data_bin, group='A',
                       sw = tmp1$ipw.weights,
                       covs=c('BO06','child','adolescent','adult','male','good'),
                       specification = 'IPTW 1')


#----------------
# SPECIFICATION 2
#----------------
tmp2 <- ipwpoint(exposure = A, family='multinomial',
                 numerator = ~ V,
                 denominator = ~ trial + age + gender + V +
                   motox_rule_pre*motox_gen_pre +  motox_rule_post*motox_gen_post,
                 data = df, trace=F)

SW_2 = tmp2$ipw.weights

spe2_cont = StdDiff.cont(data_cont, group='A', 
                         sw = tmp2$ipw.weights,
                         covs=c('motox_rule_pre','motox_rule_post','motox_gen_pre','motox_gen_post'),
                         specification = 'IPTW 2')

spe2_bin = StdDiff.bin(data_bin, group='A',
                       sw = tmp2$ipw.weights,
                       covs=c('BO06','child','adolescent','adult','male','good'),
                       specification = 'IPTW 2')


#----------------
# SPECIFICATION 3
#----------------
tmp3 <- ipwpoint(exposure = A, family='multinomial',
                 numerator = ~ V,
                 denominator = ~ trial + age + gender + V +
                   motox_rule_pre + motox_gen_pre +  motox_rule_post + motox_gen_post +
                   trial:motox_rule_pre + trial:motox_gen_pre +  trial:motox_rule_post + trial:motox_gen_post,
                 data = df, trace=F)

SW_3 = tmp3$ipw.weights


spe3_cont = StdDiff.cont(data_cont, group='A', 
                         sw = tmp3$ipw.weights,
                         covs=c('motox_rule_pre','motox_rule_post','motox_gen_pre','motox_gen_post'),
                         specification = 'IPTW 3')

spe3_bin = StdDiff.bin(data_bin, group='A',
                       sw = tmp3$ipw.weights,
                       covs=c('BO06','child','adolescent','adult','male','good'),
                       specification = 'IPTW 3')


#----------------
# SPECIFICATION 4
#----------------
tmp4 <- ipwpoint(exposure = A, family='multinomial',
                 numerator = ~ V,
                 denominator = ~ trial + age + gender + V +
                   bs(motox_rule_pre,6) + bs(motox_gen_pre,6) + 
                   bs(motox_rule_post,6) + bs(motox_gen_post,6),
                 data = df, trace=F)
SW_4 = tmp4$ipw.weights


spe4_cont = StdDiff.cont(data_cont, group='A', 
                         sw = tmp4$ipw.weights,
                         covs=c('motox_rule_pre','motox_rule_post','motox_gen_pre','motox_gen_post'),
                         specification = 'IPTW 4')

spe4_bin = StdDiff.bin(data_bin, group='A',
                       sw = tmp4$ipw.weights,
                       covs=c('BO06','child','adolescent','adult','male','good'),
                       specification = 'IPTW 4')

#----------------
# SPECIFICATION 5
#----------------
tmp5 <- ipwpoint(exposure = A, family='multinomial',
                 numerator = ~ V,
                 denominator = ~ trial + age + gender + V +
                   motox_rule_pre + motox_gen_pre +  motox_rule_post + motox_gen_post +
                   V:motox_rule_pre + V:motox_gen_pre +  V:motox_rule_post + V:motox_gen_post,
                 data = df, trace=F)

SW_5 = tmp5$ipw.weights


spe5_cont = StdDiff.cont(data_cont, group='A', 
                         sw = tmp5$ipw.weights,
                         covs=c('motox_rule_pre','motox_rule_post','motox_gen_pre','motox_gen_post'),
                         specification = 'IPTW 5')

spe5_bin = StdDiff.bin(data_bin, group='A',
                       sw = tmp5$ipw.weights,
                       covs=c('BO06','child','adolescent','adult','male','good'),
                       specification = 'IPTW 5')


## Figure 4 - right panel ##
##------------------------##
## Balance plot
love_rdi = data.table(rbind.data.frame(unadj_cont,unadj_bin,
                                       spe1_cont,spe1_bin,
                                       spe2_cont,spe2_bin,
                                       spe3_cont,spe3_bin,
                                       spe4_cont,spe4_bin,
                                       spe5_cont,spe5_bin))
love_rdi[, ASD_mean := rowMeans(.SD), .SDcols = 3:5]

ggplot(love_rdi, aes(x=ASD_mean, y=Variable, color=Method, group=Method, 
                     shape=Method, fill=Method)) +
  geom_point() + stat_summary(fun.x=sum, geom="line", size=0.5) +
  geom_vline(xintercept=0.15, linetype="dashed", color='gray10', alpha=0.5) +
  scale_color_manual(values=c('#CC0066','#FF9933','gold','darkcyan', 'dodgerblue3','gray25')) +
  labs(x='Mean Abs. Std Diff.', y='Confounder') +
  ggtitle('Confounders balance plot') + theme_bw()


## Table 3 ##
##---------##
summary_table(list(SW_1,SW_2,SW_3,SW_4,SW_5),'RDI level')

