library(data.table)
library(reshape2)

reshapeCATE <- function(cate_list, tau){
  Nt = length(tau)
  plot.data = rbind.data.frame(cbind.data.frame('tau' = tau,
                                                'V' = rep(0,Nt),
                                                'Adiff' = rep('1-0',Nt),
                                                'RD' = as.numeric(cate_list$CATE1P)),
                               cbind.data.frame('tau' = tau,
                                                'V' = rep(1,Nt),
                                                'Adiff' = rep('1-0',Nt),
                                                'RD' = as.numeric(cate_list$CATE1G)),
                               cbind.data.frame('tau' = tau,
                                                'V' = rep(0,Nt),
                                                'Adiff' = rep('2-0',Nt),
                                                'RD' = as.numeric(cate_list$CATE2P)),
                               cbind.data.frame('tau' = tau,
                                                'V' = rep(1,Nt),
                                                'Adiff' = rep('2-0',Nt),
                                                'RD' = as.numeric(cate_list$CATE2G))
  )
  
  plot.data = as.data.table(plot.data)
  plot.data$V = as.factor(plot.data$V)
  
  return(plot.data)
}
