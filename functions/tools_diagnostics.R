library(data.table)
library(ggplot2)
library(ggpubr)

mean_sd = function(x, sw){
  x_mean = sum(sw*x)/sum(sw)
  x_sd2 = (sum(sw))/(sum(sw)^2-sum(sw^2))*sum(sw*(x-x_mean)^2)
  return(list(x_mean, x_sd2))
}

p_sd = function(p, sw){
  p_hat = sum(sw*p)/sum(sw)
  p_hat_sd2 = p_hat*(1-p_hat)
  return(list(p_hat, p_hat_sd2))
}

StdD = function(x_mean, x_sd2){
  return(abs((x_mean[1]-x_mean[2])/sqrt((x_sd2[1]+x_sd2[2])/2)))
}

StdDiff.cont = function(data, sw=NULL, covs, group='group', specification = NA){
  if(is.null(sw)){
    sw = rep(1,dim(data)[1])
  }
  SD_data = NULL
  levels = as.vector(t(as.data.frame(data[,..group])))
  groups = unlist(unique(levels))
  for(var in covs){
    x = as.vector(t(as.data.frame(data[,..var])))
    X = NULL
    SD = NULL
    for(i in 1:length(groups)){
      index = which(levels==groups[i])
      sw_tmp = sw[index]
      x_tmp = x[index]
      stat = mean_sd(x_tmp, sw_tmp)
      X = c(X, stat[[1]])
      SD = c(SD, stat[[2]])
    }
    sd.matrix = cbind.data.frame('Method' = specification, 'Variable' = var)
    names = c('Method','Variable')
    for(j in 1:(length(groups)-1)){
      for(k in (j+1):length(groups)){
        sd.matrix = cbind.data.frame(sd.matrix,
                                     StdD(X[c(j,k)], SD[c(j,k)]))
        names = c(names, paste0('StdD',paste0(j,k)))
      }
      colnames(sd.matrix) = names
    }
    SD_data = rbind.data.frame(SD_data, sd.matrix)
  }
  return(SD_data)
}

StdDiff.bin = function(data, sw=NULL, covs, group='group', specification = NA){
  if(is.null(sw)){
    sw = rep(1,dim(data)[1])
  }
  SD_data = NULL
  levels = as.vector(t(as.data.frame(data[,..group])))
  groups = unlist(unique(levels))
  for(var in covs){
    x = as.vector(t(as.data.frame(data[,..var])))
    P = NULL
    P_SD = NULL
    for(i in 1:length(groups)){
      index = which(levels==groups[i])
      sw_tmp = sw[index]
      x_tmp = x[index]
      stat = p_sd(x_tmp, sw_tmp)
      P = c(P, stat[[1]])
      P_SD = c(P_SD, stat[[2]])
    }
    sd.matrix = cbind.data.frame('Method' = specification, 'Variable' = var)
    names = c('Method','Variable')
    for(j in 1:(length(groups)-1)){
      for(k in (j+1):length(groups)){
        sd.matrix = cbind.data.frame(sd.matrix,
                                     StdD(P[c(j,k)], P_SD[c(j,k)]))
        names = c(names, paste0('StdD',paste0(j,k)))
      }
      colnames(sd.matrix) = names
    }
    SD_data = rbind.data.frame(SD_data, sd.matrix)
  }
  return(SD_data)
}


summary_table = function(SW_list,exposure){
  summary_tab = NULL
  for(k in 1:length(SW_list)){
    summary_tab = rbind.data.frame(summary_tab,
                                   cbind.data.frame('Exposure'= c(exposure),
                                                    'IPTW' = c(k),
                                                    'Mean' = c(mean(SW_list[[k]])),
                                                    'Sd' = c(sd(SW_list[[k]])),
                                                    'Min' = c(range(SW_list[[k]])[1]),
                                                    'Max' = c(range(SW_list[[k]])[2])))
    
  }
  return(summary_tab)
}


iptw.balance.plot <- function(love_data, max.diff=0.15, metric.x.name=NULL, title=NULL,
                              colori = c('hotpink3','#E69F00','darkcyan', 'dodgerblue3','gray25'),
                              dim.text = 0.8, dim.axis = 0.9, dim.title = 0.95, dim.legend = 0.9,
                              leg.pos = 'right'){
  
  if(is.null(metric.x.name)){
    metric.x.name = 'Metric'
  }
  var = love_data[Method=='Unadjusted']$Variable
  values = love_data[Method=='Unadjusted']$Metric
  love_data$Variable <- factor(love_data$Variable, levels=var[order(values)])
  ymax = max(values) + 0.05
  
  p=ggplot(love_data, aes(x=Metric, y=Variable, color=Method, group=Method, shape=Method, fill=Method)) +
    geom_point(size=1) + 
    scale_shape_manual(values=c(21:25)) + 
    scale_color_manual(values=colori) +
    scale_fill_manual(values=colori) +
    stat_summary(geom="line", size=0.5) +
    scale_x_continuous(limits = c(0,ymax))+
    geom_vline(xintercept=max.diff, linetype="dashed", color='gray40', alpha=0.5) +
    xlab(metric.x.name) +
    ylab('Confounder') +
    ggtitle(title) +
    theme_bw() + 
    theme(legend.position=leg.pos, legend.title = element_text(size=rel(dim.legend)),
          axis.text=element_text(size=rel(dim.text)), axis.title=element_text(size=rel(dim.axis)),
          legend.text = element_text(size=rel(dim.legend)),
          plot.title = element_text(face="bold", size=rel(dim.title), hjust = 0.5) )

  return(p)
}

