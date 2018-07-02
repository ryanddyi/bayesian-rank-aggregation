
library(ggplot2)
library(RColorBrewer)
library(reshape)


setwd(paste0('/Users/ydd/Documents/BARC/simulation_new/','sim-0'))

param=expand.grid(seedID=c(1:100),std=c(1,2,5,10,20))
#hyperparam=expand.grid(lambda=c(5,10),s2=c(0.1,0.5))

M = 69
id_all = array(NA, dim = c(500, M))
weight_all = id_all
KDist = array(NA, dim = c(500, 3))
acc = rep(0, 500)

for (i in 1:500){
  bat = (i-1) %/% 100 + 1
  folder.char = paste('batch', bat, sep = '')
  rds = readRDS(paste(folder.char, '/output', i, '.rds', sep = ''))
  id_all[i,] = rds$idCat0
  weight_all[i,] = rds$weight_est
  KDist[i,] = rds$KDist
  acc[i] = rds$accuracy
}
KDist = 1-(KDist+1)/2

combine_std = function(acc, param){
  df_acc = c()
  for (std in unique(param$std)){
    df_acc = rbind(df_acc, c(std, colMeans(acc[which(param$std == std),])))
  }
  df_acc
}

KDist = combine_std(KDist, param)

names = c('StDev', 'Category 1', 'Category 2', 'Category 3')
colnames(KDist) = names
KDist = as.data.frame(KDist)
KDist$StDev = factor(KDist$StDev)

attributes(pdmelt$StDev)
pdmelt <- melt(KDist, id.vars = c('StDev'))

p1 <- ggplot(pdmelt, aes(x = StDev, y = value, group = variable, color = variable)) +
  geom_line(size = 1) +
  theme_bw() +
  ggtitle('Rank aggregation by category') +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  labs(x = 'Noise Level', y = parse(text = 'paste("Kendall\'s ", tau, " Distance")'), color = 'BARCM') +
  scale_colour_brewer(palette="Paired") 

require(gridExtra)

setwd('/Users/ydd/Documents/BARC/simulation_new/')
pdf('sim-mix.pdf', height = 4, width = 6)
print(p1)
dev.off()


df_weight = cbind(as.vector(id_all[201:300,]),as.vector(weight_all[201:300,]))
df_weight = as.data.frame(df_weight)

p<-ggplot(df_weight, aes(x=V1, y=V2,group = V1))+
  theme_bw() +
  ggtitle("BARCW -- Rankers' weights by category") +
  geom_boxplot()+
  labs(x = 'Category', y = 'Weight')
print(p)

setwd('/Users/ydd/Documents/BARC/simulation_new/')
pdf('sim-weights-cat.pdf', height = 4, width = 6)
print(p)
dev.off()
