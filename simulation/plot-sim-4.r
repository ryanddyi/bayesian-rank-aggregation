
library(ggplot2)
library(RColorBrewer)
library(reshape)


setwd(paste0('/Users/ydd/Documents/BARC/simulation_new/','sim-4'))

param=expand.grid(seedID=c(1:100),std=c(1,2,5,10,20))
#hyperparam=expand.grid(lambda=c(5,10),s2=c(0.1,0.5))

acc = array(NA, dim = c(500, 5))
for (i in 1:dim(acc)[1]){
  bat = (i-1) %/% 100 + 1
  folder.char = paste('batch', bat, sep = '')
  acc[i, ] = readRDS(paste(folder.char, '/acc', i, '.rds', sep = ''))
}

acc = acc[,-c(1,2)]

df_acc = c()
for (std in unique(param$std)){
  df_acc = rbind(df_acc, c(std, colMeans(acc[which(param$std == std),])))
}

names = c('StDev', 'K=3', 'K=4', 'K=5')
colnames(df_acc) = names
df_acc = as.data.frame(df_acc)
df_acc$StDev = factor(df_acc$StDev)

attributes(pdmelt$StDev)
pdmelt <- melt(df_acc, id.vars = c('StDev'))

p1 <- ggplot(pdmelt, aes(x = StDev, y = value, group = variable, color = variable)) +
  geom_line(size = 1) +
  theme_bw() +
  ggtitle('Sensitivity to the number of categories') +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  labs(x = 'Noise Level', y = 'Classification Accuracy', color = 'BARCM') +
  scale_colour_brewer(palette="Paired") 

require(gridExtra)

setwd('/Users/ydd/Documents/BARC/simulation_new/')
pdf('sim-acc.pdf', height = 4, width = 6)
print(p1)
dev.off()
