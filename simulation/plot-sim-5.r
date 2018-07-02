
library(ggplot2)
library(RColorBrewer)
library(reshape)


setwd(paste0('/Users/ydd/Documents/BARC/simulation_new/','sim-5'))


M = 69
acc = array(dim = c(100,4))
nCat = array(dim = c(100,4))

for (i in 1:100){
  bat = (i-1) %/% 100 + 1
  folder.char = paste('batch', bat, sep = '')
  accrds = readRDS(paste(folder.char, '/acc', i, '.rds', sep = ''))
  acc[i,] = accrds
  ncatrds = readRDS(paste(folder.char, '/nCat', i, '.rds', sep = ''))
  nCat[i,] = ncatrds
}

colMeans(acc)
colMeans(nCat)

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
