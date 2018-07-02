
library(ggplot2)
library(RColorBrewer)
library(reshape)

setwd('/Users/ydd/Documents/BARC/simulation_new/sim-2/')

param=expand.grid(seedID=c(1:500), m=10, n=50, std=c(5,10,20,40))
#hyperparam=expand.grid(lambda=c(5,10),s2=c(0.1,0.5))

kd = array(NA, dim = c(2000, 6))
for (i in 1:dim(kd)[1]){
  bat = (i-1) %/% 500 + 1
  folder.char = paste('batch', bat, sep = '')
  kd[i, ] = as.numeric(read.table(paste(folder.char, '/KD', i, '.txt', sep = '')))
}
kd = 1-(kd+1)/2


kdsum = c()
for (std in unique(param$std)){
    kdsum = rbind(kdsum, c(std, colMeans(kd[which(param$std == std),])))
}

names = c('StDev', 'BC', 'MC1', 'MC2', 'MC3', 'BARC', 'PL')
colnames(kdsum) = names
kdsum = as.data.frame(kdsum)
kdsum$StDev = factor(kdsum$StDev)


#display.brewer.all()
#display.brewer.pal()
attributes(pdmelt$StDev)
pdmelt <- melt(kdsum, id.vars = c('StDev'))
pdf('sim-2-kd.pdf', height = 4, width = 4)
p <- ggplot(pdmelt, aes(x = StDev, y = value, group = variable, color = variable)) +
  geom_line(size = 0.8) +
  theme_bw() +
  theme(legend.position="none") + 
  ggtitle(parse(text = 'paste("Scenario 2")')) +
  scale_y_continuous(limits = c(0, 0.5), expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  labs(x = 'Noise Level', y = NULL, color = 'Method') +
  scale_colour_brewer(palette="Paired") 
print(p)
dev.off()

#library(RColorBrewer)
#display.brewer.all()
