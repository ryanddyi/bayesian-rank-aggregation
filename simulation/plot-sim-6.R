
library(ggplot2)
library(RColorBrewer)
library(reshape)

setwd('/Users/ydd/Documents/BARC/simulation_new/sim-6/')

param=expand.grid(seedID=c(1:100),m=10,g=c(1,2,4,8,10,16),std=5)

kd = array(NA, dim = c(600, 3))
for (i in 1:dim(kd)[1]){
  bat = (i-1) %/% 100 + 1
  folder.char = paste('batch', bat, sep = '')
  kd[i, ] = as.numeric(read.table(paste(folder.char, '/KD', i, '.txt', sep = '')))
}
kd = 1-(kd+1)/2

kdsum = c()
for (std in unique(param$g)){
  kdsum = rbind(kdsum, c(std, colMeans(kd[which(param$g == std),])))
}

names = c('NumberOfEntities', 'BARC(0.5)', 'BARC(1)', 'BAR')
colnames(kdsum) = names

kdsum = as.data.frame(kdsum)
kdsum$NumberOfEntities = factor(kdsum$NumberOfEntities)

#display.brewer.all()
#display.brewer.pal()

attributes(pdmelt$StDev)
pdmelt <- melt(kdsum, id.vars = c('NumberOfEntities'))
pdf('sim-f-kd.pdf', height = 3, width = 8)
p <- ggplot(pdmelt, aes(x = NumberOfEntities, y = value, group = variable, color = variable)) +
  geom_line(size = 0.9) +
  theme_bw() +
  ggtitle("Impact of partial lists") +
  scale_y_continuous(limits = c(0, 0.25), expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  labs(x = 'Number of Non-overlapping Groups', y = parse(text = 'paste("Kendall\'s ", tau, " Distance")'), color = 'Method') +
  scale_colour_manual(values = brewer.pal(10, "Paired")[7:10])
print(p)
dev.off()

#brewer.pal(10, "Paired")[7:10]
