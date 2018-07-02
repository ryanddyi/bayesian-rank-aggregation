
library(ggplot2)
library(RColorBrewer)
library(reshape)


setwd(paste0('/Users/ydd/Documents/BARC/simulation_new/','sim-7'))


M = 69
acc = array(dim = c(100,4))
nCat = array(dim = c(100,4))
weights = array(dim = c(100,M))

for (i in 1:100){
  bat = (i-1) %/% 100 + 1
  folder.char = paste('batch', bat, sep = '')
  accrds = readRDS(paste(folder.char, '/acc', i, '.rds', sep = ''))
  acc[i,] = accrds
  ncatrds = readRDS(paste(folder.char, '/nCat', i, '.rds', sep = ''))
  nCat[i,] = ncatrds
  weight_est = readRDS(paste(folder.char, '/weight', i, '.rds', sep = ''))
  weights[i,] = weight_est
}

colMeans(acc)
colMeans(nCat)
hist(weights,nclass = 40, freq = F, main = '', xlab = 'Weight')

