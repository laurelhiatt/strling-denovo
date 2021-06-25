library(LearnBayes)
library(tidyverse)
getwd()
setwd('/Users/Quinlan/Documents/Git/strling-denovo') ### to whatever you want
list.files() ###find what we want to use

data <- read.csv("STR-denovothreshold150", stringsAsFactors = TRUE, header = TRUE, sep='\t')
###import what we want to use
head(data) ### ensure what we want to use is what we actually have
colnames(data) ### for upcoming filter

data <- data %>% filter(novel_amp == 'True') #only 'True' expansions matter
kiddeldad <- data$kiddeldad 
kiddelmom <- data$kiddelmom
###making objects of the values in these columns

plot(density(kiddeldad))
plot(density(kiddelmom))
###density plots of differences to get idea of expansion size

variantcountcontrast <- data %>% group_by(sample) %>% ###can also do sample, mutation
  summarise(all_amplifications = sum(novel_amp =='True'), 
            amplification_vs_dad = sum(kiddeldad > 150), ### can change this to whatever threshold
            amplification_vs_mom = sum(kiddelmom > 150),
            mutation = mutation[1]) %>% 
  gather(denovo_origins, counts, all_amplifications:amplification_vs_mom)
###making a new dataframe with sample, mutations, denovo origins, and counts
### can also set threshold of sum integer with something like threshold = seq(0,max_diff, by = 100)
head(variantcountcontrast) #wow, it looks perfect (or not)
ggplot(variantcountcontrast) + geom_jitter(aes(x= denovo_origins, y = counts, color = mutation), width = 0.3) +ylab('counts filtered at 150bp') + xlab('denovo amplifications')
###jitterplot with counts over threshold by origin, colored by phenotype.


