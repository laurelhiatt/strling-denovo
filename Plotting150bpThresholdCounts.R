###plot the number of expansions passing threshold for each individual. 
###Then color that by phenotype

library(LearnBayes)
library(tidyverse)
library(ggplot2)

getwd()
##setwd() if necessary
list.files()
###make sure what we want is where we want it

data <- read.csv("STR-denovothreshold150", stringsAsFactors = TRUE, header = TRUE, sep='\t')
###inputdata
head(data) ### is this our data?
colnames(data) ### we're going to work with column names so
true_data <- data %>% filter(novel_amp == 'True')
### we have now subset the data to where we believe to have de novo expansions
### or, where novel_amp is true

ind <- aggregate(novel_amp ~ sample, true_data, length)
###we are taking the count of novel_amps by sample and making a new data frame
colnames(ind) #again, want to know what we're #working with when we graph
ggplot(ind) + geom_col(aes(x = factor(sample), y = novel_amp))
### this is an easy way to see the counts of novel amplifications by sample


variantcountcontrast <- data %>% group_by(sample) %>% ###can also do sample, mutation
  summarise(all_amplifications = sum(novel_amp =='True'), 
            amplification_vs_mom = sum(kiddelmom > 150),
            amplification_vs_dad = sum(kiddeldad > 150), 
            mutation = mutation[1]) %>% 
  gather(denovo_origins, counts, all_amplifications:amplification_vs_dad)
#This is another way of getting where novel_amp is true, counting where it meets the 150bp threshold
###note: this is already the case with this data but we can alter the threshold here if we would like
### and then including the mutation data so we can color by mutation
head(variantcountcontrast) ### do we have what we want in this object
ggplot(variantcountcontrast) + geom_jitter(aes(x= denovo_origins, y = counts, color = mutation), width = 0.3) +ylab('counts filtered at 150bp') + xlab('denovo amplifications')
###now, we can see the general count by mutation for all amplifications and then those seen compared to dad versus mom.




ggplot(variantcountcontrast) + geom_jitter(aes(x= denovo_origins, y = counts, color = mutation), width = 0.3) +ylab('counts filtered at 150bp') + xlab('denovo amplifications') + facet_wrap(variantcountcontrast$sample)
###If I want this split up by sample, I can put that here