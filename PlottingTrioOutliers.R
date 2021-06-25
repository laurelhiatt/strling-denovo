library(LearnBayes)
library(tidyverse)
library(ggplot2)

getwd() #where are we
setwd('/Users/Quinlan/Documents/Git/strling-denovo') #where do we want to be
list.files() #what do we want to use

data <- read.csv("STR-denovothreshold150", stringsAsFactors = TRUE, header = TRUE, sep='\t')
head(data)
colnames(data) #to filter by
true_data <- data %>% filter(novel_amp == 'True')
### I only want the implificated amplifications here, so I filter by true
head(true_data)

ggplot(true_data,  aes(x = locus, y= kiddeldad)) + geom_point()
ggplot(true_data,  aes(x = locus, y= kiddelmom)) + geom_point()
### looking at distribution by locus

ggplot(true_data,  aes(x = kiddeldad)) + geom_histogram() + scale_y_continuous(trans = "log10") + xlab("kid allele - dad allele")
ggplot(true_data,  aes(x = kiddelmom)) + geom_histogram() + scale_y_continuous(trans = "log10") + xlab("kid allele - mom allele")
###histograms to show general count of expansion sizes  
###separate by sample with facet wrap

ggplot(true_data) + 
  geom_density(aes(x = kiddeldad, color = 'dad')) +
  geom_density(aes(x = kiddelmom, color = 'mom')) +
  scale_x_continuous(name = 'kid allele - dad/mom alleles, log10', trans = "log10") + facet_wrap(~sample)
###mom and dad density by sample
ggplot(true_data) + 
  geom_density(aes(x = kiddeldad, color = factor(sample))) +
  scale_x_continuous(name = 'kid allele - dad allele by sample, log10',trans = "log10")
### dad density, colored by sample

ggplot(true_data) + 
  geom_density(aes(x = kiddelmom, color = factor(sample))) +
  scale_x_continuous(name = 'kid allele - mom allele by sample, log10',trans = "log10")
### mom density, colored by sample

variantcountcontrast <- true_data %>% group_by(sample) %>% 
  summarise(all_amplifications = sum(novel_amp =='True'), 
            amplification_vs_dad = sum(kiddeldad > 150), 
            amplification_vs_mom = sum(kiddelmom > 150)) %>% 
            gather(denovo_origins, counts, all_amplifications:amplification_vs_mom)
###making a new dataframe based on these conditions

amplification_vs_dad = (true_data$kiddeldad)
amplification_vs_dad
### convert to long format, which is what ggplot likes

ggplot(variantcountcontrast) + 
  geom_jitter(aes(x= denovo_origins, y= counts, color = denovo_origins), height = 0, width =0.2)
###jitter plot by origins

###Below is self reference R script
###threshold =  seq(0,450, by = 50)
###all_sum_over_threshold = data.frame()
###for (this_sample in unique(data$sample)){
  ###sample_df = subset(data, sample == this_sample)
  ###sum_over_threshold_mom = sapply(threshold, function(x){sum(sample_df$kiddelmom >= x)})
  ###sum_over_threshold_dad = sapply(threshold, function(x){sum(sample_df$kiddeldad >= x)})
  ###all_sum_over_threshold = rbind(all_sum_over_threshold, 
     ###                            data.frame(sample = this_sample,threshold = threshold, 
        ###                                    sum_over_threshold_mom = sum_over_threshold_mom,
           ###                                 sum_over_threshold_dad = sum_over_threshold_dad,
              ###                              mutation = sample_df$mutation[1]))
###}
###sum_over_threshold = sapply(threshold, function(x){sum(data >= x)})

