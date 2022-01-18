#integrate DEG list from RNA-seq with differentially 'bound' genes identified via ChIP-seq

#R packages - 'ComplexUpset' and 'ggplot2'

#Reference
#https://krassowski.github.io/complex-upset/articles/Examples_R.html

#load libraries
library(ggplot2)
library(ComplexUpset)

#read in data 
data.upset <- read.csv("wm_h3k4me3_up_matrix.csv", header = TRUE)
data.upset
#define groups
groups=colnames(data.upset)[3:4]
groups
data.upset[groups] = data.upset[groups] == 1
t(head(data.upset[groups], 23000))

#generate upset plot
upset(data.upset, groups, name='groups', width_ratio=0.4, height_ratio = 0.6, 
      max_size=23000, mode = "inclusive_intersection", 
      sort_sets='ascending', n_intersections=50, min_size=250,
      set_sizes=upset_set_size()+geom_text(aes(label=..count..), stat='count',))

