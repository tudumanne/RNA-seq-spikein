#integrate DEG list from RNA-seq with differentially 'bound' genes identified via ChIP-seq

#generate a matrix 

#read in data 
#e.g. WT-Mal, up-regulated genes with H3K4me3 differentially bound genes
data = read.csv("wm_h3k4me3_up.csv", header = TRUE)

#create a data matrix, with 1 and O's
data1 = c(data$DEG.down)
data2 = c(data$H3K4me3.gain)
data3 = c(data$H3K4me3.loss)

dat <- paste0("data", 1:3)
out <- t(splitstackshape:::charMat(listOfValues = mget(dat), fill = 0L))
colnames(out) <- dat
out

#export the matrix as a .csv file
write.csv(out, file="wm_all_down_matrix.csv")


