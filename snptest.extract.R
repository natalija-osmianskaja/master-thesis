args<-commandArgs(trailingOnly=TRUE)
name <- args[1]
n <- as.double(args[2])
p <- as.double(args[3])
data<-read.table(name, header=T, sep=" ")
select = c(1:6,9,42,44:45)
write.table(data[data$frequentist_add_pvalue <= p & data$info >= n & complete.cases(data$frequentist_add_pvalue),select], file= "", quote=F, row.names = F)
