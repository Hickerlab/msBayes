obs<-scan("pods.obs")
write(obs[-c(2:37)], file=paste("pods.obs"),ncol=383)
