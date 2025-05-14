biodel <- 0.05
biom <- rep(10000,30)
set.seed(7)
for(y in 2:30){
  biom[y] <-  biom[y-1] * (1+rnorm(1,0,biodel))
}

plot(biom, type='o')

fisU <- matrix(0, ncol=30, nrow=30, dimnames = list(year=paste0('y',1:30),fisher=paste0('f',1:30)))
set.seed(7)
feff <- rnorm(30,0,1)
feff <- feff+abs(min(feff))
feff <- 1+ feff/100

set.seed(17)
fq <- rnorm(30,0,20)
fq <- fq+abs(min(fq))+4
fq <- fq/1000

for(f in 1:30){
  fisU[,f] <-  feff[f]^(1:30)*fq[f] * biom
}

dat <- reshape2::melt(fisU) %>% mutate(Year=as.numeric(gsub('y','',year)))

ggplot2::ggplot(dat,aes(x=Year,y=value,colour=fisher))+
  geom_line()
write.csv(dat,'30by30_CPUEData.csv',row.names = F)
getwd()
