library(tidyverse)
library(magrittr)
library(tictoc)

  #select dataset
test.what='Simon'
#test.what='Matias'
if(test.what=='Simon')
{
  dat <- read.csv('30by30_CPUEData.csv')
  original.efficiency=dget('origeff')
}
if(test.what=='Matias')
{
  dat <- read.csv('dat.Mat.csv')
  original.efficiency=read.csv('origeff.Mat.csv')$x
}
head(dat)

  #select method for calculating mnbio 
#get.bio='simple.mean'
get.bio='glm'

fn.glm=function(dd)
{
  mod=glm(value~Year,data=dd%>%mutate(Year=factor(Year)))
  newdata = dd%>%distinct(Year)%>%mutate(Year=factor(Year))
  newdata=newdata%>%
            mutate(mnbio=predict(mod,newdata,type='response'),
                   Year=as.integer(as.character(Year)))
  return(newdata)
}

## First loop
#Calc bad average biomass
if(get.bio=='simple.mean') dat2 <- dat %>% group_by(Year) %>% summarise(mnbio=mean(value)) 
if(get.bio=='glm') dat2 <- fn.glm(dd=dat)
p=dat%>%
  ggplot()+
  geom_line(aes(Year,value,color=fisher))+
  geom_point(data=dat%>%filter(fisher%in%c('f1','f10')),aes(Year,value,color=fisher))+
  geom_line(data=dat2%>%mutate(fisher='average'),aes(x=Year,y=mnbio),linewidth=3,alpha=.3,color='black')+
  theme_bw()


## calc q for each fisher
dat %<>% full_join(dat2) %>% mutate(baseq=value/(mnbio))

## Pull out just the first year to use for all years (calc biomass to fisher)
for(f in 1:length(unique(dat$fisher)))
{
  dat$baseq[dat$fisher==unique(dat$fisher)[f]] <- dat$baseq[dat$fisher==unique(dat$fisher)[f] & dat$Year==1]
}

## Add empty slots
dat$Anncreep <- dat$Estcreep <-  NA

## Function to estimate efficiency to match cpue to biomass 
func1 <- function(x,y)
{
  estcpue <- y$mnbio*y$baseq
  est <- estcpue * (1+x)^(0:(length(y$baseq)-1))
  SS <- sum((y$val - est)^2)
  return(SS)
} 

## Apply to each fisher
for(f in 1:length(unique(dat$fisher)))
{
  tmp <- dat[dat$fisher==unique(dat$fisher)[f],] 
  tmp$Anncreep <- optimise(func1, interval = c(0.000001, 0.10), y=tmp)$minimum  
  tmp$Estcreep <- (1+tmp$Anncreep)^(0:(nrow(tmp)-1))
  dat[dat$fisher==unique(dat$fisher)[f],] <- tmp
}


## make estimated cpue based on efficiency creep
dat %<>% mutate(estcpue=value/Estcreep)

## Do big loop to continue to get a better estimate of average biomass after taking into account he estimated efficiency creep.
## Had to iterate many times to slowly move down
#Calc biomass
tic()
for(l in 1:100)
{
  if(get.bio=='simple.mean') dat2 <- dat %>% group_by(Year) %>% summarise(mnbio=mean(estcpue))  
  if(get.bio=='glm') dat2 <- fn.glm(dd=dat%>%dplyr::select(Year,estcpue)%>%rename(value=estcpue))
  
  ## calc q for each fisher
  dat %<>%  select(-mnbio) %>% full_join(dat2) %>% mutate(baseq=value/(mnbio))

  ## Pull out just the first year
  for(f in 1:length(unique(dat$fisher)))
  {
    dat$baseq[dat$fisher==unique(dat$fisher)[f]] <- dat$baseq[dat$fisher==unique(dat$fisher)[f] & dat$Year==1]
  }
  
  for(f in 1:length(unique(dat$fisher)))
  {
    tmp <- dat[dat$fisher==unique(dat$fisher)[f],] 
    tmp$Anncreep <- optimise(func1, interval = c(0.000001, 0.10), y=tmp)$minimum  
    tmp$Estcreep <- (1+tmp$Anncreep)^(0:(nrow(tmp)-1))
    dat[dat$fisher==unique(dat$fisher)[f],] <- tmp
  }
  
  dat %<>% mutate(estcpue=value/Estcreep)
  
}
toc()
  
## Make output and compare to values used to make simulated data
out <- dat %>% mutate(fish=as.numeric(gsub('f','',fisher))) %>% group_by(fish) %>% summarise(effic=mean(Anncreep)) %>% as.data.frame()

out$orig <- original.efficiency

plot(out$fish, out$effic, type='o')  
lines(out$fish, out$orig, col=2)  
  
  
  
#out.simple=out  
  