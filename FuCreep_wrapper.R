library(emmeans)


# Simulate cpue data-------------------------------------------------------------------------
fn.sim.stand.input=function(Biom,nblock,nmonth,ndepth,nyear,nfish,dist)
{
  dat <- expand.grid(block=1:nblock,year=1:nyear,month=1:nmonth,depth=1:ndepth,fish=1:nfish)
  qblock <- runif(nblock,-0.00001,0.00001)
  qmonth <- runif(nmonth,-0.00001,0.00001)
  qdepth <- runif(ndepth,-0.00001,0.00001)
  if(dist=='unif')
  {
    qfish <- runif(nfish,0.0001,0.0003)
    fcrep <- runif(nfish,1e-4,4e-2)
  }
  if(dist=='normal')
  {
    qfish=rnorm(nfish,Q.ori,Q.ori*.1)
    fcrep=rlnorm(nfish,log(baseline.eff.creep),.5)
  }
  fcrep[2]=0.00004  #must set one to very low
  dat$block <- qblock[dat$block]
  dat$month <- qmonth[dat$month]
  dat$depth <- qdepth[dat$depth]
  dat$creep <- fcrep[dat$fish]
  dat$f.q <- qfish[dat$fish]
  
  dat$biom <- Biom[dat$year]
  dat=dat%>%
    mutate(q=f.q+block+month+depth,
           cpue.f=q*biom)%>%
    group_by(fish)%>%
    mutate(f.creep=(1+creep)^(year-1))%>%
    data.frame()%>%
    ungroup()%>%
    mutate(cpue.f.creep=cpue.f*f.creep)
  return(dat)
}

# Calculate mean cpue by fisher-------------------------------------------------------------------------
fn.glm1=function(dd)
{
  dd=dd%>%
        mutate(fish=factor(fish),
               year=factor(year),
               month=factor(month),
               block=factor(block))
   mod=glm(log(cpue)~year+month+block+depth+fish,data=dd)
  newdata=summary(emmeans (mod, ~ year+fish),type='response')%>%
      data.frame()%>%
      rename(cpue=response)%>%
      dplyr::select(year,fish,cpue)%>%
      mutate(year=as.integer(as.character(year)))
  return(newdata)
}

# Estimate efficiency to match cpue to biomass-------------------------------------------------------------------------
func1 <- function(x,y)
{
  estcpue <- y$mnbio*y$baseq
  est <- estcpue * (1+x)^(y$Year-1)
  #est <- estcpue * (1+x)^(0:(length(y$baseq)-1))
  SS <- sum((y$value - est)^2)
  #SS <- sum((y$val - est)^2)
  return(SS)
}

# Calculate mean mnbio-------------------------------------------------------------------------
fn.glm=function(dd,what.mod='Year.vessel',what.pred='emmeans')
{
  dd=dd%>%
        mutate(fisher=factor(fisher),
               Year=factor(Year))
  if(what.mod=='Year') mod=glm(value~Year,data=dd)
  if(what.mod=='Year.vessel') mod=glm(log(value)~Year+fisher,data=dd)
  
  if(what.pred=='emmeans')
  {
    newdata=summary(emmeans (mod, ~ Year),type='response')%>%
      data.frame()%>%
      rename(mnbio=response)%>%
      dplyr::select(Year, mnbio)%>%
      mutate(Year=as.integer(as.character(Year)))
  }
  if(what.pred=='predict')
  {
    newdata = dd%>%distinct(Year)%>%mutate(Year=factor(Year))
     newdata=newdata%>%
       mutate(mnbio=predict(mod,newdata,type='response'),
              Year=as.integer(as.character(Year)))
  }
  
  return(newdata)
}

# FuCreep wrapper function-------------------------------------------------------------------------
fn.FuCreep=function(dat,Sims=100,original.efficiency,get.bio,CPUE_ori=NULL)
{
  ## First loop
  #Calc bad average biomass
  if(get.bio=='simple.mean') dat2 <- dat %>% group_by(Year) %>% summarise(mnbio=mean(value)) 
  if(get.bio=='glm') dat2 <- fn.glm(dd=dat)
  dat2.ori=dat2
  
  #display fisher cpue and mean
  p1=dat%>%
    ggplot()+
    geom_line(aes(Year,value,color=fisher))+
     geom_line(data=dat2%>%mutate(fisher='average'),aes(x=Year,y=mnbio),linewidth=3,alpha=.3,color='black')+
    theme_bw()+
    ggtitle('Observed cpue by fisher and average cpue (grey)')
  
  
  ## calc q for each fisher
  dat %<>% full_join(dat2,join_by('Year')) %>% mutate(baseq=value/(mnbio))
  
  ## Pull out just the first year to use for all years (calc biomass to fisher)
  Unik.f=unique(dat$fisher)
  for(f in 1:length(Unik.f))
  {
    Min.yr=min(dat%>%filter(fisher==Unik.f[f])%>%pull(Year))
    dat$baseq[dat$fisher==Unik.f[f]] <- dat$baseq[dat$fisher==Unik.f[f] & dat$Year==Min.yr]
  }
  
  ## Add empty slots
  dat$Anncreep <- dat$Estcreep <-  NA
  

  ## Apply to each fisher
  for(f in 1:length(Unik.f))
  {
    tmp <- dat[dat$fisher==Unik.f[f],] 
    tmp$Anncreep <- optimise(func1, interval = c(1e-06, 0.10), y=tmp)$minimum  
    tmp$Estcreep <- (1+tmp$Anncreep)^(tmp$Year-1)
    #tmp$Estcreep <- (1+tmp$Anncreep)^(0:(nrow(tmp)-1))
    dat[dat$fisher==Unik.f[f],] <- tmp
  }
  
  
  ## make estimated cpue based on efficiency creep
  dat %<>% mutate(estcpue=value/Estcreep)
  
  ## Do big loop to continue to get a better estimate of average biomass after taking into account he estimated efficiency creep.
  ## Had to iterate many times to slowly move down
  #Calc biomass
  for(l in 1:Sims)
  {
    if(get.bio=='simple.mean') dat2 <- dat %>% group_by(Year) %>% summarise(mnbio=mean(estcpue))  
    if(get.bio=='glm') dat2 <- fn.glm(dd=dat%>%dplyr::select(Year,fisher,estcpue)%>%rename(value=estcpue))
    
    ## calc q for each fisher
    dat %<>%  select(-mnbio) %>% full_join(dat2,join_by('Year')) %>% mutate(baseq=value/(mnbio))
    
    ## Pull out just the first year
    for(f in 1:length(Unik.f))
    {
      Min.yr=min(dat%>%filter(fisher==Unik.f[f])%>%pull(Year))
      dat$baseq[dat$fisher==Unik.f[f]] <- dat$baseq[dat$fisher==Unik.f[f] & dat$Year==Min.yr]
    }
    
    for(f in 1:length(Unik.f))
    {
      tmp <- dat[dat$fisher==Unik.f[f],] 
      tmp$Anncreep <- optimise(func1, interval = c(0.000001, 0.10), y=tmp)$minimum  
      tmp$Estcreep <- (1+tmp$Anncreep)^(tmp$Year-1)
      #tmp$Estcreep <- (1+tmp$Anncreep)^(0:(nrow(tmp)-1))
      dat[dat$fisher==Unik.f[f],] <- tmp
    }
    
    dat %<>% mutate(estcpue=value/Estcreep)
    
  }
  
  
  ## Make output and compare to values used to make simulated data
  out <- dat %>%
            group_by(fisher) %>%
            summarise(effic=mean(Anncreep)) %>%
            as.data.frame() %>%
            left_join(original.efficiency,by='fisher')%>%
            mutate(fisher=as.numeric(gsub(".*?([0-9]+).*", "\\1", fisher))) 
  
  #out$orig <- original.efficiency
  p2=out%>%
    rename(effic_original=creep)%>%
    ggplot(aes(fisher,effic))+
    geom_point(color='red',size=3)+
    geom_line(aes(),color='red',size=1.1,linetype='longdash')+
    geom_point(aes(fisher,effic_original),size=2)+
    theme_bw()+ggtitle('Observed (black) and predicted (red) efficiency creep')
  
  #compare estimated cpue with original cpue
  p3=NULL
  if(!is.null(CPUE_ori))
  {
    if(get.bio=='simple.mean') dat2 <- dat %>% group_by(Year) %>% summarise(mnbio=mean(estcpue))  
    if(get.bio=='glm') dat2 <- fn.glm(dd=dat%>%dplyr::select(Year,fisher,estcpue)%>%rename(value=estcpue))
    p3=rbind(dat2%>%mutate(data='creep corrected'),
             data.frame(Year=dat2$Year,mnbio=CPUE_ori,data='original'),
             dat2.ori%>%mutate(data='creep not corrected'))%>%
      ggplot(aes(Year,mnbio,color=data))+
      geom_point(size=2.5)+
      geom_line(linetype='longdash')+ylab('CPUE')+
      theme_bw()+ggtitle('Observed and predicted cpue')+theme(legend.title = element_blank())
  }

  
  return(list(p1=p1,p2=p2,p3=p3))
}

