library(TropFishR)
library(tidyverse)
library(magrittr)
library(ggpubr)
Usr=Sys.getenv("USERNAME")
if(Usr=="") Usr='myb'
handl_OneDrive=function(x)paste('C:/Users',Usr,'OneDrive - Department of Primary Industries and Regional Development/Matias',x,sep='/')

hndl.dat=handl_OneDrive('Analyses/Efficiency creep simulation/Data')
setwd(handl_OneDrive('Analyses/Efficiency creep simulation'))

source('Git_Efficiency-creep/FuCreep_wrapper.r')

# PARAMETERS SECTION PSG -------------------------------------------------------------------------
#select method for calculating mnbio 
#bio.method='simple.mean'
bio.method='glm'

# 1. Simple sim data: simulate biomass and cpue trajectory based on simple dynamic model-------------------------------------------------------------------------
data(emperor)
out=prod_mod_ts(emperor, method = "Schaefer")
Years.ori=emperor$year
Q.ori=out$q
BIOMASS.ori=out$Bvec
CPUE.ori=out$CPUE
CPUE_est.ori=out$CPUE_hat


# 2. Simple sim data: simulate fisher specific cpue with efficiency creep -------------------------------------------------------------------------
set.seed(666)
n.fishers=1e1
baseline.eff.creep=0.01
n.qs=rnorm(n.fishers,Q.ori,Q.ori*.1) #;hist(n.qs)
n.eff=rlnorm(n.fishers,log(baseline.eff.creep),.5) #;hist(n.eff)
n.eff[2] <- 0.00004   #must have 1 fisher with very low efficiency

fisher.cpue=data.frame(Years=rep(Years.ori,n.fishers),
                       BIOMASS=rep(BIOMASS.ori,n.fishers),
                       CPUE=rep(CPUE_est.ori,n.fishers),
                       Q=rep(Q.ori,length(BIOMASS.ori)*n.fishers),
                       fisher=paste0('f_',rep(1:n.fishers,each=length(Years.ori))),
                       q.f=rep(n.qs,each=length(Years.ori)),
                       creep=rep(n.eff,each=length(Years.ori)))%>%
                mutate(cpue.f=q.f*BIOMASS.ori,
                       fisher=factor(fisher,levels=paste0('f_',1:n.fishers)))
fisher.cpue=fisher.cpue%>%
              group_by(fisher)%>%
              mutate(f.creep=(1+creep)^(row_number()-1))%>%
              data.frame()%>%
              ungroup()%>%
              mutate(cpue.f.creep=cpue.f*f.creep)

# 3. Simple sim data: display simulated cpues -------------------------------------------------------------------------
setwd(paste0(getwd(),'/Outputs'))
p1=fisher.cpue%>%
  ggplot()+
  geom_point(aes(x=Years,y=CPUE),size=4)+geom_line(aes(x=Years,y=CPUE),col=1,linewidth=1.5)+
  geom_point(aes(Years,cpue.f,color=fisher))+
  geom_line(aes(Years,cpue.f,color=fisher),linetype='longdash')+
  theme_bw()+ggtitle('Mean CPUE (black) & fisher cpue (without creep)')+
  theme(legend.position = 'top')+guides(color=guide_legend(nrow =2))

p2=fisher.cpue%>%
  ggplot()+
  geom_point(aes(Years,f.creep),color='steelblue')+
  geom_line(aes(Years,f.creep),color='steelblue')+
  facet_wrap(~fisher)+
  theme_bw()+ggtitle('Efficiency creep by fisher')

p3=fisher.cpue%>%
  mutate(dummy=rep(CPUE_est.ori,n.fishers))%>%
  ggplot()+
  geom_point(aes(Years,cpue.f),color=2)+
  geom_line(aes(Years,cpue.f),color=2)+
  geom_point(aes(Years,cpue.f.creep),color=3)+
  geom_line(aes(Years,cpue.f.creep),color=3)+
  geom_line(aes(x=Years,y=dummy),col=1,alpha=.6,linewidth=1.25)+
  facet_wrap(~fisher,scales='free_y')+
  theme_bw()+ggtitle('Mean CPUE (black) & fisher cpue (green, with; red, without creep)')

ggarrange(plotlist =list(p1,p2,p3), ncol=3,nrow=1)

ggsave("sim_dat.tiff",width=20,height= 8,compression="lzw")


# 4. Simple sim data: run FuCreep-------------------------------------------------------------------------

  #Matias'
dat.Mat=fisher.cpue%>%
          rename(Year=Years,
                 value=cpue.f.creep)%>%
          mutate(fisher=as.character(fisher),
                 fisher=str_remove(fisher,'_'),
                 Year=Year-(min(Year)-1),
                 year=paste0('y',Year))
original.efficiency=dat.Mat%>%
          distinct(fisher,creep)%>%
          mutate(fisher=as.character(fisher),
                 fisher=str_remove(fisher,'_'))

write.csv(dat.Mat%>%dplyr::select(year,fisher,value,Year),paste0(hndl.dat,'/dat.Mat.csv'),row.names = F)
write.csv(original.efficiency%>%pull(creep),paste0(hndl.dat,'/origeff.Mat.csv'),row.names = F)

res.Matias=fn.FuCreep(dat=dat.Mat%>%dplyr::select(year,fisher,value,Year),
                      Sims=1e2,
                      original.efficiency=original.efficiency,
                      get.bio=bio.method,
                      CPUE_ori=CPUE_est.ori)

ggarrange(plotlist =list(res.Matias$p1,res.Matias$p2,res.Matias$p3), ncol=1,nrow=3)
ggsave("sim_dat_estimates.tiff",width=8,height= 8,compression="lzw")

#Drop some years for some fishers
do.this=FALSE
if(do.this)
{
  DROP=dat.Mat%>%
    distinct(Year,fisher)%>%
    mutate(Yr.fisher=paste(Year,fisher))
  
  Prop.dropped=0.1
  drop.id=match(sample(DROP$Yr.fisher,round(nrow(DROP)*Prop.dropped),replace=FALSE),DROP$Yr.fisher)
  res.Matias.unbalanced1=fn.FuCreep(dat=dat.Mat[-drop.id,]%>%
                                      dplyr::select(year,fisher,value,Year),
                                    Sims=1e2,
                                    original.efficiency=original.efficiency,
                                    get.bio=bio.method,
                                    CPUE_ori=CPUE_est.ori)
  subtit1=paste('Matias',100*Prop.dropped,'% less observations')
  
  Prop.dropped=0.5
  drop.id=match(sample(DROP$Yr.fisher,round(nrow(DROP)*Prop.dropped),replace=FALSE),DROP$Yr.fisher)
  res.Matias.unbalanced2=fn.FuCreep(dat=dat.Mat[-drop.id,]%>%
                                      dplyr::select(year,fisher,value,Year),
                                    Sims=1e2,
                                    original.efficiency=original.efficiency,
                                    get.bio=bio.method,
                                    CPUE_ori=CPUE_est.ori)
  subtit2=paste('Matias',100*Prop.dropped,'% less observations')
  
  p.list1=ggarrange(plotlist=list(res.Matias$p2+
                                    labs(title='Observed (black) and predicted creep by fisher',
                                         subtitle = 'Matias'),
                                  res.Matias.unbalanced1$p2+labs(title='',subtitle=subtit1),
                                  res.Matias.unbalanced2$p2+labs(title='',subtitle=subtit2)),
                    nrow=3,ncol=1)
  p.list2=ggarrange(plotlist=list(res.Matias$p3+
                                    labs(title='Observed and predicted cpue',
                                         subtitle = 'Matias'),
                                  res.Matias.unbalanced1$p3+labs(title=''),
                                  res.Matias.unbalanced2$p3+labs(title='')),
                    nrow=3,ncol=1,common.legend = TRUE)
  ggarrange(p.list1,p.list2)
  ggsave("Compare Matias Matias unbalanced_estimated efficiency.tiff",width=10,height= 8,compression="lzw")
  
}

#compare Simon and Matias data
do.this=FALSE
if(do.this)
{
  #Simon's
  dat.Simon=read.csv(paste0(hndl.dat,'/30by30_CPUEData.csv'))
  res.Simon=fn.FuCreep(dat=dat.Simon,
                       Sims=100,
                       original.efficiency=data.frame(fisher=unique(dat.Simon$fisher),
                                                      creep=dget(paste0(hndl.dat,'/origeff'))),
                       get.bio=bio.method)
  
  ggarrange(plotlist =list(res.Simon$p1+
                             labs(title='Observed cpue by fisher and average cpue (grey)',subtitle = 'Simon'),
                           res.Matias$p1+labs(title='',subtitle='Matias')),ncol=1,nrow=2)
  ggsave("Compare Simon Matias_cpue by fisher.tiff",width=10,height= 8,compression="lzw")
  
  
  ggarrange(plotlist =list(res.Simon$p2+
                             labs(title='Observed (black) and predicted efficiency creep by fisher',subtitle = 'Simon'),
                           res.Matias$p2+labs(title='',subtitle='Matias')),ncol=1,nrow=2)
  ggsave("Compare Simon Matias_estimated efficiency.tiff",width=10,height= 8,compression="lzw")
  
}


# 5. Realistic sim data: run FuCreep -------------------------------------
  #5.1. Create fisher cpue
dat.sim.stand=fn.sim.stand.input(Biom=BIOMASS.ori,
                                 nblock = 2,
                                 nmonth = 12, 
                                 ndepth = 2, 
                                 nyear = length(Years.ori), 
                                 nfish = 10,
                                 dist='normal')  #dist='unif' 'normal'

  #5.2. Run glm of dat.sim.stand to get annual cpue by fisher
fisher.cpue1=fn.glm1(dd=dat.sim.stand%>%rename(cpue=cpue.f.creep))

  #5.3. Display 
display.sim=TRUE
if(display.sim)
{
  p1=fisher.cpue1%>%
    rename(Years=year,
           fisher=fish)%>%
    mutate(dummy=rep(CPUE_est.ori,length(unique(fisher.cpue1$fish))))%>%
    ggplot()+
    geom_point(aes(Years,cpue,color=fisher))+
    geom_line(aes(Years,cpue,color=fisher),linetype='longdash')+
    theme_bw()+ggtitle('Mean CPUE (black) & fisher cpue (without creep)')+
    theme(legend.position = 'top')+guides(color=guide_legend(nrow =2))+
    geom_line(aes(x=Years,y=dummy),col=1,alpha=.6,linewidth=1.25)
  
  p2=dat.sim.stand%>%
    rename(Years=year,
           fisher=fish)%>%
    group_by(Years,fisher)%>%
    summarise(f.creep=mean(f.creep),
              creep=mean(creep))%>%
    ungroup()%>%
    mutate(labl=paste0(fisher,' creep=',round(creep,3)))%>%
    ggplot()+
    geom_point(aes(Years,f.creep),color='steelblue')+
    geom_line(aes(Years,f.creep),color='steelblue')+
    facet_wrap(~labl)+
    theme_bw()+ggtitle('Efficiency creep by fisher')
  
  p3=dat.sim.stand%>%
    rename(Years=year,
           fisher=fish)%>%
    group_by(Years,fisher)%>%
    summarise(cpue.f.creep=mean(cpue.f.creep),
              cpue.f=mean(cpue.f))%>%
    ungroup()%>%
    arrange(fisher,Years)%>%
    mutate(dummy=rep(CPUE_est.ori,length(unique(dat.sim.stand$fish))))%>%
    ggplot()+
    geom_point(aes(Years,cpue.f),color=2)+
    geom_line(aes(Years,cpue.f),color=2)+
    geom_point(aes(Years,cpue.f.creep),color=3)+
    geom_line(aes(Years,cpue.f.creep),color=3)+
    geom_line(aes(x=Years,y=dummy),col=1,alpha=.6,linewidth=1.25)+
    facet_wrap(~fisher,scales='free_y')+
    theme_bw()+ggtitle('Mean CPUE (black) & fisher cpue (green, with; red, without creep)')
  
  ggarrange(plotlist =list(p1,p2,p3), ncol=3,nrow=1)
  ggsave("sim_dat_realistic.tiff",width=20,height= 8,compression="lzw")
  
}

  #5.4. run FuCreep   #ACA. Not working, not estimating creep coefficients...
dat.Mat=fisher.cpue1%>%
          rename(Year=year,
                 value=cpue,
                 fisher=fish)%>%
          mutate(fisher=as.character(fisher),
                 fisher=str_remove(fisher,'_'),
                 year=paste0('y',Year))%>%
          mutate(fisher=paste0('f',fisher))
original.efficiency=dat.sim.stand%>%
          distinct(fish,creep)%>%
          rename(fisher=fish)%>%
          mutate(fisher=paste0('f',fisher))

res.Matias=fn.FuCreep(dat=dat.Mat,
                      Sims=1e2,
                      original.efficiency=original.efficiency,
                      get.bio=bio.method,
                      CPUE_ori=CPUE_est.ori)

ggarrange(plotlist =list(res.Matias$p1,res.Matias$p2,res.Matias$p3), ncol=1,nrow=3)
ggsave("sim_dat_realistic_estimates.tiff",width=8,height= 8,compression="lzw")
