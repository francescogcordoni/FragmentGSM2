library(knitr)
library(egg)
library(MASS)
library(rgl)
library(zoo)
library(Sim.DiffProc)
library(doParallel)
library(nlme)
library(data.table)
library(plotly)
library(ggsci)
library(tidymodels)
library(ggplot2)
library(dplyr)
library(readxl)
library(tidymodels)
#library(rootSolve)

cb_a <- c("#E69F00","#000000","#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cb_b <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cb_c <- c( "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7")

cb_nob <- c("#E69F00", "#56B4E9","#6A3D9A","darkred","#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
c25 <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "black", "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown"
)

path_or<-"/home/user/Scrivania/Dottorato/RISULTATI_FIT/RBE"
setwd(path_or)
source("/home/user/Scrivania/TESI/Codici_R/utilities_doserate.R")

cell<-"H460"
#cell<-"H1437"
ion<-"He"
ion<-"H"
ion<-"C"

if(ion=="He"){
  load("/home/user/Scrivania/Dottorato/RISULTATI_FIT/RBE/CSV_RDATA/Surv_mod_He.RData")
}else if(ion=="H"){
  load("/home/user/Scrivania/Dottorato/RISULTATI_FIT/RBE/CSV_RDATA/Surv_mod_H.RData")
}else if(ion=="C"){
  load("/home/user/Scrivania/Dottorato/RISULTATI_FIT/RBE/CSV_RDATA/Surv_mod_C.RData")
}

load("/home/user/Scrivania/Dottorato/Codici_R/Survival_uniform_dose_sampling.RData")
#load("/home/user/Scrivania/Dottorato/Codici_R/Survival_uniform_dose_sampling_H1437.RData")
pide_data<-survival

#surv<-na.omit(surv_GSM2[[paste0(as.character(rd),"-",as.character(Rn))]])

# Carico dosi
setwd("/home/user/Scrivania/Dottorato/RISULTATI_FIT")

# LET values
LET_H <- c("1","2.6","4.7","7.3","8.7","11.1","13.7","15.4","16.9","18.3","20.2","21.4")
LET_He <- c("3.4","8.5","12.7","22.4","36.2","44.9","54.7","63","71.1","77.6","83.6","88.7")
LET_C <- c("20.2","39.8","63.1","70.6","84.3","100.8","126.3","157.6","196.4","242.9","285.8","308.4")

# yD values (R=0.8um)
yD_H<-c(7.36,7.40,7.68,10.6,12.2,14.9,18,19.8,21.5,22.9,25.1,26.7)
yD_He<-c(12,17.5,21.8,32.1,47.4,56.3,66.4,74.2,82.3,88.2,93,93.9)
yD_C<-c(27.4,52.9,80.2,84.6,91.1,98.8,112.2,137.4,172,235.1,292.2,319.6)

if(ion=="He"){
  LET_cycle <- LET_He[c(1,2,3,4,5,6,7,8,9,10,11,12)]
}else if(ion=="H"){
  LET_cycle <- LET_H[c(1,2,3,4,5,6,7,8,9,10,11,12)]
}else if(ion=="C"){
  LET_cycle <- LET_C[c(1,2,3,4,5,6,7,8,9,10,11,12)]
}

cell_line <-"H460"
#cell_line <-"H1437"
# Choose cell line
if(cell_line == "H1437"){
  ax <- 0.05; bx <- 0.041
}else if (cell_line == "H460"){
  ax <- 0.29; bx <- 0.083
}

# LQ H460
D10_Xray <- seq(0,10,by=0.1)[which.min(abs(exp(-0.29*seq(0,10,by=0.1)-0.083*seq(0,10,by=0.1)*seq(0,10,by=0.1))-0.1))]
D50_Xray <- seq(0,10,by=0.1)[which.min(abs(exp(-0.29*seq(0,10,by=0.1)-0.083*seq(0,10,by=0.1)*seq(0,10,by=0.1))-0.50))]
D90_Xray <- seq(0,10,by=0.1)[which.min(abs(exp(-0.29*seq(0,10,by=0.1)-0.083*seq(0,10,by=0.1)*seq(0,10,by=0.1))-0.90))]

# LQ H1437
#D10_Xray <- seq(0,10,by=0.1)[which.min(abs(exp(-ax*seq(0,10,by=0.1)-bx*seq(0,10,by=0.1)*seq(0,10,by=0.1))-0.1))]

# LQL
# seq(0,10,by=0.1)[which.min(abs(exp(-0.18*seq(0,10,by=0.1)-0.1*seq(0,10,by=0.1)*seq(0,10,by=0.1))-0.1))]

quadraticRoots <- function(a, b, c) {
  
  #print(paste0("You have chosen the quadratic equation ", a, "x^2 + ", b, "x + ", c, "."))
  
  discriminant <- (b^2) - (4*a*c)
  
  if(discriminant < 0) {
    return(paste0("This quadratic equation has no real numbered roots."))
  }
  else if(discriminant > 0) {
    x_int_plus <- (-b + sqrt(discriminant)) / (2*a)
    x_int_neg <- (-b - sqrt(discriminant)) / (2*a)
    
    return(x_int_plus)
  }
  else #discriminant = 0  case
    x_int <- (-b) / (2*a)
  return(x_int)
}

survperc<-0.1
integral_spacing<-0.001 # MID CALCULATIONS
{
Sreal<-c()
surv_mod<-c()
dose<-c()
err<-c()
LET<-c()
yD<-c()
alpha_exp<-c()
alpha_up<-c()
alpha_down<-c()
alpha_mod<-c()
beta_exp<-c()
beta_up<-c()
beta_down<-c()
beta_mod<-c()
d10_exp<-c()
d10_up<-c()
d10_down<-c()
d10_mod<-c()
MID_exp<-c()
MID_mod<-c()

ctr<-1
for(file in LET_cycle){
  if(ion=="He"){
    S_real<-pide_data$H460$He[[file]]$LQ
    dose<-pide_data$H460$He[[file]]$Dose
    err<-pide_data$H460$He[[file]]$ErrSF
    S_up<-S_real+0.5*err
    S_down<-S_real-0.5*err
  }else if(ion=="H"){
    S_real<-pide_data$H460$H[[file]]$LQ
    dose<-pide_data$H460$H[[file]]$Dose
    err<-pide_data$H460$H[[file]]$ErrSF
    S_up<-S_real+0.5*err
    S_down<-S_real-0.5*err
  }else if(ion=="C"){
    S_real<-pide_data$H460$C[[file]]$LQ
    dose<-pide_data$H460$C[[file]]$Dose
    err<-pide_data$H460$C[[file]]$ErrSF
    S_up<-S_real+0.5*err
    S_down<-S_real-0.5*err
  }
  
  surv_mod<-surv[c((10*(ctr-1)+1):(10*ctr))]
  
  # FIT for ERROR
  if(ion == "H"){
    quadratic_fit_exp <- LET_cycle
    quadratic_fit_mod <- LET_cycle
  }
  if(ion == "He"){
    quadratic_fit_exp <- LET_cycle[-c(10,11,12)]
    quadratic_fit_mod <- LET_cycle[-c(4,5,6,7,8,9,10,11,12)]
  }
  if(ion == "C"){
    quadratic_fit_exp <- LET_cycle[-c(4,5,6,7,8,9,10,11,12)]
    quadratic_fit_mod <- LET_cycle[c(1,2)]
  }
  
  ###### FIT (linear or LQ) ###########
  if(file %in% quadratic_fit_exp){
    dose2<-dose*dose
    s<-log(S_real)
    mod <- lm(s ~ dose+dose2-1)
    s_up<-log(S_up)
    mod_up <- lm(s_up ~ dose+dose2-1)
    s_down<-log(S_down)
    mod_down <- lm(s_down ~ dose+dose2-1)
  }else{
    s<-log(S_real)
    mod <- lm(s ~ dose-1, weights = 1/sqrt(err))
    s_up<-log(S_up)
    mod_up <- lm(s_up ~ dose-1)
    s_down<-log(S_down)
    mod_down <- lm(s_down ~ dose-1)
  }
  if(file %in% quadratic_fit_mod){
    dose2<-dose*dose
    s1<-log(surv_mod)
    mod1 <- lm(s1 ~ dose+dose2-1)
  }else{
    s1<-log(surv_mod)
    mod1 <- lm(s1 ~ dose-1)
  }
  alpha_exp<-c(alpha_exp,mod$coefficients[1])
  is.na(alpha_exp)<-0
  alpha_mod<-c(alpha_mod,mod1$coefficients[1])
  is.na(alpha_mod)<-0
  beta_exp<-c(beta_exp,mod$coefficients[2])
  is.na(beta_exp)<-0
  beta_mod<-c(beta_mod,mod1$coefficients[2])
  is.na(beta_mod)<-0
  alpha_up<-c(alpha_up,mod_up$coefficients[1])
  is.na(alpha_up)<-0
  beta_up<-c(beta_up,mod_up$coefficients[2])
  is.na(beta_up)<-0
  alpha_down<-c(alpha_down,mod_down$coefficients[1])
  is.na(alpha_down)<-0
  beta_down<-c(beta_down,mod_down$coefficients[2])
  is.na(beta_down)<-0
  LET<-c(LET,as.numeric(file))
  ############################
  
  # Dose corrispondente (selezionare i punti per fit solo lineare)
  if(ion == "H"){
    if(FALSE){
      d10_exp<-c(d10_exp,log(survperc)/mod$coefficients[1])
      d10_up<-c(d10_up,log(survperc)/mod_up$coefficients[1])
      d10_down<-c(d10_down,log(survperc)/mod_down$coefficients[1])
      beta_exp[ctr]<-0
      beta_up[ctr]<-0
      beta_down[ctr]<-0
    }else{
      d10_exp<-c(d10_exp,quadraticRoots(1,mod$coefficients[1]/mod$coefficients[2],-log(survperc)/mod$coefficients[2]))
      d10_up<-c(d10_up,quadraticRoots(1,mod_up$coefficients[1]/mod_up$coefficients[2],-log(survperc)/mod_up$coefficients[2]))
      d10_down<-c(d10_down,quadraticRoots(1,mod_down$coefficients[1]/mod_down$coefficients[2],-log(survperc)/mod_down$coefficients[2]))
    }
    if(FALSE){
      d10_mod<-c(d10_mod,log(survperc)/mod1$coefficients[1])
      beta_mod[ctr]<-0
    }else{
      d10_mod<-c(d10_mod,quadraticRoots(1,mod1$coefficients[1]/mod1$coefficients[2],-log(survperc)/mod1$coefficients[2]))
    }
  }
  
  if(ion == "He"){
  if(file == LET_cycle[10] || file == LET_cycle[11] || file == LET_cycle[12]){
    d10_exp<-c(d10_exp,log(survperc)/mod$coefficients[1])
    d10_up<-c(d10_up,log(survperc)/mod_up$coefficients[1])
    d10_down<-c(d10_down,log(survperc)/mod_down$coefficients[1])
    beta_exp[ctr]<-0
    beta_up[ctr]<-0
    beta_down[ctr]<-0
  }else{
    d10_exp<-c(d10_exp,quadraticRoots(1,mod$coefficients[1]/mod$coefficients[2],-log(survperc)/mod$coefficients[2]))
    d10_up<-c(d10_up,quadraticRoots(1,mod_up$coefficients[1]/mod_up$coefficients[2],-log(survperc)/mod_up$coefficients[2]))
    d10_down<-c(d10_down,quadraticRoots(1,mod_down$coefficients[1]/mod_down$coefficients[2],-log(survperc)/mod_down$coefficients[2]))
  }
  if(file == LET_cycle[4] || file == LET_cycle[5] || file == LET_cycle[6] || file == LET_cycle[7] || file == LET_cycle[8] || file == LET_cycle[9] || file == LET_cycle[10] || file == LET_cycle[11] || file == LET_cycle[12]){
    d10_mod<-c(d10_mod,log(survperc)/mod1$coefficients[1])
    beta_mod[ctr]<-0
  }else{
    d10_mod<-c(d10_mod,quadraticRoots(1,mod1$coefficients[1]/mod1$coefficients[2],-log(survperc)/mod1$coefficients[2]))
  }
  }
  
  if(ion == "C"){
    if(file == LET_cycle[4] || file == LET_cycle[5] || file == LET_cycle[6] || file == LET_cycle[7] || file == LET_cycle[8] || file == LET_cycle[9] || file == LET_cycle[10] || file == LET_cycle[11] || file == LET_cycle[12]){
      d10_exp<-c(d10_exp,log(survperc)/mod$coefficients[1])
      d10_up<-c(d10_up,log(survperc)/mod_up$coefficients[1])
      d10_down<-c(d10_down,log(survperc)/mod_down$coefficients[1])
      beta_exp[ctr]<-0
      beta_up[ctr]<-0
      beta_down[ctr]<-0
    }else{
      d10_exp<-c(d10_exp,quadraticRoots(1,mod$coefficients[1]/mod$coefficients[2],-log(survperc)/mod$coefficients[2]))
      d10_up<-c(d10_up,quadraticRoots(1,mod_up$coefficients[1]/mod_up$coefficients[2],-log(survperc)/mod_up$coefficients[2]))
      d10_down<-c(d10_down,quadraticRoots(1,mod_down$coefficients[1]/mod_down$coefficients[2],-log(survperc)/mod_down$coefficients[2]))
    }
    if(file != LET_cycle[1]){
      d10_mod<-c(d10_mod,log(survperc)/mod1$coefficients[1])
      beta_mod[ctr]<-0
    }else{
      d10_mod<-c(d10_mod,quadraticRoots(1,mod1$coefficients[1]/mod1$coefficients[2],-log(survperc)/mod1$coefficients[2]))
    }
  }
  
  # Mean Inactivation Dose
  #MID_exp_value<-0
  #MID_mod_value<-0
  #for(integral_count in seq(0, d10_exp[ctr]-integral_spacing, by=integral_spacing)){
  #  #print(integral_count)
  #  MID_exp_value<-MID_exp_value+(integral_spacing/2)*(exp(+mod$coefficients[1]*integral_count+mod$coefficients[2]*integral_count*integral_count)+exp(+mod$coefficients[1]*(integral_count+integral_spacing)+mod$coefficients[2]*(integral_count+integral_spacing)*(integral_count+integral_spacing)))
  #  MID_mod_value<-MID_mod_value+(integral_spacing/2)*(exp(+mod1$coefficients[1]*integral_count+mod1$coefficients[2]*integral_count*integral_count)+exp(+mod1$coefficients[1]*(integral_count+integral_spacing)+mod1$coefficients[2]*(integral_count+integral_spacing)*(integral_count+integral_spacing)))
  #}
  #MID_exp<-c(MID_exp,MID_exp_value)
  #MID_mod<-c(MID_mod,MID_mod_value)
  
  ctr<-ctr+1
}

#MID <- as.data.frame(c(MID_exp,MID_mod))
alpha <- as.data.frame(c(alpha_exp,alpha_mod))
beta <- as.data.frame(c(beta_exp,beta_mod))
alpha_Up <- as.data.frame(c(alpha_up,alpha_mod))
alpha_Down <- as.data.frame(c(alpha_down,alpha_mod))
beta_Up <- as.data.frame(c(beta_up,beta_mod))
beta_Down <- as.data.frame(c(beta_down,beta_mod))
#ratio <- as.data.frame(c(alpha_exp/beta_exp,alpha_mod/beta_mod))

Type<-c()
for(counter in c(1:length(LET))){
  Type<-c(Type,"MD Anderson data")
}
for(counter in c(1:length(LET))){
  Type<-c(Type,"GSM2")
}
LET<-rep(LET,2)
if(ion=="He"){
  yD<-rep(yD_He,2)
}else if(ion=="H"){
  yD<-rep(yD_H,2)
}else if(ion=="C"){
  yD<-rep(yD_C,2)
}
D10_exp<-as.data.frame(d10_exp)

D10_up<-as.data.frame(d10_up)
D10_down<-as.data.frame(as.numeric(d10_down))

D10_mod<-as.data.frame(d10_mod)
# da modificare RBE 10-50-90
RBE10_exp<-D10_Xray/D10_exp

RBE10_up <-D10_Xray/D10_down #opposite (high SF means low RBE)
RBE10_down <-D10_Xray/D10_up #opposite (high SF means low RBE)

RBE10_mod<-D10_Xray/D10_mod

RBE10_up <- as.data.frame(RBE10_up)
colnames(RBE10_up)[1]<-"Up"
errorbar_up <- rbind(RBE10_up,as.data.frame(RBE10_mod %>% rename("Up" = d10_mod)))
for(ii in c(1:12)){
  if(errorbar_up[ii,1] <= RBE10_exp[ii,1]){
    errorbar_up[ii,1]<-RBE10_exp[ii,1] + 0.1*RBE10_exp[ii,1]
  }
}

RBE10_down <- as.data.frame(RBE10_down)
colnames(RBE10_down)[1]<-"Down"
errorbar_down <- rbind(RBE10_down,as.data.frame(RBE10_mod %>% rename("Down" = d10_mod)))
for(ii in c(1:12)){
  if(errorbar_down[ii,1] >= RBE10_exp[ii,1]){
    errorbar_down[ii,1]<-RBE10_exp[ii,1] - 0.1*RBE10_exp[ii,1]
  }
}

plotD10<-rbind(as.data.frame(D10_exp %>% rename("D10" = d10_exp)),as.data.frame(D10_mod %>% rename("D10" = d10_mod)))
plotD10<-cbind(plotD10,LET,Type)
plotRBE10<-rbind(as.data.frame(RBE10_exp %>% rename("RBE10" = d10_exp)),as.data.frame(RBE10_mod %>% rename("RBE10" = d10_mod)))
plotRBE10<-cbind(plotRBE10,LET,yD,Type,alpha,beta,alpha_Up,alpha_Down,beta_Up,beta_Down,errorbar_up,errorbar_down)
#colnames(plotRBE10)[5]<-"MID"
colnames(plotRBE10)[5]<-"alpha"
colnames(plotRBE10)[6]<-"beta"
colnames(plotRBE10)[7]<-"alpha_up"
colnames(plotRBE10)[8]<-"alpha_down"
colnames(plotRBE10)[9]<-"beta_up"
colnames(plotRBE10)[10]<-"beta_down"

if(ion=="He"){
  Ion<-rep("He",24)
  RBE_He<-cbind(plotRBE10,Ion)
}else if(ion=="H"){
  Ion<-rep("H",24)
  RBE_H<-cbind(plotRBE10,Ion)
}else if(ion=="C"){
  Ion<-rep("C",24)
  RBE_C<-cbind(plotRBE10,Ion)
}

# FILTERING DATA
plot <- plotRBE10 %>% filter(plotRBE10$LET != LET_cycle[1])
#plot <- plot %>% filter(plot$LET != LET_cycle[2])
#plot <- plot %>% filter(plot$LET != LET_cycle[3])
#plot <- plot %>% filter(plot$LET != LET_cycle[4])
#plot <- plot %>% filter(plot$LET != LET_cycle[6])
#plot <- plot %>% filter(plot$LET != LET_cycle[7])
#plot <- plot %>% filter(plot$LET != LET_cycle[9])
#plot <- plot %>% filter(plot$LET != LET_cycle[10])
#plot <- plot %>% filter(plot$LET != LET_cycle[11])
#plot <- plot %>% filter(plot$LET != LET_cycle[12])

ggplot(plotRBE10,aes(yD,RBE10,color=Type,shape=Type))+
           geom_line(linewidth=0)+geom_point(size=2)+
           geom_errorbar(aes(ymax = Up, ymin = Down),width=0.1,linewidth=0.5,linetype = "solid")+
           scale_y_log10()+scale_color_manual(values=c25)+
           xlab("yD [keV/um]")+ylab("RBE10")+annotation_logticks(sides = "l")+
           labs(title = "RBE10 Helium (H460 cell line)",
                colour = "Data set")+
           guides(shape="none")+
           theme_bw()+theme(plot.title = element_text(size=10, color="black"),
                            axis.title.x = element_text(size=10, color="black"),
                            axis.title.y = element_text(size=10, color="black"),
                            axis.text.x = element_text(size=10, color="black"),
                            axis.text.y = element_text(size=10, color="black"),
                            legend.title = element_text(size=10, color="black"),
                            legend.text = element_text(size=10, color="black"),
                            legend.position = "right")
}

setwd("/home/user/Scrivania/Dottorato/RISULTATI_FIT/Global/")
save(plotRBE10, file = "FINAL_He_H460.RData")

DATA_full<-rbind(RBE_He,RBE_H,RBE_C)

# GGPLOT RBE
ggplot(DATA_full,aes(LET,MID,color=Ion,shape=Type))+
  geom_line(linewidth=0)+geom_point(size=2)+scale_color_manual(values=c25)+
  xlab("LET [keV/um]")+ylab("MID [Gy]")+
  scale_x_log10()+annotation_logticks(sides = "b")+
  scale_shape_manual(name = "Data set",labels = c("GSM2 model", "MD Anderson data"), values=c(15, 1))+
  labs(title = "MID (H460 cell line)",
       colour = "Ion type",
       shape = "Data set") +
  theme_bw()+theme(plot.title = element_text(size=10, color="black"),
                   axis.title.x = element_text(size=10, color="black"),
                   axis.title.y = element_text(size=10, color="black"),
                   axis.text.x = element_text(size=10, color="black"),
                   axis.text.y = element_text(size=10, color="black"),
                   legend.title = element_text(size=10, color="black"),
                   legend.text = element_text(size=10, color="black"),
                   legend.margin = margin(0,0,0,0, unit="cm"),
                   legend.position = "right")

# PLOT RBEalpha

################################################################################

# ERRORE Chi2 + u=0.5*Beta (DeltaBeta) per dare peso maggiore all'errore di curve
# che non predicono il corretto andamento nel termine quadratico (più sensibile)

chi_sq<-c()
sizes<-c()
rd_vec<-c()
Rn_vec<-c()
file_vec<-c()

Delta_beta<-c()
Delta_alpha<-c()

# DUE LINEE SUCCESSIVE PER SELEZIONARE VALORI rd ED Rn
domain_size <- 0.8
domain_size_cell <- 6
beta_weight <- 0.5
alpha_weight <- 0

for (file in coloumn_corrected) {
    file_vec<-c(file_vec,as.numeric(file))
    surv<-na.omit(surv_GSM2[[paste0(as.character(rd),"-",as.character(Rn))]])
    surv_model<-c()
    surv_real<-c()
    surv_model<-surv[(1+10*(as.numeric(file)-1)):(10*as.numeric(file))]
    surv_real<-S_real[(1+10*(as.numeric(file)-1)):(10*as.numeric(file))]
    Delta_beta<-c(Delta_beta,sqrt((beta_exp[as.numeric(file)]-beta_mod[as.numeric(file)])^2))
    Delta_alpha<-c(Delta_alpha,sqrt((alpha_exp[as.numeric(file)]-alpha_mod[as.numeric(file)])^2))
    chi_sq<-c(chi_sq,sum(((surv_model-surv_real)^2)/(surv_real))+beta_weight*sqrt((beta_exp[as.numeric(file)]-beta_mod[as.numeric(file)])^2)+alpha_weight*sqrt((alpha_exp[as.numeric(file)]-alpha_mod[as.numeric(file)])^2))
}

surv_ranking<-data.frame(Ion = rep(ion,length(chi_sq)),
                         Curve = file_vec,
                         Chi2 = chi_sq,
                         Delta_beta = Delta_beta) %>% View()

Data_H<-data.frame(Ion = rep(ion,length(chi_sq)),
                 Curve = file_vec,
                 Chi2 = chi_sq,
                 Delta_beta = Delta_beta)

write.csv(Data_H, "/home/user/Scrivania/Dottorato/RISULTATI_FIT/RBE/Ranking_H.csv", row.names=FALSE)
save(Data_H, file = "/home/user/Scrivania/Dottorato/RISULTATI_FIT/RBE/Ranking_H.RData")

# GLOBAL RANKING

# bind dataframes
load("~/Scrivania/Dottorato/RISULTATI_FIT/RBE/Ranking_He.RData")
load("~/Scrivania/Dottorato/RISULTATI_FIT/RBE/Ranking_H.RData")
load("~/Scrivania/Dottorato/RISULTATI_FIT/RBE/Ranking_C.RData")
#yD
yD<-c(yD_He,yD_H,yD_C)
yD_data<-as.data.frame(yD)
LET<-c(LET_He,LET_H,LET_C)
LET_vec<-c(LET)
global_ranking<-rbind(Data_He,Data_H,Data_C)
global_ranking_energy<-cbind(global_ranking,yD_data,LET_vec)
write.csv(global_ranking, "/home/user/Scrivania/Dottorato/RISULTATI_FIT/RBE/Ranking_tot.csv", row.names=FALSE)
save(global_ranking, file = "/home/user/Scrivania/Dottorato/RISULTATI_FIT/RBE/Ranking_tot.RData")
View(global_ranking)

accepted_ranking <- global_ranking_energy[which(global_ranking_energy$Chi2 <= 1),names(global_ranking_energy) %in% c("Ion","Curve","Chi2","Delta_beta","yD","LET_vec")]
accepted_ranking_yD <- global_ranking_yD[which(global_ranking_yD$Chi2 <= 1),names(global_ranking_yD) %in% c("Ion","Curve","Chi2","Delta_beta","yD","LET_vec")]


# GGPLOT RANKING
ggplot(accepted_ranking,aes(as.numeric(yD),Chi2,color=Ion,shape=Ion))+
  geom_point(size=3)+scale_color_manual(values=c25)+
  xlab("yD [keV/um]")+ylab("Chi2")+
  scale_x_log10(breaks = c(7,10,25,50,100,180), labels=c("7","10","25","50","100","180"))+
  labs(title = "Chi2 estimator",
       colour = "Data set")+
  guides(shape="none")+
  geom_hline(yintercept=0.6,color = "red",linewidth=1.5)+
  scale_y_continuous(breaks = c(0,0.25,0.5,0.6,0.75,1), labels=c("0","0.25","0.5","Threshold","0.75","1"))+
  theme_bw()+theme(plot.title = element_text(size=10, color="black"),
                   axis.title.x = element_text(size=10, color="black"),
                   axis.title.y = element_text(size=10, color="black"),
                   axis.text.x = element_text(size=10, color="black"),
                   axis.text.y = element_text(size=10, color="black"),
                   legend.title = element_text(size=10, color="black"),
                   legend.text = element_text(size=10, color="black"),
                   legend.margin = margin(0,0,0,0, unit="cm"),
                   legend.position = c(0.1,0.9))
