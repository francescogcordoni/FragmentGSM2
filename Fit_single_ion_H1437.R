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
library(readxl)
library(tidymodels)

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

#path_or<-"/home/user/Scrivania/TESI/RISULTATI_FIT/PROVA_3"
#setwd(path_or)
#data<-read.csv("Fit_parameters.csv")

path_or<-"/home/user/Scrivania/Dottorato/RISULTATI_FIT"
setwd(path_or)
source("/home/user/Scrivania/TESI/Codici_R/utilities_doserate.R")

# CELL LINE
cell<-"H460"
cell<-"H1437"

# IONE DA SELEZIONARE
ion<-"He"
ion<-"H"
ion<-"C"

# DEFINISCO ENERGY_LIST PER OGNI IONE
if(ion=="He"){
  energy_list<-c("3.4","8.5","12.7","22.4","36.2","44.9","54.7","63","71.1","77.6","83.6","88.7")  
}else if(ion=="H"){
  energy_list<-c("1","2.6","4.7","7.3","8.7","11.1","13.7","15.4","16.9","18.3","20.2","21.4")
}else if(ion=="C"){
  energy_list<-c("20.2","39.8","63.1","70.6","84.3","100.8","126.3","157.6","196.4","242.9","285.8","308.4")
}

# DEFINISCO I RAGGI TIPICI PER DOMINI E CELLULA PER GLI SPETTRI MICRO E LE DIMENSIONI CHE ENTRANO NEL CICLO
rd<-0.8; Rn<-8
domain_size<-c(0.35,0.4,0.45,0.5,0.6,0.8,1)
domain_size_cell<-c(5,6,7,8)

# DEFINISCO IL PATH PER GLI SPETTRI MICRO CON R=8um
radius<-c("SPECTRA_8um")

# DEFINISCO IL VETTORE PER CICLARE SULLE COLONNE
coloumn<-as.character(1:12)

# CARICO I DATI BIOLOGICI
#load("/home/user/Scrivania/TESI/Codici_R/Survival_CORRECTED.RData")
#load("/home/user/Scrivania/TESI/Codici_R/Survival_removed_points.RData")
#load("/home/user/Scrivania/TESI/Codici_R/Survival_all.RData")
#load("/home/user/Scrivania/Dottorato/Codici_R/Survival_uniform_dose_sampling.RData")
load("/home/user/Scrivania/Dottorato/Codici_R/Survival_uniform_dose_sampling_H1437.RData")
pide_data<-survival

# FLAG <- 0 (H460) FLAG <- 1 (H1437)
#flag<-0;
flag<-1;
# CARICO LE DISTRIBUZIONI DI DANNO E LA MULTI-EVENT
if(flag == 0){
  if(ion=="He"){
    load("/home/user/Scrivania/Dottorato/Kappa_parametrization/KAPPA_PARTRAC/Helium_x9/p0X_MD.RData")
    load("/home/user/Scrivania/Dottorato/Kappa_parametrization/KAPPA_PARTRAC/Helium_x9/p0Y_MD.RData")
    load("/home/user/Scrivania/Dottorato/Kappa_parametrization/KAPPA_PARTRAC/Helium_x9/fn_dose_all_MD.RData") 
  }else if(ion=="H"){
    load("/home/user/Scrivania/Dottorato/Kappa_parametrization/KAPPA_PARTRAC/Protons_x9/p0X_MD.RData")
    load("/home/user/Scrivania/Dottorato/Kappa_parametrization/KAPPA_PARTRAC/Protons_x9/p0Y_MD.RData")
    load("/home/user/Scrivania/Dottorato/Kappa_parametrization/KAPPA_PARTRAC/Protons_x9/fn_dose_all_MD.RData") 
  }else if(ion=="C"){
    load("/home/user/Scrivania/Dottorato/Kappa_parametrization/KAPPA_PARTRAC/Carbon_x9/p0X_MD.RData")
    load("/home/user/Scrivania/Dottorato/Kappa_parametrization/KAPPA_PARTRAC/Carbon_x9/p0Y_MD.RData")
    load("/home/user/Scrivania/Dottorato/Kappa_parametrization/KAPPA_PARTRAC/Carbon_x9/fn_dose_all_MD.RData")
  }
}else if(flag == 1){
  if(ion=="He"){
    load("/home/user/Scrivania/Dottorato/Damage_distributions/Helium_x9/p0X_MD.RData")
    load("/home/user/Scrivania/Dottorato/Damage_distributions/Helium_x9/p0Y_MD.RData")
    load("/home/user/Scrivania/Dottorato/Damage_distributions/Helium_x9/fn_dose_all_MD.RData") 
  }else if(ion=="H"){
    load("/home/user/Scrivania/Dottorato/Damage_distributions/Protons_x9/p0X_MD.RData")
    load("/home/user/Scrivania/Dottorato/Damage_distributions/Protons_x9/p0Y_MD.RData")
    load("/home/user/Scrivania/Dottorato/Damage_distributions/Protons_x9/fn_dose_all_MD.RData") 
  }else if(ion=="C"){
    load("/home/user/Scrivania/Dottorato/Damage_distributions/Carbon_x9/p0X_MD.RData")
    load("/home/user/Scrivania/Dottorato/Damage_distributions/Carbon_x9/p0Y_MD.RData")
    load("/home/user/Scrivania/Dottorato/Damage_distributions/Carbon_x9/fn_dose_all_MD.RData")
  }
}


# err^2 removed in optimization (err is a constant)
Fit_GSM2 <- function(x) {
  
  a<-x[1]
  b<-x[2]
  r<-x[3]
  
  N<-floor((Rn^3)/(rd^3))
  
  surv<-c()
  ctr_e<-1
  
  if(flag == 0){
    if(ion=="He"){
      load("/home/user/Scrivania/Dottorato/Kappa_parametrization/KAPPA_PARTRAC/Helium_x9/p0X_MD.RData")
      load("/home/user/Scrivania/Dottorato/Kappa_parametrization/KAPPA_PARTRAC/Helium_x9/p0Y_MD.RData")
      load("/home/user/Scrivania/Dottorato/Kappa_parametrization/KAPPA_PARTRAC/Helium_x9/fn_dose_all_MD.RData") 
    }else if(ion=="H"){
      load("/home/user/Scrivania/Dottorato/Kappa_parametrization/KAPPA_PARTRAC/Protons_x9/p0X_MD.RData")
      load("/home/user/Scrivania/Dottorato/Kappa_parametrization/KAPPA_PARTRAC/Protons_x9/p0Y_MD.RData")
      load("/home/user/Scrivania/Dottorato/Kappa_parametrization/KAPPA_PARTRAC/Protons_x9/fn_dose_all_MD.RData") 
    }else if(ion=="C"){
      load("/home/user/Scrivania/Dottorato/Kappa_parametrization/KAPPA_PARTRAC/Carbon_x9/p0X_MD.RData")
      load("/home/user/Scrivania/Dottorato/Kappa_parametrization/KAPPA_PARTRAC/Carbon_x9/p0Y_MD.RData")
      load("/home/user/Scrivania/Dottorato/Kappa_parametrization/KAPPA_PARTRAC/Carbon_x9/fn_dose_all_MD.RData")
    }
  }else if(flag == 1){
    if(ion=="He"){
      load("/home/user/Scrivania/Dottorato/Damage_distributions/Helium_x9/p0X_MD.RData")
      load("/home/user/Scrivania/Dottorato/Damage_distributions/Helium_x9/p0Y_MD.RData")
      load("/home/user/Scrivania/Dottorato/Damage_distributions/Helium_x9/fn_dose_all_MD.RData") 
    }else if(ion=="H"){
      load("/home/user/Scrivania/Dottorato/Damage_distributions/Protons_x9/p0X_MD.RData")
      load("/home/user/Scrivania/Dottorato/Damage_distributions/Protons_x9/p0Y_MD.RData")
      load("/home/user/Scrivania/Dottorato/Damage_distributions/Protons_x9/fn_dose_all_MD.RData") 
    }else if(ion=="C"){
      load("/home/user/Scrivania/Dottorato/Damage_distributions/Carbon_x9/p0X_MD.RData")
      load("/home/user/Scrivania/Dottorato/Damage_distributions/Carbon_x9/p0Y_MD.RData")
      load("/home/user/Scrivania/Dottorato/Damage_distributions/Carbon_x9/fn_dose_all_MD.RData")
    }
  }
  
  
  surv<-c()
  ctr_e<-1
  
  for (file in coloumn_corrected) {
    #la dose diventa la dose contenuta nel file survival rispetto alla linea cellulare, allo ione e all'energia
    if(ion=="He"){
      dose<-pide_data[[cell]]$He[[energy_list[ctr_e]]]$Dose  
    }else if(ion=="H"){
      dose<-pide_data[[cell]]$H[[energy_list[ctr_e]]]$Dose
    }else if(ion=="C"){
      dose<-pide_data[[cell]]$C[[energy_list[ctr_e]]]$Dose
    }
    
    ctr<-1
    surv_GSM<-c()
    
    for (zn in dose_zn) {
      
      dose_name<-as.character(zn)
      
      surv_GSM[ctr]<-p0X[[paste0(as.character(rd),"-",as.character(Rn))]][[file]][[dose_name]][1]*p0Y[[paste0(as.character(rd),"-",as.character(Rn))]][[file]][[dose_name]][1]
      
      l_max<-length(p0X[[paste0(as.character(rd),"-",as.character(Rn))]][[file]][[dose_name]])
      
      for (k in 2:l_max) {
        
        surv_GSM[ctr]<-surv_GSM[ctr]+
          ifelse(is.na(p0X[[paste0(as.character(rd),"-",as.character(Rn))]][[file]][[dose_name]][k]*p0Y[[paste0(as.character(rd),"-",as.character(Rn))]][[file]][[dose_name]][1]*Csur(0,0,k-1,a,b,r)),0,
                 p0X[[paste0(as.character(rd),"-",as.character(Rn))]][[file]][[dose_name]][k]*p0Y[[paste0(as.character(rd),"-",as.character(Rn))]][[file]][[dose_name]][1]*Csur(0,0,k-1,a,b,r))
      }
      ctr<-ctr+1
    }
    
    surv_GSM[is.na(surv_GSM)]<-0
    
    surv_GSM<-surv_GSM^N
    
    surv_cell<-c()
    N<-floor((Rn^3)/(rd^3))
    ctr<-1
    for (d in unique(dose)) {
      surv_cell[ctr]<-sum(fn_dose_all[[as.character(Rn)]][[file]][[ctr]]$zfz*(surv_GSM))/sum(fn_dose_all[[as.character(Rn)]][[file]][[ctr]]$zfz)
      ctr<-ctr+1
    }
    
    s_p_s<-log(surv_cell)
    dose_s<-dose
    dose2<-dose_s*dose_s
    
    lq<-lm(formula = s_p_s ~ dose_s + dose2)
    #surv<-c(surv,exp(coefficients(lq)[2]*dose_s + coefficients(lq)[3]*dose2))
    
    # CODICE PER EVITARE BETA POSITIVO
    if(coefficients(lq)[3] > 0){
      s_p_s<-log(surv_cell)
      dose_s<-dose
      dose2<-dose_s*dose_s
      
      lq<-lm(formula = s_p_s ~ dose_s)
      surv<-c(surv,exp(coefficients(lq)[2]*dose_s))
      
    }else{
      surv<-c(surv,exp(coefficients(lq)[2]*dose_s + coefficients(lq)[3]*dose2))
    }
    
    # surv<-c(surv,surv_cell)
    
    ctr_e<-ctr_e+1
  }
  
  surv<-na.omit(surv)
  
  #sum(((surv-S_real)^2)/((S_real)*err^2))
  sum(((surv-S_real)^2)/(S_real))
  
  # return(res)
}

# SOTTRAGGO ALCUNE DELLE 12 COLONNE (DA MODIFICARE SE NECESSARIO)
rem_col<-c(1,2,3,4,5,6,7,8,11,12)
no_removal<-27 # valore a caso
no_removal<-0 # set to 0 for plotting all SF curves
dose_zn<-fn_dose_all[[1]][[1]][[1]]$z
if(no_removal==0){
  coloumn_corrected<-coloumn
  energy_list<-energy_list
}else{
  coloumn_corrected<-coloumn[-rem_col]
  energy_list<-energy_list[-rem_col]
}

# CARICO LA VERA SURVIVAL SF CICLANDO SULLE ENERGIE (da mettere dentro il ciclo per togliere poi alcune colonne)
S_real<-c(); err<-c()
ctr_e<-1
for(file in coloumn_corrected){
  if(ion=="He"){
    S_real<-c(S_real,pide_data[[cell]]$He[[energy_list[ctr_e]]]$LQ) 
    err<-c(err,pide_data[[cell]]$He[[energy_list[ctr_e]]]$ErrSF) 
  }else if(ion=="H"){
    S_real<-c(S_real,pide_data[[cell]]$H[[energy_list[ctr_e]]]$LQ)
    err<-c(err,pide_data[[cell]]$H[[energy_list[ctr_e]]]$ErrSF) 
  }else if(ion=="C"){
    S_real<-c(S_real,pide_data[[cell]]$C[[energy_list[ctr_e]]]$LQ)
    err<-c(err,pide_data[[cell]]$C[[energy_list[ctr_e]]]$ErrSF)
  }
  ctr_e<-ctr_e+1
}

# DEFINISCO I PARAMETRI DA FITTARE
a_param<-c();b_param<-c();r_param<-c()
surv_GSM2<-list()
flag_linear_fit <- 0
domain_size<-c(0.35,0.4,0.45,0.5,0.6,0.8,1)
domain_size_cell<-c(5,6,7,8)
# Customized sizes
#domain_size<-c(0.35)
#domain_size_cell<-c(8)
#play with starting point of a,b and r (use expectations)
# INIZIO CICLO
for (rd in domain_size) {
  for (Rn in domain_size_cell) {
    
    print(paste0("rd=", rd," Rn=", Rn))
    
    result_<-optim(c(0.01,0.01,3),lower=c(0.0001,0.00001,1),upper=c(1,1,16), 
                   Fit_GSM2, method="L-BFGS-B")
    
    (a<-result_$par[1])
    (b<-result_$par[2])
    (r<-result_$par[3])
    
    # SCRIVO PARAMETRI OTTENUTI DA FIT PRECEDENTI PER PREDIRE (H460)
    #rd=0.8;
    #Rn=6;
    #a=0.0008997624
    #b=0.06421062
    #r=2.714594
    
    # SCRIVO PARAMETRI OTTENUTI DA FIT PRECEDENTI PER PREDIRE (H1437)
    #rd=0.45;
    #Rn=8;
    #a=0.003026402
    #b=0.147680249
    #r=2.951922
    
    # From global fit
    #rd=0.4;
    #Rn=5;
    #a=0.018889
    #b=0.0030425
    #r=3.0051
    
    # stampa parametri
    print(paste0("a=", a))
    print(paste0("b=", b))
    print(paste0("r=", r))
    print(paste0("chi2=", Fit_GSM2(c(a,b,r))))
    
    a_param<-c(a_param,a)
    b_param<-c(b_param,b)
    r_param<-c(r_param,r)
    
    rc<-Rn
    N<-floor((Rn^3)/(rd^3))
    
    if(flag == 0){
      if(ion=="He"){
        load("/home/user/Scrivania/Dottorato/Kappa_parametrization/KAPPA_PARTRAC/Helium_x9/p0X_MD.RData")
        load("/home/user/Scrivania/Dottorato/Kappa_parametrization/KAPPA_PARTRAC/Helium_x9/p0Y_MD.RData")
        load("/home/user/Scrivania/Dottorato/Kappa_parametrization/KAPPA_PARTRAC/Helium_x9/fn_dose_all_MD.RData") 
      }else if(ion=="H"){
        load("/home/user/Scrivania/Dottorato/Kappa_parametrization/KAPPA_PARTRAC/Protons_x9/p0X_MD.RData")
        load("/home/user/Scrivania/Dottorato/Kappa_parametrization/KAPPA_PARTRAC/Protons_x9/p0Y_MD.RData")
        load("/home/user/Scrivania/Dottorato/Kappa_parametrization/KAPPA_PARTRAC/Protons_x9/fn_dose_all_MD.RData") 
      }else if(ion=="C"){
        load("/home/user/Scrivania/Dottorato/Kappa_parametrization/KAPPA_PARTRAC/Carbon_x9/p0X_MD.RData")
        load("/home/user/Scrivania/Dottorato/Kappa_parametrization/KAPPA_PARTRAC/Carbon_x9/p0Y_MD.RData")
        load("/home/user/Scrivania/Dottorato/Kappa_parametrization/KAPPA_PARTRAC/Carbon_x9/fn_dose_all_MD.RData")
      }
    }else if(flag == 1){
      if(ion=="He"){
        load("/home/user/Scrivania/Dottorato/Damage_distributions/Helium_x9/p0X_MD.RData")
        load("/home/user/Scrivania/Dottorato/Damage_distributions/Helium_x9/p0Y_MD.RData")
        load("/home/user/Scrivania/Dottorato/Damage_distributions/Helium_x9/fn_dose_all_MD.RData") 
      }else if(ion=="H"){
        load("/home/user/Scrivania/Dottorato/Damage_distributions/Protons_x9/p0X_MD.RData")
        load("/home/user/Scrivania/Dottorato/Damage_distributions/Protons_x9/p0Y_MD.RData")
        load("/home/user/Scrivania/Dottorato/Damage_distributions/Protons_x9/fn_dose_all_MD.RData") 
      }else if(ion=="C"){
        load("/home/user/Scrivania/Dottorato/Damage_distributions/Carbon_x9/p0X_MD.RData")
        load("/home/user/Scrivania/Dottorato/Damage_distributions/Carbon_x9/p0Y_MD.RData")
        load("/home/user/Scrivania/Dottorato/Damage_distributions/Carbon_x9/fn_dose_all_MD.RData")
      }
    }
    
    surv<-c()
    ctr_e<-1
    for (file in coloumn_corrected) {
      print(paste0("file=", file))
      
      #la dose diventa la dose contenuta nel file survival rispetto alla linea cellulare, allo ione e all'energia
      if(ion=="He"){
        dose<-pide_data[[cell]]$He[[energy_list[ctr_e]]]$Dose  
      }else if(ion=="H"){
        dose<-pide_data[[cell]]$H[[energy_list[ctr_e]]]$Dose
      }else if(ion=="C"){
        dose<-pide_data[[cell]]$C[[energy_list[ctr_e]]]$Dose
      }
      
      ctr<-1
      surv_GSM<-c()
      
      for (zn in dose_zn) {
        
        dose_name<-as.character(zn)
        
        surv_GSM[ctr]<-p0X[[paste0(as.character(rd),"-",as.character(Rn))]][[file]][[dose_name]][1]*p0Y[[paste0(as.character(rd),"-",as.character(Rn))]][[file]][[dose_name]][1]
        
        #damage_x[ctr]<-p0X[[paste0(as.character(rd),"-",as.character(Rn))]][[file]][[dose_name]][1]
        #c_x[ctr]<-Csur(0,0,0,a,b,r)
        
        l_max<-length(p0X[[paste0(as.character(rd),"-",as.character(Rn))]][[file]][[dose_name]])
        
        for (k in 2:l_max) {
          surv_GSM[ctr]<-surv_GSM[ctr]+
            ifelse(is.na(p0X[[paste0(as.character(rd),"-",as.character(Rn))]][[file]][[dose_name]][k]*p0Y[[paste0(as.character(rd),"-",as.character(Rn))]][[file]][[dose_name]][1]*Csur(0,0,k-1,a,b,r)),0,
                   p0X[[paste0(as.character(rd),"-",as.character(Rn))]][[file]][[dose_name]][k]*p0Y[[paste0(as.character(rd),"-",as.character(Rn))]][[file]][[dose_name]][1]*Csur(0,0,k-1,a,b,r))
        }
        ctr<-ctr+1
      }
      surv_GSM[is.na(surv_GSM)]<-0
      
      surv_GSM<-surv_GSM^N
      
      surv_cell<-c()
      N<-floor((Rn^3)/(rd^3))
      ctr<-1
      for (d in unique(dose)) {
        surv_cell[ctr]<-sum(fn_dose_all[[as.character(Rn)]][[file]][[ctr]]$zfz*(surv_GSM))/sum(fn_dose_all[[as.character(Rn)]][[file]][[ctr]]$zfz)
        ctr<-ctr+1
      }
      
      s_p_s<-log(surv_cell)
      dose_s<-dose
      dose2<-dose_s*dose_s
      
      lq<-lm(formula = s_p_s ~ dose_s + dose2)
      #surv<-c(surv,exp(coefficients(lq)[2]*dose_s + coefficients(lq)[3]*dose2))
      
      # CODICE PER EVITARE BETA POSITIVO
      if(coefficients(lq)[3] > 0 || flag_linear_fit == 1){
        s_p_s<-log(surv_cell)
        dose_s<-dose
        dose2<-dose_s*dose_s
        
        lq<-lm(formula = s_p_s ~ dose_s)
        surv<-c(surv,exp(coefficients(lq)[2]*dose_s))
        
      }else{
        surv<-c(surv,exp(coefficients(lq)[2]*dose_s + coefficients(lq)[3]*dose2))
      }
      
      # surv<-c(surv,surv_cell)
      
      ctr_e<-ctr_e+1
    }
    
    surv_GSM2[[paste0(as.character(rd),"-",as.character(Rn))]]<-surv  
  }
}

surv<-na.omit(surv)

energy_list_all<-c()
dose_list_all<-c()
ctr_e<-1
for (file in coloumn_corrected) {
  print(file)
  
  #la dose diventa la dose contenuta nel file survival rispetto alla linea cellulare, allo ione e all'energia
  if(ion=="He"){
    dose_list_all<-c(dose_list_all,pide_data$H1437$He[[energy_list[ctr_e]]]$Dose)
    energy_list_all<-c(energy_list_all,rep(energy_list[ctr_e],length(pide_data$H1437$He[[energy_list[ctr_e]]]$Dose)))
  }else if(ion=="H"){
    dose_list_all<-c(dose_list_all,pide_data$H1437$H[[energy_list[ctr_e]]]$Dose)
    energy_list_all<-c(energy_list_all,rep(energy_list[ctr_e],length(pide_data$H1437$H[[energy_list[ctr_e]]]$Dose)))
  }else if(ion=="C"){
    dose_list_all<-c(dose_list_all,pide_data$H1437$C[[energy_list[ctr_e]]]$Dose)
    energy_list_all<-c(energy_list_all,rep(energy_list[ctr_e],length(pide_data$H1437$C[[energy_list[ctr_e]]]$Dose)))
  }
  ctr_e<-ctr_e+1
}

chi_sq<-c()
sizes<-c()
rd_vec<-c()
Rn_vec<-c()

# DUE LINEE SUCCESSIVE PER SELEZIONARE VALORI rd ED Rn
domain_size <- 0.4
domain_size_cell <- 5
domain_size<-c(0.35,0.4,0.45,0.5,0.6,0.8,1)
domain_size_cell<-c(5,6,7,8)
# Customized sizes
#domain_size<-c(0.45,0.5)
#domain_size_cell<-c(7,8)

for (rd in domain_size) {
  for (Rn in domain_size_cell) {
    sizes<-c(sizes,paste0(as.character(rd),"-",as.character(Rn)))
    rd_vec<-c(rd_vec,rd)
    Rn_vec<-c(Rn_vec,Rn)
    surv<-na.omit(surv_GSM2[[paste0(as.character(rd),"-",as.character(Rn))]])
    chi_sq<-c(chi_sq,sum(((surv-S_real)^2)/S_real*err))
    
    plot_df<-data.frame(Dose = rep(dose_list_all,(2*length(surv))/length(dose_list_all)),
                        Survival = c(surv,S_real),
                        Energy = rep(energy_list_all,2),
                        Type = c(rep("GSM2",length(surv)),rep("Exp",length(surv))))
  }
}

data.frame(Chi=chi_sq,Size=sizes,a=a_param,b=b_param,r=r_param,ratio=a_param/b_param) %>% View()
Fit_parameters <- data.frame(Chi=chi_sq,rd=rd_vec,Rn=Rn_vec,a=a_param,b=b_param,r=r_param,ratio=a_param/b_param)
# CAMBIARE IL PATH PER ASSEGNARLO ALLA PROVA CORRISPONDENTE
write.csv(Fit_parameters, "/home/user/Scrivania/Dottorato/RISULTATI_FIT/Par_C_1_2_4_6.csv", row.names=FALSE)

############ ERROR BAR ############
plot_s<-plot_df
#barra errore esperimento
error_up<-c()
error_down<-c()
error_up_<-c()
error_down_<-c()
for (file in as.numeric(coloumn_corrected)) {
  error_up<-c(error_up,survival$H1437[[ion]][[file]]$LQ + 0.5*survival$H1437[[ion]][[file]]$ErrSF)
  error_down<-c(error_down,survival$H1437[[ion]][[file]]$LQ - 0.5*survival$H1437[[ion]][[file]]$ErrSF)
}
#barra errore modello coincide con la survival DEL MODELLO (non ho errore da aggiungere)
error_up_<-c(error_up_,surv)
error_down_<-c(error_down_,surv)

colnames(plot_s)[3]<-"Energy [MeV/u]"

plot_s<-plot_s %>% mutate(Up = c(error_up_,error_up), Down = c(error_down_,error_down))

plot_s$Up<-ifelse(plot_s$Up > 1,1-0.001,plot_s$Up)

#setwd("/home/user/Scrivania/TESI/RISULTATI_FIT/Data_fit_He")
#save(plot_s,file="Data_fit_He_1_6.RData")

# GRAFICO CON BARRE D'ERRORE (AGGIUNTA DI SIZE IN GEOM_POINT E SCALE_SHAPE_MANUAL)
ggplotly(ggplot(plot_s,aes(Dose,Survival,color=`Energy [MeV/u]`,linetype=Type,shape=Type))+
           geom_line(linewidth=0.5)+geom_point(size=1)+
           scale_y_log10()+scale_color_manual(values=c25)+scale_linetype_manual(values=c("dashed","solid"))+
           scale_shape_manual(values=c(16, 18))+
           geom_errorbar(aes(ymax = Up, ymin = Down),width=0.01,linewidth=1,linetype = "solid")+
           xlab("Dose [Gy]")+ylab("Surviving fraction")+annotation_logticks(sides = "l")+
           theme_bw()+theme(plot.title = element_text(size=10, color="black"),
                            axis.title.x = element_text(size=10, color="black"),
                            axis.title.y = element_text(size=10, color="black"),
                            axis.text.x = element_text(size=10, color="black"),
                            axis.text.y = element_text(size=10, color="black"),
                            legend.title = element_text(size=10, color="black"),
                            legend.text = element_text(size=10, color="black"),
                            legend.position = c(0.9,.8))
)


## GGPLOT PER PLOT FINALI
ggplot(plot_s,aes(Dose,Survival,color=`Energy [MeV/u]`,linetype=Type,shape=Type))+
  geom_line(linewidth=0.5)+geom_point(size=1.5)+
  scale_y_log10()+scale_color_manual(values=c25)+scale_linetype_manual(values=c("blank","solid"))+
  scale_shape_manual(labels = c("Experiments", "GSM2 model"), values=c(16, 18))+
  geom_errorbar(aes(ymax = Up, ymin = Down),width=0.01,linewidth=0.5,linetype = "solid")+
  xlab("Dose [Gy]")+ylab("Surviving fraction")+annotation_logticks(sides = "l")+
  scale_x_continuous(breaks = c(0,0.5,1,1.5,2,2.5))+
  labs(title = "Surviving fraction: protons (prediction from C ions)",
       colour = "LET [keV/um]",
       shape = "Data set") +
  guides(linetype="none")+
  theme_bw()+theme(plot.title = element_text(size=10, color="black"),
                   axis.title.x = element_text(size=10, color="black"),
                   axis.title.y = element_text(size=10, color="black"),
                   axis.text.x = element_text(size=10, color="black"),
                   axis.text.y = element_text(size=10, color="black"),
                   legend.title = element_text(size=6, color="black"),
                   legend.text = element_text(size=5, color="black"),
                   legend.margin = margin(-0.3,0,0,0, unit="cm"),
                   legend.position = "right")

# SIMULATED LET
energy_list_He<-c("4.76","20.5","27.4","36.1","40.3","42","43.2","42.7","38.8","30.3","18.5","9.57")  
energy_list_H<-c("1.23","7.08","11.0","13.1","14.0","15.5","17.2","18.2","19.3","20.3","21.8","22.7")
energy_list_C<-c("21.3","55.3","65.1","66.0","66.4","66.0","64.2","60.5","53.7","39.6","28.3","25.0")

# SIMULATED yD FOR PLOT
yD_He<-c("7.99","14.7","18.2","29.9","45.4","54.3","64.1","71.7","79.3","84.9","88.9","89")  
yD_H<-c("3.13","4.45","5.88","9.08","10.6","13.3","16.0","17.6","19.0","20.3","22.2","23.1")
yD_C<-c("23.1","46.2","72.3","77.0","84.0","92.3","106.1","131.4","166.3","228.3","282.2","308.0")

###############################################################################
# PLOT DISTRIBUTIONS

p0 <- p0X$`0.4-6`
p <- rbind(data.frame(X = c(0:(length(p0[[9]][[94]])-1)),
                      p = as.vector(p0[[9]][[94]]),
                      Type = rep("9",length(p0[[9]][[94]]))),
           data.frame(X = c(0:(length(p0[[11]][[94]])-1)),
                      p = as.vector(p0[[11]][[94]]),
                      Type = rep("11",length(p0[[11]][[94]]))),
           data.frame(X = c(0:(length(p0[[12]][[94]])-1)),
                      p = as.vector(p0[[12]][[94]]),
                      Type = rep("12",length(p0[[12]][[94]])))) %>% 
  ggplot(aes(X,p,color=Type))+geom_line()+theme_bw()+scale_color_manual(values = c(cb_a))
ggplotly(p)


p0 <- rbind(fn_dose_all[[4]][[9]][[5]],fn_dose_all[[4]][[11]][[13]],
            fn_dose_all[[4]][[12]][[13]]) %>% 
  mutate(Type = c(rep("9",nrow(fn_dose_all[[4]][[9]][[2]])),
                  rep("11",nrow(fn_dose_all[[4]][[11]][[8]])),
                  rep("12",nrow(fn_dose_all[[4]][[12]][[13]])))) %>% 
  ggplot(aes(z,zfz,color=Type))+geom_line()+scale_x_log10()+theme_bw()+scale_color_manual(values = c(cb_a))
ggplotly(p0)