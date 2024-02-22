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

path_or<-"/home/user/Scrivania/TESI/RISULTATI_FIT"
setwd(path_or)
source("/home/user/Scrivania/TESI/Codici_R/utilities_doserate.R")

energy_list_He<-c("3.4","8.5","12.7","22.4","36.2","44.9","54.7","63","71.1","77.6","83.6","88.7")  
energy_list_H<-c("1","2.6","4.7","7.3","8.7","11.1","13.7","15.4","16.9","18.3","20.2","21.4")
energy_list_C<-c("20.2","39.8","63.1","70.6","84.3","100.8","126.3","157.6","196.4","242.9","285.8","308.4")

# DEFINISCO I RAGGI TIPICI PER DOMINI E CELLULA PER GLI SPETTRI MICRO E LE DIMENSIONI CHE ENTRANO NEL CICLO
rd<-0.8; Rn<-8
domain_size<-c(0.35,0.4,0.45,0.5,0.6,0.8,1)
domain_size_cell<-c(5,6,7,8)
#domain_size<-c(0.45,0.5)
#domain_size_cell<-c(7,8)

# DEFINISCO IL PATH PER GLI SPETTRI MICRO CON R=8um
radius<-c("SPECTRA_8um")

# DEFINISCO IL VETTORE PER CICLARE SULLE COLONNE
coloumn<-as.character(1:12)

# FLAG <- 0 (H460) FLAG <- 1 (H1437)
flag<-1;
cell<-"H1437"

# CARICO I DATI BIOLOGICI
#load("/home/user/Scrivania/TESI/Codici_R/Survival_CORRECTED.RData")
#load("/home/user/Scrivania/TESI/Codici_R/Survival_removed_points.RData")
#load("/home/user/Scrivania/TESI/Codici_R/Survival_all.RData")
#load("/home/user/Scrivania/TESI/Codici_R/Survival_high_doses.RData")
if(flag == 0){
  load("/home/user/Scrivania/Dottorato/Codici_R/Survival_uniform_dose_sampling.RData")
}else if(flag == 1){
  
  load("/home/user/Scrivania/Dottorato/Codici_R/Survival_uniform_dose_sampling_H1437.RData")
}

pide_data<-survival

# CARICO LE DISTRIBUZIONI DI DANNO E LA MULTI-EVENT DIRETTAMENTE DENTRO IL FIT
Fit_GSM2 <- function(x) {
  
  a<-x[1]
  b<-x[2]
  r<-x[3]
  
  N<-floor((Rn^3)/(rd^3))
  
  surv<-c()
  ctr_e<-1
  
  for (ion in c("He","H","C")){
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
    
    
    if(ion=="He"){
      coloumn <-coloumn_He  
    }else if(ion=="H"){
      coloumn <-coloumn_H
    }else if(ion=="C"){
      coloumn <-coloumn_C
    }
    
    surv<-c()
    ctr_e<-1
    
    for (file in coloumn) {
      #la dose diventa la dose contenuta nel file survival rispetto alla linea cellulare, allo ione e all'energia
      if(ion=="He"){
        dose<-pide_data[[cell]]$He[[energy_list_He[ctr_e]]]$Dose  
      }else if(ion=="H"){
        dose<-pide_data[[cell]]$H[[energy_list_H[ctr_e]]]$Dose
      }else if(ion=="C"){
        dose<-pide_data[[cell]]$C[[energy_list_C[ctr_e]]]$Dose
      }
      
      ctr<-1
      surv_GSM<-c()
      
      # riga aggiunta per ciclare sullo ione
      dose_zn<-fn_dose_all[[1]][[1]][[1]]$z
      
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
  }
  # sum(((surv-S_real)^2)/S_real)
  
  sum(((surv-S_real)^2)/S_real)
  
  # return(res)
}

# CONSIDERO SOLO ALCUNE DELLE 12 COLONNE
curves_He<-c(4,6,8)  #1,2,3,4,(5),6,7,8,11,12
curves_H<-c(6,8) #1,2,3,(5),7
curves_C<-c(4,8,9) #1,(2),3,4,5,(6),7,8

#curves_He<-c(1:12)  
#curves_H<-c(1:12) 
#curves_C<-c(1:12) 

coloumn_He<-coloumn[sort(curves_He,decreasing=FALSE)]
coloumn_H<-coloumn[sort(curves_H,decreasing=FALSE)]
coloumn_C<-coloumn[sort(curves_C,decreasing=FALSE)]
energy_list_He<-energy_list_He[sort(curves_He,decreasing=FALSE)]
energy_list_H<-energy_list_H[sort(curves_H,decreasing=FALSE)]
energy_list_C<-energy_list_C[sort(curves_C,decreasing=FALSE)]

{# HELIUM IONS
ion<-"He"
if(ion=="He"){
  energy_list<-energy_list_He
}else if(ion=="H"){
  energy_list<-energy_list_H
}else if(ion=="C"){
  energy_list<-energy_list_C
}

# CARICO LA VERA SURVIVAL SF CICLANDO SULLE ENERGIE (da mettere dentro il ciclo per togliere poi alcune colonne)
S_real<-c(); err<-c()
ctr_e<-1
for(file in coloumn){
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
Dati_He <- cbind(S_real,err)

# PROTONS
ion<-"H"
if(ion=="He"){
  energy_list<-energy_list_He
}else if(ion=="H"){
  energy_list<-energy_list_H
}else if(ion=="C"){
  energy_list<-energy_list_C
}

# CARICO LA VERA SURVIVAL SF CICLANDO SULLE ENERGIE (da mettere dentro il ciclo per togliere poi alcune colonne)
S_real<-c(); err<-c()
ctr_e<-1
for(file in coloumn){
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
Dati_H <- cbind(S_real,err)

# CARBON IONS
ion<-"C"
if(ion=="He"){
  energy_list<-energy_list_He
}else if(ion=="H"){
  energy_list<-energy_list_H
}else if(ion=="C"){
  energy_list<-energy_list_C
}

# CARICO LA VERA SURVIVAL SF CICLANDO SULLE ENERGIE (da mettere dentro il ciclo per togliere poi alcune colonne)
S_real<-c(); err<-c()
ctr_e<-1
for(file in coloumn){
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
Dati_C <- cbind(S_real,err)
}
# RIGHE PER FIT INSIEME
Dati_He_H_C <- rbind(Dati_He,Dati_H,Dati_C)

# DEFINISCO I PARAMETRI DA FITTARE
a_param<-c();b_param<-c();r_param<-c()
surv_GSM2<-list()
surv_GSM2_He<-list()
surv_GSM2_H<-list()
surv_GSM2_C<-list()

# INIZIO CICLO
for (rd in domain_size) {
  for (Rn in domain_size_cell) {
    
    print(paste0("rd=", rd," Rn=", Rn))
    
    # Cambiato lower bound di a e b (originale 0.001)
    #result_<-optim(c(0.0008997624,0.06421062,2.714594),lower=c(0.0001,0.0001,1),upper=c(1,1,16), 
    #               Fit_GSM2, method="L-BFGS-B")
    
    result_<-optim(c(0.01,0.01,3),lower=c(0.0001,0.0001,1),upper=c(1,1,16), 
                   Fit_GSM2, method="L-BFGS-B")
    
    (a<-result_$par[1])
    (b<-result_$par[2])
    (r<-result_$par[3])
    
    # SCRIVO PARAMETRI OTTENUTI DA FIT PRECEDENTI PER PREDIRE ALTRE CURVE (H460)
    #rd=0.8;
    #Rn=6;
    #a=0.0008997624
    #b=0.06421062
    #r=2.714594
    
    # SCRIVO PARAMETRI OTTENUTI DA FIT PRECEDENTI PER PREDIRE ALTRE CURVE (H1437)
    #rd=0.45;
    #Rn=5;
    #a=0.010715712
    #b=0.010708522
    #r=2.999760
    
    # stampa parametri
    print(paste0("a=", a))
    print(paste0("b=", b))
    print(paste0("r=", r))
    
    a_param<-c(a_param,a)
    b_param<-c(b_param,b)
    r_param<-c(r_param,r)
    
    rc<-Rn
    N<-floor((Rn^3)/(rd^3))
    
    for (ion in c("He","H","C")){
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
      
      if(ion=="He"){
        coloumn <-coloumn_He  
      }else if(ion=="H"){
        coloumn <-coloumn_H
      }else if(ion=="C"){
        coloumn <-coloumn_C
      }
      
      surv<-c()
      ctr_e<-1
      
      for (file in coloumn) {
        print(paste0("file=", file))
        
        #la dose diventa la dose contenuta nel file survival rispetto alla linea cellulare, allo ione e all'energia
        if(ion=="He"){
          dose<-pide_data[[cell]]$He[[energy_list_He[ctr_e]]]$Dose  
        }else if(ion=="H"){
          dose<-pide_data[[cell]]$H[[energy_list_H[ctr_e]]]$Dose
        }else if(ion=="C"){
          dose<-pide_data[[cell]]$C[[energy_list_C[ctr_e]]]$Dose
        }
        
        ctr<-1
        surv_GSM<-c()
        
        # riga aggiunta per ciclare sullo ione
        dose_zn<-fn_dose_all[[1]][[1]][[1]]$z
        
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
      if(ion=="He"){
        surv_He<-surv  
      }else if(ion=="H"){
        surv_H<-surv
      }else if(ion=="C"){
        surv_C<-surv
      }
      
      if(ion=="He"){
        surv_GSM2_He[[paste0(as.character(rd),"-",as.character(Rn))]]<-surv  
      }else if(ion=="H"){
        surv_GSM2_H[[paste0(as.character(rd),"-",as.character(Rn))]]<-surv
      }else if(ion=="C"){
        surv_GSM2_C[[paste0(as.character(rd),"-",as.character(Rn))]]<-surv
      }
      
    }
  }
}

surv <- c(surv_He,surv_H,surv_C)
surv<-na.omit(surv)

# DOSE_LIST_ALL ED ENERGY_LIST_ALL
dose_list_all<-c()
dose_list_all_He<-c()
dose_list_all_H<-c()
dose_list_all_C<-c()
energy_list_all_He<-c()
energy_list_all_H<-c()
energy_list_all_C<-c()

ctr_e<-1
for (file in coloumn_He) {
  print(file)
  dose_list_all_He<-c(dose_list_all_He,pide_data[[cell]]$He[[energy_list_He[ctr_e]]]$Dose)
  energy_list_all_He<-c(energy_list_all_He,rep(energy_list_He[ctr_e],length(pide_data[[cell]]$He[[energy_list_He[ctr_e]]]$Dose)))
  ctr_e<-ctr_e+1
}

ctr_e<-1
for (file in coloumn_H) {
  print(file)
  dose_list_all_H<-c(dose_list_all_H,pide_data[[cell]]$H[[energy_list_H[ctr_e]]]$Dose)
  energy_list_all_H<-c(energy_list_all_H,rep(energy_list_H[ctr_e],length(pide_data[[cell]]$H[[energy_list_H[ctr_e]]]$Dose)))
  ctr_e<-ctr_e+1
}

ctr_e<-1
for (file in coloumn_C) {
  print(file)
  dose_list_all_C<-c(dose_list_all_C,pide_data[[cell]]$C[[energy_list_C[ctr_e]]]$Dose)
  energy_list_all_C<-c(energy_list_all_C,rep(energy_list_C[ctr_e],length(pide_data[[cell]]$C[[energy_list_C[ctr_e]]]$Dose)))
  ctr_e<-ctr_e+1
}

energy_list_all_He<-paste0("He - ",energy_list_all_He)
energy_list_all_H<-paste0("p - ",energy_list_all_H)
energy_list_all_C<-paste0("C - ",energy_list_all_C)
energy_list_all<-c(energy_list_all_He,energy_list_all_H,energy_list_all_C)
dose_list_all<-c(dose_list_all_He,dose_list_all_H,dose_list_all_C)
  
# ANALISI CHI2
chi_sq<-c()
sizes<-c()
rd_vec<-c()
Rn_vec<-c()
Dati_He_H_C<-as.data.frame(Dati_He_H_C)
S_real<-Dati_He_H_C$S_real
err<-Dati_He_H_C$err

# DA USARE surv_GSM2_C ed He (DUE RIGHE PER SELEZIONARE rd ed Rn DA CHI2)
domain_size<-c(0.35,0.4,0.45,0.5,0.6,0.8,1)
domain_size_cell<-c(5,6,7,8)
domain_size <- 0.45
domain_size_cell <- 5
for (rd in domain_size) {
  for (Rn in domain_size_cell) {
    sizes<-c(sizes,paste0(as.character(rd),"-",as.character(Rn)))
    rd_vec<-c(rd_vec,rd)
    Rn_vec<-c(Rn_vec,Rn)
    survival_He<-na.omit(surv_GSM2_He[[paste0(as.character(rd),"-",as.character(Rn))]])
    survival_H<-na.omit(surv_GSM2_H[[paste0(as.character(rd),"-",as.character(Rn))]])
    survival_C<-na.omit(surv_GSM2_C[[paste0(as.character(rd),"-",as.character(Rn))]])
    surv<-c(survival_He,survival_H,survival_C)
    chi_sq<-c(chi_sq,sum(((surv-S_real)^2)/S_real))
    
    plot_df<-data.frame(Dose = rep(dose_list_all,(2*length(surv))/length(dose_list_all)),
                        Survival = c(surv,S_real),
                        Energy = rep(energy_list_all,2),
                        Type = c(rep("GSM2",length(surv)),rep("Exp",length(surv))))
  }
}

data.frame(Chi=chi_sq,Size=sizes,a=a_param,b=b_param,r=r_param,ratio=a_param/b_param) %>% View()
Fit_parameters <- data.frame(Chi=chi_sq,rd=rd_vec,Rn=Rn_vec,a=a_param,b=b_param,r=r_param,ratio=a_param/b_param)
write.csv(Fit_parameters, "/home/user/Scrivania/Dottorato/RISULTATI_FIT/H1437/Fit_global_H1437_3.csv", row.names=FALSE)

############ ERROR BAR ############
plot_s<-plot_df
#barra errore esperimento
error_up<-c()
error_down<-c()
error_up_<-c()
error_down_<-c()
for (file in as.numeric(coloumn_He)) {
  error_up<-c(error_up,survival[[cell]]$He[[file]]$LQ + 0.5*survival[[cell]]$He[[file]]$ErrSF)
  error_down<-c(error_down,survival[[cell]]$He[[file]]$LQ - 0.5*survival[[cell]]$He[[file]]$ErrSF)
}
for (file in as.numeric(coloumn_H)) {
  error_up<-c(error_up,survival[[cell]]$H[[file]]$LQ + 0.5*survival[[cell]]$H[[file]]$ErrSF)
  error_down<-c(error_down,survival[[cell]]$H[[file]]$LQ - 0.5*survival[[cell]]$H[[file]]$ErrSF)
}
for (file in as.numeric(coloumn_C)) {
  error_up<-c(error_up,survival[[cell]]$C[[file]]$LQ + 0.5*survival[[cell]]$C[[file]]$ErrSF)
  error_down<-c(error_down,survival[[cell]]$C[[file]]$LQ - 0.5*survival[[cell]]$C[[file]]$ErrSF)
}
#barra errore modello coincide con la survival DEL MODELLO (non ho errore da aggiungere)
error_up_<-c(error_up_,surv)
error_down_<-c(error_down_,surv)

colnames(plot_s)[3]<-"Ion - Energy [MeV/u]"

plot_s<-plot_s %>% mutate(Up = c(error_up_,error_up), Down = c(error_down_,error_down))

plot_s$Up<-ifelse(plot_s$Up > 1,1-0.001,plot_s$Up)

# PLOT CON BARRE D'ERRORE
ggplotly(ggplot(plot_s,aes(Dose,Survival,color=`Ion - Energy [MeV/u]`,linetype=Type,shape=Type))+
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


# NO COLOR SCALE
ggplotly(ggplot(plot_s,aes(Dose,Survival,color=`Ion - Energy [MeV/u]`,linetype=Type,shape=Type))+
           geom_line(linewidth=0.5)+geom_point(size=1)+
           scale_y_log10()+scale_linetype_manual(values=c("dashed","solid"))+
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

ggplot(plot_s,aes(Dose,Survival,color=`Ion - Energy [MeV/u]`,linetype=Type,shape=Type))+
  geom_line(linewidth=0.5)+geom_point(size=1.5)+
  scale_y_log10()+scale_color_manual(values=c25)+scale_linetype_manual(values=c("blank","solid"))+
  scale_shape_manual(labels = c("Experiments", "GSM2 model"), values=c(16, 18))+
  geom_errorbar(aes(ymax = Up, ymin = Down),width=0.01,linewidth=0.5,linetype = "solid")+
  xlab("Dose [Gy]")+ylab("Surviving fraction")+annotation_logticks(sides = "l")+
  scale_x_continuous(breaks = c(0,0.5,1,1.5,2,2.5))+
  labs(title = "Surviving fraction: global fit with different radiation quality",
       colour = "LET [keV/um]",
       shape = "Data set") +
  guides(linetype="none")+
  theme_bw()+theme(plot.title = element_text(size=10, color="black"),
                   axis.title.x = element_text(size=10, color="black"),
                   axis.title.y = element_text(size=10, color="black"),
                   axis.text.x = element_text(size=10, color="black"),
                   axis.text.y = element_text(size=10, color="black"),
                   legend.title = element_text(size=5, color="black"),
                   legend.text = element_text(size=5, color="black"),
                   legend.margin = margin(-0.1,0,0,0, unit="cm"),
                   legend.position = "right")



########### CODICE PER FIT DA DATI RADIOBIOLOGICI ESPERIMENTI ###########

# CREO FILE .RDATA TOGLIENDO ALCUNI PUNTI

###### leggo i file Excel, divisi per cellula, ione e per LET, li metto tutti in una stessa lista (lista di liste di dataframe)
survival<-list()

survival[["H460"]]<-list()
survival[["H1437"]]<-list()

survival[["H460"]][["H"]]<-list()
survival[["H1437"]][["H"]]<-list()

survival[["H460"]][["He"]]<-list()
survival[["H1437"]][["He"]]<-list()

survival[["H460"]][["C"]]<-list()
survival[["H1437"]][["C"]]<-list()

LET_H <- c("1","2.6",'4.7',"7.3","8.7","11.1","13.7","15.4","16.9","18.3","20.2","21.4")
LET_He <- c("3.4","8.5",'12.7',"22.4","36.2","44.9","54.7","63","71.1","77.6","83.6","88.7")
LET_C <- c("20.2","39.8",'63.1',"70.6","84.3","100.8","126.3","157.6","196.4","242.9","285.8","308.4")

setwd("/home/user/Scrivania/TESI/Survival_data")

# load("Survival_all.RData")
# CICLO PER PRODUZIONE SURVIVAL.RDATA (MODIFICARE SCELTA PUNTI DA TOGLIERE DA DATI RAW)
for (cell in c("H460","H1437")) {
  print(cell)
  for (ion in c("H","He","C")) {
    file_ctr <- paste0(cell,"_",ion,"_ctr.xlsx")
    df_ctr<-read_excel(file_ctr, col_names = FALSE)
    
    err_ctr<-sd(df_ctr[1,])
    
    if(ion == "H"){
      for (let in LET_H) {
        pointsToremove<-c()
        print(paste0(ion,"-",let))
        file <- paste0(cell,"_",ion,"_",let,".xlsx")
        df<-read_excel(file, col_names = FALSE)
        
        df[is.na(df)]<-0
        
        dose<-df$...1
        
        media_surv<-c()
        err_surv<-c()
        
        #points selection for fitting surviving fraction with LQ model
        all_points<-c(1:dim(df)[1])
        pointsTofit<-all_points
        
        #remove some points to improve the goodness of the LQ model fit
        # RIMOZIONE PUNTI POCO SIGNIFICATIVI
        if(let == LET_H[2]){
          pointsToremove<-c(12)
          pointsTofit<-all_points[-pointsToremove]
        }
        if(let == LET_H[3]){
          pointsToremove<-c(9,10)
          pointsTofit<-all_points[-pointsToremove]
        }
        if(let == LET_H[4]){
          pointsToremove<-c(8)
          pointsTofit<-all_points[-pointsToremove]
        }
        if(let == LET_H[5]){
          pointsToremove<-c(7,8)
          pointsTofit<-all_points[-pointsToremove]
        }
        
        for (i in pointsTofit) {
          media_surv[i]<-mean(as.numeric(df[i,1+which(df[i,-1] > 0)])) #media
          err_surv[i]<-sd(as.numeric(df[i,1+which(df[i,-1] > 0)]))+media_surv[i]*err_ctr #deviazione standard
        }
        
        #sottraggo punti se necessario
        if (let == LET_H[2] | let ==LET_H[3] | let == LET_H[4] | let == LET_H[5]) {
          media_surv<-media_surv[-pointsToremove]
          err_surv<-err_surv[-pointsToremove]
          dose<-dose[-pointsToremove]
        }
        
        ############ fitto le survival per avere una curva, i punti sono troppo su e giù altrimenti
        
        s<-log(media_surv)
        dose2<-dose*dose
        
        mod <- lm(s ~ dose+dose2-1, weights = 1/sqrt(err_surv))
        pred <- exp(mod$coefficients[1]*dose+mod$coefficients[2]*dose2)
        
        survival[[cell]][[ion]][[let]]<-data.frame(Dose = dose, SF=media_surv, ErrSF=err_surv, LQ=pred)
      }#end let H
    }else if(ion == "He"){
      for (let in LET_He) {
        pointsToremove<-c()
        print(paste0(ion,"-",let))
        file <- paste0(cell,"_",ion,"_",let,".xlsx")
        df<-read_excel(file, col_names = FALSE)
         
        df[is.na(df)]<-0
        
        dose<-df$...1
        
        media_surv<-c()
        err_surv<-c()
        
        #points selection for fitting surviving fraction with LQ model
        all_points<-c(1:dim(df)[1])
        pointsTofit<-all_points
        
        #remove some points to improve the goodness of the LQ model fit
        # RIMOZIONE PUNTI POCO SIGNIFICATIVI
        if(let == LET_He[2]){
          pointsToremove<-c(13)
          pointsTofit<-all_points[-pointsToremove]
        }
        if(let == LET_He[3]){
          pointsToremove<-c(11)
          pointsTofit<-all_points[-pointsToremove]
        }
        if(let == LET_He[4]){
          pointsToremove<-c(8,9)
          pointsTofit<-all_points[-pointsToremove]
        }else if(let == LET_He[7]){
          pointsToremove<-c(4,5)
          pointsTofit<-all_points[-pointsToremove]
        }else if(let == LET_He[8]){
          pointsToremove<-c(4)
          pointsTofit<-all_points[-pointsToremove]
        }else if(let == LET_He[9]){
          pointsToremove<-c(5)
          pointsTofit<-all_points[-pointsToremove]
        }else if(let == LET_He[10]){
          pointsToremove<-c(8,9)
          pointsTofit<-all_points[-pointsToremove]
        }else if(let == LET_He[11]){
          pointsToremove<-c(12,13)
          pointsTofit<-all_points[-pointsToremove]
        }
        
        for (i in pointsTofit) {
          media_surv[i]<-mean(as.numeric(df[i,1+which(df[i,-1] > 0)])) #media
          err_surv[i]<-sd(as.numeric(df[i,1+which(df[i,-1] > 0)]))+media_surv[i]*err_ctr #deviazione standard
        }
        
        #sottraggo punti se necessario
        if (let == LET_He[2] | let == LET_He[3] | let == LET_He[4] | let == LET_He[7] | let == LET_He[8] | let == LET_He[9] | let == LET_He[10] | let == LET_He[11]) {
          media_surv<-media_surv[-pointsToremove]
          err_surv<-err_surv[-pointsToremove]
          dose<-dose[-pointsToremove]
        }
        
        ############ fitto le survival per avere una curva, i punti sono troppo su e giù altrimenti
        s<-log(media_surv)
        dose2<-dose*dose
        
        mod <- lm(s ~ dose+dose2-1, weights = 1/sqrt(err_surv))
        pred <- exp(mod$coefficients[1]*dose+mod$coefficients[2]*dose2)
        
        survival[[cell]][[ion]][[let]]<-data.frame(Dose = dose, SF=media_surv, ErrSF=err_surv, LQ=pred)
      }
    }else if(ion == "C"){
      for (let in LET_C) {
        pointsToremove<-c()
        print(paste0(ion,"-",let))
        file <- paste0(cell,"_",ion,"_",let,".xlsx")
        df<-read_excel(file, col_names = FALSE)
        
        df[is.na(df)]<-0
        
        dose<-df$...1
        
        media_surv<-c()
        err_surv<-c()
        
        #points selection for fitting surviving fraction with LQ model
        all_points<-c(1:dim(df)[1])
        pointsTofit<-all_points
        
        #points selection for fitting surviving fraction with LQ model
        all_points<-c(1:dim(df)[1])
        pointsTofit<-all_points
        
        #remove some points to improve the goodness of the LQ model fit
        # RIMOZIONE PUNTI POCO SIGNIFICATIVI
        if(let == LET_C[2]){
          pointsToremove<-c(13,14)
          pointsTofit<-all_points[-pointsToremove]
        }else if(let == LET_C[3]){
          pointsToremove<-c(7)
          pointsTofit<-all_points[-pointsToremove]
        }else if(let == LET_C[4]){
          pointsToremove<-c(9,10)
          pointsTofit<-all_points[-pointsToremove]
        }else if(let == LET_C[5]){
          pointsToremove<-c(7,8)
          pointsTofit<-all_points[-pointsToremove]
        }else if(let == LET_C[6]){
          pointsToremove<-c(7,8)
          pointsTofit<-all_points[-pointsToremove]
        }else if(let == LET_C[7]){
          pointsToremove<-c(6,7)
          pointsTofit<-all_points[-pointsToremove]
        }else if(let == LET_C[8]){
          pointsToremove<-c(7,8)
          pointsTofit<-all_points[-pointsToremove]
        }else if(let == LET_C[9]){
          pointsToremove<-c(7,8)
          pointsTofit<-all_points[-pointsToremove]
        }else if(let == LET_C[10]){
          pointsToremove<-c(9,10,11)
          pointsTofit<-all_points[-pointsToremove]
        }else if(let == LET_C[11]){
          pointsToremove<-c(13,14,15)
          pointsTofit<-all_points[-pointsToremove]
        }
        
        for (i in pointsTofit) {
          media_surv[i]<-mean(as.numeric(df[i,1+which(df[i,-1] > 0)])) #media
          err_surv[i]<-sd(as.numeric(df[i,1+which(df[i,-1] > 0)]))+media_surv[i]*err_ctr #deviazione standard
        }
        
        #sottraggo punti se necessario
        if (let == LET_C[2] | let == LET_C[3] | let == LET_C[4] | let == LET_C[5] | let == LET_C[6] | let == LET_C[7] | let == LET_C[8] | let == LET_C[9] | let == LET_C[10] | let == LET_C[11]) {
          media_surv<-media_surv[-pointsToremove]
          err_surv<-err_surv[-pointsToremove]
          dose<-dose[-pointsToremove]
        }
        
        ############ fitto le survival per avere una curva, i punti sono troppo su e giù altrimenti
        
        s<-log(media_surv)
        dose2<-dose*dose
        
        mod <- lm(s ~ dose+dose2-1, weights = 1/sqrt(err_surv))
        pred <- exp(mod$coefficients[1]*dose+mod$coefficients[2]*dose2)
        
        survival[[cell]][[ion]][[let]]<-data.frame(Dose = dose, SF=media_surv, ErrSF=err_surv, LQ=pred)
      }#end let C
    }#end ion
  }#end ion
}#end cell

# SALVO I DATI
setwd("/home/user/Scrivania/TESI/Codici_R")
save(survival, file = "Survival_CORRECTED.RData")

# Codice per eventuale plot delle curve di Survival
ggplotly(ggplot(data.frame(Dose = c(survival$H460$He[[1]]$Dose,survival$H460$He[[2]]$Dose,survival$H460$He[[3]]$Dose,survival$H460$He[[4]]$Dose,survival$H460$He[[5]]$Dose,survival$H460$He[[6]]$Dose),
                           LQ = c(survival$H460$He[[1]]$LQ,survival$H460$He[[2]]$LQ,survival$H460$He[[3]]$LQ,survival$H460$He[[4]]$LQ,survival$H460$He[[5]]$LQ,survival$H460$He[[6]]$LQ),
                           Survival = c(survival$H460$He[[1]]$SF,survival$H460$He[[2]]$SF,survival$H460$He[[3]]$SF,survival$H460$He[[4]]$SF,survival$H460$He[[5]]$SF,survival$H460$He[[6]]$SF),
                           Type = c(rep("1",nrow(survival$H460$He[[1]])),rep("2.6",nrow(survival$H460$He[[2]])),
                                    rep("4.7",nrow(survival$H460$He[[3]])),rep("7.3",nrow(survival$H460$He[[4]])),rep("8.7",nrow(survival$H460$He[[5]])),rep("11.1",nrow(survival$H460$He[[6]])))),
                aes(Dose,Survival,color=Type))+geom_line()+scale_y_log10())


###### constraint su beta #####

Fit_LQ <- function(x) {
  
  a<-x[1]
  b<-x[2]
  
  pred<- - a*dose_ - b*dose_*dose_
  
  sum(((pred-real)^2)*(1/sqrt(err_surv)))
  
  # return(res)
}
