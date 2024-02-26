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


path_or<-"~/GitHub/FragmentGSM2"
setwd(path_or)
source("utilities_fragments.R")

# CELL LINE
cell<-"H460"
cell<-"H1437"

# IONE DA SELEZIONARE
# ion<-"He"
ion<-"H"
# ion<-"C"

energy_list <- get_energy_list(ion)

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
# load("/home/user/Scrivania/Dottorato/Codici_R/Survival_uniform_dose_sampling_H1437.RData")
load("~/GitHub/FragmentGSM2/Survival_all_noIntercept.RData")
pide_data<-survival



# CALCULATIONS
radius<-c("SPECTRA_8um")
#coloumn<-as.character(1:12)
coloumn<-c("10")
hist<-list(); yF<-list(); yD<-list(); ystar<-list()

scorer<-list()
file <- "C10"
name <- paste0("SurfaceFile.phsp")
x_ <- fread(name)
scorer[[file]]<-x_

h <- compute_histogram_TOPAS_y(x_)
hist[[file]] <- h$hist
yF[[file]]<-h$zF
yD[[file]]<-h$zD
y0<-150
ystar[[file]]<-(y0*y0*sum(((1-exp(-(h$H$x^2)/(y0*y0)))*h$H$BinWidth*h$H$fy_bw)))/(yF[[file]]*sum(h$H$BinWidth*h$H$fy_bw))

hist$C10 %>% ggplot(aes(x,ydy)) + geom_line() + scale_x_log10()

damage_type <- "DSBsites"
kappa_all <- compute_kappa(ion, damage_type)



scorer[[file]]<-x_

p0X<-list(); p0Y<-list()
rd<-0.8; Rn<-8
#domain_size<-c(0.35,0.4,0.45,0.5,0.6,0.8,1)
#domain_size_cell<-c(5,6,7,8)
domain_size<-0.8
domain_size_cell<-6
coloumn<-c("C10")

plot_hist <- FALSE

radius<-c("SPECTRA_08um")

bin_start_cell <- 10^-5
bin_end_cell <- 10^3
bin_number_cell <- abs(log10(bin_start_cell) - log10(bin_end_cell))*10

bin_start_domain <- 10^-2
bin_end_domain <- 10^5
bin_number_domain <- abs(log10(bin_start_domain) - log10(bin_end_domain))*10

eps<-0.001; rho<-0.001
sim<-10^3

h <- compute_histogram_TOPAS_cell(x_, Rn, bin_start_cell, bin_end_cell, bin_number_cell)

h$H %>% 
  ggplot(aes(x,ydy)) + geom_line() + scale_x_log10()

dose_zn<-h$H$x

ncores<-detectCores()-1
cl <- parallel::makeCluster(ncores)
doParallel::registerDoParallel(cl)

for (rd in domain_size) {
  for (Rn in domain_size_cell) {
    
    ctr<-1
    
    print(rd)
    
    p0X[[paste0(as.character(rd),"-",as.character(Rn))]]<-list()
    p0Y[[paste0(as.character(rd),"-",as.character(Rn))]]<-list()
    
    ####ciclo rispetto alle cartelle usate
    for (file in coloumn) {
      
      ctr_col <- as.numeric(strsplit(coloumn,"C")[[1]][2])
      print(file)
      
      dose <- get_dose(ion, cell)
      
      h <- compute_histogram_TOPAS_domain(scorer[[file]], rd, bin_start_domain, bin_end_domain, bin_number_domain)
      
      if(plot_hist){
        h$H %>% 
          ggplot(aes(x,ydy)) + geom_line() + scale_x_log10()
      }
      
      # Ad hoc correction for Kappa in Poisson distribution
      kappaLET <- kappa_all[[ctr_col]]    
      
      # Kappa normalisation for volume
      kappa <- kappaLET/((Rn*Rn*Rn)/(rd*rd*rd)) ; l<-kappa*10^-3
      
      H_p <- h$H
      zF_s <- h$zF
      zD <- h$zD
      C <- h$C
      
      df<-data.frame(z = H_p$x, fz = H_p$fy_bw_norm, zfz = H_p$yfy)
      
      p0x<-list(); p0y<-list()
      
      for (zn in dose_zn) {
        
        print(zn)
        
        x_sample <- parSapply(cl, 1:sim,FUN=get_p0x_single,
                             df=df, kappa=kappa, zn=zn, zF_s=zF_s)
        
        y_sample <- parSapply(cl, 1:sim,FUN=get_p0y_single,
                             df=df,l=l,zn=zn,zF_s=zF_s)
        
        if(max(x_sample) > 0){
          hist_<-hist(x_sample,breaks=seq(-0.5,max(x_sample)+0.5,by=1)) 
          p0x[[as.character(zn)]]<-hist_$counts/length(x_sample)
        }else{
          p0x[[as.character(zn)]]<-1
        }
        
        if(max(y_sample) > 0){
          hist_<-hist(y_sample,breaks=seq(-0.5,max(y_sample)+0.5,by=1)) 
          p0y[[as.character(zn)]]<-hist_$counts/length(y_sample)
        }else{
          p0y[[as.character(zn)]]<-1
        }
      }
      
      p0X[[paste0(as.character(rd),"-",as.character(Rn))]][[file]]<-p0x
      p0Y[[paste0(as.character(rd),"-",as.character(Rn))]][[file]]<-p0y
      
      ctr<-ctr+1
    }
  }
}

stopCluster(cl)
registerDoSEQ()

# save(p0X,file="p0X_MD.RData")
# save(p0Y,file="p0Y_MD.RData")

radius<-c("SPECTRA_8um")

#domain_size<-c(0.35,0.4,0.45,0.5,0.6,0.8,1)
#domain_size_cell<-c(5,6,7,8)
coloumn<-c("C10")

fn_dose_all<-list()

ncores<-detectCores()-1
cl <- parallel::makeCluster(ncores)
doParallel::registerDoParallel(cl)

for (Rn in domain_size_cell) {
  ctr<-1
  
  print(Rn)
  
  fn_dose_all[[as.character(Rn)]]<-list()
  
  for (file in coloumn) {
    print(file)
    
    dose <- get_dose(ion, cell)
    
    h <- compute_histogram_TOPAS_cell(x_, Rn, bin_start_cell, bin_end_cell, bin_number_cell)
    
    if(plot_hist){
      h$H %>% 
        ggplot(aes(x,ydy)) + geom_line() + scale_x_log10()
    }
    
    hist <- h$hist
    H_p <- h$H
    zF_s <- h$zF
    zD <- h$zD
    C <- h$C
    
    df <- data.frame(z = H_p$x, fz = H_p$fy_bw_norm, zfz = H_p$yfy)
    
    kappaLET <- kappa_all[[ctr_col]]    
    kappa <- kappaLET/((Rn*Rn*Rn)/(rd*rd*rd)) ; l<-kappa*10^-3 
    
    fn_dose<-list()
    for (zn in dose) {
      
      print(zn)
      
      fn_sample <- parSapply(cl, 1:sim, FUN=get_fn, df=df, zn=zn, zF_s=zF_s)
      
      h <- compute_histogram_TOPAS_cell_sampled(fn_sample, Rn, bin_start_cell, bin_end_cell, bin_number_cell)
      
      fn_dose[[as.character(zn)]]<-data.frame(z=h$H$x,fz=h$H$fy_bw_norm,zfz=h$H$yfy)
      
      if(plot_hist){
        h$H %>% 
          ggplot(aes(x,ydy)) + geom_line() + scale_x_log10()
      }
    }
    
    fn_dose_all[[as.character(Rn)]][[file]]<-fn_dose
    
    ctr<-ctr+1
  }
}

stopCluster(cl)
registerDoSEQ()




##compare spectra
if(FALSE){
  h <- compute_histogram_TOPAS_cell(x_, Rn, bin_start_cell, bin_end_cell, bin_number_cell)
  
  plot_H <- h$H %>% 
    mutate(Type = "Original")
  
  plot_H %>% 
    ggplot(aes(x,ydy)) + geom_line() + scale_x_log10()
  
  
  ncores<-detectCores()-1
  cl <- parallel::makeCluster(ncores)
  doParallel::registerDoParallel(cl)
  
  
  fn_sample <- parSapply(cl, 1:sim, FUN=get_fn_single, df=df, zn=zn, zF_s=zF_s)
  h <- compute_histogram_TOPAS_cell_sampled(fn_sample, Rn, bin_start_cell, bin_end_cell, bin_number_cell)
  
  h$H %>% 
    mutate(Type = "Sampled") %>% 
    rbind(plot_H) %>% 
    ggplot(aes(x, ydy, color = Type)) + geom_line() + scale_x_log10()
}


