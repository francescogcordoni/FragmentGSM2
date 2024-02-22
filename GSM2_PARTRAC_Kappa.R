library(doParallel)
library(data.table)

# Number of Spectra to attach
Nspectra_start <- 9
Nspectra_end <- 10

# SPECTRA SELECTION
path_or<-"/home/giulio.bordieri/GSM2_fit/Spectra_protons"
setwd(path_or)
source("/home/giulio.bordieri/GSM2_fit/Fragmentation_protons/SF_curves/utilities_doserate_highstatistics.R")

# ION SELECTION
ion<-"p"

# LET SELECTION FROM ION TYPE
energy_list_He<-c("4.76","20.5","27.4","36.1","40.3","42","43.2","42.7","38.8","30.3","18.5","9.57")  
energy_list_H<-c("1.23","7.08","11.0","13.1","14.0","15.5","17.2","18.2","19.3","20.3","21.8","22.7")
energy_list_C<-c("21.3","55.3","65.1","66.0","66.4","66.0","64.2","60.5","53.7","39.6","28.3","25.0")

if(ion=="He"){
  energy_list<-c("3.4","8.5","12.7","22.4","36.2","44.9","54.7","63","71.1","77.6","83.6","88.7")
}else if(ion=="p"){
  energy_list<-c("1","2.6","4.7","7.3","8.7","11.1","13.7","15.4","16.9","18.3","20.2","21.4")
}else if(ion=="C"){
  energy_list<-c("20.2","39.8","63.1","70.6","84.3","100.8","126.3","157.6","196.4","242.9","285.8","308.4")
}

# COMPUTE KAPPA VALUE
LET_list <- c(1.23,7.08,11.0,13.1,14.0,15.5,17.2,18.2,19.3,20.3,21.8,22.7)
yD_list <- c(4.99,5.53,6.71,9.75,11.4,14.2,17.2,19.0,20.7,22.1,24.2,25.7)

fit <- lm(formula = LET_list[-c(1:2)] ~ yD_list[-c(1:2)])
#coefficients(fit)

# DEFINIZIONE KAPPA
damage_type <- "DSBsites"

if(damage_type == "DSB"){
  p1 <- c(6.8,6.8,6.8)
  p2 <- c(0.1835,0.1679,0.2052)
  p3 <- c(0.9583,0.9704,1.02)
  p4 <- c(0,0.004323,0.009922)
  p5 <- c(0,1.359,1.106)
}else if(damage_type == "cDSB"){
  p1 <- c(0.07,0.07,0.07)
  p2 <- c(0.01532,0.01015,0.006858)
  p3 <- c(2.396,1.794,1.498)
  p4 <- c(0,0.003817,0.002778)
  p5 <- c(0,3.255,2.208)
}else if(damage_type == "DSBsites"){
  p1 <- c(6.8,6.8,6.8)
  p2 <- c(0.1773,0.1471,0.156)
  p3 <- c(0.9314,1.038,0.9214)
  p4 <- c(0,0.006239,0.005245)
  p5 <- c(0,1.582,1.395)
}else if(damage_type == "DSB_direct"){
  p1 <- c(2.8,2.8,2.8)
  p2 <- c(0.07011,0.08076,0.09374)
  p3 <- c(1.231,0.816,1.076)
  p4 <- c(0,0,0.01033)
  p5 <- c(0,0,1.006)
}else if(damage_type == "DSB_indirect"){
  p1 <- c(2.2,2.2,2.2)
  p2 <- c(0.03598,0.02683,0.03152)
  p3 <- c(0.5834,0.6349,0.6538)
  p4 <- c(0,0.002725,0.002669)
  p5 <- c(0,2.022,1.114)
}else if(damage_type == "DSBsites_direct"){
  p1 <- c(2.8,2.8,2.8)
  p2 <- c(0.06901,0.06555,0.06191)
  p3 <- c(1.196,1.023,0.9903)
  p4 <- c(0,0.003748,0.004525)
  p5 <- c(0,1.763,1.369)
}else if(damage_type == "DSBsites_indirect"){
  p1 <- c(2.2,2.2,2.2)
  p2 <- c(0.035,0.02656,0.02946)
  p3 <- c(0.5841,0.6415,0.6435)
  p4 <- c(0,0.002875,0.002585)
  p5 <- c(0,1.994,1.166)
}  

# CALCULATIONS
radius<-c("SPECTRA_8um")
#coloumn<-as.character(1:12)
coloumn<-c("10")
hist<-list(); yF<-list(); yD<-list(); ystar<-list()

scorer<-list()

# mi creo l'istogramma e i yf e yd per ogni ione e energia presente in folder_list
for (file in coloumn) {
  print(file)
  
  #path<-paste0(path_or,"/",radius)
  #x_<-fread(paste0(path,"/Scorer_",as.character(file),".phsp"))
  
  if(file=="10"){
  	#x_<-fread(paste0("/home/giulio.bordieri/GSM2_fit/Fragmentation_protons/PMMA_",as.character(file),"/Spectra/Scorer_",as.character(file),"_8um.phsp"))
  	#load("/home/user/Scrivania/Dottorato/Micro_spectra_fragmentation/Fragmentation_protons/Spectra/PMMA_10/Spectrum_8um.RData")
	for (i in Nspectra_start:Nspectra_end) {
	  name <- paste0("/home/giulio.bordieri/GSM2_fit/Fragmentation_protons/PMMA_10/Scorers_spectra/Scorer_",as.character(i),"_8um.phsp")
	  if(i == Nspectra_start){
	    x_8um <- fread(name)
	  }else{
	    x_1 <- fread(name)
	    x_8um <- rbind(x_8um,x_1)
	    }
	}  
  }
  #x_<-histall
  #print(dim(x_))
  #setwd("/home/giulio.bordieri/GSM2_fit/KAPPA_PARTRAC/Protons_3e8")
  #save(x_,file="Spectrum_3e8_8um.RData")
  
  x_<-x_8um
  scorer[[file]]<-x_
  
  h<-compute_histogram_TOPAS_y(x_)
  hist[[file]]<-h$hist
  yF[[file]]<-h$zF
  yD[[file]]<-h$zD
  y0<-150
  ystar[[file]]<-(y0*y0*sum(((1-exp(-(h$H$x^2)/(y0*y0)))*h$H$BinWidth*h$H$fy_bw)))/(yF[[file]]*sum(h$H$BinWidth*h$H$fy_bw))
}

# Adjust Kappa 
energy_list_H[10] <- coefficients(fit)[1] + coefficients(fit)[2]*yD[[file]]

# KAPPA PROTONS
Kappa_protons<-c()
for(counter in c(1:12)){
  Kappa_protons<-c(Kappa_protons,p1[1]+(p2[1]*as.numeric(energy_list_H[counter]))^p3[1])
}
# Adjustment
Kappa_protons<-9*Kappa_protons

Kappa<-Kappa_protons

# carico i dati biologici
#load("/home/user/Scrivania/Dottorato/Codici_R/Survival_uniform_dose_sampling.RData")
load("/home/giulio.bordieri/GSM2_fit/Survival_uniform_dose_sampling.RData")
pide_data<-survival

# creo le distribuzioni iniziali di danno
# GSM2
radius<-c("SPECTRA_08um")
file<-coloumn[1]
#path<-paste0(path_or,"/",radius)
#x_<-fread(paste0(path,"/Scorer_",as.character(file),".phsp"))

if(file=="10"){
  #x_<-fread(paste0("/home/giulio.bordieri/GSM2_fit/Fragmentation_protons/PMMA_",as.character(file),"/Spectra/Scorer_",as.character(file),"_8um.phsp"))
  #load("/home/user/Scrivania/Dottorato/Micro_spectra_fragmentation/Fragmentation_protons/Spectra/PMMA_10/Spectrum_08um.RData")
	  for (i in Nspectra_start:Nspectra_end) {
	  name <- paste0("/home/giulio.bordieri/GSM2_fit/Fragmentation_protons/PMMA_10/Scorers_spectra/Scorer_",as.character(i),"_08um.phsp")
	  if(i == Nspectra_start){
	    x_08um <- fread(name)
	  }else{
	    x_1 <- fread(name)
	    x_08um <- rbind(x_08um,x_1)
	    }
	}
}
#setwd("/home/giulio.bordieri/GSM2_fit/KAPPA_PARTRAC/Protons_3e8")
#save(x_,file="Spectrum_3e8_08um.RData")
#x_<-histall
x_<-x_08um
scorer[[file]]<-x_

p0X<-list(); p0Y<-list()
rd<-0.8; Rn<-8
#domain_size<-c(0.35,0.4,0.45,0.5,0.6,0.8,1)
#domain_size_cell<-c(5,6,7,8)
domain_size<-0.8
domain_size_cell<-6
coloumn<-c("10")

#####caricare gli scorer da 0.8 um
radius<-c("SPECTRA_08um")

h<-compute_histogram_TOPAS_cell(x_,rd)
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
    print(file)
    
    ####la dose diventa la dose contenuta nel file survival rispetto alla linea cellulare, allo ione e all'energia
    if(ion=="He"){
      dose<-pide_data$H460$He[[energy_list[as.numeric(file)]]]$Dose  
    }else if(ion=="p"){
      dose<-pide_data$H460$H[[energy_list[as.numeric(file)]]]$Dose
    }else if(ion=="C"){
      dose<-pide_data$H460$C[[energy_list[as.numeric(file)]]]$Dose
    }
    
    #path<-paste0(path_or,"/",radius)
    #x_<-fread(paste0(path,"/Scorer_",as.character(file),".phsp"))
    
    x_<-x_08um
    #x_<-histall
    scorer[[file]]<-x_
    
    h<-compute_histogram_TOPAS_domain(x_,rd)
    
    # Ad hoc correction for Kappa in Poisson distribution
    KappaLET<-Kappa[[as.numeric(file)]]    

    # Kappa normalisation for volume
    kappa<-KappaLET/((Rn*Rn*Rn)/(rd*rd*rd)) ; l<-kappa*10^-3
    
    p0x<-list(); p0y<-list()
    
    for (zn in dose_zn) {
      
      print(zn)
      
      Histogram<-compute_histogram_TOPAS_domain(x_,rd)
      # hist<-Histogram$hist
      H_p<-Histogram$H
      zF_s<-Histogram$zF
      zD<-Histogram$zD
      C<-Histogram$C
      
      df<-data.frame(z = H_p$x, fz = H_p$fy_bw_norm, zfz = H_p$yfy)
      
      eps<-0.001;rho<-0.001;sim<-10^5      
        
      x_sample <-parSapply(cl, 1:sim,FUN=get_p0x_single,
                           df=df,kappa=kappa,zn=zn,zF_s=zF_s)
      
      y_sample <-parSapply(cl, 1:sim,FUN=get_p0y_single,
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
    
    # p0X[[paste0(as.character(rd),"-",as.character(Rn))]][[energy_list[ctr]]]<-p0x
    # p0Y[[paste0(as.character(rd),"-",as.character(Rn))]][[energy_list[ctr]]]<-p0y
    
    p0X[[paste0(as.character(rd),"-",as.character(Rn))]][[file]]<-p0x
    p0Y[[paste0(as.character(rd),"-",as.character(Rn))]][[file]]<-p0y
    
    ctr<-ctr+1
  }
  }
}

stopCluster(cl)
registerDoSEQ()

setwd("/home/giulio.bordieri/GSM2_fit/Fragmentation_protons/SF_curves")
save(p0X,file="p0X_MD.RData")
save(p0Y,file="p0Y_MD.RData")

################################### SECOND CYLCE ON CELL DOMAIN ##############################################

#####caricare gli scorer da 8 um
radius<-c("SPECTRA_8um")

#domain_size<-c(0.35,0.4,0.45,0.5,0.6,0.8,1)
#domain_size_cell<-c(5,6,7,8)
coloumn<-c("10")

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
    
    ####dose esperimento
    if(ion=="He"){
      dose<-pide_data$H460$He[[energy_list[as.numeric(file)]]]$Dose  
    }else if(ion=="p"){
      dose<-pide_data$H460$H[[energy_list[as.numeric(file)]]]$Dose
    }else if(ion=="C"){
      dose<-pide_data$H460$C[[energy_list[as.numeric(file)]]]$Dose
    }
    
    #path<-paste0(path_or,"/",radius)
    #x_<-fread(paste0(path,"/Scorer_",as.character(file),".phsp"))
    
    #x_<-histall
    x_<-x_8um	
    h<-compute_histogram_TOPAS_cell(x_,Rn)

    KappaLET<-Kappa[[as.numeric(file)]]
    kappa<-KappaLET/((Rn*Rn*Rn)/(rd*rd*rd)) ; l<-kappa*10^-3    

    fn_dose<-list()
    for (zn in dose) {
      
      print(zn)
      
      Histogram<-compute_histogram_TOPAS_cell(x_,Rn)
      hist<-Histogram$hist
      H_p<-Histogram$H
      zF_s<-Histogram$zF
      zD<-Histogram$zD
      C<-Histogram$C
      
      df<-data.frame(z = H_p$x, fz = H_p$fy_bw_norm, zfz = H_p$yfy)
      
      fn_sample <-parSapply(cl, 1:sim,FUN=get_fn,df=df,zn=zn,zF_s=zF_s)
      
      bin_start<-(0.0001)
      bin_end<-100
      bin_number<-300
      
      br<-10^seq(log10(bin_start),log10(bin_end),length.out=bin_number)
      
      hist_<-hist(fn_sample[which(fn_sample > bin_start & fn_sample < bin_end)],breaks=br) #use hist() and specify your breaks
      hist<-data.frame(count=hist_$counts,x=0.5*(br[-1]+br[-length(br)]),xmin=br[-length(br)],xmax=br[-1])
      hist$BinWidth<-abs(hist$xmax - hist$xmin)
      
      H<-hist
      #colnames(H[1,])
      B<-1/diff(log10(H$BinWidth))[1]
      
      # B<-1/0.05
      C<-log(10)*diff(log10(H$BinWidth))[1]
      H$fy_bw<-H$count/H$BinWidth
      H$fy_bw_norm<-H$fy_bw/(C*sum(H$x*H$fy_bw))
      
      H$yfy<-H$fy_bw_norm*H$x
      H$yfy_norm<-H$yfy/(C*sum(H$x*H$yfy))
      H$ydy<-H$yfy_norm*H$x
      H$ydy_norm<-H$ydy/(C*sum(H$x*H$ydy))
      
      #yF value with H
      zF<-sum(H$yfy*H$x)/sum(H$yfy)
      
      #yD value withH
      zD<-sum(H$ydy*H$x)/sum(H$yfy)

      fn_dose[[as.character(zn)]]<-data.frame(z=H$x,fz=H$fy_bw_norm,zfz=H$yfy)
    }
    
    fn_dose_all[[as.character(Rn)]][[file]]<-fn_dose
    
    ctr<-ctr+1
  }
}

stopCluster(cl)
registerDoSEQ()

setwd("/home/giulio.bordieri/GSM2_fit/Fragmentation_protons/SF_curves")
save(fn_dose_all,file="fn_dose_all_MD.RData")

# SALVATE LE DISTRIBUZIONI DI DANNO E LA MULTI-EVENT

################### fine cluster ####################
