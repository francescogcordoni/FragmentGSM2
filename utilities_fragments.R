# DEFINISCO ENERGY_LIST PER OGNI IONE
get_energy_list <- function(ion){
  if(ion == "He"){
    energy_list <- c("3.4","8.5","12.7","22.4","36.2","44.9","54.7","63","71.1","77.6","83.6","88.7")  
  }else if(ion == "H"){
    energy_list <- c("1","2.6","4.7","7.3","8.7","11.1","13.7","15.4","16.9","18.3","20.2","21.4")
  }else if(ion == "C"){
    energy_list <- c("20.2","39.8","63.1","70.6","84.3","100.8","126.3","157.6","196.4","242.9","285.8","308.4")
  }else{
    energy_list <- c()
  }
  
  return(energy_list)
}

compute_histogram_TOPAS_y <- function(tp){
  
  x<-tp
  
  # z<-(hist$x*rd)/md
  
  bin_start<-(0.1)
  bin_end<-10000
  bin_number<-1000
  
  br<-10^seq(log10(bin_start),log10(bin_end),length.out=bin_number)
  
  hist_<-hist(x$V1[which(x$V1>bin_start & x$V1<bin_end)],breaks=br) #use hist() and specify your breaks
  hist<-data.frame(count=hist_$counts,x=0.5*(br[-1]+br[-length(br)]),xmin=br[-length(br)],xmax=br[-1])
  hist$BinWidth<-abs(hist$xmax - hist$xmin)
  
  H<-hist
  
  B<-1/diff(log10(H$BinWidth))[1]
  
  # B<-1/0.05
  C<-log(10)*diff(log10(H$BinWidth))[1]
  H$fy_bw<-H$count/H$BinWidth
  H$fy_bw_norm<-H$fy_bw/(C*sum(H$x*H$fy_bw))
  
  H$yfy<-H$fy_bw_norm*H$x
  H$yfy_norm<-H$yfy/(C*sum(H$x*H$yfy))
  H$ydy<-H$yfy_norm*H$x
  H$ydy_norm<-H$ydy/(C*sum(H$x*H$ydy))
  
  (zF<-sum(H$yfy*H$x)/sum(H$yfy))
  
  #yD value withH
  (zD<-sum(H$ydy*H$x)/sum(H$yfy))
  
  return(list(C = C, hist = H, zF = zF, zD = zD))
}

compute_kappa <- function(ion, damage_type, col_number){
  
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
  
  # PROTONS
  if(ion == "H"){
    # COMPUTE KAPPA VALUE
    LET_list <- c(1.23,7.08,11.0,13.1,14.0,15.5,17.2,18.2,19.3,20.3,21.8,22.7)
    yD_list <- c(4.99,5.53,6.71,9.75,11.4,14.2,17.2,19.0,20.7,22.1,24.2,25.7)
    
    fit <- lm(formula = LET_list[-c(1:2)] ~ yD_list[-c(1:2)])
    # Adjust Kappa 
    energy_list[col_number] <- coefficients(fit)[1] + coefficients(fit)[2]*yD[[file]]
    
    # KAPPA PROTONS
    Kappa_protons<-c()
    for(counter in c(1:12)){
      Kappa_protons<-c(Kappa_protons,p1[1]+(p2[1]*as.numeric(energy_list[counter]))^p3[1])
    }
    # Adjustment
    Kappa<-9*Kappa_protons
  }
  
  # HELIUM IONS
  if(ion == "He"){
    # COMPUTE KAPPA VALUE
    LET_list <- c(1.23,7.08,11.0,13.1,14.0,15.5,17.2,18.2,19.3,20.3,21.8,22.7)
    yD_list <- c(12,17.5,21.8,32.1,47.4,56.3,66.4,74.2,82.3,88.2,93,93.9)
    
    fit <- lm(formula = LET_list[-c(1:2)] ~ yD_list[-c(1:2)])
    # Adjust Kappa 
    energy_list[col_number] <- coefficients(fit)[1] + coefficients(fit)[2]*yD[[file]]
    
    # KAPPA HELIUM IONS
    Kappa_helium<-c()
    for(counter in c(1:12)){
      Kappa_helium<-c(Kappa_helium,(p1[2]+(p2[2]*as.numeric(energy_list[counter]))^p3[2])/(1+(p4[2]*as.numeric(energy_list[counter]))^p5[2]))
    }
    # Adjustment
    Kappa<-9*Kappa_helium
  }
  
  
  # CARBON IONS
  if(ion == "C"){
    # COMPUTE KAPPA VALUE
    LET_list <- c(21.3,55.3,65.1,66.0,66.4,66.0,64.2,60.5,53.7,39.6,28.3,25.0)
    yD_list<-c(27.4,52.9,80.2,84.6,91.1,98.8,112.2,137.4,172,235.1,292.2,319.6)
    
    fit <- lm(formula = LET_list[-c(1:2)] ~ yD_list[-c(1:2)])
    # Adjust Kappa 
    energy_list[col_number] <- coefficients(fit)[1] + coefficients(fit)[2]*yD[[file]]
    
    # KAPPA CARBON IONS
    Kappa_carbon<-c()
    for(counter in c(1:12)){
      Kappa_carbon<-c(Kappa_carbon,(p1[3]+(p2[3]*as.numeric(energy_list[counter]))^p3[3])/(1+(p4[3]*as.numeric(energy_list[counter]))^p5[3]))
    }
    # Adjustment
    Kappa<-9*Kappa_carbon
  }
  
  return(Kappa)
}

compute_histogram_TOPAS_cell <- function(tp, rd, bin_start, bin_end, bin_number){
  
  x<-tp
  
  x$z<-(x$V1*0.204)/(4*rd*rd)
  
  max(x$z)
  # z<-(hist$x*rd)/md
  
  br<-10^seq(log10(bin_start),log10(bin_end),length.out=bin_number)
  
  hist_<-hist(x$z[which(x$z > bin_start & x$z < bin_end)],breaks=br) #use hist() and specify your breaks
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
  (zF<-sum(H$yfy*H$x)/sum(H$yfy))
  
  #yD value withH
  (zD<-sum(H$ydy*H$x)/sum(H$yfy))
  
  return(list(H = H, C = C, zF = zF, zD = zD))
}

compute_histogram_TOPAS_cell_sampled <- function(z, rd, bin_start, bin_end, bin_number){
  
  br<-10^seq(log10(bin_start),log10(bin_end),length.out=bin_number)
  
  hist_<-hist(z[which(z > bin_start & z < bin_end)],breaks=br) #use hist() and specify your breaks
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
  (zF<-sum(H$yfy*H$x)/sum(H$yfy))
  
  #yD value withH
  (zD<-sum(H$ydy*H$x)/sum(H$yfy))
  
  return(list(H = H, C = C, zF = zF, zD = zD))
}

get_dose <- function(ion, cell){
  pide_data[[cell]][[ion]][[as.numeric(strsplit(coloumn,"C")[[1]][2])]]$Dose
}

compute_histogram_TOPAS_domain<-function(tp, rd, bin_start, bin_end, bin_number){
  
  x<-tp
  
  x$z<-(x$V1*0.204)/(4*rd*rd)
  max(x$z)
  # z<-(hist$x*rd)/md
  
  br<-10^seq(log10(bin_start),log10(bin_end),length.out=bin_number)
  
  hist_<-hist(x$z[which(x$z > bin_start & x$z < bin_end)],breaks=br) #use hist() and specify your breaks
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
  (zF<-sum(H$yfy*H$x)/sum(H$yfy))
  
  #yD value withH
  (zD<-sum(H$ydy*H$x)/sum(H$yfy))
  
  return(list(H = H, C = C, zF = zF, zD = zD))
}

get_p0x_single<-function(i,df=df,kappa=kappa,zn=zn,zF_s=zF_s){
  nu <- rpois(1,zn/zF_s)
  if(nu > 0){
    rpois(1,kappa*sum(sample(df$z, size = nu, replace = TRUE, prob = df$zfz)))
  }else{
    0
  }
}

get_p0y_single<-function(i,df=df,l=l,zn=zn,zF_s=zF_s){
  nu <- rpois(1,zn/zF_s)
  if(nu > 0){
    rpois(1,l*sum(sample(df$z, size = nu, replace = TRUE, prob = df$zfz)))
  }else{
    0
  }
}

get_fn <- function(i, df=df, zn=zn, zF_s=zF_s){
  
  nu_sample <- rpois(1,zn/zF_s)
  
  if(nu_sample>0){
    z_sample <- sum(sample(df$z, size = nu_sample, replace = TRUE, prob = df$zfz))
  }else{
    z_sample<-df$z[1]
  }
}

get_fn_single <- function(i, df=df, zn=zn, zF_s=zF_s){
  
  nu_sample <- 1
  
  if(nu_sample>0){
    z_sample <- sum(sample(df$z, size = nu_sample, replace = TRUE, prob = df$zfz))
  }else{
    z_sample<-df$z[1]
  }
}

