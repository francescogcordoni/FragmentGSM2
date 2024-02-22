
get_p0x_single_OLD<-function(i,df=df,kappa=kappa,zn=zn,zF_s=zF_s){
  
  nu_sample<-rpois(1,zn/zF_s)
  
  if(nu_sample>0){
    z_sample<-sum(sample(df$z, size = nu_sample, replace = TRUE, prob = df$zfz))
  }else{
    z_sample<-0.01
  }
  
  rpois(1,kappa*z_sample)
}


get_p0x_single<-function(i,df=df,kappa=kappa,zn=zn,zF_s=zF_s){
  rpois(1,kappa*sum(sample(df$z, size = rpois(1,zn/zF_s), replace = TRUE, prob = df$zfz)))
}

get_p0x_single_Poisson<-function(i,df=df,kappa=kappa,zn=zn,zF_s=zF_s){
  rpois(1,kappa*zn)
}

get_p0x_single_FLASH<-function(i,z=z,kappa=kappa){
  rpois(1,kappa*z)
}

get_p0y_single_OLD<-function(i,df=df,l=l,zn=zn,zF_s=zF_s){
  
  nu_sample<-rpois(1,zn/zF_s)
  
  if(nu_sample>0){
    z_sample<-sum(sample(df$z, size = nu_sample, replace = TRUE, prob = df$zfz))
  }else{
    z_sample<-0.01
  }
  
  rpois(1,l*z_sample)
}

get_p0y_single<-function(i,df=df,l=l,zn=zn,zF_s=zF_s){
  rpois(1,l*sum(sample(df$z, size = rpois(1,zn/zF_s), replace = TRUE, prob = df$zfz)))
}

get_p0y_single_Poisson<-function(i,df=df,l=l,zn=zn,zF_s=zF_s){
  rpois(1,l*zn)
}

get_p0y_single_FLASH<-function(i,z=z,l=l){
  rpois(1,l*z)
}

get_p0_single<-function(i,df=df,kappa=kappa,l=l,zn=zn,zF_s=zF_s){
  
  nu_sample<-rpois(1,zn/zF_s)
  
  if(nu_sample>0){
    z_sample<-sum(sample(df$z, size = nu_sample, replace = TRUE, prob = df$zfz))
  }else{
    z_sample<-0.01
  }
  
  x0<-rpois(1,kappa*z_sample)
  y0<-rpois(1,l*z_sample)
  
  return(list(x0=x0,y0=y0))
  
}

get_fn<-function(i,df=df,zn=zn,zF_s=zF_s){
  
  nu_sample<-rpois(1,zn/zF_s)
  
  if(nu_sample>0){
    z_sample<-sum(sample(df$z, size = nu_sample, replace = TRUE, prob = df$zfz))
  }else{
    z_sample<-0.0001001
  }
}

get_f_nu<-function(i,df=df,zn=zn,zF_s=zF_s,nu=nu){
  
  nu_sample<-nu
  
  if(nu_sample>0){
    z_sample<-sum(sample(df$z, size = nu_sample, replace = TRUE, prob = df$zfz))
  }else{
    z_sample<-0.01
  }
}

compute_histogram<-function(tp,rd){
  
  x<-tp
  
  x$z<-(x$Edep*0.204)/(4*rd*rd)
  max(x$z)
  # z<-(hist$x*rd)/md
  
  bin_start<-(0.01)
  bin_end<-1000
  bin_number<-1000
  
  br<-10^seq(log10(bin_start),log10(bin_end),length.out=bin_number)
  
  p_energy_deposition<-ggplot(x, aes(x=z)) +
    geom_histogram(position="identity", alpha=0.5,breaks=br)
  
  data<-ggplot_build(p_energy_deposition)
  
  hist<-data$data[[1]][,c(2,3,4,5,6)]
  hist$BinWidth<-abs(hist$xmax - hist$xmin)
  H<-hist[,c(2,1,6)]
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
  
  hist<-data$data[[1]][,c(2,3,4,5,6)]
  hist$BinWidth<-abs(hist$xmax - hist$xmin)
  
  B<-1/diff(log10(hist$BinWidth))[1]
  C<-log(10)*diff(log10(hist$BinWidth))[1]
  hist$fy_bw<-hist$count/hist$BinWidth
  hist$fy_bw_norm<-hist$fy_bw/(C*sum(hist$x*hist$count))
  
  hist$yfy<-hist$fy_bw*hist$x
  hist$yfy_norm<-hist$yfy/(C*sum(hist$x*hist$yfy))
  hist$ydy<-hist$fy_bw*hist$x*hist$x
  
  #yF value
  (yF_hist<-sum(hist$BinWidth*hist$yfy)/sum(hist$BinWidth*hist$fy_bw))
  
  #yD value
  (yD_hist<-sum(hist$BinWidth*hist$ydy)/sum(hist$BinWidth*hist$yfy))
  
  #yF value with H
  (zF<-sum(H$yfy*H$x)/sum(H$yfy))
  
  #yD value withH
  (zD<-sum(H$ydy*H$x)/sum(H$yfy))
  
  df<-data.frame(z=H$x,fz=H$fy_bw_norm,zfz=H$yfy_norm)
  
  return(list(H = H, C = C, zF = zF, zD = zD, hist = hist, zF_hist = yF_hist, zD_hist = yD_hist))
}

compute_histogram_TOPAS_OLD<-function(tp,rd){
  
  x<-tp
  
  x$z<-(x$V1*0.204)/(4*rd*rd)
  max(x$z)
  # z<-(hist$x*rd)/md
  
  bin_start<-(0.01)
  bin_end<-1000
  bin_number<-1000
  
  br<-10^seq(log10(bin_start),log10(bin_end),length.out=bin_number)
  
  p_energy_deposition<-ggplot(x, aes(x=z)) +
    geom_histogram(position="identity", alpha=0.5,breaks=br)
  
  data<-ggplot_build(p_energy_deposition)
  
  hist<-data$data[[1]][,c(2,3,4,5,6)]
  hist$BinWidth<-abs(hist$xmax - hist$xmin)
  H<-hist[,c(2,1,6)]
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
  
  hist<-data$data[[1]][,c(2,3,4,5,6)]
  hist$BinWidth<-abs(hist$xmax - hist$xmin)
  
  B<-1/diff(log10(hist$BinWidth))[1]
  C<-log(10)*diff(log10(hist$BinWidth))[1]
  hist$fy_bw<-hist$count/hist$BinWidth
  hist$fy_bw_norm<-hist$fy_bw/(C*sum(hist$x*hist$count))
  
  hist$yfy<-hist$fy_bw*hist$x
  hist$yfy_norm<-hist$yfy/(C*sum(hist$x*hist$yfy))
  hist$ydy<-hist$fy_bw*hist$x*hist$x
  
  #yF value
  (yF_hist<-sum(hist$BinWidth*hist$yfy)/sum(hist$BinWidth*hist$fy_bw))
  
  #yD value
  (yD_hist<-sum(hist$BinWidth*hist$ydy)/sum(hist$BinWidth*hist$yfy))
  
  #yF value with H
  (zF<-sum(H$yfy*H$x)/sum(H$yfy))
  
  #yD value withH
  (zD<-sum(H$ydy*H$x)/sum(H$yfy))
  
  df<-data.frame(z=H$x,fz=H$fy_bw_norm,zfz=H$yfy_norm)
  
  return(list(H = H, C = C, zF = zF, zD = zD, hist = hist, zF_hist = yF_hist, zD_hist = yD_hist))
}

compute_histogram_TOPAS<-function(tp,rd){
  
  x<-tp
  
  x$z<-(x$V1*0.204)/(4*rd*rd)
  
  max(x$z)
  # z<-(hist$x*rd)/md
  
  bin_start<-(0.01)
  bin_end<-1000
  bin_number<-1000
  
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

compute_histogram_TOPAS_cell<-function(tp,rd){
  
  x<-tp
  
  x$z<-(x$V1*0.204)/(4*rd*rd)
  
  max(x$z)
  # z<-(hist$x*rd)/md
  
  bin_start<-(0.0001)
  bin_end<-100
  bin_number<-300
  
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

compute_histogram_TOPAS_domain<-function(tp,rd){
  
  x<-tp
  
  x$z<-(x$V1*0.204)/(4*rd*rd)
  max(x$z)
  # z<-(hist$x*rd)/md
  
  bin_start<-(0.01)
  bin_end<-1000
  bin_number<-250
  
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

compute_histogram_TOPAS_y_OLD<-function(tp){
  
  x<-tp
  
  # z<-(hist$x*rd)/md
  
  bin_start<-(0.1)
  bin_end<-1000
  bin_number<-100
  
  br<-10^seq(log10(bin_start),log10(bin_end),length.out=bin_number)
  
  p_energy_deposition<-ggplot(x, aes(x=V1)) +
    geom_histogram(position="identity", alpha=0.5,breaks=br)
  
  data<-ggplot_build(p_energy_deposition)
  
  hist<-data$data[[1]][,c(2,3,4,5,6)]
  hist$BinWidth<-abs(hist$xmax - hist$xmin)
  H<-hist[,c(2,1,6)]
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
  
  hist<-data$data[[1]][,c(2,3,4,5,6)]
  hist$BinWidth<-abs(hist$xmax - hist$xmin)
  
  B<-1/diff(log10(hist$BinWidth))[1]
  C<-log(10)*diff(log10(hist$BinWidth))[1]
  hist$fy_bw<-hist$count/hist$BinWidth
  hist$fy_bw_norm<-hist$fy_bw/(C*sum(hist$x*hist$count))
  
  hist$yfy<-hist$fy_bw*hist$x
  hist$yfy_norm<-hist$yfy/(C*sum(hist$x*hist$yfy))
  hist$ydy<-hist$fy_bw*hist$x*hist$x
  
  #yF value
  (yF_hist<-sum(hist$BinWidth*hist$yfy)/sum(hist$BinWidth*hist$fy_bw))
  
  #yD value
  (yD_hist<-sum(hist$BinWidth*hist$ydy)/sum(hist$BinWidth*hist$yfy))
  
  #yF value with H
  (zF<-sum(H$yfy*H$x)/sum(H$yfy))
  
  #yD value withH
  (zD<-sum(H$ydy*H$x)/sum(H$yfy))
  
  df<-data.frame(z=H$x,fz=H$fy_bw_norm,zfz=H$yfy_norm)
  
  return(list(H = H, C = C, zF = zF, zD = zD, hist = hist, zF_hist = yF_hist, zD_hist = yD_hist))
}

compute_histogram_TOPAS_y<-function(tp){
  
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



compute_ystar<-function(hist,y0,yF){
  
  y0<-(y0)^2
  
  ystar<-(y0*sum(((1-exp(-(hist$x^2)/(y0)))*hist$BinWidth*hist$count)))/(yF*sum(hist$BinWidth*hist$count))
  
  return(ystar)
}

compute_initial_distribution_single_TOPAS<-function(particle,rd,kappa,l,a,b,r,nu_max,dose,tp){
  
  print(paste0("Computing particle ",particle))
  
  Histogram<-compute_histogram_TOPAS(tp,rd)
  hist<-Histogram$hist
  H_p<-Histogram$H
  zF_s<-Histogram$zF
  zD<-Histogram$zD
  C<-Histogram$C
  
  df<-data.frame(z = H_p$x, fz = H_p$fy_bw_norm, zfz = H_p$yfy)
  
  p0x<-list()
  p0y<-list()
  
  p0x[[particle]]<-list()
  p0y[[particle]]<-list()
  
  for (zn in dose) {
    print(zn)
    
    sim<-2*10^6
    
    ncores<-detectCores()-1
    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)
    x_sample <-parSapply(cl, 1:sim,FUN=get_p0x_single,df=df,kappa=kappa,
                         zn=zn,zF_s=zF_s)
    stopCluster(cl)
    registerDoSEQ()
    
    ncores<-detectCores()-1
    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)
    y_sample <-parSapply(cl, 1:sim,FUN=get_p0y_single,df=df,l=l,
                         zn=zn,zF_s=zF_s)
    stopCluster(cl)
    registerDoSEQ()
    
    p0_x_p<-c()
    p0_y_p<-c()
    for (i in 0:max(x_sample)) {
      p0_x_p[i+1]<-length(which(x_sample == i))/length(x_sample)
      p0_y_p[i+1]<-length(which(y_sample == i))/length(y_sample)
    }
    
    dose_name<-paste0("Dose = ",as.character(zn))
    
    p0x[[particle]][[dose_name]]<-p0_x_p
    p0y[[particle]][[dose_name]]<-p0_y_p
  }
  
  return(list(p0x=p0x[[particle]],p0y=p0y[[particle]]))
}

compute_all_distribution_TOPAS<-function(particle,rd1,rd2,kappa,l,a,b,r,nu_max,dose,dose_full,Nmax,tp){
  
  rd<-rd1
  nu_max<-100
  p0<-compute_initial_distribution_single_TOPAS(particle,rd,kappa,l,a,b,r,nu_max,dose_full,tp)
  
  p0x<-list()
  p0y<-list()
  p0x[[particle]]<-p0$p0x
  p0y[[particle]]<-p0$p0y
  
  p0x_fullD_rd1<-p0$p0x
  p0y_fullD_rd1<-p0$p0y
  
  rd<-rd2
  # p0<-compute_initial_distribution_single_TOPAS(particle,rd,kappa,l,a,b,r,nu_max,dose_full,tp)
  
  p0x<-list()
  p0y<-list()
  p0x[[particle]]<-p0$p0x
  p0y[[particle]]<-p0$p0y
  
  p0x_fullD_rd2<-p0$p0x
  p0y_fullD_rd2<-p0$p0y
  
  ctr<-1
  
  rd_or<-rd
  
  rd<-2.5
  
  zF_cell<-list()
  Histogram<-compute_histogram_TOPAS(tp,rd)
  hist<-Histogram$hist
  H_p<-Histogram$H
  zF_cell[[particle]]<-Histogram$zF
  zD<-Histogram$zD
  C<-Histogram$C
  
  df<-data.frame(z = H_p$x, fz = H_p$fy_bw_norm, zfz = H_p$yfy)
  
  fn<-list()
  
  ctr<-1
  fn_all<-list()
  for (zn in dose) {
    dose_name<-paste0("Dose = ",as.character(zn))
    print(dose_name)
    
    ncores<-detectCores()-1
    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)
    z_sample <-parSapply(cl, 1:(2*10^6),FUN=get_fn,df=df,zn=zn,zF_s=zF_cell[[particle]])
    stopCluster(cl)
    registerDoSEQ()
    
    z_sample<-ifelse(z_sample==0,0.01,z_sample)
    
    br<-10^seq(log10(0.01),log10(100),length.out=100)
    
    p_energy_deposition<-ggplot(data.frame(z=z_sample), aes(x=z)) +
      geom_histogram(position="identity", alpha=0.5,breaks=br)
    
    data<-ggplot_build(p_energy_deposition)
    
    hist<-data$data[[1]][,c(2,3,4,5,6)]
    hist$BinWidth<-abs(hist$xmax - hist$xmin)
    H<-hist[,c(2,1,6)]
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
    
    fn[[particle]]<-data.frame(z=H$x,fz=H$fy_bw_norm,zfz=H$yfy_norm)
    
    fn_all[[dose_name]]<-fn[[particle]]
    
    ctr<-ctr+1
  }
  
  return(list(p0x_fullD_r1 = p0x_fullD_rd1, p0y_fullD_r1 = p0y_fullD_rd1,
              p0x_fullD_r2 = p0x_fullD_rd2, p0y_fullD_r2 = p0y_fullD_rd2, fn = fn_all))
}

compute_all_distribution_TOPAS_C<-function(particle,rd1,rd2,kappa,l,a,b,r,nu_max,dose,dose_full,Nmax,tp){
  
  rd<-rd1
  nu_max<-100
  p0<-compute_initial_distribution_single_TOPAS(particle,rd,kappa,l,a,b,r,nu_max,dose_full,tp)
  
  p0x<-list()
  p0y<-list()
  p0x[[particle]]<-p0$p0x
  p0y[[particle]]<-p0$p0y
  
  p0x_fullD_rd1<-p0$p0x
  p0y_fullD_rd1<-p0$p0y
  
  rd<-rd2
  # p0<-compute_initial_distribution_single_TOPAS(particle,rd,kappa,l,a,b,r,nu_max,dose_full,tp)
  
  p0x<-list()
  p0y<-list()
  p0x[[particle]]<-p0$p0x
  p0y[[particle]]<-p0$p0y
  
  p0x_fullD_rd2<-p0$p0x
  p0y_fullD_rd2<-p0$p0y
  
  ctr<-1
  
  rd_or<-rd
  
  rd<-2.5
  
  zF_cell<-list()
  Histogram<-compute_histogram_TOPAS(tp,rd)
  hist<-Histogram$hist
  H_p<-Histogram$H
  zF_cell[[particle]]<-Histogram$zF
  zD<-Histogram$zD
  C<-Histogram$C
  
  df<-data.frame(z = H_p$x, fz = H_p$fy_bw_norm, zfz = H_p$yfy)
  
  fn<-list()
  
  ctr<-1
  fn_all<-list()
  for (zn in dose) {
    dose_name<-paste0("Dose = ",as.character(zn))
    print(dose_name)
    
    ncores<-detectCores()-1
    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)
    z_sample <-parSapply(cl, 1:(2*10^6),FUN=get_fn,df=df,zn=zn,zF_s=zF_cell[[particle]])
    stopCluster(cl)
    registerDoSEQ()
    
    z_sample<-ifelse(z_sample==0,0.01,z_sample)
    
    br<-10^seq(log10(0.01),log10(100),length.out=100)
    
    p_energy_deposition<-ggplot(data.frame(z=z_sample), aes(x=z)) +
      geom_histogram(position="identity", alpha=0.5,breaks=br)
    
    data<-ggplot_build(p_energy_deposition)
    
    hist<-data$data[[1]][,c(2,3,4,5,6)]
    hist$BinWidth<-abs(hist$xmax - hist$xmin)
    H<-hist[,c(2,1,6)]
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
    
    fn[[particle]]<-data.frame(z=H$x,fz=H$fy_bw_norm,zfz=H$yfy_norm)
    
    fn_all[[dose_name]]<-fn[[particle]]
    
    ctr<-ctr+1
  }
  
  return(list(p0x_fullD_r1 = p0x_fullD_rd1, p0y_fullD_r1 = p0y_fullD_rd1,
              p0x_fullD_r2 = p0x_fullD_rd2, p0y_fullD_r2 = p0y_fullD_rd2, fn = fn_all))
}


############################### GSM2 Simulations

gillespie_dr <- function(N, n,fz,zF,zn, ...){
  tt = 0
  x = N$M
  S = t(N$Post-N$Pre)
  u = nrow(S)
  v = ncol(S)
  tvec = vector("numeric",n)
  xmat = matrix(ncol=u,nrow=n+1)
  xmat[1,] = x
  for(i in 1:n) {
    h = N$h(x,tt, ...)
    tt = tt+rexp(1,sum(h))
    j = sample(v,1,prob=h)
    
    if(j<4){
      x = x+S[,j]
    }else{
      
    }
    
    tvec[i] = tt
    xmat[i+1,] = x
  }
  return( list(t=tvec, x=xmat))
}

N= list()
N$h=function(x,t,th=c(th1=1,th2=0.005,th3=0.6)){
  with(as.list(c(x,th)),{
    return(c(th1*x1, th2*x1*x2, th3*x2 ))
  })
}

discretise <- function(out, dt = 1, start = 0){
  events = length(out$t)
  end = out$t[events]
  len = (end - start)%/%dt + 1
  x = matrix(nrow = len, ncol = ncol(out$x))
  target = 0
  j = 1
  for(i in 1:events) {
    while (out$t[i] >= target) {
      x[j, ] = out$x[i, ]
      j = j + 1
      target = target + dt
    }
  }
  ts(x, start = 0, deltat = dt)
}

gillespied <- function(N, T=100, dt=1, ...){
  tt = 0
  n = T%/%dt
  x = N$M
  S = t(N$Post-N$Pre)
  u = nrow(S)
  v = ncol(S)
  xmat = matrix(ncol=u,nrow=n)
  i = 1
  target = 0
  repeat{
    h = N$h(x, tt, ...)
    h0 = sum(h)
    if(h0 < 1e-10)
      tt = 1e99
    else
      tt = tt+rexp(1,h0)
    while(tt >= target) {
      xmat[i,] = x
      i = i+1
      target = target+dt
      if(i > n)
        return( ts(xmat,start=0,
                   deltat=dt))
    }
    j = sample(v,1,prob=h)
    x = x+S[,j]
  }
}

StepGillespie <- function(N,p0x,p0y){
  S = t(N$Post-N$Pre)
  v = ncol(S)
  return(
    function(x0, t0, deltat, ...)
    {
      t = t0
      x = x0
      termt = t0+deltat
      repeat{
        h = N$h(x,t,...)
        h0 = sum(h)
        if(h0 < 1e-10)
          t = 1e99
        else if(h0 > 1e6) {
          t = 1e99
          warning("Hazard too big - terminating simulation!"
          )
        }
        else
          t = t+rexp(1,h0)
        if(t >= termt)
          return(x)
        j = sample(v,1,prob=h)
        if(j < 4){
          x = x+S[,j]  
        }else{
          xnew<-sample(0:(length(p0x)-1),1,prob=p0x)
          ynew<-sample(0:(length(p0y)-1),1,prob=p0y)
          
          x<-x+c(xnew,ynew)
        }
        
      }
    }
  )
}

simTs <- function(x0, t0=0, tt=100, dt=0.1, stepFun,
                  ...) {
  n = (tt-t0) %/% dt + 1
  u = length(x0)
  names = names(x0)
  mat = matrix(nrow=n,ncol=u)
  x = x0
  t = t0
  mat[1,] = x
  for(i in 2:n) {
    t = t+dt
    x = stepFun(x,t,dt,...)
    mat[i,] = x
  }
  ts(mat, start=t0, deltat=dt, names=names)
}


Csur<-function(k,x,x0,a,b,r){
  
  num<-1
  for (i_num in (x+1):x0) {
    num<-num*i_num*r
  }
  
  den1<-1
  for (i_den1 in x:(k-1)) {
    den1<-den1*((a+r)*i_den1 + b*i_den1*(i_den1-1) - (a+r)*k - b*k*(k-1))
  }
  
  den2<-1
  for (i_den2 in (k+1):x0) {
    den2<-den2*((a+r)*i_den2 + b*i_den2*(i_den2-1) - (a+r)*k - b*k*(k-1))
  }
  
  if(x0 == k){
    ris<-num/(den1)
  }else if(x == k){
    ris<-num/(den2)
  }else{
    ris<-num/(den1*den2)
  }
  
  return(ris)
}
