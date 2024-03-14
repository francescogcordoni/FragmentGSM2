library(stringr)
library(gridExtra)
library(maptools)
library(ggplot2)
library(ggalt)
library(ggthemes)
library(tibble)
library(viridis)
library(knitr)
library(plotly)
library(egg)
library(ggpubr)
library(GillespieSSA)
library(MASS)
library(rgl)
library(zoo)
library(Sim.DiffProc)
library(doParallel)
library(ggsci)
library(deSolve)
library(pracma)
library(diffeqr)
library(parallel)
library(tidymodels)

source("G:\\Other computers\\Il mio laptop\\Francesco\\Università\\Articoli\\Marta\\GSM2\\Dose Rate\\Code\\utilities_doserate.R")

theme_set(theme_bw()+theme(plot.title = element_text(size=20, color="black"),
                           axis.title.x = element_text(size=20, color="black"),
                           axis.title.y = element_text(size=20, color="black"),
                           axis.text.x = element_text(size=20, color="black"),
                           axis.text.y = element_text(size=20, color="black"),
                           legend.title = element_blank(),
                           legend.text = element_text(size=20, color="black")))

cb_a <- c("#0072B2","#000000","#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cb_b <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cb_c <- c( "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7")

cb_nob <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
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

LET_list <- c(1.23,7.08,11.0,13.1,14.0,15.5,17.2,18.2,19.3,20.3,21.8,22.7)
yD_list <- c(4.99,5.53,6.71,9.75,11.4,14.2,17.2,19.0,20.7,22.1,24.2,25.7)
ystar_list <- c(4.99,5.53,6.71,9.75,11.4,14.2,17.2,19.0,20.7,22.1,24.2,25.7)

yD_H<-c(7.36,7.40,7.68,10.6,12.2,14.9,18,19.8,21.5,22.9,25.1,26.7)
yD_He<-c(12,17.5,21.8,32.1,47.4,56.3,66.4,74.2,82.3,88.2,93,93.9)
yD_C<-c(27.4,52.9,80.2,84.6,91.1,98.8,112.2,137.4,172,235.1,292.2,319.6)

yD_list <- c(yD_H,yD_He,yD_C)

load("C:/Users/Utente/Downloads/ystar_ions.RData")

ystar_list <- c(ystar$H,ystar$He,ystar$C)

load("G:/Other computers/Il mio laptop/Francesco/Università/Articoli/Marta/GSM2/Validazione/Survival_all_noIntercept.RData")
cell_line <- "H460"
data <- survival[[cell_line]]

if(cell_line == "H1437"){
  ax <- 0.05; bx <- 0.041
}else if (cell_line == "H460"){
  ax <- 0.29; bx <- 0.083
}

calculate_RBE <- function(x, ax, bx, p){
  
  dose <- x$Dose
  dose2 <- dose*dose 
  fit <- lm(formula = log(x$SF) ~ dose + dose2 - 1)
  
  ai <- -coefficients(fit)[1]; bi <- -coefficients(fit)[2]
  
  rbe <- calculate_RBE_single(ax, bx, ai, bi, p)
  
  return(rbe)
  
}

calculate_RBE_single <- function(ax, bx, ai, bi, p){
  
  dose <- seq(from = 0, to = 10, by=0.1)
  rbe <- dose[which.min(abs(exp(-ax*dose - bx*dose*dose) - p))[1]]/dose[which.min(abs(exp(-ai*dose - bi*dose*dose) - p))[1]]
  
  return(rbe)
}


df_plot2 <- data.frame(Dose = c(0:10), SF = exp(-ax*c(0:10) - bx*c(0:10)*c(0:10)),
                       Ion = "X", LET = "0", Type = "LQ")
names(data[["H"]])[11] <- "20.1"
data_1 <- data[["H"]]
for (i in 1:length(data_1)) {
  df_plot2 <- df_plot2 %>% rbind(data.frame(Dose = rep(data_1[[i]]$Dose,2), 
                                            SF = c(data_1[[i]]$SF,data_1[[i]]$LQ),
                                            Ion = "H", LET = names(data_1)[i], 
                                            Type = c(rep("SF",nrow(data_1[[i]])),rep("LQ",nrow(data_1[[i]])))))
}

data_1 <- data[["He"]]
for (i in 1:length(data_1)) {
  df_plot2 <- df_plot2 %>% rbind(data.frame(Dose = rep(data_1[[i]]$Dose,2), 
                                            SF = c(data_1[[i]]$SF,data_1[[i]]$LQ),
                                            Ion = "He", LET = names(data_1)[i], 
                                            Type = c(rep("SF",nrow(data_1[[i]])),rep("LQ",nrow(data_1[[i]])))))
}

data_1 <- data[["C"]]
for (i in 1:length(data_1)) {
  df_plot2 <- df_plot2 %>% rbind(data.frame(Dose = rep(data_1[[i]]$Dose,2), 
                                            SF = c(data_1[[i]]$SF,data_1[[i]]$LQ),
                                            Ion = "C", LET = names(data_1)[i], 
                                            Type = c(rep("SF",nrow(data_1[[i]])),rep("LQ",nrow(data_1[[i]])))))
}

ion <- c("H","He","C")
df_ <- df_plot2 %>% filter(Type == "LQ" & Ion %in% c("X",ion))

rbe <- c(); ion <- c()
for (let in unique(df_$LET)) {
  x <- df_ %>% 
    filter(LET == let)
  
  rbe <- c(rbe,calculate_RBE(x, ax, bx, 0.1))
  ion <- c(ion,x$Ion[1])
}

df_plot2 <- df_plot2 %>% mutate(Exp = paste0(Ion,"-",LET))
# ion <- "C"
pRBE <- data.frame(LET = as.numeric(unique(df_$LET))[-1], RBE = rbe[-1], ion = ion[-1]) %>% 
  ggplot(aes(LET,RBE,color = ion)) +
  scale_x_log10() +
  geom_point(size = 1.5) +
  scale_color_manual(values = cb_b)
ggplotly(pRBE)

rbe_exp <- rbe[-1]

p <- 0.1
rbe <- c()
for (y in ystar_list) {
  ai <- ax + (bx/pi*rd*rd)*y; bi <- bx
  dose <- seq(from = 0, to = 10, by=0.1)
  rbe <- c(rbe,dose[which.min(abs(exp(-ax*dose - bx*dose*dose) - p))[1]]/dose[which.min(abs(exp(-ai*dose - bi*dose*dose) - p))[1]])
}
rbe

Fit_MKM <- function(x) {
  
  rd <- x[1]
  p <- 0.1
  
  rbe <- c()
  for (y in ystar_list) {
    ai <- ax + (bx/pi*rd*rd)*y; bi <- bx
    dose <- seq(from = 0, to = 10, by=0.1)
    rbe <- c(rbe,dose[which.min(abs(exp(-ax*dose - bx*dose*dose) - p))[1]]/dose[which.min(abs(exp(-ai*dose - bi*dose*dose) - p))[1]])
  }
  
  err <- mean((abs(rbe-rbe_exp))/rbe_exp)
  
  return(err)
}


fit_SMKM <- function(x) {
  
  rd <- x[1]
  Rn <- x[2]
  p <- 0.1
  
  rbe <- c()
  for (ctr in 1:length(ystar_list)) {
    ai <- ax + (bx/(pi*rd*rd))*ystar_list[ctr]; bi <- bx*(ystar_list[ctr]/(pi*rd*rd))*((pi*rd*rd)/yD_list[ctr])
    
    dose <- seq(from = 0, to = 10, by=0.1)
    S_mkm <- exp(-ai*dose - bi*dose*dose)*(1 + dose*(-bi + 0.5*(ai + 2*bi*dose)^2)*(yD_list[ctr]/(pi*Rn*Rn)))
    
    rbe_mkm <- dose[which.min(abs(exp(-ax*dose - bx*dose*dose) - p))[1]]/dose[which.min(abs(S_mkm - p))[1]]
    
    rbe <- c(rbe,rbe_mkm)
  }
  
  err <- sum((abs(rbe-rbe_exp))/rbe_exp)
  
  return(err)
}

result_<-optim(c(1,7),lower=c(0.00001,4),upper=c(5,15),fit_SMKM, method="L-BFGS-B")

(rd <- result_$par[1])
(Rn <- result_$par[2])

predict_SMKM <- function(rd, Rn) {
  
  p <- 0.1
  
  rbe <- c()
  for (ctr in 1:length(ystar_list)) {
    ai <- ax + (bx/(pi*rd*rd))*ystar_list[ctr]; bi <- bx*(ystar_list[ctr]/(pi*rd*rd))*((pi*rd*rd)/yD_list[ctr])
    
    dose <- seq(from = 0, to = 10, by=0.1)
    S_mkm <- exp(-ai*dose - bi*dose*dose)*(1 + dose*(-bi + 0.5*(ai + 2*bi*dose)^2)*(yD_list[ctr]/(pi*Rn*Rn)))
    
    rbe_mkm <- dose[which.min(abs(exp(-ax*dose - bx*dose*dose) - p))[1]]/dose[which.min(abs(S_mkm - p))[1]]
    
    rbe <- c(rbe,rbe_mkm)
  }
  
  return(rbe)
}

# result_<-optim(c(1),lower=c(0.00001),upper=c(5),Fit_MKM, method="L-BFGS-B")

# (rd <- result_$par[1])

predict_MKM <- function(rd) {
  p <- 0.1
  rbe <- c()
  for (y in ystar_list) {
    ai <- ax + (bx/pi*rd*rd)*y; bi <- bx
    dose <- seq(from = 0, to = 10, by=0.1)
    rbe <- c(rbe,dose[which.min(abs(exp(-ax*dose - bx*dose*dose) - p))[1]]/dose[which.min(abs(exp(-ai*dose - bi*dose*dose) - p))[1]])
  }
  
  return(rbe)
}

rbe_MKM <- predict_SMKM(rd, Rn)
# rd <- 0.8; rbe_MKM <- predict_MKM(rd)

p <- data.frame(LET = as.numeric(unique(df_$LET))[-1], 
                yD = yD_list, ystar = ystar_list, Ion = ion[-1],
                RBE = rbe_MKM, 
                Type = "MKM") %>% 
  rbind(data.frame(LET = as.numeric(unique(df_$LET))[-1], 
                   yD = yD_list, ystar = ystar_list, Ion = ion[-1], 
                   RBE = rbe_exp, 
                   Type = "Exp")) %>% 
  ggplot(aes(LET,RBE,color = Type))+geom_point()+scale_x_log10()+
  scale_color_manual(values = cb_b) +
  facet_wrap(~ Ion)
ggplotly(p) 
p
