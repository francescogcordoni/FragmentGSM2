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

data.frame(LET = LET_list, yD = yD_list) %>% 
  ggplot(aes(yD,LET))+geom_point()+geom_line()


fit <- lm(formula = LET_list[-c(1:2)] ~ yD_list[-c(1:2)])
coefficients(fit)

yD <- 25
LET_predetto <- coefficients(fit)[1] + coefficients(fit)[2]*yD

# load("G:/Other computers/Il mio laptop/Francesco/Università/Articoli/Marta/GSM2/Validazione/Survival_all.RData")
load("G:/Other computers/Il mio laptop/Francesco/Università/Articoli/Marta/GSM2/Validazione/Survival_all_noIntercept.RData")
# cell_line <- "H460"
cell_line <- "H1437"
data <- survival[[cell_line]]

if(cell_line == "H1437"){
  ax <- 0.05; bx <- 0.041
}else if (cell_line == "H460"){
  ax <- 0.29; bx <- 0.083
}

df_plot <- data.frame(Dose = c(0:10), SF = exp(-ax*c(0:10) - bx*c(0:10)*c(0:10)),
                      Ion = "X", LET = "0", Type = "LQ")
data_1 <- data[["H"]]
for (i in 1:length(data_1)) {
  df_plot <- df_plot %>% rbind(data.frame(Dose = rep(data_1[[i]]$Dose,2), 
                              SF = c(data_1[[i]]$SF,data_1[[i]]$LQ),
                              Ion = "H", LET = names(data_1)[i], 
                              Type = c(rep("SF",nrow(data_1[[i]])),rep("LQ",nrow(data_1[[i]])))))
}

data_1 <- data[["He"]]
for (i in 1:length(data_1)) {
  df_plot <- df_plot %>% rbind(data.frame(Dose = rep(data_1[[i]]$Dose,2), 
                                          SF = c(data_1[[i]]$SF,data_1[[i]]$LQ),
                                          Ion = "He", LET = names(data_1)[i], 
                                          Type = c(rep("SF",nrow(data_1[[i]])),rep("LQ",nrow(data_1[[i]])))))
}

data_1 <- data[["C"]]
for (i in 1:length(data_1)) {
  df_plot <- df_plot %>% rbind(data.frame(Dose = rep(data_1[[i]]$Dose,2), 
                                          SF = c(data_1[[i]]$SF,data_1[[i]]$LQ),
                                          Ion = "C", LET = names(data_1)[i], 
                                          Type = c(rep("SF",nrow(data_1[[i]])),rep("LQ",nrow(data_1[[i]])))))
}

df_plot <- df_plot %>% mutate(Exp = paste0(Ion,"-",LET))

p <- df_plot %>% subset(Ion %in% c("X","He") & Dose < 5) %>% 
  ggplot(aes(Dose,SF, color = Exp, linetype = Type))+
  geom_line(linewidth = 1) + scale_y_log10() + 
  scale_color_manual(values = c25)
ggplotly(p)


cell_line <- "H460"
# cell_line <- "H1437"
data <- survival[[cell_line]]

if(cell_line == "H1437"){
  ax <- 0.05; bx <- 0.041
}else if (cell_line == "H460"){
  ax <- 0.29; bx <- 0.083
}

df_plot2 <- data.frame(Dose = c(0:10), SF = exp(-ax*c(0:10) - bx*c(0:10)*c(0:10)),
                      Ion = "X", LET = "0", Type = "LQ")
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

df_plot2 <- df_plot2 %>% mutate(Exp = paste0(Ion,"-",LET))

p2 <- df_plot2 %>% subset(Ion %in% c("X","H") & Dose < 5) %>% 
  ggplot(aes(Dose,SF, color = Exp, linetype = Type))+
  geom_line(linewidth = 1) + scale_y_log10() + 
  scale_color_manual(values = c25)
ggplotly(p2)

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

ion <- "H"
df_ <- df_plot2 %>% filter(Type == "LQ" & Ion %in% c("X",ion))
  
rbe <- c()
for (let in unique(df_$LET)) {
  x <- df_ %>% 
    filter(LET == let)
  
  rbe <- c(rbe,calculate_RBE(x, ax, bx, 0.1))
}

data.frame(LET = as.numeric(unique(df_$LET)), RBE = rbe) %>% 
  ggplot(aes(LET,RBE)) + geom_point()

df_RBE <- data.frame(LET = as.numeric(unique(df_$LET)), RBE = rbe, Type = "Original")

err <- 1

ion <- "C"
df_ <- df_plot2 %>% filter(Type == "LQ" & Ion %in% c("X",ion))

df_$Dose[which(as.numeric(df_$LET) > 12)] <- df_$Dose[which(as.numeric(df_$LET) > 12)]*(1 + err)

rbe <- c()
for (let in unique(df_$LET)) {
  x <- df_ %>% 
    filter(LET == let)
  
  rbe <- c(rbe,calculate_RBE(x, ax, bx, 0.1))
}

data.frame(LET = as.numeric(unique(df_$LET)), RBE = rbe, Type = "Corrected") %>%
  rbind(df_RBE) %>% 
  ggplot(aes(LET,RBE, color = Type)) + geom_point()




df_corr <- df_plot2 %>% subset(Ion %in% c("X","H") & Dose < 5 & Type == "LQ") %>% 
  mutate(Corr = "Corrected")

df_corr$Dose[which(as.numeric(df_corr$LET) > 12)] <- df_corr$Dose[which(as.numeric(df_corr$LET) > 12)]*(1 + err)

p3 <- df_corr %>% ggplot(aes(Dose, SF, color = Exp))+
  geom_line(linewidth = 1) + scale_y_log10() + 
  scale_color_manual(values = c25)
ggplotly(p3)

p4 <- df_plot2 %>% subset(Ion %in% c("X","H") & Dose < 5 & Type == "LQ") %>%
  mutate(Corr = "Original") %>% 
  rbind(df_corr) %>% 
  ggplot(aes(Dose, SF, color = Exp, linetype = Corr))+
  geom_line(linewidth = 1) + scale_y_log10() + 
  scale_color_manual(values = c25)
ggplotly(p4)
  
