library(readr)
percent10_corner_5000_abs_biases <- read_csv("Downloads/percent10_corner_5000_abs_biases.csv")
percent20_corner_5000_abs_biases <- read_csv("Downloads/percent20_corner_5000_abs_biases.csv")
percent30_corner_5000_abs_biases <- read_csv("Downloads/percent30_corner_5000_abs_biases.csv")
percent40_corner_5000_abs_biases <- read_csv("Downloads/percent40_corner_5000_abs_biases.csv")
percent50_corner_5000_abs_biases <- read_csv("Downloads/percent50_corner_5000_abs_biases.csv")
percent64_corner_5000_abs_biases <- read_csv("Downloads/percent64_corner_5000_abs_biases.csv")

percent10_corner_5000_abs_biases

corner_abs_biases<-rbind(percent10_corner_5000_abs_biases,percent20_corner_5000_abs_biases,
                         percent30_corner_5000_abs_biases,percent40_corner_5000_abs_biases,
                         percent50_corner_5000_abs_biases,percent64_corner_5000_abs_biases)

corner_abs_biases$percent<-c(rep(10,20),rep(20,20),rep(30,20),rep(40,20),
                             rep(50,20),rep(64,20))



##draw the plot of absolute biases of percent 
#find the maxium of absolute biases
alpha_abs_max=max(corner_abs_biases$alpha_abs_bias)+0.005
alpha_sd_abs_max=max(corner_abs_biases$alpha_sd_abs_bias)+0.005
beta_abs_max=max(corner_abs_biases$beta_abs_bias)+0.05
beta_sd_abs_max=max(corner_abs_biases$beta_sd_abs_bias)+0.005


library(ggplot2)
##10 percent plots
p10_alpha<-ggplot(data=percent10_corner_5000_abs_biases ,aes(x=c(seq(1,20)),y=alpha_abs_bias))+
  geom_line()+
  ylim(0,alpha_abs_max)+
  xlab("Index")+ylab("Alpha Absolute Biases")+
  theme_classic()+
  theme(panel.border = element_rect(color = "black",
                                    fill = NA,
                                    linewidth = 0.5))
  
ggsave(path="/Users/thetnyein/Documents/MSc_Research/results_5000_pics/k_per_corner/percent10","p10_alpha_corner_5000.pdf")
p10_alpha_sd<-ggplot(data=percent10_corner_5000_abs_biases ,aes(x=c(seq(1,20)),y=alpha_sd_abs_bias))+
  geom_line()+
  ylim(0,alpha_sd_abs_max)+
xlab("Index")+ylab("Alpha SD Absolute Biases")+
  theme_classic()+
  theme(panel.border = element_rect(color = "black",
                                    fill = NA,
                                  linewidth = 0.5))

ggsave(path="/Users/thetnyein/Documents/MSc_Research/results_5000_pics/k_per_corner/percent10","p10_alpha_sd_corner_5000.pdf")

p10_beta<-ggplot(data=percent10_corner_5000_abs_biases ,aes(x=c(seq(1,20)),y=beta_abs_bias))+
  geom_line()+
  ylim(0,beta_abs_max)+
xlab("Index")+ylab("Beta Absolute Biases")+
  theme_classic()+
  theme(panel.border = element_rect(color = "black",
                                    fill = NA,
                                    linewidth = 0.5))

ggsave(path="/Users/thetnyein/Documents/MSc_Research/results_5000_pics/k_per_corner/percent10","p10_beta_corner_5000.pdf")

p10_beta_sd<-ggplot(data=percent10_corner_5000_abs_biases ,aes(x=c(seq(1,20)),y=beta_sd_abs_bias))+
  geom_line()+
  ylim(0,beta_sd_abs_max)+
xlab("Index")+ylab("Beta SD Absolute Biases")+
  theme_classic()+
  theme(panel.border = element_rect(color = "black",
                                    fill = NA,
                                    linewidth = 0.5))

ggsave(path="/Users/thetnyein/Documents/MSc_Research/results_5000_pics/k_per_corner/percent10","p10_beta_sd_corner_5000.pdf")


##20 percent plots
p20_alpha<-ggplot(data=percent20_corner_5000_abs_biases ,aes(x=c(seq(1,20)),y=alpha_abs_bias))+
  geom_line()+
  ylim(0,alpha_abs_max)+
  xlab("Index")+ylab("Alpha Absolute Biases")+
  theme_classic()+
  theme(panel.border = element_rect(color = "black",
                                    fill = NA,
                                    linewidth = 0.5))

ggsave(path="/Users/thetnyein/Documents/MSc_Research/results_5000_pics/k_per_corner/percent20","p20_alpha_corner_5000.pdf")
p20_alpha_sd<-ggplot(data=percent20_corner_5000_abs_biases ,aes(x=c(seq(1,20)),y=alpha_sd_abs_bias))+
  geom_line()+
  ylim(0,alpha_sd_abs_max)+
  xlab("Index")+ylab("Alpha SD Absolute Biases")+
  theme_classic()+
  theme(panel.border = element_rect(color = "black",
                                    fill = NA,
                                    linewidth = 0.5))

ggsave(path="/Users/thetnyein/Documents/MSc_Research/results_5000_pics/k_per_corner/percent20","p20_alpha_sd_corner_5000.pdf")

p20_beta<-ggplot(data=percent20_corner_5000_abs_biases ,aes(x=c(seq(1,20)),y=beta_abs_bias))+
  geom_line()+
  ylim(0,beta_abs_max)+
  xlab("Index")+ylab("Beta Absolute Biases")+
  theme_classic()+
  theme(panel.border = element_rect(color = "black",
                                    fill = NA,
                                    linewidth = 0.5))

ggsave(path="/Users/thetnyein/Documents/MSc_Research/results_5000_pics/k_per_corner/percent20","p20_beta_corner_5000.pdf")

p20_beta_sd<-ggplot(data=percent20_corner_5000_abs_biases ,aes(x=c(seq(1,20)),y=beta_sd_abs_bias))+
  geom_line()+
  ylim(0,beta_sd_abs_max)+
  xlab("Index")+ylab("Beta SD Absolute Biases")+
  theme_classic()+
  theme(panel.border = element_rect(color = "black",
                                    fill = NA,
                                    linewidth = 0.5))

ggsave(path="/Users/thetnyein/Documents/MSc_Research/results_5000_pics/k_per_corner/percent20","p20_beta_sd_corner_5000.pdf")

##30 percent plots
p30_alpha<-ggplot(data=percent30_corner_5000_abs_biases ,aes(x=c(seq(1,20)),y=alpha_abs_bias))+
  geom_line()+
  ylim(0,alpha_abs_max)+
  xlab("Index")+ylab("Alpha Absolute Biases")+
  theme_classic()+
  theme(panel.border = element_rect(color = "black",
                                    fill = NA,
                                    linewidth = 0.5))

ggsave(path="/Users/thetnyein/Documents/MSc_Research/results_5000_pics/k_per_corner/percent30","p30_alpha_corner_5000.pdf")
p30_alpha_sd<-ggplot(data=percent30_corner_5000_abs_biases ,aes(x=c(seq(1,20)),y=alpha_sd_abs_bias))+
  geom_line()+
  ylim(0,alpha_sd_abs_max)+
  xlab("Index")+ylab("Alpha SD Absolute Biases")+
  theme_classic()+
  theme(panel.border = element_rect(color = "black",
                                    fill = NA,
                                    linewidth = 0.5))

ggsave(path="/Users/thetnyein/Documents/MSc_Research/results_5000_pics/k_per_corner/percent30","p30_alpha_sd_corner_5000.pdf")

p30_beta<-ggplot(data=percent30_corner_5000_abs_biases ,aes(x=c(seq(1,20)),y=beta_abs_bias))+
  geom_line()+
  ylim(0,beta_abs_max)+
  xlab("Index")+ylab("Beta Absolute Biases")+
  theme_classic()+
  theme(panel.border = element_rect(color = "black",
                                    fill = NA,
                                    linewidth = 0.5))

ggsave(path="/Users/thetnyein/Documents/MSc_Research/results_5000_pics/k_per_corner/percent30","p30_beta_corner_5000.pdf")

p30_beta_sd<-ggplot(data=percent30_corner_5000_abs_biases ,aes(x=c(seq(1,20)),y=beta_sd_abs_bias))+
  geom_line()+
  ylim(0,beta_sd_abs_max)+
  xlab("Index")+ylab("Beta SD Absolute Biases")+
  theme_classic()+
  theme(panel.border = element_rect(color = "black",
                                    fill = NA,
                                    linewidth = 0.5))

ggsave(path="/Users/thetnyein/Documents/MSc_Research/results_5000_pics/k_per_corner/percent30","p30_beta_sd_corner_5000.pdf")

##40 percent plots
p40_alpha<-ggplot(data=percent40_corner_5000_abs_biases ,aes(x=c(seq(1,20)),y=alpha_abs_bias))+
  geom_line()+
  ylim(0,alpha_abs_max)+
  xlab("Index")+ylab("Alpha Absolute Biases")+
  theme_classic()+
  theme(panel.border = element_rect(color = "black",
                                    fill = NA,
                                    linewidth = 0.5))

ggsave(path="/Users/thetnyein/Documents/MSc_Research/results_5000_pics/k_per_corner/percent40","p40_alpha_corner_5000.pdf")
p40_alpha_sd<-ggplot(data=percent40_corner_5000_abs_biases ,aes(x=c(seq(1,20)),y=alpha_sd_abs_bias))+
  geom_line()+
  ylim(0,alpha_sd_abs_max)+
  xlab("Index")+ylab("Alpha SD Absolute Biases")+
  theme_classic()+
  theme(panel.border = element_rect(color = "black",
                                    fill = NA,
                                    linewidth = 0.5))

ggsave(path="/Users/thetnyein/Documents/MSc_Research/results_5000_pics/k_per_corner/percent40","p40_alpha_sd_corner_5000.pdf")

p40_beta<-ggplot(data=percent40_corner_5000_abs_biases ,aes(x=c(seq(1,20)),y=beta_abs_bias))+
  geom_line()+
  ylim(0,beta_abs_max)+
  xlab("Index")+ylab("Beta Absolute Biases")+
  theme_classic()+
  theme(panel.border = element_rect(color = "black",
                                    fill = NA,
                                    linewidth = 0.5))

ggsave(path="/Users/thetnyein/Documents/MSc_Research/results_5000_pics/k_per_corner/percent40","p40_beta_corner_5000.pdf")

p40_beta_sd<-ggplot(data=percent40_corner_5000_abs_biases ,aes(x=c(seq(1,20)),y=beta_sd_abs_bias))+
  geom_line()+
  ylim(0,beta_sd_abs_max)+
  xlab("Index")+ylab("Beta SD Absolute Biases")+
  theme_classic()+
  theme(panel.border = element_rect(color = "black",
                                    fill = NA,
                                    linewidth = 0.5))

ggsave(path="/Users/thetnyein/Documents/MSc_Research/results_5000_pics/k_per_corner/percent40","p40_beta_sd_corner_5000.pdf")


##50 percent plots

p50_alpha<-ggplot(data=percent50_corner_5000_abs_biases ,aes(x=c(seq(1,20)),y=alpha_abs_bias))+
  geom_line()+
  ylim(0,alpha_abs_max)+
  xlab("Index")+ylab("Alpha Absolute Biases")+
  theme_classic()+
  theme(panel.border = element_rect(color = "black",
                                    fill = NA,
                                    linewidth = 0.5))

ggsave(path="/Users/thetnyein/Documents/MSc_Research/results_5000_pics/k_per_corner/percent50","p50_alpha_corner_5000.pdf")
p50_alpha_sd<-ggplot(data=percent50_corner_5000_abs_biases ,aes(x=c(seq(1,20)),y=alpha_sd_abs_bias))+
  geom_line()+
  ylim(0,alpha_sd_abs_max)+
  xlab("Index")+ylab("Alpha SD Absolute Biases")+
  theme_classic()+
  theme(panel.border = element_rect(color = "black",
                                    fill = NA,
                                    linewidth = 0.5))

ggsave(path="/Users/thetnyein/Documents/MSc_Research/results_5000_pics/k_per_corner/percent50","p50_alpha_sd_corner_5000.pdf")

p50_beta<-ggplot(data=percent50_corner_5000_abs_biases ,aes(x=c(seq(1,20)),y=beta_abs_bias))+
  geom_line()+
  ylim(0,beta_abs_max)+
  xlab("Index")+ylab("Beta Absolute Biases")+
  theme_classic()+
  theme(panel.border = element_rect(color = "black",
                                    fill = NA,
                                    linewidth = 0.5))

ggsave(path="/Users/thetnyein/Documents/MSc_Research/results_5000_pics/k_per_corner/percent50","p50_beta_corner_5000.pdf")

p50_beta_sd<-ggplot(data=percent50_corner_5000_abs_biases ,aes(x=c(seq(1,20)),y=beta_sd_abs_bias))+
  geom_line()+
  ylim(0,beta_sd_abs_max)+
  xlab("Index")+ylab("Beta SD Absolute Biases")+
  theme_classic()+
  theme(panel.border = element_rect(color = "black",
                                    fill = NA,
                                    linewidth = 0.5))

ggsave(path="/Users/thetnyein/Documents/MSc_Research/results_5000_pics/k_per_corner/percent50","p50_beta_sd_corner_5000.pdf")


##64 percent plots
p64_alpha<-ggplot(data=percent64_corner_5000_abs_biases ,aes(x=c(seq(1,20)),y=alpha_abs_bias))+
  geom_line()+
  ylim(0,alpha_abs_max)+
  xlab("Index")+ylab("Alpha Absolute Biases")+
  theme_classic()+
  theme(panel.border = element_rect(color = "black",
                                    fill = NA,
                                    linewidth = 0.5))

ggsave(path="/Users/thetnyein/Documents/MSc_Research/results_5000_pics/k_per_corner/percent64","p64_alpha_corner_5000.pdf")
p64_alpha_sd<-ggplot(data=percent64_corner_5000_abs_biases ,aes(x=c(seq(1,20)),y=alpha_sd_abs_bias))+
  geom_line()+
  ylim(0,alpha_sd_abs_max)+
  xlab("Index")+ylab("Alpha SD Absolute Biases")+
  theme_classic()+
  theme(panel.border = element_rect(color = "black",
                                    fill = NA,
                                    linewidth = 0.5))

ggsave(path="/Users/thetnyein/Documents/MSc_Research/results_5000_pics/k_per_corner/percent64","p64_alpha_sd_corner_5000.pdf")

p64_beta<-ggplot(data=percent64_corner_5000_abs_biases ,aes(x=c(seq(1,20)),y=beta_abs_bias))+
  geom_line()+
  ylim(0,beta_abs_max)+
  xlab("Index")+ylab("Beta Absolute Biases")+
  theme_classic()+
  theme(panel.border = element_rect(color = "black",
                                    fill = NA,
                                    linewidth = 0.5))

ggsave(path="/Users/thetnyein/Documents/MSc_Research/results_5000_pics/k_per_corner/percent64","p64_beta_corner_5000.pdf")

p64_beta_sd<-ggplot(data=percent64_corner_5000_abs_biases ,aes(x=c(seq(1,20)),y=beta_sd_abs_bias))+
  geom_line()+
  ylim(0,beta_sd_abs_max)+
  xlab("Index")+ylab("Beta SD Absolute Biases")+
  theme_classic()+
  theme(panel.border = element_rect(color = "black",
                                    fill = NA,
                                    linewidth = 0.5))

ggsave(path="/Users/thetnyein/Documents/MSc_Research/results_5000_pics/k_per_corner/percent64","p64_beta_sd_corner_5000.pdf")


