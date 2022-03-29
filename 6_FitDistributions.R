library(fitdistrplus)
library(truncdist)
options(warn=-1)

chemicals <- read.csv("~/Desktop/Dropbox/03_Brevity_infoquimicos/scrapPherobase/chemicals.csv")

convertstring <-function(string){
  if (string == "totaluses"){return("Semiochemicals")}
  else if (string == "P"){return("Pheromones")}
  else if (string == "Alelo"){return("Alelochemicals")}
  else if (string == "A"){return("Attractans")}
  else if (string == "Al"){return("Allomones")}
  else if (string == "K"){return("Kairomones")}
  else if (string == "Sy"){return("Synomones")}
}


### FIT
goodnesoffit <- function(datos){
  fit1 <- fitdist(datos, "lnorm", method = "mle", start = list(meanlog=0, sdlog=1))
  gof1<- gofstat(fit1)
  fit2 <- fitdist(datos, "gamma", method = "mle")
  gof2<- gofstat(fit2)
  fit3 <- fitdist(datos, "weibull", method = "mle")
  gof3<- gofstat(fit3)
  return(round(c(fit1$loglik/fit1$n, gof1$ks, fit2$loglik/fit2$n, gof2$ks, fit3$loglik/fit3$n, gof3$ks), digits=4))
}

vector = c()
vector <- c(vector, "Infochemical \t& LND loglik MW \t& LND loglik ks \t& Gamma loglik MW \t& Gamma loglik ks \t& Weibull loglik MW \t& Weibull loglik ks \t \\\\ \\hline")

for (v in c("totaluses", "P", "Alelo", "A", "Al", "K", "Sy")){
  if (v == "totaluses"){
    subset <- chemicals$MW} else {
    subset <- subset(chemicals, chemicals[[v]] != 0)$MW
    }
  result <- goodnesoffit(subset)
  vector <- c(vector, paste(convertstring(v),"&", round(result[1],3), "&", round(result[2],2), "&",
                            round(result[3],3),"&", round(result[4],2),"&",
                            round(result[5],3),"&", round(result[6],2), "\\\\ \\hline", sep = "\t"))
}


lapply("", write, "~/Desktop/Dropbox/03_Brevity_infoquimicos/scrapPherobase/results/Fit_MW.txt", append=FALSE, ncolumns=1000)
lapply(vector, write, "~/Desktop/Dropbox/03_Brevity_infoquimicos/scrapPherobase/results/Fit_MW.txt", append=TRUE, ncolumns=1000)


#######################################################
# Carbons #############################################
vector = c()
vector <- c(vector, "Infochemical \t& LND loglik C \t& LND loglik ks \t& Gamma loglik C \t& Gamma loglik ks \t& Weibull loglik C \t& Weibull loglik ks \t \\\\ \\hline")

for (v in c("totaluses", "P", "Alelo", "A", "Al", "K", "Sy")){
  if (v == "totaluses"){
    subset <- chemicals$Cs} else {
      subset <- subset(chemicals, chemicals[[v]] != 0)$Cs
    }
  result <- goodnesoffit(subset)
  vector <- c(vector, paste(convertstring(v),"&", round(result[1],3), "&", round(result[2],2), "&",
                            round(result[3],3),"&", round(result[4],2),"&",
                            round(result[5],3),"&", round(result[6],2), "\\\\ \\hline", sep = "\t"))
}


lapply("", write, "~/Desktop/Dropbox/03_Brevity_infoquimicos/scrapPherobase/results/Fit_C.txt", append=FALSE, ncolumns=1000)
lapply(vector, write, "~/Desktop/Dropbox/03_Brevity_infoquimicos/scrapPherobase/results/Fit_C.txt", append=TRUE, ncolumns=1000)


##### ################################3
## FIT LND #############
vector = c()
vector <- c(vector, "Infochemical \t& LND     fit Mu MW \t& LND     fit Sigma MW \t& LND     fit Mu C \t& LND     fit Sigma C \t \\\\ \\hline")

for (v in c("totaluses", "P", "Alelo", "A", "Al", "K", "Sy")){
  if (v == "totaluses"){
    subset1 <- chemicals$MW
    subset2 <- chemicals$Cs} else {
      subset1 <- subset(chemicals, chemicals[[v]] != 0)$MW
      subset2 <- subset(chemicals, chemicals[[v]] != 0)$Cs
    }
  
  fit1 <- fitdist(subset1, "lnorm", method = "mle", start = list(meanlog=0, sdlog=1))
  fit2 <- fitdist(subset2, "lnorm", method = "mle", start = list(meanlog=0, sdlog=1))
  
  vector <- c(vector, paste(convertstring(v),"&$", round(fit1$estimate[1],3), "\\pm", round(fit1$sd[1],3), "$&$",
                            round(fit1$estimate[2],3),"\\pm", round(fit1$sd[2],3), "$&$",
                            round(fit2$estimate[1],3),"\\pm", round(fit2$sd[1],3), "$&$",
                            round(fit2$estimate[2],3),"\\pm", round(fit2$sd[2],3), "$\\\\ \\hline", sep = "\t"))
}


lapply("", write, "~/Desktop/Dropbox/03_Brevity_infoquimicos/scrapPherobase/results/Fit_LND.txt", append=FALSE, ncolumns=1000)
lapply(vector, write, "~/Desktop/Dropbox/03_Brevity_infoquimicos/scrapPherobase/results/Fit_LND.txt", append=TRUE, ncolumns=1000)


## FIT LND #############
vector = c()
vector <- c(vector, "Infochemical \t& LND     fit k MW \t& LND     fit theta MW \t& LND     fit k C \t& LND     fit theta C \t \\\\ \\hline")

for (v in c("totaluses", "P", "Alelo", "A", "Al", "K", "Sy")){
  if (v == "totaluses"){
    subset1 <- chemicals$MW
    subset2 <- chemicals$Cs} else {
      subset1 <- subset(chemicals, chemicals[[v]] != 0)$MW
      subset2 <- subset(chemicals, chemicals[[v]] != 0)$Cs
    }
  
  fit1 <- fitdist(subset1, "gamma", method = "mle")
  fit2 <- fitdist(subset2, "gamma", method = "mle")
  
  vector <- c(vector, paste(convertstring(v),"&$", round(fit1$estimate[1],3), "\\pm", round(fit1$sd[1],3), "$&$",
                            round(fit1$estimate[2],3),"\\pm", round(fit1$sd[2],3), "$&$",
                            round(fit2$estimate[1],3),"\\pm", round(fit2$sd[1],3), "$&$",
                            round(fit2$estimate[2],3),"\\pm", round(fit2$sd[2],3), "$\\\\ \\hline", sep = "\t"))
}


lapply("", write, "~/Desktop/Dropbox/03_Brevity_infoquimicos/scrapPherobase/results/Fit_gamma.txt", append=FALSE, ncolumns=1000)
lapply(vector, write, "~/Desktop/Dropbox/03_Brevity_infoquimicos/scrapPherobase/results/Fit_gamma.txt", append=TRUE, ncolumns=1000)
