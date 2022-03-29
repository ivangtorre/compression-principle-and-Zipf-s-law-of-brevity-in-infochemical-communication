library("Hmisc")

pvaluestring <- function(pvalue) {
  if (pvalue<=0.0001){
    pvalue<-"$p<10^{-4}$"
  } else if (pvalue<=0.001 && pvalue>0.0001){
    pvalue<-"$p<10^{-3}$"
  } else if (pvalue<=0.01 && pvalue>0.001){
      pvalue<-"$p<10^{-2}$"
  }else{
    pvalue<-round(pvalue,2)
    pvalue<-paste("$p=", pvalue, "$", sep="")
  }
  return(pvalue)
}

convertstring <-function(string){
  if (string == "totaluses"){return("Semiochemicals")}
  else if (string == "P"){return("Pheromones")}
  else if (string == "Alelo"){return("Alelochemicals")}
  else if (string == "A"){return("Attractans")}
  else if (string == "Al"){return("Allomones")}
  else if (string == "K"){return("Kairomones")}
  else if (string == "Sy"){return("Synomones")}
}

### CORRELATION TEST FOR THE MAIN #############################################################
chemicals <- read.csv("~/Desktop/Dropbox/03_Brevity_infoquimicos/scrapPherobase/chemicals.csv")
vector = c()
#vector <- c(vector, "Infochemical & Total & Pearson & Spearman & Kendall \\\\ \\hline")
vector <- c(vector, "Infochemical \t& Total use \t&Total molec  \t& Spearman MW \t& Spearman MW p value \t& Spearman Cs \t& Spearman Cs pvalue \t \\\\ \\hline")

for (v in c("totaluses", "P", "Alelo", "A", "Al", "K", "Sy")){
#  res1 <- cor.test(chemicals$MW, chemicals[[v]], method = "pearson")
  subset <- chemicals[chemicals[[v]] !=0, ]
  res2 <- cor.test(subset$MW, subset[[v]], method = "kendall", "less")
  res22 <- cor.test(subset$Cs, subset[[v]], method = "kendall", "less")
#  res3 <- cor.test(chemicals$MW, chemicals[[v]], method = "kendall")
                   
  vector <- c(vector, paste(convertstring(v),"&", sum(chemicals[[v]]),"&", nrow(subset), "&", 
#                            "$", round(res1$estimate[[1]],2), "$", pvaluestring(res1$p.value), "&",
                             round(res2$estimate[[1]],2),"&", pvaluestring(res2$p.value),"&",
                            round(res22$estimate[[1]],2),"&", pvaluestring(res22$p.value),
#                            "$", round(res3$estimate[[1]],2), "$", pvaluestring(res2$p.value), 
                              "\\\\ \\hline",
                            sep = "\t"))
}

lapply("", write, "~/Desktop/Dropbox/03_Brevity_infoquimicos/scrapPherobase/results/correlation_main_Spearman.txt", append=FALSE, ncolumns=1000)
lapply(vector, write, "~/Desktop/Dropbox/03_Brevity_infoquimicos/scrapPherobase/results/correlation_main_Spearman.txt", append=TRUE, ncolumns=1000)


#############################################################################################
### CORRELATION TEST FOR THE SUPPLEMENTARY #################################################
## MW ##################################################################
vectorMW = c()
vectorMW <- c(vectorMW, "Infochemical \t& Spearman MW \t& Spearman MW p value \t& Kendall MW \t& Kendall MW pvalue \t& Pearson MW \t& Pearson MW pvalue \t \\\\ \\hline")

for (v in c("totaluses", "P", "Alelo", "A", "Al", "K", "Sy")){
  subset <- chemicals[chemicals[[v]] !=0, ]
  res1 <- cor.test(subset$MW, subset[[v]], method = "spearman", "less")
  res2 <- cor.test(subset$MW, subset[[v]], method = "kendall", "less")
  res3 <- cor.test(subset$MW, subset[[v]], method = "pearson", "less")
  
  vectorMW <- c(vectorMW, paste(convertstring(v),"&", 
                                round(res1$estimate[[1]],2),"&", pvaluestring(res1$p.value),"&",
                                round(res2$estimate[[1]],2),"&", pvaluestring(res2$p.value),"&",
                                round(res3$estimate[[1]],2),"&", pvaluestring(res3$p.value),
                                "\\\\ \\hline", sep = "\t"))
}

lapply("", write, "~/Desktop/Dropbox/03_Brevity_infoquimicos/scrapPherobase/results/correlation_SI_MW.txt", append=FALSE, ncolumns=1000)
lapply(vectorMW, write, "~/Desktop/Dropbox/03_Brevity_infoquimicos/scrapPherobase/results/correlation_SI_MW.txt", append=TRUE, ncolumns=1000)


#############################################################################################
### CORRELATION TEST FOR THE SUPPLEMENTARY #################################################
## Carbons ##################################################################
vectorC = c()
vectorC <- c(vectorC, "Infochemical \t& Spearman C \t& Spearman C p value \t& Kendall C \t& Kendall C pvalue \t& Pearson C \t& Pearson C pvalue \t \\\\ \\hline")

for (v in c("totaluses", "P", "Alelo", "A", "Al", "K", "Sy")){
  subset <- chemicals[chemicals[[v]] !=0, ]
  res1 <- cor.test(subset$Cs, subset[[v]], method = "spearman", "less")
  res2 <- cor.test(subset$Cs, subset[[v]], method = "kendall", "less")
  res3 <- cor.test(subset$Cs, subset[[v]], method = "pearson", "less")
  
  vectorC <- c(vectorC, paste(convertstring(v),"&", 
                                round(res1$estimate[[1]],2),"&", pvaluestring(res1$p.value),"&",
                                round(res2$estimate[[1]],2),"&", pvaluestring(res2$p.value),"&",
                                round(res3$estimate[[1]],2),"&", pvaluestring(res3$p.value),
                                "\\\\ \\hline", sep = "\t"))
}

lapply("", write, "~/Desktop/Dropbox/03_Brevity_infoquimicos/scrapPherobase/results/correlation_SI_C.txt", append=FALSE, ncolumns=1000)
lapply(vectorC, write, "~/Desktop/Dropbox/03_Brevity_infoquimicos/scrapPherobase/results/correlation_SI_C.txt", append=TRUE, ncolumns=1000)

