setwd("E:/Snow back up/Desktop/CRC database info/CCLE")

library(ggpubr)
library(tidyr)
library(dplyr)

#data<-read.csv("Fer Expression_CRC cell line with KRAS status and drug AUC.csv",header = T, stringsAsFactors = T)
data<-read.csv("Fer Expression_CRC cell line with KRAS status and drug AUC_KRAS_status_renamed_gender_primarytumor.csv",header = T, stringsAsFactors = T)

#data<-filter(data,data$age>="55",stage <="3")
#data<-data %>% replace_with_na_all(~.x == "N/A")## replace N/A to NA
data<-data %>% drop_na(Sex)

data<-filter(data, Sex=="Female")


data[is.na(data)] = 0

Wdata<-filter(data, Kras=="WT")
Mdata<-filter(data, Kras=="MT")


##Spearman


res1 <- lapply(seq(394, 888), function(i){
  cat(i); cat("\n")
  res1 <- cor.test(as.numeric(Wdata$MBOAT2), as.numeric(Wdata[,i]))
  result <- res1$estimate
  return(result)
}) 



res1<- do.call(rbind, res1)
res1



res2 <- lapply(seq(394, 888), function(i){
  cat(i); cat("\n")
  res2 <- cor.test(as.numeric(Mdata$MBOAT2), as.numeric(Mdata[,i]))
  result <- res2$estimate
  return(result)
})

res2<- do.call(rbind, res2)
res2



res3 <- lapply(seq(394, 888), function(i){
  cat(i); cat("\n")
  res3 <- cor.test(as.numeric(Wdata$MBOAT2), as.numeric(Wdata[,i]))
  result <- res3$p.value
  return(result)
})

res3<- do.call(rbind, res3)
res3



res4 <- lapply(seq(394, 888), function(i){
  cat(i); cat("\n")
  res4 <- cor.test(as.numeric(Mdata$MBOAT2), as.numeric(Mdata[,i]))
  result <- res4$p.value
  return(result)
})

res4<- do.call(rbind, res4)
res4



final<-cbind(res1,res2,res3,res4)

colnames(final)[1] <- 'KRAS WT cor'
colnames(final)[2] <- 'KRAS MT cor'
colnames(final)[3] <- 'KRAS WT p value'
colnames(final)[4] <- 'KRAS MT p value'


write.csv(final,"Female spearman MBOAT2.csv")

    



####Plot-------------------------------------------------------------------------------------------------------

sp <- ggscatter(mMdata, x = "ACSL3", y = "os.time",
                add = "reg.line",  # Add regressin line
                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE # Add confidence interval
)
# Add correlation coefficient
sp + stat_cor(method = "spearman", label.x = 3, label.y = 12)



### import data for spearman plot
data<-read.csv("Fer Expression_CRC cell line with KRAS status and drug AUC_KRAS_status_renamed.csv",header = T, stringsAsFactors = T)

#data<-filter(data,data$age>="55",stage <="3")
#data<-data %>% replace_with_na_all(~.x == "N/A")## replace N/A to NA
data<-data %>% drop_na(KRAS)

data[is.na(data)] = 0

my_data<-data
##by group

sp <- ggscatter(my_data, x = "ACSL3", y = "QL.XII.61..GDSC1.1203.",
                color = "KRAS", palette = "uchicago", shape = "KRAS" ,
                add = "reg.line", conf.int = TRUE,xlab ="ACSL3(Log2)", ylab = "QL.XII.61 AUC")
sp<-sp+ stat_cor(aes(color = KRAS), method="spearman",position = "jitter")


sp1 <- ggscatter(my_data, x = "ACSL3", y = "temozolomide..GDSC1.1375.",
                color = "KRAS", palette = "uchicago", shape = "KRAS" ,
                add = "reg.line", conf.int = TRUE, xlab ="ACSL3(Log2)", ylab = "Temozolomide AUC")
sp1<-sp1+ stat_cor(aes(color = KRAS),method="spearman",position = "jitter" )

sp2 <- ggscatter(my_data, x = "ACSL3", y = "pevonedistat..GDSC1.1529.",
                color = "KRAS", palette = "uchicago", shape = "KRAS" ,
                add = "reg.line", conf.int = TRUE,xlab ="ACSL3(Log2)", ylab = "Pevonedistat AUC")
sp2<-sp2+ stat_cor(aes(color = KRAS), method="spearman",position = "jitter")

p<-ggarrange(sp, sp1, sp2 + rremove("x.text"), 
          labels = c("A"),
          ncol = 1, nrow = 3)

ggsave(
  filename = ("Spearman ACSL3.jpg"),
  plot = p,
  width = 8,
  height = 21,
  dpi=600,
  device = "jpg"
)
## ACSL3

sp <- ggscatter(my_data, x = "ACSL3", y = "WY.090217..GDSC1.3.",
                color = "KRAS", palette = "uchicago", shape = "KRAS" ,
                add = "reg.line", conf.int = TRUE,xlab ="ACSL3(Log2)", ylab = "WY090217 AUC")
sp<-sp+ stat_cor(aes(color = KRAS), method="spearman",position = "jitter")


sp1 <- ggscatter(my_data, x = "ACSL3", y = "imatinib..GDSC1.34.",
                 color = "KRAS", palette = "uchicago", shape = "KRAS" ,
                 add = "reg.line", conf.int = TRUE, xlab ="ACSL3(Log2)", ylab = "Imatinib AUC")
sp1<-sp1+ stat_cor(aes(color = KRAS),method="spearman",position = "jitter" )


p<-ggarrange(sp, sp1,
             labels = c("A"),
             ncol = 1, nrow = 2)

ggsave(
  filename = ("Spearman ACSL3.pdf"),
  plot = p,
  width = 8,
  height = 14,
  dpi=600,
  device = "pdf"
)


###ACSL3



## ACSL3


sp <- ggscatter(my_data, x = "ACSL3", y = "ruxolitinib..GDSC1.206.",
                color = "KRAS", palette = "uchicago", shape = "KRAS" ,
                add = "reg.line", conf.int = TRUE,xlab ="ACSL3(Log2)", ylab = "QL.XII.61 AUC")
sp<-sp+ stat_cor(aes(color = KRAS), method="spearman",position = "jitter")


sp1 <- ggscatter(my_data, x = "ACSL3", y = "temozolomide..GDSC1.1375.",
                 color = "KRAS", palette = "uchicago", shape = "KRAS" ,
                 add = "reg.line", conf.int = TRUE, xlab ="ACSL3(Log2)", ylab = "Temozolomide AUC")
sp1<-sp1+ stat_cor(aes(color = KRAS),method="spearman",position = "jitter" )

sp2 <- ggscatter(my_data, x = "ACSL3", y = "pevonedistat..GDSC1.1529.",
                 color = "KRAS", palette = "uchicago", shape = "KRAS" ,
                 add = "reg.line", conf.int = TRUE,xlab ="ACSL3(Log2)", ylab = "Pevonedistat AUC")
sp2<-sp2+ stat_cor(aes(color = KRAS), method="spearman",position = "jitter")

p<-ggarrange(sp, sp1, sp2 + rremove("x.text"), 
             labels = c("A"),
             ncol = 1, nrow = 3)

ggsave(
  filename = ("Spearman ACSL3.pdf"),
  plot = p,
  width = 8,
  height = 20,
  dpi=600,
  device = "pdf"
)

######---------------------------data summary gene cor and metabolites cor----------------------
setwd("C:/Users/hy362/Desktop/CRC database info/CCLE/gene cor")

library(ggpubr)
library(tidyr)
library(dplyr)

#data<-read.csv("Fer Expression_CRC cell line with KRAS status and drug AUC.csv",header = T, stringsAsFactors = T)
data<-read.csv("Male spearman ACSL3.csv",header = T, stringsAsFactors = T)

data<-filter(data, KRAS.MT.p.value<"0.05")


data1<-read.csv("Male spearman ACSL4.csv",header = T, stringsAsFactors = T)

data1<-filter(data1, KRAS.MT.p.value<"0.05")

data2<-read.csv("Male spearman AIFM2.csv",header = T, stringsAsFactors = T)

data2<-filter(data2, KRAS.MT.p.value<"0.05")

data3<-read.csv("Male spearman ALOX12.csv",header = T, stringsAsFactors = T)

data3<-filter(data3, KRAS.MT.p.value<"0.05")

data4<-read.csv("Male spearman ALOX15.csv",header = T, stringsAsFactors = T)

data4<-filter(data4, KRAS.MT.p.value<"0.05")

data5<-read.csv("Male spearman ALOX15B.csv",header = T, stringsAsFactors = T)

data5<-filter(data5, KRAS.MT.p.value<"0.05")

data6<-read.csv("Male spearman ATF4.csv",header = T, stringsAsFactors = T)

data6<-filter(data6, KRAS.MT.p.value<"0.05")

data7<-read.csv("Male spearman CBS.csv",header = T, stringsAsFactors = T)

data7<-filter(data7, KRAS.MT.p.value<"0.05")

data8<-read.csv("Male spearman DHODH.csv",header = T, stringsAsFactors = T)

data8<-filter(data8, KRAS.MT.p.value<"0.05")

data9<-read.csv("Male spearman ELOVL5.csv",header = T, stringsAsFactors = T)

data9<-filter(data9, KRAS.MT.p.value<"0.05")

data10<-read.csv("Male spearman FH.csv",header = T, stringsAsFactors = T)

data10<-filter(data10, KRAS.MT.p.value<"0.05")

data11<-read.csv("Male spearman FTL.csv",header = T, stringsAsFactors = T)

data11<-filter(data11, KRAS.MT.p.value<"0.05")

data12<-read.csv("Male spearman GCH1.csv",header = T, stringsAsFactors = T)

data12<-filter(data12, KRAS.MT.p.value<"0.05")

data13<-read.csv("Male spearman GCLC.csv",header = T, stringsAsFactors = T)

data13<-filter(data13, KRAS.MT.p.value<"0.05")

data14<-read.csv("Male spearman GLS2.csv",header = T, stringsAsFactors = T)

data14<-filter(data14, KRAS.MT.p.value<"0.05")

data15<-read.csv("Male spearman HMOX1.csv",header = T, stringsAsFactors = T)

data15<-filter(data15, KRAS.MT.p.value<"0.05")

data16<-read.csv("Male spearman LPCAT3.csv",header = T, stringsAsFactors = T)

data16<-filter(data16, KRAS.MT.p.value<"0.05")

data17<-read.csv("Male spearman NFE2L2.csv",header = T, stringsAsFactors = T)

data17<-filter(data17, KRAS.MT.p.value<"0.05")



data,data1,data2,data3,data4,data5,data6,data7,data8,data9,data10,data11,data12,data13,data14,data15,data16,data17
