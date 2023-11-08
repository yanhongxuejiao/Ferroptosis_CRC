# Load Data
library(readr)
#import mRNA expression data (N=592)
df <- read.delim('C:/RealDocument/Research/ferroptosis/ferr tcgacoad/TCGA COADEAD PanCancer/coadread_tcga_pan_can_atlas_2018/data_mrna_seq_v2_rsem.txt')
clin <- read.delim("C:/RealDocument/Research/ferroptosis/ferr tcgacoad/TCGA COADEAD PanCancer/coadread_tcga_pan_can_atlas_2018/data_clinical_patient.txt")


# Match with FerrDB genes
library (plyr)
ferr.gene <- read.csv(file='C:/RealDocument/Research/ferroptosis/Sur_Ferr_GEO/driver and suppressor.csv')
ferrgene <- ferr.gene$symbol

df1 <-  df[df$Hugo_Symbol %in% ferrgene, ]
df1$Hugo_Symbol
coadread0 <- as.data.frame(t(as.matrix(df1)))
colnames(coadread0) <- coadread0[1,]
coadread<-coadread0[-c(1,2),]

# table(ldply(coadread, function(c) sum(c =="    0.0000"))$V1)# There are many genes with 227 NA's
# table(ldply(coadread, function(c) sum(is.na(c)))$V1) # There are 23 genes with 227 NA's, which should be deleted

## Remove genes with too many NA's. Remaining 381 genes
na.gene <- subset(ldply(coadread, function(c) sum(is.na(c))),V1!=0)$.id
coadread1 <- coadread[, !names(coadread) %in% na.gene]
coadread1<- data.frame(lapply(coadread1, function(x) as.numeric(as.character(x))),
                       check.names=F, row.names = rownames(coadread1))
## Log2 transform
coadread2<-coadread1
ferrgene1 <- colnames(coadread1)
for(i in ferrgene1){
  coadread2[,i]<-ifelse(coadread1[,i]< 0,0,log2(coadread1[,i]+1) )
}

coadread2$PATIENT_ID <- rownames(coadread2)
coadread2<-coadread2[,c(382,1:381)]
coadread2$PATIENT_ID <- gsub("[.]","-",substr(rownames(coadread2),start = 1, stop = 12))

# write.csv(coadread2,file='C:/RealDocument/Research/ferroptosis/ferr tcgacoad/TCGA_COADREAD_mrna_seq_v2_rsem.csv', row.names=F)
# coad_gene <- read.csv(file='C:/RealDocument/Research/ferroptosis/ferr tcgacoad/TCGA_COADREAD_mrna_seq_v2_rsem.csv')


#   Select and recode clinical variables
clin1 <- clin[,c(1,3,5,6,7,20,21,26,27,28,30:37)]
colnames(clin1) <- clin1[4,]
colnames(clin1)[5] <- "STAGE"
clin2 <- clin1[-c(1:4),]

numeric.cols<-c("AGE","WEIGHT","OS_MONTHS","PFS_MONTHS")
clin2[numeric.cols] <- lapply(clin2[numeric.cols], as.numeric)
clin2$STAGE[clin2$STAGE=='STAGE IA']<-'STAGE I'
clin2$STAGE[clin2$STAGE=='STAGE IIA']<-'STAGE II'
clin2$STAGE[clin2$STAGE=='STAGE IIB']<-'STAGE II'
clin2$STAGE[clin2$STAGE=='STAGE IIC']<-'STAGE II'
clin2$STAGE[clin2$STAGE=='STAGE IIIA']<-'STAGE III'
clin2$STAGE[clin2$STAGE=='STAGE IIIB']<-'STAGE III'
clin2$STAGE[clin2$STAGE=='STAGE IIIC']<-'STAGE III'
clin2$STAGE[clin2$STAGE=='STAGE IVA']<-'STAGE IV'
clin2$STAGE[clin2$STAGE=='STAGE IVB']<-'STAGE IV'

factor.cols <- c( "SEX","STAGE", "PATH_M_STAGE","PATH_N_STAGE" ,"RACE" ,"RADIATION_THERAPY","OS_STATUS","PFS_STATUS")
clin2[factor.cols] <- lapply(clin2[factor.cols], as.factor)  


## Plug in KRAS mutation info
mutation <- read.delim('C:/RealDocument/Research/ferroptosis/ferr tcgacoad/TCGA COADEAD PanCancer/coadread_tcga_pan_can_atlas_2018/data_mutations.txt')
mut.pid <- unique(substr(mutation$Tumor_Sample_Barcode,start = 1, stop = 12))
kras.info <- subset(mutation,Hugo_Symbol=="KRAS")
kras.pid <-   unique(substr(kras.info$Tumor_Sample_Barcode,start = 1, stop = 12))
clin2$Mutation.info<-ifelse(clin2$PATIENT_ID %in% mut.pid, "Yes", "No")

library(dplyr)
clin2 <- clin2 %>%
  mutate(KRAS_mut = case_when(
    (clin2$PATIENT_ID %in% kras.pid) ~ "1: Mutant",
    (!clin2$PATIENT_ID %in% kras.pid) & (clin2$Mutation.info=="Yes") ~ "0: WT",
    (!clin2$PATIENT_ID %in% kras.pid) & (clin2$Mutation.info=="No") ~ "Unknown"
  )
  )

coadread.full<- merge(clin2,coad_gene, by = "PATIENT_ID")
table(coadread.full$KRAS_mut)
# write.csv(coadread.full,file='C:/RealDocument/Research/ferroptosis/ferr tcgacoad/TCGA-COADREAD_mRNA_log2.csv', row.names=F)


#   Select only stage I, II, III
tcga123 <- subset(coadread.full, SEX!="" & (STAGE == 'STAGE I'|STAGE == 'STAGE II'|STAGE == 'STAGE III'))
tcga123<-tcga123[,!(names(tcga123) %in% c("MIR9-3HG", "OIP5-AS1", "QSOX1.1" ))]
write.csv(tcga123,file='C:/RealDocument/Research/ferroptosis/ferr tcgacoad/TCGA-COADREAD_123_mRNA_log2.csv', row.names=F)
table(tcga123$KRAS_mut)


#   Demographic table
library(tableone)
shapiro.test(tcga123$AGE)#p-value = 2.152e-06
shapiro.test(tcga123$OS_MONTHS)#p-value < 2.2e-16
shapiro.test(tcga123$PFS_MONTHS)#p-value < 2.2e-16
myVars <- c("AGE","STAGE","PATH_M_STAGE","PATH_N_STAGE" ,"RACE","OS_STATUS","OS_MONTHS","PFS_STATUS","PFS_MONTHS","KRAS_mut")
nonvar <- c("AGE","OS_MONTHS","OS_MONTHS")
factor.cols <- c( "SEX","STAGE", "PATH_M_STAGE","PATH_N_STAGE" ,"RACE","OS_STATUS","PFS_STATUS")
table <- CreateTableOne(vars = myVars, 
                        factorVars = factor.cols,
                        strata = "SEX", 
                        data = tcga123, 
)

table1<- print(table, 
               nonnormal = nonvar,
               catDigits = 2,contDigits = 3,pDigits = 4, 
               showAllLevels=TRUE, 
               quote = FALSE, 
               noSpaces = TRUE,
               printToggle = TRUE) ;table1

write.csv(table1, file = "C:/RealDocument/Research/ferroptosis/ferr tcgacoad/TCGA123_demo.csv")

