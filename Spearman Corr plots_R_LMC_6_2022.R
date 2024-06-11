## Partial Correlations EDSS-2year Microbiome ##

setwd("Dropbox (Personal)/Brigham_Womens/Sequencing_Shared/WL-49-MS3/WL-49-MS3-Luke/11-EDSS-Change/")
list.files()

#Read in data for correlations
data_input = read.table("input_L8_EDSS2yKeep_tax_10percent_prevalent_SDMT.txt", header=T, sep="\t", row.names=1,na.strings=c("NA", "-", "?"),stringsAsFactors=F
)


#Format the data
data_input  = data_input[order(rownames(data_input)),]
rownames(data_input)
colnames(data_input)
dim(data_input)

phen = data_input[,1:68]
taxa_L8 = data_input[,69:243]


# Split data by all subjects, Increase, and Stableressive
Increase <- data_input[ which(data_input[,"EDSS.direction"]=='Increase'), ]
Stable <- data_input[ which(data_input[,"EDSS.direction"]=='Stable'), ]
Decrease <- data_input[ which(data_input[,"EDSS.direction"]=='Decrease'),]

phen_Increase = Increase[,1:68]
taxa_L8_Increase = Increase[,69:243]
phen_Stable = Stable[,1:68]
taxa_L8_Stable = Stable[,69:243]
phen_Decrease = Decrease[,1:68]
taxa_L8_Decrease= Decrease[,69:243]

################### DEFINE THE FUNCTION FOR SPEARMAND CORRELATION ################


#Make a correlation test array function
cor.test.array = function(x) {
	a = cor.test(x, analyte, method = "spearman")
	c(a$statistic, a$estimate, p.value = a$p.value)
}

##### Define correlation variable, subset the data, perform the correlation.
cor.var1 = "ChangeEDSS"
dataset1 = "all"
analyte = phen[,cor.var1]
corr_Var1 = apply(taxa_L8, 2, cor.test.array)

cor.var2 = "ChangeEDSS"
dataset2 = "Decrease"
analyte = phen_Decrease[,cor.var2]
corr_Var2 = apply(taxa_L8_Decrease, 2, cor.test.array)

cor.var3 = "ChangeEDSS"
dataset3 = "Stable"
analyte = phen_Stable[,cor.var2]
corr_Var3 = apply(taxa_L8_Stable, 2, cor.test.array)

##Organize the data for output

## Create a matrix of rho values for each taxa (columns) at each timepoint (rows)
rho.corr = rbind(corr_Var1["rho",], corr_Var2["rho",], corr_Var3["rho",])
rownames(rho.corr) = c(dataset1, dataset2,dataset3)
rho.corr = cbind("rho",rho.corr)

## Create a matrix of pvalues values for each taxa (columns) at each timepoint (rows)
p.value.corr = rbind(corr_Var1["p.value",], corr_Var2["p.value",], corr_Var3["p.value",])
rownames(p.value.corr) = c(dataset1, dataset2,dataset3)
p.value.corr = cbind("pval",p.value.corr)

## Adjust pvalues values by fdr method
fdr.corr = rbind(p.adjust(corr_Var1["p.value",], method = "fdr"),p.adjust(corr_Var2["p.value",], method = "fdr"),p.adjust(corr_Var3["p.value",], method = "fdr"))
rownames(fdr.corr) = c(dataset1, dataset2,dataset3)
fdr.corr = cbind("fdr",fdr.corr)


## Merge rho values, p valies,  write the table out.
correlation_result = rbind(rho.corr, p.value.corr,fdr.corr)
write.table(t(correlation_result), file = "correlation_resultEDSStwo_L8-stable_dec.txt", sep ="\t", quote = F, col.names = NA)



################### DEFINE THE FUNCTION FOR PARTIAL CORRELATION ################
##Partial Spearman
install.packages("ppcor")
library(ppcor)

#Make a partial correlation test array function
pcor.test.array = function(x) {
  a = pcor.test(x=x, y=analyte,z=z,method="spearman")
  c(rho=a$estimate,p.value=a$p.value)
}


################### RUN THE PARTIAL CORRELATION ##################
cor.var1 = "ChangeEDSS"
dataset1 = "all"
analyte = phen[cor.var1]
z<-phen[names(phen) %in% c("AGE_AT_SAMPLE")]
z<-z[!is.na(analyte),]
taxa_L8_nomiss<-taxa_L8[!is.na(analyte),]
analyte<-analyte[!is.na(analyte)]
corr_Var1_adj = apply(taxa_L8_nomiss, 2, pcor.test.array)

cor.var2 = "ChangeEDSS"
dataset2 = "Increase"
analyte = phen_Increase[cor.var2]
z<-phen_Increase[names(phen_Increase) %in% c("AGE_AT_SAMPLE")]
z<-z[!is.na(analyte),]
taxa_L8_nomiss<-taxa_L8_Increase[!is.na(analyte),]
analyte<-analyte[!is.na(analyte)]
corr_Var2_adj = apply(taxa_L8_nomiss, 2, pcor.test.array)

cor.var3 = "ChangeEDSS"
dataset3 = "Stable"
analyte = phen_Stable[cor.var3]
z<-phen_Stable[names(phen_Stable) %in% c("AGE_AT_SAMPLE")]
z<-z[!is.na(analyte),]
taxa_L8_nomiss<-taxa_L8_Stable[!is.na(analyte),]
analyte<-analyte[!is.na(analyte)]
corr_Var3_adj = apply(taxa_L8_nomiss, 2, pcor.test.array)



##Organize the data for output

## Create a matrix of rho values for each taxa (columns) at each timepoint (rows)
rho.corr_adj = rbind(corr_Var1_adj["rho",], corr_Var2_adj["rho",], corr_Var3_adj["rho",])
rownames(rho.corr_adj) = c(dataset1, dataset2,dataset3)
rho.corr_adj = cbind("rho_adj",rho.corr_adj)

## Create a matrix of pvalues values for each taxa (columns) at each timepoint (rows)
p.value.corr_adj = rbind(corr_Var1_adj["p.value",], corr_Var2_adj["p.value",], corr_Var3_adj["p.value",])
rownames(p.value.corr_adj) = c(dataset1, dataset2,dataset3)
p.value.corr_adj = cbind("pval_adj",p.value.corr_adj)

## Adjust pvalues values by fdr method
fdr.corr_adj = rbind(p.adjust(corr_Var1_adj["p.value",], method = "fdr"),p.adjust(corr_Var2_adj["p.value",], method = "fdr"),p.adjust(corr_Var3_adj["p.value",], method = "fdr"))
rownames(fdr.corr_adj) = c(dataset1, dataset2,dataset3)
fdr.corr_adj = cbind("fdr_adj",fdr.corr_adj)


## Merge rho values, p valies,  write the table out.
correlation_result_adj = rbind(rho.corr_adj, p.value.corr_adj,fdr.corr_adj)
write.table(t(correlation_result_adj), file = "partial_correlation_result_L8_EDSStwo-AGE-only.txt", sep ="\t", quote = F, col.names = NA)


#######################

### Heatmap

#######################
library(gplots)
library(RColorBrewer)

corr_result_select = read.table('EDSS_change_for_heatmap.txt', header=T, sep="\t", row.names=1)

#corr_result_select = read.table('pval_Increase.txt', header=T, sep="\t", row.names=1)

#corr_result_select = read.table('pval_Stable.txt', header=T, sep="\t", row.names=1)

rho = as.matrix(corr_result_select[,1:2])
pval = as.matrix(corr_result_select[,3:4])

#Convert p-values to asterisks
asterisk.corr = ifelse(pval<=0.001,"***",pval)
asterisk.corr = ifelse(asterisk.corr<=0.01 & asterisk.corr>0.001,"**",asterisk.corr)
asterisk.corr = ifelse(asterisk.corr<=0.05 & asterisk.corr>0.01,"*",asterisk.corr)
asterisk.corr = ifelse(asterisk.corr>0.05," ",asterisk.corr)



my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 255)
my_palette <- colorRampPalette(c("#4477AA", "#77AADD", "#FFFFFF", "#EE9988","#BB4444" ))


pdf("correlation_heatmap_EDSStwo_L8.pdf")        

heatmap.2(rho, cellnote = asterisk.corr, notecol="black",     
main="ChangeEDSS",
col= my_palette ,
scale="none",
Rowv = FALSE,
Colv = FALSE,
dendrogram="none",
key=TRUE,
symkey=TRUE,
density.info="none",
		  margins =c(10,25),     # widens margins around plot
trace="none",cexCol = .5,cexRow = 0.7,
keysize=1)
dev.off()




