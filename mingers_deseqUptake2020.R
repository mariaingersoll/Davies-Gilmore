#set your working directory; if your files are in folders, you have to specify the command to enter the folder with a slash
setwd("~/R/Davies Lab/Orb Fav")
library("DESeq2")
library("ggplot2")

#read in counts using read.table, of pre-filtered count data
  #reads a file in table format and creates a data frame from it (cases corresponding to lines and variables to fields in the file)
countData <- read.table("fav_uptake_cleaned_BM3.txt")
head(countData)
#all data were previously filtered for a base mean of 3
#determine the length of the data file, 
  ####from just the first column with the argument [,1]?
length(countData[,1])

#get and plot total counts data for each reinfection experiment
totalCounts=colSums(countData)
totalCounts
barplot(totalCounts, col=c("green",  "green", "green", "blue" , "blue" , "blue"), ylab="raw counts")
min(totalCounts)
max(totalCounts)

#create the vector "treat" of the strains
treat=c( "B", "B", "B", "D", "D", "D")
#turn the vector "treat" into a data frame called "g" then rename it colData
g=data.frame(treat)
g
colData<- g
colData

#run DESeqDataSetFromMatrix and name it dds with countData for matrix input, colData for columns, and "treat" to dictate how the counts for each gene depend on the variables in colData
#then run DESeq on dds
dds<-DESeqDataSetFromMatrix(countData=countData, colData=colData, design=~treat)
dds<-DESeq(dds)

head(dds)
res<- results(dds)
res
#gives log2 fold change (MLE): treat D vs B
  #####this means that pos log2FC values indicate D is upregulated compared to B
  #####The name provided in the second element (B) is the baseline

#Look at dispersions plot
  #####basically just tells you the spread of your data?
plotDispEsts(dds, main="Dispersion plot Uptake")

#######-----------get rlog data (better transformation when size factors vary across samples)
#regularized log transformation, transforms count data to log2 scale to minimize differences between samples for rows with small counts
rld <- rlogTransformation(dds, blind=TRUE)
head(assay(rld))
hist(assay(rld))
#assay allows you to access the matrix-like data, so that you can view it bc head(rld) gives you just a summary
  ######the histogram displays the frequency of each binned rld count value??

#load the library RColorBrewer
library(RColorBrewer)
#Sample distance heatmap
#as.matrix attempts to turn its argument into a matrix
#dist computes and returns the distance matrix computed using the specified distance measure to compute distances between the rows of a data matrix (assay(rld))
sampleDists <- as.matrix(dist(t(assay(rld))))
  ######what is the "t" for in the above formula?
library(gplots)
#create a sample heatmap of the relatedness based on the count distance/difference between samples
heatmap.2(as.matrix(sampleDists), key=F, trace="none",
          col=colorpanel(100, "black", "white"),
          margin=c(10, 10))
#Now look at results of dds comparing (average) B and D
#the second term is the "control" we will use B as control in this case
resBD <- results(dds, contrast=c("treat","D","B"))
#how many FDR < 10%?
table(resBD$padj<0.1)
# FALSE  TRUE 
# 3543    28 
summary(resBD)
# out of 8573 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 2, 0.023%
# LFC < 0 (down)     : 26, 0.3%
# outliers [1]       : 19, 0.22%
# low counts [2]     : 4983, 58%
# (mean count < 9)

plotMA(resBD, main="B vs D")
  ####positive Log2FC are up in D here

write.table(resBD, file="BD_DE.txt", quote=F, sep="\t")
head(read.delim("BD_DE.txt"))
#idk what this part is but Sarah had it in her script as a note, I ran it just in case
library(tidyverse)
cd %>%
filter(baseMean == 0)
#gave me an error in eval(lhs, parent, parent) : object 'cd' not found

##############################################################################
#get p-values
#create valBD that is a table from resBD of just the pvalues and adjusted p-values (p-adj)
valBD=cbind(resBD$pvalue, resBD$padj)
head(valBD)
#make the table look like you want it
colnames(valBD)=c("pval.BD", "padj.BD")
length(valBD[,1])
#return a logical vector in which both cases (pval and padj) are complete (no missing values)
table(complete.cases(valBD))
#FALSE  TRUE 
#5002  3571 

#make rlogdata (that you can visualize via assay) and pvals table
rlog=rlogTransformation(dds, blind=TRUE) 
rld=assay(rlog)
head(rld)
#make the column names for rld the same ones that were in colData$treat
colnames(rld)=paste(colData$treat)
head(rld)
length(rld[,1])

#combine rld with your table of pval and padj values (valBD)
rldpvals=cbind(rld,valBD)
head(rldpvals)
dim(rldpvals)
#I am getting this output [1] 8573 8
#But Sarah's script says [1] 11977 7
table(complete.cases(rldpvals))
#I get this output:
  #FALSE   TRUE
  #5002    3571
#But Sarah's script says:
  #FALSE   TRUE 
  #8824    3153 

#create a new csv file from rldpvals
write.csv(rldpvals, "timmy_uptake_RLDandPVALS.csv", quote=F)

#########heat map of sample distances based on pvalues
#make a new variable from your csv file
rldpvals <- read.csv(file="timmy_uptake_RLDandPVALS.csv", row.names=1)
head(rldpvals)
#make rld that just returns the first six columns with values (column 1 is the row names)
#Sarah did only the first five columns for some reason but I did all six here
rld=rldpvals[,1:6]
head(rld)

#dist computes and returns the distance matrix computed using the specified distance measure to compute distances between the rows of a data matrix (rld); determine the overall differences in expression between each sample
sampleDists <- dist(t(rld))
#turn the distances computed above into a matrix; five by five matrix
sampleDistMatrix <- as.matrix(sampleDists)

#create a vector of treatment titles for the matrix and apply them to the rows and columns
treat=c("B", "B", "B", "D", "D", "D")
colnames(sampleDistMatrix)=paste(treat)
rownames(sampleDistMatrix)=paste(treat)

install.packages("pheatmap", repos='http://cran.us.r-project.org')
library("pheatmap")
heat.colors = colorRampPalette(rev(c("blue","yellow")),bias=0.3)(100)
quartz()

#apply your matrix of overall expression distance to a heatmap
pheatmap(sampleDistMatrix,color = heat.colors,cex=0.9,border_color=NA,cluster_rows=T,cluster_cols=T)

library(vegan)
library(ggplot2)
library(ggrepel)
library(tidyverse)

#now going to perform PCA on the data to visualize overall effect of experimental covariates and batch effects
  ########t hides all the columns somehow and just shows the gene identifiers??
rld_t=t(rld)
head(rld)
head(rld_t)
#prcomp performs a pca on given data matrix and returns results as an object of class prcomp
pca <- prcomp(rld_t,center = TRUE, scale. = TRUE, na.action=na.omit)
head(pca)
#sdev is the standard deviations of the principal components
#using the sdev, calculate the proportion that each PC corresponds to the variance
li <- pca$sdev^2 / sum(pca$sdev^2)
#round PC1 and PC2 (times 100, to 1 sigfig)
pc1v <- round(li[1] * 100, 1)
pc2v <- round(li[2] * 100, 1)
  #######x from prcomp seems like it's the coordinates of each treatment in the PCA
#turn pca$x into a dataframe
pca_s <- as.data.frame(pca$x)
head(pca_s)
#just take PC1 and PC2 and add your sample and treatment names
pca_s <- pca_s[,c(1,2)]
head(pca_s)
pca_s$Samples = row.names(pca_s)
head(pca_s)
pca_s$treat=colData$treat
head(pca_s)

#creating your PCA plot
cbPalette <- c("darkorchid4","firebrick4")
#aes is aesthetic mappings, how variables are mapped to visual properties
ggplot(pca_s, aes(PC1, PC2, color = treat, pch = treat)) +
  geom_point(size=3) +
  #  geom_text_repel(aes(label=Samples)) +
  scale_colour_manual(values=cbPalette)+
  theme_bw() +
  # geom_density2d(alpha=.5)+
  geom_polygon(alpha=.2)+
  xlab(paste0("PC1: ",pc1v,"% variance")) +
  ylab(paste0("PC2: ",pc2v,"% variance")) 

#adonis is analysis of variance using distance matrices
adonis(pca_s ~ treat, data = pca_s, method='eu', na.rm = TRUE)
#Error in vegdist(lhs, method = method, ...) : input data must be numeric
  #> type(pca_s)
  #[1] "character
#redefined pca_s <- as.data.frame(pca$x)
#pca_s <- pca_s[,c(1,2)]
  #> type(pca_s)
  #[1] "double"
#> adonis(pca_s ~ treat, data = pca_s, method='eu', na.rm = TRUE)
  #'nperm' >= set of all permutations: complete enumeration.
  #Set of permutations < 'minperm'. Generating entire set.

  #Call:
  # adonis(formula = pca_s ~ treat, data = pca_s, method = "eu",      na.rm = TRUE) 

  #Permutation: free
  #Number of permutations: 719

  #Terms added sequentially (first to last)

  #Df SumsOfSqs MeanSqs F.Model     R2 Pr(>F)
  #treat      1    5594.3  5594.3  1.1229 0.2192    0.4
  #Residuals  4   19927.3  4981.8         0.7808       
  #Total      5   25521.6                 1.0000
#Above is different from Sarah's, but I'm going to move on

########plotting genes within each GO category between B and D

#######Manipulated BD_DE.txt in excel following Sarah's verbal instructions
#######Now it is a csv file called BD_DE_GO.csv
#######Following GO_MWU.R script to create three dendrograms now
#######The BD_DE_GO has the -log(p-values)
######From dds (so the timmy_uptake file):
  #gives log2 fold change (MLE): treat D vs B
  #this means that pos log2FC values indicate D is upregulated compared to B
  #The name provided in the second element (B) is the baseline

#Edit these to match your data file names: 
#Running MF first
head(read.csv("BD_DE_GO.csv"))
input="BD_DE_GO.csv" # two columns of comma-separated values: gene id, continuous measure of significance. To perform standard GO enrichment analysis based on Fisher's exact test, use binary measure (0 or 1, i.e., either significant or not).
goAnnotations="orb_fav_iso2go.txt" # two-column, tab-delimited, one line per gene, multiple GO terms separated by semicolon. If you have multiple lines per gene, use nrify_GOtable.pl prior to running this script.
goDatabase="go.obo" # download from http://www.geneontology.org/GO.downloads.ontology.shtml
goDivision="MF" # either MF, or BP, or CC
source("gomwu.functions.R")

# Calculating stats. It might take ~3 min for MF and BP. Do not rerun it if you just want to replot the data with different cutoffs, go straight to gomwuPlot. If you change any of the numeric values below, delete the files that were generated in previous runs first.
gomwuStats(input, goDatabase, goAnnotations, goDivision,
           perlPath="perl", # replace with full path to perl executable if it is not in your system's PATH already
           largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
           smallest=5,   # a GO category should contain at least this many genes to be considered
           clusterCutHeight=0.25) # threshold for merging similar (gene-sharing) terms. See README for details.
#12  GO terms at 10% FDR for MF

# Plotting results
quartz()
resultsMF=gomwuPlot(input,goAnnotations,goDivision,
                  #	absValue=-log(0.05,10),  # genes with the measure value exceeding this will be counted as "good genes". Specify absValue=0.001 if you are doing Fisher's exact test for standard GO enrichment or analyzing a WGCNA module (all non-zero genes = "good genes").
                  absValue=1,
                  level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
                  level2=0.05, # FDR cutoff to print in regular (not italic) font.
                  level3=0.01, # FDR cutoff to print in large bold font.
                  txtsize=1.2,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
                  treeHeight=0.5, # height of the hierarchical clustering tree
                  #	colors=c("dodgerblue2","firebrick1","skyblue","lightcoral") # these are default colors, un-remar and change if needed
)
resultsMF

#Now doing BP
gomwuStats(input, goDatabase, goAnnotations, goDivision,
           perlPath="perl", # replace with full path to perl executable if it is not in your system's PATH already
           largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
           smallest=5,   # a GO category should contain at least this many genes to be considered
           clusterCutHeight=0.25) # threshold for merging similar (gene-sharing) terms. See README for details.
#8  GO terms at 10% FDR for BP

# Plotting results
quartz()
resultsBP=gomwuPlot(input,goAnnotations,goDivision,
                    #	absValue=-log(0.05,10),  # genes with the measure value exceeding this will be counted as "good genes". Specify absValue=0.001 if you are doing Fisher's exact test for standard GO enrichment or analyzing a WGCNA module (all non-zero genes = "good genes").
                    absValue=1,
                    level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
                    level2=0.05, # FDR cutoff to print in regular (not italic) font.
                    level3=0.01, # FDR cutoff to print in large bold font.
                    txtsize=1.2,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
                    treeHeight=0.5, # height of the hierarchical clustering tree
                    #	colors=c("dodgerblue2","firebrick1","skyblue","lightcoral") # these are default colors, un-remar and change if needed
)
resultsBP

#Now doing CC
gomwuStats(input, goDatabase, goAnnotations, goDivision,
           perlPath="perl", # replace with full path to perl executable if it is not in your system's PATH already
           largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
           smallest=5,   # a GO category should contain at least this many genes to be considered
           clusterCutHeight=0.25) # threshold for merging similar (gene-sharing) terms. See README for details.
#19  GO terms at 10% FDR for CC

# Plotting results
quartz()
resultsCC=gomwuPlot(input,goAnnotations,goDivision,
                    #	absValue=-log(0.05,10),  # genes with the measure value exceeding this will be counted as "good genes". Specify absValue=0.001 if you are doing Fisher's exact test for standard GO enrichment or analyzing a WGCNA module (all non-zero genes = "good genes").
                    absValue=1,
                    level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
                    level2=0.05, # FDR cutoff to print in regular (not italic) font.
                    level3=0.01, # FDR cutoff to print in large bold font.
                    txtsize=1.2,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
                    treeHeight=0.5, # height of the hierarchical clustering tree
                    #	colors=c("dodgerblue2","firebrick1","skyblue","lightcoral")
)
resultsCC

####Do I save the results from MF, BP, and CC as MWU_BP_BD_GO.csv or whatever to then match to genes from Timmy Uptake?
####Using Sarah's plot fixer script to change the dimensions of the plots to make them look good
####Just do it in base R

####plotting genes within each GO category between B and D; get the GO identifiers

#First, looking at BP
pv=read.csv("BP_BD_DE_GO.csv", sep=" ")
head(pv)
rownames(pv) <- pv$term
head(pv)
gene=subset(pv,name=="positive regulation of cell proliferation")  #write GO term of interest here
head(gene)
t=rownames(gene)
t
########it keeps saying there are zero observations even though I know that's not true
########Using Dan's code now https://github.com/wuitchik/Divergent-thermal-challenges-elicit-convergent-stress-signatures-in-aposymbiotic-Astrangia-poculata/blob/master/Astrangia_HotCold.pdf

##GO:0008284 is positive regulation of cell proliferation
library(tidyverse)
library(dplyr)

#GO_0008284 
iso2go = read.delim("orb_fav_iso2go.txt")
head(iso2go)
#Get the gene symbols from the gene descriptions
gene = read.delim("orb_fav_iso2gene.tab", sep = "\t")%>%
       mutate(gene_symbol = gsub(".* GN=", "", Description)) %>%
       mutate(gene_symbol = gsub(" .*", "", gene_symbol))
       
head(gene)

#To get my rlog FC expression values, I'm going to manipulate the timmy uptake file
head(read.csv("timmy_uptake_RLDandPVALS.csv"))
BD_DE = read.csv("timmy_uptake_RLDandPVALS.csv")

head(BD_DE)
rlog_BD = BD_DE[,1:7] %>%
  rename(Gene_id = X)
head(rlog_BD)
#rlog_BD just has the iso names called Gene_ID and the log2FC for each sample
#can also use select here instead of [,1:7]

GO_0008284 = iso2go %>%
  filter(str_detect(GO.terms, "GO:0008284")) %>%
  left_join(rlog_BD) %>%
  left_join(gene) %>%
  mutate(gene_symbol = make.names(gene_symbol, unique = TRUE)) %>%
  column_to_rownames(var = "gene_symbol") %>%
  dplyr::select(-GO.terms, -Description, -Gene_id) %>%
  drop_na() %>%
  dplyr::select(sort(current_vars()))

#Yayyyyyyyyy noice lookin good

head(GO_0008284)
colnames(rlog_BD)
head(rlog_BD)
GO_0008284

#get the means of the rows
GO_0008284_means = apply(GO_0008284, 1, mean)
#subtract the means of the rows from the original rows
explc = GO_0008284-GO_0008284_means

library(gplots)
#Make a heatmap
heatmap_GO_0008284 = heatmap.2(as.matrix(GO_0008284), Rowv = TRUE, Colv = FALSE, scale = "row",
          dendrogram = "both",
          trace = "none",
          main = "GO_0008284 positive regulation of cell proliferation",
          margin = c(5,15))
#It worked!

################
#Now looking at oxidoreductase stuff: 
#oxidoreductase complex from CC (GO:1990204); oxidoreductase activity (GO:0016491); oxidoreductase activity, acting on paired donors etc. (GO:0016715)
#to filter for multiple terms, use | (the symbol that is shift back slash) and quotes only outside all the terms
oxidoreductase = iso2go %>%
  filter(str_detect(GO.terms, "GO:1990204|GO:0016491|GO:0016715")) %>%
  left_join(rlog_BD) %>%
  left_join(gene) %>%
  mutate(gene_symbol = make.names(gene_symbol, unique = TRUE)) %>%
  column_to_rownames(var = "gene_symbol") %>%
  dplyr::select(-GO.terms, -Description, -Gene_id) %>%
  drop_na() %>%
  dplyr::select(sort(current_vars()))
head(oxidoreductase)

#get the means of the rows
oxidoreductase_means = apply(oxidoreductase, 1, mean)
#subtract the means of the rows from the original rows
explc = oxidoreductase-oxidoreductase_means

heatmap_oxidoreductase = heatmap.2(as.matrix(oxidoreductase), Rowv = TRUE, Colv = TRUE, scale = "row",
                               dendrogram = "both",
                               trace = "none",
                               main = "GO Oxidoreductase",
                               margin = c(5,15))

#filtering Timmy Uptake data for all genes significant for DE
BD_DE_sig = filter(BD_DE, pval.BD <= 0.05, preserve = TRUE)
head(BD_DE_sig) 
rlog_BD_sig = BD_DE_sig[,1:7] %>%
  rename(Gene_id = X)
head(rlog_BD_sig)


#GO search for just significant cell prolif
GO_0008284_sig = iso2go %>%
  filter(str_detect(GO.terms, "GO:0008284")) %>%
  left_join(rlog_BD_sig) %>%
  left_join(gene) %>%
  mutate(gene_symbol = make.names(gene_symbol, unique = TRUE)) %>%
  column_to_rownames(var = "gene_symbol") %>%
  dplyr::select(-GO.terms, -Description, -Gene_id) %>%
  drop_na() %>%
  dplyr::select(sort(current_vars()))
GO_0008284_sig
#lmao there is only one poor lonely guy at significance of 0.05
  #> GO_0008284_sig
  #B      B.1      B.2        D      D.1      D.2
  #N.2 1.258249 1.758614 2.325618 3.621391 4.958816 3.218951
#At alpha 0.1 there are 6! So we'll use that

GO_0008284_means_sig = apply(GO_0008284_sig, 1, mean)
explc_cellprolif = GO_0008284_sig-GO_0008284_means_sig
head(explc_cellprolif)

col0=colorRampPalette(rev(c("chocolate1","#FEE090","grey10", "cyan3","cyan")))(100)

#heat map of sig cell prolif
heatmap_GO_0008284_sig = heatmap.2(as.matrix(explc_cellprolif), col = col0, Rowv = TRUE, Colv = TRUE, scale = "row",
                               dendrogram = "both",
                               trace = "none",
                               main = "GO_0008284 positive regulation of cell proliferation",
                               margin = c(5,15))

#Now just the significant genes for oxidoreductase
oxidoreductase_sig = iso2go %>%
  filter(str_detect(GO.terms, "GO:1990204|GO:0016491|GO:0016715")) %>%
  left_join(rlog_BD_sig) %>%
  left_join(gene) %>%
  mutate(gene_symbol = make.names(gene_symbol, unique = TRUE)) %>%
  column_to_rownames(var = "gene_symbol") %>%
  dplyr::select(-GO.terms, -Description, -Gene_id) %>%
  drop_na() %>%
  dplyr::select(sort(current_vars()))
oxidoreductase_sig
#At alpha = 0.1 there are only four
oxidoreductase_means_sig = apply(oxidoreductase_sig, 1, mean)
explc_oxr = oxidoreductase_sig-oxidoreductase_means_sig
head(explc_oxr)

heatmap_oxidoreductase_sig = heatmap.2(as.matrix(explc_oxr), col = col0, Rowv = TRUE, Colv = FALSE, scale = "row",
                                   dendrogram = "both",
                                   trace = "none",
                                   main = "GO Oxidoreductase",
                                   margin = c(5,15))

###next step Sarah wants is me to mine for innate immune genes, especially upstream and downstream of NF-kB that are differentially expressed between the two conditions
##Start with the list of all DEGs BD_DE_sig from Timmy Uptake file
#First look at MF (molecular function) for GO terms involved in innate immunity
#Just control f within the MF file will not give me an exhaustive list
#I think what I have to do is make a giant heat map of all the differentially expressed genes and mine from there
#I'll do regular first and then the means

sigDEG_BD = iso2go %>%
  left_join(rlog_BD_sig) %>%
  left_join(gene) %>%
  mutate(gene_symbol = make.names(gene_symbol, unique = TRUE)) %>%
  column_to_rownames(var = "gene_symbol") %>%
  dplyr::select(-GO.terms, -Description, -Gene_id) %>%
  drop_na() %>%
  dplyr::select(sort(current_vars()))
head(sigDEG_BD)
#sigDEG_BD are all my significant DEGs between B and D in the right format
heatmap_sigDEG = heatmap.2(as.matrix(sigDEG_BD), col = col0, Rowv = TRUE, Colv = TRUE, scale = "row",
                           dendrogram = "both",
                           trace = "none",
                           main = "Significant DEG",
                           margin = c(5,15))

#heatmap of the means
sigDEG_BD_means = apply(sigDEG_BD, 1, mean)
explc_sig = sigDEG_BD-sigDEG_BD_means
heatmap_sigDEG_means = heatmap.2(as.matrix(explc_sig), col = col0, Rowv = TRUE, Colv = TRUE, scale = "row",
                           dendrogram = "both",
                           trace = "none",
                           main = "Significant DEG Mean",
                           margin = c(5,15))

sigDEG_BD
explc_sig
nrow(sigDEG_BD)
head(explc_sig)
write.csv(explc_sig, "0.5DEG_BD_genesymbols.csv", quote = F)

#there are 444 DEGs
#How do I determine which ones I want
#Go back to rlog_BD_sig, that has your Gene_id and the rlog for each nubbin
head(rlog_BD_sig)
nrow(rlog_BD_sig)
#look at BD_DE_GO.csv, that has just the gene and the pval
BD_DE_GO = read.csv("BD_DE_GO.csv") %>%
       rename(Gene_id = gene)
#joining stuff
sig_BD_DE_GO = rlog_BD_sig %>%
  left_join(BD_DE_GO) %>%
  dplyr::select(-B, -B.1, -B.2, -D, -D.1, -D.2) %>%
  drop_na() %>%
  dplyr::select(sort(current_vars()))
head(sig_BD_DE_GO)
nrow(sig_BD_DE_GO)
#sig_BD_DE_GO has 691 terms, as does rlog_BD_sig so there must be some non-unique gene_symbol terms in the larger sets
head(sigDEG_BD)

#write sig_BD_DE_GO as a new csv file to input into the GO functions in base R
write.csv(sig_BD_DE_GO, "sig_BD_DE_GO.csv", quote=F)

#now that you have sig_BD_DE_GO.csv with the Gene_id and the pval in the proper format for GO analysis, perform MF, CC, BP
#Doing it in base R!! This may just give you the exact same thing, fair warning
#actually no it shouldn't because now you're only inputting the DEGs and then ranking them, so getting rid of the ones that are ranked as enriched but not actually sig 

input="sig_BD_DE_GO.csv" # two columns of comma-separated values: gene id, continuous measure of significance. To perform standard GO enrichment analysis based on Fisher's exact test, use binary measure (0 or 1, i.e., either significant or not).
goAnnotations="orb_fav_iso2go.txt" # two-column, tab-delimited, one line per gene, multiple GO terms separated by semicolon. If you have multiple lines per gene, use nrify_GOtable.pl prior to running this script.
goDatabase="go.obo" # download from http://www.geneontology.org/GO.downloads.ontology.shtml
goDivision="MF" # either MF, or BP, or CC
source("gomwu.functions.R")
gomwuStats(input, goDatabase, goAnnotations, goDivision,
           perlPath="perl", # replace with full path to perl executable if it is not in your system's PATH already
           largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
           smallest=5,   # a GO category should contain at least this many genes to be considered
           clusterCutHeight=0.25)
head(rld)
head(read.csv("BD_DE_GO.csv"))
#I think what I want to do first is make a more stringent heat map of all DEGs
#Doing an alpha of 0.05

#line 415
#Write a csv file with the Gene_id, all the rlog values, and the pvalues
#Want the -logpval to = 1 if the pval < or = 0.05; all else = 0
#Then two more files: 
  #if pval < or = 0.05 and logfc > 0, -logpval = 1 all else = 0
  #if pval < or = 0.05 and logfc < 0, -logpval = 1 all else = 0
#BD_DE.txt has the log2FC and the pvalue
#First, make a table of just gene_id, log2fc, pvalue
gene_log_pval = read.delim("BD_DE.txt") %>%
  dplyr::select(-baseMean, -lfcSE, -stat, -padj)
head(gene_log_pval)
gene_log_pval = rownames_to_column(gene_log_pval, var = "Gene_id")
head(gene_log_pval)
#then join it to the file with the -logpvalue, which is BD_DE_GO using left_join
BD_DE_GO_id = read.csv("BD_DE_GO.csv") %>%
  rename(Gene_id = gene)
gene_log_pval_logpval = gene_log_pval %>%
  left_join(BD_DE_GO_id) %>%
  rename(neglogp = pval)
head(gene_log_pval_logpval)
#First do all the significant ones regardless of direction: allsig
# <= means less than or equal to
gene_log_pval_logpval$sig = ifelse(gene_log_pval_logpval$pvalue <= 0.05, 1, 0)
head(gene_log_pval_logpval)
allsig = gene_log_pval_logpval %>%
  dplyr::select(-log2FoldChange, -pvalue, -neglogp) %>%
  column_to_rownames(var = "Gene_id")
head(allsig)
write.csv(allsig, "BD_DE_allsig.csv", quote = F)
head(read.csv("BD_DE_allsig.csv"))
#okay I think this is what Sarah wants; it won't give me direction but I'll run it through GOmwu
input="BD_DE_allsig.csv"
goAnnotations="orb_fav_iso2go.txt"
goDatabase="go.obo"
goDivision="CC"
source("gomwu.functions.R")
gomwuStats(input, goDatabase, goAnnotations, goDivision,
           perlPath="perl", # replace with full path to perl executable if it is not in your system's PATH already
           largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
           smallest=5,   # a GO category should contain at least this many genes to be considered
           clusterCutHeight=0.25)
#MF, BP, and CC all gave me 0 GO terms at 10% FDR, make sure that the rownames thing isn't messing it up
#This is just your result, because of the level of stringency
#So basically, what you're ending up with is that heatmap of all the DEGs to explore
