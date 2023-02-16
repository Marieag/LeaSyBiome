#Code written by M. R. Aggerbeck. 
#mrag@envs.au.dk
#Creative commons, 2023. 

###----------------PACKAGES--------------
#Packages used in the following code. 
library("phyloseq"); packageVersion("phyloseq")
library("devtools"); packageVersion("devtools")
library("biomformat"); packageVersion("biomformat")
library("vegan"); packageVersion("vegan")
library("ggplot2"); packageVersion("ggplot2")
library("DESeq2"); packageVersion("DESeq2")
library("mvabund"); packageVersion("mvabund")
library("metacoder"); packageVersion("metacoder")
library("taxa"); packageVersion("taxa")
library("microbiomeSeq"); packageVersion("microbiomeSeq")
library("adespatial"); packageVersion("adespatial")
library("ggpubr"); packageVersion("ggpubr")
library("data.table"); packageVersion("data.table")
library("igraph"); packageVersion("igraph")
library("tidyverse"); packageVersion("tidyverse")
library("plotrix"); packageVersion("plotrix")
library("microbiome"); packageVersion("microbiome")
library("ggplotify"); packageVersion("ggplotify")
library("ggordiplots"); packageVersion("ggordiplots")
library("pairwiseAdonis"); packageVersion("pairwiseAdonis")


### ----- LOAD DATA - Hypothesis 1 ------------

otu_tbl3 <- read.csv('Study_2/rarefied-table.tsv', header = T, sep = "\t", quote = "\"", dec = ".", fill = TRUE, comment.char = "")
tax_tbl3 <- read.csv('Study_2/taxonomy.tsv', header = TRUE, sep = "\t", quote = "\"", dec = ".", fill = TRUE, comment.char = "")
meta3 <- read.csv('Study_2/Metadata_Template_GDF_tentativo3.txt', header = TRUE, sep = "\t", quote = "\"", dec = ".", fill = TRUE, comment.char = "")

head(otu_tbl3)
head(tax_tbl3)
head(meta3)

class(otu_tbl3)
row.names(otu_tbl3) <- otu_tbl3$OTU.ID
otu_tbl3[1] <- NULL
temp <- colnames(otu_tbl3) 
head(temp)
# temp <- gsub("M.", "M_", temp)
# colnames(otu_tbl3) <- temp
head(otu_tbl3)

OTU3 = otu_table(otu_tbl3, taxa_are_rows = TRUE)



# Rumakantas kode 
# tax <- read.delim("./taxonomy.tsv", as.is = TRUE) %>% select(-Consensus) %>%  mutate(Taxon = str_remove_all(Taxon, ".__| ")) %>%
#   separate(Taxon, into = c("Kingdom", "phylum", "class", "order", "family", "genus", "species"), sep = ";")
# rownames(tax) <- tax$Feature.ID
# tax <- tax %>% select(-Feature.ID)%>% as.matrix()

temp <- tax_tbl3 %>% add_column(taxonomy = tax_tbl3$Taxon %>% str_split(';', simplify = T))
head(temp)

# temp <- as.data.frame(tax_table(physeq_st2))


unique(temp$taxonomy[,2])

temp2 <- as.matrix(temp$taxonomy)
head(temp2)
temp2[1,7]
temp2[temp2==""] <- "Unidentified"
temp2[1,7]


rownames(temp2) = tax_tbl3$Feature.ID
colnames(temp2) = c("Kingdom","Phylum","Class","Order","Family","Genus", "Species")
head(temp2)

TAX3 = tax_table(as.matrix(temp2))



temp <- as.matrix(meta3)
temp <- gsub("F-", "F.", temp)
temp <- gsub("S-", "S.", temp)
# colnames(otu_tbl3) <- temp
rownames(temp) = temp[,1]
head(temp)
metadata3 <- temp
metadata3 <- as.data.frame(metadata3)
metadata3 <- sample_data(metadata3)

head(OTU3)
head(TAX3)
head(metadata3)

sample_names(metadata3)


physeq_st2 = phyloseq(OTU3, TAX3)
physeq_st2
physeq_st2 <- merge_phyloseq(physeq_st2,metadata3)
physeq_st2

sample_names(physeq_st2)

backup <- physeq_st2
# physeq_st2 <- backup

### -----------Make unique names for taxa---------------------------

tax_table(physeq_st2)[,7] <-  gsub("Unidentified", "s__unidentified", tax_table(physeq_st2)[,7]); 
tax_table(physeq_st2)[,6] <-  gsub("Unidentified", "g__unidentified", tax_table(physeq_st2)[,6]); 
tax_table(physeq_st2)[,5] <-  gsub("Unidentified", "f__unidentified", tax_table(physeq_st2)[,5]);
tax_table(physeq_st2)[,4] <-  gsub("Unidentified", "o__unidentified", tax_table(physeq_st2)[,4]);
tax_table(physeq_st2)[,3] <-  gsub("Unidentified", "c__unidentified", tax_table(physeq_st2)[,3]);
tax_table(physeq_st2)[,2] <-  gsub("Unidentified", "p__unidentified", tax_table(physeq_st2)[,2]);

sample_data(physeq_st2)[1:10]
otu_table(physeq_st2) [1:10]
tax_table(physeq_st2)[1:10]

###Agglomeration
phy_sp <- physeq_st2
phy_sp <- tax_glom(phy_sp, taxrank="Species"); phy_sp

sample_data(phy_sp)[1:5]
otu_table(phy_sp) [1:5]
tax_table(phy_sp)[1:5]

tax_table(physeq_st2) <-  gsub("s__unidentified", "", tax_table(physeq_st2)); 
tax_table(physeq_st2) <-  gsub("g__unidentified", "", tax_table(physeq_st2)); 
tax_table(physeq_st2) <-  gsub("f__unidentified", "", tax_table(physeq_st2));
tax_table(physeq_st2) <-  gsub("o__unidentified", "", tax_table(physeq_st2));
tax_table(physeq_st2) <-  gsub("c__unidentified", "", tax_table(physeq_st2));
tax_table(physeq_st2) <-  gsub("p__unidentified", "", tax_table(physeq_st2));

mynames = NULL

for (i in 1:length(taxa_names(physeq_st2))){
  name <- makeTaxLabel(taxa_names(physeq_st2)[i],physeq_st2)
  mynames <- rbind(mynames, c(name))
  
}

#Find duplicates
n_occur <- data.frame(table(mynames))
n_occur[n_occur$Freq > 1,]


taxa_names(physeq_st2) = make.unique(mynames, sep="_")
taxa_names(physeq_st2)[1:10]

physeq_st2
sample_data(physeq_st2)[1:5]
otu_table(physeq_st2) [1:5]
tax_table(physeq_st2)[1:15]



###----- B.	Mycobiome alterations due to treatment/cultivar/location ----------



physeq_st2
# otu_table()   OTU Table:         [ 854 taxa and 117 samples ]
# sample_data() Sample Data:       [ 117 samples by 6 sample variables ]
# tax_table()   Taxonomy Table:    [ 854 taxa by 7 taxonomic ranks ]

#Don't use this - needs redoing.
phy_sp
# otu_table()   OTU Table:         [ 187 taxa and 117 samples ]
# sample_data() Sample Data:       [ 117 samples by 6 sample variables ]
# tax_table()   Taxonomy Table:    [ 187 taxa by 7 taxonomic ranks ]

###----- palettes ------

# "BDA", "Esca", "Apox", "Asy"

set20dark <- c("darkorchid4", "darkred", "bisque4", "forestgreen")
set20 <- c("darkorchid1", "red3", "bisque3", "green3")
set21dark <- c("darkorchid4", "darkred", "forestgreen")
set21 <- c("darkorchid1", "red3", "green3")

barplotcol <- c("ivory2", "#d23f56", "#66b747", "#915bc6","#b4b137","#5e79c7", "#db9334","#58acd9","#d6522c","#5bbe7d", "#d2489a","#407f40",
                "#c886ca", "#9b486b", "#4ab29c", "#a45232", "#acaf69", "#db7987", "#957130", "#e1956d", "steelblue1", "gray")

set_objCdark <- c("forestgreen", "palegreen3", "darkred", "salmon3") 
set_objC <- c("green3", "palegreen", "red3", "salmon") 

stop("End of chunk")

###------Year & cultivar subsets----------------

phy <- physeq_st2
sample_data(phy)

phy
unique(sample_data(phy)[,4]) #Location
unique(sample_data(phy)[,5]) #Cultivar
unique(sample_data(phy)[,6]) #Fungicide_treatment

phy_lis = subset_samples(phy, Location=="Lisbon")
phy_lis = prune_taxa(taxa_sums(phy_lis) > 0, phy_lis)
phy_lis

phy_lis_cab = subset_samples(phy_lis, Cultivar=="Cabernet Sauvignon")
phy_lis_cab = prune_taxa(taxa_sums(phy_lis_cab) > 0, phy_lis_cab)
phy_lis_cab

phy_lis_syr = subset_samples(phy_lis, Cultivar=="Syrah")
phy_lis_syr = prune_taxa(taxa_sums(phy_lis_syr) > 0, phy_lis_syr)
phy_lis_syr

phy_tur = subset_samples(phy, Location=="Turin")
phy_tur = prune_taxa(taxa_sums(phy_tur) > 0, phy_tur)
phy_tur

phy_tur_cab = subset_samples(phy_tur, Cultivar=="Cabernet Sauvignon")
phy_tur_cab = prune_taxa(taxa_sums(phy_tur_cab) > 0, phy_tur_cab)
phy_tur_cab

phy_tur_syr = subset_samples(phy_tur, Cultivar=="Syrah")
phy_tur_syr = prune_taxa(taxa_sums(phy_tur_syr) > 0, phy_tur_syr)
phy_tur_syr
stop("End of chunk")

###----- 1.	Alpha diversity boxplot comparing symptom types (Esca, BDA, Apox, Asy) for each cultivar (cabernet sauvignon, touriga nacional) for year 2020. ----- 
#Boxplot comparing symptom types (Esca, BDA, Apox, Asy) for each cultivar for year 2021.

#----Location-----

p1 <- plot_richness(phy_lis, x="Fungicide_treatment", color="Cultivar", measures=c("Simpson", "InvSimpson"))
p1 <- p1 + geom_boxplot(data = p1$data, aes(color = NULL), alpha = 0.1)+
  theme_bw()+
  theme(legend.title = element_blank())+
  theme(legend.position = "none")+
  theme(strip.background = element_rect(fill="white" )); p1

p1 <- plot_richness(phy_tur, x="Fungicide_treatment", color="Cultivar", measures=c("Simpson", "InvSimpson"))
p1 <- p1 + geom_boxplot(data = p1$data, aes(color = NULL), alpha = 0.1)+
  theme_bw()+
  theme(legend.title = element_blank())+
  theme(legend.position = "none")+
  theme(strip.background = element_rect(fill="white" )); p1

#----Fungicide_treatment matrix plot-----

#Lisbon

#Cabernet sauvignon
p1 <- plot_richness(phy_lis_cab, x="Fungicide_treatment", color="Cultivar", measures=c("Simpson", "InvSimpson"))
p1 <- p1 + geom_boxplot(data = p1$data, aes(color = NULL), alpha = 0.1)+
  theme_bw()+
  theme(legend.title = element_blank())+
  theme(legend.position = "none")+
  theme(strip.background = element_rect(fill="white" )); p1

# method = c("richness","simpson", "evenness")
bonf=TRUE
grouping_column ="Fungicide_treatment"
pValueCutoff=0.05
meta_table <- sample_data(phy_lis_cab)

#create a starting matrix from a previous run - check proper columns are chosen 
#get sample id, variable, value (and maybe type) - usually the
  #second-to-last three (plus a column from the metadata set).

richness_anova <- p1$data[,7:9]
richness_anova <- cbind(richness_anova, p1$data$Fungicide_treatment)
colnames(richness_anova) <- c("samples", "measure", "value", "Fungicide_treatment")
richness_anova <- richness_anova[order(richness_anova$Fungicide_treatment),]; head(richness_anova)

#perform anova of diversity measure between groups
anova_res <- perform_anova_mra(richness_anova,meta_table,grouping_column,pValueCutoff)
df_pw <- anova_res$df_pw; df_pw #get pairwise p-values

#Apply Bonferroni correction
if(bonf){  
  
  numberofsites <- length(unique(c(as.vector(levels(df_pw$from)),as.vector(levels(df_pw$to)))))[1]
  bonf.cor <- as.numeric(as.matrix(df_pw$p))* numberofsites
  temp <- as.factor(bonf.cor)
  df_pw$p <- temp
  write.table(df_pw, file='Study_2/richness_lisbon_cabernet.tsv', quote=FALSE, sep='\t')    
}

#Draw the boxplots
p<-ggplot(aes_string(x=grouping_column,y="value",color=grouping_column, fill=grouping_column),data=richness_anova)

p<-p+geom_boxplot()+geom_jitter(position = position_jitter(height = 0, width=0))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  facet_wrap(~measure,scales="free_y",nrow=1)+ylab("Alpha Diversity Measure")+xlab("")+
  theme(strip.background = element_rect(fill = "white"))+xlab("")+ 
  labs(title = "Lisbon", subtitle = "Cabernet")
print(p)

#This loop will generate the lines and significances
for(i in 1:dim(df_pw)[1]){
  p<-p+geom_path(inherit.aes=F,aes(x,y),data = data.frame(
    x = c(which(unique(richness_anova[,grouping_column])==as.character(df_pw[i,"from"])),
          which(unique(richness_anova[,grouping_column])==as.character(df_pw[i,"to"]))), 
    y = c(as.numeric(as.character(df_pw[i,"y"])),
          as.numeric(as.character(df_pw[i,"y"]))), 
    measure=c(as.character(df_pw[i,"measure"]),as.character(df_pw[i,"measure"]))), 
    color="black",
    lineend = "butt",arrow = arrow(angle = 90, ends = "both", length = unit(0.1, "inches")))
  p
  p<-p+geom_text(inherit.aes=F,aes(x=x,y=y,label=label),data=data.frame(
    x=(which(unique(richness_anova[,grouping_column])==as.character(df_pw[i,"from"]))+
         which(unique(richness_anova[,grouping_column])==as.character(df_pw[i,"to"])))/2,
    y=as.numeric(as.character(df_pw[i,"y"])),measure=as.character(df_pw[i,"measure"]),
    label=as.character(cut(as.numeric(as.character(df_pw[i,"p"])),
                           breaks=c(-Inf, 0.001, 0.01, 0.05, Inf),label=c("***", "**", "*", "")))))
}

p1 <- p
# p1 <- p1 + scale_color_manual(values = set20dark)+ 
#   scale_fill_manual(values = set20)
print(p1)
print(df_pw)

ggsave(p1, file="Study_2/richness_lisbon_cabernet.pdf", width = 20, height = 20, units = "cm", useDingbats=F)
###


#Syrah
p2 <- plot_richness(phy_lis_syr, x="Fungicide_treatment", color="Cultivar", measures=c("Simpson", "InvSimpson"))
p2 <- p2 + geom_boxplot(data = p2$data, aes(color = NULL), alpha = 0.1)+
  theme_bw()+
  theme(legend.title = element_blank())+
  theme(legend.position = "none")+
  theme(strip.background = element_rect(fill="white" )); p2

# method = c("richness","simpson", "evenness")
bonf=TRUE
grouping_column ="Fungicide_treatment"
pValueCutoff=0.05
meta_table <- sample_data(phy_lis_syr)

#create a starting matrix from a previous run - check proper columns are chosen 
#get sample id, variable, value (and maybe type) - usually the
#second-to-last three (plus a column from the metadata set).

richness_anova <- p2$data[,7:9]
richness_anova <- cbind(richness_anova, p2$data$Fungicide_treatment)
colnames(richness_anova) <- c("samples", "measure", "value", "Fungicide_treatment")
richness_anova <- richness_anova[order(richness_anova$Fungicide_treatment),]; head(richness_anova)

#perform anova of diversity measure between groups
anova_res <- perform_anova_mra(richness_anova,meta_table,grouping_column,pValueCutoff)
df_pw <- anova_res$df_pw; df_pw #get pairwise p-values

#Apply Bonferroni correction
if(bonf){  
  
  numberofsites <- length(unique(c(as.vector(levels(df_pw$from)),as.vector(levels(df_pw$to)))))[1]
  bonf.cor <- as.numeric(as.matrix(df_pw$p))* numberofsites
  temp <- as.factor(bonf.cor)
  df_pw$p <- temp
  write.table(df_pw, file='Study_2/richness_lisbon_syrah.tsv', quote=FALSE, sep='\t')
}

#Draw the boxplots
p<-ggplot(aes_string(x=grouping_column,y="value",color=grouping_column, fill=grouping_column),data=richness_anova)

p<-p+geom_boxplot()+geom_jitter(position = position_jitter(height = 0, width=0))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  facet_wrap(~measure,scales="free_y",nrow=1)+ylab("Alpha Diversity Measure")+xlab("")+
  theme(strip.background = element_rect(fill = "white"))+xlab("")+ 
  labs(title = "Lisbon", subtitle = "Syrah")
print(p)

#This loop will generate the lines and significances
for(i in 1:dim(df_pw)[1]){
  p<-p+geom_path(inherit.aes=F,aes(x,y),data = data.frame(
    x = c(which(unique(richness_anova[,grouping_column])==as.character(df_pw[i,"from"])),
          which(unique(richness_anova[,grouping_column])==as.character(df_pw[i,"to"]))), 
    y = c(as.numeric(as.character(df_pw[i,"y"])),
          as.numeric(as.character(df_pw[i,"y"]))), 
    measure=c(as.character(df_pw[i,"measure"]),as.character(df_pw[i,"measure"]))), 
    color="black",
    lineend = "butt",arrow = arrow(angle = 90, ends = "both", length = unit(0.1, "inches")))
  p
  p<-p+geom_text(inherit.aes=F,aes(x=x,y=y,label=label),data=data.frame(
    x=(which(unique(richness_anova[,grouping_column])==as.character(df_pw[i,"from"]))+
         which(unique(richness_anova[,grouping_column])==as.character(df_pw[i,"to"])))/2,
    y=as.numeric(as.character(df_pw[i,"y"])),measure=as.character(df_pw[i,"measure"]),
    label=as.character(cut(as.numeric(as.character(df_pw[i,"p"])),
                           breaks=c(-Inf, 0.001, 0.01, 0.05, Inf),label=c("***", "**", "*", "")))))
}

p2 <- p
# p2 <- p2 + scale_color_manual(values = set20dark)+ 
#   scale_fill_manual(values = set20)
print(p2)
print(df_pw)

ggsave(p2, file="Study_2/richness_lisbon_syrah.pdf", width = 20, height = 20, units = "cm", useDingbats=F)

###


#Turin

#Cabernet sauvignon
p3 <- plot_richness(phy_tur_cab, x="Fungicide_treatment", color="Cultivar", measures=c("Simpson", "InvSimpson"))
p3 <- p3 + geom_boxplot(data = p3$data, aes(color = NULL), alpha = 0.1)+
  theme_bw()+
  theme(legend.title = element_blank())+
  theme(legend.position = "none")+
  theme(strip.background = element_rect(fill="white" )); p3

# method = c("richness","simpson", "evenness")
bonf=TRUE
grouping_column ="Fungicide_treatment"
pValueCutoff=0.05
meta_table <- sample_data(phy_tur_cab)

#create a starting matrix from a previous run - check proper columns are chosen 
#get sample id, variable, value (and maybe type) - usually the
#second-to-last three (plus a column from the metadata set).

richness_anova <- p3$data[,7:9]
richness_anova <- cbind(richness_anova, p3$data$Fungicide_treatment)
colnames(richness_anova) <- c("samples", "measure", "value", "Fungicide_treatment")
richness_anova <- richness_anova[order(richness_anova$Fungicide_treatment),]; head(richness_anova)

#perform anova of diversity measure between groups
anova_res <- perform_anova_mra(richness_anova,meta_table,grouping_column,pValueCutoff)
df_pw <- anova_res$df_pw; df_pw #get pairwise p-values

#Apply Bonferroni correction
if(bonf){  
  
  numberofsites <- length(unique(c(as.vector(levels(df_pw$from)),as.vector(levels(df_pw$to)))))[1]
  bonf.cor <- as.numeric(as.matrix(df_pw$p))* numberofsites
  temp <- as.factor(bonf.cor)
  df_pw$p <- temp
  write.table(df_pw, file='Study_2/richness_turin_cabernet.tsv', quote=FALSE, sep='\t')    
}

#Draw the boxplots
p<-ggplot(aes_string(x=grouping_column,y="value",color=grouping_column, fill=grouping_column),data=richness_anova)

p<-p+geom_boxplot()+geom_jitter(position = position_jitter(height = 0, width=0))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  facet_wrap(~measure,scales="free_y",nrow=1)+ylab("Alpha Diversity Measure")+xlab("")+
  theme(strip.background = element_rect(fill = "white"))+xlab("")+ 
  labs(title = "Turin", subtitle = "Cabernet")
print(p)

#This loop will generate the lines and significances
for(i in 1:dim(df_pw)[1]){
  p<-p+geom_path(inherit.aes=F,aes(x,y),data = data.frame(
    x = c(which(unique(richness_anova[,grouping_column])==as.character(df_pw[i,"from"])),
          which(unique(richness_anova[,grouping_column])==as.character(df_pw[i,"to"]))), 
    y = c(as.numeric(as.character(df_pw[i,"y"])),
          as.numeric(as.character(df_pw[i,"y"]))), 
    measure=c(as.character(df_pw[i,"measure"]),as.character(df_pw[i,"measure"]))), 
    color="black",
    lineend = "butt",arrow = arrow(angle = 90, ends = "both", length = unit(0.1, "inches")))
  p
  p<-p+geom_text(inherit.aes=F,aes(x=x,y=y,label=label),data=data.frame(
    x=(which(unique(richness_anova[,grouping_column])==as.character(df_pw[i,"from"]))+
         which(unique(richness_anova[,grouping_column])==as.character(df_pw[i,"to"])))/2,
    y=as.numeric(as.character(df_pw[i,"y"])),measure=as.character(df_pw[i,"measure"]),
    label=as.character(cut(as.numeric(as.character(df_pw[i,"p"])),
                           breaks=c(-Inf, 0.001, 0.01, 0.05, Inf),label=c("***", "**", "*", "")))))
}

p3 <- p
# p3 <- p3 + scale_color_manual(values = set20dark)+ 
#   scale_fill_manual(values = set20)
print(p3)
print(df_pw)

ggsave(p3, file="Study_2/richness_turin_cabernet.pdf", width = 20, height = 20, units = "cm", useDingbats=F)
###


#Syrah
p4 <- plot_richness(phy_tur_syr, x="Fungicide_treatment", color="Cultivar", measures=c("Simpson", "InvSimpson"))
p4 <- p4 + geom_boxplot(data = p4$data, aes(color = NULL), alpha = 0.1)+
  theme_bw()+
  theme(legend.title = element_blank())+
  theme(legend.position = "none")+
  theme(strip.background = element_rect(fill="white" )); p4

# method = c("richness","simpson", "evenness")
bonf=TRUE
grouping_column ="Fungicide_treatment"
pValueCutoff=0.05
meta_table <- sample_data(phy_tur_syr)

#create a starting matrix from a previous run - check proper columns are chosen 
#get sample id, variable, value (and maybe type) - usually the
#second-to-last three (plus a column from the metadata set).

richness_anova <- p4$data[,7:9]
richness_anova <- cbind(richness_anova, p4$data$Fungicide_treatment)
colnames(richness_anova) <- c("samples", "measure", "value", "Fungicide_treatment")
richness_anova <- richness_anova[order(richness_anova$Fungicide_treatment),]; head(richness_anova)

#perform anova of diversity measure between groups
anova_res <- perform_anova_mra(richness_anova,meta_table,grouping_column,pValueCutoff)
df_pw <- anova_res$df_pw; df_pw #get pairwise p-values

#Apply Bonferroni correction
if(bonf){  
  
  numberofsites <- length(unique(c(as.vector(levels(df_pw$from)),as.vector(levels(df_pw$to)))))[1]
  bonf.cor <- as.numeric(as.matrix(df_pw$p))* numberofsites
  temp <- as.factor(bonf.cor)
  df_pw$p <- temp
  write.table(df_pw, file='Study_2/richness_turin_syrah.tsv', quote=FALSE, sep='\t')
}

#Draw the boxplots
p<-ggplot(aes_string(x=grouping_column,y="value",color=grouping_column, fill=grouping_column),data=richness_anova)

p<-p+geom_boxplot()+geom_jitter(position = position_jitter(height = 0, width=0))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  facet_wrap(~measure,scales="free_y",nrow=1)+ylab("Alpha Diversity Measure")+xlab("")+
  theme(strip.background = element_rect(fill = "white"))+xlab("")+ 
  labs(title = "Turin", subtitle = "Syrah")
print(p)

#This loop will generate the lines and significances
for(i in 1:dim(df_pw)[1]){
  p<-p+geom_path(inherit.aes=F,aes(x,y),data = data.frame(
    x = c(which(unique(richness_anova[,grouping_column])==as.character(df_pw[i,"from"])),
          which(unique(richness_anova[,grouping_column])==as.character(df_pw[i,"to"]))), 
    y = c(as.numeric(as.character(df_pw[i,"y"])),
          as.numeric(as.character(df_pw[i,"y"]))), 
    measure=c(as.character(df_pw[i,"measure"]),as.character(df_pw[i,"measure"]))), 
    color="black",
    lineend = "butt",arrow = arrow(angle = 90, ends = "both", length = unit(0.1, "inches")))
  p
  p<-p+geom_text(inherit.aes=F,aes(x=x,y=y,label=label),data=data.frame(
    x=(which(unique(richness_anova[,grouping_column])==as.character(df_pw[i,"from"]))+
         which(unique(richness_anova[,grouping_column])==as.character(df_pw[i,"to"])))/2,
    y=as.numeric(as.character(df_pw[i,"y"])),measure=as.character(df_pw[i,"measure"]),
    label=as.character(cut(as.numeric(as.character(df_pw[i,"p"])),
                           breaks=c(-Inf, 0.001, 0.01, 0.05, Inf),label=c("***", "**", "*", "")))))
}

p4 <- p
# p4 <- p4 + scale_color_manual(values = set20dark)+ 
#   scale_fill_manual(values = set20)
print(p4)
print(df_pw)

ggsave(p4, file="Study_2/richness_turin_syrah.pdf", width = 20, height = 20, units = "cm", useDingbats=F)

###

###------------Arrange richness plots ----------------------

pbar3 <- ggarrange(p1, p2, p3, p4,
                   widths = 1, heights = 1.2,
                   ncol = 2, nrow = 2,
                   common.legend = TRUE, legend = "bottom"); pbar3

pbar3

ggsave(pbar3, file="Study_2/richness_anova_overall.pdf", width = 20, height = 20, units = "cm", useDingbats=F)

###


p2 <- p2 + scale_color_manual(values = set20dark)+ 
  scale_fill_manual(values = set20)
print(p2)
print(df_pw)

ggsave(p2, file="Study_2/richness_location_lisbon.pdf", width = 20, height = 20, units = "cm", useDingbats=F)

p3 <- p3 + scale_color_manual(values = set20dark)+ 
  scale_fill_manual(values = set20)
print(p3)
print(df_pw)

ggsave(p3, file="Study_2/richness_location_turin.pdf", width = 20, height = 20, units = "cm", useDingbats=F)



#---Cultivar-----

p1 <- plot_richness(phy, x="Cultivar", color="Cultivar", measures=c("Simpson", "InvSimpson"))
p1 <- p1 + geom_boxplot(data = p1$data, aes(color = NULL), alpha = 0.1)+
  theme_bw()+
  theme(legend.title = element_blank())+
  theme(legend.position = "none")+
  theme(strip.background = element_rect(fill="white" )); p1

# method = c("richness","simpson", "evenness")
bonf=TRUE
grouping_column ="Cultivar"
pValueCutoff=0.05
meta_table <- sample_data(phy)

#create a starting matrix from a previous run - check proper columns are chosen 
#get sample id, variable, value (and maybe type) - usually the
#second-to-last three (plus a column from the metadata set).

richness_anova <- p1$data[,7:9]
# richness_anova <- cbind(richness_anova, p1$data$Symp_abbr_order)
colnames(richness_anova) <- c("samples", "measure", "value") #, "Symptom")
# richness_anova <- richness_anova[order(richness_anova$Symptom),]
head(richness_anova)

#perform anova of diversity measure between groups
anova_res <- perform_anova_mra(richness_anova,meta_table,grouping_column,pValueCutoff)
df_pw <- anova_res$df_pw #get pairwise p-values
# df_pw
# write.table(df_pw, file='results_october/richness_anova_cab20.tsv', quote=FALSE, sep='\t')

#Apply Bonferroni correction
if(bonf){  
  
  numberofsites <- length(unique(c(as.vector(levels(df_pw$from)),as.vector(levels(df_pw$to)))))[1]
  bonf.cor <- as.numeric(as.matrix(df_pw$p))* numberofsites
  temp <- as.factor(bonf.cor)
  df_pw$p <- temp
  write.table(df_pw, file='results_october/richness_anova_cab20_bonf.tsv', quote=FALSE, sep='\t')    
}

#Draw the boxplots
p<-ggplot(aes_string(x=grouping_column,y="value",color=grouping_column, fill=grouping_column),data=richness_anova)

p<-p+geom_boxplot()+geom_jitter(position = position_jitter(height = 0, width=0))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  facet_wrap(~measure,scales="free_y",nrow=1)+ylab("Alpha Diversity Measure")+xlab("")+
  theme(strip.background = element_rect(fill = "white"))+xlab("")+ 
  labs(title = "title", subtitle = "subtitle")
print(p)

#This loop will generate the lines and significances
for(i in 1:dim(df_pw)[1]){
  p<-p+geom_path(inherit.aes=F,aes(x,y),data = data.frame(
    x = c(which(unique(richness_anova[,grouping_column])==as.character(df_pw[i,"from"])),
          which(unique(richness_anova[,grouping_column])==as.character(df_pw[i,"to"]))), 
    y = c(as.numeric(as.character(df_pw[i,"y"])),
          as.numeric(as.character(df_pw[i,"y"]))), 
    measure=c(as.character(df_pw[i,"measure"]),as.character(df_pw[i,"measure"]))), 
    color="black",
    lineend = "butt",arrow = arrow(angle = 90, ends = "both", length = unit(0.1, "inches")))
  p
  p<-p+geom_text(inherit.aes=F,aes(x=x,y=y,label=label),data=data.frame(
    x=(which(unique(richness_anova[,grouping_column])==as.character(df_pw[i,"from"]))+
         which(unique(richness_anova[,grouping_column])==as.character(df_pw[i,"to"])))/2,
    y=as.numeric(as.character(df_pw[i,"y"])),measure=as.character(df_pw[i,"measure"]),
    label=as.character(cut(as.numeric(as.character(df_pw[i,"p"])),
                           breaks=c(-Inf, 0.001, 0.01, 0.05, Inf),label=c("***", "**", "*", "")))))
}

p1 <- p
p1 <- p1 + scale_color_manual(values = set20dark)+ 
  scale_fill_manual(values = set20)
print(p1)
print(df_pw)

ggsave(p1, file="Study_2/richness_cultivar.pdf", width = 20, height = 20, units = "cm", useDingbats=F)

#---Fungicide_treatment-----




p1 <- plot_richness(phy, x="Fungicide_treatment", color="Fungicide_treatment", measures=c("Simpson", "InvSimpson"))
p1 <- p1 + geom_boxplot(data = p1$data, aes(color = NULL), alpha = 0.1)+
  theme_bw()+
  theme(legend.title = element_blank())+
  theme(legend.position = "none")+
  theme(strip.background = element_rect(fill="white" )); p1

# method = c("richness","simpson", "evenness")
bonf=TRUE
grouping_column ="Fungicide_treatment"
pValueCutoff=0.05
meta_table <- sample_data(phy)

#create a starting matrix from a previous run - check proper columns are chosen 
#get sample id, variable, value (and maybe type) - usually the
#second-to-last three (plus a column from the metadata set).

richness_anova <- p1$data[,7:9]
# richness_anova <- cbind(richness_anova, p1$data$Symp_abbr_order)
colnames(richness_anova) <- c("samples", "measure", "value") #, "Symptom")
# richness_anova <- richness_anova[order(richness_anova$Symptom),]
head(richness_anova)

#perform anova of diversity measure between groups
anova_res <- perform_anova_mra(richness_anova,meta_table,grouping_column,pValueCutoff)
df_pw <- anova_res$df_pw #get pairwise p-values
# df_pw
# write.table(df_pw, file='results_october/richness_anova_cab20.tsv', quote=FALSE, sep='\t')

#Apply Bonferroni correction
if(bonf){  
  
  numberofsites <- length(unique(c(as.vector(levels(df_pw$from)),as.vector(levels(df_pw$to)))))[1]
  bonf.cor <- as.numeric(as.matrix(df_pw$p))* numberofsites
  temp <- as.factor(bonf.cor)
  df_pw$p <- temp
  write.table(df_pw, file='results_october/richness_anova_cab20_bonf.tsv', quote=FALSE, sep='\t')    
}

#Draw the boxplots
p<-ggplot(aes_string(x=grouping_column,y="value",color=grouping_column, fill=grouping_column),data=richness_anova)

p<-p+geom_boxplot()+geom_jitter(position = position_jitter(height = 0, width=0))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  facet_wrap(~measure,scales="free_y",nrow=1)+ylab("Alpha Diversity Measure")+xlab("")+
  theme(strip.background = element_rect(fill = "white"))+xlab("")+ 
  labs(title = "title", subtitle = "subtitle")
print(p)

#This loop will generate the lines and significances
for(i in 1:dim(df_pw)[1]){
  p<-p+geom_path(inherit.aes=F,aes(x,y),data = data.frame(
    x = c(which(unique(richness_anova[,grouping_column])==as.character(df_pw[i,"from"])),
          which(unique(richness_anova[,grouping_column])==as.character(df_pw[i,"to"]))), 
    y = c(as.numeric(as.character(df_pw[i,"y"])),
          as.numeric(as.character(df_pw[i,"y"]))), 
    measure=c(as.character(df_pw[i,"measure"]),as.character(df_pw[i,"measure"]))), 
    color="black",
    lineend = "butt",arrow = arrow(angle = 90, ends = "both", length = unit(0.1, "inches")))
  p
  p<-p+geom_text(inherit.aes=F,aes(x=x,y=y,label=label),data=data.frame(
    x=(which(unique(richness_anova[,grouping_column])==as.character(df_pw[i,"from"]))+
         which(unique(richness_anova[,grouping_column])==as.character(df_pw[i,"to"])))/2,
    y=as.numeric(as.character(df_pw[i,"y"])),measure=as.character(df_pw[i,"measure"]),
    label=as.character(cut(as.numeric(as.character(df_pw[i,"p"])),
                           breaks=c(-Inf, 0.001, 0.01, 0.05, Inf),label=c("***", "**", "*", "")))))
}

p1 <- p
p1 <- p1 + scale_color_manual(values = set20dark)+ 
  scale_fill_manual(values = set20)
print(p1)
print(df_pw)

ggsave(p1, file="Study_2/richness_cultivar.pdf", width = 20, height = 20, units = "cm", useDingbats=F)




###-------------------------2, Beta diversity -----------------------
# (i) cultivars (cab sauv vs syrah),
# (ii) location (Lisbon vs Turin),
# (iii) treatments (control, copper, esquive, tessior, bentonite), 
# for each cultivar (cabernet sauvignon, syrah), for both locations (Lisbon, Turin).

#Cultivars
phy
pseq.rel <- microbiome::transform(phy, "compositional")
otu <- microbiome::abundances(pseq.rel)
meta <- microbiome::meta(pseq.rel)

dist <- vegan::vegdist(t(otu), method="bray")
mod <- vegan::betadisper(dist, meta$Cultivar, type="centroid")
mod

betadisp <- gg_ordiplot(mod, mod$group)
bd_cult <- betadisp$plot +
  labs(title="Cultivars", subtitle = "Betadispersion, Bray-Curtis")+
  theme_bw()+
  theme(legend.title = element_blank())+ 
  # scale_color_manual(values = set20)+ 
  # scale_fill_manual(values = set20)+
  theme(strip.background = element_rect(fill="white" )); bd_cult
bd_cult

# py20 + scale_color_manual(values = set20)+ 
#   scale_fill_manual(values = set20)

ggsave(bd_cult, file="Study_2/betadisp_cultivars.pdf", width = 20, height = 20, units = "cm", useDingbats=F)

#Cultivars

dist <- vegan::vegdist(t(otu), method="bray")
mod <- vegan::betadisper(dist, meta$Location, type="centroid")
mod

betadisp <- gg_ordiplot(mod, mod$group)
bd_loc <- betadisp$plot +
  labs(title="Location", subtitle = "Betadispersion, Bray-Curtis")+
  theme_bw()+
  theme(legend.title = element_blank())+ 
  # scale_color_manual(values = set20)+ 
  # scale_fill_manual(values = set20)+
  theme(strip.background = element_rect(fill="white" )); bd_loc


# py20 + scale_color_manual(values = set20)+ 
#   scale_fill_manual(values = set20)

ggsave(bd_loc, file="Study_2/betadisp_locations.pdf", width = 20, height = 20, units = "cm", useDingbats=F)

#Treatment

dist <- vegan::vegdist(t(otu), method="bray")
mod <- vegan::betadisper(dist, meta$Fungicide_treatment, type="centroid")
mod

betadisp <- gg_ordiplot(mod, mod$group)
bd_tre <- betadisp$plot +
  labs(title="Treatment, overall", subtitle = "Betadispersion, Bray-Curtis")+
  theme_bw()+
  theme(legend.title = element_blank())+ 
  # scale_color_manual(values = set20)+ 
  # scale_fill_manual(values = set20)+
  theme(strip.background = element_rect(fill="white" )); bd_tre


# py20 + scale_color_manual(values = set20)+ 
#   scale_fill_manual(values = set20)

ggsave(bd_tre, file="Study_2/betadisp_treatmentoverall.pdf", width = 20, height = 20, units = "cm", useDingbats=F)

###---------------Betadisp + permanova, treatment / sites + cultivars -------------------------

#Lisbon Cabernet
pseq.rel <- microbiome::transform(phy_lis_cab, "compositional")
otu <- microbiome::abundances(pseq.rel)
meta <- microbiome::meta(pseq.rel)

dist <- vegan::vegdist(t(otu), method="bray")
mod <- vegan::betadisper(dist, meta$Fungicide_treatment, type="centroid"); mod

TukeyHSD(mod)

betadisp <- gg_ordiplot(mod, mod$group)
bd_lc <- betadisp$plot +
  labs(title="Lisbon Cabernet", subtitle = "Betadispersion, Bray-Curtis")+
  theme_bw()+ theme(legend.title = element_blank())+ 
  # scale_color_manual(values = set20)+ 
  # scale_fill_manual(values = set20)+
  theme(strip.background = element_rect(fill="white" )); bd_lc


ggsave(bd_lc, file="Study_2/betadisp_treatment_liscab.pdf", width = 20, height = 20, units = "cm", useDingbats=F)

post_hoc_permanova <- pairwise.adonis(t(otu), meta$Fungicide_treatment, sim.function = "vegdist",
                                      sim.method = "bray", p.adjust.m = "fdr", reduce = NULL,
                                      perm = 999); post_hoc_permanova
write.table(as.data.frame(post_hoc_permanova), file='Study_2/posthoc_treatment_liscab.tsv', quote=FALSE, sep='\t')

#Lisbon Syrah
pseq.rel <- microbiome::transform(phy_lis_syr, "compositional")
otu <- microbiome::abundances(pseq.rel)
meta <- microbiome::meta(pseq.rel)

dist <- vegan::vegdist(t(otu), method="bray")
mod <- vegan::betadisper(dist, meta$Fungicide_treatment, type="centroid"); mod

TukeyHSD(mod)

betadisp <- gg_ordiplot(mod, mod$group)
bd_ls <- betadisp$plot +
  labs(title="Lisbon, Syrah", subtitle = "Betadispersion, Bray-Curtis")+
  theme_bw()+ theme(legend.title = element_blank())+ 
  # scale_color_manual(values = set20)+ 
  # scale_fill_manual(values = set20)+
  theme(strip.background = element_rect(fill="white" )); bd_ls

ggsave(bd_ls, file="Study_2/betadisp_treatment_lissyr.pdf", width = 20, height = 20, units = "cm", useDingbats=F)

post_hoc_permanova <- pairwise.adonis(t(otu), meta$Fungicide_treatment, sim.function = "vegdist",
                                      sim.method = "bray", p.adjust.m = "fdr", reduce = NULL,
                                      perm = 999); post_hoc_permanova
write.table(as.data.frame(post_hoc_permanova), file='Study_2/posthoc_treatment_lissyr.tsv', quote=FALSE, sep='\t')

#Turin Cabernet
pseq.rel <- microbiome::transform(phy_tur_cab, "compositional")
otu <- microbiome::abundances(pseq.rel)
meta <- microbiome::meta(pseq.rel)

dist <- vegan::vegdist(t(otu), method="bray")
mod <- vegan::betadisper(dist, meta$Fungicide_treatment, type="centroid"); mod

TukeyHSD(mod)

betadisp <- gg_ordiplot(mod, mod$group)
bd_tc <- betadisp$plot +
  labs(title="Turin Cabernet", subtitle = "Betadispersion, Bray-Curtis")+
  theme_bw()+ theme(legend.title = element_blank())+ 
  # scale_color_manual(values = set20)+ 
  # scale_fill_manual(values = set20)+
  theme(strip.background = element_rect(fill="white" )); bd_tc

ggsave(bd_tc, file="Study_2/betadisp_treatment_turcab.pdf", width = 20, height = 20, units = "cm", useDingbats=F)

post_hoc_permanova <- pairwise.adonis(t(otu), meta$Fungicide_treatment, sim.function = "vegdist",
                                      sim.method = "bray", p.adjust.m = "fdr", reduce = NULL,
                                      perm = 999); post_hoc_permanova
write.table(as.data.frame(post_hoc_permanova), file='Study_2/posthoc_treatment_turcab.tsv', quote=FALSE, sep='\t')

#Turin Syrah
pseq.rel <- microbiome::transform(phy_tur_syr, "compositional")
otu <- microbiome::abundances(pseq.rel)
meta <- microbiome::meta(pseq.rel)

dist <- vegan::vegdist(t(otu), method="bray")
mod <- vegan::betadisper(dist, meta$Fungicide_treatment, type="centroid"); mod

TukeyHSD(mod)

betadisp <- gg_ordiplot(mod, mod$group)
bd_ts <- betadisp$plot +
  labs(title="Turin, Syrah", subtitle = "Betadispersion, Bray-Curtis")+
  theme_bw()+ theme(legend.title = element_blank())+ 
  # scale_color_manual(values = set20)+ 
  # scale_fill_manual(values = set20)+
  theme(strip.background = element_rect(fill="white" )); bd_ts

ggsave(bd_ts, file="Study_2/betadisp_treatment_tursyr.pdf", width = 20, height = 20, units = "cm", useDingbats=F)

post_hoc_permanova <- pairwise.adonis(t(otu), meta$Fungicide_treatment, sim.function = "vegdist",
                                      sim.method = "bray", p.adjust.m = "fdr", reduce = NULL,
                                      perm = 999); post_hoc_permanova
write.table(as.data.frame(post_hoc_permanova), file='Study_2/posthoc_treatment_tursyr.tsv', quote=FALSE, sep='\t')

###--------------Arrange betadispersion plots--------------------

pbar4 <- ggarrange(bd_lc, bd_ls, bd_tc, bd_ts,
                   widths = 1, heights = 1.2,
                   ncol = 2, nrow = 2,
                   common.legend = TRUE, legend = "bottom")

pbar4

ggsave(pbar4, file="Study_2/betadispersion_combined.pdf", width = 20, height = 20, units = "cm", useDingbats=F)


###-------------------Permanovas-----------------------

phy
pseq.rel <- microbiome::transform(phy, "compositional")
otu <- microbiome::abundances(pseq.rel)
meta <- microbiome::meta(pseq.rel)

permanova <- vegan::adonis(t(otu)~Location+Cultivar+Fungicide_treatment,
                           data = meta, permutations=999, method = "jaccard")
permanova$aov.tab
write.table(as.data.frame(permanova$aov.tab), file='Study_2/permanova_all.tsv', quote=FALSE, sep='\t')

#post hoc cultivars

post_hoc_permanova <- pairwise.adonis(t(otu), meta$Cultivar, sim.function = "vegdist",
                                      sim.method = "jaccard", p.adjust.m = "fdr", reduce = NULL,
                                      perm = 999)
post_hoc_permanova
write.table(as.data.frame(post_hoc_permanova), file='Study_2/posthoc_cultivar.tsv', quote=FALSE, sep='\t')


#post hoc location

post_hoc_permanova <- pairwise.adonis(t(otu), meta$Location, sim.function = "vegdist",
                                      sim.method = "jaccard", p.adjust.m = "fdr", reduce = NULL,
                                      perm = 999)
post_hoc_permanova
write.table(as.data.frame(post_hoc_permanova), file='Study_2/posthoc_location.tsv', quote=FALSE, sep='\t')

#post hoc treatment

post_hoc_permanova <- pairwise.adonis(t(otu), meta$Fungicide_treatment, sim.function = "vegdist",
                                      sim.method = "jaccard", p.adjust.m = "fdr", reduce = NULL,
                                      perm = 999)
post_hoc_permanova
write.table(as.data.frame(post_hoc_permanova), file='Study_2/posthoc_treatment.tsv', quote=FALSE, sep='\t')

#separate sets - double checking p-values from overall permanova + post hoc 

phy_lis_cab
pseq.rel <- microbiome::transform(phy_lis_cab, "compositional")
otu <- microbiome::abundances(pseq.rel)
meta <- microbiome::meta(pseq.rel)
permanova <- vegan::adonis(t(otu)~Fungicide_treatment,
                           data = meta, permutations=999, method = "jaccard")
permanova$aov.tab
post_hoc_permanova <- pairwise.adonis(t(otu), meta$Fungicide_treatment, sim.function = "vegdist",
                                      sim.method = "jaccard", p.adjust.m = "fdr", reduce = NULL,
                                      perm = 999)
post_hoc_permanova

#
phy_lis_syr
pseq.rel <- microbiome::transform(phy_lis_syr, "compositional")
otu <- microbiome::abundances(pseq.rel)
meta <- microbiome::meta(pseq.rel)
permanova <- vegan::adonis(t(otu)~Fungicide_treatment,
                           data = meta, permutations=999, method = "jaccard")
permanova$aov.tab
post_hoc_permanova <- pairwise.adonis(t(otu), meta$Fungicide_treatment, sim.function = "vegdist",
                                      sim.method = "jaccard", p.adjust.m = "fdr", reduce = NULL,
                                      perm = 999)
post_hoc_permanova

#
phy_tur_cab
pseq.rel <- microbiome::transform(phy_tur_cab, "compositional")
otu <- microbiome::abundances(pseq.rel)
meta <- microbiome::meta(pseq.rel)
permanova <- vegan::adonis(t(otu)~Fungicide_treatment,
                           data = meta, permutations=999, method = "jaccard")
permanova$aov.tab
post_hoc_permanova <- pairwise.adonis(t(otu), meta$Fungicide_treatment, sim.function = "vegdist",
                                      sim.method = "jaccard", p.adjust.m = "fdr", reduce = NULL,
                                      perm = 999)
post_hoc_permanova

#
phy_tur_syr
pseq.rel <- microbiome::transform(phy_tur_syr, "compositional")
otu <- microbiome::abundances(pseq.rel)
meta <- microbiome::meta(pseq.rel)
permanova <- vegan::adonis(t(otu)~Fungicide_treatment,
                           data = meta, permutations=999, method = "jaccard")
permanova$aov.tab
post_hoc_permanova <- pairwise.adonis(t(otu), meta$Fungicide_treatment, sim.function = "vegdist",
                                      sim.method = "jaccard", p.adjust.m = "fdr", reduce = NULL,
                                      perm = 999)
post_hoc_permanova


###--------Merge objects for data viz & barplots-------------------------------------

phy
phy_sp

unique(sample_data(phy)[,4]) #Location
unique(sample_data(phy)[,5]) #Cultivar
unique(sample_data(phy_sp)[,6]) #Fungicide_treatment
unique(sample_data(phy_sp)[,7]) #Combined
head(sample_data(phy))

temp <- as.data.frame(sample_data(phy_sp))
head(temp)

temp$Sample_data <- paste(temp$Location,temp$Cultivar,temp$Fungicide_treatment, sep = "_")

temp$Sample_data <- gsub("t S", "tS", temp$Sample_data)
temp$Sample_data <- gsub("r o", "rO", temp$Sample_data)
temp$Sample_data <- gsub("l w", "lW", temp$Sample_data)
temp

temp$Location_Cultivar <- paste(temp$Location,temp$Cultivar, sep = "_")
temp$Location_Treatment <- paste(temp$Location,temp$Fungicide_treatment, sep = "_")
temp$Cultivar_Treatment <- paste(temp$Cultivar,temp$Fungicide_treatment, sep = "_")
  
sample_data(phy_sp)$Sample_data <- temp$Sample_data

sample_data(phy_sp)$Cultivar_Treatment <- temp$Cultivar_Treatment
sample_data(phy_sp)$Location_Treatment <- temp$Location_Treatment
sample_data(phy_sp)$Location_Cultivar <- temp$Location_Cultivar

temp = prune_taxa(taxa_sums(phy_sp) > 0, phy_sp)
temp2 <- as.data.frame(unique(sample_data(temp)[,7]))
head(temp2)
groups_phy <- as.vector(as.character(temp2$Symptom_Year_Cultivar))
sample_data(temp)$group_phy <- get_variable(temp, "Sample_data") %in% groups_phy

merged_phy = merge_samples(temp, "Sample_data")

merged_phy
sample_names(temp)
sample_names(merged_phy)

temp <- as.data.frame(sample_data(merged_phy)) 

head(temp)

temp$Cultivar <- gsub("1", 'Cabernet Sauvignon', temp$Cultivar)
temp$Cultivar <- gsub("2", 'Syrah', temp$Cultivar)

temp$Location  <- gsub("1", 'Lisbon', temp$Location )
temp$Location  <- gsub("2", 'Turin', temp$Location )

temp$Fungicide_treatment <- gsub("1", 'Bentonite', temp$Fungicide_treatment)
temp$Fungicide_treatment <- gsub("2", 'Control water', temp$Fungicide_treatment)
temp$Fungicide_treatment <- gsub("3", 'Copper oxychloride', temp$Fungicide_treatment)
temp$Fungicide_treatment <- gsub("4", 'Esquive', temp$Fungicide_treatment)
temp$Fungicide_treatment <- gsub("5", 'Tessior', temp$Fungicide_treatment)

temp

sample_data(merged_phy) <- temp
sample_data(merged_phy) 

merged_phy
phy_sp
phy
tmp <- sample_data(merged_phy)

class(tmp)
write.table(tmp, file='merged_metadata.tsv', quote=FALSE, sep='\t')


temp <- tax_glom(merged_phy, taxrank = "Genus")
mynames = NULL

for (i in 1:length(taxa_names(temp))){
  name <- makeTaxLabel(taxa_names(temp)[i],temp)
  mynames <- rbind(mynames, c(name))
  
}

#Find duplicates
n_occur <- data.frame(table(mynames))
n_occur[n_occur$Freq > 1,]


taxa_names(temp) = make.unique(mynames, sep="_")
taxa_names(temp)[1:10]

temp
sample_data(temp)[1:5]
otu_table(temp) [1:5]
tax_table(temp)[1:15]

phy_g <- temp

phy_g


phy_lis_m = subset_samples(temp, Location=="Lisbon")
phy_lis_m = prune_taxa(taxa_sums(phy_lis_m) > 0, phy_lis_m)
phy_lis_m

phy_lis_cab_m = subset_samples(phy_lis_m, Cultivar=="Cabernet Sauvignon")
phy_lis_cab_m = prune_taxa(taxa_sums(phy_lis_cab_m) > 0, phy_lis_cab_m)
phy_lis_cab_m

phy_lis_syr_m = subset_samples(phy_lis_m, Cultivar=="Syrah")
phy_lis_syr_m = prune_taxa(taxa_sums(phy_lis_syr_m) > 0, phy_lis_syr_m)
phy_lis_syr_m

phy_tur_m = subset_samples(temp, Location=="Turin")
phy_tur_m = prune_taxa(taxa_sums(phy_tur_m) > 0, phy_tur_m)
phy_tur_m

phy_tur_cab_m = subset_samples(phy_tur_m, Cultivar=="Cabernet Sauvignon")
phy_tur_cab_m = prune_taxa(taxa_sums(phy_tur_cab_m) > 0, phy_tur_cab_m)
phy_tur_cab_m

phy_tur_syr_m = subset_samples(phy_tur_m, Cultivar=="Syrah")
phy_tur_syr_m = prune_taxa(taxa_sums(phy_tur_syr_m) > 0, phy_tur_syr_m)
phy_tur_syr_m
stop("End of chunk")

###---------------Barplots----------------

lcbd.rel <- microbiome::transform(phy_g, "compositional")
p <- plot_taxa_mra(lcbd.rel,grouping_column="Fungicide_treatment",method="ab.sorensen")
tmp <- as.data.frame(p$data$Sample)
head(tmp)
class(tmp)
tmp <- str_split_fixed(tmp$`p$data$Sample`, "_", 2)
p$data$Location <- tmp[,1]
tmp <- str_split_fixed(tmp[,2], "_", 2)
p$data$Cultivar <- tmp[,1]

head(p$data)


bp <- p + 
  labs(title="", subtitle = "")+
  theme(legend.title = element_blank())+
  # scale_color_manual(values = barplotcol)+
  scale_fill_manual(values = barplotcol)+
  scale_x_discrete(labels= c("Bentonite","Control water", "Copper oxychloride", "Esquive", "Tessior"))
# xlab(substitute())

bp_factgrid1 <- bp + facet_grid(Cultivar~Location, shrink = TRUE, drop=TRUE, scale="free")+
  scale_x_discrete(labels= c("Bentonite","Control water", "Copper oxychloride", "Esquive", "Tessior")); bp_factgrid1
bp_factgrid2 <- bp + facet_grid(~Cultivar+Location, shrink = TRUE, drop=TRUE, scale="free")+
  labs(y= "Relative Abundance", x = "")+
  scale_x_discrete(labels= c("Bentonite","Control water", "Copper oxychloride", "Esquive", "Tessior")); bp_factgrid2

bp_factgrid2
ggsave(bp_factgrid2, file="Study_2/barplot_combined.pdf", width = 20, height = 15, units = "cm", useDingbats=F)


# Lisbon

lcbd.rel <- microbiome::transform(phy_lis_m, "compositional")
p <- plot_taxa_mra(lcbd.rel,grouping_column="Fungicide_treatment",method="ab.sorensen")
bp_lis <- p + 
  labs(title="Lisbon", subtitle = "Cultivars")+
  theme(legend.title = element_blank())+
  # scale_color_manual(values = barplotcol)+
  scale_fill_manual(values = barplotcol)+
  scale_x_discrete(labels= c("Cabernet Sauvignon", "Syrah"))
# xlab(substitute())

print(bp_lis)
ggsave(bp_lis, file="Study_2/barplot_Lisbon_cultivars.pdf", width = 30, height = 20, units = "cm", useDingbats=F)

p <- plot_taxa_mra(lcbd.rel,grouping_column="Cultivar",method="ab.sorensen")
bp_lis <- p + 
  labs(title="Lisbon", subtitle = "Fungicides")+
  theme(legend.title = element_blank())+
  # scale_color_manual(values = barplotcol)+
  scale_fill_manual(values = barplotcol)+
  scale_x_discrete(labels= c("Bentonite","Control water", "Copper oxychloride", "Esquive", "Tessior"))
# xlab(substitute())

print(bp_lis)
ggsave(bp_lis, file="Study_2/barplot_Lisbon_cultivars.pdf", width = 30, height = 20, units = "cm", useDingbats=F)


# cab
lcbd.rel <- microbiome::transform(phy_lis_cab_m, "compositional")
p <- plot_taxa_mra(lcbd.rel,grouping_column="Fungicide_treatment",method="ab.sorensen")
p1 <- p + 
  labs(title="Cabernet sauvignon", subtitle = "Lisbon")+
  theme(legend.title = element_blank(), axis.title.x = element_blank(), axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+ 
  # scale_color_manual(values = barplotcol)+ 
  scale_fill_manual(values = barplotcol)#)+
  # scale_x_discrete(labels= c("Bentonite","Control water", "Copper oxychloride", "Esquive", "Tessior"))

print(p1)
ggsave(p1, file="Study_2/barplot_Lisbon_cabernet.pdf", width = 30, height = 20, units = "cm", useDingbats=F)

# syr
lcbd.rel <- microbiome::transform(phy_lis_syr_m, "compositional")
p <- plot_taxa_mra(lcbd.rel,grouping_column="Fungicide_treatment",method="ab.sorensen")
p2 <- p + 
  labs(title="Syrah", subtitle = "Lisbon")+
  theme(legend.title = element_blank(), axis.title.x = element_blank(), axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+   
  # scale_color_manual(values = barplotcol)+ 
  scale_fill_manual(values = barplotcol)+ 
  scale_x_discrete(labels= c("Bentonite","Control water", "Copper oxychloride", "Esquive", "Tessior"))

print(p2)
ggsave(p2, file="Study_2/barplot_Lisbon_cabernet.pdf", width = 30, height = 20, units = "cm", useDingbats=F)


# Turin
lcbd.rel <- microbiome::transform(phy_tur_m, "compositional")
p <- plot_taxa_mra(lcbd.rel,grouping_column="Cultivar",method="ab.sorensen")
bp_tur <- p + 
  labs(title="Turin", subtitle = "Fungicides")+
  theme(legend.title = element_blank())+ 
  # scale_color_manual(values = barplotcol)+ 
  scale_fill_manual(values = barplotcol)+ 
  scale_x_discrete(labels= c("Bentonite","Control water", "Copper oxychloride", "Esquive", "Tessior"))

print(bp_tur)
ggsave(bp_tur, file="Study_2/barplot_Turin_overall.pdf", width = 30, height = 20, units = "cm", useDingbats=F)

# cab
lcbd.rel <- microbiome::transform(phy_tur_cab_m, "compositional")
p <- plot_taxa_mra(lcbd.rel,grouping_column="Fungicide_treatment",method="ab.sorensen")
p3 <- p + 
  labs(title="Cabernet sauvignon", subtitle = "Turin")+
  theme(legend.title = element_blank(), axis.title.x = element_blank(), axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+ 
  # scale_color_manual(values = barplotcol)+ 
  scale_fill_manual(values = barplotcol)+ 
  scale_x_discrete(labels= c("Bentonite","Control water", "Copper oxychloride", "Esquive", "Tessior"))

print(p3)
ggsave(p3, file="Study_2/barplot_Turin_cabernet.pdf", width = 30, height = 20, units = "cm", useDingbats=F)

# syr
lcbd.rel <- microbiome::transform(phy_tur_syr_m, "compositional")
p <- plot_taxa_mra(lcbd.rel,grouping_column="Fungicide_treatment",method="ab.sorensen")
p4 <- p + 
  labs(title="Syrah", subtitle = "Turin")+
  theme(legend.title = element_blank(), axis.title.x = element_blank(), axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+   
  # scale_color_manual(values = barplotcol)+ 
  scale_fill_manual(values = barplotcol)+ 
  scale_x_discrete(labels= c("Bentonite","Control water", "Copper oxychloride", "Esquive", "Tessior"))

print(p4)
ggsave(p4, file="Study_2/barplot_Turin_Syrah.pdf", width = 30, height = 20, units = "cm", useDingbats=F)



stop("End of chunk")

###------------------ggarrange barplots -----------------------

pbar5 <- ggarrange(p1, p2, p3, p4,
                   widths = 1, heights = 1.5,
                   # labels = c("Cabernet Sauvignon", "Touriga Nacional"),
                   # label.x = 0.8, label.y = .90, hjust = 0, vjust = 0,
                   ncol = 2, nrow = 2,
                   common.legend = FALSE)

pbar5 <- pbar5 + labs(title="20 most common taxa")
pbar5

ggsave(pbar5, file="Study_2/barplot_combined.pdf", width = 30, height = 30, units = "cm", useDingbats=F)

###--------------------------Subset all treatments----------------

phy_sp
taxa_names(phy_sp)
head(sample_data(phy_sp))

minTotRelAbun = 0.1 
x = taxa_sums(phy_sp)
keepTaxa = rownames(as.data.frame(which(((x / sum(x))*100) > minTotRelAbun)))
phy_0_1 = prune_taxa(keepTaxa, phy_sp)
phy_0_1
taxa_names(phy_0_1)
head(sample_data(phy_0_1))


ID<- as.vector(unique(sample_data(phy_0_1)$Sample_data))
ID
for (i in 1:length(ID)){ 
  temp <- subset_samples(phy_0_1, Sample_data==ID[i])
  name <- paste(ID[i])
  assign(name, temp)
}

###--------------------Deseq2 - Make separate objects for pairwise comparison -------------------

Deseq_objlist <- NULL
Deseq_namelist <- NULL

Lisbon_Syrah_ControlWater 
#vs
Deseq2_lis_syr <- c("Lisbon_Syrah_Bentonite","Lisbon_Syrah_Esquive","Lisbon_Syrah_Tessior",
                    "Lisbon_Syrah_CopperOxychloride"); Deseq2_lis_syr

for (i in 1:length(Deseq2_lis_syr)){ 
  temp2 <- subset_samples(phy_sp, Sample_data==Deseq2_lis_syr[i])
  temp = merge_phyloseq(Lisbon_Syrah_ControlWater, temp2)
  name <- paste(Deseq2_lis_syr[i], "vs_Ctrl", sep ="_")
  diagdds = phyloseq_to_deseq2(temp, ~Fungicide_treatment)
  diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
  res = results(diagdds, cooksCutoff = FALSE)
  alpha = 0.05
  sigtab = res[which(res$padj < alpha), ]
  sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(temp)[rownames(sigtab), ], "matrix"))
  head(sigtab)
  assign(name, temp)
  Deseq_namelist <- c(Deseq_namelist, name) 
  y <- as.data.frame(c(name, sigtab)) 
  colnames(y)[1] =  "Pairwise_comparison"
  Deseq_objlist <-  rbind(Deseq_objlist, y)
}

Turin_Syrah_ControlWater
#vs
Deseq2_tur_syr <- c("Turin_Syrah_Bentonite", "Turin_Syrah_Esquive", "Turin_Syrah_Tessior",
                    "Turin_Syrah_CopperOxychloride"); Deseq2_tur_syr

for (i in 1:length(Deseq2_tur_syr)){ 
  temp2 <- subset_samples(phy_sp, Sample_data==Deseq2_tur_syr[i])
  temp = merge_phyloseq(Turin_Syrah_ControlWater, temp2)
  name <- paste(Deseq2_tur_syr[i], "vs_Ctrl", sep ="_")
  diagdds = phyloseq_to_deseq2(temp, ~Fungicide_treatment)
  diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
  res = results(diagdds, cooksCutoff = FALSE)
  alpha = 0.05
  sigtab = res[which(res$padj < alpha), ]
  sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(temp)[rownames(sigtab), ], "matrix"))
  head(sigtab)
  assign(name, temp)
  Deseq_namelist <- c(Deseq_namelist, name) 
  y <- as.data.frame(c(name, sigtab)) 
  colnames(y)[1] =  "Pairwise_comparison"
  Deseq_objlist <-  rbind(Deseq_objlist, y)
}

Lisbon_CabernetSauvignon_ControlWater
#vs
Deseq2_lis_cab <- c("Lisbon_CabernetSauvignon_Bentonite","Lisbon_CabernetSauvignon_Esquive",
                    "Lisbon_CabernetSauvignon_Tessior","Lisbon_CabernetSauvignon_CopperOxychloride"); Deseq2_lis_cab

for (i in 1:length(Deseq2_lis_cab)){ 
  temp2 <- subset_samples(phy_sp, Sample_data==Deseq2_lis_cab[i])
  temp = merge_phyloseq(Lisbon_CabernetSauvignon_ControlWater, temp2)
  name <- paste(Deseq2_lis_cab[i], "vs_Ctrl", sep ="_")
  diagdds = phyloseq_to_deseq2(temp, ~Fungicide_treatment)
  diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
  res = results(diagdds, cooksCutoff = FALSE)
  alpha = 0.05
  sigtab = res[which(res$padj < alpha), ]
  sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(temp)[rownames(sigtab), ], "matrix"))
  head(sigtab)
  assign(name, temp)
  Deseq_namelist <- c(Deseq_namelist, name) 
  y <- as.data.frame(c(name, sigtab)) 
  colnames(y)[1] =  "Pairwise_comparison"
  Deseq_objlist <-  rbind(Deseq_objlist, y)
}

Turin_CabernetSauvignon_ControlWater
#vs
Deseq2_tur_cab <- c("Turin_CabernetSauvignon_Bentonite","Turin_CabernetSauvignon_Esquive",
                    "Turin_CabernetSauvignon_Tessior", "Turin_CabernetSauvignon_CopperOxychloride"); Deseq2_tur_cab

for (i in 1:length(Deseq2_tur_cab)){ 
  temp2 <- subset_samples(phy_sp, Sample_data==Deseq2_tur_cab[i])
  temp = merge_phyloseq(Turin_CabernetSauvignon_ControlWater, temp2)
  name <- paste(Deseq2_tur_cab[i], "vs_Ctrl", sep ="_")
  diagdds = phyloseq_to_deseq2(temp, ~Fungicide_treatment)
  diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
  res = results(diagdds, cooksCutoff = FALSE)
  alpha = 0.05
  sigtab = res[which(res$padj < alpha), ]
  sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(temp)[rownames(sigtab), ], "matrix"))
  head(sigtab)
  assign(name, temp)
  Deseq_namelist <- c(Deseq_namelist, name) 
  y <- as.data.frame(c(name, sigtab)) 
  colnames(y)[1] =  "Pairwise_comparison"
  Deseq_objlist <-  rbind(Deseq_objlist, y)
}
Deseq_namelist
Deseq_objlist

unique(Deseq_objlist$Pairwise_comparison)

write.table(Deseq_objlist, file='Study_2/Deseq_ctrl_vs_treatments_signtaxa_RA0_1.tsv', quote=FALSE, sep='\t')  

# unique(sample_data(Lisbon_Syrah_Bentonite)[,4]) #Location
# unique(sample_data(Lisbon_Syrah_Bentonite)[,5]) #Cultivar
# unique(sample_data(Lisbon_Syrah_Bentonite)[,6]) #Fungicide_treatment
# unique(sample_data(Lisbon_Syrah_Bentonite)[,7]) #Combined
# head(sample_data(Lisbon_Syrah_Bentonite))

# #DESeq2

# diagdds = phyloseq_to_deseq2(temp, ~Fungicide_treatment)
# diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
# print(diagdds)
# res = results(diagdds, cooksCutoff = FALSE)
# alpha = 0.05
# sigtab = res[which(res$padj < alpha), ]
# sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(temp)[rownames(sigtab), ], "matrix"))
# head(sigtab)
# class(sigtab)
# rownames(sigtab)
# Class order
# x = tapply(sigtab$log2FoldChange, sigtab$Class, function(x) max(x))
# x = sort(x, TRUE)
# sigtab$Class = factor(as.character(sigtab$Class), levels=names(x))
# # Subclass order
# x = tapply(sigtab$log2FoldChange, sigtab$Subclass, function(x) max(x))
# x = sort(x, TRUE)
# sigtab$Subclass = factor(as.character(sigtab$Subclass), levels=names(x))
# ggplot(sigtab, aes(x=Class, y=log2FoldChange, color=Genus))+theme_bw() + geom_point(size=2) + 
#   theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

###-------------------Heat trees------------------------------

stars_file = NULL

temp <- phy_lis_cab
minTotRelAbun = 0.1 
x = taxa_sums(temp)
keepTaxa = rownames(as.data.frame(which(((x / sum(x))*100) > minTotRelAbun)))
temp = prune_taxa(keepTaxa, temp)
temp <- tax_glom(temp, taxrank="Genus")
head(tax_table(temp))
temp

#parse phyloseq object temp
obj <- parse_phyloseq(temp)
tissuegroup <- obj$data$sample_data$Fungicide_treatment
obj$data$otu_table <- calc_obs_props(obj, data = "otu_table", cols = obj$data$sample_data$sample_id)
obj$data$tax_abund <- calc_taxon_abund(obj, data = "otu_table", cols = obj$data$sample_data$sample_id)
obj$data$tax_occ <- calc_n_samples(obj, "tax_abund", groups = tissuegroup)
obj$data$diff_table <- compare_groups(obj, data = "tax_abund", 
                                      cols = obj$data$sample_data$sample_id, groups = tissuegroup)
obj$data$diff_table$wilcox_p_value <- p.adjust(obj$data$diff_table$wilcox_p_value, method = "fdr")

obj$data$diff_table

heat_tree_matrix(obj, data = "diff_table",node_size = n_obs,node_label = taxon_names,
                 node_color = log2_median_ratio, node_color_range = diverging_palette(), 
                 node_color_trans = "linear", node_color_interval = c(-3, 3), edge_color_interval = c(-3, 3),
                 layout = "davidson-harel", initial_layout = "reingold-tilford",
                 node_size_axis_label = "Number of compounds",
                 node_color_axis_label = "Log2 ratio median proportions",
                 output_file = 'Study_2/heattree_genus_lis_cab.pdf')

temp2 <- as.data.frame(obj$data$tax_data$otu_id)
colnames(temp2) = "lis_cab"
temp2
stars_file <- c(stars_file, temp2)

###
temp <- phy_tur_cab
minTotRelAbun = 0.1 
x = taxa_sums(temp)
keepTaxa = rownames(as.data.frame(which(((x / sum(x))*100) > minTotRelAbun)))
temp = prune_taxa(keepTaxa, temp)
temp <- tax_glom(temp, taxrank="Genus")
head(tax_table(temp))

#parse phyloseq object temp
obj <- parse_phyloseq(temp)
tissuegroup <- obj$data$sample_data$Fungicide_treatment
obj$data$otu_table <- calc_obs_props(obj, data = "otu_table", cols = obj$data$sample_data$sample_id)
obj$data$tax_abund <- calc_taxon_abund(obj, data = "otu_table", cols = obj$data$sample_data$sample_id)
obj$data$tax_occ <- calc_n_samples(obj, "tax_abund", groups = tissuegroup)
obj$data$diff_table <- compare_groups(obj, data = "tax_abund", 
                                      cols = obj$data$sample_data$sample_id, groups = tissuegroup)
obj$data$diff_table$wilcox_p_value <- p.adjust(obj$data$diff_table$wilcox_p_value, method = "fdr")

heat_tree_matrix(obj, data = "diff_table",node_size = n_obs,node_label = taxon_names,
                 node_color = log2_median_ratio, node_color_range = diverging_palette(), 
                 node_color_trans = "linear", node_color_interval = c(-3, 3), edge_color_interval = c(-3, 3),
                 layout = "davidson-harel", initial_layout = "reingold-tilford",
                 node_size_axis_label = "Number of compounds",
                 node_color_axis_label = "Log2 ratio median proportions",
                 output_file = 'Study_2/heattree_genus_tur_cab.pdf')

temp2 <- as.data.frame(obj$data$tax_data$otu_id)
colnames(temp2) = "tur_cab"
temp2
stars_file <- c(stars_file, temp2)

###
temp <- phy_lis_syr
minTotRelAbun = 0.1 
x = taxa_sums(temp)
keepTaxa = rownames(as.data.frame(which(((x / sum(x))*100) > minTotRelAbun)))
temp = prune_taxa(keepTaxa, temp)
temp <- tax_glom(temp, taxrank="Genus")
head(tax_table(temp))

#parse phyloseq object temp
obj <- parse_phyloseq(temp)
tissuegroup <- obj$data$sample_data$Fungicide_treatment
obj$data$otu_table <- calc_obs_props(obj, data = "otu_table", cols = obj$data$sample_data$sample_id)
obj$data$tax_abund <- calc_taxon_abund(obj, data = "otu_table", cols = obj$data$sample_data$sample_id)
obj$data$tax_occ <- calc_n_samples(obj, "tax_abund", groups = tissuegroup)
obj$data$diff_table <- compare_groups(obj, data = "tax_abund", 
                                      cols = obj$data$sample_data$sample_id, groups = tissuegroup)
obj$data$diff_table$wilcox_p_value <- p.adjust(obj$data$diff_table$wilcox_p_value, method = "fdr")

heat_tree_matrix(obj, data = "diff_table",node_size = n_obs,node_label = taxon_names,
                 node_color = log2_median_ratio, node_color_range = diverging_palette(), 
                 node_color_trans = "linear", node_color_interval = c(-3, 3), edge_color_interval = c(-3, 3),
                 layout = "davidson-harel", initial_layout = "reingold-tilford",
                 node_size_axis_label = "Number of compounds",
                 node_color_axis_label = "Log2 ratio median proportions",
                 output_file = 'Study_2/heattree_genus_lis_syr.pdf')

temp2 <- as.data.frame(obj$data$tax_data$otu_id)
colnames(temp2) = "lis_syr"
temp2
stars_file <- c(stars_file, temp2)

###
temp <- phy_tur_syr

minTotRelAbun = 0.1 
x = taxa_sums(temp)
keepTaxa = rownames(as.data.frame(which(((x / sum(x))*100) > minTotRelAbun)))
temp = prune_taxa(keepTaxa, temp)
temp <- tax_glom(temp, taxrank="Genus")
head(tax_table(temp))

#parse phyloseq object temp
obj <- parse_phyloseq(temp)
tissuegroup <- obj$data$sample_data$Fungicide_treatment
obj$data$otu_table <- calc_obs_props(obj, data = "otu_table", cols = obj$data$sample_data$sample_id)
obj$data$tax_abund <- calc_taxon_abund(obj, data = "otu_table", cols = obj$data$sample_data$sample_id)
obj$data$tax_occ <- calc_n_samples(obj, "tax_abund", groups = tissuegroup)
obj$data$diff_table <- compare_groups(obj, data = "tax_abund", 
                                      cols = obj$data$sample_data$sample_id, groups = tissuegroup)
obj$data$diff_table$wilcox_p_value <- p.adjust(obj$data$diff_table$wilcox_p_value, method = "fdr")

heat_tree_matrix(obj, data = "diff_table",node_size = n_obs,node_label = taxon_names,
                 node_color = log2_median_ratio, node_color_range = diverging_palette(), 
                 node_color_trans = "linear", node_color_interval = c(-3, 3), edge_color_interval = c(-3, 3),
                 layout = "davidson-harel", initial_layout = "reingold-tilford",
                 node_size_axis_label = "Number of compounds",
                 node_color_axis_label = "Log2 ratio median proportions",
                 output_file = 'Study_2/heattree_genus_tur_syr.pdf')

temp2 <- as.data.frame(obj$data$tax_data$otu_id)
colnames(temp2) = "tur_syr"
temp2
stars_file <- c(stars_file, temp2)

class(stars_file)
write.table(stars_file, file='Study_2/stars_file.tsv', quote=FALSE, sep='\t')  
