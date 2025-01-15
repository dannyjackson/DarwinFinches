# Gene Lists Analyses

# CRA
cra.dxy.txt
cra.fst.txt

cra.pre.raisd.txt
cra.post.raisd.txt

cra.pre.sweed.txt
cra.post.sweed.txt

cra.pre.taj.txt
cra.post.taj.txt

# FOR
for.dxy.txt
for.fst.txt

for.pre.sweed.txt
for.post.sweed.txt

for.pre.taj.txt
for.post.taj.txt

for.pre.raisd.txt
for.post.raisd.txt

# PAR
par.dxy.txt
par.fst.txt

par.pre.sweed.txt
par.post.sweed.txt

par.pre.taj.txt
par.post.taj.txt

par.pre.raisd.txt
par.post.raisd.txt



import pandas as pd
import numpy as np

pre = pd.read_csv('par.pre.taj.txt')
post = pd.read_csv('par.post.taj.txt')
np.intersect1d(pre, post)

pre_df = np.array(pre)
post_df = np.array(post)

dif1 = np.setdiff1d(pre, post)
dif2 = np.setdiff1d(post, pre)

print(list(dif1))
print(list(dif2))


ANKRD44
FHOD3
LOC115904066
LOC115904067
LOC115904488
MYT1L
PPP1R14D
SUPT3H
UBE2K
ZFYVE19




# local analysis

library(venn)
venn(5, ilab=TRUE, zcolor = "style")

# post-introduction selection
library(ggvenn)

df<-read.csv("finches_genes.csv", header=TRUE)

x <- list(
CRA = df[df$Species == 'CRA',]$Gene,
FOR = df[df$Species == 'FOR',]$Gene,
PAR = df[df$Species == 'PAR',]$Gene
)


pdf(file = "species.pdf", width = 6, height = 6, useDingbats=FALSE)

ggvenn(x,
    columns = c("CRA", "FOR", "PAR"),
    fill_color = c("#4EAFAF", "#FF817E", "#A6C965"),
  stroke_size = 0.5, set_name_size = 4
  )
dev.off()


CRA_df = df[df$Species == 'CRA',]
FOR_df = df[df$Species == 'FOR',]
PAR_df = df[df$Species == 'PAR',]

intersect(CRA_df$Gene, FOR_df$Gene)

DNAH8
GFI1
ARMC1
NKAIN2
KIF6
ZMIZ1

intersect(CRA_df$Gene, PAR_df$Gene)

FSIP1
MYT1L
ZSWIM8
LOC115910218
STARD13

intersect(FOR_df$Gene, PAR_df$Gene)

SUSD4
GOLT1B
LOC115910543
LOC115911100
LOC115917567
VWF

# Total shared list
DNAH8
GFI1
ARMC1
NKAIN2
KIF6
ZMIZ1
FSIP1
MYT1L
ZSWIM8
STARD13
SUSD4
GOLT1B
VWF

LOC115910218
LOC115910543
LOC115911100
LOC115917567

print(list(dif1))
print(list(dif2))

# CRA

x <- list(
raisd = CRA_df[CRA_df$Analysis == 'raisd',]$Gene,
sweed = CRA_df[CRA_df$Analysis == 'sweed',]$Gene,
tajima = CRA_df[CRA_df$Analysis == 'tajima',]$Gene,
dxy = CRA_df[CRA_df$Analysis == 'dxy',]$Gene,
fst = CRA_df[CRA_df$Analysis == 'fst',]$Gene
)

CRA_raisd = CRA_df[CRA_df$Analysis == 'raisd',]
CRA_sweed = CRA_df[CRA_df$Analysis == 'sweed',]
CRA_tajima = CRA_df[CRA_df$Analysis == 'tajima',]
CRA_dxy = CRA_df[CRA_df$Analysis == 'dxy',]
CRA_fst = CRA_df[CRA_df$Analysis == 'fst',]

intersect(CRA_sweed$Gene, CRA_tajima$Gene)
ANKRD44
NEIL3
PACS2
SF3B1
UVRAG
ZNF385D
intersect(CRA_fst$Gene, CRA_dxy$Gene)
COG5
DLG2
FSIP1
GABRB3
LCA5
LOC115902421
LOC115906426
NKAIN2
SLC35F1
SLC9A8
ZNF521


pdf(file = "cra_analyses.pdf", width = 6, height = 6, useDingbats=FALSE)

venn(x, ilab="counts", zcolor = "style", ilcs = 1.5, sncs=1.5)

dev.off()



pdf(file = "cra_twopops.pdf", width = 6, height = 6, useDingbats=FALSE)

ggvenn(x,
    columns = c("dxy", "fst"),
  stroke_size = 0.5, set_name_size = 4
  )
dev.off()

pdf(file = "cra_onepop.pdf", width = 6, height = 6, useDingbats=FALSE)

ggvenn(x,
    columns = c("raisd", "sweed", "tajima"),
  stroke_size = 0.5, set_name_size = 4
  )
dev.off()

CRA_df$Population <-"onepop"

CRA_df[CRA_df$Analysis == "fst", "Population"] <- "twopop"
CRA_df[CRA_df$Analysis == "dxy", "Population"] <- "twopop"

x <- list(
onepop = CRA_df[CRA_df$Population == 'onepop',]$Gene,
twopop = CRA_df[CRA_df$Population == 'twopop',]$Gene
)

pdf(file = "cra_pops.pdf", width = 6, height = 6, useDingbats=FALSE)

ggvenn(x,
    columns = c("onepop", "twopop"),
  stroke_size = 0.5, set_name_size = 4
  )

dev.off()




# FOR

x <- list(
raisd = FOR_df[FOR_df$Analysis == 'raisd',]$Gene,
sweed = FOR_df[FOR_df$Analysis == 'sweed',]$Gene,
tajima = FOR_df[FOR_df$Analysis == 'tajima',]$Gene,
dxy = FOR_df[FOR_df$Analysis == 'dxy',]$Gene,
fst = FOR_df[FOR_df$Analysis == 'fst',]$Gene
)

pdf(file = "for_analyses.pdf", width = 6, height = 6, useDingbats=FALSE)

venn(x, ilab="counts", zcolor = "style", ilcs = 1.5, sncs=1.5)

dev.off()

FOR_raisd = FOR_df[FOR_df$Analysis == 'raisd',]
FOR_sweed = FOR_df[FOR_df$Analysis == 'sweed',]
FOR_tajima = FOR_df[FOR_df$Analysis == 'tajima',]
FOR_dxy = FOR_df[FOR_df$Analysis == 'dxy',]
FOR_fst = FOR_df[FOR_df$Analysis == 'fst',]

intersect(FOR_sweed$Gene, FOR_tajima$Gene)

ADA
ALS2
CELF1
COL12A1
DPF3
EMC7
KATNBL1
MYF6
NEK1
PDE2A
PTPRQ
RAPSN
SATB2
SLC44A5

for_st <- intersect(FOR_sweed$Gene, FOR_tajima$Gene)
intersect(for_st, FOR_fst$Gene)
 "EMC7"    "KATNBL1"

intersect(FOR_tajima$Gene, FOR_fst$Gene)

intersect(FOR_dxy$Gene, FOR_fst$Gene)
CDYL
DAAM2
DYNC1I1
FAM184A
GREM2
HIVEP1
KIF6 
MAL
VSNL1

pdf(file = "for_twopops.pdf", width = 6, height = 6, useDingbats=FALSE)

ggvenn(x,
    columns = c("dxy", "fst"),
  stroke_size = 0.5, set_name_size = 4
  )
dev.off()

pdf(file = "for_onepop.pdf", width = 6, height = 6, useDingbats=FALSE)

ggvenn(x,
    columns = c("raisd", "sweed", "tajima"),
  stroke_size = 0.5, set_name_size = 4
  )
dev.off()

FOR_df$Population <-"onepop"

FOR_df[FOR_df$Analysis == "fst", "Population"] <- "twopop"
FOR_df[FOR_df$Analysis == "dxy", "Population"] <- "twopop"

x <- list(
onepop = FOR_df[FOR_df$Population == 'onepop',]$Gene,
twopop = FOR_df[FOR_df$Population == 'twopop',]$Gene
)

pdf(file = "for_pops.pdf", width = 6, height = 6, useDingbats=FALSE)

ggvenn(x,
    columns = c("onepop", "twopop"),
  stroke_size = 0.5, set_name_size = 4
  )

dev.off()




# PAR

x <- list(
raisd = PAR_df[PAR_df$Analysis == 'raisd',]$Gene,
sweed = PAR_df[PAR_df$Analysis == 'sweed',]$Gene,
tajima = PAR_df[PAR_df$Analysis == 'tajima',]$Gene,
dxy = PAR_df[PAR_df$Analysis == 'dxy',]$Gene,
fst = PAR_df[PAR_df$Analysis == 'fst',]$Gene
)

pdf(file = "par_analyses.pdf", width = 6, height = 6, useDingbats=FALSE)

venn(x, ilab="counts", zcolor = "style", ilcs = 1.5, sncs=1.5)

dev.off()


PAR_raisd = PAR_df[PAR_df$Analysis == 'raisd',]
PAR_sweed = PAR_df[PAR_df$Analysis == 'sweed',]
PAR_tajima = PAR_df[PAR_df$Analysis == 'tajima',]
PAR_dxy = PAR_df[PAR_df$Analysis == 'dxy',]
PAR_fst = PAR_df[PAR_df$Analysis == 'fst',]

intersect(PAR_sweed$Gene, PAR_tajima$Gene)
BDH2
SLC12A7
SLC9B2
SUSD4

intersect(PAR_fst$Gene, PAR_dxy$Gene)
CWC22
DNAJC1
LOC115917567
MLLT10
PDE7B  
SMAD9YPEL5

intersect(PAR_fst$Gene, PAR_tajima$Gene)
STARD13

pdf(file = "par_twopops.pdf", width = 6, height = 6, useDingbats=FALSE)

ggvenn(x,
    columns = c("dxy", "fst"),
  stroke_size = 0.5, set_name_size = 4
  )
dev.off()

pdf(file = "par_onepop.pdf", width = 6, height = 6, useDingbats=FALSE)

ggvenn(x,
    columns = c("raisd", "sweed", "tajima"),
  stroke_size = 0.5, set_name_size = 4
  )
dev.off()

PAR_df$Population <-"onepop"

PAR_df[PAR_df$Analysis == "fst", "Population"] <- "twopop"
PAR_df[PAR_df$Analysis == "dxy", "Population"] <- "twopop"

x <- list(
onepop = PAR_df[PAR_df$Population == 'onepop',]$Gene,
twopop = PAR_df[PAR_df$Population == 'twopop',]$Gene
)

pdf(file = "par_pops.pdf", width = 6, height = 6, useDingbats=FALSE)

ggvenn(x,
    columns = c("onepop", "twopop"),
  stroke_size = 0.5, set_name_size = 4
  )

dev.off()


# consistent selection
library(ggvenn)

df<-read.csv("finches_genes_intersect.csv", header=TRUE)

x <- list(
CRA = df[df$Species == 'CRA',]$Gene,
FOR = df[df$Species == 'FOR',]$Gene,
PAR = df[df$Species == 'PAR',]$Gene
)

pdf(file = "species.pdf", width = 6, height = 6, useDingbats=FALSE)

ggvenn(x,
    columns = c("CRA", "FOR", "PAR"),
    fill_color = c("#4EAFAF", "#FF817E", "#A6C965"),
  stroke_size = 0.5, set_name_size = 4
  )
dev.off()


CRA_df = df[df$Species == 'CRA',]
FOR_df = df[df$Species == 'FOR',]
PAR_df = df[df$Species == 'PAR',]


# CRA
x <- list(
raisd = CRA_df[CRA_df$Analysis == 'raisd',]$Gene,
sweed = CRA_df[CRA_df$Analysis == 'sweed',]$Gene,
tajima = CRA_df[CRA_df$Analysis == 'tajima',]$Gene,
dxy = CRA_df[CRA_df$Analysis == 'dxy',]$Gene,
fst = CRA_df[CRA_df$Analysis == 'fst',]$Gene
)

pdf(file = "cra_twopops.pdf", width = 6, height = 6, useDingbats=FALSE)

ggvenn(x,
    columns = c("dxy", "fst"),
  stroke_size = 0.5, set_name_size = 4
  )
dev.off()

pdf(file = "cra_onepop.pdf", width = 6, height = 6, useDingbats=FALSE)

ggvenn(x,
    columns = c("raisd", "sweed", "tajima"),
  stroke_size = 0.5, set_name_size = 4
  )
dev.off()

CRA_df$Population <-"onepop"

CRA_df[CRA_df$Analysis == "fst", "Population"] <- "twopop"
CRA_df[CRA_df$Analysis == "dxy", "Population"] <- "twopop"

x <- list(
onepop = CRA_df[CRA_df$Population == 'onepop',]$Gene,
twopop = CRA_df[CRA_df$Population == 'twopop',]$Gene
)

pdf(file = "cra_pops.pdf", width = 6, height = 6, useDingbats=FALSE)

ggvenn(x,
    columns = c("onepop", "twopop"),
  stroke_size = 0.5, set_name_size = 4
  )

dev.off()




# FOR

x <- list(
raisd = FOR_df[FOR_df$Analysis == 'raisd',]$Gene,
sweed = FOR_df[FOR_df$Analysis == 'sweed',]$Gene,
tajima = FOR_df[FOR_df$Analysis == 'tajima',]$Gene,
dxy = FOR_df[FOR_df$Analysis == 'dxy',]$Gene,
fst = FOR_df[FOR_df$Analysis == 'fst',]$Gene
)

pdf(file = "for_twopops.pdf", width = 6, height = 6, useDingbats=FALSE)

ggvenn(x,
    columns = c("dxy", "fst"),
  stroke_size = 0.5, set_name_size = 4
  )
dev.off()

pdf(file = "for_onepop.pdf", width = 6, height = 6, useDingbats=FALSE)

ggvenn(x,
    columns = c("raisd", "sweed", "tajima"),
  stroke_size = 0.5, set_name_size = 4
  )
dev.off()

FOR_df$Population <-"onepop"

FOR_df[FOR_df$Analysis == "fst", "Population"] <- "twopop"
FOR_df[FOR_df$Analysis == "dxy", "Population"] <- "twopop"

x <- list(
onepop = FOR_df[FOR_df$Population == 'onepop',]$Gene,
twopop = FOR_df[FOR_df$Population == 'twopop',]$Gene
)

pdf(file = "for_pops.pdf", width = 6, height = 6, useDingbats=FALSE)

ggvenn(x,
    columns = c("onepop", "twopop"),
  stroke_size = 0.5, set_name_size = 4
  )

dev.off()




# PAR

x <- list(
raisd = PAR_df[PAR_df$Analysis == 'raisd',]$Gene,
sweed = PAR_df[PAR_df$Analysis == 'sweed',]$Gene,
tajima = PAR_df[PAR_df$Analysis == 'tajima',]$Gene,
dxy = PAR_df[PAR_df$Analysis == 'dxy',]$Gene,
fst = PAR_df[PAR_df$Analysis == 'fst',]$Gene
)

pdf(file = "par_twopops.pdf", width = 6, height = 6, useDingbats=FALSE)

ggvenn(x,
    columns = c("dxy", "fst"),
  stroke_size = 0.5, set_name_size = 4
  )
dev.off()

pdf(file = "par_onepop.pdf", width = 6, height = 6, useDingbats=FALSE)

ggvenn(x,
    columns = c("raisd", "sweed", "tajima"),
  stroke_size = 0.5, set_name_size = 4
  )
dev.off()

PAR_df$Population <-"onepop"

PAR_df[PAR_df$Analysis == "fst", "Population"] <- "twopop"
PAR_df[PAR_df$Analysis == "dxy", "Population"] <- "twopop"

x <- list(
onepop = PAR_df[PAR_df$Population == 'onepop',]$Gene,
twopop = PAR_df[PAR_df$Population == 'twopop',]$Gene
)

pdf(file = "par_pops.pdf", width = 6, height = 6, useDingbats=FALSE)

ggvenn(x,
    columns = c("onepop", "twopop"),
  stroke_size = 0.5, set_name_size = 4
  )

dev.off()





# GO analyses Gene ontology



geneIDs <- c("GO:0030449","GO:0060294","GO:0009411","GO:0000271","GO:0006611","GO:0018279","GO:0001539","GO:0060285","GO:0090630","GO:0000381","GO:0045332","GO:0034204","GO:0007519","GO:0019935","GO:0070286","GO:0019933","GO:0010965","GO:0005976","GO:0006284","GO:0048024","GO:0097035","GO:1905818","GO:0030071","GO:0003341","GO:1902099","GO:0061245","GO:0071103","GO:0045197","GO:0035088","GO:0050684","GO:0015914","GO:0043547","GO:0060537","GO:0043484","GO:0030705","GO:0097553","GO:0043087","GO:0010970","GO:0035023","GO:0099111","GO:0010959","GO:0015748","GO:0007018","GO:0051345","GO:1904062","GO:0051056","GO:0010876","GO:0031503","GO:0006869","GO:1903311","GO:0044782","GO:0043269","GO:0006816","GO:0016311","GO:0060271","GO:0009101","GO:0061024","GO:0044093","GO:0043085","GO:0120031","GO:0141124","GO:0051336","GO:0030001","GO:1901137","GO:0006812","GO:0007017","GO:0033036","GO:0065009","GO:0006811","GO:0046907","GO:0070727","GO:0008104","GO:0051252","GO:0019219","GO:0006357","GO:0051649","GO:0006355","GO:2001141","GO:0051641","GO:0010556","GO:0051171","GO:0031326","GO:0009889","GO:0080090","GO:0010468","GO:0031323","GO:0051179","GO:0006810","GO:0051234","GO:0060255","GO:0019222","GO:0065007","GO:0050794","GO:0050789","GO:0008150","GO:0002376")


BiocManager::install("org.Gg.eg.db")
library(org.Gg.eg.db)
library(limma)
library(rrvgo)

df<-read.csv("finches_genes.csv", header=TRUE)
CRA_df = df[df$Species == 'CRA',]
FOR_df = df[df$Species == 'FOR',]
PAR_df = df[df$Species == 'PAR',]



# CRA

# Do limma's goana test for over-representation
gotest <- goana(CRA_df$Gene, species = "Gg", FDR = 0.05)
# gotest <- gotest[gotest$P.DE < 0.05,]
table(gotest$Ont)
  BP   CC   MF 
7366 1018 1903

> #Run rrvgo::calculateSimMatrix on only BP terms
BPterms <- rownames(gotest)[gotest$Ont %in% "BP"]
length(BPterms)
7366


simMatrix <- calculateSimMatrix(BPterms,
                                orgdb="org.Gg.eg.db",
                                ont="BP",
                                method="Rel")

reducedTerms <- reduceSimMatrix(simMatrix,
                                orgdb="org.Gg.eg.db")

heatmapPlot(simMatrix,
            annotateParent=TRUE,
            fontsize=6)

scatterPlot(simMatrix, reducedTerms)

treemapPlot(reducedTerms)


# FOR 

gotest <- goana(FOR_df$Gene, species = "Gg")

# Do limma's goana test for over-representation
gotest <- goana(FOR_df$Gene, species = "Gg")
gotest <- gotest[gotest$P.DE < 0.05,]
table(gotest$Ont)
BP CC MF 
77 26 31 

> #Run rrvgo::calculateSimMatrix on only BP terms
BPterms <- rownames(gotest)[gotest$Ont %in% "BP"]
length(BPterms)
77


simMatrix <- calculateSimMatrix(BPterms,
                                orgdb="org.Gg.eg.db",
                                ont="BP",
                                method="Rel")

reducedTerms <- reduceSimMatrix(simMatrix,
                                orgdb="org.Gg.eg.db")

heatmapPlot(simMatrix,
            annotateParent=TRUE,
            fontsize=6)

scatterPlot(simMatrix, reducedTerms)

treemapPlot(reducedTerms)


# PAR 


gotest <- goana(PAR_df$Gene, species = "Gg")

# Do limma's goana test for over-representation
gotest <- goana(r.genes, species = "Gg")
gotest <- gotest[gotest$P.DE < 0.05,]
table(gotest$Ont)
BP CC MF 
77 26 31 

> #Run rrvgo::calculateSimMatrix on only BP terms
BPterms <- rownames(gotest)[gotest$Ont %in% "BP"]
length(BPterms)
77


simMatrix <- calculateSimMatrix(BPterms,
                                orgdb="org.Gg.eg.db",
                                ont="BP",
                                method="Rel")

reducedTerms <- reduceSimMatrix(simMatrix,
                                orgdb="org.Gg.eg.db")

heatmapPlot(simMatrix,
            annotateParent=TRUE,
            fontsize=6)

scatterPlot(simMatrix, reducedTerms)

treemapPlot(reducedTerms)



# run it on all intersection genes
CRA_FOR <- intersect(CRA_df$Gene, FOR_df$Gene)

CRA_PAR <- intersect(CRA_df$Gene, PAR_df$Gene)

FOR_PAR <- intersect(FOR_df$Gene, PAR_df$Gene)

species_intersect <- c(CRA_FOR, CRA_PAR, FOR_PAR)

write.csv(species_intersect, "species_intersect_genes.csv")


CRA_raisd = CRA_df[CRA_df$Analysis == 'raisd',]
CRA_sweed = CRA_df[CRA_df$Analysis == 'sweed',]
CRA_tajima = CRA_df[CRA_df$Analysis == 'tajima',]
CRA_dxy = CRA_df[CRA_df$Analysis == 'dxy',]
CRA_fst = CRA_df[CRA_df$Analysis == 'fst',]
FOR_raisd = FOR_df[FOR_df$Analysis == 'raisd',]
FOR_sweed = FOR_df[FOR_df$Analysis == 'sweed',]
FOR_tajima = FOR_df[FOR_df$Analysis == 'tajima',]
FOR_dxy = FOR_df[FOR_df$Analysis == 'dxy',]
FOR_fst = FOR_df[FOR_df$Analysis == 'fst',]
PAR_raisd = PAR_df[PAR_df$Analysis == 'raisd',]
PAR_sweed = PAR_df[PAR_df$Analysis == 'sweed',]
PAR_tajima = PAR_df[PAR_df$Analysis == 'tajima',]
PAR_dxy = PAR_df[PAR_df$Analysis == 'dxy',]
PAR_fst = PAR_df[PAR_df$Analysis == 'fst',]

a <- intersect(CRA_sweed$Gene, CRA_tajima$Gene)
b <- intersect(CRA_fst$Gene, CRA_dxy$Gene)
c <- intersect(FOR_sweed$Gene, FOR_tajima$Gene)
d <- intersect(FOR_sweed$Gene, FOR_tajima$Gene)
e <- intersect(FOR_tajima$Gene, FOR_fst$Gene)
f <- intersect(FOR_dxy$Gene, FOR_fst$Gene)
g <- intersect(PAR_sweed$Gene, PAR_tajima$Gene)
h <- intersect(PAR_fst$Gene, PAR_dxy$Gene)
i <- intersect(PAR_fst$Gene, PAR_tajima$Gene)

analyses_intersect <- c(a, b, c, d, e, f, g, h, i)
write.csv(analyses_intersect, "analyses_intersect_genes.csv")


df<-read.csv("species_intersect_GO.csv", header=TRUE)
df<-read.csv("analyses_intersect_GO.csv", header=TRUE)

df<-read.csv("cra_sig_GO.csv", header=TRUE)
df<-read.csv("for_sig_GO.csv", header=TRUE)
df<-read.csv("par_sig_GO.csv", header=TRUE)
df<-read.csv("all_sig_GO.csv", header=TRUE)

simMatrix <- calculateSimMatrix(df$GO,
                                orgdb="org.Hs.eg.db",
                                ont="BP",
                                method="Rel")

scores <- setNames(-log10(df$raw.P.value), df$GO)

reducedTerms <- reduceSimMatrix(simMatrix,
                                scores,
                                threshold=0.7,
                                orgdb="org.Hs.eg.db")
                                
scatterPlot(simMatrix, reducedTerms)

treemapPlot(reducedTerms)




# make plots of analyses that supported overlapping genes
df<-read.csv("finches_genes.csv", header=TRUE)

CRA_FOR <- intersect(CRA_df$Gene, FOR_df$Gene)

CRA_PAR <- intersect(CRA_df$Gene, PAR_df$Gene)

FOR_PAR <- intersect(FOR_df$Gene, PAR_df$Gene)

species_intersect <- c(CRA_FOR, CRA_PAR, FOR_PAR)


library(dplyr)
df_new <- filter(df, Gene %in% species_intersect)

df2 <- aggregate(Analysis ~ Gene, data = df_new, FUN = c)

X <- aggregate(sample(c(T, F, NA), 100, r=T), list(rep(letters[1:4], 25)), summary)

X <- aggregate(sample(c(T, F, NA), 100, r=T), list(rep(letters[1:4], 25)), summary)
X <- cbind(X[-ncol(X)], X[[ncol(X)]])

 cbind(df2[-ncol(df2)], df2[[2]])

df3 <- df2 %>%
    count(Analysis)


barplot(n ~ Analysis, data = df3)
unlist(df3$Analysis)

types <- c(1,1,1,1,1,4,2,1,1,4)

names <- c("Tajima","Tajima","Tajima","Tajima, SweeD","dxy, FST","dxy ","dxy, FST ","dxy, FST","FST, Tajima","dxy")


data <- data.frame(name= c("Tajima","Tajima","Tajima","Tajima, SweeD","dxy, FST","dxy ","dxy, FST ","dxy, FST","FST, Tajima","dxy"), value= c(1,1,1,1,1,4,2,1,1,4)
)

barplot(height=data$value, names=data$name, col="darkblue",  horiz=T, las=1)
