# ******************************************************** #
# R-Ladies BCN: R for bioinformatics
# 2017/06/12
# ******************************************************** #

# Set up the packages ---------------------------------
source("https://bioconductor.org/biocLite.R")
biocLite(c("estrogen", "limma", "affy", "affyPLM",
           "hgu95av2.db", "hgu95av2cdf", "made4", "GOstats"))

# Load data -------------------------------------------
library(estrogen)
# Where is our example data? 
system.file("extdata", package="estrogen")

setwd(system.file("extdata", package="estrogen"))
# Let's load it
library(limma)
(targets <- readTargets("estrogen.txt", sep=""))

library(affy)
affybatch <- ReadAffy(filenames=targets$filename)


# Normalisation and QC --------------------------------
# Quality control (QC)
windows()
par(mfrow = c(2,4))
MAplot(object = affybatch)
par(mfrow = c(1,1))

library(made4)
overview(affybatch)

# Normalisation methodology: RMA (robust multi-array)
eset <- rma(affybatch)
eset

overview(eset) # Data is not perfectly normalised. =(
library(affyPLM)
par(mfrow=c(2,4))
MAplot(eset) # Results are good enough for today :) 
par(mfrow=c(1,1))


# Calculating differential expression -----------------

factor <- factor(paste0(targets$estrogen,targets$time.h))
factor

design <- model.matrix(~0+factor)
colnames(design) <- levels(factor)
design

fit <- lmFit(eset, design)
names(fit)
fit
# fit$coef are mean log-exp for each treatment combination.

cont.matrix <- makeContrasts(E10="present10-absent10",
                             E48="present48-absent48",
                             Time="present48-present10",
                             levels=design)
cont.matrix

fit2  <- contrasts.fit(fit, cont.matrix)
fit2  <- eBayes(fit2)

topTable(fit2) # All together
toptable(fit2, n = 5, adjust.method = 'holm') 
  # Change the n to display, and 
    # the multiple testing adjustment method


# We did not make that contrast
topTable(fit2, coef = 4, adjust.method = 'BH')

Res <- decideTests(fit2)
summary(Res) 
  # 1 = Up-regulated
  # -1 = Down-regulated
vennDiagram(Res)
vennDiagram(Res, include = c('up', 'down'),
            counts.col = c('seagreen3', 'tomato1'))

# Annotation ------------------------------------------
# 31798_at does not mean anything to anybody. 
# library(hgu95av2cdf) 
  # Only if we want to have the full probe-annotation map

library(hgu95av2.db)
# exploring the package:
  # HGNC <- unlist(as.list(hgu95av2GENENAME))

DEG.ID <- rownames(toptable(fit2,coef = 3, number = nrow(fit2), p.value = 0.05))
All.ID <- rownames(toptable(fit2, coef = 1, number = nrow(fit2)))
DEG.HGNC <- as.character(unlist(lapply(mget(DEG.ID,
                                       env=hgu95av2GENENAME),
                                  function (symbol){ 
                                    return(paste(symbol,
                                                 collapse="; "))
                                    })))

head(DEG.HGNC)

# Gene ontology annotation ----------------------------
# What are these genes doing?
library(GOstats)
DEG <- as.character(unlist(lapply(mget(DEG.ID,
                                       env=hgu95av2ENTREZID),
                                  function (symbol){ 
                                    return(paste(symbol,
                                                 collapse="; "))
                                  })))
ALL <- as.character(unlist(lapply(mget(All.ID,
                                       env=hgu95av2ENTREZID),
                                  function (symbol){ 
                                    return(paste(symbol,
                                                 collapse="; "))
                                  })))

params <- new('GOHyperGParams', annotation = 'org.Hs.eg', 
              geneIds = DEG, universeGeneIds = ALL, 
              ontology = c('BP'), 
              pvalueCutoff = 0.05, testDirection = 'over')
hg.test <- hyperGTest(params)
hg.res <- summary(hg.test)
hg.test.adjp <- p.adjust(pvalues(hg.test), 'fdr')
sigGO.ID <- names(hg.test.adjp[hg.test.adjp < 0.05])
hg.res[hg.res[,1] %in% sigGO.ID, "Term"]

# Optional
windows()
plot(makeGOGraph(DEG, chip = 'hgu95av2.db'))
