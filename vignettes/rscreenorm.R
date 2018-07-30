## ----echo=FALSE----------------------------------------------------------
#Note the various macros within the `vignette` section of the metadata block above. These are required in order to instruct R how to build the vignette. Note that you should change the `title` field and the `\VignetteIndexEntry` to match the title of your vignette.

#The `html_vignette` template includes a basic CSS theme. To override this theme you can specify your own CSS in the document metadata as follows:

#    output: 
#      rmarkdown::html_vignette:
#        css: mystyles.css

## ----load.rscreenorm-----------------------------------------------------
library(rscreenorm)

## ----simul.defineSetup---------------------------------------------------
set.seed <- 12345 # 76543
n.clines <- 6     # number of cell lines
nlibrary <- 1000  # number of siRNAs in the library
n.cont <- 200     # number of controls - of negative and positive controls, each
n.reps <- 3        # Number of replicates per cell line (a single condition per cell line)
# Starting parameters of control distributions
mu.neg <- 0
mu.pos <- 1
sd.neg <- 0.1
sd.pos <- sd.neg

## ----simul.defineParameters----------------------------------------------
par.mat <- c(2, 6)
stretch <- rep( rep(c(0.6, 1, 1.5),2), each=n.reps) 
stretch <- stretch + rnorm(length(stretch), mean=0, sd=0.05)
mu.pos <- mu.pos*stretch
pi0 <- 0.8  # prop features with the same mean for all cell lines

## ----simul.SetMeans------------------------------------------------------
# Set means
meanv <- rbeta(nlibrary, shape1=par.mat[1], shape2=par.mat[2]) # means for null features, to be used for all replicates
mean.eff <- rep(0.5, ceiling((1-pi0)*nlibrary)) 
all.means1 <- matrix( rep(meanv, n.reps*n.clines/2), nrow=nlibrary, ncol=n.reps*n.clines/2 )
all.means2 <- rbind( all.means1[ 1:(nlibrary*pi0), ] , 
                     matrix(mean.eff, nrow=ceiling((1-pi0)*nlibrary), ncol=n.reps*n.clines/2 ) )
all.means <- cbind(all.means1, all.means2)

## ----simul.generateData--------------------------------------------------
data.mat <- matrix(0, nrow=nlibrary+2*n.cont, ncol=n.clines*n.reps) 
wt <- c(rep("sample", nlibrary), rep(c("neg","pos"), each=n.cont))

for(xi in 1:(n.clines*n.reps))
{
  if(is.null(dim(all.means))) {  
  data.mat[,xi] <- c(rnorm(nlibrary, mean=all.means, sd=2*sd.neg)*stretch[xi],
                     rnorm(n.cont, mean=mu.neg, sd=sd.neg),
                     rnorm(n.cont, mean=mu.pos[xi], sd=sd.pos))
  } else {
    data.mat[,xi] <- c(rnorm(nlibrary, mean=all.means[,xi], sd=2*sd.neg)*stretch[xi],
                       rnorm(n.cont, mean=mu.neg, sd=sd.neg),
                       rnorm(n.cont, mean=mu.pos[xi], sd=sd.pos))
  }
}

## ----simul.clines.labels-------------------------------------------------
cline.vec <- rep(1:n.clines, each=n.reps)
cl.reps <- paste( rep(paste("CL", 1:n.clines, sep=""), each=n.reps), rep(paste("R", 1:n.reps, sep=""), n.clines)  )
myp <- rep(c(rep(0.95,3),rep(0.70,3)), each=3)  # proportion of features to include in the core set
colnames(data.mat) <- cl.reps

## ----simul.annot---------------------------------------------------------
wtype <- factor(c(rep("sample", nlibrary), rep("neg", n.cont), rep("pos", n.cont)))
data.ann <- data.frame(geneID = paste("Gene", 1:nrow(data.mat), sep=""), 
                       wtype = wtype )  
data.screen <- list(data_ann = data.ann, 
                    data_only = data.mat, 
                    use_plate = FALSE)
class(data.screen) <- "rscreen.object"

## ----simul.pdata---------------------------------------------------------
f.cl <- factor(paste( rep(paste("CL",1:n.clines,sep=""),each=n.reps) ))
# Colours to be used for graphs, one corresponding to each cell line
vcols <- c("blue","darkblue","darkolivegreen4","red","deeppink4","violet")
col.cl <- rep(vcols,each=n.reps)
names(col.cl) <- as.character(f.cl)

pdata <- data.frame(names = cl.reps, clines = f.cl, treat = rep(0:1, each=n.reps))
pdata$cols_vec = var.in.colour(pdata$clines)[[1]]
pdata <- list(pdata = pdata, 
              cols_list = unique(pdata$cols_vec), 
              list_vars = colnames(pdata))

## ----simul.boxplots.rawData, fig.width=7, fig.height=4-------------------
par(las=2, mar=c(7.5, 3, 5, 1))
  mydata <- matrix(scr.data(data.screen), 
         nrow= nr.screen(data.screen)*nc.screen(data.screen),
         ncol=1)[, 1]
  wtype.vec <- rep(scr.ann(data.screen)[, "wtype"], nc.screen(data.screen) )
  cline.vec <- rep(scr.names(data.screen), each=nr.screen(data.screen) )
  all.cols <- rep(cols(pdata), 
                  each=nlevels(scr.ann(data.screen)[, "wtype"]))
  boxplot(mydata ~ wtype.vec + cline.vec, col = all.cols, pch=20, cex=.7,
          main = "Raw data per replicate and well type",
          cex.axis = 0.7)

## ----simul.get.leth.scores, fig.width=7, fig.height=4--------------------
lscores <- get.leth.scores(data.screen, rescale = FALSE)
par(las=2, mar=c(7.5, 3, 5, 1))
  mydata <- matrix(scr.data(lscores), 
         nrow= nr.screen(data.screen)*nc.screen(data.screen),
         ncol=1)[, 1]
 boxplot(mydata ~ wtype.vec + cline.vec, col = all.cols, pch=20, cex=.7,
         main = "Lethality scores per replicate and well type",
         cex.axis = 0.7)

## ----get.inv.set---------------------------------------------------------
inv.set <- get.inv.set(lscores, my_gamma = 1, set_shape = "left", var_type = "mad")

## ----plot.prop.inv, fig.width=7, fig.height=3----------------------------
par(las = 2)
thres <- get.th.inv.set(lscores, my_gamma = 1, var_type = "mad") 
myprops <- get.inv.set.prop(lscores, thres, plot = TRUE, cols_list = cols(pdata) )

## ----simul.get.qnorm,  fig.width=7, fig.height=4-------------------------
qnorm.scores <- get.qnorm(lscores, inv.set)
par(las=2, mar=c(7.5, 3, 5, 1))
  mydata <- matrix(scr.data(qnorm.scores), 
         nrow= nr.screen(data.screen)*nc.screen(data.screen),
         ncol=1)[, 1]
 boxplot(mydata ~ wtype.vec + cline.vec, col = all.cols, pch=20, cex=.7,
         main = "Rscreenorm scores per replicate and well type",
         cex.axis = 0.7)

## ----.simul.rscreenorm.wrapper-------------------------------------------
qnorm.scores <- get.rscreenorm(lscores, my_gamma = 1, var_type = "mad", set_shape = "left")

## ----simul.rscreenorm.densplots, fig.width=7, fig.height=4---------------
par(mfrow = c(1, 2), mar = c(3, 2.5, 4, 0.5))
pdw(data.screen, wtype = "sample", mytitle = "Raw data", lcex = .5)
pdw(qnorm.scores, wtype = "sample", mytitle = "Rscreenorm scores", lcex = .5)

## ----read.data.cris------------------------------------------------------
mydir <- getwd()
vcl <- c("HCT116_1", "RPE1")
filesToRead <- paste("readcount-", vcl, "-lib1.gz", sep="")
scr1 <- read.screen.data(mydir = mydir, filename = filesToRead[1],
                         n_ann_cols = 2, use_plate = FALSE )
scr2 <- read.screen.data(mydir = mydir, filename = filesToRead[2],
                         n_ann_cols = 2, use_plate = FALSE )

## ----read.ann.cris-------------------------------------------------------
data.ann.all <- read.delim(file.path(mydir, "annotation_mmc2.txt"), stringsAsFactors=FALSE)
wellType <- rep("sample", nrow(data.ann.all))
tab.target <- table(data.ann.all$Target)
cont.names <- names(tab.target[tab.target>10])
wellType[data.ann.all$Target %in% cont.names[3:5]] <- "neg"
wellType[data.ann.all$Target == "chr10Promiscuous"] <- "pos"
data.ann.all$wtype <- wellType

## ----order.data----------------------------------------------------------
ann.screen <- sapply(1:nr.screen(scr1), mysplit, 
                     char.vector = as.character(scr1[["data_ann"]]$GENE_CLONE),
                     split.by = "_", result.sel = 2)
mat.ann <- scr1[["data_ann"]][order(ann.screen), ]
mat.ann$gRNA <- ann.screen[order(ann.screen)]

scr1[["data_only"]] <- scr1[["data_only"]][order(ann.screen), ]
scr2[["data_only"]] <- scr2[["data_only"]][order(ann.screen), ]
scr1[["data_ann"]] <- scr2[["data_ann"]] <- mat.ann

## ----wtype.ann.cris------------------------------------------------------
scr1[["data_ann"]]$wtype <- scr2[["data_ann"]]$wtype <- wellType

## ----labels.cris---------------------------------------------------------
scr.names(scr1)
scr.names(scr2)

## ----newColNames---------------------------------------------------------
scr1.t <- sapply(1:nc.screen(scr1), mysplit, char.vector = scr.names(scr1), 
                 split.by = "_", result.sel = 2)
scr1.lab <- sapply(1:nc.screen(scr1), mysplit, char.vector = scr.names(scr1), 
                 split.by = "_", result.sel = 3)
colnames.new <- paste(vcl[1], scr1.t, c(scr1.lab[-length(scr1.lab)], ""), sep = ".")
colnames(scr1[["data_only"]]) <- colnames.new
pdata1 <- data.frame(names = colnames.new, clines = rep(vcl[1], length(scr1.t)), treat = scr1.t, cols_vec =rep("blue", nc.screen(scr1)))
pd1 <- list(pdata = pdata1, cols_list = "blue", list_vars =colnames(pdata1))
#
scr2.t <- c(scr.names(scr2)[1], 
            substr(scr.names(scr2)[-1], start=1, stop=nchar(scr.names(scr2)[-1])-1 ) )
scr2.lab <- c("", substr(scr.names(scr2)[-1], start=nchar(scr.names(scr2)[-1]),
                          stop=nchar(scr.names(scr2)[-1])) )
colnames.new <- paste(vcl[2], scr2.t, scr2.lab, sep = ".")
colnames(scr2[["data_only"]]) <- colnames.new
pdata2 <- data.frame(names = colnames.new, clines = rep(vcl[2], length(scr2.t)), treat = scr2.t, cols_vec =rep("red", nc.screen(scr2)))
pd2 <- list(pdata = pdata2, cols_list = "red", list_vars =colnames(pdata2))

## ----combine.screens.cris------------------------------------------------
scr.data <- combine.screens(list(scr1, scr2))
pdata <- combine.pdata(list(pd1, pd2))
col.cl <- pdata[["pdata"]]$cols_vec
names(col.cl) <- as.character(pdata[["pdata"]]$clines)

## ----densityPlotsWithControls.cripsr, fig.width=7, fig.height=3----------
par(mfrow=c(1,3), mar=c(4, 3, 1, 1))
pdmw(scr.data, pdata, use.wt = c("sample", "pos", "neg"), na.rm = TRUE,
                 plot = TRUE, rescale = TRUE, scale.fun = "asinh", myxlab = "asinh values",
                 mytitle = FALSE, add.text = TRUE, add.leg = TRUE, lcex=.5)

## ----count.leth.posControls----------------------------------------------
mycut <- rep(4, nc.screen(scr.data))
mycut[ pdata[["pdata"]]$names == "RPE1" ] <- 2
pos.0 <- NULL
for(xi in 1:nc.screen(scr.data))
   pos.0 <- cbind(pos.0, scr.data(scr.data)[ scr.ann(scr.data)[["wtype_ann"]][, xi] == "pos", xi] <= mycut[xi])
rownames(pos.0) <- scr.data[["data_ann"]][["data_ann"]][ scr.data[["data_ann"]][["wtype_ann"]][, 1] == "pos", "gRNA"]
n.0s <- apply(pos.0, 1, sum)
table(n.0s)

## ----cris.newposcont-----------------------------------------------------
oldPos <- rownames(pos.0)[n.0s < 18]
wtype.new <- scr.ann(scr.data)[["wtype_ann"]]
wtype.new[ scr.ann(scr.data)[["data_ann"]]$gRNA %in% oldPos, ] <- "random"
scr.data[["data_ann"]][["wtype_ann"]] <- wtype.new

## ----cris.rscreenorm.wrapper---------------------------------------------
lscores <- get.leth.scores(scr.data, rescale = TRUE, scale.fun = "asinh")
norm.scores <- get.rscreenorm(lscores, prop = 0.95, var_type = "mad", set_shape = "left")

## ----cris.rscreenorm.densplots, fig.width=7, fig.height=4----------------
par(mfrow = c(1, 2), mar = c(5, 2.5, 5, 0.5))
pdw(lscores, wtype = "sample", mytitle = "Lethality scores", lcex=.4,
    xlim = range(scr.data(lscores)), myxlab = "", leg.side = "left")
pdw(norm.scores, wtype = "sample", mytitle = "Rscreenorm scores", lcex=.4,
    xlim = range(scr.data(norm.scores)), myxlab = "", leg.side = "left")

## ----install.libraries, eval=FALSE---------------------------------------
#  source("http://bioconductor.org/biocLite.R")
#  # Before installation, run on terminal:
#  # sudo apt-get install libmpfr-dev -y
#  biocLite("Mulder2012")
#  biocLite("Rmpfr")

## ----muld.load.data, message = FALSE, warning = FALSE--------------------
library(Mulder2012)
# Raw data
data(Mulder2012.rawScreenData)
mrd.all <- rawScreenData
# Annotation
data(Mulder2012.rawScreenAnnotation)
mann <- rawScreenAnnotation

## ----muld.cols.data------------------------------------------------------
colnames(mrd.all)

## ----def.muld.dataPerchannel---------------------------------------------
mrd <- mrd.all[mrd.all$CHANNEL == "700", ] # Check if this is the channel to be used

## ----muld.cols.ann-------------------------------------------------------
colnames(mann)

## ----muld.check.plates.reps----------------------------------------------
table(mrd$PLATE)

## ----muld.poscontrols----------------------------------------------------
table(as.character(mann$SYMBOL)[mann$ACC.NUMB=="control"], mann$PLATE[mann$ACC.NUMB=="control"])

## ----muld.annot----------------------------------------------------------
wtype <- rep("sample", nrow(mrd))
wtype[mann$ACC.NUMB == "control"] <- rep("pos", sum(mann$ACC.NUMB == "control"))
data.ann <- data.frame(geneID = as.character(mann$SYMBOL), 
                       accnum = mann$ACC.NUMB,
                       wtype = wtype,
                       plate = mann$PLATE)  
data.screen <- list(data_ann = data.ann, data_only = as.matrix(mrd[, -(1:3)]), 
                    use_plate = TRUE, hasWtype = TRUE)
class(data.screen) <- "rscreen.object"

## ----muld.def.pdata------------------------------------------------------
treat <- colnames(mrd)[seq(from=4, by=3, length.out=5)] 
treat <- substr(treat, start=1, stop=nchar(treat)-2)
pdata <- data.frame( names = colnames(mrd)[-(1:3)],
                     treat = rep(treat, each=3),
                     replicate = rep(1:3, 5),
                     cols_vec = var.in.colour(rep(treat, each=3))[[1]])
pdata <- list(pdata = pdata, 
              cols_list = unique(pdata$cols_vec), 
              list_vars = colnames(pdata))

## ----muld.plot.dens.rawData, fig.width=6, fig.height=6-------------------
pdw(data.screen, wtype = "sample", mytitle = "Raw data", mycols = cols(pdata), lcex=.7 )

## ----muld.plot.dens.log2Data, fig.width=7, fig.height=6------------------
# Create a vector of unique colours corresponding to the treatment levels
vcols <- var.in.colour(pdata[["pdata"]]$treat)[[2]]
# Attribute names to this vector representing the treatment levels
names(vcols) <- var.in.colour(pdata[["pdata"]]$treat)[[3]]
# Make the density plots
accnum <- scr.ann(data.screen)$accnum
par(mfrow=c(2, 3), mar = c(5, 2.5, 5, 0.5))
for(xi in levels(pdata[["pdata"]]$treat))
{ 
 pdw(data.screen, 
     sel.reps = scr.names(data.screen)[ pdata[["pdata"]]$treat == xi ], wtype = "sample", 
     mycols = rep(vcols[xi], sum(pdata[["pdata"]]$treat == xi)),
     rescale = TRUE, scale.fun = "log2", 
     mytitle = "log2-raw data", lcex=.5)
  mypch <- 1:2
  mcont <- log2(scr.data(data.screen)[accnum =="control", pdata[["pdata"]]$treat == xi])
  mcann <- scr.ann(data.screen)[accnum == "control", ]
  xk <- 1
  for(xi in unique(as.character(mcann$geneID)))
  { 
    points(as.matrix(mcont[mcann$geneID==xi, ]), rep(0, ncol(mcont) * nrow(mcont)/2), 
           pch = mypch[xk], col="red", cex=2)
    xk <- xk+1
  }
  legend("topleft", 
         legend=unique( scr.ann(data.screen)$geneID[accnum == "control"] ), 
         pch=1:2, col="red", cex=.8)

}

## ----muld.boxplots.rescaledrawdata, fig.width=7, fig.height=3------------
myylim <- range(log2(scr.data(data.screen)))
par(mfrow=c(1, 3), mar = c(5, 2.5, 5, 0.5))
for(xj in 1:nc.screen(data.screen))
{  
    boxplot(log2(scr.data(data.screen))[, xj] ~ 
              factor(data.screen[["data_ann"]]$plate),
            main=paste(scr.names(data.screen)[xj]),
            col=cols(pdata)[xj], ylim=myylim )
  mcont <- log2(scr.data(data.screen)[accnum =="control", xj])
  mcann <- scr.ann(data.screen)[accnum =="control", ]
    xk <- 1
    for(xi in unique(as.character(mcann$geneID)))
    { 
      points(mcann$plate[mcann$geneID==xi], mcont[mcann$geneID==xi],  
             pch = mypch[xk], col="black", cex=2)
      xk <- xk+1
  }
}

## ----muld.def.centerPerPlateTG1------------------------------------------
tgdata <- log2(scr.data(data.screen))[accnum == "control", ]
tgplate <- scr.ann(data.screen)$plate[accnum == "control"]
data.c <- scr.data(data.screen)

for(xj in 1:nc.screen(data.screen))
{  
  mtgp <- tapply(tgdata[, xj], INDEX=list(factor(tgplate)), mean)
  mdatap <- NULL
  for(xp in 1:4)
  {
      mdata <- log2(scr.data(data.screen))[scr.ann(data.screen)$plate == xp, xj] - mtgp[xp]
      mdatap <- c(mdatap, mdata)
  }
  data.c[, xj] <- mdatap
}
mylimc <- range(data.c)

## ----muld.boxplots.logCdata, fig.width=7, fig.height=3-------------------
par(mfrow=c(1, 3), mar = c(5, 2.5, 5, 0.5))
for(xj in 1:nc.screen(data.screen))
{  
  boxplot(data.c[, xj] ~ factor(scr.ann(data.screen)$plate), main=paste(scr.names(data.screen)[xj]),
            col=cols(pdata)[xj], ylim=mylimc )
  mcont <- data.c[accnum =="control", xj]
  mcann <- scr.ann(data.screen)[accnum =="control", ]
  xk <- 1
  for(xi in unique(as.character(mcann$geneID)))
    { 
      points(mcann$plate[mcann$geneID==xi], mcont[mcann$geneID==xi],  
             pch = mypch[xk], col="black", cex=2)
      xk <- xk+1
   }
}

## ----muld.newrscreen.data------------------------------------------------
data.screen.c <- data.screen
data.screen.c[["data_only"]] <- data.c

## ----muld2.plot.dens.log2Data.all, fig.width=7, fig.height=4-------------
par(mfrow=c(1, 2), mar = c(5, 2.5, 5, 0.5))
pdw(data.screen, mycols = as.character(cols(pdata)), rescale = TRUE, scale.fun = "log2",
    xlim = c(myylim[1], myylim[2]*1.2), 
    myxlab = "log2-raw values", lcex=.5)

pdw(data.screen.c, mycols = as.character(cols(pdata)), rescale = FALSE, 
    xlim = c(mylimc[1], mylimc[2]*1.2), 
    myxlab = "log2 plate-centered values", lcex=.5)

## ----muld.lscores--------------------------------------------------------
mann$wtype <- rep("sample", nrow(mann))
lscores <- list(data_only = data.c, data_ann = mann, use_plate = FALSE)
class(lscores) <- "rscreen.object" 
norm.data <- get.rscreenorm(lscores, var_type = "mad", prop = 0.67, set_shape = "left")

## ----muld.plot.dens.rscreenormData, fig.width=7, fig.height=4------------
par(mfrow=c(1, 2), mar = c(5, 2.5, 5, 0.5)) 
pdw(data.screen, mycols = as.character(pdata$pdata$cols_vec), myxlab = "raw values", lcex=.5)
pdw(norm.data, mycols = as.character(pdata$pdata$cols_vec), myxlab = "rscreenorm scores", lcex=.5)

