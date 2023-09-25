## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse=TRUE,
  comment="#>"
)

## ---- eval=FALSE--------------------------------------------------------------
rCAM3 <- CAM3Run(data, K = 3, dim.rdc = 3, thres.low = 0.30, thres.high = 0.95)

## ---- echo=TRUE---------------------------------------------------------------
library(debCAM)
library(CAM3)

data(ratMix3) 
#ratMix3$X: X matrix containing mixture expression profiles to be analyzed
#ratMix3$A: ground truth A matrix containing proportions
#ratMix3$S: ground truth S matrix containing subpopulation-specific expression profiles

data <- ratMix3$X
#10000 genes * 21 tissues
#meet the input data requirements


## ---- echo=TRUE, results = 'hide', message=FALSE------------------------------
rCAM3 <- CAM3Run(data, K = 2:5, dim.rdc = 10, thres.low = 0.30, thres.high = 0.95)

## ---- echo=TRUE, results = 'hide', message=FALSE------------------------------
# Setting hyper-parameters
K = 2:5
dim.rdc = 10
thres.low = 0.30
thres.high = 0.95
radius.thres = 0.95
sim.thres = 0.95 
cluster.num = 50 
MG.num.thres = 20 
sample.weight = NULL

# Commencing pre-processing
PrepResult <- CAM3Prep(data, dim.rdc, thres.low, thres.high,
                       cluster.method = c('Fixed-Radius' ,'K-Means'),
                       radius.thres, sim.thres, cluster.num,
                       MG.num.thres, sample.weight)
## ---- message = FALSE --------------------------------------------------------
# Setting hyper-parameters 
fast.mode = TRUE

# Commencing marker gene selection
MGResult <- CAM3MGCluster(PrepResult, fast.mode)
MGResult <- MGResult[K]
names(MGResult) <- as.character(K)


## ---- message = FALSE---------------------------------------------------------
# Setting hyper-parameters     
generalNMF = FALSE    

# Commencing A and S estimation
ASestResultF <- lapply(MGResult, CAM3ASest, PrepResult, data,
                       1, generalNMF)
ASestResultB <- lapply(MGResult, CAM3ASest, PrepResult, data,
                       2, generalNMF)

ASestResult <- vector("list", length(K))
for (k in seq_along(K)){
  if (ASestResultF[[k]]@mdl < ASestResultB[[k]]@mdl)
    ASestResult[[k]] <- ASestResultF[[k]]
  else
    ASestResult[[k]] <- ASestResultB[[k]]
}
names(ASestResult) <- as.character(K)

## ---- message = FALSE---------------------------------------------------------
rCAM3 <- new("CAMObj",PrepResult=PrepResult, MGResult=MGResult,
             ASestResult=ASestResult)


# Post-CAM Analysis 

## Model Selection: Determining the Number of Cell Types

## ---- fig.height=5, fig.width=6, fig.align='left'-----------------------------
K = 2:5
plot(K, MDL(rCAM3, 1)@mdls, ylab='MDL Value')


## ---- echo=TRUE---------------------------------------------------------------
S<-rCAM3@ASestResult$'3'@Sest
MGlist <- cotMG(data,rCAM3@ASestResult$`3`@Sest, cos.thres = 0.99)


## Scatter Simplex Plot of Mixture
## ---- fig.height=6, fig.width=6, fig.align='right'----------------------------
simplexplot(data, rCAM3@ASestResult$`3`@Aest, MGlist = MGlist$mg.list,
            col = "grey", mg.col = c("red","orange","green"))


# Case Study
## GSE 64385

## ---- message = FALSE, results = 'hide'---------------------------------------
library(GEOquery)

gsm <- getGEO("GSE64385")
data <- 2^exprs(gsm[[1]])

rCAM<-CAM3Run(data=data, K= 2:8, dim.rdc = 10, thres.low = 0.30, 
              thres.high = 0.95)


## ---- fig.height=5, fig.width=6, fig.align='left'-----------------------------
K = 2:8
plot(K, MDL(rCAM, 1)@mdls, ylab='MDL Value')


## ---- echo=TRUE---------------------------------------------------------------
# Detecting cell type-specific markers
S<-rCAM@ASestResult$'6'@Sest
MGlist <- cotMG(data,rCAM@ASestResult$`6`@Sest, cos.thres = 0.99)


## ---- fig.height=6, fig.width=6, fig.align='left'-----------------------------
simplexplot(data, rCAM@ASestResult$`6`@Aest, MGlist = MGlist$mg.list,
            col = "grey", 
            mg.col = c("red","orange","green","purple","blue","pink"))

