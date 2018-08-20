## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE,message=FALSE,warning=FALSE)

## ------------------------------------------------------------------------
library(LRC3)
SampleData[1:5,1:10]

## ------------------------------------------------------------------------
data_list <- PrepareData(SampleData)
data_list[[1]][1:5,1:5]
data_list[[2]][1:5]
data_list[[3]][1:5]

## ------------------------------------------------------------------------
head(LRpairs)

## ------------------------------------------------------------------------
LRC3_list <- LRC3_INF(data_list)

## ----echo=TRUE,fig.height=3.5, fig.width=4, fig.show='hold',fig.align='center'----
LRC3_ContactMatrixPlot(LRC3_list,IndexNumber =1)


## ----echo=TRUE,fig.height=7, fig.width=6, fig.show='hold',fig.align='center'----
LRC3_Connection(LRC3_list,IndexOfPlot = 1)

## ----echo=TRUE,fig.height=7, fig.width=6, fig.show='hold',fig.align='center'----
LRC3_Connection(LRC3_list,IndexOfPlot = 1,By_Proportion = TRUE)

## ----echo=TRUE,fig.height=7, fig.width=6.5, fig.show='hold',fig.align='center'----
LRC3_Connection(LRC3_list,IndexOfPlot = 2)

## ----echo=TRUE,fig.height=3, fig.width=6.5, fig.show='hold',fig.align='center'----
LRC3_Connection(LRC3_list,IndexOfPlot = 3)

## ----echo=TRUE,fig.height=3, fig.width=6.5, fig.show='hold',fig.align='center'----
LRC3_Connection(LRC3_list,IndexOfPlot = 3,By_Proportion = TRUE)

## ----fig.align='center'--------------------------------------------------
LRC3_CellTypeConnection(LRC3_list,CellType_Index = 1,TargetCellType_IN = TRUE)

## ----fig.align='center'--------------------------------------------------
LRC3_CellTypeConnection(LRC3_list,CellType_Index = 1,
                        TargetCellType_OUT = TRUE)

## ----fig.align='center'--------------------------------------------------
LRC3_CellTypeConnection(LRC3_list,CellType_Index = 1,Combine_InOut = TRUE)

## ----fig.align='center'--------------------------------------------------
LRC3_CellTypeConnection(LRC3_list,CellType_Index = 1,TargetCellType_IN = TRUE,
                        By_Proportion = TRUE)

## ----fig.align='center'--------------------------------------------------
LRC3_CellTypeConnection(LRC3_list,CellType_Index = 1,
                        TargetCellType_OUT = TRUE,By_Proportion = TRUE)

## ----fig.align='center'--------------------------------------------------
LRC3_CellTypeConnection(LRC3_list,CellType_Index = 1,Combine_InOut = TRUE,By_Proportion = TRUE)

## ----echo= TRUE,fig.height=3.5, fig.width=4, fig.show='hold',fig.align='center'----
LRC3_list <- LRC3_INF(data_list,PercentForPval = 0.99)
LRC3_ContactMatrixPlot(LRC3_list,IndexNumber =1)

## ----echo= TRUE,fig.height=3.5, fig.width=4, fig.show='hold',fig.align='center'----
LRC3_list <- LRC3_INF(data_list,PvalAsThreshold = FALSE, FixedThreshold = 1)
LRC3_ContactMatrixPlot(LRC3_list,IndexNumber =1)


## ------------------------------------------------------------------------
LRC3_LRPairsCellTypeSearch(LRC3_list,CellType1 = 4,CellType2 = 5)

## ------------------------------------------------------------------------
LRC3_LRPairsCellTypeSearch(LRC3_list,CellType1 = 5,CellType2 = 4)

## ------------------------------------------------------------------------
TotalListOfLRpairs <- LRC3_TotalListOfLRpairs(LRC3_list)
TotalListOfLRpairs[1:20,]

