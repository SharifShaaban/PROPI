#!/bin/Rscript

# http://www.stat.berkeley.edu/~s133/Cluster2a.html
# http://www.stat.berkeley.edu/~s133/saving.html
# http://stackoverflow.com/questions/19245291/nas-introduced-by-coercion-during-cluster-analysis-in-r
# http://www.statmethods.net/advgraphs/parameters.html
# https://stat.ethz.ch/R-manual/R-devel/library/base/html/nrow.html
# https://stat.ethz.ch/R-manual/R-devel/library/base/html/any.html
# https://stat.ethz.ch/R-manual/R-devel/library/base/html/write.html

#############################################################################
# - Obtain Euclidian distances, followed by Hierachical clustering.         #
# - A for loop is done investigating the within cluster Eucledian distances #
# when looking at the range of possible cluster number: from 0 to as many   #
# unique binary lines. The thresholds are then used to determine the best   #
# clustering fitting the guidelines.                                        #
#############################################################################

library(magrittr)

file_matrix<-(commandArgs(TRUE)[1])
pan_matrix<-read.delim(file_matrix, sep="\t", header=TRUE)
names(pan_matrix)[1] <- "names"
pan_matrix$names <- as.character(pan_matrix$names)
isolate_number <- nrow(pan_matrix)
myData <- list(names = pan_matrix[,1], mat = as.matrix(pan_matrix[,2:(ncol(pan_matrix) - 1)])) # the last column is empty

a <- Sys.time() # ~ 2 min
nDup <- sum(duplicated(myData$mat), na.rm=TRUE) # must remove duplicate lines due to hclust method

myDist <- dist(myData$mat)
myHierachy <- hclust(myDist, method = "complete") # change methods if you want to

t0 <- Sys.time()

myClustersList <- lapply(1:(isolate_number - nDup), function(clusterNumbers) {
    myClust <- cutree(myHierachy, k = clusterNumbers)
    myDists <- sapply(1:max(myClust), function(nClust) {
        myData$mat[which(myClust == nClust), , drop = FALSE] %>% dist %>% max 
        # possible optimization: do not recompute dist, but subset 'myDist'
        # failing to use drop = FALSE is a major source of bugs
        # median might be better than mean
    })
    return(list("clusters" = myClust, withinss = myDists))
}) 
Sys.time() - t0

for (maxWithinss in c(0, 1.5, 3, 4.5, 6, 9))
{
    for (x in myClustersList)
    {
    # withinss distances within a cluster / need to change withinss adequate to matrix size
        if (all(x$withinss <= maxWithinss)) {
                finalCluster <- x
                break
        }      
    }

    myResults <- data.frame("name" = myData$names, "cluster" = finalCluster$cluster)
    myResults <- myResults[order(myResults$cluster),]
    cluster_number <- max(finalCluster$clusters)

    for (y in 1:cluster_number)
    {
        z <- capture.output(print(myResults[myResults$cluster == y,][1], digits = NULL, quote = FALSE, row.names = FALSE, right = FALSE))
        filename<-paste("./temp/groupings/cluster_", y, "_list_", maxWithinss,".txt", sep="")
        write(z, file = filename, sep = " ")
    }
    print("Done")
}
