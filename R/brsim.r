#' Brainerd-Robinson similarity coefficient matrix calculation, with permutation-based p-values and optional clustering
#'
#' @description This function calculates the Brainerd-Robinson (BR) similarity coefficient for each pair of row
#' of the input table (Robinson-Brainerd 1951, 1952).
#' It also performs a permutation test to assess the significance of each BR coefficient (DeBoer-Kintigh-Rostoker 1996), and allows to carry
#' out a hierarchical agglomerative clustering. An optimal cluster solution can be established using the silhouette method
#' (see details provided below). The function produces a correlation matrix in tabular form, which is also visually plotted as an heatmap.
#' In the heatmap (which is built using the \code{corrplot} package), the size and the color of the
#' squares are proportional to the Brainerd-Robinson coefficients. Optionally, the heatmap can be reordered on the basis
#' of the hierachical clustering, with clusters enclosed by red rectangles.\cr
#' Visit this \href{https://drive.google.com/file/d/1LSC5VE_QZNM2KOCPkItn--1zsjfXnHbH/view?usp=share_link}{LINK} to access the package's vignette.\cr
#'
#' @details \strong{Permutation-based p-values}\cr
#'
#' The rationale behind the p-value calculation is as follows: for each pair of assemblages in the input data,
#' the function first calculates the observed Brainerd-Robinson (BR) coefficient. This is a measure of the similarity
#' between the two assemblages. The function then performs a certain number of permutations (default is 1000).
#' In each permutation, it generates two new assemblages (each featuring a sample size corresponding to the size of
#' each assemblage being compared) by randomly sampling from the global pool (the combined data of all assemblages),
#' and calculates the BR coefficient for this new pair of assemblages (see DeBoer-Kintigh-Rostoker 1996). This creates
#' a distribution of BR coefficients that we would expect to see by chance. The p-value is then calculated as the
#' proportion of the permuted BR coefficients that are less than or equal to the observed BR coefficient. A small
#' p-value (typically < 0.05) suggests that the observed similarity between the two assemblages is statistically
#' significant; it is unlikely to have occurred just by chance.\cr
#'
#' In simple terms, the p-value calculation is trying to answer the question: if there were no real
#' similarity between these two assemblages, what is the probability that I would observe a similarity as extreme as
#' (or more extreme than) the one I actually observed, just by chance?\cr
#'
#' The p-values are returned in two matrices: in the first, the p-values are reported as they are, whereas in the second
#' they are classified as <0.05, <0.01, <0.001, or not significant.\cr
#'
#' \strong{Hierarchical agglomerative clustering}\cr
#'
#' By setting the parameter \code{clust} to \code{TRUE}, the units (rows) for which the BR coefficients have been
#' calculated will be clustered. Note that the clustering is based on a dissimilarity matrix which
#' is internally calculated as the maximum value of the BR coefficient (200) minus the observed BR coefficient.
#' This allows a simpler reading of the dendrogram which is produced by the function, where the less dissimilar
#' (i.e., more similar) units will be placed at lower levels, while more dissimilar (i.e., less similar) units
#' will be placed at higher levels within the dendrogram. The latter depicts the hierarchical clustering based
#' (by default) on the Ward's agglomeration method; rectangles identify the selected cluster partition. Optionally,
#' by setting the \code{sort.map} to \code{TRUE}, the heatmap can be reordered on the basis of the hierarchical clustering,
#' with clusters indicated by red rectangles. The number of clusters indicated depends on what requested by the user (see
#' the next section). Note that, internally, the reordering is based on the same agglomeration method
#' selected by the user via the \code{aggl.method} parameter, which is set to \code{ward.D2} by default.\cr
#'
#' \strong{Number of clusters and silhouette method}\cr
#'
#' Besides the dendrogram, a silhouette plot is produced, which allows to measure how 'good' is the selected cluster solution.
#' If the parameter \code{part} is left empty (default), an optimal cluster solution is
#' obtained. The optimal partition is selected via an iterative procedure which identifies at which
#' cluster solution the highest average silhouette width is achieved. The cluster solution ranges from a minimum of
#' 2 to a maximum which is equal to the number of units (i.e., the rows of the input dataset) minus 1.
#' The number of units essentially represents the maximum number of clusters that could potentially be formed
#' if each row were its own cluster. However, since a cluster solution requires at least two groups,
#' the maximum number of meaningful clusters is one less than the number of rows.
#' If a user-defined partition is needed, the user can input the desired number of clusters using the parameter \code{part}.
#' In either case, an additional plot is returned besides the cluster dendrogram and the silhouette plot; it
#' displays a scatterplot in which the cluster solution (x-axis) is plotted against the average
#' silhouette width (y-axis). A black dot represents the partition selected either by the iterative
#' procedure or by the user. Note that in the silhouette plot, the labels on the left-hand side of the
#' chart show the units' names and the cluster number to which each unit is closer.\cr
#'
#' The silhouette plot is obtained from the \code{silhouette()} function out from the \code{cluster} package.
#' For a detailed description of the silhouette plot, its rationale, and its interpretation, see
#' Rousseeuw 1987.
#'
#' \strong{Descriptive by-cluster dotplots}\cr
#'
#' The function also provides a Cleveland's dotplots that represent by-cluster proportions. The
#' clustered units are grouped according to their cluster membership, the frequencies are summed, and
#' then expressed as percentages. The latter are represented by the dotplots, along with the average
#' percentage. The latter provides a frame of reference to understand which percentage is below,
#' above, or close to the average. The raw data on which the plots are based are stored in the
#' list returned by the function (see below).\cr
#'
#' @param df A table (dataframe format) where each row represents an assemblage and each column represents an item.
#' @param num.perm A numeric value indicating the number of permutations to perform in each test (default is 1000).
#' @param clust TRUE (default) or FALSE if the user does or does not want a agglomerative hierarchical clustering to be performed.
#' @param part Desired number of clusters; if NULL (default), an optimal partition is calculated (see Details).
#' @param aggl.meth Agglomeration method ("ward.D2" by default) to be used; the selected method is internally used for the reordering
#' of the heatmap if \code{order.map} is set to \code{TRUE}; for other methods see \code{\link[stats]{hclust}}.
#' @param sort.map TRUE or FALSE (default) if the user does or does not want the rendered heatmap to be ordered on the basis of
#' the selected hierachical clustering.
#' @param oneplot TRUE (default) or FALSE if the user wants or does not want the plots to be visualized in a single window.
#' @param cex.dndr.lab Set the size of the labels used in the dendrogram.
#' @param cex.sil.lab Set the size of the labels used in the silhouette plot.
#' @param cex.dot.plt.lab Set the size of the labels used in the Cleveland's dotplots representing by-cluster proportions.
#'
#'@return The function returns a list storing the following components \itemize{
##'  \item{BR.similarity.matrix: }{similarity matrix reporting the BR coefficients}
##'  \item{P-value.matrix: }{matrix reporting the permuted p-values}
##'  \item{classified.P-values.matrix: }{matrix reporting the permuted p-value classified as <0.05, <0.01, <0.001, or not significant}
##'  \item{BR.distance_matrix: }{distance matrix on which the hierarchical clustering is performed (returned if clustering is selected)}
##'  \item{avr.silh.width.by.n.of.clusters: }{average silhouette width by number of clusters (returned if clustering is selected)}
##'  \item{partition.silh.data: }{silhouette data for the selected partition (returned if clustering is selected)}
##'  \item{data.with.cluster.membership: }{copy of the input data table with an additional column storing the cluster membership for each row (returned if clustering is selected)}
##'  \item{by.cluster.proportion: }{table reporting the proportion of column categories across each cluster; rows sum to 100 percent (returned if clustering is selected)}
##' }
#'
#' @keywords similarity
#'
#' @export
#'
#' @importFrom corrplot corrplot
#' @importFrom stats as.dist hclust cutree aggregate rect.hclust
#' @importFrom cluster silhouette
#' @importFrom RcmdrMisc assignCluster
#' @importFrom grDevices colorRampPalette
#' @importFrom graphics axis layout par
#'
#' @references Robinson, W. S. (1951). A Method for Chronologically Ordering Archaeological Deposits.
#' In American Antiquity (Vol. 16, Issue 4, pp. 293–301). Cambridge University Press.
#'
#' @references Robinson, W. S., & Brainerd, G. W. (1952). Robinson’s Coefficient of Agreement – A Rejoinder.
#' In American Antiquity (Vol. 18, Issue 1, pp. 60–61). Cambridge University Press.
#'
#' @references Rousseeuw P J. (1987). Silhouettes: A graphical aid to the interpretation and validation of cluster
#' analysis. In Journal of Computational and Applied Mathematics 20, 53-65.
#'
#' @references DeBoer, W. R., Kintigh, K., & Rostoker, A. G. (1996). Ceramic Seriation and Site Reoccupation in Lowland South America.
#' In Latin American Antiquity (Vol. 7, Issue 3, pp. 263–278). Cambridge University Press.
#'
#'
#' @examples
#' # build a toy dataset (subset of the 'Nelson' dataset out of the 'archdata' package )
#'
#' mytable <- structure(list(Biscuit = c(10, 17, 2, 10, 2, 1),
#' Type_I = c(2,2, 10, 40, 118, 107),
#' Type_II_Red = c(24, 64, 68, 91, 45, 3),
#' Type_II_Yellow = c(23, 90, 18, 20, 1, 0),
#' Type_II_Gray = c(34,76, 48, 15, 5, 0)),
#' row.names = c("1", "2", "3", "7", "8", "9"),
#' class = "data.frame")
#'
#' # run the function and store the results in the 'result' object
#'
#' result <- brsim(mytable, clust=TRUE)
#'
#' # same as above, but with an user-defined cluster partition
#'
#' result <- brsim(mytable, clust=TRUE, part=3)
#'
#' # same as above, but rendering with a reordered heatmap
#'
#' result <- brsim(mytable, clust=TRUE, part=3, sort.map=TRUE)
#'
#'
#' @seealso \code{\link[corrplot]{corrplot}} , \code{\link[cluster]{silhouette}}
#'
#'
brsim <- function(df, num.perm=1000, clust=FALSE, part=NULL, aggl.meth="ward.D2", sort.map=FALSE, oneplot=TRUE, cex.dndr.lab = 0.70, cex.sil.lab = 0.70, cex.dot.plt.lab = 0.75){
  # Save current par settings
  oldpar <- par(no.readonly = TRUE)
  # Ensure settings are restored when function exits
  on.exit(par(oldpar))

  # Convert df to matrix
  matrix <- as.matrix(df)
  # Get the row names to be later attached to the rows and columns of the output matrices
  df.row.names <- rownames(df)

  # Local function for Brainerd-Robinson Similarity Coefficient
  brainard_robinson <- function(assemblage1, assemblage2){
    prop_assemblage1 <- assemblage1 / sum(assemblage1)
    prop_assemblage2 <- assemblage2 / sum(assemblage2)
    min_sum <- sum(pmin(prop_assemblage1, prop_assemblage2))
    return(200 * min_sum)
  }

  # Local function for permutation test
  perm_test_br <- function(assemblage1, assemblage2, num.perm=1000){
    observed_br <- brainard_robinson(assemblage1, assemblage2)
    permuted_br <- numeric(num.perm)

    # Calculate the global pool of all assemblages
    global_pool <- sum(matrix)

    for(i in 1:num.perm){
      # Draw two random samples from the global pool
      permuted_assemblage1 <- sample(global_pool, size = length(assemblage1), replace = TRUE)
      permuted_assemblage2 <- sample(global_pool, size = length(assemblage2), replace = TRUE)
      permuted_br[i] <- brainard_robinson(permuted_assemblage1, permuted_assemblage2)
    }

    # Calculate the p-value as the proportion of permuted BR coefficients
    # that are less than or equal to the observed BR coefficient
    p_value <- mean(permuted_br <= observed_br) # notice: smaller or equal as per https://www.mattpeeples.net/BR.html
    return(p_value)
  }

  num_assemblages <- nrow(matrix)
  br_matrix <- matrix(0, nrow=num_assemblages, ncol=num_assemblages)
  p_matrix <- matrix(0, nrow=num_assemblages, ncol=num_assemblages)
  p_class_matrix <- matrix("-", nrow=num_assemblages, ncol=num_assemblages)

  for(i in 1:(num_assemblages-1)){
    for(j in (i+1):num_assemblages){
      br_matrix[i,j] <- brainard_robinson(matrix[i,], matrix[j,])
      br_matrix[j,i] <- br_matrix[i,j]

      p_value <- perm_test_br(matrix[i,], matrix[j,], num.perm)
      p_matrix[i,j] <- p_value
      p_matrix[j,i] <- p_value

      # Classify the p-values
      if(p_value < 0.001) {
        p_class <- "<0.001"
      } else if(p_value < 0.01) {
        p_class <- "<0.01"
      } else if(p_value < 0.05) {
        p_class <- "<0.05"
      } else {
        p_class <- "n.s."
      }

      p_class_matrix[i,j] <- p_class
      p_class_matrix[j,i] <- p_class
    }
  }

  # Set the diagonal elements of br_matrix to 200
  diag(br_matrix) <- 200

  # Attach the input table's row names to the output matrices
  rownames(br_matrix) <- df.row.names
  colnames(br_matrix) <- df.row.names
  rownames(p_matrix) <- df.row.names
  colnames(p_matrix) <- df.row.names
  rownames(p_class_matrix) <- df.row.names
  colnames(p_class_matrix) <- df.row.names

  # Define the color palette to be used later on in the corrplot
  col1 <- colorRampPalette(c("lightblue", "blue", "darkblue"))

  if(clust==TRUE & oneplot==TRUE){
    m <- rbind(c(1,2), c(3,4), c(5,6))
    layout(m)
  }

  if(clust==TRUE){
    dist.matrix <- stats::as.dist(200 - br_matrix, diag=TRUE, upper=TRUE)

    fit <- stats::hclust(dist.matrix, method = aggl.meth)

    # Calculate the max number of clusters to be later used in the loop
    # which iterates the calculation of the average silhouette width
    rd <- dim(matrix)[1]
    max.ncl <- rd-1

    # Create a slot to store the values of the average silhouette value at increasing numbers of cluster solutions
    # this will be used inside the loop
    sil.width.val <- numeric(max.ncl - 1)

    # Create a slot to store the increasing number of clusters whose silhouette is iteratively computed
    sil.width.step <- c(2:max.ncl)

    # Calculate the average silhouette width at each increasing cluster solution
    for (i in min(sil.width.step):max(sil.width.step)) {
      counter <- i - 1
      clust.sol <- cluster::silhouette(cutree(fit, k = i), dist.matrix)
      sil.width.val[counter] <- mean(clust.sol[, 3])
    }

    # Create a dataframe that stores the number of cluster solutions and the corresponding average silhouette values
    sil.res <- as.data.frame(cbind(sil.width.step, sil.width.val))

    # If the user does not enter the desired partition, the latter is equal to the optimal one,
    # ie the one with the largest average silhouette; otherwise use the user-defined partition
    ifelse(is.null(part)==TRUE,
           select.clst.num <- sil.res$sil.width.step[sil.res$sil.width.val == max(sil.res$sil.width.val)],
           select.clst.num <- part)

    # Create the final silhouette data using the partition established at the preceding step
    final.sil.data <- cluster::silhouette(stats::cutree(fit, k = select.clst.num), dist.matrix)

    # Add the cluster membership to a copy of the input dataframe
    data.w.cluster.membership <- df
    data.w.cluster.membership$clust<- assignCluster(df, df, cutree(fit, k = select.clst.num))

    # Aggregate the rows of the input dataframe by summing by cluster membership
    # x[-ncol(x)] tells R to work on all the dataframe excluding the last columnn since, as per the preceding step, it contains the
    # cluster membership
    sum.by.clust <- stats::aggregate(data.w.cluster.membership[-ncol(data.w.cluster.membership)], list(data.w.cluster.membership$clust), sum)[,-1]

    # Number of rows of the preceding dataframe
    nrows <- nrow(sum.by.clust)

    # Add the columns total as a new row of the preceding dataframe
    sum.by.clust[nrows+1,] <- colSums(sum.by.clust, dims=1)

    # Turn the counts to percentages
    prop.by.clust <- as.data.frame(round(prop.table(as.matrix(sum.by.clust), 1)*100,3))

    # Add row names to the above dataframe:
    # add reference to cluster membership to all the rows but the last one
    rownames(prop.by.clust)[seq(1,nrows,1)] <- paste0("cluster ", seq(1,select.clst.num, 1))

    # Add reference to the average percentage as name of the last row of the dataframe
    rownames(prop.by.clust)[nrows+1] <- "average"

    # Plot the Brainerd-Robinson matrix, upon condition
    if(sort.map==TRUE){
      corrplot::corrplot(br_matrix,
                         method="square",
                         addCoef.col="red",
                         is.corr=FALSE,
                         col.lim = c(0, 200),
                         col = col1(100),
                         tl.col="black",
                         tl.cex=0.8,
                         order="hclust",
                         hclust.method=aggl.meth,
                         addrect = select.clst.num,
                         rect.col="red")
    } else {
      corrplot::corrplot(br_matrix,
                         method="square",
                         addCoef.col="red",
                         is.corr=FALSE,
                         col.lim = c(0, 200),
                         col = col1(100),
                         tl.col="black",
                         tl.cex=0.8)
    }

    # Plot the dendrogram
    graphics::plot(fit, main = "Clusters Dendrogram based on dissimilarity (1-BR sim. coeff.)",
                   sub = paste0("\nAgglomeration method: ", aggl.meth),
                   xlab = "",
                   cex = cex.dndr.lab,
                   cex.main = 0.85,
                   cex.sub = 0.75)

    # Add the rectangles representing the selected partitions
    solution <- stats::rect.hclust(fit, k = select.clst.num, border="black")

    # Modify the row names of the final silhouette dataframe to also store the cluster to which each feature is closest
    rownames(final.sil.data) <- paste(rownames(df), final.sil.data[, 2], sep = "_")

    # Plot the silhouette chart
    graphics::plot(final.sil.data,
                   cex.names = cex.sil.lab,
                   max.strlen = 30,
                   nmax.lab = rd + 1,
                   main = "Silhouette plot")

    # Add a line representing the average silhouette width
    graphics::abline(v = mean(final.sil.data[, 3]), lty = 2)

    # Plot the average silhouette value at each cluster solution
    graphics::plot(sil.res, xlab = "number of clusters",
                   ylab = "average silhouette width",
                   ylim = c(0, 1),
                   xaxt = "n",
                   type = "b",
                   main = "Average silhouette width vs. number of clusters",
                   sub = paste0("values on the y-axis represent the average silhouette width at each cluster solution"),
                   cex.main=0.85,
                   cex.sub = 0.75)

    axis(1, at = 0:max.ncl, cex.axis = 0.7)

    # Draw a black dot corresponding to the selected partition
    graphics::points(x=select.clst.num,
                     y=sil.res[select.clst.num-1,2],
                     pch=20)

    # Draw the label selecting the row from the sil.res dataframe corresponding to the selected partition
    graphics::text(x = select.clst.num,
                   y = sil.res[select.clst.num-1,2],
                   labels= round(sil.res[select.clst.num-1,2],3),
                   cex = 0.65,
                   pos = 3,
                   offset = 1.2,
                   srt = 90)

    # Plot the Cleveland's dotcharts of the proportions by cluster
    graphics::dotchart(as.matrix(prop.by.clust),
                       main="By-cluster proportions",
                       xlab="%",
                       xlim=c(0,100),
                       pt.cex=0.85,
                       cex=cex.dot.plt.lab,
                       pch=20)

    graphics::dotchart(as.matrix(t(prop.by.clust)),
                       main="By-cluster proportions",
                       xlab="%",
                       xlim=c(0,100),
                       pt.cex=0.85,
                       cex=cex.dot.plt.lab,
                       pch=20)

  }

  if(clust==FALSE) {
    corrplot::corrplot(br_matrix,
                       method="square",
                       addCoef.col="red",
                       is.corr=FALSE,
                       col.lim = c(0, 200),
                       col = col1(100),
                       tl.col="black",
                       tl.cex=0.8)
    dist.matrix <- NULL
    sil.res <- NULL
    final.sil.data <- NULL
    prop.by.clust <- NULL
    data.w.cluster.membership <- NULL
  }

  # Rename the columns of the sil.res dataframe to give
  # more meaningful labels before returning it
  colnames(sil.res)[1] <- "cluster solution"
  colnames(sil.res)[2] <- "average silhouette width"

  return (list("BR.similarity.matrix"=br_matrix,
               "P-value.matrix"=p_matrix,
               "classified.P-values.matrix"=p_class_matrix,
               "BR.distance_matrix"=dist.matrix,
               "avr.silh.width.by.n.of.clusters"=sil.res,
               "partition.silh.data"=final.sil.data,
               "data.with.cluster.membership"=data.w.cluster.membership,
               "by.cluster.proportion"=prop.by.clust))

}
