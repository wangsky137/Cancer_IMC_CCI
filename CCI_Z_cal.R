
#' Calculate CCI Score for Cell-Cell Interactions
#'
#' The `CCI_score_Cal` function calculates a Cell-Cell Interaction (CCI) score matrix, representing 
#' the interaction strength between different cell types in a spatial context. The function uses 
#' k-nearest neighbors (KNN) and permutation testing to generate a z-score matrix based on 
#' observed and randomized cell type distributions.
#'
#' @param CCI_data A data frame or matrix containing spatial coordinates of cells. The first two columns should be the x and y coordinates.
#' @param CCI_cell A vector of cell type labels that will be calcualted for CCI score to the cells in `CCI_data`.
#' @param cancer_cell_types A character vector specifying the names of cancer cell types in `CCI_cell`.
#' @param Nknn Integer specifying the number of nearest neighbors to consider (default is 10).
#' @param nPermu Integer specifying the number of permutations for the null distribution (default is 100).
#' 
#' @return A matrix containing z-scores that quantify the interaction strength between each pair of cell types.
#' 
#' @details
#' The function calculates cell-cell distances using Euclidean distance and counts the number of nearest neighbors 
#' for each cell. A z-score is computed for each pair of cell types based on a permutation test, where non-cancer 
#' cell types are randomly shuffled to generate a null distribution.
#' 
#'
#' @export



CCI_score_Cal<-function(CCI_data, CCI_cell=NULL, cancer_cell_types=NULL,Nknn=10,nPermu=100) {

set.seed(123)
knn = Nknn
cel.typ<- CCI_cell
myImg= CCI_data

if (length(cel.typ)==0){
  
  cel.typ<- setdiff(data$cel.typ, cancer_cell_types)
}
  #cat("\r",k)
  data = myImg
  ## calculating cell cell distance 
  dis = dist(data[,1:2], method="euclidean", upper = TRUE)
  dis = as.matrix(dis)
  ## counting knn
  for(i in 1:ncol(dis))
  {
    vec = dis[,i]
    vec[i] = Inf
    tmp = rank(vec, ties="random")
    xx = rep(0, length(tmp))
    xx[tmp<=knn] = 1
    dis[,i] = xx
  }
  
  #--------	
  obs.knn = matrix(0, length(cel.typ), length(cel.typ))
  row.names(obs.knn) = colnames(obs.knn) = cel.typ
  for(i in 1:length(cel.typ))
  {
    se = which(data$cel.typ%in%cel.typ[i])
    if(length(se)==0)
    {
      obs.knn[i,] = 0
      next
    }
    tmp = dis[,se, drop=F]
    for(j in 1:length(cel.typ))
    {
      se = which(data$cel.typ%in%cel.typ[j])
      obs.knn[i,j] = sum(tmp[se,])
    }
  }

  pmu.knn1 = pmu.knn2 = matrix(0, length(cel.typ), length(cel.typ))
  pnn = nPermu
  for(p in 1:pnn)
  {
    pm.cel.typ = data$cel.typ
    se = which( ! pm.cel.typ %in% cancer_cell_types)
    pm.cel.typ[se] = sample(pm.cel.typ[se])
    
    for(i in 1:length(cel.typ))
    {
      se = which(pm.cel.typ%in%cel.typ[i])
      if(length(se)==0)
      {
        pmu.knn1[i,] = pmu.knn1[i,] + 0
        pmu.knn2[i,] = pmu.knn2[i,] + 0
        next
      }
      tmp = dis[,se, drop=F]
      for(j in 1:length(cel.typ))
      {
        se = which(pm.cel.typ%in%cel.typ[j])
        xx = sum(tmp[se,])
        pmu.knn1[i,j] = pmu.knn1[i,j] + xx
        pmu.knn2[i,j] = pmu.knn2[i,j] + xx^2
      }
    }
    pmu.avg = pmu.knn1/pnn
    pmu.var = (pmu.knn2 - pmu.knn1^2/pnn)/(pnn-1)
    pmu.std = sqrt(pmu.var)
  }
  zscore = (obs.knn-pmu.avg)/pmu.std
  CCI.img = zscore

  return (CCI.img)

}



