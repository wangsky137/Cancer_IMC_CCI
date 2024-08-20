

CCI_score_Cal<-function(CCI_data, outfile=NULL, cel.typ, cancer_cell_name,Nknn=10,nPermu=100) {

set.seed(123)
knn = Nknn
cel.typ<- cel.typ
myImg= CCI_data


  cat("\r",k)
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
    se = which(pm.cel.typ!=cancer_cell_name)
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

save(CCI.img, file = outfile)
}



