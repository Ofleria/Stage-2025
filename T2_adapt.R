
# data : données
# l : lag
# lam : 
#################################################
############# T2 générale dynamique #############
#################################################
GT2 = function(data,l){
  p = ncol(data) # nombre de variables
  n = nrow(data) # nombre d'observations
  N = n-l+1 # nombre de données surveillées
  T2 = rep(NA,N) # initialisation de T2
  Vt = matrix(nrow = N, ncol = p*l) # matrice des données décalées
  for (i in 1:N){
    # décalage par variable
    for (j in 1:p) Vt[i,((j-1)*l+1):(j*l)] = data[(i+l-1):i,j]
  }
  S = cov(Vt) # covariance des décalages
  # T² des décalages
  for (i in 1:N) T2[i]=t(Vt[i,])%*%solve(S)%*%(Vt[i,])
  return(T2)
}

###########################################################
############# T² adaptative avec lag constant #############
###########################################################
T2_adapt = function(data,l){
  p = ncol(data) # nombre de variables
  n = nrow(data) # nombre d'observations
  N = n-l+1 # nombre de données surveillées
  T2 = rep(NA,N) # initialisation de T2
  Vt = matrix(nrow = N, ncol = p*l) # matrice des données décalées
  Mt = matrix(nrow = N, ncol = p*l) # matrice des moyennes 
  for (i in 1:N){
    # décalage par variable
    for (j in 1:p) Vt[i,((j-1)*l+1):(j*l)] = data[(i+l-1):i,j]
    # moyenne glissante des colonnes
    if (i==1) Mt[i,] = Vt[1,]
    else Mt[i,] = colMeans(Vt[1:i,])
  }
  S = cov(Vt) # covariance des décalages
  # adaptative T²
  for (i in 1:N) T2[i]=t(Mt[i,])%*%solve(S)%*%(Vt[i,]) - t(Mt[i,])%*%solve(S)%*%(Mt[i,])/2
  return(T2)
}

###########################################################
############# T² adaptative avec lag constant #############
######### Calcul récursif de la moyenne par EWMA ##########
###########################################################
T2_adapt_EWMA = function(data,l,lam){
  p = ncol(data) # nombre de variables
  n = nrow(data) # nombre d'observations
  N = n-l+1 # nombre de données surveillées
  T2 = rep(NA,N) # initialisation de T2
  Vt = matrix(nrow = N, ncol = p*l) # matrice des données décalées
  Mt = matrix(nrow = N, ncol = p*l) # matrice des moyennes
  Dt = matrix(nrow = n, ncol = p) # matrice des EWMA
  # prédiction par EWMA
  for (i in 1:n) {
    if (i==1) Dt[i,] = lam*data[1,] + (1-lam)*colMeans(data)
    else Dt[i,] = lam*data[i,] + (1-lam)*Dt[i-1,]
  }
  for (i in 1:N){
    for (j in 1:p){
      # décalage par variable
      Vt[i,((j-1)*l+1):(j*l)] = data[(i+l-1):i,j]
      # décalage des résultats de EWMA
      Mt[i,((j-1)*l+1):(j*l)] = Dt[(i+l-1):i,j]
    }
  }
  S = cov(Vt)# covariance des décalages
  # adaptative T²
  for (i in 1:N) T2[i]=t(Mt[i,])%*%solve(S)%*%(Vt[i,]) - t(Mt[i,])%*%solve(S)%*%(Mt[i,])/2
  return(T2)
}

#########################################
############# T² adaptative #############
#########################################
T2_adapt_2 = function(data,l){
  p = ncol(data) # nombre de variables
  n = nrow(data) # nombre d'observations
  N = n-max(l)+1 # nombre de données surveillées
  T2 = rep(NA,N) # initialisation de T2
  Vt = matrix(nrow = N, ncol = sum(l)) # matrice des données décalées
  Mt = matrix(nrow = N, ncol = sum(l)) # matrice des moyennes 
  for (i in 1:N){
    # décalage par variable
    col = 0
    k = max(l) + i-1
    for (j in 1:p) {
      col = col + l[j]
      Vt[i,(col-l[j]+1):col] = data[k:(k-l[j]+1),j]
    }
    # moyenne glissante des colonnes
    if (i==1) Mt[i,] = Vt[1,]
    else Mt[i,] = colMeans(Vt[1:i,])
  }
  S = cov(Vt) # covariance des décalages
  # adaptative T²
  for (i in 1:N) T2[i]=t(Mt[i,])%*%solve(S)%*%(Vt[i,]) - t(Mt[i,])%*%solve(S)%*%(Mt[i,])/2
  return(T2)
}

###########################################################
##################### T² adaptative  ######################
######### Calcul récursif de la moyenne par EWMA ##########
###########################################################
T2_adapt_EWMA2 = function(data,l,lam){
  p = ncol(data) # nombre de variables
  n = nrow(data) # nombre d'observations
  N = n-max(l)+1 # nombre de données surveillées
  T2 = rep(NA,N) # initialisation de T2
  Vt = matrix(nrow = N, ncol = sum(l)) # matrice des données décalées
  Mt = matrix(nrow = N, ncol = sum(l)) # matrice des moyennes
  Dt = matrix(nrow = n, ncol = p) # matrice des EWMA
  # prédiction par EWMA
  for (i in 1:n) {
    if (i==1) Dt[i,] = lam*data[1,] + (1-lam)*colMeans(data)
    else Dt[i,] = lam*data[i,] + (1-lam)*Dt[i-1,]
  }
  for (i in 1:N){
    # décalage par variable
    col = 0
    k = max(l) + i-1
    for (j in 1:p) {
      col = col + l[j]
      Vt[i,(col-l[j]+1):col] = data[k:(k-l[j]+1),j]
      # décalage des résultats de EWMA
      Mt[i,(col-l[j]+1):col] = Dt[k:(k-l[j]+1),j]
    }
  }
  S = cov(Vt)# covariance des décalages
  # adaptative T²
  for (i in 1:N) T2[i]=t(Mt[i,])%*%solve(S)%*%(Vt[i,]) - t(Mt[i,])%*%solve(S)%*%(Mt[i,])/2
  return(T2)
}
