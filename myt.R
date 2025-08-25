# data : données
# Var : variables étudiées
# obs : observation étudiée
# n : nombre d'observation en phase I

# fonction pour le calcul du T2 pour une observation
MYT_T2_2D <- function(data, Var, obs, n){
  # sélection des données selon les variables étudiées
  df = data[,Var]
  # nombre de variables
  p = length(Var)
  if (p==1){
    # contribution individuelle
    # moy : moyenne de la variable
    moy = mean(df)
    T2 = (df[obs]-moy)^2 / var(df)
  }else {
    # contribution multivariée
    # moyenne de chaque variable
    moy = colMeans(df)
    # covariance
    S = cov(df)
    T2 = t(df[obs,]-moy)%*%solve(S)%*%(df[obs,]-moy)
  }
  # retourne le T2 des données avec les variables considérées
  return(as.numeric(T2))
}

# alpha : niveau de significativité

MYT = function(data, obs, n, alpha = 0.01){
  # initialisation du tableau des résultats
  result = data.frame(vars=character(), T2=numeric(), LSC=numeric(),
                      exclu=character(), stringsAsFactors = FALSE)
  # toutes les variables présentes dans les données
  var_rest = colnames(data)
  # initialisation du numéro de l'étape
  p = 1
  repeat {
    # nombre de variables
    m = length(var_rest)
    # arrêt s'il n'y a pas assez de variables pour la combinaison
    if (m<p) break
    # combinaison de variables
    variables = combn(var_rest, p)
    # nombre de combinaisons
    n_comb = ncol(variables)
    # initialisation pour chaque combinaison
    T2 = rep(0,n_comb)
    # calcul du T2 pour chaque combinaison 
    for (i in 1:n_comb) T2[i] = MYT_T2_2D(data, variables[,i], obs, n)
    # calcul de la limite de contrôle
    LSC = p*(n+1)*(n-1)*qf(1-alpha,p,n-p)/(n*(n-p))
    # décision prise pour l'exclusion
    decision = ifelse(T2 > LSC, "oui", "non")
    # liste des variables ou combinaison de variables prise en compte
    list_var = apply(variables, 2, paste, collapse=" | ")
    # création d'une dataframe avec les résultats
    res = data.frame(vars=list_var, T2=T2, LSC=rep(LSC,n_comb), exclu=decision)
    # ajout à la dataframe contenant tous les résultats
    result = rbind(result, res)
    # les variables ou combinaisons de variables exclues
    exclu = which(T2 > LSC)
    # arrêt si aucune variable ou combinaison n'est exclue
    if (length(exclu)==0 && p>1) break
    # les variables présentes dans les combinaisons exclues
    var_exclu = unique(unlist(strsplit(variables[, exclu], ' ')))
    # les variables restantes
    var_rest = setdiff(var_rest, var_exclu)
    # étape suivante
    p = p+1
  }
  return(result)
}

#########################################################
############## version pour adaptative T² ###############
#########################################################
# data : données
# nvar : nombre de variables étudiées
# obs : observation étudiée
# n : nombre d'observation en phase I
# l : lag de chaque variable
# lam : lambda, paramètre de lissage
MYT_T2_adapt <- function(data, obs, l, lam, nvar){
  if (nvar==1){
    n = length(data) # nombre d'observations
    N = n-max(l)+1 # nombre de données surveillées
    Vt = matrix(nrow = N, ncol = l) # matrice des données décalées
    Mt = matrix(nrow = N, ncol = l) # matrice des moyennes
    Dt = matrix(nrow = n, ncol = nvar) # matrice des EWMA
    for (i in 1:n) {
      if (i==1) Dt[i,] = lam*data[1] + (1-lam)*mean(data)
      else Dt[i,] = lam*data[i] + (1-lam)*Dt[i-1,]
    }
    for (i in 1:N){
      # décalage par variable
      k = l + i-1
      Vt[i,] = data[k:i]
      # décalage des résultats de EWMA
      Mt[i,] = Dt[k:i]
    }
    # covariance des décalages
    S = cov(Vt)
    # adaptative T²
    T2=t(Mt[obs,])%*%solve(S)%*%(Vt[obs,]) - t(Mt[obs,])%*%solve(S)%*%(Mt[obs,])/2
    
  }  else {
    # contribution multivariée
    n = nrow(data) # nombre d'observations
    N = n-max(l)+1 # nombre de données surveillées
    Vt = matrix(nrow = N, ncol = sum(l)) # matrice des données décalées
    Mt = matrix(nrow = N, ncol = sum(l)) # matrice des moyennes
    Dt = matrix(nrow = n, ncol = nvar) # matrice des EWMA
    # prédiction par EWMA
    for (i in 1:n) {
      if (i==1) Dt[i,] = lam*data[1,] + (1-lam)*colMeans(data)
      else Dt[i,] = lam*data[i,] + (1-lam)*Dt[i-1,]
    }
    for (i in 1:N){
      # décalage par variable
      col = 0
      k = max(l) + i-1
      for (j in 1:nvar) {
        col = col + l[j]
        Vt[i,(col-l[j]+1):col] = data[k:(k-l[j]+1),j]
        # décalage des résultats de EWMA
        Mt[i,(col-l[j]+1):col] = Dt[k:(k-l[j]+1),j]
      }
    }
    # covariance des décalages
    S = cov(Vt)
    # adaptative T²
    T2=t(Mt[obs,])%*%solve(S)%*%(Vt[obs,]) - t(Mt[obs,])%*%solve(S)%*%(Mt[obs,])/2
  }
  # retourne le T2 des données avec les variables considérées
  return(as.numeric(T2))
}

MYT_adapt = function(data, obs, n, l, lam, alpha = 0.05){
  # initialisation du tableau des résultats
  data = X_val2
  result = data.frame(vars=character(), T2=numeric(), LSC=numeric(),
                      exclu=character(), stringsAsFactors = FALSE)
  # toutes les variables présentes dans les données
  var_rest = colnames(data)
  # initialisation du numéro de l'étape
  p = 1
  repeat {
    # nombre de variables
    m = length(var_rest)
    # arrêt s'il n'y a pas assez de variables pour la combinaison
    if (m<p) break
    # combinaison de variables
    variables = combn(var_rest, p)
    # nombre de combinaisons
    n_comb = ncol(variables)
    # initialisation pour chaque combinaison
    T2 = rep(0,n_comb)
    LSC = rep(0,n_comb)
    # calcul du T2 pour chaque combinaison 
    for (i in 1:n_comb) {
      df = data[, variables[,i]]
      l1 = l[, variables[,i]]
      T2[i] = MYT_T2_adapt(df, obs, l1, lam, p)
      # calcul de la limite de contrôle
      N = n-max(l1)+1
      p1 = sum(l1)
      # limite de contrôle avec les paramètres de 
      LSC[i] = p1*(N+1)*(N-1)*qf(1-alpha,p1,N-p1)/(N*(N-p1))
    }
    # décision prise pour l'exclusion
    decision = ifelse(T2 > LSC, "oui", "non")
    # liste des variables ou combinaison de variables prise en compte
    list_var = apply(variables, 2, paste, collapse=" | ")
    # création d'une dataframe avec les résultats
    res = data.frame(vars=list_var, T2=T2, LSC=LSC, exclu=decision)
    # ajout à la dataframe contenant tous les résultats
    result = rbind(result, res)
    # les variables ou combinaisons de variables exclues
    exclu = which(T2 > LSC)
    # arrêt si aucune variable ou combinaison n'est exclue
    if (length(exclu)==0 && p>1) break
    # les variables présentes dans les combinaisons exclues
    var_exclu = unique(unlist(strsplit(variables[, exclu], ' ')))
    # les variables restantes
    var_rest = setdiff(var_rest, var_exclu)
    # étape suivante
    p = p+1
  }
  return(result)
}
