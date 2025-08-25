
# X : données standardisées
# d : nombre de décalage temporelle
# r : ordre de la SVD

FDHODMD = function(X, d, r){
  n = ncol(X) # nombre d'observations
  p = nrow(X) # nombre de variables
  # étape 1
  Xd1 = X[, (d+1):(n-1)]
  XXd = X[, 1:(n-d-1)]
  if (d>1){
    for (i in 2:d) {
      XXd = rbind(XXd, X[, i:(n+i-d-2)])
    }
  }
  
  A_ho = Xd1 %*% ginv(XXd)
  
  # étape 2
  svd_result = svd(X, nu=r, nv=r)
  U = svd_result$u
  V = svd_result$v
  sig = diag(svd_result$d[1:r])
  
  A_ho_til = t(U) %*% A_ho[, 1:p] %*% U
  if (d>1){
    for (k in 2:d) {
      A_ho_til = cbind(A_ho_til, (t(U) %*% A_ho[, ((k-1)*p+1):(k*p)] %*% U))
    }
  }
  #
  return(A_ho_til)
}


############################################################
# recherche de la dimension des caractéristiques dynamiques
# X : données
# p : nombre de variables
# d : nombre de lag
detect_r = function(X,p,d){
  n = ncol(X) # nombre d'observations
  # initialisation de vecteur des résultats
  p_val = rep(0,p)
  for (r in 1:p) {
    # SVD des données
    svd_result = svd(X, nu=r, nv=r)
    U = svd_result$u
    V = svd_result$v
    # application de FDHODMD
    A_ho_til = FDHODMD(X, d, r)
    # les caractéristiques dynamiques
    X_til = t(U) %*% X
    # prédictions des caractéristiques dynamiques
    Xkd_pred = A_ho_til[,1:r] %*% X_til[,1:(n-d)]
    if (d>1){
      for (i in 2:d) {
        Xkd_pred = Xkd_pred + A_ho_til[,((i-1)*r+1):(i*r)] %*% X_til[,i:(n-d+i-1)]
      }
    }
    # résidu statique
    E = X[,(d+1):n] - U %*% Xkd_pred
    n1 = ncol(E)
    # application de la méthode avec 1 variable latente
    A_param = FDHODMD(E, d, 1)
    # SVD des résidus
    svd_result = svd(E, nu=1, nv=1)
    U = svd_result$u
    V = svd_result$v
    # caractéristiques statiques
    E_til = t(U) %*% E
    # prédictions
    Ekd_pred = A_param[,1] * E_til[,1:(n1-d)]
    if (d>1){
      for (i in 2:d) {
        Ekd_pred = Ekd_pred + A_param[,i] %*% E_til[,i:(n1-d+i-1)]
      }
    }
    # test de Ljung Box pour l'autocorrélation
    p_val[r] = Box.test(t(Ekd_pred), type = "Ljung-Box")$p.value
  }
  return(p_val)
}


######################################################################
# calcul des différentes statistiques de contrôle
# U : m
# A : opérateur dynamique obtenu après l'entrainement du FDHODMD
# X_cntrl : données pour les statistiques de contrôle
# d : nombre de décalage temporelle
# r : ordre de la SVD
# seuil : seuil de variance expliquée pour le nombre de composantes optimal
stat_cntrl = function(U, A, X_cntrl, d, r, seuil=0.90){
  n = ncol(X_cntrl)
  # caractéristiques dynamiques
  X_til = t(U) %*% X_cntrl
  # prédiction 
  Xkd_pred = A[,1:r] %*% X_til[,1:(n-d)]
  if (d>1){
    for (i in 2:d) {
      Xkd_pred = Xkd_pred + A[,((i-1)*r+1):(i*r)] %*% X_til[,i:(n-d+i-1)]
    }
  }
  # stat de contrôle dynamique
  # sur les prédictions
  Dv = rep(0,(n-d))
  for (i in 1:(n-d)) Dv[i] = t(Xkd_pred[,i]) %*% Xkd_pred[,i]
  # sur les erreurs de prédiction
  E_dyn = X_til[,(d+1):n] - Xkd_pred
  De = rep(0,(n-d))
  for (i in 1:(n-d)) De[i] = norm(E_dyn[,i], type="2")
  # erreur de reconstruction
  Es = X_cntrl[,(d+1):n] - U %*% Xkd_pred
  # stat de contrôle statique
  acp_res = prcomp(t(Es))
  # recherche du nombre de composante optimal
  val_prp = acp_res$sdev^2
  var_exp = cumsum(val_prp)/sum(val_prp)
  l = which(var_exp >= seuil)[1]
  P = acp_res$rotation[,1:l]
  # score de l'acp
  T = acp_res$x[,1:l]
  # T2
  Sv = T2_Hotelling_k1(T,S=diag(val_prp[1:l]),MoyT=rep(0,l))
  # résidus statiques
  Er = Es - P%*%t(T)
  # SPE
  Se = diag(t(Er) %*% Er)
  
  return(list(Dv=Dv, De=De, Sv=Sv, Se=Se))
}

