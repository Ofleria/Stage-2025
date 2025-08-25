
# X : données
# seuil : seuil de variance expliquée
# s : décalage temporel
# l : nombre de composante
# max_iter : nombre d'itération maximal
# eps :
# nb_essai : plusieurs essais pour éviter les maximums locaux

DiPCA = function(X, s=1, seuil=0.95, max_iter=100, eps=10^(-8), nb_essai = 10){
  # taille des données
  n = nrow(X)
  p = ncol(X)
  N = n-s # taille des sous matrice de décalage
  # initialisation
  T = matrix(nrow=n, ncol=p) # matrice des scores selon les composantes
  P = matrix(nrow=p, ncol=p) # matrice des loadings
  best_W = matrix(nrow=p, ncol=p) # matrice des meilleurs poids
  for (j in 1:p) {
    # matrices des résultats des essais
    test_W = matrix(nrow = p, ncol = nb_essai)
    test_J = rep(0,nb_essai)
    # reproduction multiple pour éviter un maximum local
    for (essai in 1:nb_essai) {
      # initialisation aléatoire du vecteur de poids W 
      W = rnorm(p)
      # vecteur unitaire (||W|| = 1)
      W = W/(sqrt(sum(W^2)))
      # extraction des variables latentes 
      for (i in 1:max_iter) {
        # score latent de la composante j
        t_score = X %*% W
        # score des sous matrices de décalage
        ti = matrix(nrow=N, ncol=(s+1))
        # calcul des scores des sous matrices
        for (k in 1:(s+1)) {
          ti[,k] = X[k:(N+k-1),] %*% W
        }
        # calcule des beta
        b = t(ti[,1:s]) %*% ti[,(s+1)]
        b = b/(sqrt(sum(b^2))) # vecteur unitaire (||b|| = 1)
        # calcul des poids
        w = rep(0,p)
        for (k in 1:s) {
          w = w + b[k] * (t(X[(s+1):n,]) %*% ti[,k] + t(X[k:(N+k-1),]) %*% ti[,(s+1)])
        }
        w = w/(sqrt(sum(w^2))) # vecteur unitaire (||W|| = 1)
        # vérification de la condition de convergence
        if (sqrt(sum((W-w)^2))<eps) break
        # mise à jour des poids
        W = w
      }
      # calcule de J, la valeur à maximiser
      test_J[essai] = sum(t(b) %*% t(ti[,1:s]) %*% ti[,(s+1)])
      test_W[,essai] = W
    }
    # choix optimal qui maximise la fonction à optimiser
    best = which.max(test_J)
    t_score = X %*% test_W[, best]
    # loadings
    Pj = (t(X) %*% t_score)/as.numeric(t(t_score) %*% t_score)
    # déflation
    X = X - t_score %*% t(Pj)
    # ajout des scores et loadings de la composante
    T[,j] = t_score
    P[,j] = Pj
    best_W[,j] = test_W[, best]
  }
  # choix du nombre de composantes
  val_prp = eigen(cov(T))$values
  var_exp = cumsum(val_prp)/sum(val_prp)
  l = which(var_exp >= seuil)[1]
  # sélection selon le meilleur nombre de composantes
  T = T[,1:l]
  P = P[,1:l]
  best_W = best_W[,1:l]
  R = best_W %*% solve(t(P)%*%best_W)
  # modélisation dynamique: VAR
  # matrice de décalage des scores
  Ts = T[1:N,]
  if (s>1){
    for (k in 2:s) {
      Ts = cbind(Ts, T[k:(N+k-1),])
    }
  }
  # estimateur des moindres carrées
  Theta = solve(t(Ts)%*%Ts) %*% t(Ts)%*%T[(s+1):n,]
  # sortie : poids, loadings, scores, 
  return(list(W=best_W, P=P, T=T, Theta=Theta, R=R))
}


######################################################
# recherche du s optimal sur les données de validation
# data : données d'entraînement
# X_val : données de validation
# seuil : seuil de variance expliquée pour le choix du nombre de composante de DiPCA
# nb : nombre de répétitions
search_s = function(data, X_val, s_values, seuil=0.95, nb=10){
  # nombre de décalage testé
  s_max = length(s_values)
  # taille des données de validation
  n = nrow(X_val)
  # initialisation des vecteurs des résultats
  Aic_score = rep(0,s_max)
  Bic_score = rep(0,s_max)
  # parcourt des différents s testés
  for (s in s_values) {
    # résultats de DiPCA
    results = DiPCA(data, s=s, seuil=seuil, nb_essai = nb)
    # nombre de composantes
    l = ncol(results$T)
    # taille des sous matrices de décalage
    N = n-s
    # scores des données de validation
    T_val = X_val %*% results$R
    # matrice de décalage des scores des données de validation
    Ts_val = T_val[1:N,]
    if (s>1){
      for (k in 2:s) {
        Ts_val = cbind(Ts_val, T_val[k:(N+k-1),])
      }
    }
    # prédiction des scores à s+1
    Ts1_val = Ts_val %*% results$Theta
    # erreur de prédiction
    E_val = X_val[(s+1):n,] - Ts1_val%*%t(results$P) 
    
    sig = (t(E_val) %*% E_val)/N
    # AIC et BIC pénalisé par le nombre de décalage
    Aic_score[s] = N * log(sum(diag(sig))) + l*s*2
    Bic_score[s] = N * log(sum(diag(sig))) + l*s*log(N)
  }
  return(list(Aic_score=Aic_score, Bic_score=Bic_score))
}

