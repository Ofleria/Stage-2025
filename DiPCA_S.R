

library(multiSPC)
library(ggplot2)
library(spc)
library(tseries)
library(forcats)
library(forecast)
library(zoo)

# X : données
# seuil : seuil de variance expliquée
# s : décalage temporel
# l : nombre de composante
# nb_essai : plusieurs essais pour éviter les maximums locaux
# S : décalage saisonnier
# t : saisonnalité

DiPCA_S = function(X, s=1, S=0, t=0, seuil=0.95, max_iter=100, eps=10^(-8), nb_essai = 10){
  n = nrow(X)
  p = ncol(X)
  N = n-(S*t) # taille des sous matrice de décalage
  # initialisation
  T = matrix(nrow=n, ncol=p) # matrice des scores selon les composantes
  P = matrix(nrow=p, ncol=p) # matrice des loadings
  best_W = matrix(nrow=p, ncol=p) # matrice des poids
  for (j in 1:p) {
    # variables des résultats des essais
    test_W = matrix(nrow = p, ncol = nb_essai)
    test_J = rep(0,nb_essai)
    # de multiple essais pour éviter un maximun local
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
        ti = matrix(nrow=N, ncol=(S+s+1))
        # calcul des scores des sous matrices
        if (S!=0) {
          # saisonnier
          for (k in 1:S) {
            r = (k-1)*t+1
            ti[,k] = X[r:(N+r-1),] %*% W
          }
        }
        for (k in 1:(s+1)) {
          r = S*t-s+k
          ti[,k+S] = X[r:(N+r-1),] %*% W
        }
        
        b = t(ti[,1:(s+S)]) %*% ti[,(s+S+1)]
        b = b/(sqrt(sum(b^2))) # vecteur unitaire (||b|| = 1)
        # calcul des nouveaux poids
        w = rep(0,p)
        if (S!=0){
          # saisonnier
          for (k in 1:S) {
            r = (k-1)*t+1
            w = w + b[k] * (t(X[(S*t+1):n,]) %*% ti[,k] + t(X[r:(N+r-1),]) %*% ti[,(S+s+1)])
          }
        }
        for (k in 1:s) {
          r = S*t-s+k
          w = w + b[k+S] * (t(X[(S*t+1):n,]) %*% ti[,k+S] + t(X[r:(N+r-1),]) %*% ti[,(S+s+1)])
        }
        w = w/(sqrt(sum(w^2))) # vecteur unitaire (||W|| = 1)
        # condition de convergence
        if (sqrt(sum((W-w)^2))<eps) break
        # mise à jour des poids
        W = w
      }
      # calcule de J
      test_J[essai] = sum(t(b) %*% t(ti[,1:(S+s)]) %*% ti[,(S+s+1)])
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
  #  
  Ts = T[1:N,]
  if (S!=0){
    # saisonnier
    if (S>1){
      for (k in 2:S) {
        r = (k-1)*t+1
        Ts = cbind(Ts, T[r:(N+r-1),])
      }
    }
    for (k in 1:s) {
      r = S*t-s+k
      Ts = cbind(Ts, T[r:(N+r-1),])
    }
  } else {
    if (s>1){
      for (k in 2:s) {
        Ts = cbind(Ts, T[k:(N+k-1),])
      }
    }
  }
  # estimateur des moindres carrées
  Theta = solve(t(Ts)%*%Ts) %*% t(Ts)%*%T[(S*t+1):n,]
  
  return(list(W=best_W, P=P, T=T, Theta=Theta, R=R))
}


search_s = function(data, X_val, s_values, S_values, t, seuil=0.95, nb=5){
  s_max = length(s_values)
  S_max = length(S_values)
  Aic_score = matrix(nrow=s_max, ncol=S_max)
  Bic_score = matrix(nrow=s_max, ncol=S_max)
  for (S in S_values) {
    for (s in s_values) {
      results = DiPCA_S(data, s=s, S=S, t=t, seuil=seuil, nb_essai = nb)
      l = ncol(results$T)
      n2 = nrow(X_val)
      N2 = n2 - (S*t)
      # scores latent dynamique
      T_val = X_val %*% results$R
      #
      Ts_val = T_val[1:N2,]
      if (S!=0){
        # saisonnier
        if (S>1){
          for (k in 2:S) {
            r = (k-1)*t+1
            Ts_val = cbind(Ts_val, T_val[r:(N2+r-1),])
          }
        }
        for (k in 1:s) {
          r = S*t-s+k
          Ts_val = cbind(Ts_val, T_val[r:(N2+r-1),])
        }
      } else {
        if (s>1){
          for (k in 2:s) {
            Ts_val = cbind(Ts_val, T_val[k:(N2+k-1),])    }
        }
      }
      # score prédit
      Ts1_val = Ts_val %*% results$Theta
      
      # erreur de prédiction
      E_val = X_val[(S*t+1):n2,] - Ts1_val%*%t(results$P) 
      
      sig = (t(E_val) %*% E_val)/N2
      #Aic_score[s] = N * log(sum(diag(sig))) + l*2
      Aic_score[s,S] = N2 * log(sum(diag(sig))) + l*s*S*2
      
      #Bic_score[s] = N * log(sum(diag(sig))) + l*log(N)
      Bic_score[s,S] = N2 * log(sum(diag(sig))) + l*s*S*log(N2)
    }
  }
  return(list(Aic=Aic_score, Bic=Bic_score))
}

