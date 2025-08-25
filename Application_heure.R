# les librairies nécessaires
library(multiSPC)
library(ggplot2)
library(spc)
library(tseries)
library(forcats)
library(forecast)
library(zoo)
library(MASS)


###################################################################
# graphique de carte de contrôle
plot_chart2 = function(data, dates, LSC = NULL,
                       xlab = "Dates",
                       ylab = "Valeurs",
                       Type = "") {
  
  df <- data.frame(t = dates, values = data)
  # Marquage des points hors contrôle
  df$out_of_control <- rep(FALSE, nrow(df))
  if (!is.null(LSC)) {
    df$out_of_control <- df$out_of_control | (df$values > LSC)
  }
  # Création du plot de base
  p <- ggplot(df, aes(x = t, y = values)) +
    geom_point(aes(color = out_of_control), size = 1) +
    geom_line(color = "darkgray", linewidth = 0.6) +
    scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red"), guide = "none") +
    labs(x = xlab,
         y = ylab,
         title = Type) +
    theme_minimal()
  # Ajouter LSC si présent
  if (!is.null(LSC)) {
    p <- p +
      geom_hline(yintercept = LSC, linetype = "dashed", color = "red") +
      annotate("label", x = dates[2], y = LSC, label = "LSC", color = "red")
  }
  return(p)
}

# limite de contrôle par KDE
library(kde1d)
# limites de contrôle
lim_cntrl = function(x, alpha=0.05){
  fit = kde1d(x)
  return(qkde1d(1-alpha,fit))
}
####################################################################
###################### Importation des données #####################
####################################################################
# dossier contenant les données
setwd("C:/Users/hp/OneDrive/Bureau/Work/ONIRIS/BAM/Codes/données")
# données horaire
data_hour_O2=readRDS('data_hour_O2.RDS')
# sélection des jours pertinents
data = as.matrix(data_hour_O2[(data_hour_O2[,"jour"]<=34 
                               & data_hour_O2[,"jour"]>5),])
# gestion des valeurs manquantes
for (col in colnames(data)) {
  data[,col] = na.approx(data[,col])
}
# séparation train/test
train_hour = data[data[,"jour"]<=12,-32]
test_hour = data[data[,"jour"]>12,-32]
# réinitialisation des numéros de ligne
rownames(train_hour) = NULL
rownames(test_hour) = NULL

# données absorbance bleu vert jaune rouge
df = train_hour[,c(8:11, 14:17, 20:23, 25:28)]
# ACP des absorbance pour créer une nouvelle variable couleur
acp_res = prcomp(df, center=TRUE, scale=TRUE)
summary(acp_res)
# centrer/réduire test avec les résultats de l'ACP sur train
ts = scale(test_hour[,c(8:11, 14:17, 20:23, 25:28)], 
           center=acp_res$center, scale=acp_res$scale)
# résultat ACP de train sur test
ts_red = ts %*% acp_res$rotation
# sélection des variables nécessaires
data_train = train_hour[,-c(4,6, 8:11, 14:17, 20:23, 25:28)]
data_test = test_hour[,-c(4,6, 8:11, 14:17, 20:23, 25:28)]
# centrer/réduire les données d'apprentissage
X_ca2 = scale(data_train)
# moyenne et écart type de train
moy = sapply(as.data.frame(data_train), mean)
ecar = sapply(as.data.frame(data_train), sd)
# centrer/réduire test avec les paramètres de train
X_val2 = as.matrix(mapply(function(col,m,e) (col-m)/e, 
                         as.data.frame(data_test), moy, ecar))
# ajout de la nouvelle variable couleur
X_ca2 = cbind(X_ca2, couleur = acp_res$x[,1])
X_val2 = cbind(X_val2, couleur = ts_red[, 1])


# standardiser les données d'entrainement
X_ca = scale(train_hour)
# standardisation du test avec les paramètres du train
moy = sapply(as.data.frame(train_hour), mean) # moyenne de train
ecar = sapply(as.data.frame(train_hour), sd) # écart type de train
# centrer/réduire test
X_val = as.matrix(mapply(function(col,m,e) (col-m)/e, 
                         as.data.frame(test_hour), moy, ecar))
# taille des tables de données
n1 = nrow(X_ca)
n2 = nrow(X_val)

# Génération d'horodatages horaires
dates = seq.POSIXt(from = as.POSIXct("2022-05-31 00:00"), 
                   by = "hour", length.out = 696)
dates1 = seq.POSIXt(from = as.POSIXct("2022-05-31 00:00"), 
                    by = "hour", length.out = nrow(train_hour))
dates2 = seq.POSIXt(from = as.POSIXct("2022-06-07 00:00"), 
                    by = "hour", length.out = nrow(test_hour))

# train et test ensemble
X_tot = rbind(X_ca, X_val)
X_tot2 = rbind(X_ca2, X_val2)
####################################################################
######################### T² de hotelling ##########################
####################################################################
# train_hour et test_hour : toutes les variables
# X_ca et X_val : données standardisée
# X_ca2 et X_val2 : données standardisée sans absorbance et avec couleur (pc1)

# T² de base sur toutes les variables
lsc = LSC_T2_Hotelling_k1(train_hour)
res = T2_Hotelling_k1(train_hour)
plot_chart2(res, dates1, LSC=lsc[1])
res = T2_Hotelling_k1(test_hour)
plot_chart2(res, dates2, LSC=lsc[2])

###################################################
# T² avec ACP sur toutes les variables
# réduction de dimension
tt <- prcomp(train_hour, center=TRUE, scale=TRUE)
plot(tt, type="l")

summary(tt)
# entraînement
tt_red = tt$x[, 1:3]
res = T2_Hotelling_k1(tt_red)
lsc = LSC_T2_Hotelling_k1(tt_red, c(0.05,0.05))
plot_chart2(res, dates1, LSC=lsc[1], Type="ACP et T2 sur données d'entrainement")
plot_chart2(res, dates1, LSC=lsc[1])
# test
ts = scale(test_hour, center=tt$center, scale=tt$scale)
ts_red = ts %*% tt$rotation
ts_red = ts_red[, 1:3]
res = T2_Hotelling_k1(ts_red)
plot_chart2(res, dates2, LSC=lsc[2], Type="ACP et T2 sur données de test")
plot_chart2(res, dates2, LSC=lsc[2])
###################################################
# sans les absorbances
# avec l'ancienne variable couleur
df1 = train_hour[,-c(4, 8:11, 14:17, 20:23, 25:28)]
df2 = test_hour[,-c(4, 8:11, 14:17, 20:23, 25:28)]
# réduction de dimension
tt <- prcomp(df1, center=TRUE, scale=TRUE)
plot(tt, type="l")
summary(tt)
# entraînement
tt_red = tt$x[, 1:3]
res = T2_Hotelling_k1(tt_red)
lsc = LSC_T2_Hotelling_k1(tt_red, c(0.05,0.05))
plot_chart2(res, dates1, LSC=lsc[1], Type="ACP et T2 sur données d'entrainementsans les absorbances")
plot_chart2(res, dates1, LSC=lsc[1])
# test
ts = scale(df2, center=tt$center, scale=tt$scale)
ts_red = ts %*% tt$rotation
ts_red = ts_red[, 1:3]
res = T2_Hotelling_k1(ts_red)
plot_chart2(res, dates2, LSC=lsc[2], Type="ACP et T2 sur données de test sans les absorbances")
plot_chart2(res, dates2, LSC=lsc[2])
###################################################
# avec la nouvelle variable couleur (ACP)
tt <- prcomp(X_ca2, center=FALSE)
plot(tt, type="l")
summary(tt)
# entraînement
tt_red = tt$x[, 1:3]
res = T2_Hotelling_k1(tt_red)
lsc = LSC_T2_Hotelling_k1(tt_red, c(0.05,0.05))
plot_chart2(res, dates1, LSC=lsc[1], 
            Type="ACP et T2 sur données d'entrainement avec nouvelle variable couleur")
plot_chart2(res, dates1, LSC=lsc[1])
# test
ts_red = X_val2 %*% tt$rotation
ts_red = ts_red[, 1:3]
res = T2_Hotelling_k1(ts_red)
plot_chart2(res, dates2, LSC=lsc[2], 
            Type="ACP et T2 sur données de test avec nouvelle variable couleur")
plot_chart2(res, dates2, LSC=lsc[2])
####################################################################
############################## DiPCA ###############################
####################################################################
source("C:/Users/hp/OneDrive/Bureau/Work/ONIRIS/BAM/Codes/DiPCA.R")
# recherche du meilleur nombre de lag
res = search_s(X_ca2, X_val2, 1:5, seuil=0.95, 5)
# résultat de la recherche du lag optimal
s_aic = which.min(res$Aic_score)
s_bic = which.min(res$Bic_score)
# le choix du nombre de lag
s = s_bic # si BIC
# s = s_aic # si AIC
  
# résultats du DiPCA
results = DiPCA(X_ca2, s=s, seuil=0.95)
# nombre de composantes
l = ncol(results$T)
# nombre d'observations restant 
N1 = n1 - s
N2 = n2 - s
# scores latent dynamique
T_ca = X_ca2 %*% results$R
T_val = X_val2 %*% results$R

# initialisation des matrices de décalage des scores 
Ts_ca = T_ca[1:N1,]
Ts_val = T_val[1:N2,]
# remplissage des matrices
if (s>1){
  for (k in 2:s) {
    Ts_ca = cbind(Ts_ca, T_ca[k:(N1+k-1),])
    Ts_val = cbind(Ts_val, T_val[k:(N2+k-1),])
  }
}
# prédiction des scores
Ts1_ca = Ts_ca %*% results$Theta
Ts1_val = Ts_val %*% results$Theta
# erreur de prédiction des scores
V_k = T_ca[(s+1):n1,] - Ts1_ca
V_test = T_val[(s+1):n2,] - Ts1_val
##################
# carte de contrôle sur les erreurs de prédictions
# train
res = T2_Hotelling_k1(V_k)
lsc = LSC_T2_Hotelling_k1(V_k)
plot_chart2(res, dates1[(s+1):n1], LSC=lsc[1])
# test
res = T2_Hotelling_k1(V_test)
plot_chart2(res,dates2[(s+1):n2], LSC=lsc[2])
###################
# carte de contrôle sur les prédictions
res = T2_Hotelling_k1(Ts1_ca)
lsc = LSC_T2_Hotelling_k1(Ts1_ca)
plot_chart2(res, dates1[(s+1):n1], LSC=lsc[1])
# test
res = T2_Hotelling_k1(Ts1_val)
plot_chart2(res,dates2[(s+1):n2], LSC=lsc[2])
# pour le MYT
T2_dipca = res
lsc_dipca = lsc[2]
s_dipca = s

######################
# résidus dynamique
E_ca = X_ca2[(s+1):n1,] - Ts1_ca%*%t(results$P)
E_val = X_val2[(s+1):n2,] - Ts1_val%*%t(results$P)

# ACP sur les erreurs
acp_ca = prcomp(E_ca)
plot(acp_ca, type="l")
summary(acp_ca)
# scores latent statique
Tr_ca = acp_ca$x

Pr_ca = acp_ca$rotation
# résidus statiques
Er_ca = E_ca - Tr_ca%*%t(Pr_ca)


# carte de contrôle score statique
# T2
# train
res = T2_Hotelling_k1(Tr_ca[,1:4])
lsc = LSC_T2_Hotelling_k1(Tr_ca[,1:4])
plot_chart2(res, dates1[(s+1):n1],LSC=lsc[1])
# test
ts = scale(E_val, center=acp_ca$center, scale=acp_ca$scale)
ts_red = ts %*% acp_ca$rotation
ts_red = ts_red[,1:4]

res = T2_Hotelling_k1(ts_red)
plot_chart2(res, dates2[(s+1):n2],LSC=lsc[2])






####################################################################
######################## DiPCA saisonnier ##########################
####################################################################
source("C:/Users/hp/OneDrive/Bureau/Work/ONIRIS/BAM/Codes/DiPCA_S.R")

# nombre de lag
sea = search_s(X_ca2, X_val2, 1:5, 1:4, 24, seuil=0.95, 7)
# minimiser AIC
S_aic = which.min(sea$Aic)%/%5 +1
s_aic = which.min(sea$Aic)%%5
# Minimiser BIC
S_bic = which.min(sea$Bic)%/%5 +1
s_bic = which.min(sea$Bic)%%5
# choix des paramètres
S = S_bic 
s = s_bic
# Sélection avec de test avec les S jours précédents
test_hour2 = data[data[,"jour"]>(12-S),-32]
rownames(test_hour2) = NULL
# transformation selon les résultats de l'ACP de train
ts = scale(test_hour2[,c(8:11, 14:17, 20:23, 25:28)], 
           center=acp_res$center, scale=acp_res$scale)
ts = ts %*% acp_res$rotation
# sélection des variables nécessaire
data_test2 = test_hour2[,-c(6, 4, 8:11, 14:17, 20:23, 25:28)]
data_train = train_hour[,-c(4,6, 8:11, 14:17, 20:23, 25:28)]
# moyenne et écart type de train
moy = sapply(as.data.frame(data_train), mean)
ecar = sapply(as.data.frame(data_train), sd)
# centrer/réduire test avec les paramètres de train
X_val3 = as.matrix(mapply(function(col,m,e) (col-m)/e, 
                         as.data.frame(data_test2), moy, ecar))
# ajout de la nouvelle variable couleur
X_val3 = cbind(X_val3, couleur = ts[, 1])

# application de DiPCA saisonnier
results = DiPCA_S(X_ca2, s=s, S=S, t=24, seuil=0.95)
n1 = nrow(X_ca2) # nombre d'observation de train
l = ncol(results$T) # nombre de variables latentes
t = 24 # périodicité
N1 = n1-(S*t) # nombre de données surveillées de train
n3 = nrow(X_val3) # nombre d'observation de test
N3 = n3 - (S*t) # nombre de données surveillées de test
n = nrow(X_tot2) # nombre total d'observation 
N = n - (S*t) # nombre total de données surveillées
# scores latent dynamique
T_ca = X_ca2 %*% results$R # train
T_val = X_val3 %*% results$R # test
T_tot = X_tot2 %*% results$R # tout
# modèle de décalage pour le VAR
Ts_ca = T_ca[1:N1,]
Ts_val = T_val[1:N3,]
Ts_tot = T_tot[1:N,]
if (S!=0){
  # saisonnier
  if (S>1){
    for (k in 2:S) {
      r = (k-1)*t+1
      Ts_ca = cbind(Ts_ca, T_ca[r:(N1+r-1),])
      Ts_val = cbind(Ts_val, T_val[r:(N3+r-1),])
      Ts_tot = cbind(Ts_tot, T_tot[r:(N+r-1),])
    }
  }
  # non saisonnier
  for (k in 1:s) {
    r = S*t-s+k
    Ts_ca = cbind(Ts_ca, T_ca[r:(N1+r-1),])
    Ts_val = cbind(Ts_val, T_val[r:(N3+r-1),])
    Ts_tot = cbind(Ts_tot, T_tot[r:(N+r-1),])
  }
} else {
  # aucune saisonnalité
  if (s>1){
    for (k in 2:s) {
      Ts_ca = cbind(Ts_ca, T_ca[k:(N1+k-1),])
      Ts_val = cbind(Ts_val, T_val[k:(N3+k-1),])    
      Ts_tot = cbind(Ts_tot, T_tot[k:(N+k-1),])    
    }
  }
}
# prédiction des scores du temps t+1
Ts1_ca = Ts_ca %*% results$Theta
Ts1_val = Ts_val %*% results$Theta
Ts1_tot = Ts_tot %*% results$Theta
# différence entre vraie donnée et prédiction
V_k = T_ca[(S*t+1):n1,] - Ts1_ca
V_test = T_val[(S*t+1):n3,] - Ts1_val
V_tot = T_tot[(S*t+1):n,] - Ts1_tot


###################
# carte de contrôle sur les prédictions
# train
res = T2_Hotelling_k1(Ts1_ca)
lsc = LSC_T2_Hotelling_k1(Ts1_ca, alpha = c(0.05,0.05))
plot_chart2(res,dates1[(S*t+1):n1], LSC=lsc[1])
# test
res = T2_Hotelling_k1(Ts1_val)
plot_chart2(res, dates2, LSC=lsc[2])
# tout
res = T2_Hotelling_k1(Ts1_tot)
plot_chart2(res, dates[(S*t+1):n], LSC=lsc[2],
            Type="DiPCA saisonnier et T2 sur toutes les données")
plot_chart2(res, dates[(S*t+1):n], LSC=lsc[2])
##################
# carte de contrôle sur les erreurs de prédictions
# train
res = T2_Hotelling_k1(V_k)
lsc = LSC_T2_Hotelling_k1(V_k)
plot_chart2(res, dates1[(S*t+1):n1], LSC=lsc[1])
# test
res = T2_Hotelling_k1(V_test)
plot_chart2(res,dates2, LSC=lsc[2])

res = T2_Hotelling_k1(V_tot)
plot_chart2(res,dates[(S*t+1):n], LSC=lsc[2])

#####################
# résidus dynamique
E_ca = X_ca2[(S*t+1):n1,] - Ts1_ca%*%t(results$P)
E_val = X_val3[(S*t+1):n3,] - Ts1_val%*%t(results$P)
E_tot = X_tot2[(S*t+1):n,] - Ts1_tot%*%t(results$P)
# ACP sur les résidus
acp_ca = prcomp(E_ca)
plot(acp_ca, type="l")
summary(acp_ca)
# scores latent statique
Tr_ca = acp_ca$x

# carte de contrôle score statique
# T2 train
res = T2_Hotelling_k1(Tr_ca[,1:5])
lsc = LSC_T2_Hotelling_k1(Tr_ca[,1:5], alpha = c(0.05,0.05))
plot_chart2(res, dates1[(S*t+1):n1], LSC=lsc[1])
# ACP de test avec les paramètres de l'ACP de train
ts = scale(E_val, center=acp_ca$center, scale=acp_ca$scale)
ts_red = ts %*% acp_ca$rotation
ts_red = ts_red[,1:5]
# T2 test
res = T2_Hotelling_k1(ts_red)
plot_chart2(res, dates2,LSC=lsc[2])
# ACP de toute les données avec les paramètres de l'ACP de train
ts = scale(E_tot, center=acp_ca$center, scale=acp_ca$scale)
ts_red = ts %*% acp_ca$rotation
ts_red = ts_red[,1:5]
# total
res = T2_Hotelling_k1(ts_red)
plot_chart2(res, dates[(S*t+1):n],LSC=lsc[2])

####################################################################
############################## FDHODMD #################################
####################################################################
source("C:/Users/hp/OneDrive/Bureau/Work/ONIRIS/BAM/Codes/FDHODMD.R")
X = t(X_ca2)
n = ncol(X)
p = nrow(X)
d = 2

p_val = detect_r(X,p,d)
round(p_val,3)
which(p_val>0.05)
r = min(which(p_val>0.05))

##############################
# paramètres du modèle
U = svd(X, nu=r, nv=r)$u
# coefficients de VAR
A = FDHODMD(X, d, r)

######################################
# carte sur les données d'entraînement
X_cntrl = t(X_ca2)
# les différentes statistiques de contrôle
res = stat_cntrl (U, A, X_cntrl, d, r, 0.90)

# sur les prédictions
lsc1 = lim_cntrl(res$Dv, 0.01)
plot_chart2(res$Dv, dates1[(d+1):n1], LSC=lsc1)
# sur les erreurs de prédiction
lsc2 = lim_cntrl(res$De, 0.01)
plot_chart2(res$De, dates1[(d+1):n1], LSC=lsc2)
# erreur de reconstruction
lsc3 = lim_cntrl(res$Sv, 0.01)
plot_chart2(res$Sv, dates1[(d+1):n1], LSC=lsc3)
# résidus statiques
lsc4 = lim_cntrl(res$Se, 0.01)
plot_chart2(res$Se, dates1[(d+1):n1], LSC=lsc4)

################################
# carte sur les données de test
X_cntrl = t(X_val2)

res = stat_cntrl (U, A, X_cntrl, d, r, 0.90)
# sur les prédictions
plot_chart2(res$Dv, dates2[(d+1):n2], LSC=lsc1)
# sur les erreurs de prédiction
plot_chart2(res$De, dates2[(d+1):n2], LSC=lsc2)
# erreur de reconstruction
plot_chart2(res$Sv, dates2[(d+1):n2], LSC=lsc3)
# résidus statiques
plot_chart2(res$Se, dates2[(d+1):n2], LSC=lsc4)
plot_chart2(res$Se, dates2[(d+1):n2])


####################################################################
########################### Adaptative #############################
####################################################################
source("C:/Users/hp/OneDrive/Bureau/Work/ONIRIS/BAM/Codes/T2_adapt.R")



##################### GT2 ###########################
# cartes de contrôle 
l = 3
p = ncol(X_ca2) # nombre de variables
n = nrow(X_ca2) # nombre d'observations de train
n2 = nrow(X_val2) # nombre d'observations de test
n3 = nrow(X_tot2) # nombre d'observations au total
N = n-l+1
# limite de contrôle avec les paramètres de la matrice de décalage
lsc = LSC_T2_Hotelling_k1(matrix(nrow = N, ncol = p*l), alpha=c(0.05,0.05))
# GT2 sur train
T2 = GT2(X_ca2,l)
plot_chart2(T2,dates1[l:n], LSC=lsc[1],
            Type="GT2 sur données d'entrainement avec nouvelle variable couleur")
plot_chart2(T2,dates1[l:n], LSC=lsc[1])
# GT2 sur test
T2 = GT2(X_val2,l)
plot_chart2(T2, dates2[l:n2], LSC=lsc[2],
            Type="GT2 sur données de test avec nouvelle variable couleur")
plot_chart2(T2, dates2[l:n2], LSC=lsc[2])
# GT2 sur les données totales
T2 = GT2(X_tot2,l)
plot_chart2(T2, dates[l:n3], LSC=lsc[2],
            Type="GT2 sur toutes les données avec nouvelle variable couleur")
plot_chart2(T2, dates[l:n3], LSC=lsc[2])


#################################################
# calcul de la limite avec KDE
# GT2 sur train
T2 = GT2(X_ca2,l)
lsc = lim_cntrl(T2)
plot_chart2(T2,dates1[l:n], LSC=lsc,
            Type="GT2 sur données d'entrainement avec nouvelle variable couleur")
plot_chart2(T2,dates1[l:n], LSC=lsc)
# GT2 sur les données totales
T2 = GT2(X_tot2,l)
plot_chart2(T2, dates[l:n3], LSC=lsc,
            Type="GT2 sur toutes les données avec nouvelle variable couleur")
plot_chart2(T2, dates[l:n3], LSC=lsc)


##################### adapt ###########################
# cartes de contrôle
l = 2
# limite de contrôle avec les paramètres de train
lsc = LSC_T2_Hotelling_k1(X_ca2, alpha=c(0.05,0.05))
# T2_adapt sur train
T2 = T2_adapt(X_ca2,l)
plot_chart2(T2,dates1[l:n], LSC=lsc[1],
            Type="Adaptative T2 sur données d'entrainement avec lag constant")
plot_chart2(T2,dates1[l:n], LSC=lsc[1])

# T2_adapt sur test
T2 = T2_adapt(X_val2,l)
plot_chart2(T2, dates2[l:n2], LSC=lsc[2],
            Type="Adaptative T2 sur données de test avec lag constant")
plot_chart2(T2, dates2[l:n2], LSC=lsc[2])


# T2_adapt sur train
T2 = T2_adapt(X_ca2,l)
lsc = lim_cntrl(T2[5:141])
plot_chart2(T2,dates1[l:n], LSC=lsc,
            Type="Adaptative T2 sur données d'entrainement avec lag constant")
plot_chart2(T2,dates1[l:n], LSC=lsc)
# T2_adapt sur les données totales
T2 = T2_adapt(X_tot2,l)
plot_chart2(T2, dates[l:n3], LSC=lsc,
            Type="Adaptative T2 sur toutes les données avec lag constant")
plot_chart2(T2[2:length(T2)], dates[(l+1):n3], LSC=lsc)



###################### adapt ewma ##########################
# cartes de contrôle
l = 1
# limite de contrôle avec les paramètres de train
lsc = LSC_T2_Hotelling_k1(X_ca2, alpha=c(0.05,0.01))
# T2_adapt_EWMA sur train
T2 = T2_adapt_EWMA(X_ca2, l,0.1)
plot_chart2(T2,dates1[l:n], LSC=lsc[1],
            Type="Adaptative T2 sur données d'entrainement avec lag constant")
plot_chart2(T2,dates1[l:n], LSC=lsc[1])
# T2_adapt_EWMA sur test
T2 = T2_adapt_EWMA(X_val2,l,0.1)
plot_chart2(T2, dates2[l:n2], LSC=lsc[2],
            Type="Adaptative T2 sur données de test avec lag constant")
plot_chart2(T2, dates2[l:n2], LSC=lsc[2])



# T2_adapt sur train
T2 = T2_adapt_EWMA(X_ca2, l,0.1)
lsc = lim_cntrl(T2, 0.01)
plot_chart2(T2,dates1[l:n], LSC=lsc,
            Type="Adaptative T2 sur données d'entrainement avec lag constant")
plot_chart2(T2,dates1[l:n], LSC=lsc)
# T2_adapt sur les données totales
T2 = T2_adapt_EWMA(X_tot2,l,0.1)
plot_chart2(T2[13:695], dates[(l+12):n3], LSC=lsc,
            Type="Adaptative T2 sur toutes les données avec lag constant")
plot_chart2(T2[13:695], dates[(l+12):n3], LSC=lsc)



################################################
# détection des lags de chaque variable
nlag = 30 # lag maximal
cor_res = matrix(nrow = nlag, ncol = p) # 
res = list() # pic au dessus de la limite
colnames(cor_res) = colnames(X_ca2) # nom des variables
clim = qnorm((1 + 0.95)/2)/sqrt(nrow(X_ca2)) # calcul de la limite
for (col in colnames(X_ca2)) {
  cor_res[,col] = pacf(X_ca2[,col], lag.max = nlag,plot = FALSE)$acf
  pacf(X_ca2[,col], lag.max = nlag, main=col)
  # point en dehors des limites
  des = (cor_res[,col] < clim) & (cor_res[,col] > -clim)
  res[[col]] = which(des==FALSE)
}

res


####################### adapt #############################
# lag personalisé
l = c(2,3,3,1,1,1,1,1,2,2,1,2,1,5,3,1)
N = n-max(l)+1
# limite de contrôle avec les paramètres de train
lsc = LSC_T2_Hotelling_k1(matrix(nrow = N, ncol = sum(l)), alpha=c(0.05,0.05))
# T2_adapt sur train
T2 = T2_adapt_2(X_ca2,l)
plot_chart2(T2,dates1[(max(l)):n], LSC=lsc[1],
            Type="Adaptative T2 sur données d'entrainement")
plot_chart2(T2,dates1[(max(l)):n], LSC=lsc[1])
# T2_adapt sur test
T2 = T2_adapt_2(X_val2,l)
plot_chart2(T2, dates2[(max(l)):n2], LSC=lsc[2],
            Type="Adaptative T2 sur données de test")
plot_chart2(T2, dates2[(max(l)):n2], LSC=lsc[2])


# T2_adapt sur train
T2 = T2_adapt_2(X_ca2,l)
lsc = lim_cntrl(T2[5:141])
plot_chart2(T2,dates1[max(l):n], LSC=lsc,
            Type="Adaptative T2 sur données d'entrainement")
plot_chart2(T2,dates1[max(l):n], LSC=lsc)
# T2_adapt sur les données totales
T2 = T2_adapt_2(X_tot2,l)
plot_chart2(T2, dates[max(l):n3], LSC=lsc,
            Type="Adaptative T2 sur toutes les données")
plot_chart2(T2, dates[max(l):n3], LSC=lsc)


###################### adapt ewma ############################
# cartes de contrôle
l = c(2,3,3,1,1,1,1,1,2,2,1,2,1,5,3,1)
N = n-max(l)+1
# limite de contrôle avec les paramètres de train
lsc = LSC_T2_Hotelling_k1(matrix(nrow = N, ncol = sum(l)), alpha=c(0.05,0.05))
# T2_adapt_EWMA sur train
T2 = T2_adapt_EWMA2(X_ca2, l,0.1)
plot_chart2(T2,dates1[max(l):n], LSC=lsc[1],
            Type="Adaptative T2 sur données d'entrainement")
plot_chart2(T2,dates1[max(l):n], LSC=lsc[1])

# T2_adapt_EWMA sur test
T2 = T2_adapt_EWMA2(X_val2,l,0.1)
plot_chart2(T2, dates2[max(l):n2], LSC=lsc[2],
            Type="Adaptative T2 sur données de test")
plot_chart2(T2, dates2[max(l):n2], LSC=lsc[2])



# T2_adapt sur train
T2 = T2_adapt_EWMA2(X_ca2, l,0.3)
lsc = lim_cntrl(T2, 0.05)
plot_chart2(T2,dates1[max(l):n], LSC=lsc,
            Type="Adaptative T2 sur données d'entrainement")
plot_chart2(T2,dates1[max(l):n], LSC=lsc)

# T2_adapt sur les données totales
T2 = T2_adapt_EWMA2(X_tot2,l,0.3)
plot_chart2(T2, dates[max(l):n3], LSC=lsc,
            Type="Adaptative T2 sur toutes les données")
plot_chart2(T2, dates[max(l):n3], LSC=lsc)





####################################################################
################### MYT des points hors contrôle ###################
####################################################################
source("C:/Users/hp/OneDrive/Bureau/Work/ONIRIS/BAM/Codes/myt.R")

l = c(2,3,3,1,1,1,1,1,2,2,1,2,1,5,3,1)
l = matrix(l, nrow = 1)
colnames(l) = colnames(X_ca2)
n = nrow(X_ca2) # nombre d'observations de train
N = n-max(l)+1
# limite de contrôle avec les paramètres de train
lsc = LSC_T2_Hotelling_k1(matrix(nrow = N, ncol = sum(l)), alpha=c(0.05,0.05))

# T2_adapt_EWMA sur test
T2 = T2_adapt_EWMA2(X_val2,l,0.1)
df <- data.frame(t = dates2[max(l):n2], values = T2)
# Marquage des points hors contrôle
df$out_of_control <- rep(FALSE, nrow(df))
df$out_of_control <- df$out_of_control | (df$values > lsc[2])
# les dates où les points sont hors contrôle
df[df['out_of_control']==TRUE,]
# application de MYT et affichage
for (i in 1:nrow(df)) {
  if (df$out_of_control[i]){
    resultat = MYT_adapt(X_val2, i, n,l,0.1, alpha = 0.05)
    print(i)
    choix = which(resultat$exclu=="oui")
    print(resultat[choix,])
  }
}

########################## DiPCA ##########################

df <- data.frame(t = dates2[(s_dipca+1):n2], values = T2_dipca)
# Marquage des points hors contrôle
df$out_of_control <- rep(FALSE, nrow(df))
df$out_of_control <- df$out_of_control | (df$values > lsc_dipca)
# les dates où les points sont hors contrôle
df[df['out_of_control']==TRUE,]
# application de MYT et affichage
for (i in 1:nrow(df)) {
  if (df$out_of_control[i]){
    resultat = MYT(X_val2, i, n, alpha = 0.05)
    print(i)
    choix = which(resultat$exclu=="oui")
    print(resultat[choix,])
  }
}
