

checkupRes = function(Res) {
  
  # Normaliser les résidus (centrer / réduire)
  Res_norm = scale(Res)
  
  # Configuration de la mise en page des graphiques
  layout(matrix(c(1,1,1,2:7), nrow=3, ncol=3, byrow=TRUE))
  
  # Graphique 1 : Evolution de la série
  plot(Res, type='l', main="Evolution de la série", xlab="Temps", ylab="Série", col='blue')
  
  # Graphique 2 : ACF 
  acf(Res, main="ACF")
  
  # Graphique 3 : PACF 
  pacf(Res, main="PACF")
  
  # Graphique 4 : Nuage de points Res[i] en fonction de Res[i-1]
  plot(Res[-length(Res)], Res[-1], main="X_i en fonction de X_i-1", xlab="X_i-1", ylab="X_i", pch=16, col='blue')
  
  # Graphique 5 : Histogramme des résidus
  hist(Res, breaks=20, main="Histogramme", xlab="Résidus", col='lightblue', border='white')
  
  # Graphique 6 : QQ plot par rapport aux quantiles gaussiens
  qqnorm(Res, main="QQ Plot")
  qqline(Res, col='red')
  
  # Graphique 7 : Nuage de points de la série renormalisée, avec lignes horizontales à ±1.96
  plot(Res_norm, type='p', main="Série renormalisée", xlab="Temps", ylab="Résidus normés", pch=16, col='blue')
  abline(h=c(-1.96, 1.96), col='red', lty=2)
  
  par(mfrow=c(1,1))
}

# pour enlever le layout
# par(mfrow=c(1,1))

# Exemple d'utilisation

#library(TSA)

# chargement des données electricity
#data("electricity")
# observation des données
#plot(electricity)
# considération de la série transformée Yt = ln(Xt)
#log_elec = log(electricity)
#plot(log_elec)

# récupération des fluctuations
#decomp_elec = decompose(log_elec) 
#fluc_elec = na.omit(decomp_elec$random)


#checkupRes(fluc_elec)

#par(mfrow=c(1,1))

# fonction de forecast
#checkresiduals(fluc_elec)


# source("TP4.R") pour utiliser la fonction dans d'autres fichiers
