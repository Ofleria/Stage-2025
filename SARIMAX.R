
library(xts)
library(dynlm)
library(zoo)
library(forecast)
library(ggplot2)
library(dplyr)
library(changepoint)
library(strucchange)



source("C:/Users/hp/OneDrive/Bureau/Cours/S3/XMS3MU120_Séries Temporelles/TP4.R")

setwd("C:/Users/hp/OneDrive/Bureau/Work/ONIRIS/BAM/Codes/données")
data_hour_O2=readRDS('data_hour_O2.RDS')
data = as.matrix(data_hour_O2[(data_hour_O2[,"jour"]<=34 & data_hour_O2[,"jour"]>5),])

for (col in colnames(data)) {
  data[,col] = na.approx(data[,col])
}
# séparation train / test
train_hour = data[data[,"jour"]<=12,-32]
test_hour = data[data[,"jour"]>12,-32]
rownames(train_hour) = NULL
rownames(test_hour) = NULL

# sélection de variables
data_train = train_hour[, c("O2", "Temperature", "ensoleil", "Turb", "Couleur", "temp_ext")]
data_test = test_hour[, c("O2", "Temperature", "ensoleil", "Turb", "Couleur", "temp_ext")]

# Génération d'horodatages horaires
dates <- seq.POSIXt(from = as.POSIXct("2022-05-31 00:00"), by = "hour", length.out = 696)
dates1 <- seq.POSIXt(from = as.POSIXct("2022-05-31 00:00"), by = "hour", length.out = 168)
dates2 <- seq.POSIXt(from = as.POSIXct("2022-06-07 00:00"), by = "hour", length.out = 528)

# Série temporelle avec horodatage
data_tot = xts(data, order.by = dates)
data1 <- xts(data_train, order.by = dates1)
data3 <- xts(data_test, order.by = dates2)

# Affichage
plot(data_tot[,"O2"],col = "purple",  main = "Oxygène", ylab = "O2", lty = 1, lwd = 2)
plot(data_tot[,"ensoleil"], col = "purple", main = "Ensoleillement", ylab = "ensoleil", lty = 1, lwd =1 )
plot(data_tot[,"Turb"],col = "purple",  main = "Turbidité", ylab = "Valeur", lty = 1, lwd = 2)

# données sous forme de série temporelle
data = ts(data, frequency = 24)
data2 = ts(data_train, frequency = 24)
data4 = ts(data_test, frequency = 24, start=c(end(data2)[1]+1, 1))

# variables exogènes
xreg = data2[, c("Temperature", "ensoleil", "temp_ext")]
xreg_future = data4[, c("Temperature", "ensoleil", "temp_ext")]

# H0 de Ljung-Box : les résidus ne sont pas autocorrélés 


#####################################################
# couleur 
mod = Arima(data2[,"Couleur"], order=c(1,0,0), seasonal=list(order=c(1,1,0),period=24), 
            xreg = xreg, include.drift = TRUE)
summary(mod)
# analyse des résidus
checkupRes(mod$residuals)
Box.test(mod$residuals, lag = 10, type = "Ljung-Box")

# prédiction
forecast_arimax = forecast(mod, xreg = xreg_future, level = c(80, 95))
pred_df <- data.frame(
  time = dates2,
  forecast = as.numeric(forecast_arimax$mean),
  lower80 = as.numeric(forecast_arimax$lower[,1]),
  upper80 = as.numeric(forecast_arimax$upper[,1]),
  lower95 = as.numeric(forecast_arimax$lower[,2]),
  upper95 = as.numeric(forecast_arimax$upper[,2])
)

# Observé
observed_df <- data.frame(
  time = dates,
  O2 = as.numeric(data[,"Couleur"])
)

# Combiner tout avec ggplot
ggplot() +
  geom_line(data = observed_df, aes(x = time, y = O2), color = "black", size = 0.6) +
  geom_line(data = pred_df, aes(x = time, y = forecast), color = "red", size = 0.6) +
  geom_ribbon(data = pred_df, aes(x = time, ymin = lower95, ymax = upper95), fill = "orange", alpha = 0.3) +
  geom_ribbon(data = pred_df, aes(x = time, ymin = lower80, ymax = upper80), fill = "brown", alpha = 0.3) +
  labs(#title = "Prévision de Couleur",
       x = "Date",
       y = "Couleur") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"))

des = data4[,'Couleur'] > forecast_arimax$lower[,'95%'] & data4[,'Couleur'] < forecast_arimax$upper[,'95%']
des_count = rollapply(des, width=5, FUN=function(w) sum(w==FALSE), align="left", partial=TRUE)
# la première fenêtre avec que des false
detect = which(des_count==5)[1:2]
# la date correspondante
dates2[detect]


###################################################
# O2
# S ARI X
mod = Arima(data2[,"O2"], order=c(2,0,0), seasonal=list(order=c(2,1,0),period=24), xreg = xreg)
summary(mod)
# analyse des résidus
checkupRes(mod$residuals)
Box.test(mod$residuals, lag = 10, type = "Ljung-Box")

# prédiction
forecast_arimax = forecast(mod, xreg = xreg_future, level = c(80, 95))
pred_df <- data.frame(
  time = dates2,
  forecast = as.numeric(forecast_arimax$mean),
  lower80 = as.numeric(forecast_arimax$lower[,1]),
  upper80 = as.numeric(forecast_arimax$upper[,1]),
  lower95 = as.numeric(forecast_arimax$lower[,2]),
  upper95 = as.numeric(forecast_arimax$upper[,2])
)

# Observé
observed_df <- data.frame(
  time = dates,
  O2 = as.numeric(data[,"O2"])
)

# Combiner tout avec ggplot
ggplot() +
  geom_line(data = observed_df, aes(x = time, y = O2), color = "black", size = 0.5) +
  geom_line(data = pred_df, aes(x = time, y = forecast), color = "red", size = 0.5) +
  geom_ribbon(data = pred_df, aes(x = time, ymin = lower95, ymax = upper95), fill = "orange", alpha = 0.3) +
  geom_ribbon(data = pred_df, aes(x = time, ymin = lower80, ymax = upper80), fill = "brown", alpha = 0.3) +
  labs(#title = "Prévision de O2",
       x = "Date",
       y = "O2") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"))


# détection des points inhabituels 
# données présentes ou non dans l'intervalle de confiance
des = data4[,'O2'] > forecast_arimax$lower[,'95%'] & data4[,'O2'] < forecast_arimax$upper[,'95%']
# le nombre de false dans une fenêtre glissante de taille 5
des_count = rollapply(des, width=5, FUN=function(w) sum(w==FALSE), align="left", partial=TRUE)
# la première fenêtre avec que des false
detect = which(des_count==5)[1:2]
# la date correspondante
dates2[detect]


# distances entre prédiction et vraie valeur
distances = abs(data4[,"O2"] - forecast_arimax$mean)
plot(distances, type='l')

# points de ruptures
cp = cpt.meanvar(distances, method = 'PELT')
plot(cp, main = "Ruptures détectées", xlab = "Jour", ylab = "Distance")
legend("topright", legend = "Ruptures", col = "red", lty = 1, lwd = 2)

cpts(cp)
dates2[cpts(cp)]

bp = breakpoints(distances ~ 1)
summary(bp)



########################################################
# turbidité
mod = Arima(data2[,"Turb"], order=c(1,0,0), seasonal=list(order=c(2,1,0),period=24), 
            xreg = xreg, include.drift = TRUE)
summary(mod)
# analyse des résidus
checkupRes(mod$residuals)
Box.test(mod$residuals, lag = 10, type = "Ljung-Box")

# prédiction
forecast_arimax = forecast(mod, xreg = xreg_future, level = c(80, 95))
pred_df <- data.frame(
  time = dates2,
  forecast = as.numeric(forecast_arimax$mean),
  lower80 = as.numeric(forecast_arimax$lower[,1]),
  upper80 = as.numeric(forecast_arimax$upper[,1]),
  lower95 = as.numeric(forecast_arimax$lower[,2]),
  upper95 = as.numeric(forecast_arimax$upper[,2])
)

# Observé
observed_df <- data.frame(
  time = dates,
  O2 = as.numeric(data[,"Turb"])
)

# Combiner tout avec ggplot
ggplot() +
  geom_line(data = observed_df, aes(x = time, y = O2), color = "black", size = 0.5) +
  geom_line(data = pred_df, aes(x = time, y = forecast), color = "red", size = 0.5) +
  geom_ribbon(data = pred_df, aes(x = time, ymin = lower95, ymax = upper95), fill = "orange", alpha = 0.3) +
  geom_ribbon(data = pred_df, aes(x = time, ymin = lower80, ymax = upper80), fill = "brown", alpha = 0.3) +
  labs(#title = "Prévision de Turbidité",
       x = "Date",
       y = "Turbidité") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"))



# détection des points inhabituels 
# données présentes ou non dans l'intervalle de confiance
des = data4[,'Turb'] > forecast_arimax$lower[,'95%'] & data4[,'Turb'] < forecast_arimax$upper[,'95%']
# le nombre de false dans une fenêtre glissante de taille 5
des_count = rollapply(des, width=5, FUN=function(w) sum(w==FALSE), align="left", partial=TRUE)
# la première fenêtre avec que des false
detect = which(des_count==5)[1:10]
# la date correspondante
dates2[detect]

# distances entre prédiction et vraie valeur
distances = abs(data4[,"Turb"] - forecast_arimax$mean)

cp = cpt.meanvar(distances, method = 'PELT')
plot(cp, main = "Ruptures détectées", xlab = "Jour", ylab = "Distance")
legend("topleft", legend = "Ruptures", col = "red", lty = 1, lwd = 2)
cpts(cp)

dates2[cpts(cp)]


bp = breakpoints(distances ~ 1)

plot(bp)
summary(bp)


########################################################
# O2 chaque jour
date = seq(from=as.Date("2022-05-26 00:00"), 
           to=as.Date("2022-07-05 00:00"), by="day")

# O2
for (i in 12:25) {
  data = as.matrix(data_hour_O2[(data_hour_O2[,"jour"]<=i+1 & data_hour_O2[,"jour"]>5),])
  for (col in colnames(data)) {
    data[,col] = na.approx(data[,col])
  }
  # séparation
  data_train = data[data[,"jour"]<=i,-32]
  data_test = data[data[,"jour"]>i,-32]
  
  rownames(data_train) = NULL
  rownames(data_test) = NULL
  
  dates <- seq.POSIXt(from = as.POSIXct(date[6]), by = "hour", length.out = nrow(data))
  dates2 <- seq.POSIXt(from = as.POSIXct(date[i+1]), by = "hour", length.out = 24)
  
  data = ts(data, frequency = 24)
  data2 = ts(data_train, frequency = 24)
  data4 = ts(data_test, frequency = 24, start=c(end(data2)[1]+1, 1))
  # données exogènes
  xreg = data2[, c("Temperature", "ensoleil", "temp_ext")]
  xreg_future = data4[, c("Temperature", "ensoleil", "temp_ext")]
  
  # S ARI X
  mod = Arima(data2[,"O2"], order=c(2,0,0), seasonal=list(order=c(2,1,0),period=24), xreg = xreg)
  # prédiction
  forecast_arimax = forecast(mod, xreg = xreg_future, level = c(80, 95))
  # données présentes ou non dans l'intervalle de confiance
  des = data4[,'O2'] > forecast_arimax$lower[,'95%'] & data4[,'O2'] < forecast_arimax$upper[,'95%']
  # nombre d’occurrence
  des_count = table(des)
  print(des_count)
  
  pred_df <- data.frame(
    time = dates2,
    forecast = as.numeric(forecast_arimax$mean),
    lower80 = as.numeric(forecast_arimax$lower[,1]),
    upper80 = as.numeric(forecast_arimax$upper[,1]),
    lower95 = as.numeric(forecast_arimax$lower[,2]),
    upper95 = as.numeric(forecast_arimax$upper[,2])
  )
  
  # Observé
  observed_df <- data.frame(
    time = dates,
    O2 = as.numeric(data[,"O2"])
  )
  
  # Combiner tout avec ggplot
  p = ggplot() +
    geom_line(data = observed_df, aes(x = time, y = O2), color = "black", size = 0.8) +
    geom_line(data = pred_df, aes(x = time, y = forecast), color = "red", size = 0.8) +
    geom_ribbon(data = pred_df, aes(x = time, ymin = lower95, ymax = upper95), fill = "orange", alpha = 0.3) +
    geom_ribbon(data = pred_df, aes(x = time, ymin = lower80, ymax = upper80), fill = "brown", alpha = 0.3) +
    labs(x = "Date",
         y = "O2") +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold"))
  
  print(p)
  
  if (des_count["TRUE"]<14) break
}



#################
# turbidité chaqur jour
for (i in 12:25) {
  data = as.matrix(data_hour_O2[(data_hour_O2[,"jour"]<=i+1 & data_hour_O2[,"jour"]>5),])
  
  for (col in colnames(data)) {
    data[,col] = na.approx(data[,col])
  }
  # séparation
  data_train = data[data[,"jour"]<=i,-32]
  data_test = data[data[,"jour"]>i,-32]
  
  rownames(data_train) = NULL
  rownames(data_test) = NULL
  
  dates <- seq.POSIXt(from = as.POSIXct(date[6]), by = "hour", length.out = nrow(data))
  dates2 <- seq.POSIXt(from = as.POSIXct(date[i+1]), by = "hour", length.out = 24)
  
  data = ts(data, frequency = 24)
  data2 = ts(data_train, frequency = 24)
  data4 = ts(data_test, frequency = 24, start=c(end(data2)[1]+1, 1))
  # données exogènes
  xreg = data2[, c("Temperature", "ensoleil", "temp_ext")]
  xreg_future = data4[, c("Temperature", "ensoleil", "temp_ext")]
  
  # S ARI X
  mod = Arima(data2[,"Turb"], order=c(1,0,0), seasonal=list(order=c(2,1,0),period=24), xreg = xreg)
  # prédiction
  forecast_arimax = forecast(mod, xreg = xreg_future, level = c(80, 95))
  # données présentes ou non dans l'intervalle de confiance
  des = data4[,'Turb'] > forecast_arimax$lower[,'95%'] & data4[,'Turb'] < forecast_arimax$upper[,'95%']
  # nombre d’occurrence
  des_count = table(des)
  print(des_count)
  
  pred_df <- data.frame(
    time = dates2,
    forecast = as.numeric(forecast_arimax$mean),
    lower80 = as.numeric(forecast_arimax$lower[,1]),
    upper80 = as.numeric(forecast_arimax$upper[,1]),
    lower95 = as.numeric(forecast_arimax$lower[,2]),
    upper95 = as.numeric(forecast_arimax$upper[,2])
  )
  # Observé
  observed_df <- data.frame(
    time = dates,
    O2 = as.numeric(data[,"Turb"])
  )
  # Combiner tout avec ggplot
  p = ggplot() +
    geom_line(data = observed_df, aes(x = time, y = O2), color = "black", size = 0.5) +
    geom_line(data = pred_df, aes(x = time, y = forecast), color = "red", size = 0.5) +
    geom_ribbon(data = pred_df, aes(x = time, ymin = lower95, ymax = upper95), fill = "orange", alpha = 0.3) +
    geom_ribbon(data = pred_df, aes(x = time, ymin = lower80, ymax = upper80), fill = "brown", alpha = 0.3) +
    labs(#title = "Prévision de turbidité",
         x = "Date",
         y = "Turbidité") +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold"))
  
  print(p)
  
  if (des_count["TRUE"]<14) break
  
}







######################################################
# couleur par ACP
data = as.matrix(data_hour_O2[(data_hour_O2[,"jour"]<=34 & data_hour_O2[,"jour"]>5),])

for (col in colnames(data)) {
  data[,col] = na.approx(data[,col])
}
# séparation
train_hour = data[data[,"jour"]<=12,-32]
test_hour = data[data[,"jour"]>12,-32]
rownames(train_hour) = NULL
rownames(test_hour) = NULL

# sélection de variables
data_train = train_hour[, c("O2", "Temperature", "ensoleil", "Turb", "Couleur", "temp_ext")]
data_test = test_hour[, c("O2", "Temperature", "ensoleil", "Turb", "Couleur", "temp_ext")]
# les absorbances
df = train_hour[,c(8:11, 14:17, 20:23, 25:28)]

# ACP des absorbance pour créer une nouvelle variable couleur
acp_res = prcomp(df, center = TRUE, scale=TRUE)
summary(acp_res)

# sur les données de test
ts = scale(test_hour[,c(8:11, 14:17, 20:23, 25:28)], center=acp_res$center, scale=acp_res$scale)
ts_red = ts %*% acp_res$rotation

# ajout aus données sélectionnées
data_train = cbind(data_train, pc1 = acp_res$x[,1])
data_test = cbind(data_test, pc1 = ts_red[, 1])
data = rbind(data_train, data_test)

# trasformation en time series
data2 = ts(data_train, frequency = 24)
data4 = ts(data_test, frequency = 24, start=c(end(data2)[1]+1, 1))
data5 = ts(data, frequency = 24)

# variables exogènes
xreg = data2[, c("Temperature", "ensoleil", "temp_ext")]
xreg_future = data4[, c("Temperature", "ensoleil", "temp_ext")]

# SARIMAX
mod = Arima(data2[,"pc1"], order=c(2,0,0), seasonal=list(order=c(1,1,0),period=24), 
            xreg = xreg, include.drift = TRUE)
forecast_arimax = forecast(mod, xreg = xreg_future, level = c(80, 95))

summary(mod)
# analyse des résidus
checkupRes(mod$residuals)
Box.test(mod$residuals, lag = 10, type = "Ljung-Box")

# prédiction
forecast_arimax = forecast(mod, xreg = xreg_future, level = c(80, 95))
pred_df <- data.frame(
  time = dates2,
  forecast = as.numeric(forecast_arimax$mean),
  lower80 = as.numeric(forecast_arimax$lower[,1]),
  upper80 = as.numeric(forecast_arimax$upper[,1]),
  lower95 = as.numeric(forecast_arimax$lower[,2]),
  upper95 = as.numeric(forecast_arimax$upper[,2])
)

# Observé
observed_df <- data.frame(
  time = dates,
  O2 = as.numeric(data[,"pc1"])
)

# Combiner tout avec ggplot
ggplot() +
  geom_line(data = observed_df, aes(x = time, y = O2), color = "black", size = 0.5) +
  geom_line(data = pred_df, aes(x = time, y = forecast), color = "red", size = 0.5) +
  geom_ribbon(data = pred_df, aes(x = time, ymin = lower95, ymax = upper95), fill = "orange", alpha = 0.3) +
  geom_ribbon(data = pred_df, aes(x = time, ymin = lower80, ymax = upper80), fill = "brown", alpha = 0.3) +
  labs(#title = "Prévision de couleur",
       x = "Date",
       y = "Couleur") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"))



# détection des points inhabituels 
# données présentes ou non dans l'intervalle de confiance
des = data4[,'pc1'] > forecast_arimax$lower[,'95%'] & data4[,'pc1'] < forecast_arimax$upper[,'95%']
# le nombre de false dans une fenêtre glissante de taille 5
des_count = rollapply(des, width=5, FUN=function(w) sum(w==FALSE), align="left", partial=TRUE)
# la première fenêtre avec que des false
detect = which(des_count==5)[1:3]
# la date correspondante
dates2[detect]

# distances entre prédiction et vraie valeur
distances = abs(data4[,"pc1"] - forecast_arimax$mean)
plot(distances, type='l')


cp = cpt.meanvar(distances, method = 'PELT')
plot(cp, main = "Ruptures détectées", xlab = "Jour", ylab = "Distance")
legend("topleft", legend = "Ruptures", col = "red", lty = 1, lwd = 2)
cpts(cp)

dates2[cpts(cp)]

bp = breakpoints(distances ~ 1)

plot(bp)
summary(bp)


#############
# chaque jour
date = seq(from=as.Date("2022-05-26 00:00"), 
           to=as.Date("2022-07-05 00:00"), by="day")
for (i in 12:25) {
  data = as.matrix(data_hour_O2[(data_hour_O2[,"jour"]<=i+1 & data_hour_O2[,"jour"]>5),])
  
  for (col in colnames(data)) {
    data[,col] = na.approx(data[,col])
  }
  # séparation
  data_train = data[data[,"jour"]<=i,-32]
  data_test = data[data[,"jour"]>i,-32]
  
  rownames(data_train) = NULL
  rownames(data_test) = NULL
  
  dates <- seq.POSIXt(from = as.POSIXct(date[6]), by = "hour", length.out = nrow(data))
  dates2 <- seq.POSIXt(from = as.POSIXct(date[i+1]), by = "hour", length.out = 24)
  
  
  acp_res = prcomp(data_train[,c(8:11, 14:17, 20:23, 25:28)], center = TRUE, scale = TRUE)
  # test
  ts = scale(data_test[,c(8:11, 14:17, 20:23, 25:28)], center=acp_res$center, scale=acp_res$scale)
  ts_red = ts %*% acp_res$rotation
  
  data_train = cbind(data_train, pc1 = acp_res$x[,1])
  data_test = cbind(data_test, pc1 = ts_red[, 1])
  data = rbind(data_train, data_test)
  
  data = ts(data, frequency = 24)
  data2 = ts(data_train, frequency = 24)
  data4 = ts(data_test, frequency = 24, start=c(end(data2)[1]+1, 1))
  # données exogènes
  xreg = data2[, c("Temperature", "ensoleil", "temp_ext")]
  xreg_future = data4[, c("Temperature", "ensoleil", "temp_ext")]
  
  # S ARI X
  mod = Arima(data2[,"pc1"], order=c(2,0,0), seasonal=list(order=c(1,1,0),period=24), 
              xreg = xreg, include.drift = TRUE)  # prédiction
  forecast_arimax = forecast(mod, xreg = xreg_future, level = c(80, 95))
  # données présentes ou non dans l'intervalle de confiance
  des = data4[,'pc1'] > forecast_arimax$lower[,'95%'] & data4[,'pc1'] < forecast_arimax$upper[,'95%']
  # nombre d’occurrence
  des_count = table(des)
  print(des_count)
  
  pred_df <- data.frame(
    time = dates2,
    forecast = as.numeric(forecast_arimax$mean),
    lower80 = as.numeric(forecast_arimax$lower[,1]),
    upper80 = as.numeric(forecast_arimax$upper[,1]),
    lower95 = as.numeric(forecast_arimax$lower[,2]),
    upper95 = as.numeric(forecast_arimax$upper[,2])
  )
  
  # Observé
  observed_df <- data.frame(
    time = dates,
    O2 = as.numeric(data[,"pc1"])
  )
  
  # Combiner tout avec ggplot
  p = ggplot() +
    geom_line(data = observed_df, aes(x = time, y = O2), color = "black", size = 0.8) +
    geom_line(data = pred_df, aes(x = time, y = forecast), color = "red", size = 0.8) +
    geom_ribbon(data = pred_df, aes(x = time, ymin = lower95, ymax = upper95), fill = "orange", alpha = 0.3) +
    geom_ribbon(data = pred_df, aes(x = time, ymin = lower80, ymax = upper80), fill = "brown", alpha = 0.3) +
    labs(x = "Date",
         y = "pc1") +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold"))
  
  print(p)
  
  if (des_count["TRUE"]<14) break
  
}



#################################################################
# couleur par interpolation
lambda <- c(465, 525, 595, 625)
# vecteur de longueurs d'onde cibles
lambda_interp <- seq(400, 700, by=1)

data = as.matrix(data_hour_O2[(data_hour_O2[,"jour"]<=34 & data_hour_O2[,"jour"]>5),])
for (col in colnames(data)) {
  data[,col] = na.approx(data[,col])
}
# séparation
train_hour = data[data[,"jour"]<=12,-32]
test_hour = data[data[,"jour"]>12,-32]

rownames(train_hour) = NULL
rownames(test_hour) = NULL

# sélection de variables
data_train = train_hour[, c("O2", "Temperature", "ensoleil", "Turb", "Couleur", "temp_ext")]
data_test = test_hour[, c("O2", "Temperature", "ensoleil", "Turb", "Couleur", "temp_ext")]

# recherche pour toutes les observations
coul = c()
for (i in 1:nrow(train_hour)) {
  abs_interp <- spline(x = lambda, y = train_hour[i,25:28], xout = lambda_interp)$y
  coul = rbind(coul, lambda_interp[which.min(abs_interp)])
}
which(coul==700)
#coul[589]=mean(c(coul[c(588,590)]))

data_train = cbind(data_train, data.frame(Coul = coul))


coul = c()
for (i in 1:nrow(test_hour)) {
  abs_interp <- spline(x = lambda, y = test_hour[i,25:28], xout = lambda_interp)$y
  coul = rbind(coul, lambda_interp[which.min(abs_interp)])
}
# remplacement du point à 700
which(coul==700)
coul[421]=mean(c(coul[c(420,422)]))
data_test = cbind(data_test, data.frame(Coul = coul))

data = rbind(data_train, data_test)


data = ts(data, frequency = 24)
data2 = ts(data_train, frequency = 24)
data4 = ts(data_test, frequency = 24, start=c(end(data2)[1]+1, 1))
# données exogènes
xreg = data2[, c("Temperature", "ensoleil", "temp_ext")]
xreg_future = data4[, c("Temperature", "ensoleil", "temp_ext")]

# S ARI X
mod = Arima(data2[,"Coul"], order=c(1,0,0), seasonal=list(order=c(3,1,0),period=24), 
            xreg = xreg, include.drift = TRUE)  # prédiction

summary(mod)
# analyse des résidus
checkupRes(mod$residuals)
Box.test(mod$residuals, lag = 10, type = "Ljung-Box")

forecast_arimax = forecast(mod, xreg = xreg_future, level = c(80, 95))

# données présentes ou non dans l'intervalle de confiance
des = data4[,'Coul'] > forecast_arimax$lower[,'95%'] & data4[,'Coul'] < forecast_arimax$upper[,'95%']
# le nombre de false dans une fenêtre glissante de taille 5
des_count = rollapply(des, width=5, FUN=function(w) sum(w==FALSE), align="left", partial=TRUE)
# la première fenêtre avec que des false
detect = which(des_count==5)[1:5]
# la date correspondante
dates2[detect]

# distances entre prédiction et vraie valeur
distances = abs(data4[,"Coul"] - forecast_arimax$mean)
plot(distances, type='l')


cp = cpt.meanvar(distances, method = 'PELT')
plot(cp, main = "Ruptures détectées", xlab = "Jour", ylab = "Distance")
legend("topleft", legend = "Ruptures", col = "red", lty = 1, lwd = 2)
cpts(cp)

dates2[cpts(cp)]


pred_df <- data.frame(
  time = dates2,
  forecast = as.numeric(forecast_arimax$mean),
  lower80 = as.numeric(forecast_arimax$lower[,1]),
  upper80 = as.numeric(forecast_arimax$upper[,1]),
  lower95 = as.numeric(forecast_arimax$lower[,2]),
  upper95 = as.numeric(forecast_arimax$upper[,2])
)

# Observé
observed_df <- data.frame(
  time = dates,
  O2 = as.numeric(data[,"Coul"])
)

# Combiner tout avec ggplot
ggplot() +
  geom_line(data = observed_df, aes(x = time, y = O2), color = "black", size = 0.5) +
  geom_line(data = pred_df, aes(x = time, y = forecast), color = "red", size = 0.5) +
  geom_ribbon(data = pred_df, aes(x = time, ymin = lower95, ymax = upper95), fill = "orange", alpha = 0.3) +
  geom_ribbon(data = pred_df, aes(x = time, ymin = lower80, ymax = upper80), fill = "brown", alpha = 0.3) +
  labs(x = "Date",
       y = "Couleur") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"))
 

