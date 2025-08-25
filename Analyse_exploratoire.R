# Charger les packages nécessaires
library(ggplot2)
library(dplyr)
library(zoo)
library(tidyr)
library(lubridate)
library(GGally)
library(colorscience)
library(patchwork)  # pour empiler facilement plusieurs ggplot

setwd("C:/Users/hp/OneDrive/Bureau/Work/ONIRIS/BAM/Codes/données")

data = readRDS('data_hour_O2.RDS')
summary(data)

data[,"dates"] = seq.POSIXt(from = as.POSIXct("2022-05-26 00:00:00"), 
                            by = "hour", length.out = 984)

data = data[(data[,"jour"]<=12& data[,"jour"]>5),]

df = data[, c("dates", "O2", "Temperature", "temp_ext", 
              "ensoleil", "capacitance", "conductivity")]
colnames(df) = c("dates", "Oxygène", "Temp_interne", "Temp_externe", 
                 "Ensoleillement", "Capacitance", "Conductivité")

df[,2:7] = scale(df[,2:7])

# Passage en format long et extraction de l’heure
df_long <- df %>%
  pivot_longer(cols = colnames(df)[-1], names_to = "variable", values_to = "value") %>%
  mutate(hour = hour(dates))

df_long$hour <- factor(df_long$hour, levels = 0:23)

# Moyenne par heure et par variable
df_agg <- df_long %>%
  group_by(variable, hour) %>%
  summarise(mean_value = mean(value), .groups = "drop")

# Heatmap : heure en X, variable en Y
ggplot(df_agg, aes(x = hour, y = variable, fill = mean_value)) +
  geom_tile() +
  scale_fill_gradient(low = "yellow", high = "red") +
  labs(#title = "Heatmap par Heure et Variable",
    x = "Heure", y = "", fill = "Valeur") +
  theme_minimal() +
  geom_hline(yintercept = seq(1.5, length(unique(df_agg$variable)) - 0.5, by = 1),
             color = "white", linewidth = 0.5)





#############################################################
data_hour_O2=readRDS('data_hour_O2.RDS')
data = as.matrix(data_hour_O2[(data_hour_O2[,"jour"]<=34 
                               & data_hour_O2[,"jour"]>5),])

for (col in colnames(data)) {
  data[,col] = na.approx(data[,col])
}

train_hour = data[data[,"jour"]<=11,-32]
test_hour = data[data[,"jour"]>11,-32]

rownames(train_hour) = NULL
rownames(test_hour) = NULL
rownames(data) = NULL

# séquence de dates
dates = seq.POSIXt(from = as.POSIXct("2022-05-31 00:00"), by = "hour", length.out = 696)

####################################################
# abs et couleur
df = train_hour[,c(6, 8:11, 14:17, 20:23, 25:28)]
# Corrélation avec la variable cible "Couleur"
cor_with_x <- cor(df, use = "complete.obs")[, "Couleur"]
# Mise en forme du tableau
cor_df <- data.frame(variable = names(cor_with_x), correlation = cor_with_x) %>%
  filter(variable != "Couleur") %>%
  arrange(desc(abs(correlation)))
# Barplot
ggplot(cor_df, aes(x = variable, y = correlation)) +
  geom_col(fill = "brown") +
  coord_flip() +
  geom_text(aes(label = round(correlation, 2)), hjust = ifelse(cor_df$correlation > 0, -0.1, 1.1)) +
  theme_minimal() +
  labs(title = "Corrélation avec la variable Couleur",
       x = "Variables",
       y = "Corrélation") +
  ylim(c(-1, 1))

# abs et ensoleil
df = train_hour[,c(8:11, 14:17, 20:23, 25:28, 31)]
# Corrélation avec la variable cible "Couleur"
cor_with_x <- cor(df, use = "complete.obs")[, "ensoleil"]
# Mise en forme du tableau
cor_df <- data.frame(variable = names(cor_with_x), correlation = cor_with_x) %>%
  filter(variable != "ensoleil") %>%
  arrange(desc(abs(correlation)))
# Barplot
ggplot(cor_df, aes(x = variable, y = correlation)) +
  geom_col(fill = "brown") +
  coord_flip() +
  geom_text(aes(label = round(correlation, 2)), hjust = ifelse(cor_df$correlation > 0, -0.1, 1.1)) +
  theme_minimal() +
  labs(title = "Corrélation avec la variable ensoleil",
       x = "Variables",
       y = "Corrélation") +
  ylim(c(-1, 1))


# abs et temp
df = train_hour[,c(3, 8:11, 14:17, 20:23, 25:28)]
# Corrélation avec la variable cible "Couleur"
cor_with_x <- cor(df, use = "complete.obs")[, "Temperature"]
# Mise en forme du tableau
cor_df <- data.frame(variable = names(cor_with_x), correlation = cor_with_x) %>%
  filter(variable != "Temperature") %>%
  arrange(desc(abs(correlation)))
# Barplot
ggplot(cor_df, aes(x = variable, y = correlation)) +
  geom_col(fill = "brown") +
  coord_flip() +
  geom_text(aes(label = round(correlation, 2)), hjust = ifelse(cor_df$correlation > 0, -0.1, 1.1)) +
  theme_minimal() +
  labs(title = "Corrélation avec la variable Temperature",
       x = "Variables",
       y = "Corrélation") +
  ylim(c(-1, 1))


# abs et temp ext
df = train_hour[,c(30, 8:11, 14:17, 20:23, 25:28)]
# Corrélation avec la variable cible "Couleur"
cor_with_x <- cor(df, use = "complete.obs")[, "temp_ext"]
# Mise en forme du tableau
cor_df <- data.frame(variable = names(cor_with_x), correlation = cor_with_x) %>%
  filter(variable != "temp_ext") %>%
  arrange(desc(abs(correlation)))
# Barplot
ggplot(cor_df, aes(x = variable, y = correlation)) +
  geom_col(fill = "brown") +
  coord_flip() +
  geom_text(aes(label = round(correlation, 2)), hjust = ifelse(cor_df$correlation > 0, -0.1, 1.1)) +
  theme_minimal() +
  labs(title = "Corrélation avec la variable Temperature externe",
       x = "Variables",
       y = "Corrélation") +
  ylim(c(-1, 1))

#############################################################
# train 
df = train_hour[,-c(8:11, 14:17, 20:23, 25:28,4)]
# Corrélation avec la variable "capacitance"
cor_with_x <- cor(df, use = "complete.obs")[, "capacitance"]
# Mise en forme du tableau
cor_df <- data.frame(variable = names(cor_with_x), correlation = cor_with_x) %>%
  filter(variable != "capacitance") %>%
  arrange(desc(abs(correlation)))
# Barplot
ggplot(cor_df, aes(x = variable, y = correlation)) +
  geom_col(fill = "purple") +
  coord_flip() +
  geom_text(aes(label = round(correlation, 2)), hjust = ifelse(cor_df$correlation > 0, -0.1, 1.1)) +
  theme_minimal() +
  labs(title = "La variable capacitance du 31/05 au 05/06",
       x = "Variables",
       y = "Corrélation") +
  ylim(c(-1, 1))

# Corrélation avec la variable cible "O2"
cor_with_x <- cor(df, use = "complete.obs")[, "O2"]
# Mise en forme du tableau
cor_df <- data.frame(variable = names(cor_with_x), correlation = cor_with_x) %>%
  filter(variable != "O2") %>%
  arrange(desc(abs(correlation)))
# Barplot
ggplot(cor_df, aes(x = variable, y = correlation)) +
  geom_col(fill = "purple") +
  coord_flip() +
  geom_text(aes(label = round(correlation, 2)), hjust = ifelse(cor_df$correlation > 0, -0.1, 1.1)) +
  theme_minimal() +
  labs(title = "La variable O2 du 31/05 au 05/06",
       x = "Variables",
       y = "Corrélation") +
  ylim(c(-1, 1))


# Corrélation avec la variable cible "Turb"
cor_with_x <- cor(df, use = "complete.obs")[, "Turb"]
# Mise en forme du tableau
cor_df <- data.frame(variable = names(cor_with_x), correlation = cor_with_x) %>%
  filter(variable != "Turb") %>%
  arrange(desc(abs(correlation)))
# Barplot
ggplot(cor_df, aes(x = variable, y = correlation)) +
  geom_col(fill = "purple") +
  coord_flip() +
  geom_text(aes(label = round(correlation, 2)), hjust = ifelse(cor_df$correlation > 0, -0.1, 1.1)) +
  theme_minimal() +
  labs(title = "La variable Turbidité du 31/05 au 05/06",
       x = "Variables",
       y = "Corrélation") +
  ylim(c(-1, 1))



#############################################################
# test 
df = test_hour[,-c(8:11, 14:17, 20:23, 25:28,4)]

# Corrélation avec la variable "capacitance"
cor_with_x <- cor(df, use = "complete.obs")[, "capacitance"]
# Mise en forme du tableau
cor_df <- data.frame(variable = names(cor_with_x), correlation = cor_with_x) %>%
  filter(variable != "capacitance") %>%
  arrange(desc(abs(correlation)))
# Barplot
ggplot(cor_df, aes(x = variable, y = correlation)) +
  geom_col(fill = "purple") +
  coord_flip() +
  geom_text(aes(label = round(correlation, 2)), hjust = ifelse(cor_df$correlation > 0, -0.1, 1.1)) +
  theme_minimal() +
  labs(title = "La variable capacitance du 06/06 au 28/06",
       x = "Variables",
       y = "Corrélation") +
  ylim(c(-1, 1))

# Corrélation avec la variable cible "O2"
cor_with_x <- cor(df, use = "complete.obs")[, "O2"]
# Mise en forme du tableau
cor_df <- data.frame(variable = names(cor_with_x), correlation = cor_with_x) %>%
  filter(variable != "O2") %>%
  arrange(desc(abs(correlation)))
# Barplot
ggplot(cor_df, aes(x = variable, y = correlation)) +
  geom_col(fill = "purple") +
  coord_flip() +
  geom_text(aes(label = round(correlation, 2)), hjust = ifelse(cor_df$correlation > 0, -0.1, 1.1)) +
  theme_minimal() +
  labs(title = "La variable O2 du 06/06 au 28/06",
       x = "Variables",
       y = "Corrélation") +
  ylim(c(-1, 1))

# Corrélation avec la variable cible "Turb"
cor_with_x <- cor(df, use = "complete.obs")[, "Turb"]
# Mise en forme du tableau
cor_df <- data.frame(variable = names(cor_with_x), correlation = cor_with_x) %>%
  filter(variable != "Turb") %>%
  arrange(desc(abs(correlation)))
# Barplot
ggplot(cor_df, aes(x = variable, y = correlation)) +
  geom_col(fill = "purple") +
  coord_flip() +
  geom_text(aes(label = round(correlation, 2)), hjust = ifelse(cor_df$correlation > 0, -0.1, 1.1)) +
  theme_minimal() +
  labs(title = "La variable Turbidité du 06/06 au 28/06",
       x = "Variables",
       y = "Corrélation") +
  ylim(c(-1, 1))


############################################################

df = train_hour[,-c(8:11, 14:17, 20:23, 25:28,4)]
# Corrélation avec la variable cible "Couleur"
cor_with_x <- cor(df, use = "complete.obs")[, "Couleur"]
# Mise en forme du tableau
cor_df <- data.frame(variable = names(cor_with_x), correlation = cor_with_x) %>%
  filter(variable != "Couleur") %>%
  arrange(desc(abs(correlation)))
# Barplot
ggplot(cor_df, aes(x = variable, y = correlation)) +
  geom_col(fill = "purple") +
  coord_flip() +
  geom_text(aes(label = round(correlation, 2)), hjust = ifelse(cor_df$correlation > 0, -0.1, 1.1)) +
  theme_minimal() +
  labs(title = "Corrélation avec la variable Couleur",
       x = "Variables",
       y = "Corrélation") +
  ylim(c(-1, 1))

df = test_hour[,-c(8:11, 14:17, 20:23, 25:28,4)]
# Corrélation avec la variable cible "Couleur"
cor_with_x <- cor(df, use = "complete.obs")[, "Couleur"]
# Mise en forme du tableau
cor_df <- data.frame(variable = names(cor_with_x), correlation = cor_with_x) %>%
  filter(variable != "Couleur") %>%
  arrange(desc(abs(correlation)))
# Barplot
ggplot(cor_df, aes(x = variable, y = correlation)) +
  geom_col(fill = "purple") +
  coord_flip() +
  geom_text(aes(label = round(correlation, 2)), hjust = ifelse(cor_df$correlation > 0, -0.1, 1.1)) +
  theme_minimal() +
  labs(title = "Corrélation avec la variable Couleur",
       x = "Variables",
       y = "Corrélation") +
  ylim(c(-1, 1))


df = train_hour[,c(8:11, 14:17, 20:23, 25:28)]
colnames(df) = c("Abs45blue" ,  "Abs45green" ,  "Abs45yellow" , "Abs45red",  "Abs90blue"  , "Abs90green",  
                 "Abs90yellow" , "Abs90red" , "Abs135blue",  "Abs135green" , "Abs135yellow", "Abs135red",
                 "Abs180blue" , "Abs180green" , "Abs180yellow" ,"Abs180red")

ggcorr(df, 
       label = TRUE, 
       label_round = 2, 
       label_size = 2.4,
       size=2.2,
       hjust = 0.7,
       low = "blue", mid = "white", high = "red")

# ACP des absorbance pour créer une nouvelle variable couleur
acp_res = prcomp(df, center = TRUE, scale=TRUE)
summary(acp_res)

library(FactoMineR)
library(factoextra)

res.pca <- PCA(df, scale.unit = TRUE, graph = FALSE)

# Cercle des corrélations
fviz_pca_var(res.pca, 
             col.var = "cos2", # Colorier selon la qualité de représentation
             gradient.cols = c("blue", "yellow", "red"),
             repel = TRUE # Évite le chevauchement des noms
)


# données de test
ts = scale(test_hour[,c(8:11, 14:17, 20:23, 25:28)], center=acp_res$center, scale=acp_res$scale)
ts_red = ts %*% acp_res$rotation

data_train = cbind(train_hour, pc1 = acp_res$x[,1])
data_test = cbind(test_hour, pc1 = ts_red[, 1])
data = rbind(data_train, data_test)

# DataFrame avec les scores de l'ACP
pc_df <- data.frame(
  Index = dates,
  PC1 = data[,'pc1']
)
# DataFrame avec les couleurs
color_df <- data.frame(
  Index = dates,
  Couleur = data[, "Couleur"]
)
# Graphique ACP
p2 <- ggplot(pc_df, aes(x = Index, y = PC1)) +
  geom_line(color = "purple") +
  labs(title = "ACP des variables abs (bleu, vert, jaune, rouge)", x = "Date", y = "PC1") +
  theme_minimal()
# Graphique Couleur
p1 <- ggplot(color_df, aes(x = Index, y = Couleur)) +
  geom_line(color = "tomato") +
  labs(title = "Variable couleur", x = "Date", y = "valeur") +
  theme_minimal()
# Affichage empilé
p1 / p2

###################
# barplot sur train
df = train_hour[,-c(8:11, 14:17, 20:23, 25:28,4)]
df <- cbind(df, pc1 = acp_res$x[, 1])
# Corrélation avec la composante principale
cor_with_x <- cor(df, use = "complete.obs")[, "pc1"]
# Mise en forme du tableau
cor_df <- data.frame(variable = names(cor_with_x), correlation = cor_with_x) %>%
  filter(variable != "pc1") %>%
  arrange(desc(abs(correlation)))
# Barplot
ggplot(cor_df, aes(x = variable, y = correlation)) +
  geom_col(fill = "purple") +
  coord_flip() +
  geom_text(aes(label = round(correlation, 2)), hjust = ifelse(cor_df$correlation > 0, -0.1, 1.1)) +
  theme_minimal() +
  labs(title = "Corrélation avec la variable pc1",
       x = "Variables",
       y = "Corrélation") +
  ylim(c(-1, 1))


# barplot sur test
df = test_hour[,-c(8:11, 14:17, 20:23, 25:28,4)]
df <- cbind(df, pc1 = ts_red[, 1])
# Corrélation avec la composante principale
cor_with_x <- cor(df, use = "complete.obs")[, "pc1"]
# Mise en forme du tableau
cor_df <- data.frame(variable = names(cor_with_x), correlation = cor_with_x) %>%
  filter(variable != "pc1") %>%
  arrange(desc(abs(correlation)))
# Barplot
ggplot(cor_df, aes(x = variable, y = correlation)) +
  geom_col(fill = "purple") +
  coord_flip() +
  geom_text(aes(label = round(correlation, 2)), hjust = ifelse(cor_df$correlation > 0, -0.1, 1.1)) +
  theme_minimal() +
  labs(title = "Corrélation avec la variable pc1",
       x = "Variables",
       y = "Corrélation") +
  ylim(c(-1, 1))

#############################################################
# interpolation de l'absorbance
# le spectre d'absorbance
lambda <- c(465, 525, 595, 625)
# vecteur de longueurs d'onde cibles
lambda_interp <- seq(400, 700, by=1)
# 180
# exemple d'interpolation
abs_interp <- spline(x = lambda, y = data[20,25:28], xout = lambda_interp)$y
# le minimum de la courbe interpolée
min_index <- which.min(abs_interp)
lambda_min <- lambda_interp[min_index]
abs_min <- abs_interp[min_index]
# data.frame long pour ggplot
df_interp <- data.frame(lambda = lambda_interp, absorbance = abs_interp)
df_points <- data.frame(lambda = lambda, absorbance = data[20,25:28])
df_min <- data.frame(lambda = lambda_min, absorbance = abs_min)
# représentation graphique
ggplot() +
  geom_line(data = df_interp, aes(x = lambda, y = absorbance, color = "Interpolé"), size = 1) +
  geom_point(data = df_points, aes(x = lambda, y = absorbance, color = "Mesuré"), size = 3) +
  geom_point(data = df_min, aes(x = lambda, y = absorbance, color = "Minimun"), size = 4, shape = 8) +  # ★ point au minimum
  scale_color_manual(values = c("Interpolé" = "purple", "Mesuré" = "red", "Minimun" = "black")) +
  labs(
    x = "Longueur d'onde (nm)",
    y = "Absorbance",
    color = "",
    title = "Spectre interpolé vs valeurs mesurées"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top")
#########################################
# recherche pour toutes les observations
coul = c()
for (i in 1:nrow(data)) {
  abs_interp <- spline(x = lambda, y = data[i,25:28], xout = lambda_interp)$y
  coul = rbind(coul, lambda_interp[which.min(abs_interp)])
}
# correction de l'erreur
which(coul==700)
coul[589]=mean(c(coul[c(588,590)]))
# conversion des longueur d'onde en couleur
get_rgb <- function(wl) {
  data("ciexyz31", package = "colorscience")
  if (!(wl %in% ciexyz31$wlnm)) return(NA)
  # Extraire XYZ pour cette longueur d’onde
  xyz <- ciexyz31[ciexyz31$wlnm == wl, 2:4]
  # Conversion en sRGB
  rgb_vals <- XYZ2RGB(as.numeric(xyz))
  # S'assurer que les valeurs sont dans [0,1]
  rgb_vals <- pmin(pmax(rgb_vals, 0), 1)
  # Convertir en couleur R
  rgb(rgb_vals[1], rgb_vals[2], rgb_vals[3])
}
# dataframe pour le graphique
df = data.frame(coul)
df$x = dates
df$couleur = sapply(coul, get_rgb)
# représentation
ggplot(df, aes(x = x, y = coul)) +
  geom_line(color = "darkgray") +
  geom_point(aes(color = couleur), size = 1) +
  scale_color_identity() +
  labs(title = "Prévision de la couleur avec abs180",
       x = "Date",
       y = "longeur d'onde") +
  theme_minimal()




df = train_hour[,-c(8:11, 14:17, 20:23, 25:28,4)]
df <- cbind(df, coul = coul[1:144])
# Corrélation avec la composante principale
cor_with_x <- cor(df, use = "complete.obs")[, "coul"]
# Mise en forme du tableau
cor_df <- data.frame(variable = names(cor_with_x), correlation = cor_with_x) %>%
  filter(variable != "coul") %>%
  arrange(desc(abs(correlation)))
# Barplot
ggplot(cor_df, aes(x = variable, y = correlation)) +
  geom_col(fill = "purple") +
  coord_flip() +
  geom_text(aes(label = round(correlation, 2)), hjust = ifelse(cor_df$correlation > 0, -0.1, 1.1)) +
  theme_minimal() +
  labs(title = "Corrélation avec la variable couleur à partir de abs180",
       x = "Variables",
       y = "Corrélation") +
  ylim(c(-1, 1))



df = test_hour[,-c(8:11, 14:17, 20:23, 25:28,4)]
df <- cbind(df, coul = coul[145:696])
# Corrélation avec la composante principale
cor_with_x <- cor(df, use = "complete.obs")[, "coul"]
# Mise en forme du tableau
cor_df <- data.frame(variable = names(cor_with_x), correlation = cor_with_x) %>%
  filter(variable != "coul") %>%
  arrange(desc(abs(correlation)))
# Barplot
ggplot(cor_df, aes(x = variable, y = correlation)) +
  geom_col(fill = "purple") +
  coord_flip() +
  geom_text(aes(label = round(correlation, 2)), hjust = ifelse(cor_df$correlation > 0, -0.1, 1.1)) +
  theme_minimal() +
  labs(title = "Corrélation avec la variable couleur à partir de abs180",
       x = "Variables",
       y = "Corrélation") +
  ylim(c(-1, 1))


###############################################
# 90
# exemple d'interpolation
abs_interp <- spline(x = lambda, y = data[20,14:17], xout = lambda_interp)$y
# le minimum de la courbe interpolée
min_index <- which.min(abs_interp)
lambda_min <- lambda_interp[min_index]
abs_min <- abs_interp[min_index]
# data.frame long pour ggplot
df_interp <- data.frame(lambda = lambda_interp, absorbance = abs_interp)
df_points <- data.frame(lambda = lambda, absorbance = data[20,14:17])
df_min <- data.frame(lambda = lambda_min, absorbance = abs_min)
# représentation graphique
ggplot() +
  geom_line(data = df_interp, aes(x = lambda, y = absorbance, color = "Interpolé"), size = 1) +
  geom_point(data = df_points, aes(x = lambda, y = absorbance, color = "Mesuré"), size = 3) +
  geom_point(data = df_min, aes(x = lambda, y = absorbance, color = "Minimun"), size = 4, shape = 8) +  # ★ point au minimum
  scale_color_manual(values = c("Interpolé" = "purple", "Mesuré" = "red", "Minimun" = "black")) +
  labs(
    x = "Longueur d'onde (nm)",
    y = "Absorbance",
    color = "",
    title = "Spectre interpolé vs valeurs mesurées"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top")
########################################
# recherche pour toutes les observations
coul = c()
for (i in 1:nrow(data)) {
  abs_interp <- spline(x = lambda, y = data[i,14:17], xout = lambda_interp)$y
  coul = rbind(coul, lambda_interp[which.min(abs_interp)])
}
# conversion des longueur d'onde en couleur
# dataframe pour le graphique
df = data.frame(coul)
df$x = seq.POSIXt(from = as.POSIXct("2022-05-31 00:00"), by = "hour", length.out = 696)
df$couleur = sapply(coul, get_rgb)
# représentation
ggplot(df, aes(x = x, y = coul)) +
  geom_line(color = "darkgray") +
  geom_point(aes(color = couleur), size = 1) +
  scale_color_identity() +
  labs(title = "Prévision de la couleur avec abs90",
       x = "Date",
       y = "longeur d'onde") +
  ylim(555,575)+
  theme_minimal()

df = train_hour[,-c(8:11, 14:17, 20:23, 25:28,4)]
df <- cbind(df, coul = coul[1:144])
# Corrélation avec la composante principale
cor_with_x <- cor(df, use = "complete.obs")[, "coul"]
# Mise en forme du tableau
cor_df <- data.frame(variable = names(cor_with_x), correlation = cor_with_x) %>%
  filter(variable != "coul") %>%
  arrange(desc(abs(correlation)))
# Barplot
ggplot(cor_df, aes(x = variable, y = correlation)) +
  geom_col(fill = "purple") +
  coord_flip() +
  geom_text(aes(label = round(correlation, 2)), hjust = ifelse(cor_df$correlation > 0, -0.1, 1.1)) +
  theme_minimal() +
  labs(title = "Corrélation avec la variable couleur à partir de abs90",
       x = "Variables",
       y = "Corrélation") +
  ylim(c(-1, 1))



df = test_hour[,-c(8:11, 14:17, 20:23, 25:28,4)]
df <- cbind(df, coul = coul[145:696])
# Corrélation avec la composante principale
cor_with_x <- cor(df, use = "complete.obs")[, "coul"]
# Mise en forme du tableau
cor_df <- data.frame(variable = names(cor_with_x), correlation = cor_with_x) %>%
  filter(variable != "coul") %>%
  arrange(desc(abs(correlation)))
# Barplot
ggplot(cor_df, aes(x = variable, y = correlation)) +
  geom_col(fill = "purple") +
  coord_flip() +
  geom_text(aes(label = round(correlation, 2)), hjust = ifelse(cor_df$correlation > 0, -0.1, 1.1)) +
  theme_minimal() +
  labs(title = "Corrélation avec la variable couleur à partir de abs90",
       x = "Variables",
       y = "Corrélation") +
  ylim(c(-1, 1))




