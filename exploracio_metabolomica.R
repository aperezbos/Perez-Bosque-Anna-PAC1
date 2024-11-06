library(SummarizedExperiment)

# Carregar el dataset
dataset <- read.csv("human_cachexia.csv", header = TRUE)

# Crear el DataFrame de metadades
metadades <- DataFrame(
  PatientID = dataset$Patient.ID,
  MuscleLoss = dataset$Muscle.loss
)

# Convertir les dades dels metabòlits en una matriu
data <- as.matrix(dataset[, -c(1, 2)])  # Exclou les columnes 1 i 2
data <- t(data)  # Transposa la matriu per a que les mostres siguin columnes

dim(metadades)  # Ha de coincidir amb el nombre de columnes de data
dim(data)       # Comprova el nombre de files i columnes

# Crea l'objecte SummarizedExperiment i guardar-lo
se <- SummarizedExperiment(assays = list(counts = data), colData = metadades)

save(se, file = "dades_metabolomica.Rda")

library(dplyr)
library(tidyr)

# Obtenir els noms dels metabòlits
metabolites <- rownames(assay(se))

# Inicialitzar un dataframe per guardar els resultats de la mitjana i la desviació estàndard
summary_results <- data.frame(metabolite = character(), group = character(), mean = numeric(), sd = numeric(), stringsAsFactors = FALSE)

# Obtenir els valors per a cada metabòlit i calcular les mitjanes i desviacions estàndards per a cada grup (cachexic vs control)
for (metabolite in metabolites) {
  values <- assay(se)[metabolite, ]
  group <- colData(se)$MuscleLoss
  
  # Calcular la mitjana i desviació estàndard per a cada grup
  mean_cachexic <- round(mean(values[group == "cachexic"], na.rm = TRUE),2)
  sd_cachexic <- round(sd(values[group == "cachexic"], na.rm = TRUE),2)
  mean_control <- round(mean(values[group == "control"], na.rm = TRUE),2)
  sd_control <- round(sd(values[group == "control"], na.rm = TRUE),2)
  
  # Emmagatzemar els resultats al dataframe
  summary_results <- rbind(summary_results, 
                           data.frame(metabolite = metabolite, group = "cachexic", mean = mean_cachexic, sd = sd_cachexic),
                           data.frame(metabolite = metabolite, group = "control", mean = mean_control, sd = sd_control))
}

# Convertir el dataframe a un format ampli per mostrar la mitjana i la desviació estàndard per cada grup en columnes separades
summary_table <- summary_results %>%
  pivot_wider(names_from = group, values_from = c(mean, sd))

# Mostrar la taula resum
print(summary_table)

library(factoextra)
library(ggplot2)

# Aplicar l'anàlisi de PCA
pca <- prcomp(t(assay(se)), center = TRUE, scale. = TRUE)

# Crear un vector de colors segons el grup (cachexic vs control)
group_labels <- colData(se)$MuscleLoss
group_colors <- ifelse(group_labels == "cachexic", "#c60000", "#08ee47")

# Visualitzar el PCA amb factoextra (components principals 1 i 2)
fviz_pca_ind(pca, geom = "point", col.ind = group_labels, 
             palette = c("#c60000", "#08ee47"),
             addEllipses = TRUE, ellipse.level = 0.95,
             legend.title = "Grup",
             title = "Anàlisi de Components Principals (PCA)") +
  labs(x = "PC1 (40.4%)", y = "PC2 (8.2%)") +
  theme_minimal()

# Calcular la contribució de cada metabòlit a les components principals
contrib_pca <- as.data.frame(pca$rotation)
contrib_pca$metabolite <- rownames(contrib_pca)

# Filtrar només els metabòlits amb contribució significativa a PC1 o PC2 (superior a 0.2)
threshold <- 0.2
significant_metabolites <- contrib_pca %>%
  select(metabolite, PC1, PC2) %>%  # Seleccionem només les columnes de PC1 i PC2
  filter(abs(PC1) > threshold | abs(PC2) > threshold) %>%  # Filtre per contribució significativa
  mutate(across(starts_with("PC"), ~ round(.x, 3)))  # Arrodonir a 3 decimals

# Visualitzar el PCA de contribució dels metabòlits amb només els més significatius
fviz_pca_var(pca, col.var = "contrib",
             select.var = list(name = significant_metabolites$metabolite),
             gradient.cols = c("#2697f5", "#ffd733", "#ff6b33"),
             repel = TRUE, title = "Contribució dels Metabòlits") +
  labs(x = "PC1 (40.4%)", y = "PC2 (8.2%)") +
  theme_minimal()

# Visualitzar la taula amb els metabòlits significatius (contribució > 0.2 en PC1 o PC2)
print(significant_metabolites)

library(ggplot2)
library(dendextend)

# Matriu de distàncies euclidianes
dist_matrix <- dist(t(assay(se)))

# Clúster jeràrquic amb mètode d'enllaç complet
hc <- hclust(dist_matrix, method = "complete")
dend <- as.dendrogram(hc)

colors <- ifelse(colData(se)$MuscleLoss == "cachexic", "#c60000", "#08ee47")

# Assignar els colors a les etiquetes del dendrograma
labels_colors(dend) <- colors

# Visualitzar el dendrograma amb els colors aplicats
plot(dend, main = "Dendrograma de les mostres",
     ylab = "Distància", xlab = "Mostres")
legend("topright", legend = c("Cachexic", "Control"),
       fill = c("#c60000", "#08ee47"), bty = "n", title = "Grups")

# Aplicar k-means amb 2 grups (cachexic i control)
library(factoextra)
library(ggplot2)

# Aplicar k-means amb 2 grups (cachexic i control)
set.seed(123)
kmeans_result <- kmeans(t(assay(se)), centers = 2)

# Crear un vector amb la informació dels grups coneguts (cachexic vs control)
group_labels <- colData(se)$MuscleLoss

# Visualitzar els resultats del k-means amb fviz_cluster
fviz_plot <- fviz_cluster(kmeans_result, data = t(assay(se)), geom = "point", ellipse.type = "convex",
                          main = "Clúster k-means de les mostres", 
                          palette = c("#ff6b33", "#2697f5"),
                          ggtheme = theme_minimal(), show.clust.cent = FALSE)

# Convertir el plot a un objecte ggplot per a una personalització posterior
fviz_data <- fviz_plot$data
fviz_data$group <- group_labels

# Modificar el plot per afegir formes plenes o buides segons el grup conegut
final_plot <- ggplot() +
  # Dibuixar les el·lipses (convex hulls) només amb el contorn
  stat_ellipse(data = fviz_data, aes(x = x, y = y, group = cluster, color = factor(cluster)),
               geom = "polygon", alpha = 1, fill = NA) +  # fill = NA per deixar sense farciment
  # Afegir els punts amb formes segons el grup conegut
  geom_point(data = fviz_data, aes(x = x, y = y, color = factor(cluster), shape = group), size = 3) +
  scale_color_manual(values = c("1" = "#ff6b33", "2" = "#2697f5")) +
  scale_shape_manual(values = c("cachexic" = 16, "control" = 1)) +  # Ple per cachexic, buit per control
  labs(title = "Clúster k-means de les mostres",
       x = "PC1", y = "PC2", color = "Clúster", shape = "Grup") +
  theme_minimal() +
  theme(legend.position = "right")

# Mostrar el gràfic final
print(final_plot)

# Inicialitzar un dataframe per guardar els resultats del Mann-Whitney U Test
mw_test_results <- data.frame(metabolite = character(), p_value = numeric(), stringsAsFactors = FALSE)

# Realitzar el Mann-Whitney U Test per a cada metabòlit
for (metabolite in metabolites) {
  # Obtenir els valors del metabòlit per a cada grup (cachexic vs control)
  values <- assay(se)[metabolite, ]
  group <- colData(se)$MuscleLoss
  
  # Realitzar el Mann-Whitney U Test amb exact = FALSE per evitar l'avís
  mw_test <- wilcox.test(values ~ group, exact = FALSE)
  
  # Emmagatzemar el nom del metabòlit i el p-valor arrodonit a 3 decimals en el dataframe
  mw_test_results <- rbind(mw_test_results, 
                           data.frame(metabolite = metabolite, p_value = round(mw_test$p.value, 3)))
}

# Visualitzar els resultats
print(mw_test_results)