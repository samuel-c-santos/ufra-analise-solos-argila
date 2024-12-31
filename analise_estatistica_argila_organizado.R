# ----------------------------------
# Configuração inicial
# ----------------------------------
# Configuração do diretório de trabalho
setwd("G:/Meu Drive/UFRA/SOL 24/")

# Carregar pacotes necessários
library(sf)
library(dplyr)
library(psych)
library(ggplot2)
library(raster)
library(gstat)
library(MASS)
library(gt)
library(ggspatial)
library(gridExtra)

# ----------------------------------
# Importação e preparação dos dados
# ----------------------------------
# Importar shapefiles
limite <- st_read("./imovel/area_imovel_utm.shp")
camada <- st_read("./imovel/Analise_solo_completo.shp")

# Verificar e alinhar projeções
if (!st_crs(camada) == st_crs(limite)) {
  limite <- st_transform(limite, st_crs(camada))
}

# Listar colunas disponíveis e primeiras linhas da tabela de atributos
colnames(camada)
head(camada)

# ----------------------------------
# Estatísticas descritivas e remoção de outliers
# ----------------------------------
# Remover geometria para manipulação tabular
variaveis <- camada %>% st_drop_geometry()

# Selecionar e renomear colunas
variaveis <- variaveis[, c("Argila.", "pHCaCl2", "pHÁgua")]
colnames(variaveis) <- c("Argila", "pHCaCl2", "pHÁgua")
variaveis$Argila_original <- variaveis$Argila  # Preservar dados originais

# Remover outliers da variável "Argila"
variaveis <- variaveis %>%
  filter(!Argila %in% boxplot.stats(Argila)$out)

# Estatísticas básicas e detalhadas
summary(variaveis)
describe(variaveis)

# Teste de normalidade antes da transformação
shapiro_original <- shapiro.test(variaveis$Argila_original)
print(shapiro_original)

# ----------------------------------
# Transformação Box-Cox
# ----------------------------------
# Calcular parâmetro lambda para Box-Cox
boxcox_model <- boxcox(lm(Argila ~ 1, data = variaveis))
lambda <- boxcox_model$x[which.max(boxcox_model$y)]

# Aplicar a transformação
variaveis$Argila_boxcox <- (variaveis$Argila^lambda - 1) / lambda

# Teste de normalidade após a transformação
shapiro_boxcox <- shapiro.test(variaveis$Argila_boxcox)
print(shapiro_boxcox)

# Atualizar transformação no objeto "camada" para usar os valores transformados
camada$Argila_boxcox <- variaveis$Argila_boxcox

# ----------------------------------
# Visualização das distribuições
# ----------------------------------
# Histogramas e boxplots antes e depois da transformação
hist_antes <- ggplot(variaveis, aes(x = Argila_original)) +
  geom_histogram(binwidth = 1, fill = "skyblue", color = "black", alpha = 0.7) +
  theme_minimal() +
  labs(title = "Histograma de Argila (Antes da Transformação)", x = "Argila (%)", y = "Frequência")

hist_depois <- ggplot(variaveis, aes(x = Argila_boxcox)) +
  geom_histogram(binwidth = 0.1, fill = "lightgreen", color = "black", alpha = 0.7) +
  theme_minimal() +
  labs(title = "Histograma de Argila (Depois da Transformação)", x = "Argila (Box-Cox)", y = "Frequência")

boxplot_antes <- ggplot(variaveis, aes(x = "Antes", y = Argila_original)) +
  geom_boxplot(fill = "skyblue") +
  theme_minimal() +
  labs(title = "Boxplot de Argila (Antes da Transformação)", y = "Argila (%)", x = "")

boxplot_depois <- ggplot(variaveis, aes(x = "Depois", y = Argila_boxcox)) +
  geom_boxplot(fill = "lightgreen") +
  theme_minimal() +
  labs(title = "Boxplot de Argila (Depois da Transformação)", y = "Argila (Box-Cox)", x = "")

grid.arrange(hist_antes, hist_depois, boxplot_antes, boxplot_depois, ncol = 2)

# ----------------------------------
# Análise variográfica
# ----------------------------------
# Converter camada para formato Spatial
camada_sp <- as(camada, "Spatial")

# Criar semivariograma experimental
variogram_model <- variogram(Argila_boxcox ~ 1, data = camada_sp)

# Ajustar modelo teórico ao semivariograma experimental
fit_model <- fit.variogram(variogram_model, model = vgm("Sph", nugget = 0.1, psill = 0.5, range = 500))

# Exportar semivariograma experimental
png(filename = "./saida/semivariograma_experimental.png", width = 800, height = 600)
plot(variogram_model, main = "Semivariograma Experimental de Argila (Box-Cox)")
dev.off()
cat("Semivariograma experimental exportado para ./saida/semivariograma_experimental.png\n")

# Exportar semivariograma com modelo teórico ajustado
png(filename = "./saida/semivariograma_ajustado.png", width = 800, height = 600)
plot(variogram_model, fit_model, main = "Semivariograma de Argila (Box-Cox)")
dev.off()
cat("Semivariograma ajustado exportado para ./saida/semivariograma_ajustado.png\n")

# ----------------------------------
# Interpolação por krigagem
# ----------------------------------
# Ajustar a grade para a extensão do polígono limite
extent_limite <- extent(st_bbox(limite))
res <- 5  # Resolução da grade
grid <- raster(extent_limite, res = res) %>% as("SpatialPixelsDataFrame")
proj4string(grid) <- st_crs(limite)$proj4string

# Realizar a krigagem usando valores transformados
krig_result <- krige(Argila_boxcox ~ 1, camada_sp, grid, model = fit_model)

# Converter resultado para raster e recortar pelo limite
krig_raster <- rasterFromXYZ(as.data.frame(krig_result, xy = TRUE))
proj4string(krig_raster) <- st_crs(limite)$proj4string
grid_cortado <- mask(krig_raster, limite)

# ----------------------------------
# Cálculo de variabilidade espacial
# ----------------------------------
# Gerar camada de variabilidade
cv_layer <- calc(krig_raster, fun = function(x) sd(x) / mean(x) * 100)

# Ajustar projeção e recortar
if (!st_crs(cv_layer) == st_crs(limite)) {
  cv_layer <- projectRaster(cv_layer, crs = st_crs(limite)$proj4string)
}
cv_layer_cortado <- mask(cv_layer, limite)

# Converter raster para data frame com nome correto para a coluna "layer"
cv_layer_df <- as.data.frame(cv_layer_cortado, xy = TRUE, na.rm = TRUE)
colnames(cv_layer_df) <- c("x", "y", "layer")

# ----------------------------------
# Ajustar limites para centralizar com margem
# ----------------------------------
bbox_limite <- st_bbox(limite)
margem <- 0.2 * max(bbox_limite["xmax"] - bbox_limite["xmin"], bbox_limite["ymax"] - bbox_limite["ymin"])

bbox_limite["xmin"] <- bbox_limite["xmin"] - margem
bbox_limite["xmax"] <- bbox_limite["xmax"] + margem
bbox_limite["ymin"] <- bbox_limite["ymin"] - margem
bbox_limite["ymax"] <- bbox_limite["ymax"] + margem

# ----------------------------------
# Visualizações
# ----------------------------------
# Mapa de variabilidade espacial
variabilidade_map <- ggplot() +
  geom_tile(data = cv_layer_df, aes(x = x, y = y, fill = layer)) +
  scale_fill_gradientn(colors = c("#ffffcc", "#41b6c4", "#0c2c84"), name = "Variabilidade (%)") +
  geom_sf(data = limite, fill = NA, color = "black", size = 1.2) +
  annotation_scale(location = "bl", width_hint = 0.5) +
  annotation_north_arrow(location = "tr", style = north_arrow_fancy_orienteering) +
  theme_minimal() +
  coord_sf(xlim = c(bbox_limite["xmin"], bbox_limite["xmax"]),
           ylim = c(bbox_limite["ymin"], bbox_limite["ymax"]), expand = FALSE) +
  labs(title = "Mapa de Variabilidade Espacial (CV)", x = "Longitude", y = "Latitude")

print(variabilidade_map)

# Exportar mapa de variabilidade espacial
ggsave(filename = "./saida/mapa_variabilidade_espacial.jpg", plot = variabilidade_map, dpi = 300, width = 10, height = 8)

# Mapa de krigagem contínua
krigagem_df <- as.data.frame(grid_cortado, xy = TRUE, na.rm = TRUE)
colnames(krigagem_df) <- c("x", "y", "layer")

# Adicionar intervalos originais à legenda
krigagem_intervals <- c("< 3.42", "3.42 - 3.44", "3.44 - 3.46", "3.46 - 3.48", "> 3.48")
krigagem_colors <- c("#feffd3", "#fdbb84", "#fc8d59", "#d7301f", "#993504")

krigagem_continuo_map <- ggplot() +
  geom_tile(data = krigagem_df, aes(x = x, y = y, fill = layer)) +
  scale_fill_gradientn(colors = krigagem_colors, name = "Interpolação Contínua",
                       breaks = seq(min(krigagem_df$layer, na.rm = TRUE), max(krigagem_df$layer, na.rm = TRUE), length.out = 5),
                       labels = krigagem_intervals) +
  geom_sf(data = limite, fill = NA, color = "black", size = 1.2) +
  annotation_scale(location = "bl", width_hint = 0.5) +
  annotation_north_arrow(location = "tr", style = north_arrow_fancy_orienteering) +
  theme_minimal() +
  coord_sf(xlim = c(bbox_limite["xmin"], bbox_limite["xmax"]),
           ylim = c(bbox_limite["ymin"], bbox_limite["ymax"]), expand = FALSE) +
  labs(title = "Mapa de Krigagem (Contínuo)", x = "Longitude", y = "Latitude")

print(krigagem_continuo_map)

# Exportar mapa de krigagem contínua
ggsave(filename = "./saida/mapa_krigagem_continuo.jpg", plot = krigagem_continuo_map, dpi = 300, width = 10, height = 8)

# Mapa de krigagem reclassificada
reclass_matrix <- matrix(c(-Inf, 3.42, 1, 3.42, 3.44, 2, 3.44, 3.46, 3, 3.46, 3.48, 4, 3.48, Inf, 5), ncol = 3, byrow = TRUE)
krig_raster_reclass <- reclassify(grid_cortado, reclass_matrix)

krig_reclass_df <- as.data.frame(krig_raster_reclass, xy = TRUE, na.rm = TRUE)
colnames(krig_reclass_df) <- c("x", "y", "class")

reclass_labels <- c("< 3.42", "3.42 - 3.44", "3.44 - 3.46", "3.46 - 3.48", "> 3.48")
reclass_colors <- c("#feffd3", "#fdbb84", "#fc8d59", "#d7301f", "#993504")

krigagem_reclass_map <- ggplot() +
  geom_tile(data = krig_reclass_df, aes(x = x, y = y, fill = factor(class))) +
  scale_fill_manual(values = reclass_colors, labels = reclass_labels, name = "Classes") +
  geom_sf(data = limite, fill = NA, color = "black", size = 1.2) +
  annotation_scale(location = "bl", width_hint = 0.5) +
  annotation_north_arrow(location = "tr", style = north_arrow_fancy_orienteering) +
  theme_minimal() +
  coord_sf(xlim = c(bbox_limite["xmin"], bbox_limite["xmax"]),
           ylim = c(bbox_limite["ymin"], bbox_limite["ymax"]), expand = FALSE) +
  labs(title = "Mapa de Krigagem Reclassificada", x = "Longitude", y = "Latitude")

print(krigagem_reclass_map)

# Exportar mapa de krigagem reclassificada
ggsave(filename = "./saida/mapa_krigagem_reclass.jpg", plot = krigagem_reclass_map, dpi = 300, width = 10, height = 8)

# Exportar raster reclassificado
writeRaster(krig_raster_reclass, filename = "./saida/raster_krigagem.tif", format = "GTiff", overwrite = TRUE)
