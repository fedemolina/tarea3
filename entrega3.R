library(dplyr)
library(tidyr)
library(readr)
library(tidyverse)
library(reshape2)
library(Rcpp)
data <- read.delim("datos_rna.txt", stringsAsFactors = FALSE)
str(data)
#data2 <-readr::read_delim("datos_rna.txt", delim = "\t")

######## EDA

data.tidy <- data %>% 
  gather(., "cod", "valor", 2:ncol(.)) %>%
  separate(data=. , col = cod, into = c("rep","geno","alelo")) %>%
  separate(data = . , col = rep, sep = 3, into = c("REP","rep")) %>% select(.,-REP) %>%
  spread(data=., alelo, valor)
data.tidy$geno[data.tidy$geno==c("B","BM","M","MB")] <- c("B73","B73xMo17","Mo17","Mo17xB73")
data.tidy[,4:6] <- lapply(data.tidy[,4:6], as.numeric) #podría usar purrr
data.tidy[2:3] <- lapply(data.tidy[,2:3], as.factor) 
str(data.tidy) #Perfeito

#Obs: Los pen-últimos dos pasos se pueden combinar en uno solo.

data.plot <- data.tidy
data.plot[,3] <- as.character(data.plot[,3])
data.plot$geno[data.plot$geno=="Mo17xB73"] <- "B73xMo17"
#data.plot <- data.plot %>% filter(!geno=="Mo17xB73")
data.plot[,3] <- as.factor(data.plot[,3])

data.plot %>% ggplot(., aes(x=log(m+1),y=log(b+1), color = geno)) +
  geom_point(alpha=1/6) +
  facet_wrap(~rep) +
  geom_abline(slope=1, intercept=0, size = 1.5) +
  # scale_x_log10() +
  # scale_y_log10() +
  labs(x = "Expresión génica de alelo M (en logs)", y = "Expresión génica de alelo B (en logs)") +
  theme_bw() +
  scale_color_manual(values = c("green4","orangered2","purple4"))

# Comentarios figura
  # Comentar las correlaciones

### Ejercicio 2
  # No queda claro si habría que filtrar por "B73xMo17"
  # Lo hago de las dos formas.

compara <- function(x, y) {
  m <- length(x)
  n <- length(y)
  # calculo el estadistico de la prueba
  sp <- sqrt(((m-1)*sd(x)^2 + (n-1)*sd(y)^2) / (m+n-2))
  tstat <- (mean(x) - mean(y)) / (sp*sqrt(1/m + 1/n))
  # calculo el p-valor
  2*(1 - pt( abs(tstat), df = n+m-2) )
  #el valor t es 8.743106
}
Rcpp::sourceCpp("test_medias.cpp")

# a)
BM <- data.tidy %>% filter(geno=="B73xMo17",GeneID == "AC155377.1_FG001") %>% select(b,m,GeneID)
compara(BM$b,BM$m) 
t.test(BM$b,BM$m)$p.value 
test_medias(BM$b,BM$m) 

BM0 <- data.tidy %>% filter(GeneID == "AC155377.1_FG001")
compara(BM$b,BM$m) 
t.test(BM$b,BM$m)$p.value 
test_medias(BM$b,BM$m) 

# b
a <- microbenchmark::microbenchmark(
  compara = compara(BM$b,BM$m),
  'R base' = t.test(BM$b,BM$m)$p.value,
  Rcpp = test_medias(BM$b,BM$m),
  times=10000L
)
a %>% ggplot(., aes(y=log(time+1), x = expr, col = expr)) +
  geom_boxplot() +
  #scale_y_log10() +
  labs(x = "", y = "Tiempo (en logs)", col = "Funciones test de medias") +
  theme(legend.position = "bottom")
  
# c
  # Explicar las diferencias en cada corrida debido a procesamieno interno de la pc

# d
# d. Pensando en la comparación entre alelos para un sólo gen (ejemplo
# el ‘AC155377.1_FG001’), cuál es el problema más grave de hacer la
# prueba t.test()?
    # El problema es que son solo 4 casos. Por lo tanto, la prueba pierde potencia? significatividad?

### Ejercicio 3
# ¿Cuál es la probabilidad de concluir que hay
# al menos un gen en que las medias son distintas (es decir, rechazar H0
#                                                  con alpha = .05) ?
BM3 <- data.tidy %>% filter(geno=="B73xMo17") 
compara(BM3$b,BM3$m)
t.test(BM3$b,BM3$m)
test_medias(BM3$b,BM3$m)
# No esta el cero en el intervalo. P value es cero.
# Por defecto alpha es 0.05
  # No se rechaza H0 con alpha = 0.05. La probabilidad es prácticamente uno?
  # Prob(x>1)
BM3 <- data.tidy %>% filter(geno=="B73xMo17") # %>% group_by(GeneID) %>% summarise(test = t.test)

compara(BM3$b,BM3$m)
t.test(BM3$b,BM3$m)
test_medias(BM3$b,BM3$m)


# b. Asumiendo que: para todos los genes, las medias de m y b son exactamente
# iguales (H0g cierta), y que la expresión de cada gen es
# independiente del resto. ¿En cuantos genes se espera concluir, erróneamente,
# que la diferencia es significativa con alpha = .05?

#   c. Para cada gen, calcula un valor-p para la prueba que compara las medias
# de cada alelo, utiliza la expresión transformada en logaritmos (puedes
# sumar el valor 1 para asegurarte que siempre existe la transformación).
BM3 <- BM3 %>% mutate(logm = log(1+.$m),
                      logb = log(1+.$b))
# Con broom 
BM3 %>% group_by(GeneID) %>% dplyr::do(broom::tidy(t.test(log(.$b+1),log(.$m+1), data=.)))
# R base
lista <- split(BM3, BM3$GeneID,drop=FALSE) 
output <- vector("list", length = length(lista))
output <- lapply(1:length(lista),function (x) output[[x]] <- t.test(lista[[x]]$logb,lista[[x]]$logm)) # NO rechazo H0, el cero incluido en el intervalo.
pvalor <- vector(length = length(lista))
pvalor <- lapply(1:length(lista), function (x) pvalor[[x]] <- output[[x]]$p.value) %>% unlist() 
  # Los NaN se explican porque en la agrupación hay genes de alelos a y b que son iguales, por lo tanto no hay variancia en ningún grupo ni standar error.
  # o sea queda indefinido.
mean(pvalor<0.05, na.rm = TRUE)
summary(pvalor)

# ¿Qué proporción de genes tiene valor-p menor al 5% ?
#   d. Utilizando la función p.adjust() se pueden corregir los valores-p
# debido a que estamos realizando 2500 pruebas a la vez. ¿Qué proporción
# de genes tiene valor-p (ajustado) menos a 5% ?
mean(p.adjust(pvalor)<0.05,na.rm = TRUE)
summary(p.adjust(pvalor))
  # Aproximadamente un 0.3%. Sin contar mas de 1000 NaN
