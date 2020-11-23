#Ejemplo práctico R-Ladies
# Referencies:
# https://www.aridhia.com/blog/r-for-researchers-8-essential-cheatsheets-for-research-data-analysis/
# https://www.datacamp.com/community/tutorials/survival-analysis-R#comments
# https://rpkgs.datanovia.com/survminer/index.html
# http://www.stat.unipg.it/luca/R/


# Cargamos las librerias
# install.packages() si no se han instalado antes las librerias
install.packages("digest") si no se han instalado antes las librerias
library(digest)
library(survival)
# devtools::install_github("kassambara/survminer", build_vignettes = FALSE)
library(survminer)
library(dplyr)

# Importamos la base de datos de cancer de ovario ("ovarian") y le echamos un vistazo
data(ovarian)
glimpse(ovarian)
# Rows: 26
# Columns: 6
# $ futime   <dbl> 59, 115, 156, 421, 431, 448, 464, 475, 477, 563, 638, 744, 769, 770, 803, 855, 1040, 1106, 1129...
# $ fustat   <dbl> 1, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0
# $ age      <dbl> 72.3315, 74.4932, 66.4658, 53.3644, 50.3397, 56.4301, 56.9370, 59.8548, 64.1753, 55.1781, 56.75...
# $ resid.ds <dbl> 2, 2, 2, 2, 2, 1, 2, 2, 2, 1, 1, 1, 2, 2, 1, 1, 2, 1, 1, 2, 1, 2, 2, 1, 2, 1
# $ rx       <dbl> 1, 1, 1, 2, 1, 1, 2, 2, 1, 2, 1, 2, 2, 2, 1, 1, 1, 1, 2, 2, 2, 1, 1, 2, 2, 2
# $ ecog.ps  <dbl> 1, 1, 2, 1, 1, 2, 2, 2, 1, 2, 2, 1, 2, 1, 1, 2, 2, 1, 1, 1, 2, 2, 1, 2, 1, 1
help(ovarian)

# Dicotomizamos la edad y cambiamos las etiquetas
ovarian$rx <- factor(ovarian$rx, 
                     levels = c("1", "2"), 
                     labels = c("A", "B"))
ovarian$resid.ds <- factor(ovarian$resid.ds, 
                           levels = c("1", "2"), 
                           labels = c("no", "si"))
ovarian$ecog.ps <- factor(ovarian$ecog.ps, 
                          levels = c("1", "2"), 
                          labels = c("bueno", "malo"))

# La edad parece ser bimodal
hist(ovarian$age, xlab= 'Años', main = 'Histograma: Edad') 


ovarian <- ovarian %>% mutate(age_group = ifelse(age >=50, "mayor", "joven"))
ovarian$age_group <- factor(ovarian$age_group)

# Kaplan-Meier
surv_object <- Surv(time = ovarian$futime, event = ovarian$fustat)
surv_object 
# [1]   59   115   156   421+  431   448+  464   475   477+  563   638   744+  769+  770+  803+  855+ 1040+ 1106+
# [19] 1129+ 1206+ 1227+  268   329   353   365   377+

fit1 <- survfit(surv_object ~ rx, data = ovarian)
class(fit1)
summary(fit1)

# La funcion ggsurvplot crea graficos ggplot2 a partir de objetos survfit (Graficar las curvas de supervivencia)
ggsurvplot(fit1, data = ovarian, pval = TRUE)
# Añadiendo el argumento "fun" podemos decidir que funcion queremos graficar.
# "event" para eventos acumulados, "cumhaz" para la funcion de riesgo acumulada
# "pct" para la probabilidad de supervivencia en porcentaje
par(mfrow = c(3,1))
ggsurvplot(fit1, data = ovarian, fun = "event")
ggsurvplot(fit1, data = ovarian, fun = "cumhaz")
ggsurvplot(fit1, data = ovarian, fun = "pct")
par(mfrow = c(1,1))

# ggsurvplot tiene muchos parametros graficos para escoger, por ejemplo, intervalos de confianza, 
# mostrar la tabla de individuos a riesgo, la posicion de la leyenda, anotaciones adicionales 
# (como añadir el p-valor del test de logrank), titulo, subtitulo, ...

ggsurv1 <- ggsurvplot(fit1, data = ovarian, conf.int = TRUE, conf.int.style = "step", pval = TRUE, risk.table = TRUE, 
           risk.table.height = 0.25, risk.table.y.text.col = TRUE, risk.table.y.text = FALSE, size = 1, 
           linetype = "strata", palette = c("#E7B800", "#2E9FDF"), legend = "bottom", legend.title = "RX",
           legend.labs = c("A", "B"), censor.shape="|", censor.size = 4, ncensor.plot = TRUE, 
           ncensor.plot.height = 0.25, ggtheme = theme_light(), xlim = c(0,1200), xlab = "Time in days", 
           break.time.by = 100, surv.median.line = "hv")
ggsurv1

# Valores predictivos del status residual de la enfermedad
fit2 <- survfit(surv_object ~ resid.ds, data = ovarian)
ggsurvplot(fit2, data = ovarian, pval = TRUE)

ggsurvplot(fit2, data = ovarian, conf.int = TRUE, pval = TRUE, fun = "pct", risk.table = TRUE, size = 1, 
           linetype = "strata", palette = c("#E7B800", "#2E9FDF"), legend = "bottom", legend.title = "Resid",
           legend.labs = c("no", "si"))

# Ajustamos un modelo de riesgos proporcionales de Cox
fit.coxph <- coxph(surv_object ~ rx + resid.ds + age_group + ecog.ps, 
                   data = ovarian)

summary(fit.coxph)
# La funcion ggforest() de la libreria survminer crea un foresplot para el modelo de Cox ajustado
# Muestra la estimacion de los hazard ratios, intervalos de confianza y p-valores para cada covariable
ggforest(fit.coxph, data = ovarian)


# La funcion cox.zph() de la libreria "survival" se puede usar para verificar la asuncion de riesgos 
# proporcionales de un modelo de Cox ajustado
ftest <- cox.zph(fit.coxph)
ftest

# La funcion ggcxzph() de la libreria "survminer" podemos hacer una verificacion grafica 
 # Produce un grafico del tiempo vs los residuos de Scheonfeld para cada covariable
ggcoxzph(ftest)


# La funcion ggcoxdiagnostics() de la libreria "survminer" grafica diferentes residuos como funciones del tiempo
# predictores lineales o las propias observaciones
# El parametro type permite escoger el tipo de residuos: "martingale", "deviance", "score", "schoenfeld", "dfbeta",
# "dfbetas", "scaledsch"
# El parametro ox.scale permite escoger que graficar en el eje OX: "linear.predictions", "observation.id" o "time" 
# El parametro hline añade una linea horizontal de referencia
# El parametro sline suaviza la curva del grafico (smooth line)
ggcoxdiagnostics(fit.coxph, type = "deviance", ox.scale = "linear.predictions")
ggcoxdiagnostics(fit.coxph, type = "schoenfeld", ox.scale = "time")

# La funcion ggadjustedcurves() de la libreria survimer grafica las curvas de supervivencia ajustadas
# por el modelo de riesgos proporcionales de Cox 
# Estas curvas difieren de las creadas con Kaplan-Meier ya que contienen la informacion del modelo de Cox ajustado
ggadjustedcurves(fit.coxph, data = ovarian)
curve <- surv_adjustedcurves(fit.coxph, data = ovarian)
curve
summary(curve)

surv_adjustedcurves(fit.coxph, data = ovarian, variable = 'age_group')
ggadjustedcurves(fit.coxph, data = ovarian, variable = 'age_group')

surv_adjustedcurves(fit.coxph, data = ovarian, variable = 'rx')
ggadjustedcurves(fit.coxph, data = ovarian, variable = 'rx')


# Ajuste de la funcion de incidencia acumulada

#Competing risks 
set.seed(2)
ss <- rexp(100)
gg <- factor(sample(1:3,100,replace=TRUE),1:3,c('BRCA','LUNG','OV'))
cc <- factor(sample(0:2,100,replace=TRUE),0:2,c('no event', 'death', 'progression'))
strt <- sample(1:2,100,replace=TRUE)

df <- data.frame(time = ss, group = gg, status = cc, strt)
fit2 <- survfit(Surv(time, status, type="mstate") ~ 1, data=df)
ggcompetingrisks(fit2)
fit3 <- survfit(Surv(time, status, type="mstate") ~ group, data=df)
ggcompetingrisks(fit3)

library(ggsci)
library(cowplot)
ggcompetingrisks(fit3) + theme_cowplot() + scale_fill_jco()


library(cmprsk)
# source ("CumIncidence.R")
fit4 <- cuminc(df$time, df$status, cencode = 0)
plot(fit4)
