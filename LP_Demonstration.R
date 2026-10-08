# VAR ESTRUCTURALES Y LOCAL PROJECTIONS CON tidyMacro ----------------------
# Script independiente. Usa solo los datos que trae el paquete, replicando
# los ejemplos de Ciftci, "tidyMacro: A Fast, Tidy Toolkit for Applied
# Macroeconometrics". Los nombres de las funciones siguen la versión 1.0.0
# de GitHub, que difiere del paper (por ejemplo fPlotIRFBQ en vez de
# fplotirf_bq).
#
# Instalación (el paquete está en GitHub, no en CRAN):
#   install.packages("remotes")
#   remotes::install_github("muhsinciftci/tidyMacro")

library(tidyverse)
library(tidyMacro)
library(zoo)       # para convertir trimestres ("1985Q1") en fechas

# Tema de los gráficos del paquete (usa theme_set, funciona con cualquier
# versión de ggplot2). Si la fuente Palatino no está instalada, puede haber
# avisos, pero los gráficos se dibujan igual.
fSetTidyMacroTheme()

# Réplicas de bootstrap para los SVAR. El paper usa 1000; con 500 corre más
# rápido en clase y las bandas cambian poco.
n_boot    <- 500
n_threads <- 3     # núcleos para el bootstrap en paralelo

###########################################################################
# PARTE A: LOCAL PROJECTIONS
###########################################################################

###########################################################################
# 0. IDEA DE LAS LOCAL PROJECTIONS
###########################################################################
# Un VAR estima la dinámica de un período y la itera hacia adelante.
# Una LP estima directamente una regresión por cada horizonte:
#
#   y_{t+h} = a_h + b_h * shock_t + controles_t + error_{t+h},  h = 0, 1, ..., H
#
# El impulso-respuesta es la secuencia b_0, b_1, ..., b_H.
# - Ventaja: no hay que especificar toda la dinámica del sistema, así que es
#   más robusta a errores de especificación.
# - Costo: estimaciones más ruidosas (bandas más anchas) que un VAR.
# - Los errores se superponen entre horizontes, por eso se usan errores
#   estándar Newey-West.

###########################################################################
# 1. LP-OLS: RESPUESTA DEL IPC A UN SHOCK MONETARIO DE ROMER-ROMER
###########################################################################
# El shock narrativo de Romer-Romer se trata como directamente observable:
# entra a la derecha de la regresión como cualquier variable.

data("JordaTaylor2025")

jt_lp <- JordaTaylor2025 |>
  mutate(Date = as.Date(as.yearqtr(Date, format = "%YQ%q"))) |>   # "1985Q1" -> fecha
  filter(Date >= "1984-10-01", Date <= "2007-10-01") |>            # Gran Moderación
  transmute(Date,
            CPI = 100 * logCPI,       # log del IPC x 100 -> respuesta en %
            RRShock,                  # shock monetario narrativo
            dLogGDP,                  # controles: crecimiento del PIB,
            dLogCPI,                  #            inflación,
            dSTIR)                    #            cambio en la tasa corta

head(jt_lp)

controls <- c("dLogGDP", "dLogCPI", "dSTIR")

lp_jt <- fLP(
  CPI ~ RRShock + l(..controls, 1:4),   # l(x, 1:4) crea 4 rezagos de cada control;
  # ..controls usa el vector definido arriba
  data       = jt_lp,
  horizons   = 17,                      # trimestres hacia adelante
  shock      = "RRShock",               # el término cuyo coeficiente es la IRF
  conf       = c(68, 95),               # niveles de confianza en porcentaje
  nw_lags    = 0,                       # ancho de banda Newey-West:
  nw_offset  = 0,                       # regla clásica de Jordà
  cumulative = TRUE
)

fPlotLP(lp_jt) +
  labs(title = "Respuesta del IPC a un shock de Romer-Romer",
       x = "Trimestres", y = "Porcentaje")

# Tabla de resultados. Se llama al método directamente porque el paquete
# exporta tidy.fLP() como función.
tidy.fLP(lp_jt)

###########################################################################
# 2. LP-IV: DESEMPLEO ANTE UN SHOCK EN LA TASA DE FONDOS FEDERALES
###########################################################################
# Ahora la pregunta es el efecto de la tasa de la Fed. Pero la tasa es
# endógena: la Fed la mueve en respuesta a la economía. Por eso el shock de
# Romer-Romer se usa como INSTRUMENTO para la tasa, y cada horizonte se
# estima por mínimos cuadrados en dos etapas (2SLS).

data("JordaTaylor2025IV")

jt_iv <- JordaTaylor2025IV |>
  filter(Date >= as.POSIXct("1985-01-01"),
         Date <= as.POSIXct("2000-01-01")) |>
  transmute(Date,
            urate = Unemployment,   # outcome: tasa de desempleo
            infl  = Inflation,
            ffr   = FFRates,        # tratamiento endógeno: tasa de fondos federales
            rz    = RRShock)        # instrumento

head(jt_iv)

iv_ctrl <- c("urate", "infl", "ffr")

lp_iv <- fLPIV(
  urate ~ ffr + l(..iv_ctrl, 1:6),   # 6 rezagos de cada control
  instruments = ~ l(rz, 0:6),        # instrumento: shock de hoy y 6 rezagos
  data        = jt_iv,
  endog       = "ffr",               # la variable instrumentada
  horizons    = 48,                  # meses
  conf        = c(68, 95),
  nw_lags_iv  = 6,                   # ancho de banda Newey-West fijo en 6
  cumulative  = TRUE
)

fPlotLP(lp_iv) +
  labs(title    = "Desempleo ante un shock de 1 p.p. en la tasa de la Fed (LP-IV)",
       subtitle = "Instrumento: shock de Romer-Romer (0 a 6 rezagos)",
       x = "Meses", y = "Puntos porcentuales")

tidy.fLP(lp_iv)   # fLPIV hereda de fLP, así que el mismo método funciona

# Resultado esperado: respuesta con forma de joroba. El desempleo sube
# gradualmente hasta cerca de 1 p.p. a los ~2.5 años y luego revierte.

###########################################################################
# 3. DIAGNÓSTICO: ¿EL INSTRUMENTO ES FUERTE?
###########################################################################
# F de la primera etapa en cada horizonte. Regla práctica: F < 10 indica un
# instrumento débil, y en ese caso las estimaciones 2SLS no son confiables.

tibble(horizon = 0:48, F = as.numeric(lp_iv$Fstat_fs)) |>
  ggplot(aes(x = horizon, y = F)) +
  geom_col(alpha = 0.75) +
  geom_hline(yintercept = 10, linetype = "dashed") +
  labs(title = "F de la primera etapa por horizonte",
       x = "Meses", y = "F")

###########################################################################
# PARTE B: VAR ESTRUCTURALES (SVAR)
###########################################################################
# Los cuatro ejemplos siguen el mismo flujo:
#   1. fVAR()           -> estima el VAR en forma reducida
#   2. fWoldIRF()       -> respuestas a las sorpresas sin identificar (Wold)
#   3. identificación   -> convierte sorpresas en shocks estructurales
#   4. fBootstrap*()    -> bandas de confianza (90 % y 68 % por defecto)
#   5. fPlotIRF*()      -> gráfico (objeto ggplot)
# Lo único que cambia entre ejemplos es el paso 3: CÓMO se identifica el shock.

###########################################################################
# 4. CHOLESKY: SHOCKS DE INCERTIDUMBRE (BLOOM, 2009)
###########################################################################
# Pregunta: ¿qué le pasa a la economía cuando sube la incertidumbre?
# VAR mensual de 9 variables (1962-2013), 12 rezagos.
# Identificación: Cholesky con la incertidumbre ORDENADA PRIMERO, es decir,
# la incertidumbre no reacciona a las demás variables dentro del mes.

data("Bloom2009")
head(Bloom2009)

y_bloom     <- Bloom2009 |> select(-Date) |> as.matrix()   # una columna por variable
names_bloom <- colnames(y_bloom)
shock_bloom <- match("UNCERT", names_bloom)                # posición del shock

var_bloom <- fVAR(y_bloom, p = 12, c = 1)    # c = 1: incluye constante
fCheckStability(var_model = var_bloom)       # equivalente a roots() en vars

wold_bloom <- fWoldIRF(var_bloom, horizon = 48)   # 48 meses

# Factor de Cholesky: chol() devuelve la triangular SUPERIOR, por eso se
# transpone con t() para obtener la INFERIOR que usa la identificación.
S_bloom   <- t(chol(var_bloom$sigma))
irf_bloom <- fCholeskyIRF(wold_bloom, S_bloom)

boot_bloom <- fBootstrapChol(
  y          = y_bloom,
  var_result = var_bloom,
  nboot      = n_boot,
  horizon    = 48,
  bootscheme = "wild",       # bootstrap "wild": robusto a heterocedasticidad
  n_threads  = n_threads
)

fPlotIRFChol(
  point       = irf_bloom,
  boot_result = boot_bloom,
  shock       = shock_bloom,
  varnames    = names_bloom,
  facet_ncol  = 3
) +
  labs(title = "Respuestas a un shock de incertidumbre (Bloom, 2009)", y = NULL)

# FEVD: qué parte de la varianza de cada variable explica la incertidumbre
fevd_bloom <- fFEVDChol(irf_bloom, shock = shock_bloom)
fPlotVarDec(fevd       = fevd_bloom$fevd,
            varnames   = names_bloom,
            shocknames = "UNCERT")

# Resultado esperado: caída rápida del empleo y la producción industrial,
# seguida de una recuperación con "rebote" (forma de V con sobrerreacción).

###########################################################################
# 5. RESTRICCIONES DE LARGO PLAZO: SHOCKS TECNOLÓGICOS (GALÍ, 1999)
###########################################################################
# Pregunta: ¿qué le pasa a las horas trabajadas ante un shock tecnológico?
# VAR bivariado: crecimiento de la productividad y de las horas.
# Identificación (Blanchard-Quah): solo el shock tecnológico puede tener un
# efecto PERMANENTE sobre el nivel de la productividad. En vez de restringir
# el efecto del mismo período (Cholesky), se restringe el efecto de largo plazo.

data("Gali1999")
head(Gali1999)

y_gali <- Gali1999 |> select(Productivity, Hours) |> as.matrix()

var_gali  <- fVAR(y_gali, p = 4, c = 1)
wold_gali <- fWoldIRF(var_gali, horizon = 40)

# Matriz de impacto de Blanchard-Quah, en tres pasos:
C1 <- apply(wold_gali, c(1, 2), sum)                    # efecto acumulado de largo plazo
D1 <- t(chol(C1 %*% var_gali$sigma %*% t(C1)))          # Cholesky en el LARGO plazo
K  <- solve(C1, D1)                                     # impacto en el mismo período

irf_gali <- fBQIRF(wold_gali, K)

boot_gali <- fBootstrapBQ(
  y          = y_gali,
  var_result = var_gali,
  nboot      = n_boot,
  horizon    = 40,
  bootscheme = "residual",
  cumulate   = c(1, 2)      # acumular: de tasas de crecimiento a niveles
)

fPlotIRFBQ(
  point       = irf_gali,
  boot_result = boot_gali,
  varnames    = c("Productividad", "Horas"),
  shocknames  = c("Shock tecnológico", "Shock no tecnológico"),
  cumulate    = c(1, 2),    # igual que en el bootstrap
  shocks      = c(1, 2)     # qué shocks graficar (1 = tecnológico, 2 = no tecnológico)
)

# Resultado esperado: ante un shock tecnológico positivo la productividad sube
# de forma persistente y las horas CAEN en el corto plazo, en contra de lo que
# predice el modelo clásico de ciclo real.

###########################################################################
# 6. INSTRUMENTO EXTERNO (PROXY-SVAR): NOTICIAS DE OFERTA DE PETRÓLEO
#    (KÄNZIG, 2021)
###########################################################################
# Pregunta: ¿qué efecto tiene un shock de noticias sobre la oferta de petróleo?
# Identificación: un instrumento externo, el cambio en los futuros del petróleo
# en una ventana estrecha alrededor de los anuncios de la OPEP. Esa variación
# se debe a la noticia y no a otras cosas que pasan en la economía.
# El instrumento debe cumplir: relevancia (correlacionado con el shock de
# petróleo) y exogeneidad (no correlacionado con los demás shocks).

data("Kaenzig2021")

y_vars <- c("Oil_Price", "World_Oil_Prod", "World_Oil_Inven",
            "World_IP", "US_IP", "US_CPI")
names_kan <- c("Precio del petróleo", "Producción mundial de petróleo",
               "Inventarios mundiales", "Producción industrial mundial",
               "Producción industrial EE. UU.", "IPC EE. UU.")

y_kan  <- Kaenzig2021 |> select(all_of(y_vars)) |> as.matrix()
iv_oil <- Kaenzig2021 |> select(iv_kanzig_final) |> drop_na() |> as.matrix()

var_kan <- fVAR(y_kan, p = 12, c = 1)

# El instrumento es más corto que la muestra del VAR. Hay que indicar qué
# filas del instrumento (adjustZ) se alinean con qué filas de los residuos
# (adjustu): el instrumento corresponde al FINAL de la muestra.
n_res   <- nrow(var_kan$residuals)
adjustZ <- c(1, nrow(iv_oil))
adjustu <- c(n_res - nrow(iv_oil) + 1, n_res)

# Bootstrap por bloques (moving-block): remuestrea bloques de 24 meses de
# residuos e instrumento juntos, para conservar su dependencia temporal.
boot_kan <- fBootstrapIVMBB(
  y          = y_kan,
  var_result = var_kan,
  Z          = iv_oil,
  nboot      = n_boot,
  blocksize  = 24,
  adjustZ    = adjustZ,
  adjustu    = adjustu,
  policyvar  = 1,            # el shock se normaliza sobre la variable 1 (precio del petróleo)
  horizon    = 48,
  n_threads  = n_threads
)

fPlotIRFIV(
  result_mbb = boot_kan,
  varnames   = names_kan,
  shockname  = "Shock de noticias de petróleo",
  scale      = 10,           # shock que sube el precio del petróleo 10 %
  facet_ncol = 3
) +
  labs(x = "Meses", y = NULL)

# Resultado esperado: un recorte de producción de la OPEP sube el precio del
# petróleo y reduce la actividad mundial y de EE. UU.

###########################################################################
# 7. MÁXIMO EFECTO DE LARGO PLAZO: SHOCKS DE NOTICIAS
#    (BEAUDRY Y PORTIER, 2014)
###########################################################################
# Pregunta: ¿las noticias sobre la tecnología futura mueven la economía hoy?
# VAR trimestral: PTF (productividad total), S&P 500, consumo y horas.
# Identificación: el shock de noticias es el que
#   (a) NO mueve la PTF en el mismo período (es una noticia, todavía no ocurrió), y
#   (b) tiene el MAYOR efecto sobre la PTF en el largo plazo.

data("BeaudryPortier2014")

y_bp     <- BeaudryPortier2014 |> as.matrix()
names_bp <- colnames(y_bp)
names_bp

var_bp  <- fVAR(y_bp, p = 2, c = 1)
wold_bp <- fWoldIRF(var_bp, horizon = 40)
S_bp    <- t(chol(var_bp$sigma))                    # Cholesky inferior

irf_bp <- fMaxIRF(wold_bp, S_bp, var_idx = 1)       # var_idx = 1: la PTF es la variable objetivo

boot_bp <- fBootstrapMax(
  y          = y_bp,
  var_result = var_bp,
  nboot      = n_boot,
  horizon    = 40,
  var_idx    = 1,
  cumulate   = integer(0),   # sin acumular: los datos ya están en niveles
  n_threads  = n_threads
)

S_bp   <- t(chol(var_bp$sigma))                  # Cholesky inferior, no sigma
irf_bp <- fMaxIRF(wold_bp, S_bp, var_idx = 1)

fPlotIRFLR(
  point       = irf_bp,
  boot_result = boot_bp,
  varnames    = colnames(y_bp),
  facet_ncol  = 2
) +
  labs(title = "Respuestas a un shock de noticias", y = NULL)
# Resultado esperado: la bolsa y el consumo suben de inmediato, mientras que
# la PTF sube recién gradualmente. Los mercados anticipan la tecnología futura.
# Variante: fUhligIRF() elige el shock que más explica la VARIANZA de la PTF
# (en vez del efecto de largo plazo), como en Barsky y Sims (2012).

###########################################################################
# RESUMEN: CINCO FORMAS DE IDENTIFICAR UN SHOCK
###########################################################################
# Cholesky (Bloom)          -> orden: la primera variable no reacciona en el mes
# Largo plazo (Galí)        -> solo un shock tiene efectos permanentes
# Instrumento (Känzig)      -> una variable externa que solo refleja el shock
# Máximo efecto (B&P)       -> el shock que más mueve una variable a futuro
# Local projections (J&T)   -> el shock se observa (o se instrumenta) directamente