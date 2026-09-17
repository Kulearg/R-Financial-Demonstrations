library(tidyverse) #Manejo de Datos
library(tstests) #Tests Varios de Series de Tiempo
library(rugarch) #Paquete principal de este codigo, para modelos GARCH, hay alternativas
library(ggplot2) #Graficos
library(quantmod) #Datos financieros
library(tseries) #Paquete Principal de series de tiempo
library(strucchange) #Paquete para quiebres estructurales
library(forecast) #Para modelos ARMA, principalmente para auto.arima()
library(greeks) #Modelado de Opciones, para aplicaciones financieras
library(patchwork) #Juntar varios graficos en una imagen

# 0. BUSCAR DATOS ---------------------------------------------------------
  #Traemos los datos usando el Ticker, con el paquete quantmod
    #Pide al usuario que seleccione un ticker (accion/activo) de finance.yahoo.com
    #Este codigo fue hecho con el Oro (GLD) en mente
    #Otros ejemplos son Futuro del Petroleo (BZ=F), Nvidia (NVDA), etc
  ticker <- as.character(readline(prompt = "Enter Ticker: "))
  getSymbols(ticker, src = "yahoo", auto.assign = FALSE) -> GLD
  #Hacemos un Grafico simple, de la evolucion de precio del activo
  candleChart(GLD,multi.col=TRUE,theme="white")
  #Al grafico de velas le agregamos unas Bandas de Bollinguer, como primera aproximacion a la volatilidad
  addBBands()


# 1. PREPARACION DE DATOS -------------------------------------------------
  #Sacamos valores vacios, de haberlos y convertimos a retornos logaritmicos para el analisis
  gdx_returns <- na.omit(dailyReturn(GLD, type = "log"))
  #Multiplicamos por 100 por temas matematicos en el analisis
  gdx_returns <- gdx_returns * 100
  #Tests de Estacionariedad, por si hay que integrar
  adf.test(gdx_returns)
  kpss.test(gdx_returns)
  pp.test(gdx_returns)
  #ACF/PACF
  acf(gdx_returns)
  pacf(gdx_returns)
  #Graficos de la Serie, Retornos en el tiempo y Distribucion de los retornos (Recordemos las colas pesadas y la ineficiencia de la distribucion normal)
    #Retornos en el Tiempo
  plot(gdx_returns, type = "c", main = "GLD Returns")
    #Distribucion de los resultados
      #Caracteristicas de la serie por distribucion a evaluar, media, desvio, varianza
  retornos <- as.numeric(gdx_returns)
  mu_norm <- mean(retornos)
  sd_norm <- sd(retornos)
  ajuste_std <- fitdist(distribution = "std", x = retornos)
  mu_std <- ajuste_std$pars["mu"]
  sig_std <- ajuste_std$pars["sigma"]
  shape_std <- ajuste_std$pars["shape"] # Grados de libertad
  ajuste_sstd <- fitdist(distribution = "sstd", x = retornos)
  mu_sstd <- ajuste_sstd$pars["mu"]
  sig_sstd <- ajuste_sstd$pars["sigma"]
  skew_sstd <- ajuste_sstd$pars["skew"] # Sesgo
  shape_sstd <- ajuste_sstd$pars["shape"] 
  df_plot <- data.frame(Retornos = retornos)
  
  ggplot(df_plot, aes(x = Retornos)) +
    # Distribución Empírica
    geom_density(aes(color = "Empírica"), size = 1.2, fill = "gray", alpha = 0.3) +
    
    # Distribución Normal
    stat_function(fun = dnorm, 
                  args = list(mean = mu_norm, sd = sd_norm),
                  aes(color = "Normal"), size = 1, linetype = "dashed") +
    
    # Distribución t-Student
    stat_function(fun = ddist, 
                  args = list(distribution = "std", mu = mu_std, 
                              sigma = sig_std, shape = shape_std),
                  aes(color = "t-Student"), size = 1) +
    
    # Distribución Skew-t-Student
    stat_function(fun = ddist, 
                  args = list(distribution = "sstd", mu = mu_sstd, 
                              sigma = sig_sstd, skew = skew_sstd, shape = shape_sstd),
                  aes(color = "Skew-t-Student"), size = 1) +
    
    coord_cartesian(xlim = c(-5, 5)) + # Foco en el área de interés
    scale_color_manual(values = c("Empírica" = "black", 
                                  "Normal" = "blue", 
                                  "t-Student" = "orange", 
                                  "Skew-t-Student" = "red")) +
    labs(title = "GLD: Distribución Incondicional de Retornos",
         x = "Retornos Logarítmicos (%)", y = "Densidad") +
    theme_minimal() +
    theme(legend.position = "bottom", legend.title = element_blank())
  #Aunque el codigo para este grafico es complejo, el punto es simple, ver la asimetria de la serie y lo inadecuado que es el uso de la distribucion normal en la estimacion de los modelos GARCH en Series Financieras
  #Podemos captar gran parte de esa asimetria simplemente cambiando la distribucion utilizada, en particular, la Skew T-Student (sstd)
  #Retorno Promedio (Para contrastar con las salidas de auto.arima)
  avg_return <- mean(gdx_returns)
  cat("Average Return: ", avg_return)

# 2. ESTIMACION ARMA ------------------------------------------------------
  #Buscamos un modelo ARMA(p,q) con auto.arima, usamos BIC de IC
  modelo_arima <- auto.arima(gdx_returns, stationary = TRUE, trace = TRUE, ic = "bic")
  summary(modelo_arima)
  
  #Ordenes del Modelo ARIMA
  orders <- arimaorder(modelo_arima)
  p <- as.numeric(orders["p"])
  q <- as.numeric(orders["q"])
  #Graficamos p,q
  autoplot(modelo_arima)
  #polyroot(orders) #Para verificar sobreparametrizacion, por ahi un modelo mas chico, no aplica con un ARMA(0,0)
  #ARMA(0,0) es esperado para una serie financiera
  #Utilizamos BIC, penaliza mas en series con numeros mas grandes
  #Si utilizamos AIC o AICc, da un ARMA(1,1), con ar1 0.8420, y ma1  -0.8576
  #Eso es un problema, por que se cancelan entre si, verificamos esto usando polyroot
  #Nos darian raices caracteristicas casi exactas, significa que se cancelan, no estaban aportando nada real, por ende, BIC es la eleccion correcta como IC.
  
# 3. ESTIMACION GARCH -----------------------------------------------------

    #Tipos y Subtipos de modelos GARCH
    #Hacemos esta lista para ilustrar los distintos modelos disponibles
      #Existen muchos otros: ZDGARCH, MCSGARCH, etc
    model_specs <- list(
    list(model = "sGARCH"),
    list(model = "eGARCH"),
    list(model = "gjrGARCH"),
    list(model = "iGARCH"),
    list(model = "csGARCH"),
    list(model = "apARCH"),                             
    list(model = "fGARCH", submodel = "AVGARCH"),     
    list(model = "fGARCH", submodel = "NGARCH"),      
    list(model = "fGARCH", submodel = "NAGARCH"),      
    list(model = "fGARCH", submodel = "ALLGARCH")           
  )
  
  spec_results <- data.frame() #Estructura para luego mostrar los resultados
  
  #Aca un codigo complejo, en resumen, es utilizado para ver las especificaciones de cada uno de los modelos posibles de GARCH, con el orden ARMA visto anteriormente
  #Es un bucle iterativo, el cual prueba cada tipo de modelo (Con la especificacion ugarchspec), ve como se ajusta a los datos (ugarchfit) y le hace ciertas pruebas
  #Las pruebas son las siguientes
    #L Jung Box: Buscamos ver si el modelo dejo heterocedasticidad sin explicar luego que el modelo sea aplicado, es decir, si logran capturar el componente estocastico (Temporal) de la varianza
    #Sign Bias Test: Evaluamos si el modelo reacciona distinto a shocks negativos que a shocks positivos. La asimetria de eso es algo real en las series financieras, y debemos tenerlo en cuenta cuando queremos especificar correctamente
    #Nyblom Test: Vemos si los parametros son estables en el tiempo (Constancia de Parametros), para modelos de maxima verosimilitud, instabilidad general, no quiebres espeficos como Chow o Bai-Perron
      #Importante, solo de demostracion, no hicimos dummies para ciertas crisis, como 2008 o COVID
      #Buen primer filtro de quiebres estructurales 
    #Criterio de Informacion: Parsimonia de los modelos, penalizando cantidad de parametros, recordemos el uso de cada uno para forecasting y modelo "verdadero" visto en clases anteriores.
  #No se preocupen por esto, la estructura es irrelevante, nomas para evitar hacer esto manualmente y poder mostrar de forma concisa todos los modelos disponibles
  alpha_diag <- 0.05 #Nivel de Confianza para los Test, por alguna razon el default es 0.1
  
  for (spec_def in model_specs) {
    label <- if (!is.null(spec_def$submodel)) {
      paste0(spec_def$model, "(", spec_def$submodel, ")")
    } else spec_def$model
    tryCatch({
      vm <- list(model = spec_def$model, garchOrder = c(1, 1))
      if (!is.null(spec_def$submodel)) vm$submodel <- spec_def$submodel
      spec_temp <- ugarchspec(
        variance.model = vm,
        mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),
        distribution.model = "sstd"
      )
      fit_temp <- ugarchfit(spec = spec_temp, data = gdx_returns, solver = "hybrid")
      ic <- infocriteria(fit_temp)
      std_resid  <- as.numeric(residuals(fit_temp, standardize = TRUE))
      lb_sq_p    <- Box.test(std_resid^2, lag = 10, type = "Ljung")$p.value
      lb_p    <- Box.test(std_resid, lag = 10, type = "Ljung")$p.value
      sb         <- tryCatch(signbias(fit_temp), error = function(e) NULL)
      sb_joint_p <- if (!is.null(sb)) sb["Joint Effect", "prob"] else NA
      ny <- tryCatch(nyblom(fit_temp), error = function(e) NULL)
      ny_joint <- ny_crit <- ny_beta1 <- ny_beta1_crit <- NA
      if (!is.null(ny)) {
        col_name      <- paste0(alpha_diag * 100, "%")
        ny_joint      <- tryCatch(as.numeric(ny$JointStat), error = function(e) NA)
        ny_crit       <- tryCatch(as.numeric(ny$JointCritical[col_name]), error = function(e) NA)
        ny_beta1      <- tryCatch(as.numeric(ny$IndividualStat["beta1", ]), error = function(e) NA)
        ny_beta1_crit <- tryCatch(as.numeric(ny$IndividualCritical[col_name]), error = function(e) NA)
      }
      spec_results <- rbind(spec_results, data.frame(
        Modelo = label,
        AIC = round(ic[1], 4),
        BIC = round(ic[2], 4),
        LjungBox_p = round(lb_p, 4),
        LjungBox_pasa = lb_p > alpha_diag,
        LjungBoxSq_p = round(lb_sq_p, 4),
        LjungBoxSq_pasa = lb_sq_p > alpha_diag,
        SignBiasJoint_p = round(sb_joint_p, 4),
        SignBias_pasa = sb_joint_p > alpha_diag,
        Nyblom_Joint_ratio = round(ny_joint / ny_crit, 3),
        Nyblom_Joint_pasa = (ny_joint / ny_crit) < 1,
        Nyblom_beta1 = round(ny_beta1, 4),
        Nyblom_beta1_pasa = if (is.na(ny_beta1)) NA else (ny_beta1 < ny_beta1_crit)
      ))
    }, error = function(e) cat("Error en modelo:", label, "-", conditionMessage(e), "\n"))
  }
  
  #Resultados del Bucle
  spec_results <- spec_results[order(spec_results$BIC), ]
  
  #Tabla con los Resultados
    #Criterios de Informacion - AIC y BIC, lo vimos anteriormente
    #Pasa: TRUE/FALSE, es decir, si rechazamos o no la hipotesis nula (TRUE no rechazamos)
    #LjungBoxSq_p:  p-valor de Ljung-Box sobre residuos estandarizados al cuadrado
      #Buscamos TRUE: No queda varianza sin explicar
    #SignBiasJoint_p: p-valor del Sign Bias test
      #Buscamos TRUE: No hay asimetria en el modelo sin capturar (No en la serie)
    #Nyblom_Joint_ratio: Estadistico conjunto de Nyblom  
      #Buscamos TRUE: El modelo es estable en el tiempo en su conjunto
    #Nyblom_beta1: Persistencia de la Varianza en el Tiempo
      #Buscamos TRUE: La persistencia en el corto plazo es estable en el tiempo
  print(spec_results)
  #Resultados
    #Modelos hechos para capturar asimetrias, no pasan la SignBias Test
      #Esto se debe a que logramos capturar las asimetrias con la distribucion de los modelos, sstd
  
  
  #Funcion que elige automaticamente el mejor modelo si paso todos los tests con el BIC, sino, elige el que mas tests alla pasado.
  #IMPORTANTE: La Seleccion del modelo deberia ser manual, segun los resultados de la tabla anterior, el que esta arriba no es automaticamente el ganador, por el contrario, es probable que no lo sea
    #De nuevo, no se preocupen por la complejidad del codigo, nomas por los resultados de los modelos y tests.
  select_bestmodel <- function(results) {
    logicos <- results %>% select(where(is.logical))
    rows <- nrow(logicos)
    filas_a_mantener <- logical(rows) 
    for (row in 1:rows) {
      tryCatch({
        filas_a_mantener[row] <- all(as.logical(logicos[row, ]), na.rm = TRUE)
      }, error = function(e) {
        cat("Error en fila", row, "-", conditionMessage(e), "\n")
        filas_a_mantener[row] <<- FALSE   # <<- para que el cambio salga del tryCatch
      })
    }
    resultado <- results[filas_a_mantener, ]
    if (nrow(resultado) > 0) {
      return(as.character(resultado$Modelo[1]))
    } else {
      cat("No hay modelos apropiados para los datos\n")
      results <- results %>% mutate(tests_aprobados = rowSums(logicos)) %>% arrange(desc(tests_aprobados))
      return(results$Modelo[1])
    }
  }
  #Seleccionar el Mejor Modelo segun la funcion (Lo pueden hacer manualmente editando la variable mejor modelo con el que prefieran)
  mejor_modelo <- select_bestmodel(spec_results)
  print(mejor_modelo)
  
  if(is.null(mejor_modelo)){
    mejor_modelo <- "sGARCH"
  }
  #mejor_modelo <- "csGARCH" #Por si quieren hacerlo manual
  
  #Especificamos Explicitamente el Mejor Modelo
  #"Gano" el csGARCH, es un GARCH estandar que separa la Varianza entre componentes Permanentes y Estacionarios
    #Podemos ver que es el unico que tiene los resultados que queremos al p-value de 0.05
    #Notamos dos cosas:
      #Aunque no usamos dummies, csGARCH pasa el test de estabilidad de parametros, dada la separacion de la varianza
      #El modelo ganador tiene los ICs mas altos, aunque muy cercanos entre si (Tiene mas parametros)
  best_spec <- ugarchspec(
    variance.model = list(model = mejor_modelo, garchOrder = c(1, 1)),
    mean.model = list(armaOrder = c(p, q), include.mean = FALSE),
    distribution.model = "sstd"
  )
  #Seria un csGARCH(1,1)
    #Aca podemos ver algunas de las caracteristicas del modelo
  show(best_spec)

# 4. FITTING Y TESTEO -----------------------------------------------------
  
  #Fitting del Modelo, llevamos la especificacion del modelo a los datos, encontramos el modelo concreto
  best_fit <- ugarchfit(spec = best_spec, data = gdx_returns)
  #Caracteristicas del modelo con los datos y su ajuste al mismo, incluye varias cosas que vimos antes
    #omega: Constante (Piso de la Volatilidad)
    #alpha: Magnitud de shock pasado en varianza presente
    #beta: Varianza pasada en el presente
    #csGARCH muestra parametros demas eta11 y eta21
    #skew: Medicion de la asimetria
    #shape: Forma de la Distribucion
    #Tests de Fitting
      #LM
      #Pearson GoF: Si la forma de la distribucion se ajusta bien a los residuos
  show(best_fit)
  
  #Retorno contra Volatilidad, vemos como la misma cambia y como el modelo captura la heterocedasticidad
  plot(best_fit, which = 1)
  
  #Clusters de Volatilidad, consistente con lo que sabemos de series financieras
  plot(best_fit, which = 3)
  
  #Vemos si la distribucion que elegimos fue la correcta
    #Como se puede ver, lo es
  plot(best_fit, which = 8)
  
  #Podemos ver las colas de la distribucion, el QQ plot es una herramienta comun para estas cosas
    #Mencionado en el libro
    #Si hay desvio, es prueba que nos equivocamos de distribucion
  plot(best_fit, which = 9)
  
  #Este ultimo grafico nos muestra una parabola, en resumen, la asimetria de la serie la capturamos con la distribucion, no era innata a la varianza de la misma
  #Muestra la varianza presente en base a un shock pasado
  plot(best_fit, which = 12)
  

# 5. BACKTESTING Y FORECASTING --------------------------------------------

  #Forecasting: Mismo concepto que ARMA
    #Out of  (Backtesting): Separamos parte de la muestra y hacemos un testeo comparativo
    #Muestra para comparar, retornos de 2026
  gdx_comparison <- GLD['2026']
  gdx_comparison <- (na.omit(dailyReturn(gdx_comparison, type = "log"))) * 100
    #Muestra para el modelo
  gdx_sample <- GLD['2007::2025'] #No utilizada, ilustrativa
  gdx_sample <- (na.omit(dailyReturn(gdx_sample, type = "log"))) * 100
  n_oos <- nrow(GLD['2026'])
  
    #Hacemos un nuevo modelo, con las mismas caracteristicas que el anterior, especificando el metodo out of sample
  fit_oos <- ugarchfit(spec = best_spec, data = gdx_returns, out.sample = n_oos)
    #Forecasting con las funciones del paquete rugarch
      #22 dias, o un mes de trading (n.ahead)
      #n.roll son los forecasts desde cada uno de los dias reservados
  forecast_oos <- ugarchforecast(fit_oos, n.ahead = 22, n.roll = n_oos - 1)
    
    #Mostramos el Forecast Out of Sample
      #La Serie va a 0 dado el componente ARMA(0,0) sin media
      #Nos importa mas el forecast de la volatilidad
      #No se tiende a comparar los valores, comparar volatilidad contra un retorno de un dia es un proxy ruidoso (Tsay)
        #Volatilidad Permanente y Transitoria (Por csGARCH)
  show(forecast_oos)
    #Rolling Forecast vs Actual Series: los 174 pronosticos superpuestos contra lo que realmente paso en 2026
  plot(forecast_oos, which = 2)
    #Volatilidad, cuan bien la predijo el 2026
  plot(forecast_oos, which = 4)
  #Forecast
    #Lo mismo, un mes de Trading
  forecast <- ugarchforecast(best_fit, n.ahead = 22)
  show(forecast)
    #Valor 0, dado el ARMA(0,0), nos interesa mas el bloque amarillo que la linea roja
  plot(forecast, which = 1)
    #Trayectoria de la volatilidad futura
  plot(forecast, which = 3)

# 6. APLICACIONES FINANCIERAS ---------------------------------------------

#Volatilidad
  #Series Volatility
  #Recordar que csGARCH divide entre los componentes de la volatilidad
    #Total Combinado
  volatility <- sigma(best_fit)
  #Historical Volatility (mean)
  avg_vol <- mean(volatility)
  #Current Volatility
  current_vol <- as.numeric(tail(volatility, 1)) #Ultimo valor de la serie, asumimos que es el presente
  current_vol <- current_vol %>% as.numeric(coredata(current_vol))
  #Ratio IV/HV, is it rising or falling
  vol_ratio <- current_vol/avg_vol
  cat(sprintf("Volatilidad actual: %.3f | Promedio historico: %.3f | Ratio: %.2fx\n",
              current_vol, avg_vol, vol_ratio))
  if (current_vol > avg_vol) {
    cat("La volatilidad actual esta por ENCIMA del promedio historico.\n") #Opciones mas Caras
  } else {
    cat("La volatilidad actual esta por DEBAJO del promedio historico.\n") #Opciones mas Baratas
  }
#Riesgo
  #Value at Risk
    #VaR Dinamico, dia por dia, basado en una distribucion
      #Cual es la peor perdida esperable en un dia normal, dado una probabilidad
  VaR_95 <- as.numeric(quantile(best_fit, probs = 0.05))
      # En el 95% de los dias no deberiamos perder mas de esto
  retornos_reales <- as.numeric(gdx_returns)
  var_test <- VaRTest(alpha = 0.05, actual = retornos_reales, VaR = VaR_95)
    #Cuantos dias deberian haber roto el VaR
  cat("Fallas VaR Esperadas:", var_test$expected.exceed, "\n")
    #Cuantos efectivamente lo hicieron
  cat("Fallas VaR Reales:", var_test$actual.exceed, "\n")
    #263 vs 247, modelo bien especificado
  
  #Expected Shortfall
    #Habiendo roto el VaR (umbral), cuanto podes perder?
      #Promedio de perdida en el peor de los casos
  dias_exceso <- retornos_reales[retornos_reales < VaR_95]
  expected_shortfall <- mean(dias_exceso)
  cat("Expected Shortfall (95%):", round(expected_shortfall, 2), "%\n")
  #Estos dos metodos son utilizados para Composicion de Portfolio y Tolerancia al Riesgo
  
  #Half Life - Vida Media de un Shock
  halflife(best_fit)
    #Transitorio y Permanente
      #"Permanente", no es absoluto
  #ln(0,5)/ln(alpha1+beta1) -> Formula Transitoria
  
  #Greeks
    #Greeks are derivatives of the option value with respect to underlying parameters
    #Tenemos distintos tipos de opciones (Segun los derechos y obligaciones)
      #Americana: Are like European options, but allow the holder to exercise at any time until expiration
      #Europea: Give to option holder the right but not the obligation to buy or sell, at expiry time
      #Asiatica: Have a similar payoff to European put- and call options but differ from European options in that they are path dependent
    #Extraemos ciertos valores
      # Delta: the derivative with respect to the price of the underlying asset
      # Vega: the derivative with respect to the volatility
      # Theta: the negative derivative with respect to the time until expiration of the option
      # rho: the derivative with respect to the risk-free interest rate
      # Epsilon: the derivative with respect to the dividend yield of the underlying asset
      # Lambda: the percentage change in an option contract's price to the percentage change in the price of the underlying security.
      # Gamma: the second derivative with respect to the price of the underlying asset
      # Vanna: the derivative of Delta with respect to the volatility
      # Vomma: the second derivative with respect to the volatility
      # Veta: the derivative of Vega with respect to the time until expiration of the option
      # Vera: the derivative of Vega with respect to the risk-free interest rate
      # Speed: the third derivative of the option value with respect to the price of the underlying asset
      # Zomma: the derivative of Gamma with respect to the volatility
      # Color: the derivative of Gamma with respect to the time
      # Ultima: the third derivative with respect to the volatility
        #Delta, Gamma, Vega, Theta, Rhoson las mas usadas en la practica.
    #Reconocemos tambien por que movimiento consideramos en el instrumento
      #Put
      #Call
    #Ademas de su periodo (Long y Short)
    #Su ubicacion en relacion al precio (At the Money/Out of the money)
      #Greeks_UI()
  vol_anualizada_decimal <- (current_vol / 100) * sqrt(252)
  spot_price <- as.numeric(tail(Cl(GLD), 1))
    #Esto Seria una opcion call at the money (Lo mas comun para mostrar griegas)
      #At the money: No perdes si se ejecuta ahora
    #Basado en Black-Scholes como motor matematico para extraer estos valores
  opcion_call <- as.list(Greeks(
    initial_price = spot_price, #Precio Actual
    exercise_price = round(spot_price), #Precio del contrato; Usamos el precio de hoy redondeado, como referencia; En general, el precio al que predecimos que se va a mover el instrumento
    r = 0.05, #Tasa libre de riesgo (Proxy)
    time_to_maturity = 30 / 365, #Tiempo de duracion del contrato, en años (Por B-S)
    volatility = vol_anualizada_decimal, #Volatilidad que calculamos en el modelo GARCH, es lo que mueve el precio de la accion
    payoff = "call" #Derecho a comprar la accion al precio de ejercicio; Put derecho a vender
  ))
    #Las griegas van a ser los comportamientos matematicos del instrumento en base a los valores que le damos ahora a esta funcion
  cat("\n--- PRECIO JUSTO Y GRIEGAS (CALL 30 DÍAS) ---\n")
  print(opcion_call)
  
  greeks_claros <- function(g) {
    data.frame(
      Griega = c("Precio", "Delta", "Gamma", "Vega", "Theta", "Rho"),
      Valor = c(round(g$fair_value, 2), round(g$delta, 4), round(g$gamma, 4),
                round(g$vega/100, 4), round(g$theta/365, 4), round(g$rho/100, 4)),
      Interpretacion = c(
        sprintf("La opcion vale $%.2f hoy", g$fair_value),
        sprintf("Si GLD sube $1, la opcion sube ~$%.2f", g$delta),
        sprintf("El delta se mueve ~%.4f por cada $1 que cambia GLD", g$gamma),
        sprintf("Si la vol anualizada sube 1 punto (ej. 26%%->27%%), la opcion sube ~$%.2f", g$vega/100),
        sprintf("Cada dia que pasa sin que nada mas cambie, la opcion pierde ~$%.2f", abs(g$theta/365)),
        sprintf("Si la tasa sube 1 punto, la opcion sube ~$%.2f", g$rho/100)
      )
    )
  }
  print(greeks_claros(opcion_call), right = FALSE)
  
  #Graficos Griegas
    #Recalcula la opcion 100 veces, una por cada precio hipotetico de GLD -- misma
    #volatilidad y mismo strike, solo varia el spot
    #Codigo complejo, importa mas el grafico que otra cosa
  K <- round(spot_price)
  spot_grid <- seq(spot_price * 0.7, spot_price * 1.3, length.out = 100)
  
  greeks_grid <- lapply(spot_grid, function(S) {
    as.list(Greeks(initial_price = S, exercise_price = K, r = 0.05,
                   time_to_maturity = 30/365, volatility = vol_anualizada_decimal,
                   payoff = "call"))
  })
  
  df_grid <- data.frame(
    spot = spot_grid,
    precio_hoy = sapply(greeks_grid, function(g) g$fair_value),
    delta = sapply(greeks_grid, function(g) g$delta),
    gamma = sapply(greeks_grid, function(g) g$gamma),
    vega  = sapply(greeks_grid, function(g) g$vega),
    theta = sapply(greeks_grid, function(g) g$theta)
  )
  df_grid$payoff_vencimiento <- pmax(df_grid$spot - K, 0)
  
  tangente <- opcion_call[["fair_value"]] + opcion_call[["delta"]] * (spot_grid - spot_price)
  
  p3 <- ggplot(df_grid, aes(x = spot)) +
    geom_line(aes(y = payoff_vencimiento, color = "Payoff al vencimiento"), linewidth = 1) +
    geom_line(aes(y = precio_hoy, color = "Valor de la opcion hoy"), linewidth = 1) +
    geom_line(aes(y = tangente), color = "gray40", linetype = "dashed") +
    geom_point(data = data.frame(x = spot_price, y = opcion_call[["fair_value"]]),
               aes(x = x, y = y), size = 3, color = "black") +
    geom_vline(xintercept = K, linetype = "dotted") +
    scale_color_manual(values = c("Payoff al vencimiento" = "gray50",
                                  "Valor de la opcion hoy" = "steelblue")) +
    labs(title = "Call sobre GLD: payoff vs. valor actual",
         subtitle = sprintf("Strike = $%d | Vol anualizada = %.1f%% | 30 dias al vencimiento",
                            K, vol_anualizada_decimal*100),
         x = "Precio de GLD", y = "Valor (USD)", color = NULL) +
    theme_minimal() + theme(legend.position = "bottom")
  
  df_long <- data.frame(
    spot = rep(df_grid$spot, 4),
    griega = rep(c("Delta", "Gamma", "Vega", "Theta"), each = nrow(df_grid)),
    valor = c(df_grid$delta, df_grid$gamma, df_grid$vega/100, df_grid$theta/365)
  )
  
  p4 <- ggplot(df_long, aes(spot, valor)) +
    geom_line(color = "firebrick") +
    geom_vline(xintercept = spot_price, linetype = "dotted") +
    facet_wrap(~griega, scales = "free_y") +
    labs(title = "Como cambian las Griegas segun el precio de GLD",
         x = "Precio de GLD", y = NULL) +
    theme_minimal()
  
  p3 / p4 #Cuanto vale la opcion el dia de vencimiento, comportamiento de las griegas en relacion a esa fecha
  
  #Risk Premium - Teorico
  risk_premium_teorico <- mean(retornos_reales) / mean(volatility^2)
  cat("\n--- RISK PREMIUM ---\n")
  cat("Retorno por unidad de varianza:", round(risk_premium_teorico, 4), "\n")
  
  #Riego Vizualizado
  df_plot <- data.frame(
    date = index(gdx_returns),
    cum_return = cumsum(as.numeric(gdx_returns)),
    ann_vol = as.numeric(sigma(best_fit)) * sqrt(252)
  )
  
  p1 <- ggplot(df_plot, aes(date, cum_return)) +
    geom_line(color = "steelblue") +
    labs(title = "Retorno acumulado de GLD", x = NULL, y = "Retorno acum. (%)") +
    theme_minimal()
  
  p2 <- ggplot(df_plot, aes(date, ann_vol)) +
    geom_line(color = "firebrick") +
    labs(title = "Volatilidad condicional anualizada", x = NULL, y = "Vol (%)") +
    theme_minimal()
  p1 / p2 #Volatilidad y Riesgo Coinciden Visualmente
  