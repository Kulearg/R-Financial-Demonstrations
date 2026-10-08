library(tidyverse) #Manejo de Datos
library(tstests) #Tests Varios de Series de Tiempo
library(ggplot2) #Graficos
library(quantmod) #Datos financieros
library(tseries) #Paquete Principal de series de tiempo
library(strucchange) #Paquete para quiebres estructurales
library(forecast) #Para modelos ARMA, principalmente para auto.arima()
library(vars)

# 0. CARGA DE DATOS -------------------------------------------------------

  #Traemos los datos usando el Ticker, con el paquete quantmod
    #Pide al usuario que seleccione un ticker (accion/activo) de finance.yahoo.com
    #Este codigo fue hecho con SPY en Mente
  #ticker <- as.character(readline(prompt = "Enter Ticker: "))
  #getSymbols(ticker, src = "yahoo", auto.assign = FALSE) -> ticker_data_raw
  getSymbols("SPY", src = "yahoo", auto.assign = FALSE) -> ticker_data_raw
  #Data Macro FRED
  getSymbols("DGS1",  src = "FRED", auto.assign = FALSE) -> rate_raw

# PREPARACION DE DATOS -----------------------------------------------
  #Hacemos que los dos datos esten con las mismas fechas
  #Sacamos los Retornos
  ticker_returns = dailyReturn(ticker_data_raw, type = "log")
  #Omitimos 0s y NaNs
  ticker_returns = na.omit(ticker_returns)
  rate = na.omit(rate_raw)
  rate = rate["2016/2019"]
  ticker_returns = ticker_returns["2016/2019"]
  #VAR Requiere estacionariedad
  # Stationarity tests on a series and its first difference
  # Requires: install.packages("tseries")
  stationarity_tests <- function(ts, alpha = 0.05) {
    ts      <- as.numeric(na.omit(ts))   # tests fail on NAs; also strips xts/zoo dates
    diff_ts <- diff(ts)
    results <- data.frame()
    for (level in c("Level", "First difference")) {
      series <- if (level == "Level") ts else diff_ts
      # tseries warns when p-values fall outside its lookup tables; silence that
      adf  <- suppressWarnings(adf.test(series))
      pp   <- suppressWarnings(pp.test(series))
      kpss <- suppressWarnings(kpss.test(series))
      p_values <- c(adf$p.value, pp$p.value, kpss$p.value)
      # ADF and PP: H0 = unit root  -> reject (p < alpha) means stationary
      # KPSS:       H0 = stationary -> reject (p < alpha) means non-stationary
      stationary <- c(p_values[1] < alpha,
                      p_values[2] < alpha,
                      p_values[3] >= alpha)
      results <- rbind(results, data.frame(
        Series     = level,
        Test       = c("ADF", "PP", "KPSS"),
        Null       = c("Unit root", "Unit root", "Stationary"),
        Statistic  = round(unname(c(adf$statistic, pp$statistic, kpss$statistic)), 4),
        P_value    = round(p_values, 4),
        Conclusion = ifelse(stationary, "Stationary", "Non-stationary")
      ))
    }
    return(results)
  }
  rate_tests = stationarity_tests(rate)
  ticker_tests = stationarity_tests(ticker_returns)
  print(rate_tests)
  print(ticker_tests)
  #Unimos ambas cosas en un objeto para luego ver el modelo VAR
    #Cholesky, ordenar mas inmediato a menos
    #Vamos a usar la Tasa diferenciada
  model_data = merge(diff(rate), ticker_returns)
  model_data = na.omit(model_data)

###########################################################################
  #Diferenciar series no estacionarias saca informacion de la tendencia de largo plazo
  #Es un problema si son series cointegradas, en ese caso utilizamos un VECM (Vector ECM)
    #Esto preservaria la tendencia a largo plazo
  #Vamos a Asumir que no hay Cointegracion - Por simplicidad del ejemplo
###########################################################################
# 2. MODELO VAR -------------------------------------------------------
  #VAR nos permite capturar como una variable afecta a la otra en una regresion
  #Todas las Vairables Dinamicamente Conectadas
  #Recordemos las dependencias lineales de cada variable
  #Seleccion de Lags
  model_selection = VARselect(model_data, lag.max = 12, type = 'const')
  print(model_selection)
  lag_bic = model_selection$selection["SC(n)"]
  #Estimacion del Modelo, Basado en los Lags Anteriores
  var_model = VAR(model_data, p = lag_bic, type = 'const')
  summary(var_model)
  #Raices
  roots(var_model)
  #Diagnostico Errores, Volatilidad, Normalidad
  serial.test(var_model)
  arch.test(var_model)
  normality.test(var_model)
  #Nombres de las Variables
  data_names = names(var_model$varresult)
  #Causalidad de Granger
  causality(var_model)
  for (v in data_names) {
    cat("Causalidad de:", v)
    print(causality(var_model, cause = v))
    print(causality(var_model, cause = v, vcov. = vcovHC(var_model))$Granger)
  }
    #Mide si los valores historicos de una variable dan informacion predicitiva de la otra
  #Sacar Parametros no Significativos - Mas parsimonia
  try(restrict(var_model))
  #Forma Estructural
    #Verificamos los resultados de Cholesky
    #Debemos darle Amat - B0
  #El modelo A impone restricciones sobre B0. El modelo B asumeque las restricciones aplican directamente a los errores. 
  #Para indicar que algunos par ́ametros no se van a estimar se usa NA.
  amat <- diag(2)
  diag(amat) <- NA      # NA = estimate this value
  amat[2, 1] <- NA      # returns can react to the rate within the same day
  svar_model <- suppressWarnings(SVAR(var_model, estmethod = "direct", Amat = amat))
  #supressWarnings, Si el modelo esta justamente identificado
    #In SVAR(var_model, estmethod = "direct", Amat = amat) :  The A-model is just identified. No test possible.
  summary(svar_model)
  solve(svar_model$A)
  
###########################################################################

# 3. BACKTESTING Y FORECASTING ----------------------------------------
  #Vamos a hacer Forecasting con 22 Dias (Un mes de trading)
  forecasting_foward = 22
  #Impulso Respuesta
    #Probar todas con un Bucle
  for (v in data_names) {
    for (f in data_names) {
      if(v != f){
        cat("Funcion Impulso Respuesta de", v, "sobre", f)
        temp_irf = irf(
          var_model,
          impulse = v,
          response = f,
          n.ahead = forecasting_foward
        )
        print(temp_irf)
        plot(temp_irf)
        break
      }
    }
  }
  backtest = fitted(var_model)
  plot.ts(backtest, main = "Fitted values")
  #Forecast
  forecast <- predict(var_model, n.ahead = forecasting_foward, ci = 0.95)
  plot(forecast)
  #Forecast Error Variance Decomposition
  model_forecasterror = fevd(var_model, n.ahead = forecasting_foward)
  plot(model_forecasterror)
###########################################################################

# 4. APLICACIONES FINANCIERAS -----------------------------------------

  
  #Sensitivity to Rates by Sector
    #Los Tickes X son ETFs Sectoriales
    #Cuales tienen mayor sensibilidad a un cambio de Tasa
  tickers <- c("SPY", "XLF", "XLU", "XLK", "XLRE")   # market, banks, utilities, tech, real estate
  rate_share <- sapply(tickers, function(t) {
    r <- diff(log(Ad(getSymbols(t, from = "2016-01-01", to = "2019-12-31", auto.assign = FALSE))))
    d <- na.omit(merge(diff(rate), r))
    fevd(VAR(coredata(d), p = 1, type = "const"), n.ahead = 10)[[2]][10, 1] * 100
  })
  rate_share <- sort(rate_share)
  names(rate_share) <- sub("\\..*", "", names(rate_share))   # "XLF.DGS1" -> "XLF"
  cat("\n--- 1. Rate exposure by sector ---\n")
  for (t in names(rate_share)) {
    cat(t, ":", round(rate_share[t], 2), "% of its return variance comes from rate shocks\n")
  }
  cat("Most rate-sensitive:", names(which.max(rate_share)),
      "| Least rate-sensitive:", names(which.min(rate_share)), "\n")
  
  #Rate Exposure over Time
    #Estabilidad de la Sensibilidad
  win  <- 250
  ends <- seq(win, nrow(model_data), by = 10)
  roll <- sapply(ends, function(e) {
    m <- VAR(coredata(model_data[(e - win + 1):e, ]), p = 1, type = "const")
    fevd(m, n.ahead = 10)$daily.returns[10, "DGS1"] * 100
  })
  plot(xts(roll, index(model_data)[ends]), main = "% of return variance from rate shocks")
  cat("\n--- 2. Rate exposure over time ---\n")
  cat("Rate shocks explained between", round(min(roll), 1), "% and", round(max(roll), 1),
      "% of return variance (average", round(mean(roll), 1), "%).\n")
  cat("The relationship between rates and stocks changes over time.\n")
  
  #Rate Shocks (IRF)
    #Shock de 25 bps
  ir <- irf(var_model, impulse = "DGS1", response = "daily.returns",
            n.ahead = 5, ortho = TRUE, boot = FALSE)
  s  <- sd(residuals(var_model)[, "DGS1"])       
  stress <- ir$irf$DGS1 * (0.25 / s)           
  cbind(daily = stress, cumulative = cumsum(stress)) 
  
  cat("\n--- 3. Stress test ---\n")
  cat("A 25 bp surprise rise in the 1-year yield implies a return of",
      round(100 * stress[1], 2), "% on the same day and", round(100 * sum(stress), 2), "% in total.\n")
  cat("25 bp is about", round(0.25 / s, 1),
      "standard deviations of a typical daily rate shock, so it is an extreme scenario.\n")
  
  #VAR Value at Risk
    #Cuanto perdes en un mal dia
    #Version con Distribucion Normal y Empirica
  fc <- predict(var_model, n.ahead = 1, ci = 0.98)
  var_normal <- fc$fcst$daily.returns[1, "lower"]                          
  var_emp    <- fc$fcst$daily.returns[1, "fcst"] +
    quantile(residuals(var_model)[, "daily.returns"], 0.01)           
  cat("\n--- 4. Value at Risk (99%, one day) ---\n")
  cat("Normal VaR:", round(100 * var_normal, 2),
      "% -> losses should exceed this on only 1 in 100 days.\n")
  cat("Empirical VaR:", round(100 * var_emp, 2),
      "% -> larger, because real returns have fat tails.\n")
  
  #Trading Signal
    #Prediccion Out of Sample
    #La tasa de ayer puede predecir el retorno de hoy
  n <- nrow(model_data); train <- 500
  signal <- sapply(train:(n - 1), function(t) {
    m <- VAR(coredata(model_data[1:t, ]), p = 1, type = "const")
    predict(m, n.ahead = 1)$fcst$daily.returns[1, "fcst"]
  })
  real  <- as.numeric(model_data$daily.returns)[(train + 1):n]
  strat <- ifelse(signal > 0, 1, 0) * real
  hit   <- mean(sign(signal) == sign(real))
  
  cat("\n--- 5. Trading signal ---\n")
  cat("Strategy return:", round(100 * sum(strat), 1), "% | Buy and hold:",
      round(100 * sum(real), 1), "%\n")
  cat("Hit rate:", round(100 * hit, 1),
      "% of days the model predicted the direction correctly (50% = coin flip).\n")
  cat(ifelse(sum(strat) > sum(real),
             "The strategy beat buy-and-hold, before transaction costs.",
             "The strategy did not beat buy-and-hold."), "\n")
###########################################################################


# 5. GRAFICOS -------------------------------------------------------------

  labels <- c(DGS1 = "1-year yield change", daily.returns = "Stock return")
  units  <- c(DGS1 = "basis points", daily.returns = "%")
  cols   <- c("steelblue", "grey80")
  
  # Rolling rate exposure
  roll_dates <- index(model_data)[ends]
  plot(roll_dates, roll, type = "l", lwd = 2, col = "steelblue",
       xlab = "", ylab = "% of return variance",
       main = "Share of stock return variance explained by rate shocks",
       sub  = "Rolling 250-day window, FEVD at 10 days")
  abline(h = mean(roll), lty = 2, col = "grey40")
  grid()
  legend("topleft", legend = c("Rolling share", "Average"),
         col = c("steelblue", "grey40"), lty = c(1, 2), lwd = c(2, 1), bty = "n")
  
  #FEVD
  fe <- fevd(var_model, n.ahead = 10)
  par(mfrow = c(1, 2), oma = c(3, 0, 0, 0))
  for (v in data_names) {
    barplot(t(fe[[v]]) * 100, names.arg = 1:10, col = cols, border = NA,
            ylim = c(0, 100), main = paste("FEVD:", labels[v]),
            xlab = "Horizon (days)", ylab = "% of forecast error variance")
  }
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot.new()
  legend("bottom", legend = paste("Shock to", labels[data_names]),
         fill = cols, border = NA, horiz = TRUE, bty = "n")
  par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)
  
  # Forecast (last 60 days + forecast with band)
  last_n <- 60
  par(mfrow = c(2, 1), mar = c(4, 4, 3, 1))
  for (v in data_names) {
    hist_v <- tail(as.numeric(model_data[, v]), last_n) * 100
    fc_v   <- forecast$fcst[[v]][, c("fcst", "lower", "upper")] * 100
    x_hist <- 1:last_n
    x_fc   <- last_n + 1:nrow(fc_v)
    
    plot(x_hist, hist_v, type = "l", col = "grey40",
         xlim = c(1, max(x_fc)), ylim = range(hist_v, fc_v),
         xlab = "Days", ylab = units[v],
         main = paste("Forecast:", labels[v], "(last", last_n, "days + forecast)"))
    polygon(c(x_fc, rev(x_fc)), c(fc_v[, "lower"], rev(fc_v[, "upper"])),
            col = adjustcolor("steelblue", 0.25), border = NA)
    lines(c(last_n, x_fc), c(hist_v[last_n], fc_v[, "fcst"]), col = "steelblue", lwd = 2)
    abline(v = last_n + 0.5, lty = 3)
  }
  par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)
  
  # Actual vs fitted (with dates)
  p      <- var_model$p
  dates  <- index(model_data)[-(1:p)]
  actual <- var_model$y[-(1:p), ]
  par(mfrow = c(2, 1), mar = c(3, 4, 3, 1))
  for (v in data_names) {
    plot(dates, actual[, v] * 100, type = "l", col = "grey75",
         xlab = "", ylab = units[v], main = paste("Actual vs fitted:", labels[v]))
    lines(dates, backtest[, v] * 100, col = "firebrick", lwd = 1.5)
    legend("topleft", legend = c("Actual", "Fitted"), col = c("grey75", "firebrick"),
           lwd = c(1, 1.5), bty = "n", horiz = TRUE)
  }
  par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)
  
  # Impulse responses (10 days, readable units, shaded band)
  set.seed(123)
  for (v in data_names) {
    for (f in data_names) {
      if (v != f) {
        ir  <- irf(var_model, impulse = v, response = f, n.ahead = 10,
                   boot = TRUE, runs = 500, ci = 0.95)
        est <- ir$irf[[v]][, f]   * 100
        lo  <- ir$Lower[[v]][, f] * 100
        up  <- ir$Upper[[v]][, f] * 100
        h   <- 0:10
        
        plot(h, est, type = "n", ylim = range(lo, up, 0),
             xlab = "Days after the shock", ylab = paste0("Response (", units[f], ")"),
             main = paste("Response of", labels[f], "to a shock in", labels[v]),
             sub  = "One-standard-deviation shock, 95% bootstrap interval (500 runs)")
        polygon(c(h, rev(h)), c(lo, rev(up)), col = adjustcolor("steelblue", 0.25), border = NA)
        abline(h = 0, col = "grey40")
        lines(h, est, type = "b", pch = 19, lwd = 2, col = "steelblue")
      }
    }
  }

  