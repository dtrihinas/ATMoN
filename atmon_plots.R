library(plotly)

#error
metrics <- c('mit.maxclique', 'sg09.diameter', 'sg09.components', 'sg09.degdistro', 'vanet.effdiameter', 'vanet.gcRatio')
adam <-          c(18, 20, 16, 26, 09, 15)
adamtopo <-      c(19, 23, 24, 27, 12, 19)
atmon <-         c(07, 08, 05, 12, 05, 07)
avgrandomwalk <- c(49, 38, 43, 46, 28, 36)

pdata <- data.frame(metrics, adam, adamtopo, atmon, avgrandomwalk)

#The default order will be alphabetized unless specified as below:
pdata$x <- factor(pdata$metrics, levels = pdata[["metrics"]])

plot_ly(pdata, x = ~metrics, y = ~adam, type = 'bar', name = 'AdaM', marker = list(color = 'rgb(0,0,0)')) %>%
  add_trace(y = ~adamtopo, name = 'AdaM-topo', marker = list(color = 'rgb(80,80,80)')) %>%
  add_trace(y = ~atmon, name = 'ATMoN', marker = list(color = 'rgb(144,144,144)')) %>%
  add_trace(y = ~avgrandomwalk, name = 'AvgRandomWalk', marker = list(color = 'rgb(216,216,216)')) %>%
  
  layout(
         xaxis = list(title = "metrics", tickangle = -45),
         yaxis = list(title = "MAPE (%)"),
         margin = list(b = 100),
         barmode = 'group',
         legend = list(orientation = "h", xanchor = "center", x = 0.5, yanchor = "top", y=1.1)           
        )

#cc
edatasets <- c('mit', 'sg09', 'vanet')

ccbaseline <- c(7300, 12770, 31023)
ccadam <- c(1095, 2128, 7065)
ccatmon <- c(1314, 2689, 7986)
ccavgrandomwalk <- c(3900, 5618, 13411)

ccdata <- data.frame(edatasets, ccbaseline, ccadam, ccatmon, ccavgrandomwalk)
ccdata$x <- factor(ccdata$edatasets, levels = ccdata[["edatasets"]])

plot_ly(ccdata, x = ~edatasets, y = ~ccbaseline, type = 'bar', name = 'Baseline', marker = list(color = 'rgb(0,0,0)')) %>%
  add_trace(y = ~ccadam, name = 'AdaM', marker = list(color = 'rgb(80,80,80)')) %>%
  add_trace(y = ~ccatmon, name = 'ATMoN', marker = list(color = 'rgb(144,144,144)')) %>%
  add_trace(y = ~ccavgrandomwalk, name = 'AvgRandomWalk', marker = list(color = 'rgb(216,216,216)')) %>%
  
  layout(
    xaxis = list(title = "datasets", tickangle = -45),
    yaxis = list(title = "CPU Cycles (x10^9)", range=c(0,36000)),
    margin = list(b = 100),
    barmode = 'group',
    legend = list(orientation = "h", xanchor = "center", x = 0.5, yanchor = "top", y=1.1)           
  )

#cost

prbaseline <- c(5.85, 19.02, 94.15)
pradam <- c(2.11, 4.45, 26.30)
pratmon <- c(2.49, 5.12, 28.80)
pravgrandomwalk <- c(4.10, 9.13, 50.41)

prdata <- data.frame(edatasets, prbaseline, pradam, pratmon, pravgrandomwalk)
prdata$x <- factor(prdata$edatasets, levels = prdata[["edatasets"]])

plot_ly(prdata, x = ~edatasets, y = ~prbaseline, type = 'bar', name = 'Baseline', marker = list(color = 'rgb(0,0,0)')) %>%
  add_trace(y = ~pradam, name = 'AdaM', marker = list(color = 'rgb(80,80,80)')) %>%
  add_trace(y = ~pratmon, name = 'ATMoN', marker = list(color = 'rgb(144,144,144)')) %>%
  add_trace(y = ~pravgrandomwalk, name = 'AvgRandomWalk', marker = list(color = 'rgb(216,216,216)')) %>%
  
  layout(
    xaxis = list(title = "datasets", tickangle = -45),
    yaxis = list(title = "Cost ($)", range = c(0, 101)),
    margin = list(b = 100),
    barmode = 'group',
    legend = list(orientation = "h", xanchor = "center", x = 0.5, yanchor = "top", y=1.1)
  )

#data volume

dvbaseline <- c(16.49, 45.24, 108)
dvadam <- c(6.10, 11.05, 28.24)
dvatmon <- c(7.03, 12.23, 32.32)
dvavgrandomwalk <- c(12.03, 31.3, 60.7)

dvdata <- data.frame(edatasets, dvbaseline, dvadam, dvatmon, dvavgrandomwalk)
dvdata$x <- factor(dvdata$edatasets, levels = dvdata[["edatasets"]])

plot_ly(dvdata, x = ~edatasets, y = ~dvbaseline, type = 'bar', name = 'Baseline', marker = list(color = 'rgb(0,0,0)')) %>%
  add_trace(y = ~dvadam, name = 'AdaM', marker = list(color = 'rgb(80,80,80)')) %>%
  add_trace(y = ~dvatmon, name = 'ATMoN', marker = list(color = 'rgb(144,144,144)')) %>%
  add_trace(y = ~dvavgrandomwalk, name = 'AvgRandomWalk', marker = list(color = 'rgb(216,216,216)')) %>%
  
  layout(
    xaxis = list(title = "datasets", tickangle = -45),
    yaxis = list(title = "GB", range = c(0,121)),
    margin = list(b = 100),
    barmode = 'group',
    legend = list(orientation = "h", xanchor = "center", x = 0.5, yanchor = "top", y=1.1)
  )