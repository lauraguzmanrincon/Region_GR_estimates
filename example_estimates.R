
# Code for estimating growth rate of a region based on finer-scale data
# Example of how to run model (example: Midlands/Scotland/England)
# Output: tableEstimates

library(data.table)
library(INLA)
library(ggplot2)

source("functions.R")

## ---- Settings ----
# Run one of these blocks:

# -------------------------------- #
# Example Scotland
load("Datasets/Scotland_councilAreas.RData", verbose = T)
load("Datasets/Scotland_exampleData.RData", verbose = T)

typeModel <- "GaussianProcessRaw" # Choose "GaussianProcessRaw" or "GaussianProcessProp"
lastDateAvailable <- as.Date("2022-09-11")
regionToModel <- "Scotland"
ltlaInRegion <- areasScotland[order(Code), Code]
populationInRegion <- areasScotland[order(Code), population_estim]
# -------------------------------- #

# -------------------------------- #
# Example Midlands
#load("Datasets/England_LTLA_NHSER.RData", verbose = T)
#load("Datasets/England_exampleData.RData", verbose = T)

#typeModel <- "GaussianProcessRaw" # Choose "GaussianProcessRaw" or "GaussianProcessProp"
#lastDateAvailable <- as.Date("2022-09-22")
#regionToModel <- "Midlands"
#ltlaInRegion <- LTLAtoNHSER[NHSER_name == "Midlands"][order(LTLA_code), LTLA_code]
#populationInRegion <- ltlaToNHSER[NHSER_name == "Midlands"][order(LTLA_code), population_estim]
# -------------------------------- #

# -------------------------------- #
# Example England
#load("Datasets/England_LTLA_NHSER.RData", verbose = T)
#load("Datasets/England_exampleData.RData", verbose = T)

#typeModel <- "GaussianProcessRaw" # Choose "GaussianProcessRaw" or "GaussianProcessProp"
#lastDateAvailable <- as.Date("2022-09-22")
#regionToModel <- "England"
#ltlaInRegion <- LTLAtoNHSER[order(LTLA_code), LTLA_code]
#populationInRegion <- ltlaToNHSER[order(LTLA_code), population_estim]
# -------------------------------- #

## ---- Run model for each local authority ----
daysToModel <- 121
minDateModel <- lastDateAvailable - daysToModel + 1
maxDateModel <- lastDateAvailable

parametersModel <- setParametersFn(linkType = ifelse(typeModel == "GaussianProcessProp", "BB", "NB"),
                                   prior2.sigma0 = ifelse(typeModel == "GaussianProcessProp", 2.26, 5.97),
                                   prior2.lengthscale0 = ifelse(typeModel == "GaussianProcessProp", 59.34, 120.97),
                                   theta.prior2.mean = c(0,0),
                                   theta.prior2.prec = c(1,100),
                                   BB.prior.rho = list(overdispersion = list(prior = "gaussian", param = c(0, 0.4))),
                                   NB.prior.rho = list(size = list(prior = "pc.mgamma", param = c(7))),
                                   dayWeek = T,
                                   dayWeek.prior.prec = list(theta = list(prior = 'loggamma', param = c(1, 0.01))),
                                   sizeSample = 1000)

cubeSampleGP <- array(0, dim = c(daysToModel, length(ltlaInRegion), parametersModel$config$sizeSample))
cubeSampleDerivatives <- array(0, dim = c(daysToModel, length(ltlaInRegion), parametersModel$config$sizeSample))
removeLTLA <- rep(1, length(ltlaInRegion))
for(i in 1:length(ltlaInRegion)){
  tryCatch({
    cat("\nModel for ", ltlaInRegion[i], ":\n", sep = "")
    countTable <- counTable_region[LTLA_code == ltlaInRegion[i] & date >= minDateModel & date <= maxDateModel, .(date, positiveResults, numberTest)]
    (outputSimulation <- runModelGrowthRate(countTable = countTable, parametersModel = parametersModel, saveSamples = T, minDate = minDateModel, maxDate = maxDateModel))
    cubeSampleGP[,i,] <- outputSimulation$matrixSampleDays
    cubeSampleDerivatives[,i,] <- outputSimulation$sampleDerivatives
  }, error = function(e) {
    removeLTLA[i] <<- 0
    })
}

## ---- Get samples of rR and RR ----
if(parametersModel$params$linkType == "BB"){
  weigthMatrix <- aperm(array(populationInRegion/10000, dim = c(length(populationInRegion), daysToModel, parametersModel$config$sizeSample)), c(2,1,3))
}else if(parametersModel$params$linkType == "NB"){
  weigthMatrix <- exp(cubeSampleGP)
}else{
  error("Wrong choice of link type from setParametersFn()")
}

INC_P <- 0.4162
GAMMA <- 0.7309
tttNum <- (weigthMatrix[,removeLTLA,]*cubeSampleDerivatives[,removeLTLA,])
tttDen <- weigthMatrix[,removeLTLA,]

samples_rR <- apply(tttNum, c(1,3), sum)/apply(tttDen, c(1,3), sum)
samples_RR <- ( (1+samples_rR/(3*INC_P))^3 )*(1+samples_rR/GAMMA)
tableSamples_rR <- data.table(date = seq.Date(from = minDateModel, to = maxDateModel, by = "day"),
                              mean = apply(samples_rR, 1, mean, na.rm = T),
                              median = apply(samples_rR, 1, quantile, 0.5, na.rm = T),
                              q0.025 = apply(samples_rR, 1, quantile, 0.025, na.rm = T),
                              q0.975 = apply(samples_rR, 1, quantile, 0.975, na.rm = T),
                              q0.05 = apply(samples_rR, 1, quantile, 0.05, na.rm = T),
                              q0.10 = apply(samples_rR, 1, quantile, 0.10, na.rm = T),
                              q0.15 = apply(samples_rR, 1, quantile, 0.15, na.rm = T),
                              q0.20 = apply(samples_rR, 1, quantile, 0.20, na.rm = T),
                              q0.25 = apply(samples_rR, 1, quantile, 0.25, na.rm = T),
                              q0.30 = apply(samples_rR, 1, quantile, 0.30, na.rm = T),
                              q0.35 = apply(samples_rR, 1, quantile, 0.35, na.rm = T),
                              q0.40 = apply(samples_rR, 1, quantile, 0.40, na.rm = T),
                              q0.45 = apply(samples_rR, 1, quantile, 0.45, na.rm = T),
                              q0.50 = apply(samples_rR, 1, quantile, 0.50, na.rm = T),
                              q0.55 = apply(samples_rR, 1, quantile, 0.55, na.rm = T),
                              q0.60 = apply(samples_rR, 1, quantile, 0.60, na.rm = T),
                              q0.65 = apply(samples_rR, 1, quantile, 0.65, na.rm = T),
                              q0.70 = apply(samples_rR, 1, quantile, 0.70, na.rm = T),
                              q0.75 = apply(samples_rR, 1, quantile, 0.75, na.rm = T),
                              q0.80 = apply(samples_rR, 1, quantile, 0.80, na.rm = T),
                              q0.85 = apply(samples_rR, 1, quantile, 0.85, na.rm = T),
                              q0.90 = apply(samples_rR, 1, quantile, 0.90, na.rm = T),
                              q0.95 = apply(samples_rR, 1, quantile, 0.95, na.rm = T),
                              #
                              meanRR = apply(samples_RR, 1, mean, na.rm = T),
                              medianRR = apply(samples_RR, 1, quantile, 0.5, na.rm = T),
                              q0.025RR = apply(samples_RR, 1, quantile, 0.025, na.rm = T),
                              q0.975RR = apply(samples_RR, 1, quantile, 0.975, na.rm = T),
                              q0.05RR = apply(samples_RR, 1, quantile, 0.05, na.rm = T),
                              q0.10RR = apply(samples_RR, 1, quantile, 0.10, na.rm = T),
                              q0.15RR = apply(samples_RR, 1, quantile, 0.15, na.rm = T),
                              q0.20RR = apply(samples_RR, 1, quantile, 0.20, na.rm = T),
                              q0.25RR = apply(samples_RR, 1, quantile, 0.25, na.rm = T),
                              q0.30RR = apply(samples_RR, 1, quantile, 0.30, na.rm = T),
                              q0.35RR = apply(samples_RR, 1, quantile, 0.35, na.rm = T),
                              q0.40RR = apply(samples_RR, 1, quantile, 0.40, na.rm = T),
                              q0.45RR = apply(samples_RR, 1, quantile, 0.45, na.rm = T),
                              q0.50RR = apply(samples_RR, 1, quantile, 0.50, na.rm = T),
                              q0.55RR = apply(samples_RR, 1, quantile, 0.55, na.rm = T),
                              q0.60RR = apply(samples_RR, 1, quantile, 0.60, na.rm = T),
                              q0.65RR = apply(samples_RR, 1, quantile, 0.65, na.rm = T),
                              q0.70RR = apply(samples_RR, 1, quantile, 0.70, na.rm = T),
                              q0.75RR = apply(samples_RR, 1, quantile, 0.75, na.rm = T),
                              q0.80RR = apply(samples_RR, 1, quantile, 0.80, na.rm = T),
                              q0.85RR = apply(samples_RR, 1, quantile, 0.85, na.rm = T),
                              q0.90RR = apply(samples_RR, 1, quantile, 0.90, na.rm = T),
                              q0.95RR = apply(samples_RR, 1, quantile, 0.95, na.rm = T))

## ---- Create table to export ----
subsetSamples <- tableSamples_rR[date >= maxDateModel - 28 & date < maxDateModel]
table_r <- subsetSamples[order(date),
                         .(Group = "Warwick",
                           Model = ifelse(parametersModel$params$linkType == "BB", "GaussianProcessesProp", "GaussianProcessesRaw"),
                           Scenario = "Nowcast",
                           ModelType = "Multiple",
                           Version = "1.0",
                           `Creation Day` = format(Sys.Date(), "%d"), `Creation Month` = format(Sys.Date(), "%m"),`Creation Year` = format(Sys.Date(), "%Y"),
                           `Day of Value` = format(date,"%d"),`Month of Value` = format(date,"%m"),`Year of Value` = format(date,"%Y"),
                           AgeBand = "All",
                           Geography = regionToModel,
                           ValueType = "growth_rate",
                           Value = q0.50,
                           `Quantile 0.05` = q0.05,
                           `Quantile 0.1` = q0.10, `Quantile 0.15` = q0.15,
                           `Quantile 0.2` = q0.20, `Quantile 0.25` = q0.25,
                           `Quantile 0.3` = q0.30, `Quantile 0.35` = q0.35,
                           `Quantile 0.4` = q0.40, `Quantile 0.45` = q0.45,
                           `Quantile 0.5` = q0.50, `Quantile 0.55` = q0.55,
                           `Quantile 0.6` = q0.60, `Quantile 0.65` = q0.65,
                           `Quantile 0.7` = q0.70, `Quantile 0.75` = q0.75,
                           `Quantile 0.8` = q0.80, `Quantile 0.85` = q0.85,
                           `Quantile 0.9` = q0.90, `Quantile 0.95` = q0.95)]
# Table R
table_R <- subsetSamples[order(date),
                         .(Group = "Warwick",
                           Model = ifelse(parametersModel$params$linkType == "BB", "GaussianProcessesProp", "GaussianProcessesRaw"),
                           Scenario = "Nowcast",
                           ModelType = "Multiple",
                           Version = "1.0",
                           `Creation Day` = format(Sys.Date(), "%d"), `Creation Month` = format(Sys.Date(), "%m"),`Creation Year` = format(Sys.Date(), "%Y"),
                           `Day of Value` = format(date,"%d"),`Month of Value` = format(date,"%m"),`Year of Value` = format(date,"%Y"),
                           AgeBand = "All",
                           Geography = regionToModel,
                           ValueType = "R",
                           Value = q0.50RR,
                           `Quantile 0.05` = q0.05RR,
                           `Quantile 0.1` = q0.10RR, `Quantile 0.15` = q0.15RR,
                           `Quantile 0.2` = q0.20RR, `Quantile 0.25` = q0.25RR,
                           `Quantile 0.3` = q0.30RR, `Quantile 0.35` = q0.35RR,
                           `Quantile 0.4` = q0.40RR, `Quantile 0.45` = q0.45RR,
                           `Quantile 0.5` = q0.50RR, `Quantile 0.55` = q0.55RR,
                           `Quantile 0.6` = q0.60RR, `Quantile 0.65` = q0.65RR,
                           `Quantile 0.7` = q0.70RR, `Quantile 0.75` = q0.75RR,
                           `Quantile 0.8` = q0.80RR, `Quantile 0.85` = q0.85RR,
                           `Quantile 0.9` = q0.90RR, `Quantile 0.95` = q0.95RR)]
tableEstimates <- rbind(table_r, table_R)

## ---- Output ----
print(tableEstimates)
ggplot(tableEstimates, aes(x = as.Date(paste0(`Day of Value`, `Month of Value`, `Year of Value`), format = "%d%m%Y"), y = Value)) + theme_minimal() +
  facet_grid(ValueType ~ ., scales = "free_y") +
  geom_hline(data = data.table(ValueType = c("growth_rate", "R"), yi = 0:1), aes(yintercept = yi), linetype = 2, colour = "gray20") +
  geom_ribbon(aes(ymin = `Quantile 0.05`, ymax = `Quantile 0.95`), colour = NA, fill = "gray70", alpha = 0.2) +
  geom_path() +
  #scale_x_date(labels = scales::label_date(formatBreaks), breaks = dateBreaks, expand = c(0,0)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        panel.grid.major = element_line(linetype = 2, colour = "gray90")) +
  labs(x = "day", y = "median (90% CI)")
