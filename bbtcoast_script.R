library(tidyverse) # general data handling
library(spaMM) # fitting spatial mixed models
library(DHARMa) # model diagnostics

setwd() # set working directory
data <- read_csv("Final/GitHub/New/bbtcoast_data.csv")

## Choose "PIPNAT" for Pipistrellus nathusii and "EPTNIL" for Eptesicus nilssonii
taxon <- "PIPNAT"

## Choose TRUE to fit interaction of northerness and coastalness
interact <- FALSE


#### Model fitting with spaMM ####
set.seed(0000)
if (interact) {
  modelformula <- bats ~ north*log(coast+1) + Matern(1|x+y) + (1|year/rec.period) + (1|site) 
} else {
  modelformula <- bats ~ north+log(coast+1) + Matern(1|x+y) + (1|year/rec.period) + (1|site) 
}
fits <- list()
test.data <- list()
seasons <- list(1:3, 4:6, 7:9) # recording periods to be included in the 3 subseasons (early, mid and late)
for (i in 1:3){
  test <- filter(data, is.na(obs) == FALSE, species == taxon, rec.period2 %in% seasons[[i]]) %>% 
    mutate(obs.binom = as.numeric(as.logical(obs))) %>%
    group_by(rec.period, site, north, coast, year, rec.period2, x, y) %>%
    summarize(bats = mean(obs.min))%>%
    mutate(bats = as.numeric(as.logical(bats)))
  test.data[[i]] <- test
  fit <- fitme(formula = modelformula, data = test, family = binomial(link = "logit"))
  print(paste0("#### Season ", i, " / ", length(seasons)))
  summary(fit)
  fits[[i]] <- fit
}

#### Model diagnostics ####
i <- 1 # Early season
i <- 2 # Mid season
i <- 3 # Late season

test <- test.data[[i]]
fit <- fits[[i]]
dd <- dist(test[,c("x","y")])
mm <- MaternCorr(dd, nu = fit$corrPars$'1'$nu ,rho= fit$corrPars$'1'$rho)
plot(as.numeric(dd)/1000, as.numeric(mm), xlab = "Distance between pairs of location [in km]", ylab = "Estimated correlation")
sims <- simulateResiduals(fit)
plot(sims)

