#### Setup ####

# Get required packages
#install.packages("tidyverse", "spaMM", "DHARMa", "effects", "sf", "ggpubr", "Cairo")
library(tidyverse)
library(spaMM)
library(DHARMa)
library(effects)
library(sf)
library(ggpubr)
library(Cairo)


# Set working directory and load data (and make a groupng variable for subseasons)
setwd("C:/Data/R Projects/Lukiolaiset/Revision_2024_08")
data <- read_csv("bbtcoast_data.csv") %>%
  mutate(subseason = case_when(rec.period2 %in% c(1,2,3) ~ "Early",
                               rec.period2 %in% c(4,5,6) ~ "Mid",
                               rec.period2 %in% c(7,8,9) ~ "Late"))

# Create sf-object for device locations (for visualisation)
sites <- filter(data, is.na(obs) == FALSE) %>%
  group_by(site, coast, north, year, x, y) %>%
  summarise() %>%
  st_as_sf(coords = c("x", "y"), crs = 3067)

# Measure distances between the sites
dist_all <- as.vector(st_distance(sites))/1000
dist19 <- as.vector(st_distance(filter(sites, year == 2019)))/1000
dist20 <- as.vector(st_distance(filter(sites, year == 2020)))/1000
tibble(mean = c(mean(dist_all), mean(dist19), mean(dist20)), 
       median = c(median(dist_all), median(dist19), median(dist20)), 
       sd = c(sd(dist_all), sd(dist19), sd(dist20)),
       min = c(min(dist_all), min(dist19), min(dist20)),
       max = c(max(dist_all), max(dist19), max(dist20)),
       n = c(length(dist_all), length(dist19), length(dist20)))

# Create spatial object of the study sites (for plotting only)
sites <- group_by(data, site) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  dplyr::select(site, x, y, coast, north) %>%
  st_as_sf(coords = c("x", "y"), crs = 3067)

#### Model fitting ####

# Custom function for fitting task
## Parameters
### df: data frame with all observations
### separating_var: categorical variable in df to split analysis into separate fits using with same parameters
### model_formula: a formula object for fitme() in spaMM. See: ?fitme
### model_family: a family object for fitme() in spaMM. See: ?fitme
### print_summaries: if TRUE, model summaries are printed on the go

fit_models <- function(df, separating_var, model_formula, model_family = binomial(link = "logit"), print_summaries = FALSE) {
  # Check if separating_var is in the data frame
  if (!separating_var %in% names(df)) {
    stop("The separating variable is not present in the data frame.")
  }
  # Check if any variables in the model formula are missing from the data frame
  formula_vars <- all.vars(model_formula)
  missing_vars <- setdiff(formula_vars, names(df))
  if (length(missing_vars) > 0) {
    stop("The following variables from the model formula are missing in the data frame: ", paste(missing_vars, collapse = ", "))
  }
  # Check if the separating variable is included in the model formula
  if (separating_var %in% formula_vars) {
    warning("The separating variable is included in the model formula. It will be excluded from the model fitting.")
  }
  # Get unique levels of the separating variable
  levels <- unique(df[[separating_var]])
  # Print initial message about the number of models to be fitted
  cat("Fitting", length(levels), "models for the levels of", separating_var, "\n")
  # Fit models for each level of the separating variable
  models <- list()
  for (level in levels) {
    cat("Fitting model for", separating_var, "=", level, "\n")
    data_subset <- df[df[[separating_var]] == level, ]
    fit <- fitme(model_formula, family = model_family, data = data_subset)
    models[[as.character(level)]] <- fit
    if (print_summaries) {
      cat("Model summary for", separating_var, "=", level, ":\n")
      print(summary(fit))
    }
  }
  return(models)
}

# Remove NA:s, summarize data by species, site and rec.period, and create a
## binary variable for  bat occurrence
model_data <- filter(data, is.na(obs) == FALSE, species %in% c("EPTNIL", "PIPNAT")) %>%
  mutate(obs.binom = as.numeric(as.logical(obs)), year = as.character(year)) %>%
  group_by(species, subseason, site, north, coast, rec.period, year, x, y) %>%
  summarize(bats = mean(obs.min))%>%
  mutate(bats = as.numeric(as.logical(bats)))
# Set model formula for fit_models()
model_formula = bats ~ year + subseason * north * log(coast + 1) + Matern(1 | x + y) + (1 | site)
# Run fit_models() NOTE: fit_models() could include time calculation
start_time <- Sys.time()
fits <- fit_models(df = model_data, separating_var = "species", model_family  = binomial(link = "logit"), model_formula = model_formula)
print(Sys.time()-start_time)

#write_rds(fits, "fits.rds")
#### Diagnostics ####


##### Null models #####

# Save model formulas

# Null Model 1 (Intercept Only)
null1_model_formula <- bats ~ 1 + Matern(1 | x + y) + (1 | site)

# Null Model 2 (No Interaction)
null2_model_formula <- bats ~ year + subseason + north + log(coast + 1) + Matern(1 | x + y) + (1 | site)

# Null Model 3 (No Subseason)
null3_model_formula <- bats ~ year + north * log(coast + 1) + Matern(1 | x + y) + (1 | site)

# Null Model 4 (No Year)
null4_model_formula <- bats ~ subseason * north * log(coast + 1) + Matern(1 | x + y) + (1 | site)

# Null Model 5 (No Coast)
null5_model_formula <- bats ~ subseason * north + Matern(1 | x + y) + (1 | site)

# Null Model 6 (No North)
null6_model_formula <- bats ~ subseason * log(coast + 1) + Matern(1 | x + y) + (1 | site)


# Run fit_models() NOTE: fit_models() could include time calculation
start_time <- Sys.time()
fits_null1 <- fit_models(df = model_data, separating_var = "species", model_family  = binomial(link = "logit"), model_formula = null1_model_formula)
fits_null2 <- fit_models(df = model_data, separating_var = "species", model_family  = binomial(link = "logit"), model_formula = null2_model_formula)
fits_null3 <- fit_models(df = model_data, separating_var = "species", model_family  = binomial(link = "logit"), model_formula = null3_model_formula)
fits_null4 <- fit_models(df = model_data, separating_var = "species", model_family  = binomial(link = "logit"), model_formula = null4_model_formula)
fits_null5 <- fit_models(df = model_data, separating_var = "species", model_family  = binomial(link = "logit"), model_formula = null5_model_formula)
fits_null6 <- fit_models(df = model_data, separating_var = "species", model_family  = binomial(link = "logit"), model_formula = null6_model_formula)
print(Sys.time()-start_time)


aic_table <- tibble(model = rep(c("full", "null1", "null2", "null3", "null4", "null5", "null6"), times = 2),
                    species = c(rep("PIPNAT", 7), rep("EPTNIL", 7)),
                    aicc = c(AIC(fits[["PIPNAT"]])[2], AIC(fits_null1[["PIPNAT"]])[2], AIC(fits_null2[["PIPNAT"]])[2], AIC(fits_null3[["PIPNAT"]])[2], AIC(fits_null4[["PIPNAT"]])[2], AIC(fits_null5[["PIPNAT"]])[2], AIC(fits_null6[["PIPNAT"]])[2],
                             AIC(fits[["EPTNIL"]])[2], AIC(fits_null1[["EPTNIL"]])[2], AIC(fits_null2[["EPTNIL"]])[2], AIC(fits_null3[["EPTNIL"]])[2], AIC(fits_null4[["EPTNIL"]])[2], AIC(fits_null5[["EPTNIL"]])[2], AIC(fits_null6[["EPTNIL"]])[2]))

# Replace Full model for P. nathusii in fits-object with the "No Year" -null model
fits[["PIPNAT"]] <- fits_null4[["PIPNAT"]]


##### DHARMa #####
#See vignette("DHARMa", package = "DHARMa") for reference

# For P. nathusii
fit <- fits[["PIPNAT"]]
sims <- simulateResiduals(fit)
plot(sims) 

# For E. nilssonii
fit <- fits[["EPTNIL"]]
sims <- simulateResiduals(fit)
plot(sims)


##### Spatial autocorrelation ######

# For P. nathusii
fit <- fits[["PIPNAT"]]
dd <- dist(fit$data[,c("x","y")])
mm <- MaternCorr(dd, nu = fit$corrPars$'1'$nu ,rho= fit$corrPars$'1'$rho) # Matérn spatial correlation
plot(as.numeric(dd)/1000, as.numeric(mm), xlab = "Distance between pairs of location [in km]", ylab = "Estimated correlation") #Plot spatial correlation

# For E. nilssonii
fit <- fits[["EPTNIL"]]
dd <- dist(fit$data[,c("x","y")])
mm <- MaternCorr(dd, nu = fit$corrPars$'1'$nu ,rho= fit$corrPars$'1'$rho) # Matérn spatial correlation
plot(as.numeric(dd)/1000, as.numeric(mm), xlab = "Distance between pairs of location [in km]", ylab = "Estimated correlation") #Plot spatial correlation



##### Confidence intervals #####
# Warning. Calculation of Confidence intervals is very slow. Executing the code 
# below took a moderately powerful laptop approximately 12 hours. Consider using
# a server or excluding parameters with t-values between -1.96 and 1.96 as they 
# will include zero (selected effects are set with the parameter "parm").

#ci_pipnat <- read_rds("new_ci_pipnat.rds")
#ci_eptnil <- read_rds("ci_eptnil.rds")

results_pipnat <- read_rds("results_pipnat.rds") %>%
  mutate(Lower_Cl = ci_pipnat[,"2.5 %"], Upper_Cl = ci_pipnat[, "97.5 %"])
results_eptnil <- read_rds("results_eptnil.rds")


start_time <- Sys.time()

# For P. nathusii
model_family <- fits[["PIPNAT"]]$family
# set parm = c(4,6) to exclude parameters with t-values between -1.96 and 1.96
ci_pipnat <- confint.HLfit(fits[["PIPNAT"]], parm = 9:12, format = "stats")
write_rds(ci_pipnat, "new_ci_pipnat.rds")

print(Sys.time()-start_time)

# For E. nilssonii
model_family <- fits[["EPTNIL"]]$family
# set parm = c(1:3, 5, 7:8, 11:12) to exclude parameters with t-values between -1.96 and 1.96
ci_eptnil <- confint.HLfit(fits[["EPTNIL"]], parm = 1:13, format = "stats") 

##### Results as table #####
results_pipnat <- tibble(Fixed_effect = names(fits[["PIPNAT"]]$fixef), 
                         Estimate = fits[["PIPNAT"]]$fixef, 
                         T_value = summary(fits[["PIPNAT"]])$beta_table[,"t-value"], 
                         Lower_Cl = ci_pipnat[,"2.5 %"], 
                         Upper_Cl = ci_pipnat[, "97.5 %"])

results_eptnil <- tibble(Fixed_effect = names(fit_en$fixef), 
                  Estimate = fit_en$fixef, 
                  T_value = summary(fit_en)$beta_table[,"t-value"], 
                  Lower_Cl = ci_eptnil[,"2.5 %"], 
                  Upper_Cl = ci_eptnil[, "97.5 %"])

#write_csv(results_pipnat, "results_pipnat.csv")
#write_csv(results_eptnil, "results_eptnil.csv")


#### Visualisation ####


##### Estimate maps #####

finland <- st_read("mml_maakunnat_2021/maakunnat_2021_milj.shp") %>%
  st_union()
finland <- st_simplify(finland, dTolerance = 2000)

baltic <- list(rbind(c(74601.95,6632492.46), c(543717,6632492.46), c(543717,7304430), c(74601.95,7304430), c(74601.95,6632492.46)))
baltic <- st_sfc(st_polygon(baltic), crs = st_crs(finland))
baltic <- st_difference(baltic, finland)

# Create line for for the 60 deg latitude
lat60 <- st_linestring(matrix(c(seq(19,31,0.1),rep(60,121)),121,2), dim = "XY") %>%
  st_sfc(crs = 4326, dim = "XY") %>%
  st_transform(crs = 3067)

# Create a shorter version for visualization
lat60s <- st_linestring(matrix(c(seq(20.7,30,0.1),rep(60,94)),94,2), dim = "XY") %>%
  st_sfc(crs = 4326, dim = "XY") %>%
  st_transform(crs = 3067)

# Function to modify the y coordinates by adding a specified value
modify_y_coords <- function(lat60, y_offset) {
  # Extract the coordinates from lat60
  lat60_coords <- st_coordinates(lat60)
  
  # Add the offset to the y coordinates (second column)
  lat60_coords[, 2] <- lat60_coords[, 2] + y_offset
  
  # Recreate the LINESTRING with the modified coordinates
  lat60_modified <- st_sfc(st_linestring(lat60_coords), crs = st_crs(lat60))
  
  # Create an sf object with the same structure as lat60
  return(st_sf(geometry = lat60_modified))
}

# Example vector of values to add to the y coordinates
y_offsets <- c(9000, 200000, 300000, 500000, 600000)  # Replace with desired offsets

# Loop over the vector of values and create a new line for each
modified_lines <- lapply(y_offsets, function(offset) modify_y_coords(lat60s, offset))

# Optional: Combine all the lines into one sf object (if needed)
combined_sf <- do.call(rbind, modified_lines)


map_grid <- st_make_grid(sites, n = c(25,25))
map_grid <- map_grid [st_intersects(map_grid , finland, sparse = FALSE)]
grid_centroids <- st_centroid(map_grid)
map_grid  <- st_sf(map_grid , tibble(coast = as.vector(st_distance(grid_centroids, baltic))/1000,
                                     coast.log = log(coast+1),
                                     north = as.vector(st_distance(grid_centroids,lat60))/1000,
                                     cell_id = seq(1, length(map_grid), 1)))


# Create a data frame of estimates for the map_grid (without year for PIPNAT)
subseasons <- c("Early", "Mid", "Late")
grid_estimates <- list()
for (i in "PIPNAT"){
  fit <- fits[[i]]
  # Create a long-format data frame for all combinations of year and sub-season
  estimates <- expand.grid(subseason = subseasons) %>%
    left_join(map_grid, by = character()) %>%  # Cartesian join to get all combinations with grid
    rowwise() %>%
    mutate(est = fit$fixef["(Intercept)"] +
             fit$fixef["subseasonMid"] * (subseason == "Mid") +
             fit$fixef["subseasonLate"] * (subseason == "Late") +
             fit$fixef["north"] * north +
             fit$fixef["log(coast + 1)"] * coast.log +
             fit$fixef["subseasonMid:north"] * (subseason == "Mid") * north +
             fit$fixef["subseasonLate:north"] * (subseason == "Late") * north +
             fit$fixef["subseasonMid:log(coast + 1)"] * (subseason == "Mid") * coast.log +
             fit$fixef["subseasonLate:log(coast + 1)"] * (subseason == "Late") * coast.log +
             fit$fixef["north:log(coast + 1)"] * north * coast.log +
             fit$fixef["subseasonMid:north:log(coast + 1)"] * (subseason == "Mid") * north * coast.log +
             fit$fixef["subseasonLate:north:log(coast + 1)"] * (subseason == "Late") * north * coast.log,
           species = i,
           prob.est = exp(est)/(1+exp(est))) %>%
    ungroup()
  grid_estimates[[i]] <- estimates
}
estimate_grid_pnat <- bind_rows(grid_estimates) %>%
  st_as_sf(sf_column_name = "map_grid")

# Plot spatial estimates of P. nathusii occurrence
map_pipnat <-ggplot() +
  geom_sf(data = modified_lines[[1]], color = "#d59c00", size = 1, linetype = "solid") +
  geom_sf(data = modified_lines[[2]], color = "#ffcd41", size = 1, linetype = "solid") +
  geom_sf(data = modified_lines[[3]], color = "#ffe68b", size = 1, linetype = "solid") +
  geom_sf(data = modified_lines[[4]], color = "#d583cd", size = 1, linetype = "solid") +
  geom_sf(data = modified_lines[[5]], color = "#83416a", size = 1, linetype = "solid") +
  geom_sf(aes(fill = prob.est), data = filter(estimate_grid_pnat, species == "PIPNAT"), color = "NA", alpha = 1) +
  geom_sf(data = finland, alpha = 0, size = 0.7, color = "grey50", fill = "NA") +
  geom_sf(data = sites, shape = 16,  size = 1, color = "black") +
  scale_fill_viridis_c(begin = 0.2, direction = -1, name = "Prob. of\noccurrence") +
  coord_sf(ylim = c(6670412.4,7277589.0)) +
  facet_grid(.~subseason)+
  theme_void() + theme(panel.spacing = unit(-2, "lines"),
                       legend.position = "right",
                       strip.text = element_text(size = 10, margin = margin(2,0,2,0, "mm")),
                       strip.background = element_rect(colour = "white", fill = "white", size = 7),
                       legend.title = element_text(size = 10))


# Create a data frame of estimates for the map_grid
years <- c("2019", "2020")
subseasons <- c("Early", "Mid", "Late")
grid_estimates <- list()
for (i in "EPTNIL"){
  fit <- fits[[i]]
  # Create a long-format data frame for all combinations of year and sub-season
  estimates <- expand.grid(year = years, subseason = subseasons) %>%
    left_join(map_grid, by = character()) %>%  # Cartesian join to get all combinations with grid
    rowwise() %>%
    mutate(est = fit$fixef["(Intercept)"] +
             fit$fixef["year2020"] * (year == "2020") +
             fit$fixef["subseasonMid"] * (subseason == "Mid") +
             fit$fixef["subseasonLate"] * (subseason == "Late") +
             fit$fixef["north"] * north +
             fit$fixef["log(coast + 1)"] * coast.log +
             fit$fixef["subseasonMid:north"] * (subseason == "Mid") * north +
             fit$fixef["subseasonLate:north"] * (subseason == "Late") * north +
             fit$fixef["subseasonMid:log(coast + 1)"] * (subseason == "Mid") * coast.log +
             fit$fixef["subseasonLate:log(coast + 1)"] * (subseason == "Late") * coast.log +
             fit$fixef["north:log(coast + 1)"] * north * coast.log +
             fit$fixef["subseasonMid:north:log(coast + 1)"] * (subseason == "Mid") * north * coast.log +
             fit$fixef["subseasonLate:north:log(coast + 1)"] * (subseason == "Late") * north * coast.log,
           species = i,
           prob.est = exp(est)/(1+exp(est))) %>%
    ungroup()
  grid_estimates[[i]] <- estimates
}

average_years <- bind_rows(grid_estimates) %>%
  group_by(subseason, species, cell_id) %>%
  summarize(ave.prob.est = mean(prob.est))

estimate_grid_enil <- bind_rows(grid_estimates) %>%
  st_as_sf(sf_column_name = "map_grid") %>%
  left_join(average_years, by = c("subseason", "species", "cell_id"))



# Plot spatial estimates of E. nilssonii occurrence
map_eptnil <- ggplot() +
  geom_sf(data = modified_lines[[1]], color = "#d59c00", size = 1, linetype = "solid") +
  geom_sf(data = modified_lines[[2]], color = "#ffcd41", size = 1, linetype = "solid") +
  geom_sf(data = modified_lines[[3]], color = "#ffe68b", size = 1, linetype = "solid") +
  geom_sf(data = modified_lines[[4]], color = "#d583cd", size = 1, linetype = "solid") +
  geom_sf(data = modified_lines[[5]], color = "#83416a", size = 1, linetype = "solid") +
  geom_sf(aes(fill = ave.prob.est), data = filter(estimate_grid_enil, species == "EPTNIL", year == 2019), color = NA, alpha = 1) +
  geom_sf(data = finland, alpha = 0, size = 0.4, color = "grey50", fill = "NA") +
  geom_sf(data = sites, shape = 16, size = 1, color = "black") +
  scale_fill_viridis_c(begin = 0.2, direction = -1, name = "Prob. of\noccurrence") +
  coord_sf(ylim = c(6670412.4,7277589.0)) +
  facet_grid(.~subseason)+
  theme_void() + theme(panel.spacing = unit(-2, "lines"),
                       legend.position = "right",
                       strip.text = element_text(size = 10, margin = margin(2,0,2,0, "mm")),
                       strip.background = element_rect(colour = "white", fill = "white", size = 7),
                       legend.title = element_text(size = 10))


##### Line graphs #####
eff <- list()
for (i in names(fits)){
  fit <- fits[[i]]
  data_subset <- fit$data
  eff[[i]] <- as_tibble(Effect(c("coast", "north", "subseason"), fits[[i]])) %>%
  mutate(north = as.factor(north), species = i)
}
eff <- bind_rows(eff) %>%
  mutate(subseason = factor(subseason, levels = c("Early", "Mid", "Late")))



lineplot_pipnat <- ggplot(data=filter(eff, species == "PIPNAT"))+#, north %in% c(9,300,600))) + 
  geom_errorbar(aes(y=fit, x = coast, ymin=lower, ymax=upper, color = north), size = 0.75, position=position_dodge(10)) +
  geom_line(aes(coast, fit, color=north), size = 1) + 
  geom_point(aes(coast, fit, color=north), size = 2.5) + 
  scale_x_continuous(breaks = c(0, 100, 200), name = "Distance to coast (km)") +
  scale_y_continuous(breaks = c(0.0, 0.2, 0.4, 0.6), name = "Prob. of occurrence") +
  scale_color_manual(values = c("#d59c00", "#ffcd41", "#ffe68b", "#d583cd", "#83416a"), name = "Distance\nto 60°N\n(km)")+
  coord_cartesian(ylim = c(0, 0.6)) +
  facet_grid(.~subseason)+
  theme_minimal() + theme(panel.grid.major = element_line(linetype = "solid", size = 0.7, colour = "gray85"),
                          panel.grid.minor = element_line(linetype = "dotted", size = 0.7, colour = "gray85"),
                          panel.spacing = unit(1, "lines"),
                          legend.position = "right",
                          legend.text = element_text(size = 9),
                          strip.text = element_text(size = 10, margin = margin(2,0,2,0, "mm")),
                          strip.background = element_rect(colour = "white", fill = "white", size = 7),
                          legend.title = element_text(size = 10),
                          axis.title = element_text(size = 10),
                          axis.text = element_text(size = 9))

lineplot_eptnil <- ggplot(data=filter(eff, species == "EPTNIL"))+#, north %in% c(9,300,600))) + 
  geom_errorbar(aes(y=fit, x = coast, ymin=lower, ymax=upper, color = north), size = 0.75, position=position_dodge(10)) +
  geom_line(aes(coast, fit, color=north), size = 1) + 
  geom_point(aes(coast, fit, color=north), size = 2.5) +
  scale_x_continuous(breaks = c(0, 100, 200), labels = c("0", "100", "200"), name = "Distance to coast (km)") +
  scale_y_continuous(name = "Prob. of occurrence") +
  scale_color_manual(values = c("#d59c00", "#ffcd41", "#ffe68b", "#d583cd", "#83416a"), name = "Distance\nto 60°N\n(km)")+
  facet_grid(.~subseason)+
  theme_minimal() + theme(panel.grid.major = element_line(linetype = "solid", size = 0.7, colour = "gray85"),
                          panel.grid.minor = element_line(linetype = "dotted", size = 0.7, colour = "gray85"),
                          panel.spacing = unit(1, "lines"),
                          legend.position = "right",
                          strip.text = element_text(size = 10, margin = margin(2,0,2,0, "mm")),
                          strip.background = element_rect(colour = "white", fill = "white", size = 7),
                          legend.title = element_text(size = 10),
                          legend.text = element_text(size = 9),
                          axis.title = element_text(size = 10),
                          axis.text = element_text(size = 9))


##### Arrange and save figures #####

# P. nathusii
ggpubr::ggarrange(map_pipnat, lineplot_pipnat, ncol = 1, nrow = 2, labels = "auto")
ggsave("new_results_pipnat.png", device = "png", type = "cairo-png", width = 170, height = 140, units = "mm", dpi = 600, bg = "white")

# E. nilssonii
ggpubr::ggarrange(map_eptnil, lineplot_eptnil, ncol = 1, nrow = 2, labels = "auto")
ggsave("new_results_eptnil.png", device = "png", type = "cairo-png", width = 170, height = 140, units = "mm", dpi = 600, bg = "white")

