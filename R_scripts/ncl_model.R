library(ggplot2)
library(MASS)

#_______________________________________________________________________________
# Looking at each sample ----
# check distribution ----

methods_pa <- read.csv('all_week_pa.csv', row.names <- 1)
ca_vegmet <- read.csv('all_week_meta.csv', row.names <- 1)
row.names(ca_vegmet) <- ca_vegmet$sample

methods_pa$species_richness <- rowSums(methods_pa)
ca_vegmet <- ca_vegmet[match(row.names(methods_pa), row.names(ca_vegmet)), ]
merged_data <- cbind(ca_vegmet, methods_pa[,80])
colnames(merged_data)[5] <- "species_richness"


# Visualize species richness
boxplot(species_richness ~ type, data <- merged_data, 
        main <- "Richness by Type")
boxplot(species_richness ~ habitat, data <- merged_data, 
        main <- "Richness by Habitat")
boxplot(species_richness ~ location, data <- merged_data, 
        main <- "Richness by location")

hist(merged_data$species_richness, breaks <- 20, main <- "Species Richness", 
     xlab <- "Species Richness")
# poisson distribution

mean(merged_data$species_richness)
var(merged_data$species_richness)
#_______________________________________________________________________________
# check model fits ----

# Additive model
model_additive <- glm(species_richness ~ type + habitat, 
                      data <- merged_data, family <- poisson())
# Interaction model
model_interaction <- glm(species_richness ~ type * habitat, 
                            data <- merged_data, , family <- poisson())
# Null model
model_null <- glm(species_richness ~ 1, data <- merged_data, family <- poisson())
# Without habitat
model_type <- glm(species_richness ~ type, data <- merged_data, 
                       family <- poisson())
# Without type
model_habitat <- glm(species_richness ~ habitat, data <- merged_data, 
                    family <- poisson())

model_location_interaction <- glm(species_richness ~ type * location, 
                                 data <- merged_data, , family <- poisson())
model_location_addition <- glm(species_richness ~ type + location, 
                                 data <- merged_data, , family <- poisson())
# Likelihood ratio tests
anova(model_null, model_interaction, test <- "Chisq")  # Additive vs Interaction
anova(model_null, model_additive, test <- "Chisq")      # Full model vs null
anova(model_null, model_type, test <- "Chisq")    # Full vs no habitat
anova(model_null, model_habitat, test <- "Chisq")       # Full vs no type
anova(model_null, model_type, model_habitat, model_additive, model_interaction, test <- "Chisq")

aic_values <- data.frame(
  Model <- c("Null", "Type only", "Habitat only", "Additive", "Interaction"),
  AIC <- c(
    AIC(model_null),
    AIC(model_type),
    AIC(model_habitat),
    AIC(model_additive),
    AIC(model_interaction)
  )
)
print(aic_values)

summary(model_type)
#_______________________________________________________________________________
# diagnostic plots ----
par(mfrow <- c(2,2))
plot(model_type)
par(mfrow <- c(1,1))

# Histogram of residuals
hist(residuals(model_type), breaks <- 10, 
     main <- "Histogram of Residuals", xlab <- "Residuals", 
     xlim <- c(-3, 3))

# Q-Q plot for residuals
qqnorm(residuals(model_type))
qqline(residuals(model_type), col <- "red")

deviance(model_type) / df.residual(model_type)

# % Deviance explained
null_deviance <- deviance(model_null)
residual_deviance <- deviance(model_type)
deviance_explained <- 1 - (residual_deviance / null_deviance)
deviance_explained * 100  # Convert to percentage

# Add predicted values to the data
merged_data$predicted <- predict(model_type, type <- "response")

ggplot(merged_data, aes(x <- type, y <- species_richness)) +
  geom_boxplot() +
  geom_point(aes(y <- predicted), color <- "red")

#_______________________________________________________________________________

# Looking at environmental vairables and variation within method types ----

## factors only relevant to water samples ----
merged_data_W <- subset(merged_data, type == "W")
hist(merged_data_W$species_richness)

flow <- c('low', 'moderate', 'fast')
merged_data_W$water_flow[merged_data_W$water_flow %in% flow] <- 'flowing'
size_l <- c('L', 'M')
merged_data_W$water_size[merged_data_W$water_size %in% size_l] <- 'M/L'
merged_data_W$water_size <- factor(merged_data_W$water_size, 
                        levels <- c("XS", "S", "M/L"), 
                        ordered <- TRUE)

model_waternull <- glm(species_richness ~ 1, data <- merged_data_W, family <- poisson())
model_waterflow <- glm(species_richness ~ water_flow, 
                         data <- merged_data_W, family <- poisson())
model_watersize <- glm(species_richness ~ water_size, 
                       data <- merged_data_W, family <- poisson())

anova(model_waternull, model_waterflow, test <- "Chisq") # no difference, p <- 0.123
anova(model_waternull, model_watersize, test <- "Chisq") # no difference, p <- 0.497

## soil cover data ----

merged_data_no_soilcover <- subset(merged_data, !is.na(soil_cover))
hist(merged_data_no_soilcover$species_richness)

model_soilnull <- glm(species_richness ~ type, data <- merged_data_no_soilcover, family <- poisson())
model_groundcover <- glm(species_richness ~ soil_cover, 
                                     data <- merged_data_no_soilcover, family <- poisson())
model_groundcover_additive <- glm(species_richness ~ type + soil_cover, 
                                     data <- merged_data_no_soilcover, family <- poisson())
model_groundcover_interaction <- glm(species_richness ~ type * soil_cover, 
                       data <- merged_data_no_soilcover, family <- poisson())
AIC(model_soilnull)
AIC(model_groundcover)
AIC(model_groundcover_additive)
AIC(model_groundcover_interaction)
anova(model_soilnull, model_groundcover_additive, test <- "Chisq") # no sig differences

## tree cover data ----
merged_data_no_treecover <- subset(merged_data, !is.na(tree_cover))
hist(merged_data_no_treecover$species_richness)

model_treenull <- glm(species_richness ~ type, data <- merged_data_no_treecover, family <- poisson())
model_tree <- glm(species_richness ~ tree_cover, 
                         data <- merged_data_no_treecover, family <- poisson())
model_tree_additive <- glm(species_richness ~ type + tree_cover, 
                                  data <- merged_data_no_treecover, family <- poisson())
model_tree_interaction <- glm(species_richness ~ type * tree_cover, 
                                     data <- merged_data_no_treecover, family <- poisson())
AIC(model_treenull)
AIC(model_tree)
AIC(model_tree_additive)
AIC(model_tree_interaction)
anova(model_treenull, model_tree_additive, test <- "Chisq")
anova(model_treenull, model_tree_interaction, test <- "Chisq") 
summary(model_tree_additive)
summary(model_tree_interaction)
# no differences BUT close to significant negative influence of low tree cover in additive model (p <- 0.063)
# Significant negative influence of open tree cover on water samples in additive model (p <- 0.007)

#_______________________________________________________________________________
