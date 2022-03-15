# Load libraries
library(ggplot2)
library(ggrepel)
library(ggforce)
library(cowplot) 
library(dplyr)
library(brms)
library(rstan)
library(bayestestR)
library(rotl)

# Functions
mytheme <- function() {
  theme_bw() + 
    theme(panel.border          = element_rect(fill = NA, colour = "black"), # set border around plot.
          panel.grid.major      = element_blank(), # remove major grid lines
          panel.grid.minor      = element_blank(), # remove minor grid lines
          axis.line             = element_blank(), # remove axis lines
          axis.ticks            = element_line(colour = "black"),
          axis.text             = element_text(size = 10, colour = "black"), # axis text size
          axis.title            = element_text(size = 10), # axis title size
          axis.title.y          = element_text(vjust = 3), # increase distance from the y-axis
          axis.title.x          = element_text(vjust = -1), # increase distance from the x-axis
          panel.background      = element_rect(fill = NA),
          plot.background       = element_rect(fill = NA, color = NA), # remove background colour
          plot.margin           = unit(c(1, 1, 1, 1), units = , "cm"), 
          legend.background     = element_rect(fill = NA, color = NA), # get rid of legend bg
          legend.box.background = element_rect(fill = NA, color = NA), # get rid of legend panel bg
          strip.text.x          = element_text(size = 10, color = "black", face = "bold"), # for facet plots
          strip.background      = element_rect(fill = NA, color = NA)
          )
} # set up plot theme

# Set options in Rstan
rstan::rstan_options(auto_write = TRUE) # translate to STAN platform for running Bayesian model
options(mc.cores = parallel::detectCores()) # detects how many cores available to use

# Set directory
setwd('/Users/nicholaswu/Dropbox/Chytrid meta-analysis') # Macbook

# Load data
raw_data <- read.csv("trait_raw_data.csv") 
es_data  <- read.csv("trait_corr_data.csv") 

# Cleaning raw data
traits_data_clean <- raw_data %>% 
  dplyr::filter(trait != "Survival") %>% 
  dplyr::mutate(species_OTL = as.factor(species_OTL),
                trait_value = as.numeric(trait_value),
                resistance  = as.factor(resistance),
                trait       = as.factor(trait),
                response    = as.factor(response),
                lnBd        = log(Bd + 1),
                bio_hier    = factor(bio_hier, levels = c("Cellular","Tissue","Organism")),
                temp_K      = trt_temp + 273.15, # convert temperature Celcius to Kelvin
                inv_temp    = 1 / 8.62e-5  * (1 / mean(temp_K, na.rm = TRUE) - 1 / temp_K) # standardise to mean temp, Boltzmann constant as eV/K
  ) %>%
  dplyr::select(-one_of(c("notes", "title", "ref")))

# Clean effect size data
es_data_clean <- es_data %>% 
  dplyr::mutate(study_ID    = ifelse(notes == "repeated", study_ID, study_ID + max(traits_data_clean$study_ID)),
                species_OTL = as.factor(species_OTL),
                lnBd        = log(Bd + 1),
                resistance  = as.factor(resistance),
                trait       = as.factor(trait),
                response    = as.factor(response),
                bio_hier    = factor(bio_hier, levels = c("Cellular","Tissue","Organism")),
                Zr          = (1 / 2) * log((1 + corr_coeff) / (1 - corr_coeff)), # Fisher-transformed correlation coefficient
                Zr_v        = 1 / (sample_size - 3), # sampling variance (v)
                Zr_sei      = 1 / sqrt(sample_size - 3), # standard error (SE)
                Zr_inv      = 1 / Zr_sei, # precision (inverse of SE) 
                Zr_z        = Zr / Zr_sei, # Egger - z score 
                Zr_w        = 1 / Zr_v) %>% # weight (inverse of Zr_v))
  dplyr::select(-one_of(c("notes", "title", "ref")))

# Clean survival data
surv_data_clean <- raw_data %>% 
  dplyr::filter(trait == "Survival") %>% 
  dplyr::mutate(species_OTL = as.factor(species_OTL),
                resistance  = as.factor(resistance),
                life_stage  = as.factor(life_stage),
                lnBd        = log(Bd + 1),
                lnMass      = log(mass_g),
                lnTime      = log(time_post_days),
                temp_K      = trt_temp + 273.15, # convert temperature Celcius to Kelvin
                inv_temp    = 1 / 8.62e-5  * (1 / mean(temp_K, na.rm = TRUE) - 1 / temp_K), # standardise to mean temp, Boltzmann constant as eV/K
                alive       = ifelse(trait_value == "Alive", 1, 0)) %>%
  tidyr::unite("species_lifestage", c("species", "life_stage"), sep = "_", remove = FALSE) %>%
  dplyr::select(-one_of(c("notes", "title", "ref"))) 

# Infection intensity distribution
ggplot(traits_data_clean %>% dplyr::filter(exposed != 0), aes(x = lnBd, colour = resistance, fill = resistance)) +
  geom_density(alpha = 0.2) +
  scale_colour_manual(values = c("#FB8D46", "#CD5A7E", "#159291")) +
  scale_fill_manual(values = c("#FB8D46", "#CD5A7E", "#159291")) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  xlab("Bd load (log ZE + 1)") +
  ylab("Response density") +
  mytheme()

## EFFECT SIZE CALCULATION ## ----------------------------------------------------------------------
effect_size_data <- as.data.frame(traits_data_clean %>%
                                    dplyr::group_by(study_ID, bio_hier, trait, species, response) %>%
                                    dplyr::summarise(corr_coeff     = cor(Bd, trait_value, method = "pearson"),
                                                     sample_size    = length(trait_value),
                                                     author_first   = unique(author_first),
                                                     year           = median(year),
                                                     extracted      = unique(extracted),
                                                     Order          = unique(Order),
                                                     Family         = unique(Family),
                                                     species        = unique(species),
                                                     species_OTL    = unique(species_OTL),
                                                     resistance     = unique(resistance),
                                                     origin         = unique(origin),
                                                     setting        = unique(setting),
                                                     life_stage     = unique(life_stage),
                                                     mean_mass_g    = mean(mass_g),
                                                     trt_temp       = mean(trt_temp),
                                                     Bd             = mean(Bd[Bd != 0]),
                                                     time_post_days = median(time_post_days[Bd != 0]),
                                                     lnBd           = mean(lnBd[lnBd != 0]))) %>%
  dplyr::mutate(Zr     = (1 / 2) * log((1 + corr_coeff) / (1 - corr_coeff)), # Fisher-transformed correlation coefficient
                Zr_v   = 1 / (sample_size - 3), # sampling variance (v)
                Zr_sei = 1 / sqrt(sample_size - 3), # standard error (SE)
                Zr_inv = 1 / Zr_sei, # precision (inverse of SE) 0
                Zr_z   = Zr / Zr_sei, # Egger - z score 
                Zr_w   = 1 / Zr_v, # weight (inverse of Zr_v)
                Zr     = ifelse(response %in% c("Evaporative water loss", "AQP activity", "Skin permeability", 
                                                'Cutaneous ion loss', 'Muscle water content', "Caspase activity"), -abs(Zr), Zr)) 

# Combine effect size dataset
effect_size_full <- bind_rows(effect_size_data, es_data_clean) %>%
  tibble::rowid_to_column("es_ID") %>% # add effect size ID
  dplyr::mutate(year_centre = year - mean(year)) # mean-centring year of publication)

## DATA SUMMARY ## ---------------------------------------
unique(effect_size_full$trait)
length(unique(effect_size_full$study_ID))
nrow(effect_size_full)

table_1 <- data.frame(effect_size_full %>% 
  dplyr::group_by(trait, response) %>% 
  dplyr::summarise(ef_n = n(),
                   study_n = length(unique(study_ID)),
                   species_n = length(unique(species))))

as.data.frame(effect_size_full %>% 
                dplyr::group_by(life_stage) %>% 
                dplyr::summarise(ef_n = n())) %>%
  dplyr::arrange(-ef_n)

as.data.frame(effect_size_full %>% 
                dplyr::group_by(species) %>% 
                dplyr::summarise(study_n = length(unique(study_ID)))) %>%
  dplyr::arrange(-study_n)

as.data.frame(effect_size_full %>% 
                dplyr::group_by(trait) %>% 
                dplyr::summarise(ef_n = n()))

## PHYLOGENETIC TREE ## -------------------------------------------------------------------------
species_comb <- rbind(effect_size_full %>% dplyr::select(Order, Family, species, species_OTL, resistance), 
                      surv_data_clean %>% dplyr::select(Order, Family, species, species_OTL, resistance)) # 'surv_data_clean' created below
species_all  <- sort(unique(as.character(species_comb$species_OTL))) # generate list of species (as character format)
taxa         <- rotl::tnrs_match_names(names = species_all) # match taxonomic names to the OTL

# check if species list match OT identifier
taxa[taxa$approximate_match == TRUE,] # none so far

# retrieving phylogenetic relationships among taxa in the form of a trimmed sub-tree
tree <- rotl::tol_induced_subtree(ott_ids = rotl::ott_id(taxa), label_format = "name")
plot(tree, cex = 0.6, label.offset = 0.1, no.margin = TRUE) # plot tree

# Compute branch lengths
set.seed(1) 
tree <- ape::compute.brlen(tree, method = "Grafen", power = 1)
tree <- ape::multi2di(tree, random = TRUE) # use a randomization approach to deal with polytomies

# Check tree is ultrametric
ape::is.ultrametric(tree) # TRUE

# Create correlation matrix for analysis
phylo_cor <- ape::vcv(tree, cor = T)

# Fig production
tree_tip_label <- tree$tip.label # extract tree tip names
species_list   <- levels(species_comb$species_OTL) # extract species name 

# Check if lengths match for both data and tree
length(unique(tree_tip_label)) # 57
length(unique(species_list)) # 57

effect_size_full$species <- sub(" ", "_", effect_size_full$species_OTL)
effect_size_full$species <- factor(effect_size_full$species, levels = tree_tip_label) # relevel order by tree
species_comb$species <- sub(" ", "_", species_comb$species_OTL)
species_comb$species <- factor(species_comb$species, levels = tree_tip_label) # relevel order by tree

species_trait_data       <- species_comb %>% 
  dplyr::distinct(species) %>% 
  dplyr::mutate(Family     = as.factor(species_comb$Family[match(species, species_comb$species)]),
                resistance = as.factor(species_comb$resistance[match(species, species_comb$species)]),
                Order      = as.factor(species_comb$Order[match(species, species_comb$species)])) %>% 
  tibble::column_to_rownames(var = 'species')

levels(species_trait_data$Family)
mycol <- viridis::viridis(11) # set 9 discrete colours
diversitree::trait.plot(tree, species_trait_data,
                        cols = list(Family = mycol, resistance = c("#FB8D46", "#CD5A7E", "#159291")),
                        type = 'p', cex.lab = 0.7, w = 0.05)

## META-ANALYSIS ## ----------------------------------------------------------------------
# OVERALL EFFECT #
meta_prior    <- c(set_prior("cauchy(0, 1)", class = "sd")) # Cauchy on tau (random effect variance), normal on fixed effect
overall_model <- brms::brm(Zr | se(Zr_sei) ~  time_post_days + life_stage + origin + resistance + trt_temp + year_centre + 
                             (1 | es_ID) + (1 | study_ID) + (1 | species_OTL) + (1 | gr(species, cov = phylo)),
                           data    = effect_size_full,
                           family  = gaussian,
                           data2   = list(phylo = phylo_cor),
                           prior   = meta_prior,
                           iter    = 1e4, warmup = 5e3, cores = 4, chains = 4,
                           control = list(adapt_delta = 0.999, max_treedepth = 18))

# Check convergence
brms::pp_check(overall_model)

# Model summary
summary(overall_model)

# Heterogeneity
overall_post <- brms::posterior_samples(overall_model) # extracting the posterior distributions from our models

# WI = weight
Zr_WI <- na.omit(effect_size_full$Zr_w)

# s2I = measurement error variance = sigma2m
s2I_Zr <- sum(Zr_WI[is.finite(Zr_WI)] * (length(Zr_WI) - 1)) / (sum(Zr_WI[is.finite(Zr_WI)]) ^ 2 - sum(Zr_WI[is.finite(Zr_WI)] ^ 2))

# total variance, including measurement error variance
total_var_Zr <- overall_post$sd_es_ID__Intercept +
  overall_post$sd_species__Intercept +
  overall_post$sd_species_OTL__Intercept +
  overall_post$sd_study_ID__Intercept +
  s2I_Zr

# total heterogeneity I2
I2_total_Zr     <- (total_var_Zr - s2I_Zr) / total_var_Zr
I2_total_Zr_est <- round(c(MCMCglmm::posterior.mode(I2_total_Zr),
                           bayestestR::hdi(I2_total_Zr, ci = 0.95)$CI_low,
                           bayestestR::hdi(I2_total_Zr, ci = 0.95)$CI_high), 3) * 100

# observational level I2
I2_esID_Zr     <- overall_post$sd_es_ID__Intercept / total_var_Zr
I2_esID_Zr_est <- round(c(MCMCglmm::posterior.mode(I2_esID_Zr),
                          bayestestR::hdi(I2_esID_Zr, ci = 0.95)$CI_low,
                          bayestestR::hdi(I2_esID_Zr, ci = 0.95)$CI_high), 3) * 100

# studyID I2
I2_studyID_Zr     <- overall_post$sd_study_ID__Intercept / total_var_Zr
I2_studyID_Zr_est <- round(c(MCMCglmm::posterior.mode(I2_studyID_Zr),
                             bayestestR::hdi(I2_studyID_Zr, ci = 0.95)$CI_low,
                             bayestestR::hdi(I2_studyID_Zr, ci = 0.95)$CI_high), 3) * 100

# phylogeny I2: notice that s2I is substracted from this calculation as phylogenetic
# relatedness is a "fixed random effect"
I2_phylo_Zr     <- overall_post$sd_species__Intercept / (total_var_Zr - s2I_Zr)
I2_phylo_Zr_est <- round(c(MCMCglmm::posterior.mode(I2_phylo_Zr),
                           bayestestR::hdi(I2_phylo_Zr, ci = 0.95)$CI_low,
                           bayestestR::hdi(I2_phylo_Zr, ci = 0.95)$CI_high), 3) * 100

# speciesID I2
I2_species_Zr     <- overall_post$sd_species_OTL__Intercept / total_var_Zr
I2_species_Zr_est <- round(c(MCMCglmm::posterior.mode(I2_species_Zr),
                             bayestestR::hdi(I2_species_Zr, ci = 0.95)$CI_low,
                             bayestestR::hdi(I2_species_Zr, ci = 0.95)$CI_high), 3) * 100

# Funnel plot
metafor::funnel(x = effect_size_full$Zr, sei = effect_size_full$Zr_sei, pch = 1)

# TRAIT MODEL
trait_model <- brms::brm(Zr | se(Zr_sei) ~ -1 + trait + resistance + life_stage + year_centre + 
                           (1 | es_ID) + (1 | study_ID) + (1 | species_OTL) + (1 | gr(species, cov = phylo)),
                         data    = effect_size_full %>%
                           dplyr::group_by(trait) %>% 
                           dplyr::filter(n() >= 5),
                         family  = gaussian,
                         data2   = list(phylo = phylo_cor),
                         prior   = meta_prior,
                         iter    = 1e4, warmup = 5e3, cores = 4, chains = 4,
                         control = list(adapt_delta = 0.999, max_treedepth = 18))

brms::fixef(trait_model)
summary(trait_model)

# Extract marginal effects from model
model_me      <- brms::conditional_effects(trait_model, c("trait", "resistance", "life_stage"))
trait_me      <- as.data.frame(model_me[[1]]) %>% dplyr::rename(estimate = estimate__, ci.lb = lower__, ci.ub = upper__)
resistance_me <- as.data.frame(model_me[[2]]) %>% dplyr::rename(estimate = estimate__, ci.lb = lower__, ci.ub = upper__)
age_me        <- as.data.frame(model_me[[3]]) %>% dplyr::rename(estimate = estimate__, ci.lb = lower__, ci.ub = upper__)

# Plot model output
trait_plot <- trait_me %>% 
  ggplot(aes(x = trait, y = estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  ggforce::geom_sina(data = effect_size_full, aes(x = trait, y = Zr, size = Zr_inv), colour = "#cfcfcf") +
  geom_point(aes(colour = trait), size = 4, show.legend = FALSE) +
  geom_errorbar(aes(ymin = ci.lb, ymax = ci.ub, colour = trait), size = 0.8, width = 0.1, show.legend = FALSE) +
  colorspace::scale_colour_discrete_sequential(palette = "RedOr", nmax = 14, order = 4:14) +
  xlab(NULL) + ylab(expression("Effect size"~(italic(Z)[r]))) +
  coord_flip() +
  mytheme() + theme(legend.position = "bottom")

resist_plot <- resistance_me %>% 
  ggplot(aes(x = resistance, y = estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  ggforce::geom_sina(data = effect_size_full, aes(x = resistance, y = Zr, size = Zr_inv), colour = "#cfcfcf") +
  geom_point(aes(colour = resistance), size = 4, show.legend = FALSE) +
  geom_errorbar(aes(ymin = ci.lb, ymax = ci.ub, colour = resistance), size = 0.8, width = 0.1, show.legend = FALSE) +
  scale_colour_manual(values = c("#FB8D46", "#CD5A7E", "#159291")) +
  xlab(NULL) + ylab(expression("Effect size"~(italic(Z)[r]))) +
  coord_flip() +
  mytheme() + theme(legend.position = "bottom")

age_plot <- age_me %>% 
  ggplot(aes(x = life_stage, y = estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  ggforce::geom_sina(data = effect_size_full, aes(x = life_stage, y = Zr, size = Zr_inv), colour = "#cfcfcf") +
  geom_point(aes(colour = life_stage), size = 4, show.legend = FALSE) +
  geom_errorbar(aes(ymin = ci.lb, ymax = ci.ub, colour = life_stage), size = 0.8, width = 0.1, show.legend = FALSE) +
  colorspace::scale_colour_discrete_sequential(palette = "Blues 3", nmax = 5, order = c(3,5)) +
  xlab(NULL) + ylab(expression("Effect size"~(italic(Z)[r]))) +
  coord_flip() +
  mytheme() + theme(legend.position = "bottom")

# Combine plots
right_row <- cowplot::plot_grid(resist_plot + theme(legend.position = "none"), age_plot +  theme(legend.position = "none"), labels = c('b', 'c'), ncol = 1)
cowplot::plot_grid(trait_plot, right_row, labels = 'a', ncol = 2, rel_widths = c(1, 0.85))

## ESTIMATE BD SENSITIVITY ## ---------------------------------------------------------------------------------- 
# Extract mean, 10th and 90th quantile control values by species and responses
control_mean <- as.data.frame(traits_data_clean %>%
  dplyr::filter(Bd == 0) %>% 
  dplyr::group_by(species, response) %>% 
  dplyr::summarise(mean            = mean(trait_value, na.rm = TRUE),
                   quantile_10     = quantile(trait_value, .10),
                   quantile_90     = quantile(trait_value, .90)) %>%
  tidyr::unite("specific_response", c(species, response), sep = "_",  remove = FALSE)) # unique identifier by species and response

# Combine mean controls and calculate relative change of response
traits_exposed <- traits_data_clean %>% 
  tidyr::unite("specific_response", c(species, response), sep = "_",  remove = FALSE) %>%
  dplyr::mutate(control           = control_mean$mean[match(specific_response, control_mean$specific_response)], 
                change            = ifelse(trait_unit == "relative", trait_value, (trait_value - control) / control), # calculate relative change
                specific_response = as.factor(specific_response),
                quantile_90       = control_mean$quantile_90[match(specific_response, control_mean$specific_response)],
                quantile_10       = control_mean$quantile_10[match(specific_response, control_mean$specific_response)],
                control_90        = ifelse(trait_unit == "relative", trait_value, (quantile_90 - control) / control), # calculate relative 90th quantile
                control_10        = ifelse(trait_unit == "relative", trait_value, (quantile_10 - control) / control)) # calculate relative 10th quantile

# Run lm models to find ones with good R2 
R_sq_model <- as.data.frame(traits_exposed %>% 
                              dplyr::filter(!is.na(change)) %>% # remove NA's
                              dplyr::group_by(specific_response) %>%
                              dplyr::do(model = broom::glance(lm(change ~ lnBd, data = .))) %>%
                              tidyr::unnest(model)) %>%
  dplyr::select(specific_response, r.squared, adj.r.squared, p.value)

# extract equation (intercept and slope)
lm_eq <- as.data.frame(traits_exposed %>% 
                         dplyr::filter(!is.na(change)) %>%
                         dplyr::group_by(specific_response) %>%
                         dplyr::do(model = broom::tidy(lm(change ~ lnBd, data = .))) %>% 
                         tidyr::unnest(model)) %>%
  dplyr::select(specific_response, term, estimate) %>%
  reshape(idvar = "specific_response", timevar = "term", direction = "wide") %>% # from long to wide for estimate column
  dplyr::rename(intercept = "estimate.(Intercept)",
                slope     = "estimate.lnBd")

# Deal with relative units that were not converted previously (ifelse(trait_unit == "relative", trait_value, (quantile_90 - control) / control))
control_range <- as.data.frame(traits_exposed %>%
                                 dplyr::filter(Bd == 0) %>%
                                 dplyr::group_by(specific_response) %>%
                                 dplyr::summarise(control_90 = max(control_90),
                                                  control_10 = min(control_10)))

# Generate Bd load sequence from min and max and predict response
xseq       <- seq(0, max(traits_exposed$lnBd[!is.na(traits_exposed$lnBd)]), length.out = 50)
eq_bd_pred <- data.frame(lm_eq %>%
                           dplyr::group_by(specific_response) %>%
                           dplyr::summarise(pred = intercept + slope * xseq) %>%
                           dplyr::mutate(pred_bd    = xseq,
                                         trait      = traits_exposed$trait[match(specific_response, traits_exposed$specific_response)],
                                         control    = traits_exposed$control[match(specific_response, traits_exposed$specific_response)],
                                         control_90 = control_range$control_90[match(specific_response, control_range$specific_response)],
                                         control_10 = control_range$control_10[match(specific_response, control_range$specific_response)]))

# Keep specific responses with R2 >= 0.1 and p value <= 0.05
responses_kept <- R_sq_model %>% dplyr::filter(adj.r.squared >= 0.1 | p.value <= 0.05) %>% droplevels()
kept_pred      <- eq_bd_pred %>% dplyr::filter(specific_response %in% c(levels(responses_kept$specific_response))) %>% droplevels()

#kept_pred %>% filter(specific_response == "Litoria verreauxii alpina_Serotonin")

# Filter Bd load lower than 10th and higher than 90th quantile control values 
sensitivity_pos <- data.frame(kept_pred %>%
                                dplyr::group_by(specific_response) %>%
                                dplyr::filter(pred_bd != 0 & pred > 0) %>% # in positive direction
                                dplyr::filter(pred > control_90) %>%
                                dplyr::summarize(minBd = min(pred_bd)))

sensitivity_neg <- data.frame(kept_pred %>%
                                dplyr::group_by(specific_response) %>%
                                dplyr::filter(pred_bd != 0 & pred < 0) %>%  # in negative direction
                                dplyr::filter(if (control < 0) pred < control_90 else pred < control_10) %>% # negative controls are flipped
                                dplyr::summarize(minBd = min(pred_bd)))

# Combine sensitivity_pos and sensitivity_neg
sensitivity <- rbind(sensitivity_pos, sensitivity_neg) %>%
  dplyr::mutate(lnminBd     = minBd,
                response    = traits_exposed$response[match(specific_response, traits_exposed$specific_response)],
                trait       = traits_exposed$trait[match(specific_response, traits_exposed$specific_response)],
                bio_hier    = traits_exposed$bio_hier[match(specific_response, traits_exposed$specific_response)],
                resistance  = traits_exposed$resistance[match(specific_response, traits_exposed$specific_response)],
                species     = traits_exposed$species[match(specific_response, traits_exposed$specific_response)],
                species_OTL = traits_exposed$species_OTL[match(specific_response, traits_exposed$specific_response)],
                life_stage  = traits_exposed$life_stage[match(specific_response, traits_exposed$specific_response)],
                inv_temp    = traits_exposed$inv_temp[match(specific_response, traits_exposed$specific_response)])


# Due to low number of traits in resilient and tolerant category with strong correlation, only sensitive species were included.
sensitivity$species <- sub(" ", "_", sensitivity$species_OTL)
sensitivity$species <- factor(sensitivity$species, levels = tree_tip_label) # relevel order by tree

sensitive_model <- brms::brm(lnminBd ~ -1 + trait + inv_temp + (1 | species_OTL) + (1 | gr(species, cov = phylo)),
                             data    = sensitivity %>%
                               dplyr::group_by(trait) %>% 
                               dplyr::filter(resistance == "Sensitive"),
                             family  = gaussian,
                             data2   = list(phylo = phylo_cor),
                             prior   = c(prior(normal(0, 5), lb = 0), # mean of 0 and SD of 5
                                       prior(student_t(3, 0, 5), "sd"), # class of random effect deviation
                                       prior(student_t(3, 0, 5), "sigma")), # residual SD parameter
                             iter    = 1e4, warmup = 5e3, cores = 4, chains = 4,
                             control = list(adapt_delta = 0.99, max_treedepth = 18))

summary(sensitive_model)

# Extract marginal effects
sensitive_me <- brms::conditional_effects(sensitive_model)
sensitive_me <- as.data.frame(sensitive_me[[1]]) %>% 
  dplyr::rename(estimate = estimate__, ci.lb = lower__, ci.ub = upper__) %>%
  dplyr::mutate(estimate = exp(estimate),
                ci.lb    = exp(ci.lb),
                ci.ub    = exp(ci.ub))

sensitive_post <- brms::posterior_samples(sensitive_model) %>%
  dplyr::sample_n(4000) %>%
  dplyr::select(b_traitBehaviour:b_traitSkinintegrity) %>%
  reshape2::melt(value.name = "sample")
  
sensitive_post$variable <- substring(sensitive_post$variable, 8)
source("https://raw.githubusercontent.com/datavizpyr/data/master/half_flat_violinplot.R")

sensitive_post %>%
  dplyr::group_by(variable) %>%
  dplyr::summarise(mean = mean(sample)) %>%
  dplyr::mutate(variable = forcats::fct_reorder(variable, mean)) %>%
  ggplot() +
  geom_point(aes(x = variable, y = mean), size = 2) +
  geom_flat_violin(data = sensitive_post, aes(x = variable, y = sample, fill = variable), colour = NA) +
  geom_point(aes(x = variable, y = mean), size = 3) +
  geom_point(aes(x = variable, y = mean, colour = variable), size = 2) +
  colorspace::scale_fill_discrete_sequential(palette = "RedOr", nmax = 12, order = 4:12) +
  colorspace::scale_colour_discrete_sequential(palette = "RedOr", nmax = 12, order = 4:12) +
  xlab(NULL) +
  ylab("Infection intensity (log ZE + 1)") +
  coord_flip() +
  mytheme()
  

## SURVIVAL ## ----------------------------------------------------------------------------------
as.data.frame(surv_data_clean %>% 
                dplyr::summarise(study_n = length(unique(study_ID)),
                                 species_n = length(unique(species))))

# Rename species based on species_OTL
surv_data_clean$species <- sub(" ", "_", surv_data_clean$species_OTL)
surv_data_clean$species <- factor(surv_data_clean$species, levels = tree_tip_label) # relevel order by tree

# Binary logistic regression model
surv_data_clean <- surv_data_clean %>% dplyr::mutate(species = factor(species))
surv_model <- brms::brm(alive ~ 1 + lnBd + time_post_days + resistance + life_stage + inv_temp + species_OTL + (1 | study_ID:species_OTL) + (1 | gr(species, cov = phylo)),  
                        data   = surv_data_clean,
                        family = bernoulli(link = "logit"), 
                        data2   = list(phylo = phylo_cor),
                        prior = c(prior(normal(0, 2), class = Intercept),
                                  prior(normal(0, 2), class = b)),
                        iter = 5e3, warmup = 2e3, chains = 4, cores = 4,
                        control = list(adapt_delta = 0.99, max_treedepth = 18))

summary(surv_model)
#plot(conditional_effects(surv_model), points = TRUE)

# Extract marginal effects
surv_marg_eff <- as.data.frame(ggeffects::ggpredict(surv_model, terms = c("lnBd[sample = 300]", "species_OTL")))
surv_marg_eff <- surv_marg_eff %>%
  dplyr::mutate(resistance = surv_data_clean$resistance[match(group, surv_data_clean$species_OTL)],
                species    = surv_data_clean$species[match(group, surv_data_clean$species_OTL)],
                family     = surv_data_clean$Family[match(group, surv_data_clean$species_OTL)],
                life_stage = surv_data_clean$life_stage[match(group, surv_data_clean$species_OTL)])
  
# Plot
surv_plot <- surv_marg_eff %>%
  ggplot(aes(x = x, y = predicted, colour = resistance)) +
  geom_vline(xintercept = log(1e4), linetype = "dashed", alpha = 0.5) +
  geom_line(aes(group = group, linetype = life_stage), size = 0.5) +
  geom_point(data = surv_data_clean, aes(x = lnBd, y = alive, colour = resistance), size = 2, alpha = 0.5) +
  scale_colour_manual(values = c("#FB8D46", "#CD5A7E", "#159291")) +
  ylab("Survival probabilty") +
  xlab("Infection intensity (ln ZE + 1)") +
  facet_wrap(vars(resistance), ncol = 1) +
  mytheme() +
  theme(legend.position = "top")

# Predict 50% infliction point
inflict_est <- as.data.frame(surv_marg_eff %>%
                                dplyr::group_by(group, resistance, life_stage) %>% 
                                dplyr::filter(predicted < 0.53 & predicted > 0.47) %>% 
                                dplyr::summarise(threshold = mean(x),
                                                 resistance = unique(resistance),
                                                 life_stage = unique(life_stage),
                                                 family     = unique(family)) %>%
                               dplyr::mutate(bd = exp(threshold)))

inflict_90_est <- as.data.frame(surv_marg_eff %>%
                               dplyr::group_by(group, resistance, life_stage) %>% 
                               dplyr::filter(predicted < 0.12 & predicted > 0.08) %>% 
                               dplyr::summarise(threshold = mean(x),
                                                resistance = unique(resistance),
                                                life_stage = unique(life_stage),
                                                family     = unique(family)) %>%
                               dplyr::mutate(bd = exp(threshold)))

inflict_plot <- inflict_est %>%
  dplyr::mutate(group = forcats::fct_reorder(group, -threshold)) %>%
  ggplot(aes(x = group, y = threshold, colour = resistance)) +
  geom_point(aes(shape = life_stage), size = 3) +
  scale_colour_manual(values = c("#FB8D46", "#CD5A7E", "#159291")) +
  ylab(expression(atop(italic("Bd")~"load (ln ZE + 1)", "with 50% mortality probabilty"))) +
  xlab(NULL) +
  scale_x_discrete(position = "top") +
  coord_flip() +
  mytheme() + 
  theme(legend.position = "top",
        axis.text.y = element_text(size = rel(0.8)))

cowplot::plot_grid(surv_plot, inflict_plot, 
                   labels = c('a', 'b'), 
                   ncol = 2, 
                   align = "h", 
                   axis = "bt",
                   rel_widths = c(0.9, 1))

## BD & TIME-DEPENDENT MORTALITY ## -----------------------------------------------------------
# (1 | study_ID:species) = different studies might cause different effects on a given species. 
time_bd_model <- brms::brm(lnBd ~ lnTime * life_stage + resistance + inv_temp + (1 | study_ID:species_OTL),  
                           data    = surv_data_clean %>% dplyr::filter(alive == 0),
                           family  = gaussian(),
                           prior   = c(prior(normal(0, 2), class = Intercept),
                                       prior(normal(0, 2), class = b)),
                           iter    = 5e3, warmup = 2e3, chains = 4, cores = 4,
                           control = list(adapt_delta = 0.99, max_treedepth = 18))

summary(time_bd_model)
#plot(conditional_effects(time_bd_model), points = TRUE)
parnames(time_bd_model)

# Life stage only as resistance showed no difference
time_marg_eff <- as.data.frame(ggeffects::ggpredict(time_bd_model, terms = c("lnTime[sample = 50]", "life_stage")))

# Create matrix of min and max values per group
time_range <- surv_data_clean %>%
  dplyr::filter(alive == 0) %>%
  dplyr::group_by(life_stage) %>%
  dplyr::summarise(min = min(lnTime[!is.na(lnTime)]),
                   max = max(lnTime[!is.na(lnTime)])) %>%
  as.data.frame()

# Add min and max values to model df and keep predictions within data range
time_marg_eff$min <- time_range$min[match(time_marg_eff$group, time_range$life_stage)]
time_marg_eff$max <- time_range$max[match(time_marg_eff$group, time_range$life_stage)]
time_marg_eff     <- time_marg_eff %>% 
  dplyr::group_by(group) %>%
  dplyr::filter(x >= min & x < max) %>%
  dplyr::mutate(life_stage = group,
                days       = exp(x),
                Bd         = exp(predicted))

# Plot 
time_bd_plot <- ggplot(data = time_marg_eff, aes(x = x, y = predicted, group = life_stage)) +
  geom_point(data = surv_data_clean %>% dplyr::filter(alive == 0), aes(x = lnTime, y = lnBd, colour = life_stage, shape = life_stage), size = 2, alpha = 0.4) +
  geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high, fill = life_stage), alpha = 0.1) +
  geom_line(aes(colour = life_stage, linetype = life_stage), size = 1) +
  colorspace::scale_colour_discrete_sequential(palette = "Blues 3", nmax = 5, order = c(3,5)) +
  colorspace::scale_fill_discrete_sequential(palette = "Blues 3", nmax = 5, order = c(3,5)) +
  ylab("Infection intensity (ln ZE + 1)") +
  xlab("Time post exposure (ln days)") +
  mytheme() +
  theme(legend.position = "top")

time_plot <- surv_data_clean %>% dplyr::filter(alive == 0) %>% 
  ggplot(aes(x = life_stage, y = lnTime)) +
  geom_flat_violin(aes(fill = life_stage), colour = NA, alpha = 0.5) +
  geom_point(aes(shape = life_stage, colour = life_stage), size = 1.5) +
  colorspace::scale_colour_discrete_sequential(palette = "Blues 3", nmax = 5, order = c(3,5)) +
  colorspace::scale_fill_discrete_sequential(palette = "Blues 3", nmax = 5, order = c(3,5)) +
  ylab("Time post exposure (ln days)") + xlab(NULL) +
  mytheme() +
  theme(legend.position = "top")

cowplot::plot_grid(time_bd_plot, time_plot, 
                   labels = c('a', 'b'), 
                   ncol = 2, 
                   align = "h", 
                   axis = "bt",
                   rel_widths = c(1, 0.5))

## SUPPLEMENTARY FIGURE ## ---------------------------------------------------------------------
effect_size_full %>%
  ggplot(aes(x = response, y = Zr, colour = trait)) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_point(aes(shape = resistance), show.legend = FALSE) +
  viridis::scale_colour_viridis(discrete = TRUE) +
  xlab(NULL) + ylab(expression("Effect size"~(italic(Z)[r]))) +
  coord_flip() +
  facet_wrap(~ trait, ncol = 3, scales = "free_y") +
  mytheme() + theme(axis.text = element_text(size = 5),
                    legend.position = "bottom")

unique(effect_size_full$trait)
effect_size_full %>%
  filter(trait == "Skin integrity") %>%
  summarise(author = unique(author_first))

# Fig 
bg <- surv_data_clean %>% dplyr::filter(alive == 0) %>% dplyr::select(-resistance)

ggplot() +
  geom_point(data = bg, aes(x = lnTime, y = lnBd, shape = life_stage), colour = "#ebebeb", size = 2) +
  geom_point(data = surv_data_clean %>% dplyr::filter(alive == 0), aes(x = lnTime, y = lnBd, colour = life_stage, shape = life_stage), size = 2) +
  colorspace::scale_colour_discrete_sequential(palette = "Blues 3", nmax = 5, order = c(3,5)) +
  ylab("Infection intensity (ln ZE + 1)") +
  xlab("Time post exposure (ln days)") +
  facet_wrap(~ resistance, ncol = 1) +
  mytheme() +
  theme(legend.position = "top")


## METAFOR TEST ##--------------------------------------------------------------------------------------------
library(metafor)

# use Zr-v, not Zr_sei as V (https://www.metafor-project.org/doku.php/tips:input_to_rma_function?fbclid=IwAR1C8q2n0qp8TEHmMWWY1BsY7qIqfXjxi9loiNslPLXwAMzvPLOeKuUalD0)
overall_model_2 <- metafor::rma.mv(yi = Zr, V = Zr_v, 
                                   mod    = ~ 1 + time_post_days + life_stage + origin + resistance + Zr_sei + year_centre, 
                                   random = list(~1 | species, ~1 | study_ID, ~1 | es_ID, ~1 | species_OTL),
                                   R      = list(species = phylo_cor),
                                   method = "REML", test = "t",
                                   data   = effect_size_full)
summary(overall_model_2)

trait_model_2 <- metafor::rma.mv(yi = Zr, V = Zr_v, 
                                 mod    = ~ -1 + trait + resistance + Zr_sei + year_centre, 
                                 random = list(~1 | species, ~1 | study_ID, ~1 | es_ID, ~1 | species_OTL),
                                 R      = list(species = phylo_cor),
                                 method = "ML", test = "t",
                                 data   = effect_size_full %>%
                                   dplyr::group_by(trait) %>% 
                                   dplyr::filter(n() >= 5))
summary(trait_model_2)

# R2
fix <- var(as.numeric(as.vector(trait_model_2$b) %*% t(as.matrix(trait_model_2$X))))
R2m <- fix / (fix + sum(trait_model_2$sigma2))
100 * R2m

trait_est <- data.frame(
  trait      = substr(row.names(trait_model_2$b), 6, 100), #remove the word trait
  estimate   = trait_model_2$b, 
  ci.lb      = trait_model_2$ci.lb, 
  ci.ub      = trait_model_2$ci.ub)
trait_est <- head(trait_est, -4)

trait_est %>% 
  #dplyr::mutate(trait = forcats::fct_reorder(trait, estimate)) %>%
  ggplot(aes(x = trait, y = estimate, group = trait)) +
  geom_vline(aes(xintercept = trait), linetype = "dotted", colour = "grey30", alpha = 0.4) +
  geom_hline(yintercept = 0, linetype = "dashed") + 
  geom_point(data = effect_size_full %>% dplyr::filter(resistance != "Sensitive"), aes(x = trait, y = Zr, size = Zr_inv), position = position_nudge(x = -0.05), colour = "#cfcfcf") +
  geom_point(data = effect_size_full %>% dplyr::filter(resistance == "Sensitive"), aes(x = trait, y = Zr, size = Zr_inv), position = position_nudge(x = 0.05), colour = "#FB8D46", alpha = 0.2) +
  geom_point(aes(colour = trait), size = 4, show.legend = FALSE) +
  geom_errorbar(aes(ymin = ci.lb, ymax = ci.ub, colour = trait), size = 0.8, width = 0.1, show.legend = FALSE) +
  colorspace::scale_colour_discrete_sequential(palette = "RedOr") +
  xlab(NULL) + ylab(expression("Effect size"~(italic(Z)[r]))) +
  coord_flip() +
  mytheme() + 
  theme(legend.position = "bottom")

trait_est %>% 
  ggplot(aes(x = trait, y = estimate, group = trait)) +
  geom_hline(yintercept = 0, linetype = "dashed") + 
  ggforce::geom_sina(data = effect_size_full, aes(x = trait, y = Zr, size = Zr_inv), colour = "#cfcfcf") + 
  geom_point(aes(colour = trait), size = 4, show.legend = FALSE) +
  geom_errorbar(aes(ymin = ci.lb, ymax = ci.ub, colour = trait), size = 0.8, width = 0.1, show.legend = FALSE) +
  colorspace::scale_colour_discrete_sequential(palette = "RedOr") +
  xlab(NULL) + ylab(expression("Effect size"~(italic(Z)[r]))) +
  coord_flip() +
  mytheme() + 
  theme(legend.position = "bottom")


## EXTRA CODE ## --------------------------------------------------------
# to delete
effect_size_full %>% 
  ggplot(aes(x = trait, y = Zr)) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_jitter(aes(shape = life_stage, colour = trait), size = 2, position = position_dodge(0.5)) +
  colorspace::scale_colour_discrete_sequential(palette = "RedOr") +
  xlab(NULL) + ylab(expression("Effect size"~(italic(Z)[r]))) +
  facet_grid(. ~ resistance) +
  coord_flip() +
  mytheme() + theme(legend.position = "bottom")


