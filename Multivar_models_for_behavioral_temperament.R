################################################################################
# Microbial exposure at birth and the development of behavioral 
# temperament during the first three years of childhood
# 
# This R-script fits the cross-sectional and longitudinal models
# investigating associations between neonate microbial exposure and 
# the emerging behavioral temperament measures at the ages of 1, 2, and 3 years
################################################################################


# libraries
library(haven)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(lubridate)
library(tidyr)
library(lme4)
library(emmeans)
library(margins)
library(gridExtra)
library(haven)
library(broom)
library(lmerTest)
library(mice)


###### DATA
dat <- read_sav("~/Data/full_dat.sav")


# Missing values, mean imputation for covariates
# Missing values in outcome (behavioral temperament, will be handled either with multiple imputation or in a mixed model setting)
dat$mat.age.imp <- ifelse(is.na(dat$mat_age), mean(dat$mat_age, na.rm = TRUE), dat$mat_age)
dat$epds.tr3.imp <- ifelse(is.na(dat$epds.tr3), mean(dat$epds.tr3, na.rm = TRUE), dat$epds.tr3)

# prepare all data for all temeperament subscales
# - long formats
# - year 1 timepoint only

### EFF
eff.df <- dat %>%
          select(id, gender, mat.age.imp, epds.tr3.imp, extr.premat, mat_age, epds.tr3, 
                 alpha.pc1:beta.pco2, Observed, Lactobacillus:Anaerococcus, starts_with("eff"))

eff.long <- pivot_longer(eff.df, cols = eff_1y:eff_3y, values_to = "scale.score", names_to = "year")
eff.long <- eff.long %>%
            group_by(id) %>%
            mutate(year.num = row_number()-1) %>%
            ungroup()
# year 1 only
eff1 <- eff.long %>%
  filter(year.num == 0)

# AB use during birth, variable from a separate file
ab.dat <- read_sav("~/Data/NID5.sav")

ab.dat <- ab.dat %>%
        select(J‰rjestelm‰ID, birth_AB) %>%
        rename(id = J‰rjestelm‰ID)
ab.eff.df <- left_join(eff.df, ab.dat, by = "id")

# EFF long formats with and without AB use
ab.yes.eff.df <- ab.eff.df %>%
                filter(birth_AB == 1)
ab.no.eff.df <- ab.eff.df %>%
                filter(birth_AB == 0)

ab.yes.eff.long <- pivot_longer(ab.yes.eff.df, cols = eff_1y:eff_3y, values_to = "scale.score", names_to = "year")
ab.yes.eff.long <- ab.yes.eff.long %>%
                    group_by(id) %>%
                    mutate(year.num = row_number()-1) %>%
                    ungroup()
# year 1 only
ab.yes.eff1 <- ab.yes.eff.long %>%
                filter(year.num == 0)

# no AB long
ab.no.eff.long <- pivot_longer(ab.no.eff.df, cols = eff_1y:eff_3y, values_to = "scale.score", names_to = "year")
ab.no.eff.long <- ab.no.eff.long %>%
                  group_by(id) %>%
                  mutate(year.num = row_number()-1) %>%
                  ungroup()
# year 1 only
ab.no.eff1 <- ab.no.eff.long %>%
              filter(year.num == 0)


### NEG
neg.df <- dat %>%
          select(id, gender, mat.age.imp, epds.tr3.imp, extr.premat, mat_age, epds.tr3, 
               alpha.pc1:beta.pco2, Observed, Lactobacillus:Anaerococcus, starts_with("neg"))

neg.long <- pivot_longer(neg.df, cols = neg_1y:neg_3y, values_to = "scale.score", names_to = "year")
neg.long <- neg.long %>%
            group_by(id) %>%
            mutate(year.num = row_number()-1) %>%
            ungroup()
# year 1 only
neg1 <- neg.long %>%
        filter(year.num == 0)

# NEG long formats with and without AB use
ab.neg.df <- left_join(neg.df, ab.dat, by = "id")

ab.yes.neg.df <- ab.neg.df %>%
                filter(birth_AB == 1)
ab.no.neg.df <- ab.neg.df %>%
                filter(birth_AB == 0)

ab.yes.neg.long <- pivot_longer(ab.yes.neg.df, cols = neg_1y:neg_3y, values_to = "scale.score", names_to = "year")
ab.yes.neg.long <- ab.yes.neg.long %>%
                  group_by(id) %>%
                  mutate(year.num = row_number()-1) %>%
                  ungroup()
# year 1 only
ab.yes.neg1 <- ab.yes.neg.long %>%
                filter(year.num == 0)

# no AB long
ab.no.neg.long <- pivot_longer(ab.no.neg.df, cols = neg_1y:neg_3y, values_to = "scale.score", names_to = "year")
ab.no.neg.long <- ab.no.neg.long %>%
                  group_by(id) %>%
                  mutate(year.num = row_number()-1) %>%
                  ungroup()
# year 1 only
ab.no.neg1 <- ab.no.neg.long %>%
              filter(year.num == 0)

### SUR
sur.df <- dat %>%
          select(id, gender, mat.age.imp, epds.tr3.imp, extr.premat, mat_age, epds.tr3, 
                 alpha.pc1:beta.pco2, Observed, Lactobacillus:Anaerococcus, starts_with("sur"))

sur.long <- pivot_longer(sur.df, cols = sur_1y:sur_3y, values_to = "scale.score", names_to = "year")
sur.long <- sur.long %>%
            group_by(id) %>%
            mutate(year.num = row_number()-1) %>%
            ungroup()
# year 1 only
sur1 <- sur.long %>%
        filter(year.num == 0)

# SUR long formats with and without AB use
ab.sur.df <- left_join(sur.df, ab.dat, by = "id")

ab.yes.sur.df <- ab.sur.df %>%
                filter(birth_AB == 1)
ab.no.sur.df <- ab.sur.df %>%
                filter(birth_AB == 0)

ab.yes.sur.long <- pivot_longer(ab.yes.sur.df, cols = sur_1y:sur_3y, values_to = "scale.score", names_to = "year")
ab.yes.sur.long <- ab.yes.sur.long %>%
                  group_by(id) %>%
                  mutate(year.num = row_number()-1) %>%
                  ungroup()
# year 1 only
ab.yes.sur1 <- ab.yes.sur.long %>%
                filter(year.num == 0)

# no AB long
ab.no.sur.long <- pivot_longer(ab.no.sur.df, cols = sur_1y:sur_3y, values_to = "scale.score", names_to = "year")
ab.no.sur.long <- ab.no.sur.long %>%
                  group_by(id) %>%
                  mutate(year.num = row_number()-1) %>%
                  ungroup()
# year 1 only
ab.no.sur1 <- ab.no.sur.long %>%
              filter(year.num == 0)


# plot Figure 2

plot.all <- rbind(eff.long, neg.long, sur.long) %>%
            mutate(scale = as.factor(substr(year, 1, 3)))

p <- ggplot(plot.all, aes(x=year.f, y=scale.score, fill = gender)) + 
  geom_boxplot() + 
  scale_fill_manual(values=c("#6699ff", "#ffcc99")) +
  theme( 
    strip.text.x = element_text(size = 12, face = "bold"), 
    panel.border=element_blank(),
    panel.spacing.x = unit(0,"line")) +
  xlab("Year") +
  ylab("Scale score") +
  scale_y_continuous(breaks=seq(1, 7, 1), limits=c(1, 7)) +
  labs(fill = " ") +
  facet_wrap(~scale, ncol = 3, labeller = labeller(scale = 
                                                     c("eff" = "Effortful control",
                                                       "neg" = "Negative affect",
                                                       "sur" = "Surgency")))
p
ggsave("~/Figures/Scale_scores_by_sex.png", p, width = 3000, height =1500, units = "px")


####### MODELS


#### Subscale EFF

### EFF: Year 1 only vs.  alpha and beta diversity score PC(O)s, multiple imputed year 1 values (n=335)


# multiple imputation
set.seed(123)

imp.eff <- mice(ab.eff.df, maxit=0)
pred.m <- imp.eff$predictorMatrix

# Setting values of variables we'd like to leave out to 0 in the predictor matrix
pred.m[, c("id")] <- 0
pred.m[, c("extr.premat")] <- 0
pred.m[, c("mat_age")] <- 0
pred.m[, c("epds.tr3")] <- 0
pred.m[, c("alpha.pc1")] <- 0
pred.m[, c("alpha.pc2")] <- 0
pred.m[, c("beta.pco1")] <- 0
pred.m[, c("beta.pco2")] <- 0
pred.m[, c("Observed")] <- 0
pred.m[, c("Lactobacillus")] <- 0
pred.m[, c("Streptococcus")] <- 0
pred.m[, c("Staphylococcus")] <- 0
pred.m[, c("Corynebacterium")] <- 0
pred.m[, c("Stenotrophomonas")] <- 0
pred.m[, c("Bacteroides")] <- 0
pred.m[, c("Prevotella")] <- 0
pred.m[, c("Anaerococcus")] <- 0

imp.eff2 <- mice(eff.imp.data, maxit = 5, predictorMatrix = pred.m, print =  FALSE)

# Regression with MI data
fitimp.eff <- with(imp.eff2,
                   lm(eff_1y ~ gender + mat.age.imp + epds.tr3.imp +  
                        alpha.pc1 + alpha.pc2 +
                        beta.pco1 + beta.pco2))
sum <- summary(pool(fitimp.eff)) %>% mutate_if(is.numeric, ~round(., 3))
sum


### EFF: Year 1 only vs. microbiome alpha and beta diversity score PC(O)s, no ABs used during birth, multiple imputed year 1 values (n=198)

# Analysis separately in AB and non-AB groups -->
# transform imputed data into long format, remove ids with AB use/ non-AB use

# First, turn the datasets into long format
imp_long <- mice::complete(imp.eff2, action="long", include = TRUE)

# combine with AB information
imp_long_ab <- left_join(imp_long, ab.dat, by = "id")
imp_long_ab0 <- imp_long_ab %>% filter(birth_AB == 0)
imp_long_ab1 <- imp_long_ab %>% filter(birth_AB == 1)

# Convert back to mids type
imp_long_ab0_mids<-as.mids(imp_long_ab0)
imp_long_ab1_mids<-as.mids(imp_long_ab1)


# regression for MI y1, no AB
fitimp.eff.noAB <- with(imp_long_ab0_mids,
                        lm(eff_1y ~ gender + mat.age.imp + epds.tr3.imp +  
                             alpha.pc1 + alpha.pc2 +
                             beta.pco1 + beta.pco2))

sum <- summary(pool(fitimp.eff.noAB)) %>% mutate_if(is.numeric, ~round(., 3))
sum


### EFF: Year 1 only vs. microbiome alpha and beta diversity score PC(O)s, ABs used during birth, multiple imputed year 1 values (n=137)

# regression for MI y1, no AB
fitimp.eff.yesAB <- with(imp_long_ab1_mids,
                         lm(eff_1y ~ gender + mat.age.imp + epds.tr3.imp +  
                              alpha.pc1 + alpha.pc2 +
                              beta.pco1 + beta.pco2))

sum <- summary(pool(fitimp.eff.yesAB)) %>% mutate_if(is.numeric, ~round(., 3))
sum



### EFF: Year 1 only vs. microbiome alpha and beta diversity score PC(O)s, complete case analysis (n=253)

model.eff1.pco <- lm(scale.score ~ gender + mat.age.imp + epds.tr3.imp +  
                       alpha.pc1 + alpha.pc2 +
                       beta.pco1 + beta.pco2, data = eff1)
round(summary(model.eff1.pco)$coefficients, 3)


### EFF: Year 1 only vs. top genera, multiple imputed year 1 values (n=335)

# Regression 
fitimp.eff.genera <- with(imp.eff2,
                          lm(eff_1y ~ gender + mat.age.imp + epds.tr3.imp + Lactobacillus + 
                               Streptococcus + Staphylococcus + Corynebacterium + Stenotrophomonas + 
                               Bacteroides + Prevotella + Anaerococcus))

sum <- summary(pool(fitimp.eff.genera)) %>% mutate_if(is.numeric, ~round(., 3))
sum


### EFF: Year 1 only vs. top genera, no ABs used during birth, multiple imputed year 1 values (n=198)

# regression for MI y1, no AB
fitimp.eff.genera.noAB <- with(imp_long_ab0_mids,
                               lm(eff_1y ~ gender + mat.age.imp + epds.tr3.imp + Lactobacillus + 
                                    Streptococcus + Staphylococcus + Corynebacterium + Stenotrophomonas + 
                                    Bacteroides + Prevotella + Anaerococcus))

sum <- summary(pool(fitimp.eff.genera.noAB)) %>% mutate_if(is.numeric, ~round(., 3))
sum

### EFF: Year 1 only vs. genera, ABs used during birth, multiple imputed year 1 values (n=137)


# regression for MI y1, no AB
fitimp.eff.genera.yesAB <- with(imp_long_ab1_mids,
                                lm(eff_1y ~ gender + mat.age.imp + epds.tr3.imp + Lactobacillus + 
                                     Streptococcus + Staphylococcus + Corynebacterium + Stenotrophomonas + 
                                     Bacteroides + Prevotella + Anaerococcus))

sum <- summary(pool(fitimp.eff.genera.yesAB)) %>% mutate_if(is.numeric, ~round(., 3))
sum


### EFF: Year 1 only vs. top genera, complete case analysis (n=253)


model.eff1.genera <- lm(scale.score ~ gender + mat.age.imp + epds.tr3.imp +  
                          Lactobacillus + Streptococcus + Staphylococcus + Corynebacterium + 
                          Stenotrophomonas + Bacteroides + Prevotella + Anaerococcus, data = eff1)
round(summary(model.eff1.genera)$coefficients, 3)


## Microbiome vs. repeatedly measured EFF over years 1, 2, 3
## Missing values handled within the mixed model framework

### EFF: Repeatedly measured outcome at years 1, 2, and 3 vs. microbiome alpha and beta diversity score PC(O)s (n=335)

model.eff123.pco <- lmer(scale.score ~ gender + mat.age.imp + epds.tr3.imp + year.num + 
                           alpha.pc1 + alpha.pc2 +
                           beta.pco1 + beta.pco2 + (1|id),
                         data = eff.long, REML=TRUE)
sum <- coef(summary(model.eff123.pco))
round(sum, 3)


### EFF: Repeatedly measured outcome at years 1, 2, and 3 vs. microbiome alpha and beta diversity score PC(O)s, no ABs during birth (n=198)

model.eff123.pco.noAB <- lmer(scale.score ~ gender + mat.age.imp + epds.tr3.imp + year.num + 
                                alpha.pc1 + alpha.pc2 +
                                beta.pco1 + beta.pco2 + (1|id),
                              data = ab.no.eff.long, REML=TRUE)
sum <- coef(summary(model.eff123.pco.noAB))
round(sum, 3)


### EFF: Repeatedly measured outcome at years 1, 2, and 3 vs. microbiome alpha and beta diversity score PC(O)s, ABs used during birth (n=137)

model.eff123.pco.AB <- lmer(scale.score ~ gender + mat.age.imp + epds.tr3.imp + year.num + 
                              alpha.pc1 + alpha.pc2 +
                              beta.pco1 + beta.pco2 + (1|id),
                            data = ab.yes.eff.long, REML=TRUE)
sum <- coef(summary(model.eff123.pco.AB))
round(sum, 3)


### EFF: Repeatedly measured outcome at years 1, 2, and 3 vs. top genera (n=335)

model.eff123.genera <- lmer(scale.score ~ gender + mat.age.imp + epds.tr3.imp + year.num + 
                              Lactobacillus + Streptococcus + Staphylococcus + Corynebacterium + 
                              Stenotrophomonas + Bacteroides + Prevotella + Anaerococcus + (1|id),
                            data = eff.long, REML=TRUE)
sum <- coef(summary(model.eff123.genera))
round(sum, 3)

### EFF: Repeatedly measured outcome at years 1, 2, and 3 vs. top genera, no AB during birth (n=198)

model.eff123.genera.noAB <- lmer(scale.score ~ gender + mat.age.imp + epds.tr3.imp + year.num + 
                                   Lactobacillus + Streptococcus + Staphylococcus + Corynebacterium + 
                                   Stenotrophomonas + Bacteroides + Prevotella + Anaerococcus + (1|id),
                                 data = ab.no.eff.long, REML=TRUE)
sum <- coef(summary(model.eff123.genera.noAB))
round(sum, 3)


### EFF: Repeatedly measured outcome at years 1, 2, and 3 vs. top genera, ABs used during birth (n=198)


model.eff123.genera.AB <- lmer(scale.score ~ gender + mat.age.imp + epds.tr3.imp + year.num + 
                                 Lactobacillus + Streptococcus + Staphylococcus + Corynebacterium + 
                                 Stenotrophomonas + Bacteroides + Prevotella + Anaerococcus + (1|id),
                               data = ab.yes.eff.long, REML=TRUE)
sum <- coef(summary(model.eff123.genera.AB))
round(sum, 3)




#### Subscale NEG

### NEG: Year 1 only vs.  microbiome alpha and beta diversity score PC(O)s, multiple imputed year 1 values (n=335)

# Multiple imputation
imp.neg <- mice(ab.neg.df, maxit=0)
pred.m <- imp.neg$predictorMatrix

# Setting values of variables we'd like to leave out to 0 in the predictor matrix
pred.m[, c("id")] <- 0
pred.m[, c("extr.premat")] <- 0
pred.m[, c("mat_age")] <- 0
pred.m[, c("epds.tr3")] <- 0
pred.m[, c("alpha.pc1")] <- 0
pred.m[, c("alpha.pc2")] <- 0
pred.m[, c("beta.pco1")] <- 0
pred.m[, c("beta.pco2")] <- 0
pred.m[, c("Observed")] <- 0
pred.m[, c("Lactobacillus")] <- 0
pred.m[, c("Streptococcus")] <- 0
pred.m[, c("Staphylococcus")] <- 0
pred.m[, c("Corynebacterium")] <- 0
pred.m[, c("Stenotrophomonas")] <- 0
pred.m[, c("Bacteroides")] <- 0
pred.m[, c("Prevotella")] <- 0
pred.m[, c("Anaerococcus")] <- 0

imp.neg2 <- mice(neg.imp.data, maxit = 5, predictorMatrix = pred.m, print =  FALSE)

# Regression 
fitimp.neg <- with(imp.neg2,
                   lm(neg_1y ~ gender + mat.age.imp + epds.tr3.imp +  
                        alpha.pc1 + alpha.pc2 +
                        beta.pco1 + beta.pco2))

fitimp.neg <- with(imp.neg2,
                   lm(neg_1y ~ gender + mat.age.imp + epds.tr3.imp + birth_AB + 
                        alpha.pc1 + alpha.pc2 +
                        beta.pco1 + beta.pco2 +
                        birth_AB:alpha.pc1 + birth_AB:alpha.pc2 +
                        birth_AB:beta.pco1 + birth_AB:beta.pco2))

sum <- summary(pool(fitimp.neg)) %>% mutate_if(is.numeric, ~round(., 3))
sum


### NEG: Year 1 only vs.  microbiome alpha and beta diversity score PC(O)s, no ABs used during birth, multiple imputed year 1 values (n=335)

# Analysis separately in AB and non-AB groups -->
# transform imputed data into long format, remove ids with AB use/ non-AB use

# First, turn the datasets into long format
imp_long <- mice::complete(imp.neg2, action="long", include = TRUE)

# combine with AB information
imp_long_ab <- left_join(imp_long, ab.dat, by = "id")
imp_long_ab0 <- imp_long_ab %>% filter(birth_AB == 0)
imp_long_ab1 <- imp_long_ab %>% filter(birth_AB == 1)

# Convert back to mids type
imp_long_ab0_mids<-as.mids(imp_long_ab0)
imp_long_ab1_mids<-as.mids(imp_long_ab1)

# regression for MI y1, no AB
fitimp.neg.noAB <- with(imp_long_ab0_mids,
                        lm(neg_1y ~ gender + mat.age.imp + epds.tr3.imp +  
                             alpha.pc1 + alpha.pc2 +
                             beta.pco1 + beta.pco2))
sum <- summary(pool(fitimp.neg.noAB)) %>% mutate_if(is.numeric, ~round(., 3))
sum



### NEG: Year 1 only vs.  microbiome alpha and beta diversity score PC(O)s, ABs used during birth, multiple imputed year 1 values (n=335)


# regression for MI y1, no AB
fitimp.neg.yesAB <- with(imp_long_ab1_mids,
                         lm(neg_1y ~ gender + mat.age.imp + epds.tr3.imp +  
                              alpha.pc1 + alpha.pc2 +
                              beta.pco1 + beta.pco2))
sum <- summary(pool(fitimp.neg.yesAB)) %>% mutate_if(is.numeric, ~round(., 3))
sum



### NEG: Year 1 only vs. microbiome alpha and beta diversity score PC(O)s, complete case analysis (n=253)

model.neg1.pco <- lm(scale.score ~ gender + mat.age.imp + epds.tr3.imp +  
                       alpha.pc1 + alpha.pc2 +
                       beta.pco1 + beta.pco2, data = neg1)
round(summary(model.neg1.pco)$coefficients, 3)


## NEG: Year 1 only vs.  top genera, multiple imputed year 1 values (n=335)

fitimp.neg.genera <- with(imp.neg2,
                          lm(neg_1y ~ gender + mat.age.imp + epds.tr3.imp + Lactobacillus + Streptococcus + Staphylococcus + Corynebacterium + 
                               Stenotrophomonas + Bacteroides + Prevotella + Anaerococcus))
sum <- summary(pool(fitimp.neg.genera)) %>% mutate_if(is.numeric, ~round(., 3))
sum


## NEG: Year 1 only vs.  top genera, no ABs used during birth, multiple imputed year 1 values (n=198)

fitimp.neg.genera.noAB <- with(imp_long_ab0_mids,
                               lm(neg_1y ~ gender + mat.age.imp + epds.tr3.imp + Lactobacillus + Streptococcus + Staphylococcus + Corynebacterium + 
                                    Stenotrophomonas + Bacteroides + Prevotella + Anaerococcus))
sum <- summary(pool(fitimp.neg.genera.noAB)) %>% mutate_if(is.numeric, ~round(., 3))
sum

## NEG: Year 1 only vs.  top genera, ABs used during birth, multiple imputed year 1 values (n=137)

fitimp.neg.genera.yesAB <- with(imp_long_ab1_mids,
                                lm(neg_1y ~ gender + mat.age.imp + epds.tr3.imp + Lactobacillus + Streptococcus + Staphylococcus + Corynebacterium + 
                                     Stenotrophomonas + Bacteroides + Prevotella + Anaerococcus))
sum <- summary(pool(fitimp.neg.genera.yesAB)) %>% mutate_if(is.numeric, ~round(., 3))
sum



## NEG: Year 1 only vs. top genera, complete case analysis (n=253)

model.neg1.genera <- lm(scale.score ~ gender + mat.age.imp + epds.tr3.imp +  
                          Lactobacillus + Streptococcus + Staphylococcus + Corynebacterium + 
                          Stenotrophomonas + Bacteroides + Prevotella + Anaerococcus, data = neg1)
round(summary(model.neg1.genera)$coefficients, 3)


## Microbiome vs. repeatedly measured EFF over years 1, 2, 3
## Missing values handled within the mixed model framework


### NEG: Repeatedly measured outcome at years 1, 2, and 3 vs. microbiome alpha and beta diversity score PC(O)s (n=335)

model.neg123.pco <- lmer(scale.score ~ gender + mat.age.imp + epds.tr3.imp + year.num + 
                           alpha.pc1 + alpha.pc2 +
                           beta.pco1 + beta.pco2 + (1|id),
                         data = neg.long, REML=TRUE)
sum <- coef(summary(model.neg123.pco))
round(sum, 3)


### NEG: Repeatedly measured outcome at years 1, 2, and 3 vs. microbiome alpha and beta diversity score PC(O)s,  no ABs used during birth (n=198)

model.neg123.pco.noAB <- lmer(scale.score ~ gender + mat.age.imp + epds.tr3.imp + year.num + 
                                alpha.pc1 + alpha.pc2 +
                                beta.pco1 + beta.pco2 + (1|id),
                              data = ab.no.neg.long, REML=TRUE)
sum <- coef(summary(model.neg123.pco.noAB))
round(sum, 3)


### NEG: Repeatedly measured outcome at years 1, 2, and 3 vs. microbiome alpha and beta diversity score PC(O)s, ABs used during birth (n=137)

model.neg123.pco.AB <- lmer(scale.score ~ gender + mat.age.imp + epds.tr3.imp + year.num + 
                              alpha.pc1 + alpha.pc2 +
                              beta.pco1 + beta.pco2 + (1|id),
                            data = ab.yes.neg.long, REML=TRUE)
sum <- coef(summary(model.neg123.pco.AB))
round(sum, 3)


### NEG: Repeatedly measured outcome at years 1, 2, and 3 vs. top genera (n=335)

model.neg123.genera <- lmer(scale.score ~ gender + mat.age.imp + epds.tr3.imp + year.num + 
                              Lactobacillus + Streptococcus + Staphylococcus + Corynebacterium + 
                              Stenotrophomonas + Bacteroides + Prevotella + Anaerococcus + (1|id),
                            data = neg.long, REML=TRUE)
sum <- coef(summary(model.neg123.genera))
round(sum, 3)


### NEG: Repeatedly measured outcome at years 1, 2, and 3 vs. top genera, no AB during birth (n=198)

model.neg123.genera.noAB <- lmer(scale.score ~ gender + mat.age.imp + epds.tr3.imp + year.num + 
                                   Lactobacillus + Streptococcus + Staphylococcus + Corynebacterium + 
                                   Stenotrophomonas + Bacteroides + Prevotella + Anaerococcus + (1|id),
                                 data = ab.no.neg.long, REML=TRUE)
sum <- coef(summary(model.neg123.genera.noAB))
round(sum, 3)


### NEG: Repeatedly measured outcome at years 1, 2, and 3 vs. top genera, ABs used during birth (n=198)

model.neg123.genera.AB <- lmer(scale.score ~ gender + mat.age.imp + epds.tr3.imp + year.num + 
                                 Lactobacillus + Streptococcus + Staphylococcus + Corynebacterium + 
                                 Stenotrophomonas + Bacteroides + Prevotella + Anaerococcus + (1|id),
                               data = ab.yes.neg.long, REML=TRUE)
sum <- coef(summary(model.neg123.genera.AB))
round(sum, 3)



## Subscale SUR


### SUR: Year 1 only vs. microbiome alpha and beta diversity score PC(O)s, multiple imputed year 1 values (n=335)

imp.sur <- mice(ab.sur.df, maxit=0)
pred.m <- imp.sur$predictorMatrix

# Setting values of variables we'd like to leave out to 0 in the predictor matrix
pred.m[, c("id")] <- 0
pred.m[, c("extr.premat")] <- 0
pred.m[, c("mat_age")] <- 0
pred.m[, c("epds.tr3")] <- 0
pred.m[, c("alpha.pc1")] <- 0
pred.m[, c("alpha.pc2")] <- 0
pred.m[, c("beta.pco1")] <- 0
pred.m[, c("beta.pco2")] <- 0
pred.m[, c("Observed")] <- 0
pred.m[, c("Lactobacillus")] <- 0
pred.m[, c("Streptococcus")] <- 0
pred.m[, c("Staphylococcus")] <- 0
pred.m[, c("Corynebacterium")] <- 0
pred.m[, c("Stenotrophomonas")] <- 0
pred.m[, c("Bacteroides")] <- 0
pred.m[, c("Prevotella")] <- 0
pred.m[, c("Anaerococcus")] <- 0

imp.sur2 <- mice(sur.imp.data, maxit = 5, predictorMatrix = pred.m, print =  FALSE)

# Regression with MI data
fitimp.sur <- with(imp.sur2,
                   lm(sur_1y ~ gender + mat.age.imp + epds.tr3.imp +  
                        alpha.pc1 + alpha.pc2 +
                        beta.pco1 + beta.pco2))

sum <- summary(pool(fitimp.sur)) %>% mutate_if(is.numeric, ~round(., 3))
sum


### SUR: Year 1 only vs. microbiome alpha and beta diversity score PC(O)s, no ABs used during birth, multiple imputed year 1 values (n=198)

# Analysis separately in AB and non-AB groups -->
# transform imputed data into long format, remove ids with AB use/ non-AB use

# First, turn the datasets into long format
imp_long <- mice::complete(imp.sur2, action="long", include = TRUE)

# combine with AB information
imp_long_ab <- left_join(imp_long, ab.dat, by = "id")
imp_long_ab0 <- imp_long_ab %>% filter(birth_AB == 0)
imp_long_ab1 <- imp_long_ab %>% filter(birth_AB == 1)

# Convert back to mids type
imp_long_ab0_mids<-as.mids(imp_long_ab0)
imp_long_ab1_mids<-as.mids(imp_long_ab1)


# regression for MI y1, no AB
fitimp.sur.noAB <- with(imp_long_ab0_mids,
                        lm(sur_1y ~ gender + mat.age.imp + epds.tr3.imp +  
                             alpha.pc1 + alpha.pc2 +
                             beta.pco1 + beta.pco2))

sum <- summary(pool(fitimp.sur.noAB)) %>% mutate_if(is.numeric, ~round(., 3))
sum


### SUR: Year 1 only vs. microbiome alpha and beta diversity score PC(O)s, ABs used during birth, multiple imputed year 1 values (n=137)

# regression for MI y1, no AB
fitimp.sur.yesAB <- with(imp_long_ab1_mids,
                         lm(sur_1y ~ gender + mat.age.imp + epds.tr3.imp +  
                              alpha.pc1 + alpha.pc2 +
                              beta.pco1 + beta.pco2))

sum <- summary(pool(fitimp.sur.yesAB)) %>% mutate_if(is.numeric, ~round(., 3))
sum


### SUR: Year 1 only vs. microbiome alpha and beta diversity score PC(O)s, complete case analysis (n=253)

model.sur1.pco <- lm(scale.score ~ gender + mat.age.imp + epds.tr3.imp +  
                       alpha.pc1 + alpha.pc2 +
                       beta.pco1 + beta.pco2, data = sur1)
round(summary(model.sur1.pco)$coefficients, 3)


### SUR: Year 1 only vs. top genera, multiple imputed year 1 values (n=335)

# Regression 
fitimp.sur.genera <- with(imp.sur2,
                          lm(sur_1y ~ gender + mat.age.imp + epds.tr3.imp + Lactobacillus + Streptococcus + Staphylococcus + Corynebacterium + 
                               Stenotrophomonas + Bacteroides + Prevotella + Anaerococcus))

sum <- summary(pool(fitimp.sur.genera)) %>% mutate_if(is.numeric, ~round(., 3))
sum


### SUR: Year 1 only vs. top genera, no ABs used during birth, multiple imputed year 1 values (n=198)

# regression for MI y1, no AB
fitimp.sur.genera.noAB <- with(imp_long_ab0_mids,
                               lm(sur_1y ~ gender + mat.age.imp + epds.tr3.imp + Lactobacillus + Streptococcus + Staphylococcus + Corynebacterium + 
                                    Stenotrophomonas + Bacteroides + Prevotella + Anaerococcus))

sum <- summary(pool(fitimp.sur.genera.noAB)) %>% mutate_if(is.numeric, ~round(., 3))
sum

### SUR: Year 1 only vs. genera, ABs used during birth, multiple imputed year 1 values (n=137)

# regression for MI y1, no AB
fitimp.sur.genera.yesAB <- with(imp_long_ab1_mids,
                                lm(sur_1y ~ gender + mat.age.imp + epds.tr3.imp + Lactobacillus + Streptococcus + Staphylococcus + Corynebacterium + 
                                     Stenotrophomonas + Bacteroides + Prevotella + Anaerococcus))

sum <- summary(pool(fitimp.sur.genera.yesAB)) %>% mutate_if(is.numeric, ~round(., 3))
sum


### SUR: Year 1 only vs. top genera (n=253)


model.sur1.genera <- lm(scale.score ~ gender + mat.age.imp + epds.tr3.imp +  
                          Lactobacillus + Streptococcus + Staphylococcus + Corynebacterium + 
                          Stenotrophomonas + Bacteroides + Prevotella + Anaerococcus, data = sur1)
round(summary(model.sur1.genera)$coefficients, 3)


### Microbiome vs. repeatedly measured SUR over years 1, 2, 3

### SUR: Repeatedly measured outcome at years 1, 2, and 3 vs. microbiome alpha and beta diversity score PC(O)s (n=335)

model.sur123.pco <- lmer(scale.score ~ gender + mat.age.imp + epds.tr3.imp + year.num + 
                           alpha.pc1 + alpha.pc2 +
                           beta.pco1 + beta.pco2 + (1|id),
                         data = sur.long, REML=TRUE)

sum <- coef(summary(model.sur123.pco))
round(sum, 3)

### SUR: Repeatedly measured outcome at years 1, 2, and 3 vs. microbiome alpha and beta diversity score PC(O)s, no ABs during birth (n=198)

model.sur123.pco.noAB <- lmer(scale.score ~ gender + mat.age.imp + epds.tr3.imp + year.num + 
                                alpha.pc1 + alpha.pc2 +
                                beta.pco1 + beta.pco2 + (1|id),
                              data = ab.no.sur.long, REML=TRUE)

sum <- coef(summary(model.sur123.pco.noAB))
round(sum, 3)


### SUR: Repeatedly measured outcome at years 1, 2, and 3 vs. microbiome alpha and beta diversity score PC(O)s, ABs used during birth (n=137)

model.sur123.pco.AB <- lmer(scale.score ~ gender + mat.age.imp + epds.tr3.imp + year.num + 
                              alpha.pc1 + alpha.pc2 +
                              beta.pco1 + beta.pco2 + (1|id),
                            data = ab.yes.sur.long, REML=TRUE)

sum <- coef(summary(model.sur123.pco.AB))
round(sum, 3)

### SUR: Repeatedly measured outcome at years 1, 2, and 3 vs. top genera (n=335)

model.sur123.genera <- lmer(scale.score ~ gender + mat.age.imp + epds.tr3.imp + year.num + 
                              Lactobacillus + Streptococcus + Staphylococcus + Corynebacterium + 
                              Stenotrophomonas + Bacteroides + Prevotella + Anaerococcus + (1|id),
                            data = sur.long, REML=TRUE)
sum <- coef(summary(model.sur123.genera))
round(sum, 3)


### SUR: Repeatedly measured outcome at years 1, 2, and 3 vs. top genera, no AB during birth (n=198)

model.sur123.genera.noAB <- lmer(scale.score ~ gender + mat.age.imp + epds.tr3.imp + year.num + 
                                   Lactobacillus + Streptococcus + Staphylococcus + Corynebacterium + 
                                   Stenotrophomonas + Bacteroides + Prevotella + Anaerococcus + (1|id),
                                 data = ab.no.sur.long, REML=TRUE)

sum <- coef(summary(model.sur123.genera.noAB))
round(sum, 3)

### SUR: Repeatedly measured outcome at years 1, 2, and 3 vs. top genera, ABs used during birth (n=198)

model.sur123.genera.AB <- lmer(scale.score ~ gender + mat.age.imp + epds.tr3.imp + year.num + 
                                 Lactobacillus + Streptococcus + Staphylococcus + Corynebacterium + 
                                 Stenotrophomonas + Bacteroides + Prevotella + Anaerococcus + (1|id),
                               data = ab.yes.sur.long, REML=TRUE)

sum <- coef(summary(model.sur123.genera.AB))
round(sum, 3)

