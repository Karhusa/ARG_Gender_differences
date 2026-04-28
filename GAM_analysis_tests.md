## 1. Prepare the data
### 1.1 Load packages
```r
library(SummarizedExperiment)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggsci)
library(mgcv)
library(ggpubr)
library(tidyverse)

#dir.create("~/R/library", recursive = TRUE)
.libPaths(c("~/R/library", .libPaths()))
.libPaths() #[1] "/users/saankarh/R/library"
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
# BiocManager::install("mia")
# install.packages("gt")
# install.packages("gratia")    

library(gt)
library(gratia)
library(mia)

```
### 1.2 Load TSE
```r
tse_path <- "/scratch/project_2008149/USER_WORKSPACES/karhula/DATA/TSE.rds"
TSE <- readRDS(tse_path)
```

### 1.3 Extract colData
```r
colData_df <- as.data.frame(colData(TSE))

Subset <- colData_df %>%
  select(
    log10_ARG_load,
    ARG_div_shan,
    region,
    city,
    sex,
    country,
    age_years,
    age_category,
    World_Bank_Income_Group,
    country,
    continent,
    parent,
    pregnant,
    bioproject,
    biosample)

Subset[Subset == ""] <- NA
```
### 1.4 Remove samples from Bioproject PRJNA544527
```r
dim(Subset)

remove_ids <- c(
  paste0("SAMN", 11950082:11950286),
  paste0("SAMN", 11950288:11950349),
  paste0("SAMN", 11950351:11950423)
)

Subset <- Subset[!Subset$biosample %in% remove_ids, ]
dim(Subset)
```
1. dim: 60996    14
2. dim: 60666    14
--> 330 samples removed

---

## 2. Sample distribution

```r
plot_df <- Subset %>%
  filter(!is.na(log10_ARG_load), !is.na(sex), !is.na(age_years)
  ) %>%
  mutate(
    sex = recode(sex, "female" = "Female", "male" = "Male"),
    sex = factor(sex),
    Income_group = case_when(
      World_Bank_Income_Group == "High income" ~ "HIC",
      World_Bank_Income_Group %in% c(
        "Low income",
        "Lower middle income",
        "Upper middle income"
      ) ~ "LMIC",
      TRUE ~ NA_character_
    ),
    Income_group = factor(Income_group, levels = c("HIC", "LMIC"))
  ) %>% filter(!is.na(Income_group))

table_df <- plot_df %>%
  count(Income_group, sex) %>%
  pivot_wider(
    names_from = sex,
    values_from = n,
    values_fill = 0
  ) %>% mutate(Total = Female + Male)

table_image <- table_df %>%
  gt() %>%
  cols_label(
    Income_group = "Income Group",
    Female = "Female (n)",
    Male = "Male (n)",
    Total = "Total (n)"
  ) %>% tab_header(title = "Total Sample Counts by Sex and Income Group")

table_image
```
## 3. Original Gam model 

* How ARG abundance changes
* Difference in baseline ARG load between males and females
* separate nonlinear age curves for each sex
* Does aging affect ARG load differently in males vs females in a nonlinear way?
* there is NO shared age trend, all age structure is sex-specific
* Does the effect of income group differ between males and females”

Question: Does ARG load vary with age differently in males and females, and is this modified by socioeconomic status?

### 3.1 Model
```r
fit <- gam(log10_ARG_load ~ sex + s(age_years, by = sex, k = 10) + sex*Income_group,
  method = "REML", data= plot_df)
```
### 3.2 Model check
```r
summary(fit)
draw(fit)
gam.check(fit)
```
### 3.3 Cross validation
```r
set.seed(1)

K <- 5
fold_id <- sample(rep(1:K, length.out = nrow(plot_df)))
rmse <- numeric(K)

for (k in 1:K) {

  train <- plot_df[fold_id != k, ]
  test  <- plot_df[fold_id == k, ]

  fit_test <- gam(
    log10_ARG_load ~ sex +
      s(age_years, by = sex, k = 10) +
      sex * Income_group,
    method = "REML",
    data = train
  )

  pred <- predict(fit_test, newdata = test)

  rmse[k] <- sqrt(mean((test$log10_ARG_load - pred)^2))
}

mean(rmse)
#[1] 0.3028779

sd(plot_df$log10_ARG_load)
#[1] 0.3121431
```
--> Gam effects are real but weak in predictive power

--> Your model explains population-level trends, but not individual-level variation


## 4 GAM analyses separated by income groups

* Within HIC, what is the age/sex pattern?
* Within LMIC, what is the age/sex pattern?

### 4.1 Model
```r
HIC_df <- plot_df %>% filter(Income_group == "HIC")
LMIC_df <- plot_df %>% filter(Income_group == "LMIC")

fit_HIC <- gam(log10_ARG_load ~ sex + s(age_years, by = sex, k = 10),
  method = "REML", data= HIC_df)

fit_LMIC <- gam(log10_ARG_load ~ sex + s(age_years, by = sex, k = 10),
  method = "REML", data= LMIC_df)
```
### 4.1 Model check
```
summary(fit_HIC)
draw(fit_HIC)
gam.check(fit_HIC)

summary(fit_LMIC)
draw(fit_LMIC)
gam.check(fit_LMIC)
```
### 4.2. Cross validation (Root Mean Squared Error)
HIC
```r
set.seed(1)

K <- 5
folds <- sample(rep(1:K, length.out = nrow(HIC_df)))

rmse <- numeric(K)

for (k in 1:K) {

  train <- HIC_df[folds != k, ]
  test  <- HIC_df[folds == k, ]

  model <- gam(
    log10_ARG_load ~ sex + s(age_years, by = sex, k = 10),
    method = "REML",
    data = train
  )

  pred <- predict(model, newdata = test)

  rmse[k] <- sqrt(mean((test$log10_ARG_load - pred)^2))
}

mean(rmse)
sd(rmse)

mean(rmse)
sd(rmse)
```
LMIC:
```r
set.seed(1)

K <- 5
folds <- sample(rep(1:K, length.out = nrow(LMIC_df)))

rmse <- numeric(K)

for (k in 1:K) {

  train <- LMIC_df[folds != k, ]
  test  <- LMIC_df[folds == k, ]

  model <- gam(
    log10_ARG_load ~ sex + s(age_years, by = sex, k = 10),
    method = "REML",
    data = train
  )

  pred <- predict(model, newdata = test)

  rmse[k] <- sqrt(mean((test$log10_ARG_load - pred)^2))
}

mean(rmse)
sd(rmse)

mean(rmse)
sd(rmse)

```

* Mean CV RMSE ≈ 0.303	weak predictive power
* SD of RMSE ≈ 0.007	very stable model
* different train/test splits give almost identical performance
* no evidence of instability or overfitting fluctuations

---

## 5 GAM analyses separated by sexes
### 5.1 Model
```r
Male_df   <- plot_df %>% filter(sex == "Male")
Female_df <- plot_df %>% filter(sex == "Female")

fit_female <- gam(
  log10_ARG_load ~ s(age_years, k = 10) + Income_group,
  method = "REML",
  data = Female_df
)

fit_male <- gam(
  log10_ARG_load ~ s(age_years, k = 10) + Income_group,
  method = "REML",
  data = Male_df
)
```

