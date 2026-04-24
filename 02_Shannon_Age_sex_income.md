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
```
--> 330 samples removed

--

## 2. Age sex and and income analyses

### 2.1 Sample distriburions

Age categogy
```r
Subset$sex <- recode(Subset$sex,"female" = "Female", "male"   = "Male")

plot_df <- Subset %>%
  filter(
    !is.na(ARG_div_shan), !is.na(sex), !is.na(age_category)) %>%
  mutate(
    Income_group = case_when(
      World_Bank_Income_Group == "High income" ~ "HIC",
      World_Bank_Income_Group %in% c("Low income", "Lower middle income", "Upper middle income") ~ "LMIC",
      TRUE ~ NA_character_), Income_group = factor(Income_group, levels = c("HIC", "LMIC")),
    age_category = factor(age_category, levels = c(
      "Infant", "Toddler", "Child", "Teenage",
      "Young adult", "Middle-Age Adult",
      "Older Adult", "Oldest Adult"))) %>%
  filter(!is.na(Income_group))

table_df <- plot_df %>%
  count(Income_group, age_category, sex) %>%
  pivot_wider(
    names_from = sex,
    values_from = n,
    values_fill = 0
  ) %>% mutate(Total = Female + Male)

table_image <- table_df %>%
  gt(groupname_col = "Income_group") %>%
  cols_label(
    age_category = "Age Category",
    Female = "Female (n)",
    Male = "Male (n)",
    Total = "Total (n)"
  ) %>%
  tab_header(title = "Sample Distribution by Age, Sex, and Income Group (HIC vs LMIC)")

table_image
```

Numeric Age
```r
plot_df <- Subset %>% filter(
    !is.na(ARG_div_shan),
    !is.na(sex),
    !is.na(age_years)) %>%
  mutate(
    sex = recode(sex, "female" = "Female", "male" = "Male"),
    Income_group = case_when(
      World_Bank_Income_Group == "High income" ~ "HIC",
      World_Bank_Income_Group %in% c("Low income", "Lower middle income", "Upper middle income") ~ "LMIC",
      TRUE ~ NA_character_),
    Income_group = factor(Income_group, levels = c("HIC", "LMIC")),
    age_category = factor(age_category, levels = c(
      "Infant", "Toddler", "Child", "Teenage",
      "Young adult", "Middle-Age Adult",
      "Older Adult", "Oldest Adult"))) %>%filter(!is.na(Income_group))

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

### 2.2 Linear model
With numeric age
```r
HIC_df <- plot_df %>% filter(Income_group == "HIC")
LMIC_df <- plot_df %>% filter(Income_group == "LMIC")

lm_HIC <- lm(ARG_div_shan ~ age_years + sex, data = HIC_df)
summary(lm_HIC)

lm_LMIC <- lm(ARG_div_shan ~ age_years + sex, data = LMIC_df)
summary(lm_LMIC)
```
Call:
lm(formula = ARG_div_shan ~ age_years + sex, data = HIC_df)

Residuals:
     Min       1Q   Median       3Q      Max 
-1.73511 -0.28980  0.02956  0.32344  1.43123 

Coefficients:
              Estimate Std. Error t value Pr(>|t|)    
(Intercept)  1.6989980  0.0106261 159.889   <2e-16 ***
age_years    0.0058535  0.0002085  28.077   <2e-16 ***
sexMale     -0.0084761  0.0104798  -0.809    0.419    
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.4752 on 8254 degrees of freedom
Multiple R-squared:  0.08723,	Adjusted R-squared:  0.08701 
F-statistic: 394.4 on 2 and 8254 DF,  p-value: < 2.2e-16


Call:
lm(formula = ARG_div_shan ~ age_years + sex, data = LMIC_df)

Residuals:
    Min      1Q  Median      3Q     Max 
-1.2185 -0.2862  0.0118  0.2661  1.3388 

Coefficients:
              Estimate Std. Error t value Pr(>|t|)    
(Intercept)  1.9532166  0.0279786  69.811  < 2e-16 ***
age_years    0.0037365  0.0005602   6.670 4.25e-11 ***
sexMale     -0.0276851  0.0265875  -1.041    0.298    
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.4174 on 984 degrees of freedom
Multiple R-squared:  0.04401,	Adjusted R-squared:  0.04207 
F-statistic: 22.65 on 2 and 984 DF,  p-value: 2.412e-10



### 2.3 Boxplot
With Categorised age
```r
plot_df <- Subset %>% filter(
    !is.na(ARG_div_shan),
    !is.na(sex),
    !is.na(age_category)) %>%
  mutate(
    sex = recode(sex, "female" = "Female", "male" = "Male"),
    Income_group = case_when(
      World_Bank_Income_Group == "High income" ~ "HIC",
      World_Bank_Income_Group %in% c("Low income", "Lower middle income", "Upper middle income") ~ "LMIC",
      TRUE ~ NA_character_),
    Income_group = factor(Income_group, levels = c("HIC", "LMIC")),
    age_category = factor(age_category, levels = c(
      "Infant", "Toddler", "Child", "Teenage",
      "Young adult", "Middle-Age Adult",
      "Older Adult", "Oldest Adult"))) %>%filter(!is.na(Income_group))

plot_df <- plot_df %>%
  mutate(
    age_category = factor(age_category, levels = c(
      "Infant", "Toddler", "Child", "Teenage",
      "Young adult", "Middle-Age Adult",
      "Older Adult", "Oldest Adult"
    )),
    Income_group = factor(Income_group, levels = c("HIC", "LMIC")),
    sex = recode(sex, "female" = "Female", "male" = "Male"))

ggplot(plot_df, aes(x = age_category, y = ARG_div_shan, fill = sex)) +
  geom_jitter(
    position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.6),
    size = 1.2, alpha = 0.25, color = "grey30"
  ) +
  geom_boxplot(
    width = 0.55, outlier.shape = NA, alpha = 0.8,
    position = position_dodge(width = 0.6)
  ) +
  facet_wrap(~Income_group) +  # separate HIC vs LMIC
  scale_fill_npg() +
  scale_color_npg() +
  labs(
    title = "ARG Shannon Diversity by Age Category and Sex",
    x = "Age Category",
    y = expression("Shannon Diversity"),
    fill = "Sex"
  ) +
  stat_compare_means(
    aes(group = sex), 
    label = "p.signif",    # shows significance symbols
    method = "wilcox.test" # non-parametric comparison
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "plain"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold")
  )
```

### 2.4 LOESS

--> Not statistically best, but nice for visualization. Not in article.
```
ggplot(plot_df, aes(x = age_years, y = ARG_div_shan, color = sex)) +
  geom_point(alpha = 0.2, size = 1.0) +
  geom_smooth(method = "loess", se = TRUE, span = 0.6) +
  facet_wrap(~Income_group) +
  scale_color_npg() +
  labs(
    x = "Age (years)",
    y = "Shannon diversity",
    color = "Sex",
    title = "LOESS-smoothed ARG Shannon diversity by age, sex, and income group"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "lightgrey", color = NA))
```

### 2.5 GAM 

```r
df_temp <- Subset |>
  dplyr::select(ARG_div_shan, age_years, sex, World_Bank_Income_Group) |>
  dplyr::mutate(
    sex = factor(sex, levels = c("Male", "Female")),
    age_years = as.numeric(age_years),
    Income_group = dplyr::case_when(
      World_Bank_Income_Group == "High income" ~ "HIC",
      World_Bank_Income_Group %in% c(
        "Low income",
        "Lower middle income",
        "Upper middle income"
      ) ~ "LMIC",
      TRUE ~ NA_character_
    ),
    Income_group = factor(Income_group, levels = c("HIC", "LMIC"))
  ) |>
  dplyr::filter(
    sex %in% c("Male", "Female"),
    !is.na(Income_group),
    !is.na(age_years))


fit <- gam(ARG_div_shan ~ sex + s(age_years, by = sex, k = 10) + sex*Income_group,
  method = "REML", data= df_temp)

summary(fit)
draw(fit)
```
| Model | Family   | Link     | Adjusted R² | Deviance explained | REML   | Scale   | N     |
| ----- | -------- | -------- | ----------- | ------------------ | ------ | ------- | ----- |
| GAM   | Gaussian | Identity | 10.9%       | **11%**            | 6076.7 | 0.21631 | 9,244 |


| Term                | Estimate | Std. Error | t value |  p-value | Significance |
| ------------------- | -------: | ---------: | ------: | -------: | ------------ |
| Intercept           |   1.9138 |    0.00742 |  258.11 |   <2e-16 | ***          |
| Sex (Female)        |  0.00451 |    0.01034 |   0.436 |    0.663 |              |
| Income group (LMIC) |  0.15471 |    0.02204 |   7.020 | 2.37e-12 | ***          |
| Sex (Female) × LMIC |  0.02834 |    0.03161 |   0.896 |    0.370 |              |


| Smooth term          |   edf | Ref.df | F value | p-value | Significance |
| -------------------- | ----: | -----: | ------: | ------: | ------------ |
| s(age_years): Male   | 4.795 |  5.839 |  110.76 |  <2e-16 | ***          |
| s(age_years): Female | 7.646 |  8.374 |   45.72 |  <2e-16 | ***          |


# Model check:
--> G through with Mahkameh
```r

# Some exterme values (above 1.5)
# Over 1.0 large resuduals
# 0 to 0.5 normal deviation
```

### 2.6 GAM differences

```r
diff <- difference_smooths(
  fit,
  smooth = "s(age_years)",
  n = 200)

draw(diff)
```


### 2.7 Visualize GAM differences
```r
newdata <- expand.grid(age_years = seq(min(df_temp$age_years), max(df_temp$age_years), length = 200),
  sex = c("Female","Male"),
  Income_group = levels(df_temp$Income_group))  # changed from World_Bank_Income_Group)

predictions <- predict(fit, newdata, se.fit = TRUE)
newdata$pred <- predictions$fit
newdata$se <- predictions$se.fit

diff_data <- newdata %>%
  tidyr::pivot_wider(
    names_from = sex,
    values_from = c(pred, se))

diff_data$diff <- diff_data$pred_Male - diff_data$pred_Female
```

### Options for visualizations

--> Options
```
# 1. Curves

ggplot(diff_data, aes(age_years, diff)) +
  geom_smooth() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~Income_group) +  # changed from World_Bank_Income_Group
  labs(
    y = "Male - Female predicted Shannon Diversity",
    x = "Age (years)",
    title = "Sex Differences in Predicted Shannon Diversity Across Age by Income Group"
  ) +
  ggsci::scale_color_npg() +
  ggsci::scale_fill_npg() +
  theme_minimal()

# 2. Heatmap

ggplot(diff_data, aes(age_years, Income_group, fill = diff)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "#b40426",   # red (Female > Male)
    mid = "white",
    high = "#3b4cc0",  # blue (Male > Female)
    midpoint = 0,
    limits = c(-max(abs(diff_data$diff)),
               max(abs(diff_data$diff)))
  ) +
  theme_minimal()

# Blue = Male > Female
# Red = Female > Male
# White = no difference

# 3.Trajectories of predicted means: 

ggplot(newdata, aes(age_years, pred, color = sex, fill = sex)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = pred - 1.96*se,
                  ymax = pred + 1.96*se),
              alpha = 0.2, color = NA) +
  facet_wrap(~Income_group) +
  labs(
    y = "Predicted log10(ARG load)",
    x = "Age (years)",
    title = "Age trajectories of ARG diversity by sex and income group"
  ) +
  theme_minimal()

```

