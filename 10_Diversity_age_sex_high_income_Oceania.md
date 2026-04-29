
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

## 2.High income countries in Oceania

### 2.1 Sample distribution

```r
hic_Os <- Subset %>%
  filter(
    World_Bank_Income_Group == "High income",
    continent == "Oceania")

unique(hic_Os$country)

```
Age category (for boxplot)

```r
hic_Os_clean <- hic_Os %>%
  filter(!is.na(age_category), !is.na(sex))
table(hic_Os_clean$sex, useNA = "ifany")

table_df <- hic_Os_clean %>%
  count(age_category, sex) %>%
  tidyr::pivot_wider(
    names_from = sex,
    values_from = n,
    values_fill = 0
  ) %>% mutate(Total = female + male)

table_df <- table_df %>%
  mutate(
    age_category = factor(age_category, levels = c(
      "Infant", "Toddler", "Child", "Teenage",
      "Young adult", "Middle-Age Adult",
      "Older Adult", "Oldest Adult"))
  ) %>% arrange(age_category)

table_image <- table_df %>%
  gt() %>%
  cols_label(
    age_category = "Age Category",
    female = "Female (n)",
    male = "Male (n)",
    Total = "Total (n)"
  ) %>%
  tab_header(title = "High-Income Oceania: Sample Distribution by Age and Sex")
```
Numeric age (for linear regression)

```r
hic_Os_clean_num <- hic_Os %>% filter(!is.na(age_years), !is.na(sex))
sex_table <- hic_Os_clean_num %>% count(sex)

sex_gt <- sex_table %>%
  gt() %>%
  cols_label(
    sex = "Sex",
    n = "Count (n)"
  ) %>%tab_header(title = "Sex distribution in the High-Income Oceanian cohort")
```
## 2.2 Boxplot

```r
plot_df <- hic_Os_clean %>%
  filter(
    !is.na(ARG_div_shan),
    !is.na(sex),
    !is.na(age_category)
  ) %>%
  mutate(
    sex = case_when(
      sex %in% c("female", "Female") ~ "Female",
      sex %in% c("male", "Male") ~ "Male",
      TRUE ~ NA_character_
    ),
    sex = factor(sex, levels = c("Female", "Male")),
    
    age_category = factor(age_category, levels = c(
      "Infant", "Toddler", "Child", "Teenage",
      "Young adult", "Middle-Age Adult",
      "Older Adult", "Oldest Adult"
    ))
  ) %>%
  filter(!is.na(sex))

ggplot(plot_df, aes(x = age_category, y = ARG_div_shan, fill = sex)) +
  geom_jitter(
    position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.6),
    size = 1.2, alpha = 0.25, color = "grey30"
  ) +
  geom_boxplot(
    width = 0.55,
    outlier.shape = NA,
    alpha = 0.8,
    position = position_dodge(width = 0.6)
  ) +
  scale_fill_manual(values = c(
    "Female" = "#8B0000",
    "Male"   = "#003366"
  )) +
  labs(
    title = "ARG Diversity by Age Category and Sex (High-Income Oceania)",
    x = "Age Category",
    y = expression(" ARG Diversity, shannon index"),
    fill = "Sex"
  ) +
  stat_compare_means(
    aes(group = sex),
    method = "wilcox.test",
    label = "p.signif"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold")
  )

```

### 2.3 Linear regression
```r
model <- lm(ARG_div_shan ~ sex + age_years, data = hic_Os_clean_num)
summary(model)
```

