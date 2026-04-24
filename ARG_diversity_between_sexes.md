## 1. ARG diversity between sexes

### 1.2 Sample distribution (ARG div x sex)

```r
Subset_shan <- Subset %>%
  filter(!is.na(ARG_div_shan), !is.na(sex))

npg_cols <- pal_npg("nrc")(4)[3:4]

ggplot(Subset_shan, aes(x = ARG_div_shan, fill = sex)) +
  geom_histogram(bins = 30, color = "black", alpha = 0.8, position = "identity") +
  facet_wrap(~ sex, scales = "fixed") +  # separate panels for each sex
    scale_fill_manual(values = npg_cols) +
  labs(
    title = "Distribution of Shannon Diversity Index by Sex",
    x = "Shannon Diversity Index",
    y = "Sample Count"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

Shan_summary <- Subset_clean %>%
  group_by(sex) %>%
  summarise(
    n = n(),                                # sample count
    mean_shan = mean(ARG_div_shan),        # mean
    sd_shan = sd(ARG_div_shan),            # standard deviation
    median_shan = median(ARG_div_shan),    # median
    min_shan = min(ARG_div_shan),          # minimum
    max_shan = max(ARG_div_shan)           # maximum
  )
Shan_summary
```
![Sample distribution ARG diversity](https://github.com/Karhusa/F_AMR_project/blob/main/Results/Thesis/Sample_distribution_ARG_diversity.png)

### 2.2 Boxplot of shannon index and sex
```r

Subset$sex[Subset$sex == "" | Subset$sex == "NA"] <- NA
Subset$ARG_div_shan[Subset$ARG_div_shan == "" | Subset$ARG_div_shan == "NA"] <- NA

Subset$sex <- recode(Subset$sex,
                     "female" = "Female",
                     "male"   = "Male")

plot_df <- Subset %>% filter(!is.na(sex), !is.na(ARG_div_shan))

n_df <- plot_df %>%
  group_by(sex) %>%
  summarise(
    n = n(),
    y_max = max(ARG_div_shan, na.rm = TRUE)
  )

npg_cols <- pal_npg("nrc")(4)[3:4]

ggplot(plot_df, aes(x = sex, y = ARG_div_shan, fill = sex)) +
  geom_jitter(
    width = 0.15,
    size = 1.2,
    alpha = 0.25,
    color = "grey30"
  ) +
  geom_boxplot(
    width = 0.55,
    outlier.shape = NA,
    alpha = 0.8,
    color = "grey25"
  ) +
  geom_text(
    data = n_df,
    aes(
      x = sex,
      y = y_max + 0.05,
      label = paste0("N = ", n)
    ),
    inherit.aes = FALSE,
    size = 4
  ) +
  scale_fill_manual(values = npg_cols) +
  labs(
    title = "ARG Shannon diversity by Sex",
    x = "Sex",
    y = "ARG Shannon diversity"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold"),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10))
  )

ggsave("Boxplot_Shannon_diversity_by_sex_age.png", width = 8, height = 6, dpi = 300)
```

![Shannon diversity by sex by Sex](https://github.com/Karhusa/F_AMR_project/blob/main/Results/Shannon_Analyses/Boxplot_Shannon_diversity_by_sex_age.png)

### 3.1.2 Regression analysis of shannon index and sex
```r
plot_df <- Subset %>% filter(!is.na(sex), !is.na(ARG_div_shan))
plot_df$sex <- factor(plot_df$sex, levels = c("Female", "Male"))

lm_sex <- lm(ARG_div_shan ~ sex, data = plot_df)
summary(lm_sex)
```





