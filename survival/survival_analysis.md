
# Load packages

``` r
library(here)       # For relative file paths
library(tidyverse)  # Includes ggplot2, dplyr, purrr, readr, tidyr, stringr
library(survival)   # For core survival analysis models
library(survminer)  # For ggplot-based survival curves (ggsurvplot)
```

# Preprocess the data

``` r
file_2024 <- here::here("survival", "data", "survival_2024_with_pd_number.csv")
file_2025 <- here::here("survival", "data", "survival_2025_with_pd_number.csv")

table_2024 <- read_delim(file_2024, delim = ";", locale = locale(encoding = "latin1"))
table_2025 <- read_delim(file_2025, delim=";",locale = locale(encoding = "latin1"))


transform_survival <- function(df, year_val) {
  
  # Fix encoding issues in column names and character columns
  colnames(df) <- iconv(colnames(df), to = "UTF-8", sub = "byte")
  df <- df %>% mutate(across(where(is.character), ~iconv(., to = "UTF-8", sub = "byte")))
  
  colnames(df)[1] <- "Category"
  colnames(df)[2] <- "Petridish_number"
  colnames(df)[3] <- "Start"
  
  time_row <- df %>% 
    filter(str_detect(tolower(Category), "time"))
  
  time_dict <- tibble(
    column_name = colnames(df)[3:ncol(df)],
    time_hours  = as.numeric(time_row[1, 3:ncol(df)])
  ) %>% filter(!is.na(time_hours))
  
  data_clean <- df %>%
    filter(!is.na(Category), !str_detect(tolower(Category), "time")) %>%
    mutate(across(3:ncol(df), ~as.numeric(as.character(.))))
  
  data_long <- data_clean %>%
    pivot_longer(
      cols      = any_of(time_dict$column_name),
      names_to  = "column_name",
      values_to = "alive_count"
    ) %>%
    left_join(time_dict, by = "column_name") %>%
    filter(!is.na(alive_count)) %>%
    arrange(Category, Petridish_number, time_hours)
  
  survival_events <- data_long %>%
    group_by(Category, Petridish_number) %>%
    mutate(
      new_dead = lag(alive_count, default = dplyr::first(alive_count)) - alive_count,
      new_dead = if_else(new_dead < 0, 0, new_dead)
    ) %>%
    ungroup()
  
  morts <- survival_events %>%
    filter(new_dead > 0) %>%
    uncount(new_dead) %>%
    mutate(status = 1, time = time_hours)
  
  censures <- survival_events %>%
    group_by(Category, Petridish_number) %>%
    slice_tail(n = 1) %>%
    ungroup() %>%
    filter(alive_count > 0) %>%
    uncount(alive_count) %>%
    mutate(status = 0, time = time_hours)
  
  final_df <- bind_rows(morts, censures) %>%
    mutate(
      Category_base = str_trim(tolower(Category)),
      Antibiotic    = case_when(
        str_detect(Category_base, "high") ~ "high",
        str_detect(Category_base, "low")  ~ "low",
        TRUE                              ~ "no"
      ),
      Reinoculation = if_else(str_detect(Category_base, "gut|homogenate"), "yes", "no"),
      Treatment     = paste0(Antibiotic, "_", Reinoculation),
      Year          = year_val
    ) %>%
    dplyr::select(Category_base, Antibiotic, Reinoculation, Year, 
                  Petridish_number, time, status, Treatment)
  
  return(final_df)
}

survivial_2024 <- transform_survival(table_2024, 2024)
survivial_2025 <- transform_survival(table_2025, 2025)

dataset_final <- bind_rows(survivial_2024, survivial_2025)

surv_individual <- dataset_final %>%
  mutate(
    Treatment = factor(Treatment),
    Treatment = relevel(Treatment, ref = "no_no"),
    Antibiotic = factor(Antibiotic, levels = c("no", "low", "high")),
    Reinoculation = factor(Reinoculation),
    Reinoculation = relevel(Reinoculation, ref = "no"),
    Year = factor(Year),
    Petridish_number = factor(Petridish_number),
    dish_id = interaction(Year, Category_base, Petridish_number, drop = TRUE))
```

# Fit a survival curve using the Kaplan-Meier method (descriptive, not stratified by year)

``` r
fit_ant <- survfit(Surv(time, status) ~ Antibiotic, type = "kaplan-meier", data=surv_individual)
fit_reino <- survfit(Surv(time, status) ~ Reinoculation, type = "kaplan-meier", data=surv_individual)
fit_inter <- survfit(Surv(time, status) ~ interaction(Antibiotic, Reinoculation), type = "kaplan-meier", data=surv_individual)
```

# Log-rank tests (stratified by year), individual factors

# NOTE: these are descriptive / sensitivity checks. The primary inferential

# test for the manuscript is the Cox robust Wald test below (m_add), since

# it accounts for dish-level clustering; the log-rank test does not.

``` r
surv_diff_ant   <- survdiff(Surv(time, status) ~ Antibiotic + strata(Year), data = surv_individual)
surv_diff_reino <- survdiff(Surv(time, status) ~ Reinoculation + strata(Year), data = surv_individual)
```

# Time-dependent re-inoculation (switches on at 72 h)

``` r
# Split follow-up at 72 h: re-inoculation is a time-dependent covariate
cp <- survSplit(Surv(time, status) ~ ., data = surv_individual,
                cut = 72, episode = "period")   # creates tstart / time / status

cp$reinoc_td <- factor(ifelse(cp$period == 2 & cp$Reinoculation == "yes", "yes", "no"),
                       levels = c("no", "yes"))

# Sanity check: deaths before / after 72 h
table(period = cp$period, status = cp$status)
table(period = cp$period, reinoc_td = cp$reinoc_td, status = cp$status)
```

# Cox models: strata = year, robust SE clustered by dish

# This is the PRIMARY model for the manuscript (dose effect, re-inoculation

# effect, and interaction should be reported from here, not from the

# unstratified-by-dish log-rank tests above).

``` r
# Additive model
m_add <- coxph(Surv(tstart, time, status) ~ Antibiotic + reinoc_td +
                 strata(Year) + cluster(dish_id), data = cp)

# Interaction model
m_int <- coxph(Surv(tstart, time, status) ~ Antibiotic * reinoc_td +
                 strata(Year) + cluster(dish_id), data = cp)

# Antibiotic dose effect: robust Wald test (2 df)
idx        <- grep("^Antibiotic", names(coef(m_add)))
b          <- coef(m_add)[idx]
V          <- vcov(m_add)[idx, idx]          # robust variance (cluster)
wald       <- as.numeric(t(b) %*% solve(V) %*% b)
p_ant_wald <- pchisq(wald, df = length(idx), lower.tail = FALSE)
c(chisq = wald, df = length(idx), p = p_ant_wald)

# Interaction test 1: likelihood ratio (anova() does not work with cluster(),
# so refit the same models without it)
m_add_nc  <- coxph(Surv(tstart, time, status) ~ Antibiotic + reinoc_td + strata(Year), data = cp)
m_int_nc  <- coxph(Surv(tstart, time, status) ~ Antibiotic * reinoc_td + strata(Year), data = cp)
lrt_inter <- anova(m_add_nc, m_int_nc)
lrt_inter

# Proportional hazards check
cox.zph(m_add)
```

# Create the plots

``` r
time_2024 <- c(0, 4, 21, 24, 26.5, 28.5, 45, 48, 50.5, 52.5, 69, 72, 115.5, 117, 118)
time_2025 <- c(0, 15, 18, 21, 24, 39, 42, 45, 48, 63, 66, 69, 72, 135, 138)

p_ant <- ggsurvplot(
  fit_ant,
  pval = FALSE,
  pval.size = 5,
  pval.method = TRUE,
  conf.int = FALSE,
  censor = FALSE,     
  surv.scale = "percent",
  xlab = "Time [h]",
  ylab = "Survival probability",
  

  ggtheme = theme_classic(), 
  legend.title = "Antibiotic dose",
  palette = c("#fcca46", "#fe7f2d", "#233d4d"), 
  linewidth = 1,
  alpha = 0.65
)

heures_observation <- sort(unique(c(time_2024, time_2025)))
df_surv <- surv_summary(fit_ant)

points_artificiels <- df_surv %>%
  group_by(strata) %>%
  do({
    strata_data <- .
    map_df(heures_observation, function(ch) {
      val_surv <- strata_data %>% filter(time <= ch) %>% slice_tail(n = 1)
      if(nrow(val_surv) == 0) {
        return(tibble(time = ch, surv = 1))
      } else {
        return(tibble(time = ch, surv = val_surv$surv))
      }
    })
  }) %>%
  ungroup()

p_ant$plot <- p_ant$plot +
  geom_point(
    data = points_artificiels, 
    aes(x = time, y = surv, color = strata), 
    shape = 4,       
    size = 1,        
    stroke = 0.8    
  ) +
  geom_vline(xintercept = 72, linetype = "dotted", linewidth = 0.7) +
  annotate(
    "text",
    x = 72, y = 0.20,
    label = "Re-inoculation",
    angle = 90,
    vjust = -0.5,
    color = "black",
    size = 4
  ) +

  
  theme(
    axis.line = element_line(color = "black"),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10, color = "black"),
    panel.grid.major.y = element_line(color = "grey95"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11)
  )

print(p_ant)
```

``` r
p_reino <- ggsurvplot(
  fit_reino,
  title = "Impact of Re-inoculation on survival",
  pval = FALSE,   # unstratified log-rank p-value removed; report the Cox result instead
  conf.int = FALSE,
  surv.scale = "percent",
  xlab = "Time [h]",
  ylab = "Survival probability",
  ggtheme = theme_survminer(base_size = 8, base_family = "Helvetica"),
  legend.title = "Re-inoculation (yes/no)",
  palette = "Dark2",
  linewidth = 1.0,
  alpha = 0.8
)
p_reino$plot <- p_reino$plot +
  geom_vline(xintercept = 72, linetype = "dotted", linewidth = 0.7) +
  annotate(
    "text",
    x = 72, y = 0.20,
    label = "Re-inoculation",
    angle = 90,
    vjust = -0.5,
    color = "black",
    size = 3
  )

p_reino
```

``` r
p_inter <- ggsurvplot(
  fit_inter,
  title = "Impact of Antibiotics and Re-inoculation on survival",
  pval = FALSE,
  conf.int = FALSE,
  surv.scale = "percent",
  xlab = "Time [h]",
  ylab = "Survival probability",
  ggtheme = theme_survminer(base_size = 8, base_family = "Helvetica"),
  legend.title = "Treatment",
  palette = "Dark2",
  linewidth = 1.0,
  alpha = 0.8
)
p_inter$plot <- p_inter$plot +
  geom_vline(xintercept = 72, linetype = "dotted", linewidth = 0.7) +
  annotate(
    "text",
    x = 72, y = 0.20,
    label = "Re-inoculation",
    angle = 90,
    vjust = -0.5,
    color = "black",
    size = 3
  )
p_inter
```

# Justifying the joining of zero and low

``` r
survdiff(Surv(time, status) ~ Antibiotic + strata(Year), data = subset(surv_individual, Antibiotic %in% c("no", "low")))
```

# Pairwise comparisons (stratified by year)

# NOTE: three comparisons are run on the same family of hypotheses

# (no vs low, no vs high, low vs high), so the raw p-values are corrected

# for multiple comparisons with a Holm step-down procedure below.

``` r
pair_no_low   <- survdiff(Surv(time, status) ~ Antibiotic + strata(Year),
                          data = subset(surv_individual, Antibiotic %in% c("no", "low")))
pair_no_high  <- survdiff(Surv(time, status) ~ Antibiotic + strata(Year),
                          data = subset(surv_individual, Antibiotic %in% c("no", "high")))
pair_low_high <- survdiff(Surv(time, status) ~ Antibiotic + strata(Year),
                          data = subset(surv_individual, Antibiotic %in% c("low", "high")))

# Raw p-values (1 df chi-square each)
pairwise_p_raw <- c(
  no_vs_low   = pair_no_low$pvalue,
  no_vs_high  = pair_no_high$pvalue,
  low_vs_high = pair_low_high$pvalue
)

# Holm correction across the 3 pairwise comparisons
pairwise_p_holm <- p.adjust(pairwise_p_raw, method = "holm")

pairwise_table <- data.frame(
  comparison = names(pairwise_p_raw),
  p_raw      = signif(pairwise_p_raw, 3),
  p_holm     = signif(pairwise_p_holm, 3)
)
pairwise_table
```

# Save the plot and stats

``` r
results_dir <- here::here("survival", "result")

if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}

ggsave(file.path(results_dir, "p_ant.pdf"),   plot = p_ant$plot,   dpi = 300, width = 10, height = 5)
ggsave(file.path(results_dir, "p_reino.pdf"), plot = p_reino$plot, dpi = 300, width = 8,  height = 4)
ggsave(file.path(results_dir, "p_inter.pdf"), plot = p_inter$plot, dpi = 300, width = 8,  height = 4)

output_file <- file.path(results_dir, "survival_result.txt")
sink(output_file)

cat("=========================================================\n")
cat("1. SURVIVAL FIT OBJECTS\n")
cat("=========================================================\n\n")
print(fit_ant)
cat("\n---------------------------------------------------------\n")
print(fit_reino)
cat("\n---------------------------------------------------------\n")
print(fit_inter)

cat("\n\n=========================================================\n")
cat("2. SURVIVAL DIFFERENCES (descriptive / sensitivity checks only)\n")
cat("=========================================================\n\n")
cat("The log-rank tests below do not account for dish-level clustering and\n")
cat("are reported here only as a sensitivity check.\n\n")
cat("--> Log-rank (stratified by year): Antibiotic\n")
print(surv_diff_ant)
cat("P-value:", surv_diff_ant$pvalue, "\n")

cat("\n---------------------------------------------------------\n")
cat("--> Log-rank (stratified by year): Reinoculation (from time 0, descriptive)\n")
print(surv_diff_reino)
cat("P-value:", surv_diff_reino$pvalue, "\n")

cat("\n---------------------------------------------------------\n")
cat("--> Deaths before / after 72 h\n")
print(table(period = cp$period, status = cp$status))
print(table(period = cp$period, reinoc_td = cp$reinoc_td, status = cp$status))

cat("\n---------------------------------------------------------\n")
cat("--> Cox additive model (PRIMARY TEST; time-dependent re-inoculation,\n")
cat("    strata(Year), cluster(dish)) -- report this for dose & reinoc effects\n")
print(summary(m_add))
cat("Antibiotic Wald test (robust): chisq =", wald, " df =", length(idx), " p =", p_ant_wald, "\n")

cat("\n---------------------------------------------------------\n")
cat("--> Interaction test: likelihood ratio (models without cluster)\n")
print(lrt_inter)
cat("(Robust Wald not reported: interaction model unstable, groups with no deaths after 72 h)\n")

cat("\n---------------------------------------------------------\n")
cat("--> Proportional hazards check\n")
print(cox.zph(m_add))

cat("\n\n=========================================================\n")
cat("3. PAIRWISE POST-HOC COMPARISONS (ANTIBIOTIC EFFECTS, stratified by year)\n")
cat("   Holm-corrected for 3 comparisons -- report p_holm, not raw p\n")
cat("=========================================================\n\n")

cat("--> Pairwise: Control (no) vs Low Dose (low)\n")
print(pair_no_low)
cat("P-value (raw):", pair_no_low$pvalue, "\n\n")

cat("---------------------------------------------------------\n")
cat("--> Pairwise: Control (no) vs High Dose (high)\n")
print(pair_no_high)
cat("P-value (raw):", pair_no_high$pvalue, "\n\n")

cat("---------------------------------------------------------\n")
cat("--> Pairwise: Low Dose (low) vs High Dose (high)\n")
print(pair_low_high)
cat("P-value (raw):", pair_low_high$pvalue, "\n")

cat("\n---------------------------------------------------------\n")
cat("--> Holm-adjusted pairwise p-values (report these)\n")
print(pairwise_table)

cat("\n\n=========================================================\n")
cat("4. VALUES REPORTED (Results text)\n")
cat("=========================================================\n\n")

# Antibiotic dose: robust Wald, 2 df (computed above)
chi_ab <- wald
df_ab  <- length(idx)
p_ab   <- p_ant_wald

# Re-inoculation: robust Wald, 1 df (chi2 = robust z^2)
cf <- summary(m_add)$coefficients
reinoc_row <- grep("^reinoc_td", rownames(cf))
stopifnot("Could not find a unique reinoc_td coefficient - check contrasts()" =
            length(reinoc_row) == 1)
chi_re <- unname(cf[reinoc_row, "z"]^2)
df_re  <- 1
p_re   <- unname(cf[reinoc_row, "Pr(>|z|)"])

# Interaction: likelihood ratio (second row of the anova table)
chi_in <- lrt_inter$Chisq[2]
df_in  <- lrt_inter$Df[2]
p_in   <- lrt_inter$`Pr(>|Chi|)`[2]

fmt_p <- function(p) ifelse(p < 0.001, "p < 0.001", sprintf("p = %.3f", p))

# Deaths / N
mort_ab     <- surv_individual %>% group_by(Antibiotic)    %>% summarise(d = sum(status), n = n())
mort_re_all <- surv_individual %>% group_by(Reinoculation) %>% summarise(d = sum(status), n = n())
mort_re_72  <- cp %>% filter(period == 2) %>% group_by(Reinoculation) %>% summarise(d = sum(status), n = n())

cat("--> Test summary \n")
print(data.frame(
  test  = c("Antibiotic (robust Wald)", "Re-inoculation (robust Wald)", "Interaction (LRT)"),
  chisq = round(c(chi_ab, chi_re, chi_in), 2),
  df    = c(df_ab, df_re, df_in),
  p     = signif(c(p_ab, p_re, p_in), 3)
))

cat("\n--> Deaths / N by antibiotic (whole follow-up)\n");              print(as.data.frame(mort_ab))
cat("\n--> Deaths / N by re-inoculation (whole follow-up)\n");          print(as.data.frame(mort_re_all))
cat("\n--> Deaths / N by re-inoculation, after 72 h only (N = bees alive at 72 h)\n"); print(as.data.frame(mort_re_72))

ab <- setNames(paste0(mort_ab$d, "/", mort_ab$n), mort_ab$Antibiotic)
rp <- setNames(paste0(mort_re_72$d, "/", mort_re_72$n), mort_re_72$Reinoculation)

sink()

file.show(output_file)
```

## Session info

``` r
sessionInfo()
```
