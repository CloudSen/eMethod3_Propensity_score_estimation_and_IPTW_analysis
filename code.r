# eMethod3. Propensity score estimation and IPTW analysis (R code)
library(survey)
library(ResourceSelection)
library(pROC)
library(car)

## 1) Read data
dat <- read.csv("complete_data.csv", header = TRUE, stringsAsFactors = FALSE)

## 2) Set variable types (baseline/preoperative covariates ONLY)
cat_vars <- c(
  "Sex","History_of_hypertension","Diabetes_mellitus",
  "Chronic_obstructive_pulmonary_disease","Coronary_heart_disease",
  "History_of_tumor","History_of_tuberculosis","Preoperative_cough",
  "Preoperative_pulmonary_infection","Preoperative_hypoxemia",
  "History_of_smoking","Alcohol_consumption"
)
cat_vars <- intersect(cat_vars, names(dat))
dat[cat_vars] <- lapply(dat[cat_vars], factor)

# exposure: must be 0/1 numeric
dat$exposure <- as.integer(dat$exposure)
stopifnot(all(is.na(dat$exposure) | dat$exposure %in% c(0L, 1L)))

# outcome: must be 0/1 numeric
dat$antitussive <- as.integer(dat$antitussive)
stopifnot(all(is.na(dat$antitussive) | dat$antitussive %in% c(0L, 1L)))

# optional: set ordered factors if applicable
if ("ASA" %in% names(dat)) dat$ASA <- ordered(dat$ASA)
if ("Type_of_lung_resection" %in% names(dat)) dat$Type_of_lung_resection <- factor(dat$Type_of_lung_resection)

## 3) Propensity score model (logistic regression; pre-treatment covariates only)
ps_model <- glm(
  exposure ~ Sex + Age + BMI + Chronic_obstructive_pulmonary_disease +
    Preoperative_cough + CHARLSON + ASA + Preoperative_hypoxemia +
    History_of_tuberculosis + Preoperative_pulmonary_infection +
    Coronary_heart_disease + Type_of_lung_resection + FEV1,
  data   = dat,
  family = binomial()
)

ps <- predict(ps_model, type = "response")

## Optional: check multicollinearity
vif(ps_model)

## Model discrimination (AUC)
cc <- complete.cases(dat$exposure, ps)
roc_obj <- roc(dat$exposure[cc], ps[cc], quiet = TRUE)
auc(roc_obj)

## 4) Stabilized IPTW
pt <- mean(dat$exposure, na.rm = TRUE)
eps <- 1e-6
ps  <- pmin(pmax(ps, eps), 1 - eps)

w1 <- ifelse(dat$exposure == 1, pt / ps, (1 - pt) / (1 - ps))

## 5) Weight truncation at 99th percentile
cut <- quantile(w1, 0.99, na.rm = TRUE)
dat$w_iptw <- pmin(w1, cut)

## 6) Effective Sample Size (ESS)
ESS <- (sum(dat$w_iptw, na.rm = TRUE)^2) / sum(dat$w_iptw^2, na.rm = TRUE)
ESS

## 7) Survey design
des_iptw <- svydesign(ids = ~1, weights = ~w_iptw, data = dat)

## 8) IPTW-weighted logistic regression
fit <- svyglm(antitussive ~ exposure, design = des_iptw, family = quasibinomial())
beta <- coef(fit)["exposure"]
se   <- sqrt(vcov(fit)["exposure","exposure"])

IPTW_logit_results <- data.frame(
  outcome = "antitussive",
  OR  = exp(beta),
  LCL = exp(beta - 1.96*se),
  UCL = exp(beta + 1.96*se),
  p   = summary(fit)$coefficients["exposure", 4],
  stringsAsFactors = FALSE
)

IPTW_logit_results
