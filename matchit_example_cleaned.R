library(MatchIt)
library(dplyr)
library(ggplot2)

setwd('PropensityScores/Scripts')
sink("sink-examp.txt")



#-----------------
#settings
#-----------------
path_to_loan_data <- "Scripts/datasets/loan.csv"
path_to_catholic_vs_math_data <- "Scripts/datasets/catholic_school_vs_math.csv"

#-----------------
#raw_data
#-----------------
#dependent on parameters
path_to_data <- path_to_catholic_vs_math_data
ds <- read.csv(path_to_data)

print(path_to_data)

if (path_to_data == path_to_catholic_vs_math_data)
{
  ds_cov <- c('race_white', 'p5hmage', 'w3income', 'p5numpla', 'w3momed_hsb')
  trt <- 'catholic'
  target_metric <- 'c5r2mtsc_std'
  ds<-ds[which(p5numpla < 4),]
}

if (path_to_data == path_to_loan_data)
{
  ds <- ds %>% 
    mutate(loan_status_bin = if_else(loan_status == 'Fully Paid', 1, 0))
  
  ds_cov <- c('loan_amnt', 'funded_amnt', 'funded_amnt_inv', 'term', 'int_rate', 'installment', 'grade', 'sub_grade')
  trt <- 'loan_status_bin'
  target_metric <- 'None'
}

#--
ds['trt']<-ds[trt]
ds['target_metric'] <- ds[target_metric]

ds_nomiss <- ds %>%  # MatchIt does not allow missing values
  select(target_metric, trt, one_of(ds_cov)) %>%
  na.omit()

ds_nomiss_id <- ds %>%  # MatchIt does not allow missing values
  select(target_metric,'childid', trt, one_of(ds_cov)) %>%
  na.omit()
write.csv(ds_nomiss_id,'tmp.csv')

if (path_to_data == path_to_loan_data)
{
  ds_nomiss <- ds %>%  # MatchIt does not allow missing values
  select(trt, one_of(ds_cov)) %>%
  na.omit()
}

#-----------------
#Descriptive BEFORE the Matching
#-----------------
print('####################################')
print('#Descriptive BEFORE the Matching')
print('####################################')
ds_nomiss  %>%
  group_by(trt) %>%
  select(one_of(ds_cov)) %>%
  summarise_all(funs(mean(., na.rm = T)))

lapply(ds_cov, function(v) {
  t.test(ds_nomiss [, v] ~ ds_nomiss [, 'trt'])
})

#-----------------
#Analysis of target metric BEFORE matching
#-----------------
print('####################################')
print('#Analysis of target metric BEFORE matching')
print('####################################')
with(ds_nomiss , t.test(target_metric ~ trt))

lm_treat1_before <- lm(target_metric ~ trt, data = ds_nomiss)
summary(lm_treat1_before)

#dependent on parameters
if (path_to_data == path_to_catholic_vs_math_data)
{
  lm_treat2_before <- lm(target_metric ~ trt + race_white + p5hmage +
                    I(w3income / 10^3) + p5numpla + w3momed_hsb, data = ds_nomiss)
}

if (path_to_data == path_to_loan_data)
{
  print("NONE")
}
#--
summary(lm_treat2_before)


#-----------------
#Matching/Modeling
#-----------------

#dependent on parameters
if (path_to_data == path_to_catholic_vs_math_data)
{
  ds <- ds %>% mutate(w3income_1k = w3income / 1000)
  m_ps <- glm(trt ~ race_white + w3income_1k + p5hmage + p5numpla + w3momed_hsb,
              family = binomial(), data = ds)
}

if (path_to_data == path_to_loan_data)
{
  m_ps <- glm(trt ~ loan_amnt+funded_amnt+funded_amnt_inv+term+int_rate+installment+grade+sub_grade,
              family = binomial(), data = ds)
}

#--
print('####################################')
print('#Summary on PS')
print('####################################')
summary(m_ps)
prs_df <- data.frame(pr_score = predict(m_ps, type = "response"),
                     trt = m_ps$model$trt)
head(prs_df)

#--
labs <- paste("Treatment Type", c(1, 0))
prs_df %>%
  mutate(trt = ifelse(trt == 1, labs[1], labs[2])) %>%
  ggplot(aes(x = pr_score)) +
  geom_histogram(color = "white") +
  facet_wrap(~trt) +
  xlab("Probability of treatment") +
  theme_bw()

#dependent on parameters
if (path_to_data == path_to_catholic_vs_math_data)
{
  set.seed(20200429)
  mod_match <- matchit(trt ~ race_white + w3income + p5hmage + p5numpla + w3momed_hsb,
                       method = "nearest", data = ds_nomiss)
}

if (path_to_data == path_to_loan_data)
{
  mod_match <- matchit(trt ~ loan_amnt+funded_amnt+funded_amnt_inv+term+int_rate+installment+grade+sub_grade,
                       method = "nearest", data = ds_nomiss)
}

dta_m <- match.data(mod_match)
dim(dta_m)

summary(mod_match)
plot(mod_match)

#-----------------
#Some Visuals of Results after Matching
#-----------------
#dependent on parameters
if (path_to_data == path_to_catholic_vs_math_data)
{
  fn_bal <- function(dta, variable) {
    dta$variable <- dta[, variable]
    if (variable == 'w3income') dta$variable <- dta$variable / 10^3
    dta$trt <- as.factor(dta$trt)
    support <- c(min(dta$variable), max(dta$variable))
    ggplot(dta, aes(x = distance, y = variable, color = trt)) +
      geom_point(alpha = 0.2, size = 1.3) +
      geom_smooth(method = "loess", se = F) +
      xlab("Propensity score") +
      ylab(variable) +
      theme_bw() +
      ylim(support)
  }
  
  library(gridExtra)
  grid.arrange(
    fn_bal(dta_m, "w3income"),
    fn_bal(dta_m, "p5numpla") + theme(legend.position = "none"),
    fn_bal(dta_m, "p5hmage"),
    fn_bal(dta_m, "w3momed_hsb") + theme(legend.position = "none"),
    fn_bal(dta_m, "race_white"),
    nrow = 3, widths = c(1, 0.8)
  )
}

if (path_to_data == path_to_loan_data)
{
  print("NONE")
}

#-----------------
#Analysis of feature AFTER matching
#-----------------
#--
print('####################################')
print('#Analysis of feature AFTER matching')
print('####################################')
dta_m %>%
  group_by(trt) %>%
  select(one_of(ds_cov)) %>%
  summarise_all(funs(mean))

lapply(ds_cov, function(v) {
  t.test(dta_m[, v] ~ dta_m$trt)
})

#-----------------
#Analysis of target metric AFTER matching
#-----------------

with(dta_m, t.test(target_metric ~ trt))

lm_treat1 <- lm(target_metric ~ trt, data = dta_m)
summary(lm_treat1)

#dependent on parameters
if (path_to_data == path_to_catholic_vs_math_data)
{
  lm_treat2 <- lm(target_metric ~ trt + race_white + p5hmage +
                    I(w3income / 10^3) + p5numpla + w3momed_hsb, data = dta_m)
}

if (path_to_data == path_to_loan_data)
{
  print("NONE")
}
#--
summary(lm_treat2)

#-----------------
#summary_to the file
#-----------------
#dependent on parameters
if (path_to_data == path_to_catholic_vs_math_data)
{
  write.csv(mod_match$match.matrix, 'catholic_vs_math_data_match.matrix.csv')
  write.csv(row.names(mod_match$match.matrix), 'ind_catholic_vs_math_data_match.matrix.csv')
}
if (path_to_data == path_to_loan_data)
{
  write(mod_match$match.matrix, 'loan_data_match.matrix.csv')
}
sink()
print("end")
#--------


