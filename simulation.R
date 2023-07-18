## simulation

library(simex)
library(data.table)
library(ggplot2)
library(xtable)
library(nonprobsvy)
library(sampling)
library(RANN)
library(Rcpp)
library(survey)
library(simputation)
library(scales)
library(ggh4x)

pop_data <- fread("pop-example.csv", dec = ",")
pop_data[, short_n := pop*short]
pop_data[, long_n := pop - pop*short]
pop_data_1 <- pop_data[rep(1:nrow(pop_data), pop_data$short_n), .(age, gender, country)][, short:=1]
pop_data_2 <- pop_data[rep(1:nrow(pop_data), pop_data$long_n), .(age, gender, country)][, short:=0]
pop_data <- rbind(pop_data_1, pop_data_2)
nonprob_df <- pop_data[, .N, .(short)] ## prob on age
nonprob_df[, prob:=c(0.5, 0.7)]
pop_data[nonprob_df, prob:=i.prob, on = c("short")]
pop_data[, flag:=rbinom(nrow(pop_data), 1, prob)]

## age 18-24 25-29 30-39   40+ 
ages <- c("18-24", "25-29", "30-39", "40+")
pop_data[age == ages[1], age_m := rep(x = ages, times = round(sum(age == ages[1])*c(0.75, 0.18, 0.06, 0.01)))]
pop_data[age == ages[2], age_m := rep(x = ages, times = round(sum(age == ages[2])*c(0.12, 0.80, 0.07, 0.01)))]
pop_data[age == ages[3], age_m := rep(x = ages, times = round(sum(age == ages[3])*c(0.01, 0.04, 0.85, 0.10)))]
pop_data[age == ages[4], age_m := rep(x = ages, times = round(sum(age == ages[4])*c(0.01, 0.01, 0.08, 0.90)))]

## gender
pop_data[gender == "Male", gender_m := rep(x = c("Male", "Female"), times = sum(gender == "Male")*c(0.95, 0.05))]
pop_data[gender == "Female", gender_m := rep(x = c("Male", "Female"), times = sum(gender == "Female")*c(0.20, 0.80))]

## country
pop_data[country == "Belarus", country_m := rep(x = c("Belarus", "Ukraine"), times = sum(country == "Belarus")*c(0.90, 0.10))]
pop_data[country == "Ukraine", country_m := rep(x = c("Belarus", "Ukraine"), times = sum(country == "Ukraine")*c(0.05, 0.95))]

## stay
pop_data[short == 0, short_m := rep(x = 0:1, times = round(sum(short==0)*c(0.70, 0.30)))]
pop_data[short == 1, short_m := rep(x = 0:1, times = round(sum(short==1)*c(0.05, 0.95)))]

N <- nrow(pop_data)
n <- 1000
n_valid <- 1000
R <- 500
result_case <- list()

for (r in 1:R) {
  set.seed(23+r)
  print(r)
  
  if (r %% 25 == 0) {
    saveRDS(result_case, file = "results-whole-1000-list.rds")
  }
  
  sample_nonprob <- pop_data[, flag:=rbinom(nrow(pop_data), 1, prob)][flag == 1]
  
  ## validation sample
  ## imputation
  sample_nonprob[, id:=1:.N]
  sample_nonprob_simex <- copy(sample_nonprob)
  sample_nonprob_simex[, short_m:=as.factor(short_m)]
  sample_nonprob_simex[, age_m:=as.factor(age_m)]
  sample_nonprob_simex[, gender_m:=as.factor(gender_m)]
  sample_nonprob_simex[, country_m:=as.factor(country_m)]
  sample_nonprob_valid <- sample(1:nrow(sample_nonprob), n_valid)
  
  ## known
  short_mat <- with(pop_data, prop.table(table(short_m, short), margin = 1))
  age_mat <- with(pop_data, prop.table(table(age_m, age), margin = 1))
  gender_mat <- with(pop_data, prop.table(table(gender_m, gender), margin = 1))
  country_mat <- with(pop_data, prop.table(table(country_m, country), margin = 1))
  
  ## estimated error rates -- data is sparse ...
  # short_mat <- with(sample_nonprob[id %in% sample_nonprob_valid], prop.table(table(short_m, short), margin = 1))
  # age_mat <- with(sample_nonprob[id %in% sample_nonprob_valid], prop.table(table(age_m, age), margin = 1))
  # gender_mat <- with(sample_nonprob[id %in% sample_nonprob_valid], prop.table(table(gender_m, gender), margin = 1))
  # country_mat <- with(sample_nonprob[id %in% sample_nonprob_valid], prop.table(table(country_m, country), margin = 1))
  
  ## imputation of age variable only
  sample_nonprob_age <- copy(sample_nonprob)
  sample_nonprob_age[, ":="(short_t=short, age_t=age, gender_t=gender, country_t = country)]
  sample_nonprob_age[!id %in% sample_nonprob_valid, ":="(age=NA)]
  sample_nonprob_age <- sample_nonprob_age |> impute_shd(age ~ age_m + gender + country + short, k = 1) 
  
  ## imputation of all variables
  sample_nonprob_all <- copy(sample_nonprob)
  sample_nonprob_all[, ":="(short_t=short, age_t=age, gender_t=gender, country_t = country)]
  sample_nonprob_all[!id %in% sample_nonprob_valid, ":="(short=NA, age=NA, gender=NA, country = NA)]
  sample_nonprob_all <- sample_nonprob_all |> impute_shd(age + gender + country + short ~ age_m + gender_m + country_m + short_m, k = 1) 
  sample_nonprob_all[, short:=as.numeric(short)]
  
  sample_prob <- pop_data[sample(1:nrow(pop_data), size = n), ]
  sample_prob[, w:=N/n]
  ## in sample survey we observe it without error
  sample_prob[, age_m := age]
  sample_prob[, gender_m := gender]
  sample_prob[, country_m := country]
  sample_prob_svy <- svydesign(ids=~1, weights =~w, data = sample_prob)
  
  ############################################
  ## without error
  ipw_y1 <- tryCatch(expr = {nonprob(selection = ~ age + gender + country,
                                     target = ~ short,
                                     svydesign = sample_prob_svy,
                                     data = sample_nonprob)},
                     error = function(e) 'blad')
  
  if (class(ipw_y1)[1] !="nonprobsvy") {
    print("error")
    next
  }
  
  
  mi_y1 <- nonprob(outcome = short ~ age + gender + country,
                   svydesign = sample_prob_svy,
                   data = sample_nonprob,
                   method_outcome = "nn")
  
  dr_y1 <- nonprob(selection = ~  age + gender + country,
                   outcome = short ~ age + gender + country,
                   svydesign = sample_prob_svy,
                   data = sample_nonprob,
                   family_outcome = "binomial")
  
  ############################################
  
  ############################################
  ## age only
  ipw_y1_a <- nonprob(selection = ~ age,
                      target = ~ short,
                      svydesign = sample_prob_svy,
                      data = sample_nonprob)
  
  mi_y1_a <- nonprob(outcome = short ~ age,
                     svydesign = sample_prob_svy,
                     data = sample_nonprob,
                     method_outcome = "nn")
  
  dr_y1_a <- nonprob(selection = ~  age,
                     outcome = short ~ age,
                     svydesign = sample_prob_svy,
                     data = sample_nonprob,
                     family_outcome = "binomial")
  
  ############################################
  
  ############################################
  ## age only measurement -- no correction
  ipw_y1_a_m <- tryCatch(expr = {nonprob(selection = ~ age_m,
                                         target = ~ short_m,
                                         svydesign = sample_prob_svy,
                                         data = sample_nonprob)},
                         error = function(e) 'blad')
  
  if (class(ipw_y1_a_m)[1] !="nonprobsvy") {
    print("error")
    next
  }
  
  ## age only - correction
  ipw_y1_a_m_corr <- tryCatch(expr = {nonprob(selection = ~ age,
                                              target = ~ short,
                                              svydesign = sample_prob_svy,
                                              data = sample_nonprob_age)},
                              error = function(e) 'blad')
  
  ############################################
  ## age and short measurement -- no correction
  mi_y1_a_m <- nonprob(outcome = short_m ~ age_m,
                       svydesign = sample_prob_svy,
                       data = sample_nonprob,
                       method_outcome = "nn")
  
  ## age only measurement -- correction
  mi_y1_a_m_corr <- nonprob(outcome = short ~ age,
                            svydesign = sample_prob_svy,
                            data = sample_nonprob_age,
                            method_outcome = "nn")
  
  ## correction using simex
  mi_y1_a_m_simex <- mcsimex(model = glm(short_m ~ age_m, # + gender_m + country_m, 
                                         data = sample_nonprob_simex, 
                                         family = binomial, 
                                         x = TRUE),
                             mc.matrix = list(short_m = short_mat,
                                              age_m = age_mat),
                             #gender_m = gender_mat,
                             #country_m = country_mat),
                             SIMEXvariable = c("short_m", "age_m"), #, "gender_m", "country_m"),
                             jackknife.estimation=FALSE,
                             fitting.method = "quadratic",
                             asymptotic = F,
                             B = 1)
  
  sample_prob$short_pred <- predict(mi_y1_a_m_simex, newdata = sample_prob, type = "response")
  sample_nonprob$short_pred <- predict(mi_y1_a_m_simex, newdata = sample_nonprob, type = "response")
  inds <- nn2(data = sample_nonprob$short_pred, query = sample_prob$short_pred, k=1)
  mi_y1_a_m_simex_est_nn <- mean(as.numeric(as.character(sample_nonprob$short_m[inds$nn.idx[,1]])))
  mi_y1_a_m_simex_est_imp <- mean(sample_prob$short_pred)
  
  ## age only measurement -- no correction
  dr_y1_a_m <- tryCatch(expr = {nonprob(selection = ~  age_m,
                                        outcome = short_m ~ age_m,
                                        svydesign = sample_prob_svy,
                                        data = sample_nonprob,
                                        family_outcome = "binomial")},
                        error = function(e) 'blad'
  )
  
  dr_y1_a_m_corr <- tryCatch(expr = {nonprob(selection = ~  age,
                                             outcome = short ~ age,
                                             svydesign = sample_prob_svy,
                                             data = sample_nonprob_age,
                                             family_outcome = "binomial")},
                             error = function(e) 'blad'
  )
  ############################################
  ## in all
  ############################################
  ipw_y1_all <- tryCatch(expr = {nonprob(selection = ~ age_m + gender_m + country_m,
                                         target = ~ short_m,
                                         svydesign = sample_prob_svy,
                                         data = sample_nonprob)},
                         error = function(e) 'blad')
  if (class(ipw_y1_all)[1] !="nonprobsvy") {
    print("error")
    next
  }
  # correction
  ipw_y1_all_corr <- tryCatch(expr = {nonprob(selection = ~ age + gender + country,
                                              target = ~ short,
                                              svydesign = sample_prob_svy,
                                              data = sample_nonprob_all)},
                              error = function(e) 'blad')
  
  
  mi_y1_all <- nonprob(outcome = short_m ~ age_m + gender_m + country_m,
                       svydesign = sample_prob_svy,
                       data = sample_nonprob,
                       method_outcome = "nn")
  
  ## correction
  mi_y1_all_corr <- nonprob(outcome = short ~ age + gender + country,
                            svydesign = sample_prob_svy,
                            data = sample_nonprob_all,
                            method_outcome = "nn")
  
  ## correction using simex
  mi_y1_all_simex <- mcsimex(model = glm(short_m ~ age_m + gender_m + country_m, 
                                         data = sample_nonprob_simex, 
                                         family = binomial, 
                                         x = TRUE),
                             mc.matrix = list(short_m = short_mat,
                                              age_m = age_mat,
                                              gender_m = gender_mat,
                                              country_m = country_mat),
                             SIMEXvariable = c("short_m", "age_m", "gender_m", "country_m"),
                             jackknife.estimation=FALSE,
                             fitting.method = "quadratic",
                             asymptotic = F,
                             B = 1)
  
  sample_prob$short_pred <- predict(mi_y1_all_simex, newdata = sample_prob, type = "response")
  sample_nonprob$short_pred <- predict(mi_y1_all_simex, newdata = sample_nonprob, type = "response")
  inds <- nn2(data = sample_nonprob$short_pred, query = sample_prob$short_pred, k=1)
  mi_y1_all_simex_est_nn <- mean(as.numeric(as.character(sample_nonprob$short_m[inds$nn.idx[,1]])))
  mi_y1_all_simex_est_imp <- mean(sample_prob$short_pred)
  
  dr_y1_all <- tryCatch(expr = {nonprob(selection = ~  age_m + gender_m + country_m,
                                        outcome = short_m ~ age_m + gender_m + country_m,
                                        svydesign = sample_prob_svy,
                                        data = sample_nonprob,
                                        family_outcome = "binomial")},
                        error = function(e) 'blad')
  
  ## corr
  dr_y1_all_corr <- tryCatch(expr = {nonprob(selection = ~  age + gender + country,
                                             outcome = short ~ age + gender + country,
                                             svydesign = sample_prob_svy,
                                             data = sample_nonprob_all,
                                             family_outcome = "binomial")},
                             error = function(e) 'blad')
  
  
  result_case[[r]] <- data.table(iter= r,
                                 naive = c(mean(sample_nonprob$short)), 
                                 ipw_all_errno = ipw_y1$output$mean,
                                 mi_all_errno = mi_y1$output$mean,
                                 dr_all_errno = dr_y1$output$mean,
                                 ipw_age_errno = ipw_y1_a$output$mean,
                                 mi_age_errno = mi_y1_a$output$mean,
                                 dr_age_errno = dr_y1_a$output$mean,
                                 ipw_age_errage = ipw_y1_a_m$output$mean,
                                 ipw_age_errage_corr = ipw_y1_a_m_corr$output$mean,
                                 mi_age_errage = mi_y1_a_m$output$mean,
                                 mi_age_errage_corr = mi_y1_a_m_corr$output$mean,
                                 mi_age_errage_simexmi =mi_y1_a_m_simex_est_imp,
                                 mi_age_errage_simexnn =mi_y1_a_m_simex_est_nn,
                                 dr_age_errage = dr_y1_a_m$output$mean,
                                 dr_age_errage_corr = dr_y1_a_m_corr$output$mean,
                                 ipw_all_errall = ipw_y1_all$output$mean,
                                 ipw_all_errall_corr = ipw_y1_all_corr$output$mean,
                                 mi_all_errall = mi_y1_all$output$mean,
                                 mi_all_errall_corr = mi_y1_all_corr$output$mean,
                                 mi_all_errall_simexmi =mi_y1_all_simex_est_imp,
                                 mi_all_errall_simexnn = mi_y1_all_simex_est_nn,
                                 dr_all_errall = dr_y1_all$output$mean,
                                 dr_all_errall_corr = dr_y1_all_corr$output$mean)

}

result_case_df <- rbindlist(result_case)
saveRDS(result_case_df, file = "results-whole-1000-new.rds")
