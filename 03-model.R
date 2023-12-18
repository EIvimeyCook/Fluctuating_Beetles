
#load the source######
source("code/01-load_packages.R")
source("code/02-load_clean_data.R")

# Body Mass_Four #####
#four-way
m1a <- glmmTMB(beet_mass_mg ~
                 regime +
                 env +
                 pop +
                 sex +
                 regime:env +
                 regime:pop +
                 regime:sex +
                 env:pop +
                 env:sex +
                 pop:sex +
                 regime:env:pop +
                 regime:env:sex +
                 regime:pop:sex +
                 env:pop:sex +
                 regime:env:pop:sex +
                 (1 | jar),
               data = body_mass,
               REML = TRUE)

summary(m1a)
s1a <- simulateResiduals(m1a, plot = TRUE)
check_predictions(m1a)
#summary of four-way results

tab_model(m1a,
          transform = NULL,
          show.ci = FALSE,
          show.se = TRUE,
          show.r2 = FALSE,
          show.stat = TRUE,
          show.icc = FALSE)


#anova of full results
tidy(car::Anova(glmmTMB(beet_mass_mg ~
                          regime +
                          env +
                          pop +
                          sex +
                          regime:env +
                          regime:pop +
                          regime:sex +
                          env:pop +
                          env:sex +
                          pop:sex +
                          regime:env:pop +
                          regime:env:sex +
                          regime:pop:sex +
                          env:pop:sex +
                          regime:env:pop:sex +
                          (1 | jar),
                        data = body_mass,
                        family = "gaussian",
                        REML = TRUE, 
                        contrasts=list(regime = "contr.sum",
                                       env = "contr.sum",
                                       pop = "contr.sum",
                                       sex = "contr.sum")), type = "III")) %>%
  mutate_if(is.numeric, round, 5) %>%
  print()

#model means - no bias adjustment needed
em1a <- emmeans(m1a,
                specs = ~ sex * regime * env * pop)

confint(em1a, type = "response", calc = c(n = ~.wgt.), adjust = "mvt")

#model contrasts
summary(contrast(em1a, "pairwise",
                 combine = TRUE,
                 simple = list("pop","sex", "regime", "env", c("regime", "env")),
                 adjust = "mvt"),
        type = "response", infer = c(TRUE, TRUE)) %>%
  print()


#contrast of sex contrast
pairs(contrast(em1a, "consec",
                 combine = TRUE,
                 simple = list("sex"),
                 adjust = "mvt"), by = NULL) %>%
  print()

# DT-Four ######

#four-way
m1a <- glmmTMB(log(dt) ~
                 regime +
                 env +
                 pop +
                 sex +
                 regime:env +
                 regime:pop +
                 regime:sex +
                 env:pop +
                 env:sex +
                 pop:sex +
                 regime:env:pop +
                 regime:env:sex +
                 regime:pop:sex +
                 env:pop:sex +
                 regime:env:pop:sex +
                 (1 | jar / mat_id),
               data = devtime,
               family = "gaussian",
               REML = TRUE)

summary(m1a)
s1a <- simulateResiduals(m1a, plot = TRUE)

#summary of four-way results

tab_model(m1a,
          transform = NULL,
          show.ci = FALSE,
          show.se = TRUE,
          show.r2 = FALSE,
          show.stat = TRUE,
          show.icc = FALSE)


#anova of full results
tidy(car::Anova(glmmTMB(log(dt) ~
                          regime +
                          env +
                          pop +
                          sex +
                          regime:env +
                          regime:pop +
                          regime:sex +
                          env:pop +
                          env:sex +
                          pop:sex +
                          regime:env:pop +
                          regime:env:sex +
                          regime:pop:sex +
                          env:pop:sex +
                          regime:env:pop:sex +
                          (1 | jar / mat_id),
                        data = devtime,
                        family = "gaussian",
                        REML = TRUE, 
                        contrasts=list(regime = "contr.sum",
                                       env = "contr.sum",
                                       pop = "contr.sum",
                                       sex = "contr.sum")), type = "III")) %>%
  mutate_if(is.numeric, round, 5) %>%
  print()

#model means - no bias adjustment needed
em1a <- emmeans(m1a,
                specs = ~ sex * regime * env * pop)

confint(em1a, calc = c(n = ~.wgt.), adjust = "mvt",
        type = "response")

#model contrasts
summary(contrast(regrid(em1a), "pairwise",
                 combine = TRUE,
                 simple = list("pop","sex", "regime", "env",c("regime", "env")),
                 adjust = "mvt"),
        type = "response", infer = c(TRUE, TRUE)) %>%
  print()


#contrast of sex contrast
pairs(contrast(em1a, "consec",
               combine = TRUE,
               simple = list("sex"),
               adjust = "mvt"), by = NULL) %>%
  print()

# LRS_Three  #######

#create list for LRS models
list_lrs <- list()

#first model (poisson), detect overdispersion and zi
list_lrs[[1]] <- glmmTMB(lrs ~
                           regime +
                           env +
                           pop +
                           regime:env +
                           regime:pop +
                           env:pop +
                           regime:env:pop +
                           (1 | jar),
                           data = repro_wide,
                           family = "poisson")

#check
s1a <- simulateResiduals(list_lrs[[1]], plot = TRUE)
testZeroInflation(s1a)
testDispersion(s1a)

#second model (poisson + olre), detect overdispersion and zi
list_lrs[[2]] <- glmmTMB(lrs ~
                           regime +
                           env +
                           pop +
                           regime:env +
                           regime:pop +
                           env:pop +
                           regime:env:pop +
                           (1 | jar) +
                           (1 | obs),
                         data = repro_wide,
                         family = "poisson")

summary(list_lrs[[2]])

#check
s1a <- simulateResiduals(list_lrs[[2]], plot = TRUE)
testZeroInflation(s1a)
testDispersion(s1a)

#third model (nbinom), detect overdispersion and zi
list_lrs[[3]] <- glmmTMB(lrs ~
                           regime +
                           env +
                           pop +
                           regime:env +
                           regime:pop +
                           env:pop +
                           regime:env:pop +
                           (1 | jar),
                         data = repro_wide,
                         family = "nbinom2")

summary(list_lrs[[3]])

#check
s1a <- simulateResiduals(list_lrs[[3]], plot = TRUE)
testZeroInflation(s1a)
testDispersion(s1a)

#modelchecker code
modelchecker(list_lrs)

#fit with zero-inflated parameter - zi detected

list_lrs[[4]] <-glmmTMB(lrs ~ regime + env + pop + regime:env + regime:pop + env:pop + regime:env:pop + (1 | jar), zi = ~pop, family = poisson, data = repro_wide)
list_lrs[[5]] <-glmmTMB(lrs ~ regime + env + pop + regime:env + regime:pop + env:pop + regime:env:pop + (1 | jar), zi = ~env, family = poisson, data = repro_wide)
list_lrs[[6]] <-glmmTMB(lrs ~ regime + env + pop + regime:env + regime:pop + env:pop + regime:env:pop + (1 | jar), zi = ~regime, family = poisson, data = repro_wide)
list_lrs[[7]] <-glmmTMB(lrs ~ regime + env + pop + regime:env + regime:pop + env:pop + regime:env:pop + (1 | jar), zi = ~pop + env, family = poisson, data = repro_wide)
list_lrs[[8]] <-glmmTMB(lrs ~ regime + env + pop + regime:env + regime:pop + env:pop + regime:env:pop + (1 | jar), zi = ~pop + regime, family = poisson, data = repro_wide)
list_lrs[[9]] <-glmmTMB(lrs ~ regime + env + pop + regime:env + regime:pop + env:pop + regime:env:pop + (1 | jar), zi = ~env + regime, family = poisson, data = repro_wide)
list_lrs[[10]] <-glmmTMB(lrs ~ regime + env + pop + regime:env + regime:pop + env:pop + regime:env:pop + (1 | jar), zi = ~pop + env + regime, family = poisson, data = repro_wide)

list_lrs[[11]] <-glmmTMB(lrs ~ regime + env + pop + regime:env + regime:pop + env:pop + regime:env:pop + (1 | jar), zi = ~pop, family = nbinom2, data = repro_wide)
list_lrs[[12]] <-glmmTMB(lrs ~ regime + env + pop + regime:env + regime:pop + env:pop + regime:env:pop + (1 | jar), zi = ~env, family = nbinom2, data = repro_wide)
list_lrs[[13]] <-glmmTMB(lrs ~ regime + env + pop + regime:env + regime:pop + env:pop + regime:env:pop + (1 | jar), zi = ~regime, family = nbinom2, data = repro_wide)
list_lrs[[14]] <-glmmTMB(lrs ~ regime + env + pop + regime:env + regime:pop + env:pop + regime:env:pop + (1 | jar), zi = ~pop + env, family = nbinom2, data = repro_wide)
list_lrs[[15]] <-glmmTMB(lrs ~ regime + env + pop + regime:env + regime:pop + env:pop + regime:env:pop + (1 | jar), zi = ~pop + regime, family = nbinom2, data = repro_wide)
list_lrs[[16]] <-glmmTMB(lrs ~ regime + env + pop + regime:env + regime:pop + env:pop + regime:env:pop + (1 | jar), zi = ~env + regime, family = nbinom2, data = repro_wide)
list_lrs[[17]] <-glmmTMB(lrs ~ regime + env + pop + regime:env + regime:pop + env:pop + regime:env:pop + (1 | jar), zi = ~pop + env + regime, family = nbinom2, data = repro_wide)

list_lrs[[18]] <-glmmTMB(lrs ~ regime + env + pop + regime:env + regime:pop + env:pop + regime:env:pop + (1 | jar), zi = ~pop, family = genpois, data = repro_wide)
list_lrs[[19]] <-glmmTMB(lrs ~ regime + env + pop + regime:env + regime:pop + env:pop + regime:env:pop + (1 | jar), zi = ~env, family = genpois, data = repro_wide)
list_lrs[[20]] <-glmmTMB(lrs ~ regime + env + pop + regime:env + regime:pop + env:pop + regime:env:pop + (1 | jar), zi = ~regime, family = genpois, data = repro_wide)
list_lrs[[21]] <-glmmTMB(lrs ~ regime + env + pop + regime:env + regime:pop + env:pop + regime:env:pop + (1 | jar), zi = ~pop + env, family = genpois, data = repro_wide)
list_lrs[[22]] <-glmmTMB(lrs ~ regime + env + pop + regime:env + regime:pop + env:pop + regime:env:pop + (1 | jar), zi = ~pop + regime, family = genpois, data = repro_wide)
list_lrs[[23]] <-glmmTMB(lrs ~ regime + env + pop + regime:env + regime:pop + env:pop + regime:env:pop + (1 | jar), zi = ~env + regime, family = genpois, data = repro_wide)
list_lrs[[24]] <-glmmTMB(lrs ~ regime + env + pop + regime:env + regime:pop + env:pop + regime:env:pop + (1 | jar), zi = ~pop + env + regime, family = genpois, data = repro_wide)

list_lrs[[25]] <-glmmTMB(lrs ~ regime + env + pop + regime:env + regime:pop + env:pop + regime:env:pop + (1 | jar), zi = ~pop, family = compois, data = repro_wide)
list_lrs[[26]] <-glmmTMB(lrs ~ regime + env + pop + regime:env + regime:pop + env:pop + regime:env:pop + (1 | jar), zi = ~env, family = compois, data = repro_wide)
list_lrs[[27]] <-glmmTMB(lrs ~ regime + env + pop + regime:env + regime:pop + env:pop + regime:env:pop + (1 | jar), zi = ~regime, family = compois, data = repro_wide)
list_lrs[[28]] <-glmmTMB(lrs ~ regime + env + pop + regime:env + regime:pop + env:pop + regime:env:pop + (1 | jar), zi = ~pop + env, family = compois, data = repro_wide)
list_lrs[[29]] <-glmmTMB(lrs ~ regime + env + pop + regime:env + regime:pop + env:pop + regime:env:pop + (1 | jar), zi = ~pop + regime, family = compois, data = repro_wide)
list_lrs[[30]] <-glmmTMB(lrs ~ regime + env + pop + regime:env + regime:pop + env:pop + regime:env:pop + (1 | jar), zi = ~env + regime, family = compois, data = repro_wide)
list_lrs[[31]] <-glmmTMB(lrs ~ regime + env + pop + regime:env + regime:pop + env:pop + regime:env:pop + (1 | jar), zi = ~pop + env + regime, family = compois, data = repro_wide)

list_lrs[[32]] <-glmmTMB(lrs ~ regime + env + pop + regime:env + regime:pop + env:pop + regime:env:pop + (1 | jar), zi = ~1, family = poisson, data = repro_wide)
list_lrs[[33]] <-glmmTMB(lrs ~ regime + env + pop + regime:env + regime:pop + env:pop + regime:env:pop + (1 | jar), zi = ~1, family = nbinom2, data = repro_wide)
list_lrs[[34]] <-glmmTMB(lrs ~ regime + env + pop + regime:env + regime:pop + env:pop + regime:env:pop + (1 | jar), zi = ~1, family = genpois, data = repro_wide)
list_lrs[[35]] <-glmmTMB(lrs ~ regime + env + pop + regime:env + regime:pop + env:pop + regime:env:pop + (1 | jar), zi = ~1, family = compois, data = repro_wide)

list_lrs[[36]] <-glmmTMB(lrs ~ regime + env + pop + regime:env + regime:pop + env:pop + regime:env:pop + (1 | jar), zi = ~0, family = poisson, data = repro_wide)
list_lrs[[37]] <-glmmTMB(lrs ~ regime + env + pop + regime:env + regime:pop + env:pop + regime:env:pop + (1 | jar), zi = ~0, family = nbinom2, data = repro_wide)
list_lrs[[38]] <-glmmTMB(lrs ~ regime + env + pop + regime:env + regime:pop + env:pop + regime:env:pop + (1 | jar), zi = ~0, family = genpois, data = repro_wide)
list_lrs[[39]] <-glmmTMB(lrs ~ regime + env + pop + regime:env + regime:pop + env:pop + regime:env:pop + (1 | jar), zi = ~0, family = compois, data = repro_wide)

#model ranking via aic
model.sel(list_lrs, rank = "AIC")
modelchecker(list_lrs)

#best model is compois w/ env zi
summary(list_lrs[[26]])

#complete seperation, add as a random effect
list_lrs_44 <- glmmTMB(lrs ~
                         regime +
                         env +
                         pop +
                         regime:env +
                         regime:pop +
                         env:pop +
                         regime:env:pop +
                         (1 | jar),
                       data = repro_wide,
                       zi = ~(1 | env),
                       family = poisson)

summary(list_lrs_44)

#simulate residuals
s1a <- simulateResiduals(list_lrs_44)
plot(s1a)
testZeroInflation(s1a)

#slight heteroskedastic residuals
DHARMa::testCategorical(s1a,
                        repro_wide$pop:repro_wide$regime:repro_wide$env)

#change dispersion
list_lrs_disp <- list()

list_lrs_disp[[1]] <-glmmTMB(lrs ~ regime + env + pop + regime:env + regime:pop + env:pop + regime:env:pop + (1 | jar), zi = ~(1 | env), disp = ~pop, family = compois, data = repro_wide)
list_lrs_disp[[2]] <-glmmTMB(lrs ~ regime + env + pop + regime:env + regime:pop + env:pop + regime:env:pop + (1 | jar), zi = ~(1 | env), disp = ~env, family = compois, data = repro_wide)
list_lrs_disp[[3]] <-glmmTMB(lrs ~ regime + env + pop + regime:env + regime:pop + env:pop + regime:env:pop + (1 | jar), zi = ~(1 | env), disp = ~regime, family = compois, data = repro_wide)
list_lrs_disp[[4]] <-glmmTMB(lrs ~ regime + env + pop + regime:env + regime:pop + env:pop + regime:env:pop + (1 | jar), zi = ~(1 | env), disp = ~pop + env, family = compois, data = repro_wide)
list_lrs_disp[[5]] <-glmmTMB(lrs ~ regime + env + pop + regime:env + regime:pop + env:pop + regime:env:pop + (1 | jar), zi = ~(1 | env), disp = ~pop + regime, family = compois, data = repro_wide)
list_lrs_disp[[6]] <-glmmTMB(lrs ~ regime + env + pop + regime:env + regime:pop + env:pop + regime:env:pop + (1 | jar), zi = ~(1 | env), disp = ~env + regime, family = compois, data = repro_wide)
list_lrs_disp[[7]] <-glmmTMB(lrs ~ regime + env + pop + regime:env + regime:pop + env:pop + regime:env:pop + (1 | jar), zi = ~(1 | env), disp = ~pop + env + regime, family = compois, data = repro_wide)

#model ranking via aic
model.sel(list_lrs_disp, rank = "AIC")
modelchecker(list_lrs_disp)

#best model is compois w/ env zi
summary(list_lrs_disp[[7]])

#simulate residuals
s1a <- simulateResiduals(list_lrs_disp[[7]])
plot(s1a)
testZeroInflation(s1a)
check_predictions(list_lrs_disp[[7]])
check_convergence(list_lrs_disp[[7]])

#summary
tab_model(list_lrs_disp[[7]],
          transform = NULL,
          show.ci = FALSE,
          show.se = TRUE,
          show.r2 = FALSE,
          show.stat = TRUE,
          show.icc = FALSE)

#anova of full results
tidy(car::Anova(glmmTMB(lrs ~
                          regime +
                          env +
                          pop +
                          regime:env +
                          regime:pop +
                          env:pop +
                          regime:env:pop +
                          (1 | jar),
                        data = repro_wide,
                        zi = ~(1 | env),
                        disp = ~pop + env + regime,
                        family = compois, 
                        contrasts=list(regime = "contr.sum",
                                       env = "contr.sum",
                                       pop = "contr.sum")), type = "III")) %>%
  mutate_if(is.numeric, round, 5) %>%
  print()

#model means + bias ajdutmsent
lme4::VarCorr(list_lrs_disp[[7]])
sigma <- sqrt(0.048^2)

em1a <- emmeans(list_lrs_disp[[7]],
                specs = ~ regime * env * pop,
                bias.adjust = T,
                sigma = sigma)

confint(em1a, calc = c(n = ~.wgt.), adjust = "mvt",
        type = "response")

#model contrasts
summary(contrast(regrid(em1a), "pairwise",
                 combine = TRUE,
                 simple = list("pop", "regime", "env", c("regime", "env")),
                 adjust = "mvt"),
        type = "response", infer = c(TRUE, TRUE))

#save
summary(contrast(regrid(em1a), "pairwise",
                 combine = TRUE,
                 simple = list("pop", "regime", "env", c("regime", "env")),
                 adjust = "mvt"),
        type = "response", infer = c(TRUE, TRUE)) %>%
  print()

# ASR_Four########

#create list for asr models
list_asr <- list()

#first model (poisson), detect overdispersion and zi
list_asr[[1]] <- glmmTMB(offspring ~
                           regime +
                           env +
                           pop +
                           day +
                           regime:env +
                           regime:pop +
                           env:pop +
                           regime:day +
                           env:day +
                           pop:day +
                           regime:env:pop +
                           regime:env:day +
                           regime:pop:day +
                           env:pop:day +
                           regime:env:pop:day +
                           (1 | jar / mat_id),
                         data = repro_long,
                         family = "poisson")

#check
s1a <- simulateResiduals(list_asr[[1]], plot = TRUE)
testZeroInflation(s1a)
testDispersion(s1a)

#second model (poisson + olre), detect overdispersion and zi
list_asr[[2]] <- glmmTMB(offspring ~
                           regime +
                           env +
                           pop +
                           day +
                           regime:env +
                           regime:pop +
                           env:pop +
                           regime:day +
                           env:day +
                           pop:day +
                           regime:env:pop +
                           regime:env:day +
                           regime:pop:day +
                           env:pop:day +
                           regime:env:pop:day +
                           (1 | jar / mat_id) +
                           (1 | obs),
                         data = repro_long,
                         family = "poisson")

summary(list_asr[[2]])

#check
s1a <- simulateResiduals(list_asr[[2]], plot = TRUE)
testZeroInflation(s1a)
testDispersion(s1a)

#third model (nbinom), detect overdispersion and zi
list_asr[[3]] <- glmmTMB(offspring ~
                           regime +
                           env +
                           pop +
                           day +
                           regime:env +
                           regime:pop +
                           env:pop +
                           regime:day +
                           env:day +
                           pop:day +
                           regime:env:pop +
                           regime:env:day +
                           regime:pop:day +
                           env:pop:day +
                           regime:env:pop:day +
                           (1 | jar / mat_id),
                         data = repro_long,
                         family = "nbinom2")

summary(list_asr[[3]])

#check
s1a <- simulateResiduals(list_asr[[3]], plot = TRUE)
testZeroInflation(s1a)
testDispersion(s1a)

#modelchecker code
modelchecker(list_asr)

#model selection on zeroinflation parameters

list_asr[[4]] <-glmmTMB(offspring ~ regime + env + pop + day + regime:env + regime:pop + env:pop + regime:day + env:day + pop:day + regime:env:pop + regime:env:day + regime:pop:day + env:pop:day + regime:env:pop:day + (1 | jar/mat_id), zi = ~pop, family = poisson, data = repro_long)
list_asr[[5]] <-glmmTMB(offspring ~ regime + env + pop + day + regime:env + regime:pop + env:pop + regime:day + env:day + pop:day + regime:env:pop + regime:env:day + regime:pop:day + env:pop:day + regime:env:pop:day + (1 | jar/mat_id), zi = ~env, family = poisson, data = repro_long)
list_asr[[6]] <-glmmTMB(offspring ~ regime + env + pop + day + regime:env + regime:pop + env:pop + regime:day + env:day + pop:day + regime:env:pop + regime:env:day + regime:pop:day + env:pop:day + regime:env:pop:day + (1 | jar/mat_id), zi = ~regime, family = poisson, data = repro_long)
list_asr[[7]] <-glmmTMB(offspring ~ regime + env + pop + day + regime:env + regime:pop + env:pop + regime:day + env:day + pop:day + regime:env:pop + regime:env:day + regime:pop:day + env:pop:day + regime:env:pop:day + (1 | jar/mat_id), zi = ~day, family = poisson, data = repro_long)
list_asr[[8]] <-glmmTMB(offspring ~ regime + env + pop + day + regime:env + regime:pop + env:pop + regime:day + env:day + pop:day + regime:env:pop + regime:env:day + regime:pop:day + env:pop:day + regime:env:pop:day + (1 | jar/mat_id), zi = ~pop + env, family = poisson, data = repro_long)
list_asr[[9]] <-glmmTMB(offspring ~ regime + env + pop + day + regime:env + regime:pop + env:pop + regime:day + env:day + pop:day + regime:env:pop + regime:env:day + regime:pop:day + env:pop:day + regime:env:pop:day + (1 | jar/mat_id), zi = ~pop + regime, family = poisson, data = repro_long)
list_asr[[10]] <-glmmTMB(offspring ~ regime + env + pop + day + regime:env + regime:pop + env:pop + regime:day + env:day + pop:day + regime:env:pop + regime:env:day + regime:pop:day + env:pop:day + regime:env:pop:day + (1 | jar/mat_id), zi = ~pop + day, family = poisson, data = repro_long)
list_asr[[11]] <-glmmTMB(offspring ~ regime + env + pop + day + regime:env + regime:pop + env:pop + regime:day + env:day + pop:day + regime:env:pop + regime:env:day + regime:pop:day + env:pop:day + regime:env:pop:day + (1 | jar/mat_id), zi = ~env + regime, family = poisson, data = repro_long)
list_asr[[12]] <-glmmTMB(offspring ~ regime + env + pop + day + regime:env + regime:pop + env:pop + regime:day + env:day + pop:day + regime:env:pop + regime:env:day + regime:pop:day + env:pop:day + regime:env:pop:day + (1 | jar/mat_id), zi = ~env + day, family = poisson, data = repro_long)
list_asr[[13]] <-glmmTMB(offspring ~ regime + env + pop + day + regime:env + regime:pop + env:pop + regime:day + env:day + pop:day + regime:env:pop + regime:env:day + regime:pop:day + env:pop:day + regime:env:pop:day + (1 | jar/mat_id), zi = ~regime + day, family = poisson, data = repro_long)
list_asr[[14]] <-glmmTMB(offspring ~ regime + env + pop + day + regime:env + regime:pop + env:pop + regime:day + env:day + pop:day + regime:env:pop + regime:env:day + regime:pop:day + env:pop:day + regime:env:pop:day + (1 | jar/mat_id), zi = ~pop + env + regime, family = poisson, data = repro_long)
list_asr[[15]] <-glmmTMB(offspring ~ regime + env + pop + day + regime:env + regime:pop + env:pop + regime:day + env:day + pop:day + regime:env:pop + regime:env:day + regime:pop:day + env:pop:day + regime:env:pop:day + (1 | jar/mat_id), zi = ~pop + env + day, family = poisson, data = repro_long)
list_asr[[16]] <-glmmTMB(offspring ~ regime + env + pop + day + regime:env + regime:pop + env:pop + regime:day + env:day + pop:day + regime:env:pop + regime:env:day + regime:pop:day + env:pop:day + regime:env:pop:day + (1 | jar/mat_id), zi = ~pop + regime + day, family = poisson, data = repro_long)
list_asr[[17]] <-glmmTMB(offspring ~ regime + env + pop + day + regime:env + regime:pop + env:pop + regime:day + env:day + pop:day + regime:env:pop + regime:env:day + regime:pop:day + env:pop:day + regime:env:pop:day + (1 | jar/mat_id), zi = ~env + regime + day, family = poisson, data = repro_long)
list_asr[[18]] <-glmmTMB(offspring ~ regime + env + pop + day + regime:env + regime:pop + env:pop + regime:day + env:day + pop:day + regime:env:pop + regime:env:day + regime:pop:day + env:pop:day + regime:env:pop:day + (1 | jar/mat_id), zi = ~pop + env + regime + day, family = poisson, data = repro_long)

list_asr[[19]] <-glmmTMB(offspring ~ regime + env + pop + day + regime:env + regime:pop + env:pop + regime:day + env:day + pop:day + regime:env:pop + regime:env:day + regime:pop:day + env:pop:day + regime:env:pop:day + (1 | jar/mat_id), zi = ~pop, family = nbinom2, data = repro_long)
list_asr[[20]] <-glmmTMB(offspring ~ regime + env + pop + day + regime:env + regime:pop + env:pop + regime:day + env:day + pop:day + regime:env:pop + regime:env:day + regime:pop:day + env:pop:day + regime:env:pop:day + (1 | jar/mat_id), zi = ~env, family = nbinom2, data = repro_long)
list_asr[[21]] <-glmmTMB(offspring ~ regime + env + pop + day + regime:env + regime:pop + env:pop + regime:day + env:day + pop:day + regime:env:pop + regime:env:day + regime:pop:day + env:pop:day + regime:env:pop:day + (1 | jar/mat_id), zi = ~regime, family = nbinom2, data = repro_long)
list_asr[[22]] <-glmmTMB(offspring ~ regime + env + pop + day + regime:env + regime:pop + env:pop + regime:day + env:day + pop:day + regime:env:pop + regime:env:day + regime:pop:day + env:pop:day + regime:env:pop:day + (1 | jar/mat_id), zi = ~day, family = nbinom2, data = repro_long)
list_asr[[23]] <-glmmTMB(offspring ~ regime + env + pop + day + regime:env + regime:pop + env:pop + regime:day + env:day + pop:day + regime:env:pop + regime:env:day + regime:pop:day + env:pop:day + regime:env:pop:day + (1 | jar/mat_id), zi = ~pop + env, family = nbinom2, data = repro_long)
list_asr[[24]] <-glmmTMB(offspring ~ regime + env + pop + day + regime:env + regime:pop + env:pop + regime:day + env:day + pop:day + regime:env:pop + regime:env:day + regime:pop:day + env:pop:day + regime:env:pop:day + (1 | jar/mat_id), zi = ~pop + regime, family = nbinom2, data = repro_long)
list_asr[[25]] <-glmmTMB(offspring ~ regime + env + pop + day + regime:env + regime:pop + env:pop + regime:day + env:day + pop:day + regime:env:pop + regime:env:day + regime:pop:day + env:pop:day + regime:env:pop:day + (1 | jar/mat_id), zi = ~pop + day, family = nbinom2, data = repro_long)
list_asr[[26]] <-glmmTMB(offspring ~ regime + env + pop + day + regime:env + regime:pop + env:pop + regime:day + env:day + pop:day + regime:env:pop + regime:env:day + regime:pop:day + env:pop:day + regime:env:pop:day + (1 | jar/mat_id), zi = ~env + regime, family = nbinom2, data = repro_long)
list_asr[[27]] <-glmmTMB(offspring ~ regime + env + pop + day + regime:env + regime:pop + env:pop + regime:day + env:day + pop:day + regime:env:pop + regime:env:day + regime:pop:day + env:pop:day + regime:env:pop:day + (1 | jar/mat_id), zi = ~env + day, family = nbinom2, data = repro_long)
list_asr[[28]] <-glmmTMB(offspring ~ regime + env + pop + day + regime:env + regime:pop + env:pop + regime:day + env:day + pop:day + regime:env:pop + regime:env:day + regime:pop:day + env:pop:day + regime:env:pop:day + (1 | jar/mat_id), zi = ~regime + day, family = nbinom2, data = repro_long)
list_asr[[29]] <-glmmTMB(offspring ~ regime + env + pop + day + regime:env + regime:pop + env:pop + regime:day + env:day + pop:day + regime:env:pop + regime:env:day + regime:pop:day + env:pop:day + regime:env:pop:day + (1 | jar/mat_id), zi = ~pop + env + regime, family = nbinom2, data = repro_long)
list_asr[[30]] <-glmmTMB(offspring ~ regime + env + pop + day + regime:env + regime:pop + env:pop + regime:day + env:day + pop:day + regime:env:pop + regime:env:day + regime:pop:day + env:pop:day + regime:env:pop:day + (1 | jar/mat_id), zi = ~pop + env + day, family = nbinom2, data = repro_long)
list_asr[[31]] <-glmmTMB(offspring ~ regime + env + pop + day + regime:env + regime:pop + env:pop + regime:day + env:day + pop:day + regime:env:pop + regime:env:day + regime:pop:day + env:pop:day + regime:env:pop:day + (1 | jar/mat_id), zi = ~pop + regime + day, family = nbinom2, data = repro_long)
list_asr[[32]] <-glmmTMB(offspring ~ regime + env + pop + day + regime:env + regime:pop + env:pop + regime:day + env:day + pop:day + regime:env:pop + regime:env:day + regime:pop:day + env:pop:day + regime:env:pop:day + (1 | jar/mat_id), zi = ~env + regime + day, family = nbinom2, data = repro_long)
list_asr[[33]] <-glmmTMB(offspring ~ regime + env + pop + day + regime:env + regime:pop + env:pop + regime:day + env:day + pop:day + regime:env:pop + regime:env:day + regime:pop:day + env:pop:day + regime:env:pop:day + (1 | jar/mat_id), zi = ~pop + env + regime + day, family = nbinom2, data = repro_long)

list_asr[[34]] <-glmmTMB(offspring ~ regime + env + pop + day + regime:env + regime:pop + env:pop + regime:day + env:day + pop:day + regime:env:pop + regime:env:day + regime:pop:day + env:pop:day + regime:env:pop:day + (1 | jar/mat_id), zi = ~pop, family = genpois, data = repro_long)
list_asr[[35]] <-glmmTMB(offspring ~ regime + env + pop + day + regime:env + regime:pop + env:pop + regime:day + env:day + pop:day + regime:env:pop + regime:env:day + regime:pop:day + env:pop:day + regime:env:pop:day + (1 | jar/mat_id), zi = ~env, family = genpois, data = repro_long)
list_asr[[36]] <-glmmTMB(offspring ~ regime + env + pop + day + regime:env + regime:pop + env:pop + regime:day + env:day + pop:day + regime:env:pop + regime:env:day + regime:pop:day + env:pop:day + regime:env:pop:day + (1 | jar/mat_id), zi = ~regime, family = genpois, data = repro_long)
list_asr[[37]] <-glmmTMB(offspring ~ regime + env + pop + day + regime:env + regime:pop + env:pop + regime:day + env:day + pop:day + regime:env:pop + regime:env:day + regime:pop:day + env:pop:day + regime:env:pop:day + (1 | jar/mat_id), zi = ~day, family = genpois, data = repro_long)
list_asr[[38]] <-glmmTMB(offspring ~ regime + env + pop + day + regime:env + regime:pop + env:pop + regime:day + env:day + pop:day + regime:env:pop + regime:env:day + regime:pop:day + env:pop:day + regime:env:pop:day + (1 | jar/mat_id), zi = ~pop + env, family = genpois, data = repro_long)
list_asr[[39]] <-glmmTMB(offspring ~ regime + env + pop + day + regime:env + regime:pop + env:pop + regime:day + env:day + pop:day + regime:env:pop + regime:env:day + regime:pop:day + env:pop:day + regime:env:pop:day + (1 | jar/mat_id), zi = ~pop + regime, family = genpois, data = repro_long)
list_asr[[40]] <-glmmTMB(offspring ~ regime + env + pop + day + regime:env + regime:pop + env:pop + regime:day + env:day + pop:day + regime:env:pop + regime:env:day + regime:pop:day + env:pop:day + regime:env:pop:day + (1 | jar/mat_id), zi = ~pop + day, family = genpois, data = repro_long)
list_asr[[41]] <-glmmTMB(offspring ~ regime + env + pop + day + regime:env + regime:pop + env:pop + regime:day + env:day + pop:day + regime:env:pop + regime:env:day + regime:pop:day + env:pop:day + regime:env:pop:day + (1 | jar/mat_id), zi = ~env + regime, family = genpois, data = repro_long)
list_asr[[42]] <-glmmTMB(offspring ~ regime + env + pop + day + regime:env + regime:pop + env:pop + regime:day + env:day + pop:day + regime:env:pop + regime:env:day + regime:pop:day + env:pop:day + regime:env:pop:day + (1 | jar/mat_id), zi = ~env + day, family = genpois, data = repro_long)
list_asr[[43]] <-glmmTMB(offspring ~ regime + env + pop + day + regime:env + regime:pop + env:pop + regime:day + env:day + pop:day + regime:env:pop + regime:env:day + regime:pop:day + env:pop:day + regime:env:pop:day + (1 | jar/mat_id), zi = ~regime + day, family = genpois, data = repro_long)
list_asr[[44]] <-glmmTMB(offspring ~ regime + env + pop + day + regime:env + regime:pop + env:pop + regime:day + env:day + pop:day + regime:env:pop + regime:env:day + regime:pop:day + env:pop:day + regime:env:pop:day + (1 | jar/mat_id), zi = ~pop + env + regime, family = genpois, data = repro_long)
list_asr[[45]] <-glmmTMB(offspring ~ regime + env + pop + day + regime:env + regime:pop + env:pop + regime:day + env:day + pop:day + regime:env:pop + regime:env:day + regime:pop:day + env:pop:day + regime:env:pop:day + (1 | jar/mat_id), zi = ~pop + env + day, family = genpois, data = repro_long)
list_asr[[46]] <-glmmTMB(offspring ~ regime + env + pop + day + regime:env + regime:pop + env:pop + regime:day + env:day + pop:day + regime:env:pop + regime:env:day + regime:pop:day + env:pop:day + regime:env:pop:day + (1 | jar/mat_id), zi = ~env + regime + day, family = genpois, data = repro_long)
list_asr[[47]] <-glmmTMB(offspring ~ regime + env + pop + day + regime:env + regime:pop + env:pop + regime:day + env:day + pop:day + regime:env:pop + regime:env:day + regime:pop:day + env:pop:day + regime:env:pop:day + (1 | jar/mat_id), zi = ~pop + env + regime + day, family = genpois, data = repro_long)

list_asr[[48]] <-glmmTMB(offspring ~ regime + env + pop + day + regime:env + regime:pop + env:pop + regime:day + env:day + pop:day + regime:env:pop + regime:env:day + regime:pop:day + env:pop:day + regime:env:pop:day + (1 | jar/mat_id), zi = ~pop, family = compois, data = repro_long)
list_asr[[49]] <-glmmTMB(offspring ~ regime + env + pop + day + regime:env + regime:pop + env:pop + regime:day + env:day + pop:day + regime:env:pop + regime:env:day + regime:pop:day + env:pop:day + regime:env:pop:day + (1 | jar/mat_id), zi = ~env, family = compois, data = repro_long)
list_asr[[50]] <-glmmTMB(offspring ~ regime + env + pop + day + regime:env + regime:pop + env:pop + regime:day + env:day + pop:day + regime:env:pop + regime:env:day + regime:pop:day + env:pop:day + regime:env:pop:day + (1 | jar/mat_id), zi = ~regime, family = compois, data = repro_long)
list_asr[[51]] <-glmmTMB(offspring ~ regime + env + pop + day + regime:env + regime:pop + env:pop + regime:day + env:day + pop:day + regime:env:pop + regime:env:day + regime:pop:day + env:pop:day + regime:env:pop:day + (1 | jar/mat_id), zi = ~day, family = compois, data = repro_long)
list_asr[[52]] <-glmmTMB(offspring ~ regime + env + pop + day + regime:env + regime:pop + env:pop + regime:day + env:day + pop:day + regime:env:pop + regime:env:day + regime:pop:day + env:pop:day + regime:env:pop:day + (1 | jar/mat_id), zi = ~pop + env, family = compois, data = repro_long)
list_asr[[53]] <-glmmTMB(offspring ~ regime + env + pop + day + regime:env + regime:pop + env:pop + regime:day + env:day + pop:day + regime:env:pop + regime:env:day + regime:pop:day + env:pop:day + regime:env:pop:day + (1 | jar/mat_id), zi = ~pop + regime, family = compois, data = repro_long)
list_asr[[54]] <-glmmTMB(offspring ~ regime + env + pop + day + regime:env + regime:pop + env:pop + regime:day + env:day + pop:day + regime:env:pop + regime:env:day + regime:pop:day + env:pop:day + regime:env:pop:day + (1 | jar/mat_id), zi = ~pop + day, family = compois, data = repro_long)
list_asr[[55]] <-glmmTMB(offspring ~ regime + env + pop + day + regime:env + regime:pop + env:pop + regime:day + env:day + pop:day + regime:env:pop + regime:env:day + regime:pop:day + env:pop:day + regime:env:pop:day + (1 | jar/mat_id), zi = ~env + regime, family = compois, data = repro_long)
list_asr[[56]] <-glmmTMB(offspring ~ regime + env + pop + day + regime:env + regime:pop + env:pop + regime:day + env:day + pop:day + regime:env:pop + regime:env:day + regime:pop:day + env:pop:day + regime:env:pop:day + (1 | jar/mat_id), zi = ~env + day, family = compois, data = repro_long)
list_asr[[57]] <-glmmTMB(offspring ~ regime + env + pop + day + regime:env + regime:pop + env:pop + regime:day + env:day + pop:day + regime:env:pop + regime:env:day + regime:pop:day + env:pop:day + regime:env:pop:day + (1 | jar/mat_id), zi = ~regime + day, family = compois, data = repro_long)
list_asr[[58]] <-glmmTMB(offspring ~ regime + env + pop + day + regime:env + regime:pop + env:pop + regime:day + env:day + pop:day + regime:env:pop + regime:env:day + regime:pop:day + env:pop:day + regime:env:pop:day + (1 | jar/mat_id), zi = ~pop + env + regime, family = compois, data = repro_long)
list_asr[[59]] <-glmmTMB(offspring ~ regime + env + pop + day + regime:env + regime:pop + env:pop + regime:day + env:day + pop:day + regime:env:pop + regime:env:day + regime:pop:day + env:pop:day + regime:env:pop:day + (1 | jar/mat_id), zi = ~pop + env + day, family = compois, data = repro_long)
list_asr[[60]] <-glmmTMB(offspring ~ regime + env + pop + day + regime:env + regime:pop + env:pop + regime:day + env:day + pop:day + regime:env:pop + regime:env:day + regime:pop:day + env:pop:day + regime:env:pop:day + (1 | jar/mat_id), zi = ~pop + regime + day, family = compois, data = repro_long)
list_asr[[61]] <-glmmTMB(offspring ~ regime + env + pop + day + regime:env + regime:pop + env:pop + regime:day + env:day + pop:day + regime:env:pop + regime:env:day + regime:pop:day + env:pop:day + regime:env:pop:day + (1 | jar/mat_id), zi = ~env + regime + day, family = compois, data = repro_long)
list_asr[[62]] <-glmmTMB(offspring ~ regime + env + pop + day + regime:env + regime:pop + env:pop + regime:day + env:day + pop:day + regime:env:pop + regime:env:day + regime:pop:day + env:pop:day + regime:env:pop:day + (1 | jar/mat_id), zi = ~pop + env + regime + day, family = compois, data = repro_long)

list_asr[[63]] <-glmmTMB(offspring ~ regime + env + pop + day + regime:env + regime:pop + env:pop + regime:day + env:day + pop:day + regime:env:pop + regime:env:day + regime:pop:day + env:pop:day + regime:env:pop:day + (1 | jar/mat_id), zi = ~1, family = poisson, data = repro_long)
list_asr[[64]] <-glmmTMB(offspring ~ regime + env + pop + day + regime:env + regime:pop + env:pop + regime:day + env:day + pop:day + regime:env:pop + regime:env:day + regime:pop:day + env:pop:day + regime:env:pop:day + (1 | jar/mat_id), zi = ~1, family = nbinom2, data = repro_long)
list_asr[[65]] <-glmmTMB(offspring ~ regime + env + pop + day + regime:env + regime:pop + env:pop + regime:day + env:day + pop:day + regime:env:pop + regime:env:day + regime:pop:day + env:pop:day + regime:env:pop:day + (1 | jar/mat_id), zi = ~1, family = genpois, data = repro_long)
list_asr[[66]] <-glmmTMB(offspring ~ regime + env + pop + day + regime:env + regime:pop + env:pop + regime:day + env:day + pop:day + regime:env:pop + regime:env:day + regime:pop:day + env:pop:day + regime:env:pop:day + (1 | jar/mat_id), zi = ~1, family = compois, data = repro_long)

list_asr[[67]] <-glmmTMB(offspring ~ regime + env + pop + day + regime:env + regime:pop + env:pop + regime:day + env:day + pop:day + regime:env:pop + regime:env:day + regime:pop:day + env:pop:day + regime:env:pop:day + (1 | jar/mat_id), zi = ~0, family = poisson, data = repro_long)
list_asr[[68]] <-glmmTMB(offspring ~ regime + env + pop + day + regime:env + regime:pop + env:pop + regime:day + env:day + pop:day + regime:env:pop + regime:env:day + regime:pop:day + env:pop:day + regime:env:pop:day + (1 | jar/mat_id), zi = ~0, family = nbinom2, data = repro_long)
list_asr[[69]] <-glmmTMB(offspring ~ regime + env + pop + day + regime:env + regime:pop + env:pop + regime:day + env:day + pop:day + regime:env:pop + regime:env:day + regime:pop:day + env:pop:day + regime:env:pop:day + (1 | jar/mat_id), zi = ~0, family = genpois, data = repro_long)
list_asr[[70]] <-glmmTMB(offspring ~ regime + env + pop + day + regime:env + regime:pop + env:pop + regime:day + env:day + pop:day + regime:env:pop + regime:env:day + regime:pop:day + env:pop:day + regime:env:pop:day + (1 | jar/mat_id), zi = ~0, family = compois, data = repro_long)

#model ranking via aic
model.sel(list_asr, rank = "AIC")
modelchecker(list_asr, filename = "outputs/zi.asr.four_")

#simulate residuals to test for dispersion issues
s1a <- simulateResiduals(list_asr[[49]], plot = T)
testZeroInflation(s1a)
testDispersion(s1a)
check_predictions(list_asr[[49]])
check_convergence(list_asr[[49]])
summary(list_asr[[49]])

#tab_model
tab_model(list_asr[[49]],
          transform = NULL,
          show.ci = FALSE,
          show.se = TRUE,
          show.r2 = FALSE,
          show.stat = TRUE,
          show.icc = FALSE)

#anova of full results
tidy(car::Anova(glmmTMB(offspring ~ regime + env + pop + day + regime:env + regime:pop +  
                          env:pop + regime:day + env:day + pop:day + regime:env:pop +  
                          regime:env:day + regime:pop:day + env:pop:day + regime:env:pop:day +  
                          (1 | jar/mat_id),
                        data = repro_long,
                        zi = ~env,
                        family = compois), type = "III")) %>%
  mutate_if(is.numeric, round, 5) %>%
  print()

#model means + bias adjustment
lme4::VarCorr(list_asr[[46]])
sigma <- sqrt(0.191^2 + 0.047^2)

em1a <- emtrends(list_asr[[46]],
                specs = ~ regime * env * pop,
                var = "day",
                bias.adjust = T,
                sigma = sigma)

confint(em1a, calc = c(n = ~.wgt.), adjust = "mvt",
        type = "response")

#model contrasts
summary(contrast(em1a, "consec",
                 combine = TRUE,
                 simple = list("pop", "regime", "env"),
                 adjust = "mvt"),
        type = "response", infer = c(TRUE, TRUE)) %>%
  print()

