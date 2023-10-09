#load in packages
source("code/01-load_packages.R")

##data for all the traits##

##Load body mass data###
body_mass <- read_csv("data/processed/bodymass.csv") %>%
  convert(fct(Pop, Treat, Rep, Sex)) %>%
  mutate(Jar = paste(Pop, Treat, Rep, sep = "-")) %>%
  convert(fct(Jar)) %>%
  mutate(PopTrt = paste(Pop, Treat, sep = "-")) %>%
  convert(fct(PopTrt)) %>%
  as.data.frame() %>%
  mutate(Regime = case_when(
    Treat == "C" ~ "Constant",
    Treat == "I" ~ "Fluctuating",
    Treat == "CIA" ~ "Constant",
    Treat == "ICA" ~ "Fluctuating")) %>%
  mutate(Env = case_when(
    Treat == "C" ~ "Constant",
    Treat == "I" ~ "Fluctuating",
    Treat == "CIA" ~ "Fluctuating",
    Treat == "ICA" ~ "Constant")) %>%
  convert(fct(Regime, Env)) %>%
  mutate(mass_12 = (Beet.Mass / 12)*1000) %>%
  mutate(beet_mass_mg = Beet.Mass*1000) %>%
  janitor::clean_names()

###Load and create LRS data###
#Change the data to create LRS - no sex as sex not in development time file
development_time <- read_csv("data/processed/lifehistory.csv") %>%
  convert(fct(Pop, Treat, Rep)) %>%
  mutate(MatID = paste(Pop, Treat, Rep, VC, sep = "-")) %>%
  mutate(Jar = paste(Pop, Treat, Rep, sep = "-")) %>%
  convert(fct(Jar)) %>%
  mutate(matprob = ifelse(Total != 0, 1, 0)) %>%
  mutate(PopTrt = paste(Pop, Treat, sep = "-")) %>%
  convert(fct(PopTrt)) %>%
  mutate(TempTreat = case_when(
    Treat == "C" ~ "C",
    Treat == "I" ~ "S",
    Treat == "CIA" ~ "C-to-S",
    Treat == "ICA" ~ "S-to-C")) %>%
  convert(fct(TempTreat)) %>%
  mutate(Regime = case_when(
    Treat == "C" ~ "Constant",
    Treat == "I" ~ "Fluctuating",
    Treat == "CIA" ~ "Constant",
    Treat == "ICA" ~ "Fluctuating")) %>%
  mutate(Env = case_when(
    Treat == "C" ~ "Constant",
    Treat == "I" ~ "Fluctuating",
    Treat == "CIA" ~ "Fluctuating",
    Treat == "ICA" ~ "Constant")) %>%
  convert(fct(Regime, Env)) %>%
  janitor::clean_names()

#create total sum of data for LRS
repro_wide <- development_time %>%
  group_by(pop, env, regime, mat_id, rep, jar, day_egg) %>%
  summarise(totrep = sum(total)) %>%
  spread(day_egg, totrep) %>%
  replace(is.na(.), 0) %>%
  rename(day_0 = "0", day_1 = "1", day_2 = "2") %>%
  mutate(lrs = sum(c_across(day_0:day_2))) %>%
  as.data.frame() %>%
  mutate(obs = seq_len(nrow(.))) %>%
  convert(fct(obs)) %>%
  mutate(treatment_group = paste(pop, env, regime)) %>%
  janitor::clean_names()

###Age-specific reproduction###
repro_long <- development_time %>%
group_by(pop, env, regime, mat_id, rep, jar, day_egg) %>%
  summarise(offspring = sum(total)) %>%
  as.data.frame() %>%
  mutate(day = case_when(
    day_egg == "0" ~ "1",
    day_egg == "1" ~ "2",
   day_egg == "2" ~ "3")) %>%
  mutate(obs = seq_len(nrow(.))) %>%
  convert(fct(obs), num(day_egg, day)) %>%
  janitor::clean_names()


###Development time###
#convert the data and remove non zeros
development_time <- read_csv("data/processed/lifehistory.csv") %>%
  convert(fct(Pop, Treat, Rep)) %>%
  mutate(MatID = paste(Pop, Treat, Rep, VC, sep = "-")) %>%
  mutate(Jar = paste(Pop, Treat, Rep, sep = "-")) %>%
  convert(fct(Jar)) %>%
  mutate(matprob = ifelse(Total != 0, 1, 0)) %>%
  filter(Total != 0) %>%
  mutate(PopTrt = paste(Pop, Treat, sep = "-")) %>%
  convert(fct(PopTrt)) %>%
  mutate(TempTreat = case_when(
    Treat == "C" ~ "C",
    Treat == "I" ~ "S",
    Treat == "CIA" ~ "C-to-S",
    Treat == "ICA" ~ "S-to-C")) %>%
  convert(fct(TempTreat)) %>%
  mutate(Regime = case_when(
    Treat == "C" ~ "Constant",
    Treat == "I" ~ "Fluctuating",
    Treat == "CIA" ~ "Constant",
    Treat == "ICA" ~ "Fluctuating")) %>%
  mutate(Env = case_when(
    Treat == "C" ~ "Constant",
    Treat == "I" ~ "Fluctuating",
    Treat == "CIA" ~ "Fluctuating",
    Treat == "ICA" ~ "Constant")) %>%
  convert(fct(Regime, Env)) %>%
  janitor::clean_names()


#loop over each row to create a development_time for each individual.
dtlist <- list()

for (x in seq_len(nrow(development_time))) {

  store <- development_time[x, ] %>%
    slice(rep(seq_along(n()), each = development_time$total[x]))

  if (development_time[x, ]$males[1] == "0"  &&
     development_time[x, ]$females[1] == "0") {
    sex <- rep("NA", length = development_time[x, ]$total[1])
  } else {
    sex <- c(rep("M", length = development_time[x, ]$males[1]),
             rep("F", length = development_time[x, ]$females[1]))
  }

  dtlist[[x]] <- cbind(store, sex)

}

devtime <- data.table::rbindlist(dtlist)

#remove individuals where sex is not known as this is probably an error?
devtime %<>% convert(fct(sex)) %>%
  filter(sex != "NA")