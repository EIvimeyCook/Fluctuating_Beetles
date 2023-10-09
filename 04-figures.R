#source the other code####
source("code/01-load_packages.R")
source("code/02-load_clean_data.R")

###Body mass#####
distribeaut_shape(data = body_mass,
 xvar = "env",
  yvar = "beet_mass",
   fillvar = "regime",
    shapevar = "sex") +
  labs(x = "Assay Environment",
   y = "Individual Body Mass (g)",
    fill = "Thermal Regime",
     colour = "Thermal Regime",
      shape = "Sex") +
  facet_grid(.~pop)
  scale_shape_manual(values = c(25, 21))
  
#emmeans plot
  em1a <- emmeans(m1a,
                  specs = ~ sex * pop * regime * env)
  
  plot_means <- as_tibble(confint(em1a, type = "response", calc = c(n = ~.wgt.),
                        adjust = "mvt")) %>%
    mutate(env_short = case_when(env == "Fluctuating" ~ "F",
                                 env == "Constant" ~ "C")) %>%
    mutate(regime_short = case_when(regime == "Fluctuating" ~ "F",
                                 regime == "Constant" ~ "C")) %>%
    mutate(sex_long = case_when(sex == "M" ~ "Male",
                                    sex == "F" ~ "Female")) %>%
    mutate(Treatment_relevel = paste0(pop,"-",regime_short, "-", env_short, sex)) %>%
    mutate(Treatment_relevel = fct_relevel(Treatment_relevel,
                                   "LEIC-C-CF",
                                   "LEIC-C-CM",
                                   "USA-C-CF",
                                   "USA-C-CM",
                                   "LEIC-F-CF",
                                   "LEIC-F-CM",
                                   "USA-F-CF",
                                   "USA-F-CM",
                                   "LEIC-C-FF",
                                   "LEIC-C-FM",
                                   "USA-C-FF",
                                   "USA-C-FM",
                                   "LEIC-F-FF",
                                   "LEIC-F-FM",
                                   "USA-F-FF",
                                   "USA-F-FM")) %>%
    mutate(treatment_plot = paste0(regime_short, "-",env_short, "\n", "(",pop,")"))
  
  
  g1 <- ggplot(plot_means,
         aes(x = treatment_plot,
             y = emmean, 
             colour = interaction(regime, env),
             shape = sex)) +
    geom_point(size = 4,  position = position_dodge(width = 0.25)) + 
    scale_colour_manual(values = c("blue", "darkblue", "red", "darkred")) + 
    geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.3, 
                  linewidth = 1.2,
                  position = position_dodge(width = 0.25)) + 
    theme(axis.text.x = ggtext::element_markdown()) +
    #facet_grid(.~sex_long, scales = "free") + 
    theme_publication() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none") +
    labs(y = "Average Mass (mg)", x = NULL)
  
  contrast_means <- as_tibble(summary(contrast(em1a, "pairwise",
                   combine = TRUE,
                   simple = list("pop","sex", "regime", "env", c("regime", "env")),
                   adjust = "mvt"),
          type = "response", infer = c(TRUE, TRUE))) %>%
    slice(9:16) %>%
    mutate(env_short = case_when(env == "Fluctuating" ~ "F",
                                 env == "Constant" ~ "C")) %>%
    mutate(regime_short = case_when(regime == "Fluctuating" ~ "F",
                                    regime == "Constant" ~ "C")) %>%
    mutate(treatment_plot = paste0(regime_short, "-",env_short, "\n", "(",pop,")"))
  
  
  g2 <- ggplot(contrast_means,
               aes(x = treatment_plot,
                   y = estimate, 
                   colour = interaction(regime,env),
                   shape = sex)) +
    geom_point(size = 4,  position = position_dodge(width = 0.25)) + 
    scale_colour_manual(values = c("blue", "darkblue", "red", "darkred")) + 
    geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.3, 
                  linewidth = 1.2,
                  position = position_dodge(width = 0.25)) + 
    theme(axis.text.x = ggtext::element_markdown()) +
    #facet_grid(.~sex_long, scales = "free") + 
    theme_publication() +
    geom_hline(yintercept = 0) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none") +
    labs(y = "Sexual Dimorphism (Female-Male)", x = "Treatment: Regime-Assay (Background)")
  
g1/g2 + plot_annotation(tag_levels = 'A')
  
ggsave(path = "figures/",
 filename = "bodymass.png",
  width = 250,
   height = 270,
   units = "mm",
    dpi = 300)

###LRS#####
lrs <-
distribeaut(data = repro_wide,
 xvar = "env",
  yvar = "lrs",
   fillvar = "regime") +
  labs(x = "Assay Environment",
   y = "Lifetime Reproductive Success",
    fill = "Thermal Regime",
    colour = "Thermal Regime")

#emmeans plot
  em1a <- emmeans(list_lrs_disp[[7]],
                  specs = ~ pop * regime * env,
                  type = "response")
  
  plot_means <- as_tibble(confint(em1a, type = "response", calc = c(n = ~.wgt.),
                        adjust = "mvt")) %>%
    mutate(env_short = case_when(env == "Fluctuating" ~ "F",
                                 env == "Constant" ~ "C")) %>%
    mutate(regime_short = case_when(regime == "Fluctuating" ~ "F",
                                 regime == "Constant" ~ "C")) %>%
    mutate(Treatment_relevel = paste0(pop,"-",regime_short, "-", env_short)) %>%
    mutate(treatment_plot = paste0(regime_short, "-",env_short, "\n", "(",pop,")"))
  
  
  lrs <- ggplot(plot_means,
         aes(x = treatment_plot,
             y = response, 
             colour = interaction(regime,env))) +
    geom_point(size = 4,  position = position_dodge(width = 0.25)) + 
    scale_colour_manual(values = c("blue", "darkblue", "red", "darkred")) + 
    geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.3, 
                  linewidth = 1.2,
                  position = position_dodge(width = 0.25)) + 
    theme(axis.text.x = ggtext::element_markdown()) +
    #facet_grid(.~sex_long, scales = "free") + 
    theme_publication() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none") +
    labs(y = "Total offspring", x = "Treatment: Regime-Assay (Background)")

ggsave(path = "figures/",
 filename = "lrs.png",
  width = 297, height = 250,
   units = "mm", dpi = 300)

####Age-specific reproduction####

dodge <- position_dodge(0.2)

asr <- ggplot(repro_long,
 aes(x = day,
  y = offspring,
   colour = regime,
   linetype = env,
    group = interaction(regime, env, pop))) +
  stat_summary(fun.data = "mean_cl_boot",
   geom = "errorbar",
    width = 0.2,
      position = dodge) +
  stat_summary(fun.data = "mean_cl_boot",
   geom = "point",
    size = 2,
      position = dodge) +
  stat_summary(fun.data = "mean_cl_boot",
   geom = "line",
    size = 1.1,
      position = dodge) +
  theme_publication() +
  ylab("Offspring Count") +
  xlab("Day") +
  scale_linetype_manual(values = c("dotdash", "solid")) +
  guides(linetype = guide_legend(title = "Assay Environment"),
         colour = "none") +
  facet_grid(. ~ pop) +
  coord_cartesian(y = c(0, 40))

#emmeans

em1a <- emmeans(list_asr[[46]],
                specs = ~ env * pop * regime*day,
                at = list(day = c(1,2,3)),
                type = "response")

plot_means <- as_tibble(confint(em1a, type = "response", calc = c(n = ~.wgt.),
                                adjust = "mvt")) %>%
  mutate(env_short = case_when(env == "Fluctuating" ~ "F",
                               env == "Constant" ~ "C")) %>%
  mutate(regime_short = case_when(regime == "Fluctuating" ~ "F",
                                  regime == "Constant" ~ "C")) %>%
  mutate(Treatment_relevel = paste0(pop,"-",regime_short, "-", env_short)) %>%
  mutate(treatment_plot = paste0(regime_short, "-",env_short, "\n", "(",pop,")"))


asr <- ggplot(plot_means,
              aes(x = day,
                  y = response, 
                  colour = interaction(regime,env))) +
  geom_point(aes(group = treatment_plot), size = 4,  position = position_dodge(width = 0.25)) + 
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL,
                    group = treatment_plot), width = 0.3, 
                linewidth = 1.2,
                position = position_dodge(width = 0.25)) + 
  geom_line(aes(group = treatment_plot),
            position = position_dodge(width = 0.25)) +
  theme(axis.text.x = ggtext::element_markdown()) +
  facet_grid(.~pop, scales = "free") + 
  theme_publication() +
  labs(y = "Age-specific reproduction", x = "Day") +
  scale_x_continuous(breaks = seq(from = 1, to = 3, by = 1)) +
  scale_colour_manual(labels=c('C-C', 'F-C', "C-F", "F-F"), name = "Treatment: Regime-Assay",
                        values = c("blue", "darkblue", "red","darkred")) +
  guides(shape = "none")
  
lrs / asr +
plot_annotation(tag_levels = "A")

ggsave(path = "figures/",
 filename = "fitness.png",
  width = 297,
   height = 260,
    units = "mm",
     dpi = 300)


###Development Time####
distribeaut_shape(data = devtime,
 xvar = "env",
  yvar = "dt",
   fillvar = "pop",
    shapevar = "sex") +
  labs(x = "Assay Environment",
   y = "Development Time (Days)",
    fill = "Genetic Background",
     colour = "Genetic Background",
      shape = "Sex")

#emmeans plot
em1a <- emmeans(m1a,
                specs = ~ sex * pop * regime * env)

plot_means <- as_tibble(confint(em1a, type = "response", calc = c(n = ~.wgt.),
                                adjust = "mvt")) %>%
  mutate(env_short = case_when(env == "Fluctuating" ~ "F",
                               env == "Constant" ~ "C")) %>%
  mutate(regime_short = case_when(regime == "Fluctuating" ~ "F",
                                  regime == "Constant" ~ "C")) %>%
  mutate(sex_long = case_when(sex == "M" ~ "Male",
                              sex == "F" ~ "Female")) %>%
  mutate(Treatment_relevel = paste0(pop,"-",regime_short, "-", env_short, sex)) %>%
  mutate(Treatment_relevel = fct_relevel(Treatment_relevel,
                                         "LEIC-C-CF",
                                         "LEIC-C-CM",
                                         "USA-C-CF",
                                         "USA-C-CM",
                                         "LEIC-F-CF",
                                         "LEIC-F-CM",
                                         "USA-F-CF",
                                         "USA-F-CM",
                                         "LEIC-C-FF",
                                         "LEIC-C-FM",
                                         "USA-C-FF",
                                         "USA-C-FM",
                                         "LEIC-F-FF",
                                         "LEIC-F-FM",
                                         "USA-F-FF",
                                         "USA-F-FM")) %>%
  mutate(treatment_plot = paste0(regime_short, "-",env_short, "\n", "(",pop,")"))

g1 <- ggplot(plot_means,
             aes(x = treatment_plot,
                 y = response, 
                 colour = interaction(regime,env),
                 shape = sex)) +
  geom_point(size = 4,  position = position_dodge(width = 0.25),
             show.legend = F) + 
  scale_colour_manual(values = c("blue", "darkblue", "red", "darkred")) + 
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.3, 
                linewidth = 1.2,
                position = position_dodge(width = 0.25),
                show.legend = F) + 
  theme(axis.text.x = ggtext::element_markdown(),
        legend.position = "none") +
  #facet_grid(.~sex_long, scales = "free") + 
  theme_publication() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(y = "Development Time (Days)", x = NULL)

contrast_means <- as_tibble(summary(contrast(regrid(em1a), "pairwise",
                                             combine = TRUE,
                                             simple = list("pop","sex", "regime", "env", c("regime", "env")),
                                             adjust = "mvt"),
                                    type = "response", infer = c(TRUE, TRUE))) %>%
  slice(9:16) %>%
  mutate(env_short = case_when(env == "Fluctuating" ~ "F",
                               env == "Constant" ~ "C")) %>%
  mutate(regime_short = case_when(regime == "Fluctuating" ~ "F",
                                  regime == "Constant" ~ "C")) %>%
  mutate(treatment_plot = paste0(regime_short, "-",env_short, "\n", "(",pop,")"))


g2 <- ggplot(contrast_means,
             aes(x = treatment_plot,
                 y = estimate, 
                 colour = interaction(regime,env),
                 shape = sex)) +
  geom_point(size = 4,  position = position_dodge(width = 0.25)) + 
  scale_colour_manual(values = c("blue", "darkblue", "red", "darkred")) + 
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.3, 
                linewidth = 1.2,
                position = position_dodge(width = 0.25)) + 
  theme(axis.text.x = ggtext::element_markdown()) +
  #facet_grid(.~sex_long, scales = "free") + 
  theme_publication() +
  geom_hline(yintercept = 0) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  labs(y = "Sexual Dimorphism (Female-Male)", x = "Treatment: Regime-Assay (Background)")

g1/g2 + plot_annotation(tag_levels = 'A')

ggsave(path = "figures/",
 filename = "development_time.png",
 width = 250,
 height = 270,
   units = "mm",
    dpi = 300)
