library(dplyr)
library(ggplot2)
library(ggpattern)
library(geofacet)
library(ggh4x)
library(here)

# Main figure: Prospective impact -------------------------------------------------------------
txt_size = 10
theme_script <- theme_minimal() +
  theme(axis.text = element_text(colour="black", size = txt_size), 
        axis.ticks = element_line(colour="black"),
        legend.text = element_text(colour="black", size = txt_size),
        legend.title = element_text(colour="black", size = txt_size),
        axis.title = element_text(colour="black", size = txt_size),
        strip.text = element_text(colour = "black", size = txt_size))
theme_set(theme_script)

vacc_impact_main_plot <- read.csv(here("3_Output_Data", "Vaccine_impact", "vacc_impact_main_plot_observed_cases.csv")) %>%
  filter(!(name %in% c("YLLs averted", "YLDs averted", "Mild cases averted",
                       "Severe cases averted", "Severe cases averted per 10,000 doses",
                       "YLLs averted per 10,000 doses", 
                       "YLDs averted per 10,000 doses", 
                       "Mild cases averted per 10,000 doses"))) %>%
  mutate(name = gsub("All cases", "Cases", name))

cols <- c("#C39328FF", "#F4A3A8FF", "#E38191FF", "#AD466CFF","#8B3058FF", "#672044FF", "#00496FFF")


design <- matrix(c(
  1, NA,
  2, 3,
  4, 5,
  6, 7,
  8, 9,
  10, 11, 
  12, 13
), byrow = TRUE, ncol = 2)

col_pal <- c(cols[1], as.character(sapply(cols[2:length(cols)], rep, 2)))

p <- vacc_impact_main_plot %>%
  mutate(name = factor(name, levels = c("Doses required",
                                        "Infections averted", "Infections averted per 10,000 doses",
                                        "Detected cases averted", "Detected cases averted per 10,000 doses",
                                        "Cases averted", "Cases averted per 10,000 doses",
                                        "Chronic cases averted",
                                        "Chronic cases averted per 10,000 doses",
                                        "Deaths averted", "Deaths averted per 10,000 doses",
                                        "DALYs averted", "DALYs averted per 10,000 doses"))) %>%
  ggplot(aes(x = x, y = med)) +
  facet_manual(~ name, design = design, scales = "free_y") +
  geom_bar_pattern(aes(pattern = Subset, fill = name), stat = "identity",
                   pattern_fill = "white", pattern_color = "white", pattern_spacing = 0.1, pattern_density = 0.1) +
  geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.5) +
  scale_y_continuous(labels = scales::comma) +
  scale_pattern_manual("none", values = c("none", "stripe")) +
  scale_fill_manual(values = col_pal) +
  labs(y = NULL) +
  theme(axis.title.x=element_blank(), 
        axis.text.x=element_text(angle = 45, hjust = 1, vjust = 1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.position = "none")

ggsave(here("4_Output_Figures", "Fig_4_vaccine_impact_top_5_wo_severe_cases.pdf"), p, width=6.5, height = 8.4, units="in")
ggsave(here("4_Output_Figures", "Fig_4_vaccine_impact_top_5_wo_severe_cases.png"), p, width=6.5, height = 8.4, units="in")
ggsave(here("4_Output_Figures", "Fig_4_vaccine_impact_top_5_wo_severe_cases.svg"), p, width=6.5, height = 8.4, units="in")
ggsave(here("4_Output_Figures", "Fig_4_vaccine_impact_top_5_wo_severe_cases.eps"), p, width=6.5, height = 8.4, units="in")

# Sensitivity for potential impact -------------------------------------------------------------

theme_script <- theme_minimal() +
  theme(axis.text = element_text(colour="black", size = txt_size), 
        axis.ticks = element_line(colour="black"),
        legend.text = element_text(colour="black", size = txt_size),
        legend.title = element_text(colour="black", size = txt_size),
        axis.title = element_text(colour="black", size = txt_size + 1),
        strip.text = element_text(colour = "black", size = txt_size + 1),
        strip.background = element_rect(colour = "gray", fill = "white"))
theme_set(theme_script)

transform_names <- function(x) {
  replacements <- c(
    "Doses" = "Doses used",
    "C_all_averted" = "Cases averted",
    "D_averted" = "Deaths averted",
    "DALYs_averted" = "DALYs averted"
  )
  result <- replacements[x]
  return(result)
}

all_vacc_effect_df<- read.csv(here("3_Output_Data", "Vaccine_impact", "all_vacc_effect_forward_summary_df.csv")) %>%
  filter(Case !=  "Target > 50yo") %>%
  select(c(Type, Category, Case, Bound, Positive, Label,x_axis_lab,  Doses, ends_with("_med"))) %>%
  rename_with(.cols = ends_with("_med"), ~ gsub("_med", "", .x)) %>%
  mutate(across(where(is.numeric), ~ (.x - .x[which(Type == "Base case")])/.x[which(Type == "Base case")])) %>%
  tidyr::pivot_longer(cols = where(is.numeric), names_to = "quantity") %>%
  mutate(quantity = transform_names(quantity)) %>%
  filter(!is.na(quantity), Case != "Base case", x_axis_lab != "Prob. chronic\ndisease\nBC: 2.4%")%>%
  mutate(y_label = if_else(Positive, 0.8, -0.8),
         y_label = if_else(Positive & Label == "0.008%", 0.5, y_label),
         y_label = if_else(!Positive & Label == "18-59yo", -0.97, y_label)) %>%
  mutate(Category = if_else(Category == "Campaign", "Rollout campaign", Category),
         Category = factor(Category, levels = c("Vaccine", "Rollout campaign", "Clinical chikungunya", "Transmission dynamics")),
         quantity = factor(quantity, levels = c("Doses used", "Cases averted", "Deaths averted", "DALYs averted")),
         Bound = paste(Bound, "Bound"),
         x_axis_lab = factor(x_axis_lab, levels = c("Efficacy vs\ninfection\nBC: 70%", 
                                                    "Efficacy vs\ndisease\nBC: 95%",
                                                    "Number of doses\nBC: 1 dose", "Duration of \nprotection\nBC: Life-long",
                                                    "Coverage\nBC: 40%",
                                                    "Target age\nBC: >12yo",
                                                    "Routine vacc. in 12yo\nBC: Without",
                                                    "Target pop.\nBC: Brazil",
                                                    "Prob. disease\nBC: 50%",
                                                    "IFR\nBC: 0.004%",
                                                    "Annual FOI",
                                                    "Representative period\nBC: 2020-2024"))) 

cols <- c("#C39328FF", "#AD466CFF", "#672044FF", "#00496FFF")

design <- matrix(c(
  1, 1,
  2, 2,
  3, 4
), byrow = TRUE, ncol = 2)

p <- all_vacc_effect_df %>%
  ggplot(aes(x = x_axis_lab, y = value)) +
  facet_manual(~ Category, design = design, scales = "free_x") +
  geom_bar_pattern(data = subset(all_vacc_effect_df, Positive), 
                   aes(pattern = Bound, fill = quantity), stat = "identity", position = "dodge",
                   pattern_fill = "white", pattern_color = "white", pattern_spacing = 0.05, pattern_density = 0.2,
                   pattern_key_scale_factor = 0.3) +
  geom_bar_pattern(data = subset(all_vacc_effect_df, !Positive), 
                   aes(pattern = Bound, fill = quantity), stat = "identity", position = "dodge",
                   pattern_fill = "white", pattern_color = "white", pattern_spacing = 0.05, pattern_density = 0.2,
                   pattern_key_scale_factor = 0.3) +
  geom_hline(aes(yintercept = 0), color = "black") +
  scale_pattern_manual(breaks = c("Upper Bound", "Lower Bound"), values = c("none", "stripe")) +
  geom_label(aes(label = Label, y = y_label), label.r = unit(0, "pt"),
             vjust = "center", hjust = "center", label.padding = unit(3.5, "pt"),
             size = 10, size.unit = "pt", color = "black") +
  scale_fill_manual(values = cols) +
  scale_y_continuous(labels = scales::percent, limits = c(-1,1.001)) +
  scale_x_discrete(position = "top") +
  labs(y = "Percentage change", x = NULL, pattern = NULL, fill = NULL) +
  theme(legend.position = "bottom",
        axis.title.x.top = element_text(),
        axis.title.x.bottom = element_blank(),
        axis.ticks.x = element_blank(),
        strip.placement = "outside",
        strip.switch.pad.grid = unit(0.5, "cm"),
        panel.grid.major.x = element_blank(),
        strip.background = element_rect(color = "black")) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(pattern = "none")))

ggsave(here("4_Output_Figures", "Fig_5_sensitivity_vacc_impact_forward.pdf"), width = 6.5, height = 8.4, units = "in")
ggsave(here("4_Output_Figures", "Fig_5_sensitivity_vacc_impact_forward.png"), width = 6.5, height = 8.4, units = "in")
ggsave(here("4_Output_Figures", "Fig_5_sensitivity_vacc_impact_forward.svg"), width = 6.5, height = 8.4, units = "in")
ggsave(here("4_Output_Figures", "Fig_5_sensitivity_vacc_impact_forward.eps"), width = 6.5, height = 8.4, units = "in")
























