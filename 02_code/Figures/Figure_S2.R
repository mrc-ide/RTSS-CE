# Figure S2
# Insecticide-treated net (ITN) use by country, 2015
# Data from the Malaria Atlas Project

# set-up
source("./02_code/Figures/data_and_libraries.R")


# https://malariaatlas.org/project-resources/modelling-coverage-of-insecticide-treated-nets-itns/
# import ITN data from MAP

nets_data <- read_csv("./01_data/Intervention_ITN.csv")

SSA_ITN <- read_csv("./01_data/00_ITN_table_Africa_admin1_2000-2020.csv") |>
  filter(Year == 2019) |>
  mutate(SSA = case_when(Name_0 %in% nets_data$Name ~ 1,
                         grepl("Cape|Comor|Ivoire|Lesotho|Maurit|Tome|Seyc|Tanz", Name_0) ~ 1)) |>
  filter(SSA == 1) |>
  group_by(Name_0) |>
  summarize(median = median(Mean_ITN_coverage_rate),
            min = min(Mean_ITN_coverage_rate),
            max = max(Mean_ITN_coverage_rate)) |>
  arrange(median) |>
  filter(!is.na(median)) |>
  mutate(Name_0 = factor(Name_0, levels=Name_0))

# plot
ggplot(SSA_ITN) +
  geom_rect(xmin=1, xmax=48, ymin=-0.01, ymax=.20, fill="#F6A1A5", alpha=0.009) +
  geom_rect(xmin=1, xmax=48, ymin=.20, ymax=.40, fill="#F8CD9C", alpha=0.01) +
  geom_rect(xmin=1, xmax=48, ymin=.40, ymax=.60, fill="#1BB6AF", alpha=0.005) +
  geom_rect(xmin=1, xmax=48, ymin=.60, ymax=.90, fill="#088BBE", alpha=0.006) +
  geom_pointrange(aes(x=Name_0, y=median, ymin=min, ymax=max), lwd=.3, size=0.3) +
  labs(x="Country", y="ITN use by country, 2019") +
  scale_x_discrete(labels = scales::wrap_format(20)) +
  scale_y_continuous(breaks = c(0, .20, .40, .60)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


ggsave("./03_output/figureS2.pdf", width=10, height=4)
