# Figure S1
# Median P. falciparum prevalence (2-10 year olds) by country, 2019
# Data from the Malaria Atlas Project

# set-up
source("./02_code/Figures/data_and_libraries.R")


nets_data <- read_csv('./01_data/Intervention_ITN.csv')

PR <- read_csv('./01_data/00_PfPR_table_Global_admin0_2000-2019.csv') %>%
  filter(Year == 2019) %>% arrange(PfPR_median) %>%
  mutate(SSA = case_when(Name_0 %in% nets_data$Name ~ 1,
                         grepl('Cape|Comor|Ivoire|Lesotho|Maurit|Tome|Seyc|Tanz', Name_0) ~ 1)) %>%
  filter(SSA == 1) %>%
  mutate(Name_f = factor(Name_0, levels=Name_0))

ggplot(PR) +
  geom_rect(xmin=1, xmax=49, ymin=0.09, ymax=.11, fill="#F6A1A5", alpha=0.009) +
  geom_rect(xmin=1, xmax=49, ymin=.19, ymax=.21, fill="#F8CD9C", alpha=0.01) +
  geom_rect(xmin=1, xmax=49, ymin=.39, ymax=.41, fill="#1BB6AF", alpha=0.005) +
  geom_pointrange(aes(x=Name_f, y = PfPR_median, ymin = PfPR_LCI, ymax = PfPR_UCI), color='#088BBE') +
  labs(x='Country', y=expression(italic(Pf)~PR[2-10]~', 2019')) +
  scale_y_continuous(breaks = c(0,0.1,0.2,0.3,0.4)) +
  scale_x_discrete(labels = scales::wrap_format(20)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave('./03_output/figureS1.pdf', width=10, height=4)


# PfPR and mortality by country
nets_data <- read_csv('./01_data/Intervention_ITN.csv')

PR <- read_csv('./01_data/00_PfPR_table_Global_admin0_2000-2019.csv') %>%
  filter(Year == 2019) %>% dplyr::select(ISO, Name_0, starts_with('PfPR'))

mortality <- read_csv('./01_data/00_Pf_mortality_rate_table_Global_admin0_2000-2019.csv') %>%
  filter(Year == 2019) %>% dplyr::select(ISO, Name_0, starts_with('mortality'))

PRmort <- PR %>% left_join(mortality) %>%
  arrange(PfPR_median) %>%
  mutate(SSA = case_when(Name_0 %in% nets_data$Name ~ 1,
                         grepl('Cape|Comor|Ivoire|Lesotho|Maurit|Tome|Seyc|Tanz', Name_0) ~ 1)) %>%
  filter(SSA == 1) %>%
  mutate(Name_f = factor(Name_0, levels=Name_0))

ggplot(PRmort) +
  geom_rect(xmin=1, xmax=49, ymin=0.09, ymax=.11, fill="#F6A1A5", alpha=0.009) +
  geom_rect(xmin=1, xmax=49, ymin=.19, ymax=.21, fill="#F8CD9C", alpha=0.01) +
  geom_rect(xmin=1, xmax=49, ymin=.39, ymax=.41, fill="#1BB6AF", alpha=0.005) +
  geom_pointrange(aes(x=Name_f, y = PfPR_median, ymin = PfPR_LCI, ymax = PfPR_UCI, color='PfPR'), alpha=0.4) +
  geom_pointrange(aes(x=Name_f, y = mortality_rate_median*200, ymin = mortality_rate_LCI*200, ymax = mortality_rate_UCI*200, color='mortality'), alpha=0.4) +
  labs(x='Country', y=expression(italic(Pf)~PR[2-10]~', 2019'), color='') +
  scale_y_continuous(sec.axis = sec_axis(~ . / 200,
                                         name = expression(italic(Pf)~'mortality rate (all ages)')),
                     breaks = c(0,0.1,0.2,0.3,0.4)) +
  scale_x_discrete(labels = scales::wrap_format(20)) +
  scale_color_manual(values = c('#088BBE','blue')) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

