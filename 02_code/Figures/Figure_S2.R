# Figure S2
# Insecticide-treated net (ITN) use by country, 2015
# Data from the Malaria Atlas Project

# set-up
source("./02_code/Figures/data_and_libraries.R")


# https://malariaatlas.org/research-project/the-impact-of-malaria-control-on-plasmodium-falciparum-in-africa-2000-2015/
# import ITN data from MAP
SSA_ITN <- readRDS('./03_output/MAP_ITN.rds') %>% arrange(median) %>%
  filter(!is.na(median)) %>%
  mutate(name_0 = factor(name_0, levels=name_0))

ggplot(SSA_ITN) +
  geom_rect(xmin=1, xmax=47, ymin=-0.01, ymax=.25, fill="#F6A1A5", alpha=0.009) +
  geom_rect(xmin=1, xmax=47, ymin=.25, ymax=.50, fill="#F8CD9C", alpha=0.01) +
  geom_rect(xmin=1, xmax=47, ymin=.50, ymax=.75, fill="#1BB6AF", alpha=0.005) +
  geom_rect(xmin=1, xmax=47, ymin=.75, ymax=.90, fill="#088BBE", alpha=0.006) +
  geom_pointrange(aes(x=name_0, y=median, ymin=min, ymax=max), lwd=.3) +
  labs(x='Country', y='ITN use by country, 2019') +
  scale_x_discrete(labels = scales::wrap_format(20)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


ggsave('./03_output/figureS2.pdf', width=10, height=4)
