# Q: what prevalence values should we use in our runs?
# A: look at the global distribution of PfPRs and what ITN coverages these countries use

library(tidyverse)

# pull in PfPR data from MAP layer
# https://malariaatlas.org/
PR <- read_csv('./01_data/00_PfPR_table_Global_admin0_2000-2019.csv') %>%
  filter(Year==2019) %>%
  select(ISO, PfPR_rmean)

# pull in ITN coverage by country
# https://github.com/bertozzivill/map-itn-cube/blob/publication-2021/paper_figures/figure_data/fig_4_access_npc.csv
ITN <- read.csv('https://raw.githubusercontent.com/bertozzivill/map-itn-cube/publication-2021/paper_figures/figure_data/fig_4_access_npc.csv') %>%
  group_by(iso3) %>%
  summarize(access_mean = mean(access_mean, na.rm=T))

# combine
output <- PR %>% left_join(ITN, by=c("ISO"="iso3")) %>% filter(!is.na(access_mean))

# plot
ggplot(output, aes(x=PfPR_rmean, y=access_mean)) +
  geom_point(color = "cornflowerblue") +
  labs(title = "PfPR and ITN access from MAP",
       x = "mean PfPR, 2019",
       y = "mean ITN access, 2020") +
  geom_text(label=output$ISO, nudge_x = 0, nudge_y = .03, check_overlap = T, size = 2) +
  theme_classic()

# save
ggsave("./03_output/PR_ITN_dist.pdf", width=5, height=5)

