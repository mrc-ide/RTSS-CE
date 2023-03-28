# Tabular results

# set-up
source("./02_code/Figures/data_and_libraries.R")


# < the impact of RTSS on top of other interventions --------------
output <- scenarios |>
  filter(cost_per_dose == 12.01 & delivery_cost == 1.62) |>
  mutate(ID = paste(pfpr, seasonality, ITNuse, resistance, treatment, ITN, sep = "_")) |>
  filter(ITNuse == 0.60)

none <- output |>
  mutate(set = case_when(seasonality %in% c("perennial", "highly seasonal") & intervention %in%
                           c("ITN 10% increase","ITN PBO") ~ 1,
                         seasonality == "seasonal" & intervention %in%
                           c("ITN 10% increase + SMC","ITN PBO + SMC") ~ 1)) |>
  filter(set == 1) |>
  dplyr::select(file, ID, drawID, daly, cases, cost_total, u5_dalys, n_0_1825) |>
  rename(daly_baseline = daly,
         cases_baseline = cases,
         cost_total_baseline = cost_total,
         u5_daly_baseline = u5_dalys) |>
  dplyr::select(file, ID, drawID, daly_baseline, cases_baseline, cost_total_baseline, u5_daly_baseline)

base_IDs <- none$file

output2 <- output |> filter(!(file %in% base_IDs)) |>
  mutate(set = case_when(seasonality %in% c("perennial", "highly seasonal") & intervention %in%
                           c("ITN 10% increase + RTS,S","ITN PBO + RTS,S") ~ 1,
                         seasonality=="seasonal" & intervention %in%
                           c("ITN 10% increase + RTS,S + SMC","ITN PBO + RTS,S + SMC") ~ 1)) |>
  filter(set == 1) |>
  dplyr::select(file, ID, drawID, pfpr, seasonality, intervention, daly, cases, cost_total, u5_dalys, dose3) |>
  left_join(none |> dplyr::select(-file), by=c("ID", "drawID")) |>
  mutate(CE = (cost_total - cost_total_baseline) / (daly_baseline - daly),
         deltadaly = daly_baseline - daly,
         deltacases = cases_baseline - cases,
         CE_u5 = (cost_total - cost_total_baseline) / (u5_daly_baseline - u5_dalys))


# adjusting for 15 year simulation period and 200,000 population arguments
summary(output2$deltadaly / (2*15)) # additional dalys averted per year in a population of 100,000 people
summary(output2$deltacases / (2*15)) # additional cases averted per year in a population of 100,000 people
summary(output2$deltacases / output2$dose3 * 100000) # additional cases averted per year per fully vaccinated child (dose3)
summary(output2$CE_u5) # additional cases averted per year in a population of 100,000 people


# < rankings table for guidance note --------------
output <- scenarios |>
  filter(cost_per_dose == 12.01 & delivery_cost == 1.62) |>
  filter(intervention %in% c("ITN 10% increase","ITN PBO","SMC","RTS,S age-based","RTS,S seasonal")) |>
  mutate(seasonality = factor(seasonality, levels = c("perennial", "seasonal", "highly seasonal"))) |>
  group_by(ID, drawID) |> arrange(ID, drawID, CE) |>
  slice(1L)

table <- output |>
  group_by(seasonality, pfpr, ITNuse, SMC, resistance, intervention) |>
  summarize(n = n()) |>
  rename(SMCuse = SMC) |>
  mutate(SMCuse = ifelse(seasonality == "seasonal", 0, SMCuse)) |>
  group_by(seasonality, pfpr, ITNuse, SMCuse, resistance) |>
  mutate(rank = rank(desc(n))) |>
  select(-n) |>
  pivot_wider(names_from = "rank", values_from = "intervention") |>
  mutate(across(`1`:`4`, ~replace_na(.x, "-")))

head(table)

table_n <- output |>
  group_by(seasonality, pfpr, ITNuse, SMC, resistance, intervention) |>
  summarize(n = n()) |>
  rename(SMCuse = SMC) |>
  mutate(SMCuse = ifelse(seasonality == "seasonal", 0, SMCuse)) |>
  group_by(seasonality, pfpr, ITNuse, SMCuse, resistance) |>
  mutate(rank = rank(desc(n))) |>
  select(-intervention) |>
  pivot_wider(names_from = "rank", values_from = "n")


head(table_n)

# copy out to clipboard
table |> print(noSpaces = T) |> write.table("clipboard", sep = "\t")
table_n |> print(noSpaces = T) |> write.table("clipboard", sep = "\t")





