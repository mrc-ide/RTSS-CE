# RTSS-CE  :earth_africa: :syringe: :heavy_dollar_sign:
Strategic resource allocation to maximize the impact of malaria response.

## Objectives
:one: To evaluate the relative cost-effectiveness and related uncertainty of introducing the RTS,S malaria vaccine compared to scale-up of existing malaria interventions. 

   - As part of a threshold analysis: what is the minimum cost per fully vaccinated child to make the vaccine cost-effective compared to scaling existing malaria interventions for a given country setting?  

:two: What is the optimal sequencing of malaria interventions in order to maximize impact in resource constrained settings? 

   - What are the key determinants of the cost-effectiveness rankings? 

   - How do the uncertainties/assumptions on costs and effectiveness of malaria interventions play a role in this ranking? 

   - Is this ranking consistent across regions or countries with similar epidemiological profiles? 

   - What is the equity implication of alternative scale up scenarios? 


## Directory

```
.
├── 01_data                                # Data files
├── 02_code                                # R code file
|   ├── HPC_draws                          # HPC runs for parameter draws
|   |   ├── A_parameter_draws.R              # List parameters for each draw
|   |   ├── B_PfPR_EIR_match.R               # Match PfPR to EIR
|   |   ├── C_HPC_runs.R                     # HPC runs
|   |   ├── function_draws.R                 # Functions for HPC
|   |   ├── Processing                     # Data processing
|   |   |   ├── cost_effectiveness.R         # Process CE
|   |   |   ├── HPC_processing.R             # Helper-function HPC output
|   |   |   ├── add_costs.R                  # Helper-function add costs
|   |   |   ├── deaths_dalys.R               # Helper-function calc deaths & DALYs
|   |   |   ├── netz_dist.R                  # Helper-function calc nets distributed
|   |   |   ├── outcome_averted.R            # Helper-function calc DALYs and cases averted
|   |   ├── Processing_casestudy           # Data processing - case study
|   |   |   ├── cost_effectiveness.R         # Process CE
|   |   |   ├── HPC_processing.R             # Helper-function HPC output
|   ├── Figures                            # Re-create paper figures
|   ├── MISC_calibrate_test.R              # Test calibration output function
|   ├── MISC_ITN_usage_distribution.R      # Determine model input ITN coverage values to reach population usages
|   ├── MISC_SMCdosetest.R                 # Determine which SMC profiles should be used
|   ├── MISC_warmup_time_test.R            # Determine the warmup time needed
├── 03_output                              # PDF figures and .rds files
├── GF_RTSS_CE.Rproj                       # R.Studio project file
└── README.md                              # Project overview

```

## R package versions

malariasimulation: 1.3.0

malariaequilibrium: 1.0.1

individual: 0.1.9

cali: https://github.com/mrc-ide/cali/tree/adbc370450716a2231347094c68cf7b31419896c

netz: https://github.com/mrc-ide/netz/tree/0ce236c34f36cb4d533260f88563dc4e09d3ba65

treasure: https://github.com/mrc-ide/treasure/tree/36a5498f4a3bea07e4a6cb4390e955dac359409a
