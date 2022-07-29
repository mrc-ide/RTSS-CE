# GF-RTSS-CE  :earth_africa: :syringe: :heavy_dollar_sign:
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
|   |   |   ├── add_costs.R                  # Helper-function add costs
|   |   |   ├── deaths_dalys.R               # Helper-function calc deaths & DALYs
|   |   |   ├── netz_dist.R                  # Helper-function calc nets distributed
|   |   |   ├── outcome_averted.R            # Helper-function calc DALYs and cases averted
|   ├── HPC_median                         # HPC runs for median parameter values
|   |   ├── A_PfPR_EIR_match.R               # Match PfPR to EIR
|   |   ├── B_HPC_runs.R                     # Run malariasimulation scenarios on HPC
|   |   ├── C_HPC_processing.R               # Process model output
|   |   ├── D_HPC_case_studies.R             # Process equity analysis output
|   |   ├── E_cost_effectiveness.R           # Process model output, cost-effectiveness variables
|   |   ├── functions.R                      # Functions for HPC
|   ├── Figures                            # Re-create paper figures
|   ├── MISC_ITN_PR_distributions.R        # Visualize world ITN use and PfPR patterns
|   ├── MISC_seasonality_script.R          # Visualize timing of seasonal interventions
|   ├── MISC_SMCdosetest.R                 # Determine which SMC profiles should be used
|   ├── MISC_warmup_time_test.R            # Determine the warmup time needed
├── 03_output                              # PDF figures and .rds files
├── GF_RTSS_CE.Rproj                       # R.Studio project file
└── README.md                              # Project overview

```
