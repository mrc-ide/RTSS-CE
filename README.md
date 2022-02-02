# GF-RTSS-CE  :earth_africa: :syringe: :heavy_dollar_sign:
Strategic resource allocation to maximize the impact of malaria response.

## Objectives
:one: To evaluate the relative cost-effectiveness and related uncertainty of introducing the RTS,S compared to scale-up of existing malaria interventions especially in settings where scale up of conventional interventions has not been achieved or where there is a risk of coverage slipping as the vaccine is being introduced, including addressing parameters such as variation by vaccine delivery mode across different settings (EPI, seasonal, hybrid).
  
   - As part of a threshold analysis: what is the minimum cost per fully vaccinated child to make the vaccine cost-effective compared to scaling existing malaria interventions for a given country setting?  

:two: What is the optimal sequencing of malaria interventions in order to maximize impact in resource constrained settings? 

   - What are the key determinants of the ranking on cost-effectiveness of malaria interventions? 

   - How do the uncertainties/assumptions on costs and effectiveness of malaria interventions play a role in this ranking? 

   - Is this ranking consistent across regions or countries with similar epidemiological profiles? 

   - What is the equity implication of alternative scale up scenarios? 


## Directory

```
.
├── 01_data                            # Data files
├── 02_code                            # R code file
|   ├── HPC                            # Copy of Q: drive folder for HPC runs
|   |   ├── functions.R                # Functions
|   |   ├── PfPR_EIR_match.R           # Matching model PfPR outputs to EIRs
|   |   ├── EIRestimates.rds           # Stored EIR estimates for select PfPR settings
|   |   ├── HPC_runs.R                 # Running malariasimulation scenarios on HPC
|   |   ├── HPC_processing.R           # Processing model output
├── 03_output                          # Figures and .rds files
├── GF_RTSS_CE.Rproj                   # R.Studio project file
└── README.md                          # Project overview

```
