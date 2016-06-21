# Tw2
Multivariate Welch t-test on distances

## File listing
1. `code/` -- R source code for simulations, figure and table generation
  * `Tw2.R` -- implementation to Tw2, distance Cohen's d
  * `simulateDists.R` -- simulation of distance matrices with prescribed effect size and heteroscedasticity
  * `simulateAll.R` -- simulation of the main data set for the manuscript
  * `simulate10.R` -- simulation of data with sample sizes in the vicinity of 10 per group
  * `simulate50.R` -- simulation of data with sample sizes in the vicinity of 50 per group
  * `heatmaps.R` -- generation of figures and tables for main simulation study
  * `plot10.R` -- generation of figures and tables for simulation study with sample sizes in the vicinity of 10 per group
  * `plot50.R` -- generation of figures and tables for simulation study with sample sizes in the vicinity of 50 per group
  * `STAT.R` -- analysis of sub-therapeutic antibiotic treatment data
  * `Psoriasis.R`	-- analysis of psoriasis skin microbiome data
2. `data.zip` -- Data for applications
  * `obesity_mapping.txt` -- metadata for STAT dataset
  * `seqs_otu_table.txt` -- abundace table for STAT dataset
  * `PsoriasisMetadata.txt`	-- metadata for psoriasis dataset
  * `otu_table.txt` -- abundace table for psoriasis dataset
3. `results/`
  * `simulation_results2.txt.zip` -- main simulation results
  * `simulation10_results.txt.zip` -- simulation results for sample sizes in the vicinity of 10 per group
  * `simulation50_results.txt.zip` -- simulation results for sample sizes in the vicinity of 50 per group

