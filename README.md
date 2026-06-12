# BePhyNE_ms_analyses
includes data and scripts for reporudcing major analyses and figures from the BePhyNE mansuscript

the repo is divided into two main directories: 

sim: Directory including all simulation analyses simulations with branch length transformations.


emp 
- data: includes data used in all analyses on Eastern North American Plethodontidae (Note on the Dryad location occurrence data has been removed due to copyright, for a full completely reproducible pipeline got to the github link: )
- scripts: includes all necessary files for performing:
	- data collection and inspection 
		BePhyNE_Get_PA_dryad.R
	    Inspecting_presence_data.R
	- BePhyNE model fitting with all data:
		Batch_BePhyNE_pletho_full.sh
		BePhyNE_pletho_full.R
		BePhyNE_pletho_full.sh
	- BePhyNE model fitting with missing data:
		Batch_BePhyNE_pletho_miss.sh
		Batch_process_miss_sp_logs.sh
		BePhyNE_pletho_miss.R
		BePhyNE_pletho_miss.sh
	-Plotting BePhyNE outputs:
		ridge_plots_pletho.R
		AUC_summarizing_plotting.R
		ecdf_mahalanobis_new_logdf_breadth_exp.R
- plots
	- folder where all plots are stored
- outfiles
	-folder where all outputs from runs are stored, currently stores many summarized outputs from larger log files from the large number of MCMC analyses performed
		
		
sim
- scripts: includes all necessary files for performing:
	- Kappa simulations:n 
		Kappa_BePhyNE_sim.R
		BePhyNE_Kappa.sh
		BatchBePhyNE_Kappa.sh
		plot_sim_comps_final.R
	- Background point simulations:
		BatchBePhyNE_Background.sh
		BePhyNE_Background.sh
		BePhyNE_collect.sh
		Background_BephyNE_sim_collect.R
		Background_BePhyNE_sim.R
		plot_sim_comps_final.R
	-Lambda Simulation
		BatchBePhyNE_Lambda.sh
		BePhyNE_Lambda.sh
		lambda_BePhyNE_sim.R
		plot_sim_comps_final.R
	- Height simulations and summarization
		BatchBePhyNE_Height_priors.sh
		BePhyNE_Height_priors.sh
		Height_BePhyNE_sim.R
		BePhyNE_Height_collect_priors.sh
		BePhyNE_Height_collect_priors.R
		BatchBePhyNE_Height_collect_priors.sh
		Batch_posterior_corr.sh
		posterior_correlation_summary.sh
		posterior_correlation_summary.R
		plot_correlation_histograms_with_representatives.R
		plot_sim_comps_final.R
		height_sim_accuracy_error_plots.R		

Most empirical and simulation analyses were performed on an HPC, with separate simulation experiment replicates/ Empirical MCMCs executed in parallel over separate CPUs using the shell scripts (.sh) scripts included.