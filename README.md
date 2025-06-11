# Fetal_MEG_frequency_tagging
Code for the analyses of fetal and neonatal MEG in Frohlich et al. 2025 BioRxiv
https://www.biorxiv.org/content/10.1101/2025.06.06.658307v1

Last updated 2025-06-11

Experiment 2 of this project uses archived data from Moser et al. 2020 (https://zenodo.org/records/4018827) and Moser et al. 2021 (https://zenodo.org/records/4541463). Please visit those archives on Zenodo for the corresponding data. Data from Experiment 1 will be uploaded at the time of publication. 

NOTE: I have tried to include as many external dependencies here as possible (and please note that the repo license only applies to my code and not these dependencies), but if you run into problems, try also adding Fieldtrip to your path before running my code (https://www.fieldtriptoolbox.org/). 

# Principal scripts

## Main analysis scripts

- **stat_learn_freq_tag_fetal.m:** runs the main analysis for fetal data from Experiment 1
- **local_global_freq_tag_fetal.m:** runs the main analysis for fetal data from Experiment 2
- **local_global_freq_tag_newborn.m:** runs the main analysis for neonatal data from Experiment 2

## Statistical analysis

- **freq_tag_FDR_correct.m:** runs the FDR correction across all the main tests
- **power_analysis.m:** Monte Carlo subsampling-based statistical power analysis, fetal data from Experiment 1

## Supplemental Material 

- **deHeering_valid_nonparametric_test.m:** validates the nonparametric approach with public data from de Heering and Rossion 2015 eLife
- **validate_boostrapped_correlation_test.m:** validates the boostrapped corrleation approach 
  
# Dependencies
## The above scripts have a lot of dependencies! I've tried my best to list and as many as possible below ... 

### Scripts (These include helper functions, my custom mods of native MATLAB functions, external code, etc.) 
### Note that many of these are publicly available external scripts! See code comments for proper attribution. 


- **customcolormap.m:** Creates alternative colormaps (downloaded from [MathWorks File Exchange](https://www.mathworks.com/matlabcentral/fileexchange/69470-custom-colormap)). **external code**
- **customcolormap_preset.m:** Works with the above function. **external code**
- **makefighandsome.m:** Makes aesthetic changes to figures.
- **makefigpretty.m:** Another one of my scripts to make aesthetic changes to figures.
- **mycolorbar.m:** Nicer-looking colorbar.
- **myfigure.m:** Bigger window for the figure.
- **myfigure2.m:** Even bigger window for the figure.
- **mylsline.m:** Nicer-looking least squares line.
- **sortref.m:** Sorts one list with respect to another.
- **ro_freq_wavelet_TFT2.m:** Computes the Morlet wavelet transform, written by Joerg Hipp with slight modifications. **external code**
- **surrogate.m:** Compute surrogate data. **external code**
- **table2latex.m:** External code for writing tables in LaTeX format. **external code**
- **DataHash.m:** Generates hash values for MATLAB data structures. **external code**
- **compare_correlation_coefficients.m:** Compares correlation coefficients statistically.  **external code**
- **powernoise.m:** Generates power law noise.
- **ft_topoplotER_JF.m:** Modified Fieldtrip function for topographic plotting.
- **topoplot_common_JF.m:** Modified Fieldtrip function for topographic plotting.

### Data Files (.mat, .csv, .xlsx)
#### Experiment 2 data from Zenodo archives 
- **Overview_Participants.xlsx:** 
- **dataset_differencesT3T4.csv:**
- **Data_traces.mat**
- **time_traces_data_manuscript.mat**

#### Statistical learning frequency tagging results (Experiment 1)
- **SL_fetal_stats_out.csv:** Experimental frequency output for fetal statistical learning frequency tagging analysis
- **SL_fetal_CONTROL_stats_out.csv:** Control frequency output for fetal statistical learning frequency tagging analysis
 
#### Local-global frequency tagging results (Experiment 2)
- **LG_fetal_stats_out.csv:** Experimental frequency output for fetal local-global frequency tagging analysis  
- **LG_fetal_CONTROL_stats_out.csv:** Control frequency output for fetal local-global frequency tagging analysis  
- **LG_newborn_stats_out.csv:** Experimental frequency output for neonatal local-global frequency tagging analysis  
- **LG_newborn_CONTROL_stats_out.csv:** Control frequency output for neonatal local-global frequency tagging analysis  

#### MATLAB workspace data

- **Alltable.mat**
- **stat_learning_output_RMS_27-May-2025.mat**
- **stat_learning_output_RMS_embedded_pattern_27-May-2025.mat**
- **stat_learn_hrv.mat:** 
- **stat_learning_demo_vars.mat:** 
- **stat_learning_hrv_26-Mar-2025.mat:** 

