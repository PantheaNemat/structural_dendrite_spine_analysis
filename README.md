## structural_dendrite_spine_analysis

#### This repository contains resources to generate and analyse data on reconstructed dendrites and spines as presented in Nemat et al. (2025) - doi: 10.1016/j.nlm.2025.108033 - presented step by step in the three branches.



#### The "reconstruction_dendrites_spines" branch entails a step by step protocol to reconstruct dendrites and spines based on z-stacks in Imaris - Nemat_etal_2025_NLM_Appendix1_reconstruction_protocol.pdf.

#### The "data_compilation" branch contains a Python script that can be openeed in Juypter notebook, for instance, to compile individual data files obtained in Imaris into one big data file - Spine_Parameter_Data_Compiler.ipynb
#### FolderOrganisation_before_Running_Script.png shows how files should be ordered before running the Python script.
#### FolderOrganisation_after_Running_Script.png shows how the folder should look like after successful application of the script. 
#### Data_Compiler_Example.zip contains a folder with the raw data of this publication and aids as an example on how to run the Spine_Parameter_Data_Compiler.ipynb
#### The obtained data files can subsequently be used for statistical analyses. 
#### As an alternative to the Python script, the .zip file "Data_Compiler_Alternative_R.zip" contains a different and more flexible approach to compile data in R. The folder also contains images of the necessary folder organisation. 

#### In the "spine_parameter_analysis" branch, an R script is provided that can be used to perform spine type classification and to calculate spine clusters. 
#### The parameters can be adjusted to the users wishes. 
#### As we advise to perform multilevel modeling due to the nested data structure of dendrite and spine parameters, 
#### we also included an example on a generalised linear model taking the nested data structure into account and how to perform this in R.
#### The R script also contains usefule packages to perform these analyses. 
#### Further resources on multilevel modeling can be found in the supplementary material of this publication - Nemat_etal_2025_NLM_Appendix2_multilevelanalyses_resources.pdf. 
