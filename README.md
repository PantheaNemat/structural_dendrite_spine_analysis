# structural_dendrite_spine_analysis

# This repository contains resources to generate and analyse data on reconstructed dendrites and spines. 

# It entails a step by step protocol to reconstruct dendrites and spines based on z-stacks in Imaris - Nemat_etal_Appendix1_protocol_reconstruction_dendrites_spines.docx

# Exported data on reconstructed dendrites and spines can then be compiled based on the following script: Spine_Parameter_Data_Compiler.ipynb
# Step1.JPG and Step2.JPG contain screenshots of the necessary folder organisation. 
# Data_Compiler_Example.zip contains a folder with the raw data of this publication and aids as an example on how to run the Spine_Parameter_Data_Compiler.ipynb
# The obtained data files can subsequently be used for statistical analyses. 

# The provided R script can be used to perform spine type classification and to calculate spine clusters. 
# As we advise to perform multilevel modeling due to the nested data structure of dendrite and spine parameters, 
# we also included an exmample on a generalised linear model taking the nested data structure into account and how to perform this in R.
# The R script also contains usefule packages to perform these analyses. 
# Further resources on multilevel modeling can be found in the supplementary material of this publication. 
