# MINA
This repository contains data analysis and simulation code used in paper titled "Multiplexed imaging of nucleome architectures in single cells of mammalian tissue".

Instructions for use:

With new imaging data (single dataset), run the runM1-runM7 programs in the “RawDataProcessingCode” folder sequentially to extract the RNA MERFISH results. Run the runT1-runT7 programs in the “RawDataProcessingCode” folder sequentially to extract the chromatin tracing results. Then run the runMT8 program in the “RawDataProcessingCode” folder to extract the lamina and nucleolar profiles, and analyze the nucleome architectures by the cell types identified in this dataset. 

In the "SummaryAnalysisCode" folder, analyze the combined data from multiple datasets by running the run1-run8 programs sequentially. 

The "SimulationCode" folder contains code for polymer simulation of chromatin organization.


