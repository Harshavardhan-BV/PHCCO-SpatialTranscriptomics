# PHCCO-SpatialTranscriptomics
Material for Summer Training Programme on Spatial Transcriptomics

## Pre-requisites
- It is necessary to create a virtual environment with python 3.11. [uv](https://docs.astral.sh/uv/) is recommended for creating the virtual environment and handling the dependencies. One can also use poetry+pyenv or conda for this. Install the packages using:
    ```bash
    uv sync
    ```
- Download the required dataset files from 10x Genomics:
    - [`filtered_feature_bc_matrix.h5`](https://cf.10xgenomics.com/samples/spatial-exp/2.1.0/CytAssist_FFPE_Protein_Expression_Human_Breast_Cancer/CytAssist_FFPE_Protein_Expression_Human_Breast_Cancer_filtered_feature_bc_matrix.h5)
    - [`spatial.tar.gz`](https://cf.10xgenomics.com/samples/spatial-exp/2.1.0/CytAssist_FFPE_Protein_Expression_Human_Breast_Cancer/CytAssist_FFPE_Protein_Expression_Human_Breast_Cancer_spatial.tar.gz)
- Place both files inside the `Data/CytAssist_FFPE_Protein_Expression_Human_Breast_Cancer` directory.
- Extract the `spatial.tar.gz` file in the same directory:
    ```bash
    tar -xzf CytAssist_FFPE_Protein_Expression_Human_Breast_Cancer_spatial.tar.gz
    ```