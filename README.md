# Transcriptional Regulator Identification Using Prize-Collecting Steiner Trees (TRIPS)

TRIPS is a two-step workflow for identifying disease-associated candidate regulators that combines 1) an active module identification step on a protein-protein interaction (PPI) network and 2) transcriptional module identification by solving the prize-collecting Steiner tree (PCST) problem on a template gene regulatory network (GRN). The input to TRIPS are:

1) Differential expression analysis results. See folder 'data' for examples.
2) PPI network (a filtered BIOGRID network is provided as default).
3) Gene regulatory network (a filtered DOROTHEA is provided as default).

# Installation

Install conda environment as follows (there also exists a environment.yml but it contains more packages than necessary):
```bash
conda create --name trips_env python=3.8
conda activate trips_env
conda install  --channel conda-forge numpy pandas networkx matplotlib pip
pip install mygene
```

## Install DOMINO
Next, follow instructions for installing DOMINO here: https://github.com/Shamir-Lab/DOMINO. Please pay attention if pcst-fast package was succesfully installed, which normally is the case for Python3.7 or lower. Otherwise install it with python setup.py install from here https://github.com/fraenkel-lab/pcst_fast.

## Install DAPCSTP
Finally, follow the instructions for installation here: https://github.com/mluipersbeck/dapcstp.

# Running TRIPS
You can run TRIPS by calling
```bash
python run_trips.py file_ppi file_grn file_degs path_domino path_dapcstp output_folder 1.0 50
```
The positional arguments are:
```
[1] Path to PPI network file. Should be a tab-separated file with two columns: "node1" and "node2".
[2] Path to gene regulatory network file. Should be a tab-separated file with two columns: "node1" and "node2".
[3] Path to the differential expression results. Should be a tab-separated file with columns "Gene_symbol", "Log_FoldChange", "AdjPValue".
[4] Path to DOMINO.
[5] Path to DAPCSTP.
[6] Path to output folder.
[7] The logFoldChange threshold.
[8] Percentile value to determine the edge cost.
```
