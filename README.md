# Pharmacomaps


Pharmacomaps tool allows to generate a heatmap of pharmacophore features between a receptor (protein) and a ligand (peptide) from a molecular dynamics trajectory.


Use `environment.yml` to create the environment:
```bash
conda env create -f environment.yml
conda activate pharmacomaps
```

**Adapt path to LigandScout:** set the variable Pharmacophore_generator_path in `config.ini` file.


## Usage

:warning: Before providing your input files, be careful to prepare a xtc file with a centered trajectory.

Run the following command to run pharmacophore analysis through MD trajectory and generate a pharmacomap:

```bash
python pharmacomaps.py [-h] -xtc XTC -tpr TPR [-n N] [-o OUTPUT]

options:
  -h, --help                        show this help message and exit
  -xtc XTC                          XTC trajectory file with centered system
  -tpr TPR                          TPR file used for MD
  -n N                              number of processes to perform the analysis, default to 1
  -o OUTPUT, --output OUTPUT        PNG output file with heatmap of pharmacophore features of the whole MD, default to pharmacomap.png
```