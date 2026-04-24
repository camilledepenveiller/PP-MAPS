# PP-MAPS: Protein-Peptide Molecular dynamics Assisted Pharmacophore Signatures


PP-MAPS tool allows to generate a heatmap of pharmacophore features of a protein-peptide complex from a molecular dynamics trajectory.


## Create environment

After cloning the PP-MAPS repository, create and activate a new environment with the following commands:
```bash
conda create -n pp_maps python=3.13
conda activate pp_maps
```

Use the `requirements.txt` file to install dependencies:
```bash
pip install -r requirements.txt
```


## Software dependency

You need to have GROMACS installed on your system and `gmx` command accessible through the `PATH`.


## Choice of pharmacophore generator

PP-MAPS tool is able to use either LigandScout (under license) or CDPKit (open source) in the worflow. You can choose your pharmacophore generator by specifying it when running PP-MAPS tool. Before that, you have to install LigandScout or CDPKit.

- For LigandScout, you have to adapt path to the pharmacophore generator: set the variable Pharmacophore_generator_path in `config.ini` file.
- For CDPKit, you have to get at least version 1.3 and to be careful to install also Python bindings (see https://cdpkit.org/installation.html). If version 1.3 is not yet available, consider compiling CDPKit from source from the GitHub master branch.


## Usage

:warning: Before providing your input files, be careful to prepare an XTC file with a centered trajectory, provided that in the original PDB file of the complex, the peptide is at the end of the file.

Run the following command to run pharmacophore analysis through MD trajectory and generate a pharmacomap:

```bash
python pp_maps.py [-h] -xtc XTC -tpr TPR [-n N] [-o OUTPUT]

options:
  -h, --help                        show this help message and exit
  -xtc XTC                          XTC trajectory file with centered system
  -tpr TPR                          TPR file used for MD
  -ligandscout                      add this argument to use LigandScout as pharmacophore generator
  -cdpkit                           add this argument to use CDPKit as pharmacophore generator
  -n N                              number of processes to perform the analysis, default to 1
  -o OUTPUT, --output OUTPUT        PNG output file with heatmap of pharmacophore features of the whole MD, default to pharmacomap.png
```

Note that when using PP-MAPS with LigandScout, it will automatically produce two additional pharmacomaps, based on the same statistics, but allowing to differentiate between backbone (`pharmacomap_BB.png`) and side-chain (`pharmacomap_SC.png`) interactions with the protein amino acids.


## Example

In the `example` folder, you can find a short trajectory XTC file and a TPR file for tool testing.

Here is the command to run for pharmacomap generation from these files:

- Using LigandScout

```bash
python pp_maps.py -xtc example/md.xtc -tpr example/md.tpr -ligandscout -n 8
```

- Using CDPKit

```bash
python pp_maps.py -xtc example/md.xtc -tpr example/md.tpr -cdpkit -n 8
```


## References

PP-MAPS tool is based on `ipharmgen` tool from LigandScout and on `gen_ia_ph4s.py` script (modified) from CDPKit.

- LigandScout

*Installation process:* https://docs.inteligand.com/ligandscout/

*Article:* G. Wolber and T. Langer, ‘LigandScout: 3-D Pharmacophores Derived from Protein-Bound Ligands and Their Use as Virtual Screening Filters’, J. Chem. Inf. Model., vol. 45, no. 1, pp. 160–169, Jan. 2005, doi: 10.1021/ci049885e.

- CDPKit

*Source code:* Thomas Seidel, Chemical Data Processing Toolkit source code repository, https://github.com/molinfo-vienna/CDPKit

*Documentation:* Thomas Seidel, Oliver Wieder, Chemical Data Processing Toolkit documentation pages, https://cdpkit.org


## Formatting

`isort .; black -l 79 .`