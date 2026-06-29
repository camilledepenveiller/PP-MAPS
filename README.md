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

For trajectory conversion from XTC to PDBs, you have two possibilities.
Either you need to have GROMACS installed on your system and `gmx` command accessible through the `PATH`.
Otherwise, you have to install MDTraj (1.11.1) in your environment through this command:
```bash
pip install mdtraj==1.11.1
```


## Choice of pharmacophore generator

PP-MAPS tool is able to use either LigandScout (under license) or CDPKit (open source) in the workflow. You can choose your pharmacophore generator by specifying it when running PP-MAPS tool. Before that, you have to install LigandScout or CDPKit.

- For LigandScout, you have to adapt path to the pharmacophore generator: set the variable Pharmacophore_generator_path in `config.ini` file.
- For CDPKit:
  * first install CDPKit (1.3.0) with your package manager (see https://cdpkit.org/installation.html);
  * then install Python bindings with `pip install cdpkit==1.3.0`.


## Usage

You can choose to provide as input either an XTC trajectory, or directly a directory of PDB files.

:warning: Before providing your input files, be careful to prepare an XTC file with a centered trajectory, provided that in the original PDB file of the complex, the peptide is at the end of the file.

Run the following command to run pharmacophore analysis through MD trajectory and generate a pharmacomap:

```bash
pp_maps.py [-h] [-xtc XTC] [-pdb PDB] -topol TOPOL [-ligandscout] [-cdpkit] [-gromacs] [-mdtraj] [-n N] [-o OUTPUT]

options:
  -h, --help           show this help message and exit
  -xtc XTC             Path to XTC trajectory file with centered system.
  -pdb PDB             Path to PDB files directory.
  -topol TOPOL         Topology file (TPR with GROMACS or PDB with MDTraj). Required when providing an XTC as input.
  -ligandscout         Add this argument to use LigandScout as pharmacophore generator.
  -cdpkit              Add this argument to use CDPKit as pharmacophore generator.
  -gromacs             Add this argument to use GROMACS to convert XTC to PDBs.
  -mdtraj              Add this argument to use MDTraj (mdconvert) to convert XTC to PDBs.
  -n N                 Number of processors/CPUs to perform the analysis. Default to 1.
  -o, --output OUTPUT  PNG output file with heatmap of pharmacophore features of the whole MD. Default to pharmacomap.png
```

Note that when using PP-MAPS with LigandScout, it will automatically produce two additional pharmacomaps, based on the same statistics, but allowing to differentiate between backbone (`pharmacomap_BB.png`) and side-chain (`pharmacomap_SC.png`) interactions with the protein amino acids.


## Example

In the `example` folder, you can find a short trajectory XTC file and a TPR file for tool testing.

Here is the command to run for pharmacomap generation from these files:

- Using LigandScout

```bash
python pp_maps.py -xtc example/md.xtc -topol example/md.tpr -ligandscout -gromacs -n 8
```

- Using CDPKit

```bash
python pp_maps.py -xtc example/md.xtc -topol example/md.tpr -cdpkit -gromacs -n 8
```


## References

PP-MAPS tool is based on `ipharmgen` tool from LigandScout and on `gen_ia_ph4s.py` script (modified) from CDPKit.

- LigandScout

*Installation process:* https://docs.inteligand.com/ligandscout/.

*Article:* G. Wolber and T. Langer, ‘LigandScout: 3-D Pharmacophores Derived from Protein-Bound Ligands and Their Use as Virtual Screening Filters’, J. Chem. Inf. Model., vol. 45, no. 1, pp. 160–169, Jan. 2005, doi: 10.1021/ci049885e.

- CDPKit

*Source code:* Thomas Seidel, Chemical Data Processing Toolkit source code repository, https://github.com/molinfo-vienna/CDPKit.

*Documentation:* Thomas Seidel, Oliver Wieder, Chemical Data Processing Toolkit documentation pages, https://cdpkit.org.

- GROMACS

*Installation process:* https://manual.gromacs.org/current/install-guide/index.html.

*Article:* Berendsen, H. J. C.; van der Spoel, D.; van Drunen, R. GROMACS: A Message-Passing Parallel Molecular Dynamics Implementation. Computer Physics Communications 1995, 91 (1–3), 43–56. https://doi.org/10.1016/0010-4655(95)00042-E.

- MDTraj

*Article:* McGibbon, R. T.; Beauchamp, K. A.; Harrigan, M. P.; Klein, C.; Swails, J. M.; Hernández, C. X.; Schwantes, C. R.; Wang, L.-P.; Lane, T. J.; Pande, V. S. MDTraj: A Modern Open Library for the Analysis of Molecular Dynamics Trajectories. Biophysical Journal 2015, 109 (8), 1528–1532. https://doi.org/10.1016/j.bpj.2015.08.015.


## Formatting

`isort .; black -l 79 .`