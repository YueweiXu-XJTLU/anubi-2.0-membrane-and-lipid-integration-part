# ANUBI 2.0 - Membrane and lipid integration part

This repository archives research work on integrating protein mutation with membrane lipids, water, and ions for Amber/GROMACS molecular dynamics. The work may become an iterative branch of ANUBI 2.0; integration into ANUBI 2.0 and end-to-end validation are not established by this archive.

## Relationship to ANUBI

[ANUBI (Buratto et al., *Journal of Chemical Theory and Computation*, 2026)](https://doi.org/10.1021/acs.jctc.5c01640) explores protein and peptide mutations and evaluates candidates using binding-affinity calculations. This repository addresses a different preparation task: carrying a selected protein mutation into a membrane system while adapting its lipid, solvent, and ion components. It does not include ANUBI's sequence-search or binding free-energy workflow.

## Tutorials and paired scripts

| Group | Retained material | Scope |
|---|---|---|
| [Membrane preparation](membrane/) | [Membrane.md](membrane/Membrane.md) | Earlier manual notes for residue naming, lipid conversion, topology, and MD preparation; no separate matching `.sh` file is retained. |
| [Membrane with protein mutation](membrane-with-protein-mutation/) | [Tutorial](membrane-with-protein-mutation/Membrane_With_Protein_Mutation.md), [Mutation.sh](membrane-with-protein-mutation/Mutation.sh), [MD_Simulation.sh](membrane-with-protein-mutation/MD_Simulation.sh) | The route explicitly paired in the tutorial; the MD script uses GROMACS `genion`. |
| [Later mutation route](membrane-with-protein-mutation/) | [Mutation_2.sh](membrane-with-protein-mutation/Mutation_2.sh), [MD_Simulation_no_ion_equilibrium.sh](membrane-with-protein-mutation/MD_Simulation_no_ion_equilibrium.sh), [neutralize_mut_system.py](membrane-with-protein-mutation/neutralize_mut_system.py) | Option-based mutation script followed by pre-MD neutralization and an MD script that starts from the neutralized system. |

The mutation group also retains [clash checking](membrane-with-protein-mutation/check_mutant_clash.py), the [project presentation](membrane-with-protein-mutation/membrane-with-protein-mutation.pptx), and the [MODELLER-derived mutation helper](membrane-with-protein-mutation/MakeNewMutant_Modeller.py). The [lipid converter](shared/charmmlipid2amber.py) is shared by the tutorials.

## External sources and workflow boundary

[MODELLER](https://salilab.org/modeller/) is the protein-structure modeling software used for the mutation step. The included `MakeNewMutant_Modeller.py` was adapted from MODELLER's official [`mutate_model.py` example](https://salilab.org/modeller/resources/mutate_model.html); it is not original code written from scratch for this project. The included `charmmlipid2amber.py` comes from [gmx_lipid17.ff](https://github.com/xiki-tempula/gmx_lipid17.ff).

The scripts refer to helpers by local filename and were written for a particular HPC environment. Run a selected route from a separate working directory containing its required helpers and input files. The original membrane-protein coordinates, MDP files, POPC template, force-field distribution, and complete validation records are not included. See [dependencies and provenance](DEPENDENCIES.md) before reuse.
