# Dependencies and provenance

## Environment and inputs

The retained workflows use Bash, GROMACS, Python, [MODELLER](https://salilab.org/modeller/), NumPy, and MDAnalysis. The archived Slurm scripts load GROMACS 2023.2; their partition, GPU, module, and mail settings are specific to the original HPC environment.

Stage a selected script pair and its helpers in a separate working directory. The scripts refer to `MakeNewMutant_Modeller.py`, `charmmlipid2amber.py`, `neutralize_mut_system.py`, and their input files by local filename; they do not discover those files from this repository's directories. Required external inputs include an appropriate membrane-protein PDB, MDP files, `gro/PA_PC_OL.gro`, `lipid17.itp`, the `amber14sb_OL15.ff` force field, and compatible lipid parameters. These inputs are not bundled here.

## Tutorial and script relationships

| Tutorial or route | Scripts and helpers |
|---|---|
| [Membrane.md](membrane/Membrane.md) | Manual precursor; no dedicated standalone `.sh` counterpart is retained. It uses the [lipid converter](shared/charmmlipid2amber.py). |
| [Membrane_With_Protein_Mutation.md](membrane-with-protein-mutation/Membrane_With_Protein_Mutation.md) | [Mutation.sh](membrane-with-protein-mutation/Mutation.sh) followed by [MD_Simulation.sh](membrane-with-protein-mutation/MD_Simulation.sh), as stated in the tutorial. |
| Later option-based route | [Mutation_2.sh](membrane-with-protein-mutation/Mutation_2.sh) calls [neutralize_mut_system.py](membrane-with-protein-mutation/neutralize_mut_system.py); [MD_Simulation_no_ion_equilibrium.sh](membrane-with-protein-mutation/MD_Simulation_no_ion_equilibrium.sh) starts from its neutralized coordinates. |

The optional [clash-check helper](membrane-with-protein-mutation/check_mutant_clash.py) is retained beside the mutation scripts. The tutorial records a residual-charge error, and this archive does not establish that either complete route has been rerun successfully after reorganization.

## Third-party code

| Retained file | Source |
|---|---|
| [MakeNewMutant_Modeller.py](membrane-with-protein-mutation/MakeNewMutant_Modeller.py) | Adapted from the official [MODELLER `mutate_model.py` example](https://salilab.org/modeller/resources/mutate_model.html); its local header also cites the MODELLER mutation page. |
| [charmmlipid2amber.py](shared/charmmlipid2amber.py) | Retained from [gmx_lipid17.ff](https://github.com/xiki-tempula/gmx_lipid17.ff). |

These links identify upstream sources; attribution does not by itself establish redistribution rights for third-party code. No blanket license is assigned to this assembled research archive.
