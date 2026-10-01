import argparse
import subprocess
import random
import shutil
import os
from pathlib import Path
import re
import logging
import MDAnalysis as mda
from MDAnalysis.core.groups import AtomGroup


out_file = "neutralize_system.out"

def neutralize_system(grofile, topfile, mdp_file, wgro_1mol, gmx_path, out_gro="membrane_amber_lipid_neutral.gro", pname="K", nname="CL"):
    """
    Neutralize the system according to charge
    - starting_system: GRO file from last step
    - topology: Topology file for the system
    - mdp_file: MDP file
    - wgro_1mol: GRO file for 1 molecule of water
    - out_gro: output GRO file
    - pname: K
    - nname: CL
    """

    # grofile = f"{starting_system}.gro"
    # topfile = f"{topology}.top"
    tmp_tpr = "system_tmp.tpr"

    # Check file readability
    if not os.path.exists(grofile):
        print(f"ERROR: I cannot read SYSTEM_file({grofile})")
        return 
    if not os.path.exists(topfile):
        print(f"ERROR: I cannot read TOPOLOGY_file({topfile})")
        return 

    # --- helpers --------------------------------------------------------------
    def parse_charge_from_grompp(stdout, stderr):
        combined = (stdout or "") + "\n" + (stderr or "")
        for line in combined.splitlines():
            line = line.strip()
            if line.startswith("Total charge") or "System has non-zero total charge" in line:
                m = re.search(r"([-+]?\d*\.\d+)", line)
                if m:
                    return float(m.group(1))
        return 0.0  # assume neutral if not found

    def load_gro(fname):
        with open(fname) as f:
            lines = f.readlines()
            header = lines[0]
            atoms = lines[2:-1]
            box = lines[-1]
            return header, atoms, box

    def write_gro(fname, header, atoms, box):
        with open(fname, "w") as f:
            f.write(header)
            f.write(f"{len(atoms)}\n")
            f.writelines(atoms)
            f.write(box)

    def update_topology_changes(changes):
        out = []
        with open(topfile) as f:
            for L in f:
                words = L.split()
                if len(words) == 2 and not L.lstrip().startswith(";") and words[0] in changes:
                    old = words[1]
                    new = str(int(old) + changes[words[0]])
                    i = L.rfind(old)
                    L = L[:i] + new + L[i+len(old):]
                out.append(L)
        with open(topfile, "w") as f:
            f.writelines(out)  


    def get_water_label():
        with open(topfile) as f:
            for L in f:
                if "TIP3" in L.split(): return "TIP3"
                if "TP3" in L.split(): return "TP3"
                if "SOL"  in L.split(): return "SOL"
        raise RuntimeError("No water found in topology.")


    # --- 1) Compute system charge ----------------------------------------
    cmd = f"{gmx_path} grompp -f {mdp_file} -c {grofile} -p {topfile} -o {tmp_tpr} -maxwarn 2"

    try:
        # Capture output for parsing
        proc = subprocess.run(cmd, shell=True, check=False, capture_output=True, text=True)
        # Write output to log file
        with open(out_file, 'a') as log:
            log.write(proc.stdout)
            log.write(proc.stderr)
        # Parse the charge from captured output
        charge = int(round(parse_charge_from_grompp(proc.stdout, proc.stderr)))
        logging.info(f"Detected system charge: {charge}")
    except subprocess.CalledProcessError:
        print("Something went wrong on gmx grompp!")
        return


    # --- 2) act based on charge -----------------------------------------------
    tmp_gro = "system_tmp.gro"

    # case 0: neutral
    if charge == 0:
        logging.info("System already neutral. Copying input to output.")
        shutil.copyfile(grofile, tmp_gro)
        Path(out_gro).unlink(missing_ok=True)  # remove old if exists
        os.rename(tmp_gro, out_gro)
        Path(tmp_tpr).unlink(missing_ok=True)
        return
    
    absq = abs(charge)

    # case ±2: replace one ion with opposite
    if absq == 2:
        ion_out = pname if charge > 0 else nname
        ion_in  = nname if charge > 0 else pname

        # Load the original gro to preserve formatting
        header, atoms, box = load_gro(grofile)
        u = mda.Universe(grofile)
        # Select residues of the ion to remove
        ions_out = u.select_atoms(f"resname {ion_out}").residues
        if len(ions_out) == 0:
            raise RuntimeError(f"No {ion_out} ions found to replace.")
        # Choose the residue: last K+ or first Cl-
        if charge > 0:
            residue_to_replace = ions_out[-1]  # last K+
        else:
            residue_to_replace = ions_out[0]   # first Cl-
        atom_numbers = [atom.index + 1 for atom in residue_to_replace.atoms]
        for atom in residue_to_replace.atoms:
            atom.residue.resname = ion_in
            atom.name = ion_in
        # Write GRO
        u.atoms.write(tmp_gro)
        header2, atoms2, box2 = load_gro(tmp_gro)
        write_gro(out_gro, header2, atoms2, box)

        # Update topology counts
        update_topology_changes({ion_out: -1, ion_in: +1})
        
        # cleanup
        Path(tmp_gro).unlink(missing_ok=True)
        Path(tmp_tpr).unlink(missing_ok=True)
        logging.info(
            f"Replaced 1 {ion_out} ion with 1 {ion_in}. "
            f"Atom numbers: {list(atom_numbers)}. "
            f"Output GRO: {out_gro}, topology: {topfile}"
        )
        return
    

    # case ±1: random method
    if absq == 1:
        method = random.choice(["genion", "genwater"])
        # water info
        wlabel = get_water_label()

        # method A: genion
        header, atoms, box = load_gro(grofile)
        if method == "genion":
            gcmd = f'echo "{wlabel}" | {gmx_path} genion -s {tmp_tpr} -o {tmp_gro} -pname {pname} -nname {nname} -neutral'
            try:
                with open(out_file, 'a') as log:
                    subprocess.run(gcmd, shell=True, check=True, stdout=log, stderr=subprocess.STDOUT)
            except subprocess.CalledProcessError:
                print("Something went wrong on genion!")
                return

            # --- move last line before box to its ion block ---
            header1, atoms1, box1 = load_gro(tmp_gro)

            new_ion_line = atoms1[-1]
            new_ion_resname = new_ion_line[5:10].strip()
            atoms1 = atoms1[:-1]

            block_indices = [i for i, line in enumerate(atoms1) if line[5:10].strip() == new_ion_resname]
            insert_index = block_indices[-1] + 1 if block_indices else len(atoms1)

            atoms1.insert(insert_index, new_ion_line)

            write_gro(tmp_gro, header1, atoms1, box)

            # --- final formatting using editconf ---
            editconf_cmd = (f"{gmx_path} editconf -f {tmp_gro} -o {out_gro}")

            try:
                with open(out_file, 'a') as log:
                    subprocess.run(editconf_cmd, shell=True, check=True, stdout=log, stderr=subprocess.STDOUT)
            except subprocess.CalledProcessError:
                print("Something went wrong on gmx editconf!")
                return

            header2, atoms2, box2 = load_gro(out_gro)
            write_gro(out_gro, header2, atoms2, box)

            update_topology_changes({wlabel: -1, new_ion_resname: +1})

            # cleanup
            Path(tmp_gro).unlink(missing_ok=True)
            Path(tmp_tpr).unlink(missing_ok=True)
            logging.info(f"Genion finished. Removed 1 {wlabel} water and added 1 {new_ion_resname}. "
                        f"Output GRO: {out_gro}, topology: {topfile}")
            return

        # method B: remove ion + add water
        ion_label = pname if charge > 0 else nname
        header, atoms, box = load_gro(grofile)

        # remove one ion
        idx = [i for i, L in enumerate(atoms) if L[5:10].strip() == ion_label]
        if not idx:
            raise RuntimeError(f"No {ion_label} ions found to remove.")
        i = random.choice(idx)
        atoms1 = [L for j, L in enumerate(atoms) if j != i]

        # temp file without ion
        write_gro(tmp_gro, header, atoms1, box)


        # insert 1 water
        icmd = f"{gmx_path} insert-molecules -f {tmp_gro} -ci {wgro_1mol} -nmol 1 -try 500 -o {out_gro}"
        
        try:
            with open(out_file, 'a') as log:
                subprocess.run(icmd, shell=True, check=True, stdout=log, stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError:
            print("Something went wrong on gmx insert-molecules!")
            return 

        header2, atoms2, box2 = load_gro(out_gro)
        write_gro(out_gro, header2, atoms2, box)

        update_topology_changes({ion_label: -1, wlabel: +1})
        Path(tmp_gro).unlink(missing_ok=True)
        Path(tmp_tpr).unlink(missing_ok=True)
        logging.info(f"Genwater finished. Removed 1 {ion_label} ion and added 1 {wlabel} water. "
            f"Output GRO: {out_gro}, topology: {topfile}")
        return


    # other charges
    raise RuntimeError(f"Unsupported system charge {charge}. Expected 0, ±1, ±2.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Neutralize a membrane-protein mutant system by adding/removing/changing ions."
    )

    parser.add_argument(
        "-c", "--grofile",
        required=True,
        help="Input GRO file, e.g. membrane_amber_lipid.gro"
    )

    parser.add_argument(
        "-p", "--topfile",
        required=True,
        help="Input/output topology file, e.g. topol.top"
    )

    parser.add_argument(
        "-f", "--mdp",
        required=True,
        help="MDP file used by grompp to compute system charge, e.g. minim.mdp"
    )

    parser.add_argument(
        "-w", "--water",
        required=True,
        help="GRO file containing one water molecule, e.g. one_SOL.gro"
    )

    parser.add_argument(
        "-g", "--gmx",
        default="gmx_mpi",
        help="GROMACS command, default: gmx_mpi"
    )

    parser.add_argument(
        "-o", "--output",
        default="membrane_amber_lipid_neutral.gro",
        help="Output neutralized GRO file"
    )

    parser.add_argument(
        "--pname",
        default="K",
        help="Positive ion residue/atom name, default: K"
    )

    parser.add_argument(
        "--nname",
        default="CL",
        help="Negative ion residue/atom name, default: CL"
    )

    args = parser.parse_args()

    neutralize_system(
        grofile=args.grofile,
        topfile=args.topfile,
        mdp_file=args.mdp,
        wgro_1mol=args.water,
        gmx_path=args.gmx,
        out_gro=args.output,
        pname=args.pname,
        nname=args.nname
    )