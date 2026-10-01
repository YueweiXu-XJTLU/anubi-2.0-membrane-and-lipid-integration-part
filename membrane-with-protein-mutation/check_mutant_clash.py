import sys
import math

if len(sys.argv) != 2:
    print("Usage: python check_mutant_clash.py input.pdb")
    sys.exit(1)

pdb = sys.argv[1]

target_chain = "H"
target_resid = "2"

nonprotein_res = {
    "POPC", "PA", "PC", "OL",
    "TP3", "TIP3", "SOL", "WAT", "HOH",
    "K", "K+", "POT",
    "CL", "Cl-", "CLA"
}

target = []
nonprotein = []

with open(pdb) as f:
    for line in f:
        if not line.startswith(("ATOM", "HETATM")):
            continue

        atom = line[12:16].strip()
        resn = line[17:21].strip()
        chain = line[21].strip()
        resid = line[22:26].strip()

        if atom.startswith("H"):
            continue

        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])

        rec = (atom, resn, chain, resid, x, y, z)

        if chain == target_chain and resid == target_resid:
            target.append(rec)
        elif resn in nonprotein_res:
            nonprotein.append(rec)

def dist(a, b):
    dx = a[4] - b[4]
    dy = a[5] - b[5]
    dz = a[6] - b[6]
    return math.sqrt(dx*dx + dy*dy + dz*dz)

clash_count = 0

for t in target:
    for n in nonprotein:
        if dist(t, n) < 2.2:
            clash_count += 1

print(f"Number of mutant-vs-nonprotein heavy-atom clashes < 2.2 Å: {clash_count}")
