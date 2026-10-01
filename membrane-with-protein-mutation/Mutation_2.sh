# Input: 
# membrane.pdb
# MakeNewMutant_Modeller.py, charmmlipid2amber.py, gro/PA_PC_OL.gro
# minim.mdp, NVT.mdp. NPT.mdp, md.mdp
# ipid17.itp, amber14sb_OL15.ff
# Mutation.sh, MD_Simulation.sh

# Example:
# sh Mutation_2.sh -i membrane.pdb -c H -r 2 -t LEU
# -i: input PDB file
# -c: chain id
# -r: residue id
# -t: target amino acid

set -euo pipefail

while getopts "i:c:r:t:" opt; do
  case "$opt" in
    i) INPUT_PDB="$OPTARG" ;;
    c) CHAIN_ID="$OPTARG" ;;
    r) RESIDUE_ID="$OPTARG" ;;
    t) NEW_RES="$OPTARG" ;;
  esac
done

# 1. Extract protein_only.pdb and no_protein_only.pdb
awk '
BEGIN{
  protein_file = "protein_only.pdb"
  nonprot_file = "no_protein_only.pdb"

  nonprot["POPC"]
  nonprot["TP3"]
  nonprot["K+"]
  nonprot["Cl-"]
}

/^(ATOM|HETATM)/{
  res = substr($0,18,4)
  gsub(/ /,"",res)

  if(res in nonprot)
    print > nonprot_file
  else
    print > protein_file
}

END{
  print "END" > protein_file
  print "END" > nonprot_file
}
' "$INPUT_PDB"

echo "=================================================="
echo "File Extraction: Done!"
echo "=================================================="

# 2.1 NXXX/CXXX --> XXX
awk '
BEGIN{FS=OFS=""}

/^(ATOM|HETATM)/{
  res = substr($0,18,4)
  r = res
  gsub(/ /,"",r)

  if(length(r)==4 && (substr(r,1,1)=="N" || substr(r,1,1)=="C")){
    new = substr(r,2,3)
    $0 = substr($0,1,17) sprintf("%3s",new) " " substr($0,22)
  }
}

{print}
' protein_only.pdb > protein_only_modeller_ready.pdb

# 2.2 Mutation
python3 MakeNewMutant_Modeller.py protein_only_modeller_ready.pdb -rl "${RESIDUE_ID}:${CHAIN_ID}" -nl "${NEW_RES}" -o "protein_mutated_modeller_finished" -v

# 2.3 XXX --> NXXX/CXXX
awk '
BEGIN{FS=OFS=""}

FNR==NR{
  if($0 ~ /^(ATOM|HETATM)/){
    res = substr($0,18,4)
    r = res
    gsub(/ /,"",r)

    if(r ~ /^[NC][A-Z][A-Z][A-Z]$/){
      chain = substr($0,22,1)
      resid = substr($0,23,4)
      icode = substr($0,27,1)
      key = chain resid icode

      term[key] = substr(r,1,1)
    }
  }
  next
}

/^(ATOM|HETATM)/{
  chain = substr($0,22,1)
  resid = substr($0,23,4)
  icode = substr($0,27,1)
  key = chain resid icode

  if(key in term){
    res = substr($0,18,4)
    r = res
    gsub(/ /,"",r)

    if(length(r)==3){
      new = term[key] r
      $0 = substr($0,1,17) sprintf("%4s",new) substr($0,22)
    }
  }
}

{print}
' protein_only.pdb protein_mutated_modeller_finished.pdb > protein_mutated.pdb

echo "=================================================="
echo "Mutation setting"
echo "Chain ID     : ${CHAIN_ID}"
echo "Residue ID   : ${RESIDUE_ID}"
echo "New residue  : ${NEW_RES}"
echo "=================================================="

echo "=================================================="
echo "Protein Mutation: Done!"
echo "=================================================="

# 3.1 Rename water/ion
awk '
BEGIN{FS=OFS=""}

/^(ATOM|HETATM)/{
  res  = substr($0,18,4)
  atom = substr($0,13,4)

  r = res
  a = atom
  gsub(/ /,"",r)
  gsub(/ /,"",a)

  if(r=="TP3"){
    $0 = substr($0,1,17) sprintf("%3s","SOL") substr($0,21)

    if(a=="O")
      $0 = substr($0,1,12) sprintf("%4s","OW") substr($0,17)
    else if(a=="H1")
      $0 = substr($0,1,12) sprintf("%4s","HW1") substr($0,17)
    else if(a=="H2")
      $0 = substr($0,1,12) sprintf("%4s","HW2") substr($0,17)
  }

  else if(r=="K+"){
    $0 = substr($0,1,17) sprintf("%3s","K") substr($0,21)
    $0 = substr($0,1,12) sprintf("%4s","K") substr($0,17)
  }

  else if(r=="Cl-"){
    $0 = substr($0,1,17) sprintf("%3s","CL") substr($0,21)
    $0 = substr($0,1,12) sprintf("%4s","CL") substr($0,17)
  }
}

{print}
' no_protein_only.pdb > no_protein_only_tmp.pdb && mv no_protein_only_tmp.pdb no_protein_only.pdb

# 3.2 Amber lipid construction
python3 charmmlipid2amber.py -f no_protein_only.pdb -o no_protein_only_amber_lipid.gro -lipids POPC

# 3.3 pdb2gmx
ml load gromacs/2023.2-gcc-9.5.0-jzxesel
gmx_mpi pdb2gmx -f protein_mutated.pdb -o protein_mutated.gro -p topol.top -ff amber14sb_OL15 -water tip3p -ignh

# 3.4 Combine structures
awk -v pdb="$INPUT_PDB" '
BEGIN{
  protein = "protein_mutated.gro"
  nonprot = "no_protein_only_amber_lipid.gro"
  out = "membrane_amber_lipid.gro"

  while((getline line < pdb) > 0){
    if(substr(line,1,6)=="CRYST1"){
      a = substr(line,7,9) + 0
      b = substr(line,16,9) + 0
      c = substr(line,25,9) + 0
      box = sprintf("%10.5f%10.5f%10.5f", a/10, b/10, c/10)
      break
    }
  }
  close(pdb)

  read_gro(protein)
  read_gro(nonprot)

  print "Mutated protein + AMBER lipid17 membrane system" > out
  printf "%5d\n", total >> out

  for(i=1; i<=total; i++){
    printf "%s%5d%s\n", substr(atoms[i],1,15), i%100000, substr(atoms[i],21) >> out
  }

  print box >> out

  print "Total atoms: " total
  print "Box: " box
}

function read_gro(file, title, n, atom, i){
  getline title < file
  getline n < file
  n += 0

  for(i=1; i<=n; i++){
    getline atom < file
    atoms[++total] = atom
  }

  close(file)
}
'

# 3.5 Add water/ion/lipid into topol.top
awk -v gro="membrane_amber_lipid.gro" '
BEGIN{
  while((getline line < gro) > 0){
    nline++
    gro_line[nline] = line
  }
  close(gro)

  natoms = gro_line[2] + 0

  for(i=3; i<=2+natoms; i++){
    res = substr(gro_line[i],6,5)
    gsub(/ /,"",res)
    count[res]++
  }

  lipid = count["PA"] / 46
  sol   = count["SOL"] / 3
  kion  = count["K"] + 0
  clion = count["CL"] + 0

  printf "%-18s %d\n", "PA_PC_OL", lipid
  printf "%-18s %d\n", "K",        kion
  printf "%-18s %d\n", "CL",       clion
  printf "%-18s %d\n", "SOL",      sol
}
' >> topol.top

# 3.6 Add lipid17.itp into topol.top
awk '/forcefield.itp/{print; print "#include \"lipid17.itp\""; next} {print}' topol.top > tmp.top && mv tmp.top topol.top

# 3.7 Generate one_SOL.gro for neutralize_mut_system.py
awk '
NR==1 {
next
}
NR==2 {
  natoms = $1
  next
}
NR>=3 && NR<3+natoms {
  res = substr($0,6,5)
  gsub(/ /,"",res)
  if(res=="SOL" && count < 3){
    water[++count] = $0
  }
}
END {
  if(count != 3){
    print "ERROR: cannot extract one complete SOL molecule from membrane_amber_lipid.gro" > "/dev/stderr"
    exit 1
  }
  print "One SOL water molecule"
  print 3
  for(i=1; i<=3; i++){
    print water[i]
  }
  print "   1.00000   1.00000   1.00000"
}
' membrane_amber_lipid.gro > one_SOL.gro

echo "=================================================="
echo "MD Preparation: Done!"
echo "=================================================="

# 4 MD simulation
# ion equilibrium
python3 neutralize_mut_system.py \
  -c membrane_amber_lipid.gro \
  -p topol.top \
  -f minim.mdp \
  -w one_SOL.gro \
  -g gmx_mpi \
  -o membrane_amber_lipid_neutral.gro \
  --pname K \
  --nname CL

# sbatch MD_Simulation_no_ion_equilibrium.sh

echo "=================================================="
echo "MD Simulation: Start!"
echo "=================================================="