# input: 
# membrane.pdb
# MakeNewMutant_Modeller.py, charmmlipid2amber.py, gro/PA_PC_OL.gro
# minim.mdp, NVT.mdp. NPT.mdp, md.mdp
# ipid17.itp, amber14sb_OL15.ff
# Mutation.sh, MD_Simulation.sh

set -euo pipefail

if [ "$#" -ne 3 ]; then
  echo "Usage: bash Mutation.sh <CHAIN_ID> <RESIDUE_ID> <NEW_RES>"
  echo "Example: bash Mutation.sh H 2 LEU"
  exit 1
fi

CHAIN_ID="$1"
RESIDUE_ID="$2"
NEW_RES="$3"

# 1.1 extract protein_only.pdb
awk '
/^(ATOM|HETATM)/{
  res=substr($0,18,4)
  gsub(/ /,"",res)

  if(res!="POPC" &&
     res!="PA" && res!="PC" && res!="OL" &&
     res!="TP3" && res!="TIP3" && res!="WAT" && res!="HOH" && res!="SOL" &&
     res!="K" && res!="K+" && res!="POT" &&
     res!="CL" && res!="Cl-" && res!="CLA")
    print
}
END{print "END"}
' membrane.pdb > protein_only.pdb

# 1.2 extract no_protein_only.pdb
awk '
/^(ATOM|HETATM)/{
  res=substr($0,18,4)
  gsub(/ /,"",res)

  if(res=="POPC" ||
     res=="TP3" || res=="TIP3" || res=="WAT" || res=="HOH" || res=="SOL" ||
     res=="K" || res=="K+" || res=="POT" ||
     res=="CL" || res=="Cl-" || res=="CLA")
    print
}
' membrane.pdb > no_protein_only.pdb

echo "=================================================="
echo "File Extraction: Done!"
echo "=================================================="

# 2.1 NXXX/CXXX --> XXX
awk 'BEGIN{FS=OFS=""}
/^(ATOM|HETATM)/{
  res=substr($0,18,4)
  r=res
  gsub(/ /,"",r)

  if(r=="NASP") new="ASP"
  else if(r=="NGLU") new="GLU"
  else if(r=="NMET") new="MET"
  else if(r=="NALA") new="ALA"
  else if(r=="NARG") new="ARG"
  else if(r=="NASN") new="ASN"
  else if(r=="NCYS") new="CYS"
  else if(r=="NGLN") new="GLN"
  else if(r=="NGLY") new="GLY"
  else if(r=="NHID" || r=="NHIE" || r=="NHIP") new=substr(r,2,3)
  else if(r=="NILE") new="ILE"
  else if(r=="NLEU") new="LEU"
  else if(r=="NLYS") new="LYS"
  else if(r=="NPHE") new="PHE"
  else if(r=="NPRO") new="PRO"
  else if(r=="NSER") new="SER"
  else if(r=="NTHR") new="THR"
  else if(r=="NTRP") new="TRP"
  else if(r=="NTYR") new="TYR"
  else if(r=="NVAL") new="VAL"

  else if(r=="CALA") new="ALA"
  else if(r=="CARG") new="ARG"
  else if(r=="CASN") new="ASN"
  else if(r=="CASP") new="ASP"
  else if(r=="CCYS") new="CYS"
  else if(r=="CGLN") new="GLN"
  else if(r=="CGLU") new="GLU"
  else if(r=="CGLY") new="GLY"
  else if(r=="CHID" || r=="CHIE" || r=="CHIP") new=substr(r,2,3)
  else if(r=="CILE") new="ILE"
  else if(r=="CLEU") new="LEU"
  else if(r=="CLYS") new="LYS"
  else if(r=="CMET") new="MET"
  else if(r=="CPHE") new="PHE"
  else if(r=="CPRO") new="PRO"
  else if(r=="CSER") new="SER"
  else if(r=="CTHR") new="THR"
  else if(r=="CTRP") new="TRP"
  else if(r=="CTYR") new="TYR"
  else if(r=="CVAL") new="VAL"
  else new=""

  if(new!=""){
    $0=substr($0,1,17) sprintf("%3s",new) " " substr($0,22)
  }
}
{print}
' protein_only.pdb > protein_only_modeller_ready.pdb

# 2.2 Mutation
python3 MakeNewMutant_Modeller.py protein_only_modeller_ready.pdb -rl "${RESIDUE_ID}:${CHAIN_ID}" -nl "${NEW_RES}" -o "protein_mutated_modeller_finished" -v

# 2.3 NXXX/CXXX --> XXX
awk '
BEGIN{FS=OFS=""}

FNR==NR{
  if($0 ~ /^(ATOM|HETATM)/){
    res=substr($0,18,4)
    r=res
    gsub(/ /,"",r)

    if(r=="NASP" || r=="NGLU" || r=="NMET" ||
       r=="NALA" || r=="NARG" || r=="NASN" || r=="NCYS" || r=="NGLN" ||
       r=="NGLY" || r=="NHID" || r=="NHIE" || r=="NHIP" || r=="NILE" ||
       r=="NLEU" || r=="NLYS" || r=="NPHE" || r=="NPRO" || r=="NSER" ||
       r=="NTHR" || r=="NTRP" || r=="NTYR" || r=="NVAL" ||

       r=="CALA" || r=="CARG" || r=="CASN" || r=="CASP" || r=="CCYS" ||
       r=="CGLN" || r=="CGLU" || r=="CGLY" || r=="CHID" || r=="CHIE" ||
       r=="CHIP" || r=="CILE" || r=="CLEU" || r=="CLYS" || r=="CMET" ||
       r=="CPHE" || r=="CPRO" || r=="CSER" || r=="CTHR" || r=="CTRP" ||
       r=="CTYR" || r=="CVAL"){
       
       chain=substr($0,22,1)
       resid=substr($0,23,4)
       icode=substr($0,27,1)
       key=chain resid icode
       term[key]=r
    }
  }
  next
}

{
  if($0 ~ /^(ATOM|HETATM)/){
    chain=substr($0,22,1)
    resid=substr($0,23,4)
    icode=substr($0,27,1)
    key=chain resid icode

    if(key in term){
      $0=substr($0,1,17) sprintf("%4s",term[key]) substr($0,22)
    }
  }
  print
}
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

# 3.1 Rename protein
awk 'BEGIN{FS=OFS=""} /^(ATOM|HETATM)/{new=""; res=substr($0,18,4); gsub(/ /,"",res); atom=substr($0,13,4); a=atom; gsub(/ /,"",a); if(res=="ARG"){if(a=="HB2") new="HB1"; else if(a=="HB3") new="HB2"; else if(a=="HG2") new="HG1"; else if(a=="HG3") new="HG2"; else if(a=="HD2") new="HD1"; else if(a=="HD3") new="HD2"; if(new!="") $0=substr($0,1,12) sprintf("%4s",new) substr($0,17)}} {print}' protein_mutated.pdb > protein_mutated_tmp.pdb && mv protein_mutated_tmp.pdb protein_mutated.pdb

awk 'BEGIN{FS=OFS=""} /^(ATOM|HETATM)/{new=""; res=substr($0,18,4); gsub(/ /,"",res); atom=substr($0,13,4); a=atom; gsub(/ /,"",a); if(res=="ASN"){if(a=="HB2") new="HB1"; else if(a=="HB3") new="HB2"; if(new!="") $0=substr($0,1,12) sprintf("%4s",new) substr($0,17)}} {print}' protein_mutated.pdb > protein_mutated_tmp.pdb && mv protein_mutated_tmp.pdb protein_mutated.pdb

awk 'BEGIN{FS=OFS=""} /^(ATOM|HETATM)/{new=""; res=substr($0,18,4); gsub(/ /,"",res); atom=substr($0,13,4); a=atom; gsub(/ /,"",a); if(res=="ASP"){if(a=="HB2") new="HB1"; else if(a=="HB3") new="HB2"; if(new!="") $0=substr($0,1,12) sprintf("%4s",new) substr($0,17)}} {print}' protein_mutated.pdb > protein_mutated_tmp.pdb && mv protein_mutated_tmp.pdb protein_mutated.pdb

awk 'BEGIN{FS=OFS=""} /^(ATOM|HETATM)/{new=""; res=substr($0,18,4); gsub(/ /,"",res); atom=substr($0,13,4); a=atom; gsub(/ /,"",a); if(res=="CALA"){if(a=="O") new="OC1"; else if(a=="OXT") new="OC2"; if(new!="") $0=substr($0,1,12) sprintf("%4s",new) substr($0,17)}} {print}' protein_mutated.pdb > protein_mutated_tmp.pdb && mv protein_mutated_tmp.pdb protein_mutated.pdb

awk 'BEGIN{FS=OFS=""} /^(ATOM|HETATM)/{new=""; res=substr($0,18,4); gsub(/ /,"",res); atom=substr($0,13,4); a=atom; gsub(/ /,"",a); if(res=="CLYS"){if(a=="HB2") new="HB1"; else if(a=="HB3") new="HB2"; else if(a=="HG2") new="HG1"; else if(a=="HG3") new="HG2"; else if(a=="HD2") new="HD1"; else if(a=="HD3") new="HD2"; else if(a=="HE2") new="HE1"; else if(a=="HE3") new="HE2"; else if(a=="O") new="OC1"; else if(a=="OXT") new="OC2"; if(new!="") $0=substr($0,1,12) sprintf("%4s",new) substr($0,17)}} {print}' protein_mutated.pdb > protein_mutated_tmp.pdb && mv protein_mutated_tmp.pdb protein_mutated.pdb

awk 'BEGIN{FS=OFS=""} /^(ATOM|HETATM)/{new=""; res=substr($0,18,4); gsub(/ /,"",res); atom=substr($0,13,4); a=atom; gsub(/ /,"",a); if(res=="CSER"){if(a=="HB2") new="HB1"; else if(a=="HB3") new="HB2"; else if(a=="O") new="OC1"; else if(a=="OXT") new="OC2"; if(new!="") $0=substr($0,1,12) sprintf("%4s",new) substr($0,17)}} {print}' protein_mutated.pdb > protein_mutated_tmp.pdb && mv protein_mutated_tmp.pdb protein_mutated.pdb

awk 'BEGIN{FS=OFS=""} /^(ATOM|HETATM)/{new=""; res=substr($0,18,4); gsub(/ /,"",res); atom=substr($0,13,4); a=atom; gsub(/ /,"",a); if(res=="CYS"){if(a=="HB2") new="HB1"; else if(a=="HB3") new="HB2"; if(new!="") $0=substr($0,1,12) sprintf("%4s",new) substr($0,17)}} {print}' protein_mutated.pdb > protein_mutated_tmp.pdb && mv protein_mutated_tmp.pdb protein_mutated.pdb

awk 'BEGIN{FS=OFS=""} /^(ATOM|HETATM)/{new=""; res=substr($0,18,4); gsub(/ /,"",res); atom=substr($0,13,4); a=atom; gsub(/ /,"",a); if(res=="GLN"){if(a=="HB2") new="HB1"; else if(a=="HB3") new="HB2"; else if(a=="HG2") new="HG1"; else if(a=="HG3") new="HG2"; if(new!="") $0=substr($0,1,12) sprintf("%4s",new) substr($0,17)}} {print}' protein_mutated.pdb > protein_mutated_tmp.pdb && mv protein_mutated_tmp.pdb protein_mutated.pdb

awk 'BEGIN{FS=OFS=""} /^(ATOM|HETATM)/{new=""; res=substr($0,18,4); gsub(/ /,"",res); atom=substr($0,13,4); a=atom; gsub(/ /,"",a); if(res=="GLU"){if(a=="HB2") new="HB1"; else if(a=="HB3") new="HB2"; else if(a=="HG2") new="HG1"; else if(a=="HG3") new="HG2"; if(new!="") $0=substr($0,1,12) sprintf("%4s",new) substr($0,17)}} {print}' protein_mutated.pdb > protein_mutated_tmp.pdb && mv protein_mutated_tmp.pdb protein_mutated.pdb

awk 'BEGIN{FS=OFS=""} /^(ATOM|HETATM)/{new=""; res=substr($0,18,4); gsub(/ /,"",res); atom=substr($0,13,4); a=atom; gsub(/ /,"",a); if(res=="GLY"){if(a=="HA2") new="HA1"; else if(a=="HA3") new="HA2"; if(new!="") $0=substr($0,1,12) sprintf("%4s",new) substr($0,17)}} {print}' protein_mutated.pdb > protein_mutated_tmp.pdb && mv protein_mutated_tmp.pdb protein_mutated.pdb

awk 'BEGIN{FS=OFS=""} /^(ATOM|HETATM)/{new=""; res=substr($0,18,4); gsub(/ /,"",res); atom=substr($0,13,4); a=atom; gsub(/ /,"",a); if(res=="HID"){if(a=="HB2") new="HB1"; else if(a=="HB3") new="HB2"; if(new!="") $0=substr($0,1,12) sprintf("%4s",new) substr($0,17)}} {print}' protein_mutated.pdb > protein_mutated_tmp.pdb && mv protein_mutated_tmp.pdb protein_mutated.pdb

awk 'BEGIN{FS=OFS=""} /^(ATOM|HETATM)/{new=""; res=substr($0,18,4); gsub(/ /,"",res); atom=substr($0,13,4); a=atom; gsub(/ /,"",a); if(res=="ILE"){if(a=="CD1") new="CD"; else if(a=="HD11") new="HD1"; else if(a=="HD12") new="HD2"; else if(a=="HD13") new="HD3"; else if(a=="HG12") new="HG11"; else if(a=="HG13") new="HG12"; if(new!="") $0=substr($0,1,12) sprintf("%4s",new) substr($0,17)}} {print}' protein_mutated.pdb > protein_mutated_tmp.pdb && mv protein_mutated_tmp.pdb protein_mutated.pdb

awk 'BEGIN{FS=OFS=""} /^(ATOM|HETATM)/{new=""; res=substr($0,18,4); gsub(/ /,"",res); atom=substr($0,13,4); a=atom; gsub(/ /,"",a); if(res=="LEU"){if(a=="HB2") new="HB1"; else if(a=="HB3") new="HB2"; if(new!="") $0=substr($0,1,12) sprintf("%4s",new) substr($0,17)}} {print}' protein_mutated.pdb > protein_mutated_tmp.pdb && mv protein_mutated_tmp.pdb protein_mutated.pdb

awk 'BEGIN{FS=OFS=""} /^(ATOM|HETATM)/{new=""; res=substr($0,18,4); gsub(/ /,"",res); atom=substr($0,13,4); a=atom; gsub(/ /,"",a); if(res=="LYS"){if(a=="HB2") new="HB1"; else if(a=="HB3") new="HB2"; else if(a=="HG2") new="HG1"; else if(a=="HG3") new="HG2"; else if(a=="HD2") new="HD1"; else if(a=="HD3") new="HD2"; else if(a=="HE2") new="HE1"; else if(a=="HE3") new="HE2"; if(new!="") $0=substr($0,1,12) sprintf("%4s",new) substr($0,17)}} {print}' protein_mutated.pdb > protein_mutated_tmp.pdb && mv protein_mutated_tmp.pdb protein_mutated.pdb

awk 'BEGIN{FS=OFS=""} /^(ATOM|HETATM)/{new=""; res=substr($0,18,4); gsub(/ /,"",res); atom=substr($0,13,4); a=atom; gsub(/ /,"",a); if(res=="MET"){if(a=="HB2") new="HB1"; else if(a=="HB3") new="HB2"; else if(a=="HG2") new="HG1"; else if(a=="HG3") new="HG2"; if(new!="") $0=substr($0,1,12) sprintf("%4s",new) substr($0,17)}} {print}' protein_mutated.pdb > protein_mutated_tmp.pdb && mv protein_mutated_tmp.pdb protein_mutated.pdb

awk 'BEGIN{FS=OFS=""} /^(ATOM|HETATM)/{new=""; res=substr($0,18,4); gsub(/ /,"",res); atom=substr($0,13,4); a=atom; gsub(/ /,"",a); if(res=="NASP"){if(a=="HB2") new="HB1"; else if(a=="HB3") new="HB2"; if(new!="") $0=substr($0,1,12) sprintf("%4s",new) substr($0,17)}} {print}' protein_mutated.pdb > protein_mutated_tmp.pdb && mv protein_mutated_tmp.pdb protein_mutated.pdb

awk 'BEGIN{FS=OFS=""} /^(ATOM|HETATM)/{new=""; res=substr($0,18,4); gsub(/ /,"",res); atom=substr($0,13,4); a=atom; gsub(/ /,"",a); if(res=="NGLU"){if(a=="HB2") new="HB1"; else if(a=="HB3") new="HB2"; else if(a=="HG2") new="HG1"; else if(a=="HG3") new="HG2"; if(new!="") $0=substr($0,1,12) sprintf("%4s",new) substr($0,17)}} {print}' protein_mutated.pdb > protein_mutated_tmp.pdb && mv protein_mutated_tmp.pdb protein_mutated.pdb

awk 'BEGIN{FS=OFS=""} /^(ATOM|HETATM)/{new=""; res=substr($0,18,4); gsub(/ /,"",res); atom=substr($0,13,4); a=atom; gsub(/ /,"",a); if(res=="NMET"){if(a=="HB2") new="HB1"; else if(a=="HB3") new="HB2"; else if(a=="HG2") new="HG1"; else if(a=="HG3") new="HG2"; if(new!="") $0=substr($0,1,12) sprintf("%4s",new) substr($0,17)}} {print}' protein_mutated.pdb > protein_mutated_tmp.pdb && mv protein_mutated_tmp.pdb protein_mutated.pdb

awk 'BEGIN{FS=OFS=""} /^(ATOM|HETATM)/{new=""; res=substr($0,18,4); gsub(/ /,"",res); atom=substr($0,13,4); a=atom; gsub(/ /,"",a); if(res=="PHE"){if(a=="HB2") new="HB1"; else if(a=="HB3") new="HB2"; if(new!="") $0=substr($0,1,12) sprintf("%4s",new) substr($0,17)}} {print}' protein_mutated.pdb > protein_mutated_tmp.pdb && mv protein_mutated_tmp.pdb protein_mutated.pdb

awk 'BEGIN{FS=OFS=""} /^(ATOM|HETATM)/{new=""; res=substr($0,18,4); gsub(/ /,"",res); atom=substr($0,13,4); a=atom; gsub(/ /,"",a); if(res=="PRO"){if(a=="HB2") new="HB1"; else if(a=="HB3") new="HB2"; else if(a=="HG2") new="HG1"; else if(a=="HG3") new="HG2"; else if(a=="HD2") new="HD1"; else if(a=="HD3") new="HD2"; if(new!="") $0=substr($0,1,12) sprintf("%4s",new) substr($0,17)}} {print}' protein_mutated.pdb > protein_mutated_tmp.pdb && mv protein_mutated_tmp.pdb protein_mutated.pdb

awk 'BEGIN{FS=OFS=""} /^(ATOM|HETATM)/{new=""; res=substr($0,18,4); gsub(/ /,"",res); atom=substr($0,13,4); a=atom; gsub(/ /,"",a); if(res=="SER"){if(a=="HB2") new="HB1"; else if(a=="HB3") new="HB2"; if(new!="") $0=substr($0,1,12) sprintf("%4s",new) substr($0,17)}} {print}' protein_mutated.pdb > protein_mutated_tmp.pdb && mv protein_mutated_tmp.pdb protein_mutated.pdb

awk 'BEGIN{FS=OFS=""} /^(ATOM|HETATM)/{new=""; res=substr($0,18,4); gsub(/ /,"",res); atom=substr($0,13,4); a=atom; gsub(/ /,"",a); if(res=="TRP"){if(a=="HB2") new="HB1"; else if(a=="HB3") new="HB2"; if(new!="") $0=substr($0,1,12) sprintf("%4s",new) substr($0,17)}} {print}' protein_mutated.pdb > protein_mutated_tmp.pdb && mv protein_mutated_tmp.pdb protein_mutated.pdb

awk 'BEGIN{FS=OFS=""} /^(ATOM|HETATM)/{new=""; res=substr($0,18,4); gsub(/ /,"",res); atom=substr($0,13,4); a=atom; gsub(/ /,"",a); if(res=="TYR"){if(a=="HB2") new="HB1"; else if(a=="HB3") new="HB2"; if(new!="") $0=substr($0,1,12) sprintf("%4s",new) substr($0,17)}} {print}' protein_mutated.pdb > protein_mutated_tmp.pdb && mv protein_mutated_tmp.pdb protein_mutated.pdb

# 3.2 Rename water/ion
awk 'BEGIN{FS=OFS=""}
/^(ATOM|HETATM)/{
  res=substr($0,18,3)
  atom=substr($0,13,4)
  r=res
  a=atom
  gsub(/ /,"",r)
  gsub(/ /,"",a)
  new_res=""
  new_atom=""
  if(r=="TP3" || r=="TIP3" || r=="WAT" || r=="HOH"){
    new_res="SOL"
    if(a=="O" || a=="OH2") new_atom="OW"
    else if(a=="H1" || a=="H") new_atom="HW1"
    else if(a=="H2") new_atom="HW2"
  }
  else if(r=="K+" || r=="POT" || r=="K"){
    new_res="K"
    if(a=="K+" || a=="POT" || a=="K") new_atom="K"
  }
  else if(r=="Cl-" || r=="CLA" || r=="CL"){
    new_res="CL"
    if(a=="Cl-" || a=="CLA" || a=="CL") new_atom="CL"
  }
  if(new_res!=""){
    $0 = substr($0,1,17) sprintf("%3s",new_res) substr($0,21)
  }
  if(new_atom!=""){
    $0 = substr($0,1,12) sprintf("%4s",new_atom) substr($0,17)
  }
}
{print}' no_protein_only.pdb > no_protein_only_tmp.pdb && mv no_protein_only_tmp.pdb no_protein_only.pdb

# 3.3 Amber lipid construction
python3 charmmlipid2amber.py -f no_protein_only.pdb -o no_protein_only_amber_lipid.gro -lipids POPC

# 3.4 pdb2gmx
ml load gromacs/2023.2-gcc-9.5.0-jzxesel
gmx_mpi pdb2gmx -f protein_mutated.pdb -o protein_mutated.gro -p topol.top -ff amber14sb_OL15 -water tip3p -ignh

# 3.5 Combine structures
awk '
BEGIN{
  protein="protein_mutated.gro"
  nonprot="no_protein_only_amber_lipid.gro"
  pdb="membrane.pdb"
  out="membrane_amber_lipid.gro"

  # read CRYST1 box from membrane.pdb, Å -> nm
  while((getline line < pdb) > 0){
    if(substr(line,1,6)=="CRYST1"){
      a=substr(line,7,9)+0
      b=substr(line,16,9)+0
      c=substr(line,25,9)+0
      box=sprintf("%10.5f%10.5f%10.5f", a/10.0, b/10.0, c/10.0)
      break
    }
  }
  close(pdb)

  if(box==""){
    print "ERROR: no CRYST1 found in membrane.pdb" > "/dev/stderr"
    exit 1
  }

  # read protein atoms
  getline title < protein
  getline n < protein
  n += 0
  for(i=1;i<=n;i++){
    getline atom < protein
    atoms[++total]=atom
  }
  close(protein)

  # read non-protein atoms
  getline title < nonprot
  getline n < nonprot
  n += 0
  for(i=1;i<=n;i++){
    getline atom < nonprot
    atoms[++total]=atom
  }
  close(nonprot)

  # write merged gro
  print "Mutated protein + AMBER lipid17 membrane system" > out
  printf "%5d\n", total >> out

  for(i=1;i<=total;i++){
    # renumber atom serial columns 16-20
    printf "%s%5d%s\n", substr(atoms[i],1,15), i%100000, substr(atoms[i],21) >> out
  }

  print box >> out

  print "Total atoms: " total
  print "Box: " box
}
'

# 3.6 Add water/ion/lipid into topol.top
awk -v gro="membrane_amber_lipid.gro" '
BEGIN{
  while((getline line < gro) > 0){
    nline++
    lines[nline]=line
  }

  natoms=lines[2]+0

  for(i=3; i<=2+natoms; i++){
    resname=substr(lines[i],6,5)
    gsub(/ /,"",resname)
    atom_count[resname]++
  }

  pa_atoms=atom_count["PA"]+0
  pc_atoms=atom_count["PC"]+0
  ol_atoms=atom_count["OL"]+0
  sol_atoms=atom_count["SOL"]+0
  k_atoms=atom_count["K"]+0
  cl_atoms=atom_count["CL"]+0

  lipid_from_pa=pa_atoms/46
  lipid_from_pc=pc_atoms/38
  lipid_from_ol=ol_atoms/50

  if(lipid_from_pa != int(lipid_from_pa) ||
     lipid_from_pc != int(lipid_from_pc) ||
     lipid_from_ol != int(lipid_from_ol) ||
     lipid_from_pa != lipid_from_pc ||
     lipid_from_pa != lipid_from_ol){
    print "ERROR: inconsistent lipid atom counts." > "/dev/stderr"
    printf "PA atoms=%d -> %.3f lipids\n", pa_atoms, lipid_from_pa > "/dev/stderr"
    printf "PC atoms=%d -> %.3f lipids\n", pc_atoms, lipid_from_pc > "/dev/stderr"
    printf "OL atoms=%d -> %.3f lipids\n", ol_atoms, lipid_from_ol > "/dev/stderr"
    exit 1
  }

  if(sol_atoms % 3 != 0){
    print "ERROR: SOL atom count is not divisible by 3." > "/dev/stderr"
    printf "SOL atoms=%d\n", sol_atoms > "/dev/stderr"
    exit 1
  }

  lipid=lipid_from_pa
  sol=sol_atoms/3
  kion=k_atoms
  clion=cl_atoms

  print "Will insert into [ molecules ]:" > "/dev/stderr"
  printf "%-18s %d\n", "PA_PC_OL", lipid > "/dev/stderr"
  printf "%-18s %d\n", "K", kion > "/dev/stderr"
  printf "%-18s %d\n", "CL", clion > "/dev/stderr"
  printf "%-18s %d\n", "SOL", sol > "/dev/stderr"
}

/^\[ *molecules *\]/{
  in_mol=1
  print
  next
}

{
  if(in_mol && $1 ~ /^Protein_chain_/ && $2 == 1){
    pending_insert=1
  }
  else if(in_mol && pending_insert){
    printf "%-18s %d\n", "PA_PC_OL", lipid
    printf "%-18s %d\n", "K", kion
    printf "%-18s %d\n", "CL", clion
    printf "%-18s %d\n", "SOL", sol
    pending_insert=0
    in_mol=0
  }

  print
}

END{
  if(in_mol && pending_insert){
    printf "%-18s %d\n", "PA_PC_OL", lipid
    printf "%-18s %d\n", "K", kion
    printf "%-18s %d\n", "CL", clion
    printf "%-18s %d\n", "SOL", sol
  }
}
' topol.top > topol_tmp.top && mv topol_tmp.top topol.top

# 3.7 Add lipid17.itp into topol.top
awk '/forcefield.itp/{print; print "#include \"lipid17.itp\""; next} {print}' topol.top > tmp.top && mv tmp.top topol.top

echo "=================================================="
echo "MD Preparation: Done!"
echo "=================================================="

# 4 MD simulation
# sbatch MD_Simulation.sh

echo "=================================================="
echo "MD Simulation: Start!"
echo "=================================================="