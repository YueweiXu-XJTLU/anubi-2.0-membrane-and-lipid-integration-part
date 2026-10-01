# Membrane

## 1. Fix incompatible residue names and atom names (protein + water + ion)

1. Input: `membrane.pdb` 

2. Protein:

   - ARG HB23/HG23/HD23 → HB12/HG12/HD12
   
   
      - ASP HB23 → HB12
   
   
      - CLYS HB/HG/HD/HE 23 → 12, O/OXT → OC1/OC2
   
   
      - CYS HB23 → HB12
   
   
      - GLU HB23/HG23 → HB12/HG12
   
   
      - HID HB23 → HB12
   
   
      - ILE CD1 → CD, HD111213 → HD123, HG1213 → HG1112
   
   
      - LEU HB23 → HB12
   
   
      - MET HB23/HG23 → HB12/HG12
   
   
      - NASP HB23 → HB12
   
   
      - NGLU HB23/HG23 → HB12/HG12
   
   
      - NMET HB23/HG23 → HB12/HG12
   
   
      - PHE HB23 → HB12
   
   
      - PRO HB/HG/HD 23 → 12
   
   
      - SER HB23 → HB12
   
   
      - TRP HB23 → HB12
   
   
      - TYR HB23 → HB12
   
3. Water + Ion:

   - TP3  → SO
     - O   → OW
     - H1  → HW1
     - H2  → HW2
   - K+  → K
   - Cl- → CL
   
   

```{bash}
# 1. backup
cp -p membrane.pdb membrane.pdb.bak

# 2. protein
awk 'BEGIN{FS=OFS=""} /^(ATOM|HETATM)/{new=""; res=substr($0,18,4); gsub(/ /,"",res); atom=substr($0,13,4); a=atom; gsub(/ /,"",a); if(res=="ARG"){if(a=="HB2") new="HB1"; else if(a=="HB3") new="HB2"; else if(a=="HG2") new="HG1"; else if(a=="HG3") new="HG2"; else if(a=="HD2") new="HD1"; else if(a=="HD3") new="HD2"; if(new!="") $0=substr($0,1,12) sprintf("%4s",new) substr($0,17)}} {print}' membrane.pdb > membrane_tmp.pdb && mv membrane_tmp.pdb membrane.pdb

awk 'BEGIN{FS=OFS=""} /^(ATOM|HETATM)/{new=""; res=substr($0,18,4); gsub(/ /,"",res); atom=substr($0,13,4); a=atom; gsub(/ /,"",a); if(res=="ASN"){if(a=="HB2") new="HB1"; else if(a=="HB3") new="HB2"; if(new!="") $0=substr($0,1,12) sprintf("%4s",new) substr($0,17)}} {print}' membrane.pdb > membrane_tmp.pdb && mv membrane_tmp.pdb membrane.pdb

awk 'BEGIN{FS=OFS=""} /^(ATOM|HETATM)/{new=""; res=substr($0,18,4); gsub(/ /,"",res); atom=substr($0,13,4); a=atom; gsub(/ /,"",a); if(res=="ASP"){if(a=="HB2") new="HB1"; else if(a=="HB3") new="HB2"; if(new!="") $0=substr($0,1,12) sprintf("%4s",new) substr($0,17)}} {print}' membrane.pdb > membrane_tmp.pdb && mv membrane_tmp.pdb membrane.pdb

awk 'BEGIN{FS=OFS=""} /^(ATOM|HETATM)/{new=""; res=substr($0,18,4); gsub(/ /,"",res); atom=substr($0,13,4); a=atom; gsub(/ /,"",a); if(res=="CALA"){if(a=="O") new="OC1"; else if(a=="OXT") new="OC2"; if(new!="") $0=substr($0,1,12) sprintf("%4s",new) substr($0,17)}} {print}' membrane.pdb > membrane_tmp.pdb && mv membrane_tmp.pdb membrane.pdb

awk 'BEGIN{FS=OFS=""} /^(ATOM|HETATM)/{new=""; res=substr($0,18,4); gsub(/ /,"",res); atom=substr($0,13,4); a=atom; gsub(/ /,"",a); if(res=="CLYS"){if(a=="HB2") new="HB1"; else if(a=="HB3") new="HB2"; else if(a=="HG2") new="HG1"; else if(a=="HG3") new="HG2"; else if(a=="HD2") new="HD1"; else if(a=="HD3") new="HD2"; else if(a=="HE2") new="HE1"; else if(a=="HE3") new="HE2"; else if(a=="O") new="OC1"; else if(a=="OXT") new="OC2"; if(new!="") $0=substr($0,1,12) sprintf("%4s",new) substr($0,17)}} {print}' membrane.pdb > membrane_tmp.pdb && mv membrane_tmp.pdb membrane.pdb

awk 'BEGIN{FS=OFS=""} /^(ATOM|HETATM)/{new=""; res=substr($0,18,4); gsub(/ /,"",res); atom=substr($0,13,4); a=atom; gsub(/ /,"",a); if(res=="CSER"){if(a=="HB2") new="HB1"; else if(a=="HB3") new="HB2"; else if(a=="O") new="OC1"; else if(a=="OXT") new="OC2"; if(new!="") $0=substr($0,1,12) sprintf("%4s",new) substr($0,17)}} {print}' membrane.pdb > membrane_tmp.pdb && mv membrane_tmp.pdb membrane.pdb

awk 'BEGIN{FS=OFS=""} /^(ATOM|HETATM)/{new=""; res=substr($0,18,4); gsub(/ /,"",res); atom=substr($0,13,4); a=atom; gsub(/ /,"",a); if(res=="CYS"){if(a=="HB2") new="HB1"; else if(a=="HB3") new="HB2"; if(new!="") $0=substr($0,1,12) sprintf("%4s",new) substr($0,17)}} {print}' membrane.pdb > membrane_tmp.pdb && mv membrane_tmp.pdb membrane.pdb

awk 'BEGIN{FS=OFS=""} /^(ATOM|HETATM)/{new=""; res=substr($0,18,4); gsub(/ /,"",res); atom=substr($0,13,4); a=atom; gsub(/ /,"",a); if(res=="GLN"){if(a=="HB2") new="HB1"; else if(a=="HB3") new="HB2"; else if(a=="HG2") new="HG1"; else if(a=="HG3") new="HG2"; if(new!="") $0=substr($0,1,12) sprintf("%4s",new) substr($0,17)}} {print}' membrane.pdb > membrane_tmp.pdb && mv membrane_tmp.pdb membrane.pdb

awk 'BEGIN{FS=OFS=""} /^(ATOM|HETATM)/{new=""; res=substr($0,18,4); gsub(/ /,"",res); atom=substr($0,13,4); a=atom; gsub(/ /,"",a); if(res=="GLU"){if(a=="HB2") new="HB1"; else if(a=="HB3") new="HB2"; else if(a=="HG2") new="HG1"; else if(a=="HG3") new="HG2"; if(new!="") $0=substr($0,1,12) sprintf("%4s",new) substr($0,17)}} {print}' membrane.pdb > membrane_tmp.pdb && mv membrane_tmp.pdb membrane.pdb

awk 'BEGIN{FS=OFS=""} /^(ATOM|HETATM)/{new=""; res=substr($0,18,4); gsub(/ /,"",res); atom=substr($0,13,4); a=atom; gsub(/ /,"",a); if(res=="GLY"){if(a=="HA2") new="HA1"; else if(a=="HA3") new="HA2"; if(new!="") $0=substr($0,1,12) sprintf("%4s",new) substr($0,17)}} {print}' membrane.pdb > membrane_tmp.pdb && mv membrane_tmp.pdb membrane.pdb

awk 'BEGIN{FS=OFS=""} /^(ATOM|HETATM)/{new=""; res=substr($0,18,4); gsub(/ /,"",res); atom=substr($0,13,4); a=atom; gsub(/ /,"",a); if(res=="HID"){if(a=="HB2") new="HB1"; else if(a=="HB3") new="HB2"; if(new!="") $0=substr($0,1,12) sprintf("%4s",new) substr($0,17)}} {print}' membrane.pdb > membrane_tmp.pdb && mv membrane_tmp.pdb membrane.pdb

awk 'BEGIN{FS=OFS=""} /^(ATOM|HETATM)/{new=""; res=substr($0,18,4); gsub(/ /,"",res); atom=substr($0,13,4); a=atom; gsub(/ /,"",a); if(res=="ILE"){if(a=="CD1") new="CD"; else if(a=="HD11") new="HD1"; else if(a=="HD12") new="HD2"; else if(a=="HD13") new="HD3"; else if(a=="HG12") new="HG11"; else if(a=="HG13") new="HG12"; if(new!="") $0=substr($0,1,12) sprintf("%4s",new) substr($0,17)}} {print}' membrane.pdb > membrane_tmp.pdb && mv membrane_tmp.pdb membrane.pdb

awk 'BEGIN{FS=OFS=""} /^(ATOM|HETATM)/{new=""; res=substr($0,18,4); gsub(/ /,"",res); atom=substr($0,13,4); a=atom; gsub(/ /,"",a); if(res=="LEU"){if(a=="HB2") new="HB1"; else if(a=="HB3") new="HB2"; if(new!="") $0=substr($0,1,12) sprintf("%4s",new) substr($0,17)}} {print}' membrane.pdb > membrane_tmp.pdb && mv membrane_tmp.pdb membrane.pdb

awk 'BEGIN{FS=OFS=""} /^(ATOM|HETATM)/{new=""; res=substr($0,18,4); gsub(/ /,"",res); atom=substr($0,13,4); a=atom; gsub(/ /,"",a); if(res=="LYS"){if(a=="HB2") new="HB1"; else if(a=="HB3") new="HB2"; else if(a=="HG2") new="HG1"; else if(a=="HG3") new="HG2"; else if(a=="HD2") new="HD1"; else if(a=="HD3") new="HD2"; else if(a=="HE2") new="HE1"; else if(a=="HE3") new="HE2"; if(new!="") $0=substr($0,1,12) sprintf("%4s",new) substr($0,17)}} {print}' membrane.pdb > membrane_tmp.pdb && mv membrane_tmp.pdb membrane.pdb

awk 'BEGIN{FS=OFS=""} /^(ATOM|HETATM)/{new=""; res=substr($0,18,4); gsub(/ /,"",res); atom=substr($0,13,4); a=atom; gsub(/ /,"",a); if(res=="MET"){if(a=="HB2") new="HB1"; else if(a=="HB3") new="HB2"; else if(a=="HG2") new="HG1"; else if(a=="HG3") new="HG2"; if(new!="") $0=substr($0,1,12) sprintf("%4s",new) substr($0,17)}} {print}' membrane.pdb > membrane_tmp.pdb && mv membrane_tmp.pdb membrane.pdb

awk 'BEGIN{FS=OFS=""} /^(ATOM|HETATM)/{new=""; res=substr($0,18,4); gsub(/ /,"",res); atom=substr($0,13,4); a=atom; gsub(/ /,"",a); if(res=="NASP"){if(a=="HB2") new="HB1"; else if(a=="HB3") new="HB2"; if(new!="") $0=substr($0,1,12) sprintf("%4s",new) substr($0,17)}} {print}' membrane.pdb > membrane_tmp.pdb && mv membrane_tmp.pdb membrane.pdb

awk 'BEGIN{FS=OFS=""} /^(ATOM|HETATM)/{new=""; res=substr($0,18,4); gsub(/ /,"",res); atom=substr($0,13,4); a=atom; gsub(/ /,"",a); if(res=="NGLU"){if(a=="HB2") new="HB1"; else if(a=="HB3") new="HB2"; else if(a=="HG2") new="HG1"; else if(a=="HG3") new="HG2"; if(new!="") $0=substr($0,1,12) sprintf("%4s",new) substr($0,17)}} {print}' membrane.pdb > membrane_tmp.pdb && mv membrane_tmp.pdb membrane.pdb

awk 'BEGIN{FS=OFS=""} /^(ATOM|HETATM)/{new=""; res=substr($0,18,4); gsub(/ /,"",res); atom=substr($0,13,4); a=atom; gsub(/ /,"",a); if(res=="NMET"){if(a=="HB2") new="HB1"; else if(a=="HB3") new="HB2"; else if(a=="HG2") new="HG1"; else if(a=="HG3") new="HG2"; if(new!="") $0=substr($0,1,12) sprintf("%4s",new) substr($0,17)}} {print}' membrane.pdb > membrane_tmp.pdb && mv membrane_tmp.pdb membrane.pdb

awk 'BEGIN{FS=OFS=""} /^(ATOM|HETATM)/{new=""; res=substr($0,18,4); gsub(/ /,"",res); atom=substr($0,13,4); a=atom; gsub(/ /,"",a); if(res=="PHE"){if(a=="HB2") new="HB1"; else if(a=="HB3") new="HB2"; if(new!="") $0=substr($0,1,12) sprintf("%4s",new) substr($0,17)}} {print}' membrane.pdb > membrane_tmp.pdb && mv membrane_tmp.pdb membrane.pdb

awk 'BEGIN{FS=OFS=""} /^(ATOM|HETATM)/{new=""; res=substr($0,18,4); gsub(/ /,"",res); atom=substr($0,13,4); a=atom; gsub(/ /,"",a); if(res=="PRO"){if(a=="HB2") new="HB1"; else if(a=="HB3") new="HB2"; else if(a=="HG2") new="HG1"; else if(a=="HG3") new="HG2"; else if(a=="HD2") new="HD1"; else if(a=="HD3") new="HD2"; if(new!="") $0=substr($0,1,12) sprintf("%4s",new) substr($0,17)}} {print}' membrane.pdb > membrane_tmp.pdb && mv membrane_tmp.pdb membrane.pdb

awk 'BEGIN{FS=OFS=""} /^(ATOM|HETATM)/{new=""; res=substr($0,18,4); gsub(/ /,"",res); atom=substr($0,13,4); a=atom; gsub(/ /,"",a); if(res=="SER"){if(a=="HB2") new="HB1"; else if(a=="HB3") new="HB2"; if(new!="") $0=substr($0,1,12) sprintf("%4s",new) substr($0,17)}} {print}' membrane.pdb > membrane_tmp.pdb && mv membrane_tmp.pdb membrane.pdb

awk 'BEGIN{FS=OFS=""} /^(ATOM|HETATM)/{new=""; res=substr($0,18,4); gsub(/ /,"",res); atom=substr($0,13,4); a=atom; gsub(/ /,"",a); if(res=="TRP"){if(a=="HB2") new="HB1"; else if(a=="HB3") new="HB2"; if(new!="") $0=substr($0,1,12) sprintf("%4s",new) substr($0,17)}} {print}' membrane.pdb > membrane_tmp.pdb && mv membrane_tmp.pdb membrane.pdb

awk 'BEGIN{FS=OFS=""} /^(ATOM|HETATM)/{new=""; res=substr($0,18,4); gsub(/ /,"",res); atom=substr($0,13,4); a=atom; gsub(/ /,"",a); if(res=="TYR"){if(a=="HB2") new="HB1"; else if(a=="HB3") new="HB2"; if(new!="") $0=substr($0,1,12) sprintf("%4s",new) substr($0,17)}} {print}' membrane.pdb > membrane_tmp.pdb && mv membrane_tmp.pdb membrane.pdb

# 3. water and ions
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
{print}' membrane.pdb > membrane_tmp.pdb && mv membrane_tmp.pdb membrane.pdb
```

## 2. Charmm Lipid → Amber Lipid

- Use files from [here](https://github.com/xiki-tempula/gmx_lipid17.ff/tree/master), make sure the work directory contains: `membrane.pdb` & `charmmlipid2amber.py` & `gro/PA_PC_OL.gro`

- Update MDAnalysis: `pip install --upgrade MDAnalysis -i https://pypi.tuna.tsinghua.edu.cn/simple --trusted-host pypi.tuna.tsinghua.edu.cn`

- Convert: 

  ```{bash}
  python charmmlipid2amber.py -f membrane.pdb -o membrane_amber_lipid.gro -lipids POPC
  ```

- Output: `PA_PC_OL 417` means 417 POPC converted into PA + PC + OL in `membrane_amber_lipid.gro`

## 3. topol preparation

1. Make a topology file with only protein under orginal force field

   1. Extract solely protein from `membrane_amber_lipid.gro`

      ```{bash}
      awk '
      NR==1{title=$0; next}
      NR==2{next}
      {line[NR]=$0}
      END{
        print title
        n=0
        for(i=3;i<NR;i++){
          res=substr(line[i],6,5)
          gsub(/ /,"",res)
          if(res!="PA" && res!="PC" && res!="OL" && res!="POPC" &&
             res!="SOL" &&
             res!="K" &&
             res!="CL"){
            n++
            keep[n]=line[i]
          }
        }
        printf "%5d\n", n
        for(i=1;i<=n;i++){
          printf "%s%5d%s\n", substr(keep[i],1,15), i, substr(keep[i],21)
        }
        print line[NR]
      }' membrane_amber_lipid.gro > protein_only.gro
      ```

   2. Generate `topol.top` with only protein:

      **make sure to using the *modified* `amber14sb_OL15`**

      ```{bash}
      gmx_mpi pdb2gmx -f protein_only.gro -ff amber14sb_OL15 -water tip3p
      ```

1. Add water & ions into `topol.top`

   1. Calculate total number of water/K/CL

      ```{bash}
      awk -v last=$(wc -l < membrane_amber_lipid.gro) '
      NR>2 && NR<last{
        resid=substr($0,1,5)
        res=substr($0,6,5)
        gsub(/ /,"",resid)
        gsub(/ /,"",res)
        key=resid":"res
        seen[key]=res
      }
      END{
        for(k in seen) count[seen[k]]++
        print "K",   count["K"]+0
        print "CL",  count["CL"]+0
        print "SOL", count["SOL"]+0
      }' membrane_amber_lipid.gro
      ```

      output: 

      ​	`K 207`
      ​	`CL 251`
      ​	`SOL 75810`
      
   1. Write into `topol.top`:

      ```{bash}
      cat >> topol.top <<'EOF'
      K        207
      CL       251
      SOL      75810
      EOF      
      ```
      

1. Add Lipid into `topol.top`

   1. Include `lipid17.ff`&`lipid17.itp` **(make sure they are in work directory)**

      ```{bash}
      awk '/forcefield.itp/{print; print "#include \"lipid17.itp\""; next} {print}' topol.top > tmp.top && mv tmp.top topol.top
      ```

      if use `amber14sb`, consider use this:

      ```{bash}
      awk '/forcefield.itp/{print; print "#include \"lipid17.ff/ffnonbonded.itp\""; print "#include \"lipid17.ff/ffbonded.itp\""; print "#include \"lipid17.itp\""; next} {print}' topol.top > topol_tmp.top && mv topol_tmp.top topol.top
      ```
   
   2. Add `PA_PC_OL      417` into `[ molecules ]`
   
      ```{bash}
      awk '{print} $1=="Protein" && $2=="1"{print "PA_PC_OL      417"}' topol.top > topol_tmp.top && mv topol_tmp.top topol.top
      ```

## 4. EM + NVT + NPT ***(NET CHARGE = -0.000851, causing an error)***

- Better to use `constraints		=	h-bonds` in `NVT.mdp` & `NPT.mdp` 

- Remove definition for RNA in `md.mdp`

  > tc-grps                  = System
  > ; Time constant (ps) and reference temperature (K)
  > tau_t                    = 1.0
  > ref_t                    = 300 

- Better to use `nstcomm                  = 1000` in `md.mdp`

- Process:

  ```{bash}
  #!/bin/bash
  #SBATCH -J  membrane
  #SBATCH --partition=gpu4090         # Partition name
  #SBATCH -n 1                         # Number of tasks (typically 1 if not using MPI)
  #SBATCH --qos=1gpu               # Quality of Service (QOS) setting
  #SBATCH --cpus-per-task=4       # 1:4 GPU:CPU ratio
  #SBATCH --gres=gpu:1
  #SBATCH -e  job.%j.err
  #SBATCH -o  job.%j.out
  #SBATCH --mail-user=Yuewei.Xu23@student.xjtlu.edu.cn
  #SBATCH --mail-type=ALL
  
  set -euo pipefail
  
  ml load gromacs/2023.2-gcc-9.5.0-jzxesel
  
  # EM
  gmx_mpi grompp -f minim.mdp -c membrane_amber_lipid.gro -p topol.top -o minim.tpr -maxwarn 1
  gmx_mpi mdrun -deffnm minim -v
  
  # NVT
  gmx_mpi grompp -f NVT.mdp -c minim.gro -r minim.gro -p topol.top -o nvt.tpr -maxwarn 1
  gmx_mpi mdrun -deffnm nvt -v
  
  # NPT
  gmx_mpi grompp -f NPT.mdp -c nvt.gro -r nvt.gro -p topol.top -o npt.tpr -maxwarn 1
  gmx_mpi mdrun -deffnm npt -v
  
  # Production
  gmx_mpi grompp -f md.mdp -c npt.gro -r npt.gro -p topol.top -o md.tpr -maxwarn 1
  gmx_mpi mdrun -deffnm md -v
  ```
  
  

