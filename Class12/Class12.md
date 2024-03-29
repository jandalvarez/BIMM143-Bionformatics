---
title: "Class 12: Structural Bioninformatics 2"
output:
  html_document:
    keep_md: yes
---



## Prep for Docking

We want to produce a pdb file that's protein only and a drug only pdb file


```r
# Load Bio3d; 

library(bio3d)

file.name <- get.pdb("1hsg")
```

```
## Warning in get.pdb("1hsg"): ./1hsg.pdb exists. Skipping download
```

```r
# 

pdb <- read.pdb(file.name)

pdb
```

```
## 
##  Call:  read.pdb(file = file.name)
## 
##    Total Models#: 1
##      Total Atoms#: 1686,  XYZs#: 5058  Chains#: 2  (values: A B)
## 
##      Protein Atoms#: 1514  (residues/Calpha atoms#: 198)
##      Nucleic acid Atoms#: 0  (residues/phosphate atoms#: 0)
## 
##      Non-protein/nucleic Atoms#: 172  (residues: 128)
##      Non-protein/nucleic resid values: [ HOH (127), MK1 (1) ]
## 
##    Protein sequence:
##       PQITLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYD
##       QILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNFPQITLWQRPLVTIKIGGQLKE
##       ALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTP
##       VNIIGRNLLTQIGCTLNF
## 
## + attr: atom, xyz, seqres, helix, sheet,
##         calpha, remark, call
```

>Q1) What is the name of the two non protein resid values in this structure? What does resid
correspond to and how would you get a listing of all reside values in this structure? 

- HOH, Water; MK1
- 


```r
# Make protein only pdb file

 pdb <- read.pdb("1hsg.pdb")
prot <- atom.select(pdb, "protein", value = TRUE)
write.pdb(prot, file = "1hsg_protein.pdb")
```



```r
# Make ligand(drug) only pdb file

 pdb <- read.pdb("1hsg.pdb")
lig <- atom.select(pdb, "ligand", value = TRUE)
write.pdb(lig, file = "1hsg_ligand.pdb")
```

*******************************************************
*******************************************************

## Analyzing Docking Results


```r
# 

library(bio3d)
res <- read.pdb("all.pdbqt", multi=TRUE)
write.pdb(res, "results.pdb")
```

















