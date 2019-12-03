---
title: "Class 11: Structural Bioinformatics 1"
output:
  html_document:
    keep_md: yes
---




## The PDB Database for Biomolecular Structure Data

>Q1) Determine % of structures solved by X-Ray and Electron Microscopy.



```r
# Read csv file

PDBstats <- read.csv("Data Export Summary.csv")

PDBstats
```

```
##   Experimental.Method Proteins Nucleic.Acids Protein.NA.Complex Other
## 1               X-Ray   131278          2059               6759     8
## 2                 NMR    11235          1303                261     8
## 3 Electron Microscopy     2899            32                999     0
## 4               Other      280             4                  6    13
## 5        Multi Method      144             5                  2     1
##    Total
## 1 140104
## 2  12807
## 3   3930
## 4    303
## 5    152
```

```r
# Determining % of structures for X-Ray and Electron Microscopy

# total <- sum(PDBstats$Total)

# x_per <- prop.table(PDBstats)
```



>What proportion are from each method?


```r
(PDBstats$Total/sum(PDBstats$Total))*100
```

```
## [1] 89.0702879  8.1419744  2.4984742  0.1926305  0.0966331
```

Proportion that are protein


```r
# (PDBstats$Proteins/(sum(PDBstats$Total))*100

# use round function to change sig figs

round((sum(PDBstats$Proteins)/(sum(PDBstats$Total)))*100, 2)
```

```
## [1] 92.71
```



## Atom Selection



```r
# Read pdb file
library(bio3d)

pdb <- read.pdb("1hsg.pdb")


a.inds <- atom.select(pdb, chain = "A")

# Select C-alphas of chain A

ca.inds <- atom.select(pdb, "calpha", chain="A")

# We can combine multiple selection criteria to return their
# intersection

cab.inds <- atom.select(pdb, elety=c("CA","CB"), chain="A",
resno=10:20)



# only look at proteins

atom.select(pdb, "protein", value = TRUE)
```

```
## 
##  Call:  trim.pdb(pdb = pdb, sele)
## 
##    Total Models#: 1
##      Total Atoms#: 1514,  XYZs#: 4542  Chains#: 2  (values: A B)
## 
##      Protein Atoms#: 1514  (residues/Calpha atoms#: 198)
##      Nucleic acid Atoms#: 0  (residues/phosphate atoms#: 0)
## 
##      Non-protein/nucleic Atoms#: 0  (residues: 0)
##      Non-protein/nucleic resid values: [ none ]
## 
##    Protein sequence:
##       PQITLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYD
##       QILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNFPQITLWQRPLVTIKIGGQLKE
##       ALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTP
##       VNIIGRNLLTQIGCTLNF
## 
## + attr: atom, helix, sheet, seqres, xyz,
##         calpha, call
```

```r
protein <- atom.select(pdb, "protein", value = TRUE)
# write.pdb(ligand, file = "1hsg_protein.pdb")
# ligands

ligand <- atom.select(pdb, "ligand", value = TRUE)
write.pdb(ligand, file = "1hsg_ligand.pdb")
```









