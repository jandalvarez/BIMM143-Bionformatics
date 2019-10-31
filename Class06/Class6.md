Class6 R Funxns
================

# **This is H1**

This is my work for **BIMM143**.

## This is H2

### This is H3

``` r
# This is to show how to use R in an R markdown file; essentially, this is R within a different "program"
plot(1:10)
```

![](Class6_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

## Practice reading files (again…)

read.delim(“filename.txt”) **TAB**

read.csv(“filename.txt”) **COMMA**

read.csv2(“filename.txt”) **SEMICOLON**

read.table(“filename.txt”) **SPACE**

``` r
# This was the hard way bc read.csv2 wouldve given the same output

read.table("test1.txt", header = TRUE, sep = ",")
```

    ##   Col1 Col2 Col3
    ## 1    1    2    3
    ## 2    4    5    6
    ## 3    7    8    9
    ## 4    a    b    c

``` r
# Same as last one but read.csv2

read.csv("test1.txt")
```

    ##   Col1 Col2 Col3
    ## 1    1    2    3
    ## 2    4    5    6
    ## 3    7    8    9
    ## 4    a    b    c

``` r
# Import test2; sep by $ so put sep="$"; header to look nice

read.table("test2.txt", header = TRUE, sep = "$")
```

    ##   Col1 Col2 Col3
    ## 1    1    2    3
    ## 2    4    5    6
    ## 3    7    8    9
    ## 4    a    b    c

``` r
#Import test3; 
read.table("test3.txt")
```

    ##   V1 V2 V3
    ## 1  1  6  a
    ## 2  2  7  b
    ## 3  3  8  c
    ## 4  4  9  d
    ## 5  5 10  e

# **Writing Functions**

format:

name.of.function \<- **function**(arg1, arg2) {statements
**return**(something)}

## Writing add funxn

``` r
add <- function(x, y=1)
  { # Sum the input x and y
 x + y}
```

``` r
add(1)
```

    ## [1] 2

``` r
# This will give 10 bc we're overriding the default value of y=1 to make it y=5

add(5,5)
```

    ## [1] 10

``` r
# These are only for the x values, this allows u to write multiple x values at once; same as writing c(1:3); or whatever the number range is

add( c(1, 2, 3))
```

    ## [1] 2 3 4

``` r
add(c(1:3))
```

    ## [1] 2 3 4

``` r
# X values are 1,2,3; y value overides default and y = 4
add( c(1, 2, 3), 4)
```

    ## [1] 5 6 7

## Rescale Funxn

``` r
# Rescale Funxn

rescale <- function(x) {
 rng <-range(x)
 (x - rng[1]) / (rng[2] - rng[1])
}



# Testing it out

rescale(1:10)
```

    ##  [1] 0.0000000 0.1111111 0.2222222 0.3333333 0.4444444 0.5555556 0.6666667
    ##  [8] 0.7777778 0.8888889 1.0000000

``` r
# NA screws everything over

rescale(c(1,2,NA,3,10))
```

    ## [1] NA NA NA NA NA

``` r
# Gives error code bc non number

# rescale(c(1,10,"string"))
```

``` r
#Remove NA's by setting na.rm = TRUE 

x <- c(1,2,NA,3,10)
rng <-range(x,na.rm = TRUE)
rng
```

    ## [1]  1 10

``` r
# Change original rescale funxn to be able to omit NA's by adding na.rm = TRUE new funxn rescale2

rescale2 <- function(x) {
 rng <-range(x, na.rm = TRUE)
 (x - rng[1]) / (rng[2] - rng[1])
}

# Test it out

rescale2(c(1,2,NA,3,10))
```

    ## [1] 0.0000000 0.1111111        NA 0.2222222 1.0000000

``` r
# Prof Grant's code for rescale

BestRescale <- function(x, na.rm=TRUE, plot=FALSE) {
 rng <-range(x, na.rm=na.rm)
 print("Hello")
 answer <- (x - rng[1]) / (rng[2] - rng[1])
 return(answer)
 print("is it me you are looking for?")
 if(plot) {
 plot(answer, typ="b", lwd=4)
 }
 print("I can see it in ...")
 #return(answer)
}

# Test it 

BestRescale(1:10)
```

    ## [1] "Hello"

    ##  [1] 0.0000000 0.1111111 0.2222222 0.3333333 0.4444444 0.5555556 0.6666667
    ##  [8] 0.7777778 0.8888889 1.0000000

``` r
# Add a plot

BestRescale(1:10, plot = TRUE)
```

    ## [1] "Hello"

    ##  [1] 0.0000000 0.1111111 0.2222222 0.3333333 0.4444444 0.5555556 0.6666667
    ##  [8] 0.7777778 0.8888889 1.0000000

``` r
# Changed the code to stop after returning the answers; anything after return "technically" no longer exists

BestRescale2 <- function(x, na.rm=TRUE, plot=FALSE) {
 rng <-range(x, na.rm=na.rm)
 print("Hello")
 answer <- (x - rng[1]) / (rng[2] - rng[1])
 return(answer)
 print("is it me you are looking for?")
 if(plot) {
 plot(answer, typ="b", lwd=4)
 }
 print("I can see it in ...")
 #return(answer)
}

#  Test it out

BestRescale2(1:10)
```

    ## [1] "Hello"

    ##  [1] 0.0000000 0.1111111 0.2222222 0.3333333 0.4444444 0.5555556 0.6666667
    ##  [8] 0.7777778 0.8888889 1.0000000

``` r
# do the thing, don't want to install in RMD, install in console

# install.packages("bio3d")

# Running Chonky code

library(bio3d)
s1 <- read.pdb("4AKE") # kinase with drug
```

    ##   Note: Accessing on-line PDB file

``` r
s2 <- read.pdb("1AKE") # kinase no drug
```

    ##   Note: Accessing on-line PDB file
    ##    PDB has ALT records, taking A only, rm.alt=TRUE

``` r
s3 <- read.pdb("1E4Y") # kinase with drug
```

    ##   Note: Accessing on-line PDB file

``` r
s1.chainA <- trim.pdb(s1, chain="A", elety="CA")
s2.chainA <- trim.pdb(s2, chain="A", elety="CA")
s3.chainA <- trim.pdb(s3, chain="A", elety="CA")
s1.b <- s1.chainA$atom$b
s2.b <- s2.chainA$atom$b
s3.b <- s3.chainA$atom$b
plotb3(s1.b, sse=s1.chainA, typ="l", ylab="Bfactor")
```

![](Class6_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

``` r
plotb3(s2.b, sse=s2.chainA, typ="l", ylab="Bfactor")
```

![](Class6_files/figure-gfm/unnamed-chunk-10-2.png)<!-- -->

``` r
plotb3(s3.b, sse=s3.chainA, typ="l", ylab="Bfactor")
```

![](Class6_files/figure-gfm/unnamed-chunk-10-3.png)<!-- -->

``` r
trim.pdb(s1, chain="A",elety="CA")
```

    ## 
    ##  Call:  trim.pdb(pdb = s1, chain = "A", elety = "CA")
    ## 
    ##    Total Models#: 1
    ##      Total Atoms#: 214,  XYZs#: 642  Chains#: 1  (values: A)
    ## 
    ##      Protein Atoms#: 214  (residues/Calpha atoms#: 214)
    ##      Nucleic acid Atoms#: 0  (residues/phosphate atoms#: 0)
    ## 
    ##      Non-protein/nucleic Atoms#: 0  (residues: 0)
    ##      Non-protein/nucleic resid values: [ none ]
    ## 
    ##    Protein sequence:
    ##       MRIILLGAPGAGKGTQAQFIMEKYGIPQISTGDMLRAAVKSGSELGKQAKDIMDAGKLVT
    ##       DELVIALVKERIAQEDCRNGFLLDGFPRTIPQADAMKEAGINVDYVLEFDVPDELIVDRI
    ##       VGRRVHAPSGRVYHVKFNPPKVEGKDDVTGEELTTRKDDQEETVRKRLVEYHQMTAPLIG
    ##       YYSKEAEAGNTKYAKVDGTKPVAEVRADLEKILG
    ## 
    ## + attr: atom, helix, sheet, seqres, xyz,
    ##         calpha, call

``` r
# call library so all packages are available to u

library(bio3d)

# Make funxn to trim in NEW FILE
```
