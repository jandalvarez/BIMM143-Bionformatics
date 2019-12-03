---
title: "Class 7 R functions/practice"
author: "Justin"
date: "10/22/2019"
output:
  html_document:
    keep_md: yes
---



**#R funxns Revisited**

Source functions from last class, in order to use them for this class.


```r
# source() runs another file, allows u to pull something from another locxn or sumn; Pushing play button sources the HTML below and adds them to Global Environment

source("http://tinyurl.com/rescale-R")
```



```r
# Test funxn out

rescale(1:10)
```

```
##  [1] 0.0000000 0.1111111 0.2222222 0.3333333 0.4444444 0.5555556 0.6666667
##  [8] 0.7777778 0.8888889 1.0000000
```

```r
# Test more; works with NA bc rescale default is na.rm = TRUE

rescale(c(1,10,5,NA,6))
```

```
## [1] 0.0000000 1.0000000 0.4444444        NA 0.5555556
```


```r
# Working with is.numeric and ! to flip logical

is.numeric(c(5:10))
```

```
## [1] TRUE
```

```r
# is.numeric only gives TRUE with numbers

is.numeric(c(5:10,"chocolate"))
```

```
## [1] FALSE
```

```r
# flip using ! to make it say TRUE with not numbers

!is.numeric(c(5:10,"chocolate"))
```

```
## [1] TRUE
```



```r
# Test rescale2 with characters; input isn't numeric, rescale2 tells you it should be (author wrote that themselves)

# rescale2(c(5:10, "barry"))
```


**##Write Funxn to Find NAs**

1st make simple input where answer is known


```r
x <- c(1, 2, NA, 3, NA)
y <- c(NA, 3, NA, 3, 4)

# Test is.na

is.na(x)
```

```
## [1] FALSE FALSE  TRUE FALSE  TRUE
```

```r
# Using which() on is.na will give the positions where NA = TRUE

which(is.na(x))
```

```
## [1] 3 5
```

```r
which(is.na(y))
```

```
## [1] 1 3
```

```r
# Do more stuff

is.na(x)
```

```
## [1] FALSE FALSE  TRUE FALSE  TRUE
```

```r
is.na(y)
```

```
## [1]  TRUE FALSE  TRUE FALSE FALSE
```

```r
is.na(x) & is.na(y)
```

```
## [1] FALSE FALSE  TRUE FALSE FALSE
```

```r
# use sum() to sum values; FALSEs(=0) and TRUEs (=1) are numbers

sum(is.na(x) & is.na(y))
```

```
## [1] 1
```


**#Create Funxn
Make funxn using **sum()** and **is.na()** from before

```r
both_na <- function(x,y){
  sum( is.na(x) & is.na(y) )
}
```


Test the funxn we just created.

```r
both_na(x,y)
```

```
## [1] 1
```

Add new variables from slides; and test funxn.

```r
x <- c(NA, NA, NA)
y1 <- c( 1, NA, NA)
y2 <- c( 1, NA, NA, NA)

# Use both_na on x and y1

both_na( x, y1 )
```

```
## [1] 2
```

```r
# use both_na on x and y2

both_na( x, y2)
```

```
## Warning in is.na(x) & is.na(y): longer object length is not a multiple of
## shorter object length
```

```
## [1] 3
```



```r
x2 <- c( NA, NA )

both_na( x2, y2 )
```

```
## [1] 3
```

Showing how stuff gets recycled or something

```r
plot(1:10, col = c("red", "blue", "green"))
```

![](Class07_files/figure-html/unnamed-chunk-10-1.png)<!-- -->

Mess around w it

```r
x <- c( NA, NA, NA)
y3 <- c(1, NA, NA, NA, NA, NA, NA)

both_na(x,y3)
```

```
## Warning in is.na(x) & is.na(y): longer object length is not a multiple of
## shorter object length
```

```
## [1] 6
```

What happens is the first (shorter vector, will recycle itself to become as long as the @nd, longer vector)


```r
x3 <- c(1, NA, NA)
#what happens is
#x3<- c(1, NA, NA, 1,  NA, NA, 1)
y3 <- c(1, NA, NA, NA, NA, NA, NA)
```




```r
length(x)
```

```
## [1] 3
```

```r
length(y3)
```

```
## [1] 7
```

Make **both_na2** in conjunction with the **stop()** funxn to allow for an error message

```r
# Create both_na2

both_na2 <- function(x,y){
  if (length(x) != length(y) ) {
    stop("Inputs x and y should be the same length")
  }
}

# Test funxn

#both_na2(x, y3)
```
##Write grading funxn


```r
# Student 1 scores

a <- c(100, 100, 100, 100, 100, 100, 100, 90)

# Student 2 scores

b <- c(100, NA, 90, 90, 90, 90, 97, 80)

# using - gives everything BUT the value the original funxn gives
a[-which.min(a)]
```

```
## [1] 100 100 100 100 100 100 100
```

```r
b[-which.min(b)]
```

```
## [1] 100  NA  90  90  90  90  97
```

```r
# Do mean of it
mean(a[ -which.min(a) ])
```

```
## [1] 100
```

```r
mean(b[ -which.min(b) ], na.rm = TRUE)
```

```
## [1] 92.83333
```
check how many NAs there are


```r
any(is.na(b))
```

```
## [1] TRUE
```


**#Create the GODLIEST funxn EVER**


```r
a <- c(100, 100, 100, 100, 100, 100, 100, 90)
grade <- function(x){
  if(any(is.na(x))) {
    warning("Student missing HW")
    
  }
  mean(x[ -which.min(x) ], na.rm = TRUE)
  
}

grade(a)
```

```
## [1] 100
```

```r
grade(b)
```

```
## Warning in grade(b): Student missing HW
```

```
## [1] 92.83333
```

Grade whole class

```r
class <- read.csv("https://tinyurl.com/gradeinput", row.names = 1)

# Use apply, 1 = row, 2 would be col; apply(x, margin, FUN)
apply(class, 1, grade)
```

```
## Warning in FUN(newX[, i], ...): Student missing HW

## Warning in FUN(newX[, i], ...): Student missing HW

## Warning in FUN(newX[, i], ...): Student missing HW

## Warning in FUN(newX[, i], ...): Student missing HW
```

```
##  student-1  student-2  student-3  student-4  student-5  student-6 
##   91.75000   82.50000   84.25000   88.00000   88.25000   89.00000 
##  student-7  student-8  student-9 student-10 student-11 student-12 
##   94.00000   93.75000   87.75000   81.33333   86.00000   91.75000 
## student-13 student-14 student-15 student-16 student-17 student-18 
##   92.25000   87.75000   83.33333   89.50000   88.00000   97.00000 
## student-19 student-20 
##   82.75000   82.75000
```


datapasta Rcade


```r
# Rcade
# install.packages("devtools")
# devtools::install_github('RLesur/Rcade')

# USAGE
#  Rcade::games
#  Rcade::games$Pacman
```






















