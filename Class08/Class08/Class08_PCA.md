---
title: "PCA Practice"
author: "Justin"
output:
  html_document:
    keep_md: yes
---



Using PCA

use `prcomp()` to do PCA


```r
# Download file

mydata <- read.csv("https://tinyurl.com/expression-CSV",
 row.names=1)

# Shows first 6 lines of data plot

head(mydata)
```

```
##        wt1 wt2  wt3  wt4 wt5 ko1 ko2 ko3 ko4 ko5
## gene1  439 458  408  429 420  90  88  86  90  93
## gene2  219 200  204  210 187 427 423 434 433 426
## gene3 1006 989 1030 1017 973 252 237 238 226 210
## gene4  783 792  829  856 760 849 856 835 885 894
## gene5  181 249  204  244 225 277 305 272 270 279
## gene6  460 502  491  491 493 612 594 577 618 638
```

```r
# Find dimensions with dim(); outputs Row, Col

dim(mydata)
```

```
## [1] 100  10
```

```r
# Do PCA

pca <- prcomp(t(mydata), scale=TRUE) 


# See what is returned by the prcomp() function
attributes(pca) 
```

```
## $names
## [1] "sdev"     "rotation" "center"   "scale"    "x"       
## 
## $class
## [1] "prcomp"
```

```r
# Basic PC1 v PC2 plot

plot(pca$x[,1], pca$x[,2])
```

![](Class08_PCA_files/figure-html/unnamed-chunk-1-1.png)<!-- -->

```r
# Variance captured per PC
pca.var <- pca$sdev^2 

# Last digit specifies number of decimal places the variance goes out to

pca.var.per <- round(pca.var/sum(pca.var)*100, 1) 

# Print, tells you how much variance there is per column; EX: 92.6% variance in column 1, etc

pca.var.per
```

```
##  [1] 92.6  2.3  1.1  1.1  0.8  0.7  0.6  0.4  0.4  0.0
```

```r
barplot(pca.var.per, main="Scree Plot",
 xlab="Principal Component", ylab="Percent Variation")
```

![](Class08_PCA_files/figure-html/unnamed-chunk-1-2.png)<!-- -->

```r
## A vector of colors for wt and ko samples

colvec <- colnames(mydata)
colvec[grep("wt", colvec)] <- "red"
colvec[grep("ko", colvec)] <- "blue"

plot(pca$x[,1], pca$x[,2], col=colvec, pch=16,
 xlab=paste0("PC1 (", pca.var.per[1], "%)"),
 ylab=paste0("PC2 (", pca.var.per[2], "%)")) 


# Click to identify which sample is which
identify(pca$x[,1], pca$x[,2], labels=colnames(mydata)) 
```

![](Class08_PCA_files/figure-html/unnamed-chunk-1-3.png)<!-- -->

```
## integer(0)
```



```r
# Import and Read the file

uk <- read.csv("UK_foods.csv", row.names = 1)
```



```r
# Using View(), head(), and tail()

View(uk)

head(uk)
```

```
##                England Wales Scotland N.Ireland
## Cheese             105   103      103        66
## Carcass_meat       245   227      242       267
## Other_meat         685   803      750       586
## Fish               147   160      122        93
## Fats_and_oils      193   235      184       209
## Sugars             156   175      147       139
```

```r
tail(uk)
```

```
##                   England Wales Scotland N.Ireland
## Fresh_fruit          1102  1137      957       674
## Cereals              1472  1582     1462      1494
## Beverages              57    73       53        47
## Soft_drinks          1374  1256     1572      1506
## Alcoholic_drinks      375   475      458       135
## Confectionery          54    64       62        41
```

```r
dim(uk)
```

```
## [1] 17  4
```

```r
# Barplot work; beside=T makes bars side by side; beside=F makes bars stacked

barplot(as.matrix(uk), beside=T, col=rainbow(nrow(uk)))
```

![](Class08_PCA_files/figure-html/unnamed-chunk-3-1.png)<!-- -->

```r
barplot(as.matrix(uk), beside=F, col=rainbow(nrow(uk)))
```

![](Class08_PCA_files/figure-html/unnamed-chunk-3-2.png)<!-- -->

```r
# Pairwise plot

pairs(uk, col=rainbow(10), pch=16)
```

![](Class08_PCA_files/figure-html/unnamed-chunk-3-3.png)<!-- -->

```r
pca <- prcomp(t(uk))
summary(pca)
```

```
## Importance of components:
##                             PC1      PC2      PC3       PC4
## Standard deviation     324.1502 212.7478 73.87622 4.189e-14
## Proportion of Variance   0.6744   0.2905  0.03503 0.000e+00
## Cumulative Proportion    0.6744   0.9650  1.00000 1.000e+00
```


```r
# plot

pca <- prcomp(t(uk))
summary(pca)
```

```
## Importance of components:
##                             PC1      PC2      PC3       PC4
## Standard deviation     324.1502 212.7478 73.87622 4.189e-14
## Proportion of Variance   0.6744   0.2905  0.03503 0.000e+00
## Cumulative Proportion    0.6744   0.9650  1.00000 1.000e+00
```

```r
plot((pca$x[,1]), (pca$x[,2]))
```

![](Class08_PCA_files/figure-html/unnamed-chunk-4-1.png)<!-- -->

