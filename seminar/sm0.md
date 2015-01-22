# sem00
Eva  
January 21, 2015  


**Read dataset in R**

```r
prDat <- read.table("GSE4051_MINI.tsv", header = TRUE, row.names = 1)
str(prDat)
```

```
## 'data.frame':	39 obs. of  6 variables:
##  $ sidNum    : int  20 21 22 23 16 17 6 24 25 26 ...
##  $ devStage  : Factor w/ 5 levels "4_weeks","E16",..: 2 2 2 2 2 2 2 4 4 4 ...
##  $ gType     : Factor w/ 2 levels "NrlKO","wt": 2 2 2 2 1 1 1 2 2 2 ...
##  $ crabHammer: num  10.22 10.02 9.64 9.65 8.58 ...
##  $ eggBomb   : num  7.46 6.89 6.72 6.53 6.47 ...
##  $ poisonFang: num  7.37 7.18 7.35 7.04 7.49 ...
```

**How many rows are there? Hint: nrow(), dim().**

```r
nrow(prDat)
```

```
## [1] 39
```

**How many columns or variables are there? Hint: ncol(), length(), dim().**

```r
ncol(prDat)
```

```
## [1] 6
```

**Inspect the first few observations or the last few or a random sample. Hint: head(), tail(), x[i, j] combined with sample().**

```r
head(prDat)
```

```
##           sidNum devStage gType crabHammer eggBomb poisonFang
## Sample_20     20      E16    wt     10.220   7.462      7.370
## Sample_21     21      E16    wt     10.020   6.890      7.177
## Sample_22     22      E16    wt      9.642   6.720      7.350
## Sample_23     23      E16    wt      9.652   6.529      7.040
## Sample_16     16      E16 NrlKO      8.583   6.470      7.494
## Sample_17     17      E16 NrlKO     10.140   7.065      7.005
```

```r
tail(prDat)
```

```
##           sidNum devStage gType crabHammer eggBomb poisonFang
## Sample_38     38  4_weeks    wt      9.767   6.608      7.329
## Sample_39     39  4_weeks    wt     10.200   7.003      7.320
## Sample_11     11  4_weeks NrlKO      9.677   7.204      6.981
## Sample_12     12  4_weeks NrlKO      9.129   7.165      7.350
## Sample_2       2  4_weeks NrlKO      9.744   7.107      7.075
## Sample_9       9  4_weeks NrlKO      9.822   6.558      7.043
```

```r
sample(prDat[2, ]) # all variables for Sample_21
```

```
##           sidNum poisonFang gType eggBomb crabHammer devStage
## Sample_21     21      7.177    wt    6.89      10.02      E16
```

**What does row correspond to – different genes or different mice?**

```r
rownames(prDat)
```

```
##  [1] "Sample_20" "Sample_21" "Sample_22" "Sample_23" "Sample_16"
##  [6] "Sample_17" "Sample_6"  "Sample_24" "Sample_25" "Sample_26"
## [11] "Sample_27" "Sample_14" "Sample_3"  "Sample_5"  "Sample_8" 
## [16] "Sample_28" "Sample_29" "Sample_30" "Sample_31" "Sample_1" 
## [21] "Sample_10" "Sample_4"  "Sample_7"  "Sample_32" "Sample_33"
## [26] "Sample_34" "Sample_35" "Sample_13" "Sample_15" "Sample_18"
## [31] "Sample_19" "Sample_36" "Sample_37" "Sample_38" "Sample_39"
## [36] "Sample_11" "Sample_12" "Sample_2"  "Sample_9"
```
Answer: Different mice

**What are the variable names? Hint: names(), dimnames().**

```r
names(prDat)
```

```
## [1] "sidNum"     "devStage"   "gType"      "crabHammer" "eggBomb"   
## [6] "poisonFang"
```

```r
dimnames(prDat)
```

```
## [[1]]
##  [1] "Sample_20" "Sample_21" "Sample_22" "Sample_23" "Sample_16"
##  [6] "Sample_17" "Sample_6"  "Sample_24" "Sample_25" "Sample_26"
## [11] "Sample_27" "Sample_14" "Sample_3"  "Sample_5"  "Sample_8" 
## [16] "Sample_28" "Sample_29" "Sample_30" "Sample_31" "Sample_1" 
## [21] "Sample_10" "Sample_4"  "Sample_7"  "Sample_32" "Sample_33"
## [26] "Sample_34" "Sample_35" "Sample_13" "Sample_15" "Sample_18"
## [31] "Sample_19" "Sample_36" "Sample_37" "Sample_38" "Sample_39"
## [36] "Sample_11" "Sample_12" "Sample_2"  "Sample_9" 
## 
## [[2]]
## [1] "sidNum"     "devStage"   "gType"      "crabHammer" "eggBomb"   
## [6] "poisonFang"
```

**What “flavor” is each variable, i.e. numeric, character, factor? Hint: str().**

```r
str(prDat)
```

```
## 'data.frame':	39 obs. of  6 variables:
##  $ sidNum    : int  20 21 22 23 16 17 6 24 25 26 ...
##  $ devStage  : Factor w/ 5 levels "4_weeks","E16",..: 2 2 2 2 2 2 2 4 4 4 ...
##  $ gType     : Factor w/ 2 levels "NrlKO","wt": 2 2 2 2 1 1 1 2 2 2 ...
##  $ crabHammer: num  10.22 10.02 9.64 9.65 8.58 ...
##  $ eggBomb   : num  7.46 6.89 6.72 6.53 6.47 ...
##  $ poisonFang: num  7.37 7.18 7.35 7.04 7.49 ...
```

**For sample, do a sanity check that each integer between 1 and the number of rows in the dataset occurs exactly once. Hint: a:b, seq(), seq_len(), sort(), table(), ==, all(), all.equal(), identical().**

```r
sort(prDat$sidNum) # sort ascending sample number
```

```
##  [1]  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
## [24] 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39
```

```r
seq_len(nrow(prDat)) # number of rows in prDat
```

```
##  [1]  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
## [24] 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39
```

```r
all(sort(prDat$sidNum) == seq_len(nrow(prDat))) # is sample number equal to number of rows?
```

```
## [1] TRUE
```

```r
identical(sort(prDat$sidNum), seq_len(nrow(prDat)))
```

```
## [1] TRUE
```

**For each factor variable, what are the levels? Hint: levels(), str().**

```r
levels(prDat$devStage) # levels for devStage
```

```
## [1] "4_weeks" "E16"     "P10"     "P2"      "P6"
```

```r
levels(prDat$gType) # levels for gType
```

```
## [1] "NrlKO" "wt"
```

**How many observations do we have for each level of devStage? For gType? Hint: summary(), table().**

```r
summary(prDat$devStage)
```

```
## 4_weeks     E16     P10      P2      P6 
##       8       7       8       8       8
```

```r
summary(prDat$gType)
```

```
## NrlKO    wt 
##    19    20
```

**Perform a cross-tabulation of devStage and gType. Hint: table().**

```r
table(prDat$devStage, prDat$gType)
```

```
##          
##           NrlKO wt
##   4_weeks     4  4
##   E16         3  4
##   P10         4  4
##   P2          4  4
##   P6          4  4
```

```r
addmargins(with(prDat, table(devStage, gType)))
```

```
##          gType
## devStage  NrlKO wt Sum
##   4_weeks     4  4   8
##   E16         3  4   7
##   P10         4  4   8
##   P2          4  4   8
##   P6          4  4   8
##   Sum        19 20  39
```

**If you had to take a wild guess, what do you think the intended experimental design was? What actually happened in real life?**  Four mice with each genotype were sacrificed at five different developmental stages to evaluate expression of three different genes. One NrlKO mice may have died or ran away before E16 developmental stage. 

**For each quantitative variable, what are the extremes? How about average or median? Hint: min(), max(), range(), summary(), fivenum(), mean(), median(), quantile().**

For crabHammer: 

```r
min(prDat$crabHammer) 
```

```
## [1] 8.214
```

```r
max(prDat$crabHammer) 
```

```
## [1] 10.34
```

```r
range(prDat$crabHammer) 
```

```
## [1]  8.214 10.340
```

```r
summary(prDat$crabHammer) 
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##   8.214   8.938   9.611   9.428   9.830  10.340
```

```r
fivenum(prDat$crabHammer)
```

```
## [1]  8.214  8.938  9.611  9.830 10.340
```

```r
mean(prDat$crabHammer)
```

```
## [1] 9.427821
```

```r
median(prDat$crabHammer)
```

```
## [1] 9.611
```

```r
quantile(prDat$crabHammer)
```

```
##     0%    25%    50%    75%   100% 
##  8.214  8.938  9.611  9.830 10.340
```

**Create a new data.frame called weeDat only containing observations for which expression of poisonFang is above 7.5.**

```r
(weedat <-subset(prDat, poisonFang > 7.5))
```

```
##           sidNum devStage gType crabHammer eggBomb poisonFang
## Sample_24     24       P2    wt      8.869   6.587      7.508
## Sample_26     26       P2    wt      9.611   6.870      7.511
## Sample_27     27       P2    wt      8.613   6.800      7.843
## Sample_8       8       P2 NrlKO      9.116   6.264      8.016
## Sample_7       7       P6 NrlKO      8.803   6.188      7.754
## Sample_33     33      P10    wt      9.004   7.082      8.086
## Sample_34     34      P10    wt      8.519   6.757      8.584
## Sample_15     15      P10 NrlKO      9.746   7.226      7.786
## Sample_19     19      P10 NrlKO      9.771   7.081      7.586
```

**For how many observations poisonFang > 7.5? How do they break down by genotype and developmental stage?**

```r
addmargins(with(weedat, table(devStage, gType)))
```

```
##          gType
## devStage  NrlKO wt Sum
##   4_weeks     0  0   0
##   E16         0  0   0
##   P10         2  2   4
##   P2          1  3   4
##   P6          1  0   1
##   Sum         4  5   9
```

**Print the observations with row names “Sample_16” and “Sample_38” to screen, showing only the 3 gene expression variables.**

```r
prDat[c("Sample_16", "Sample_38"), c("crabHammer", "eggBomb", "poisonFang")]
```

```
##           crabHammer eggBomb poisonFang
## Sample_16      8.583   6.470      7.494
## Sample_38      9.767   6.608      7.329
```

**Which samples have expression of eggBomb less than the 0.10 quantile?**

```r
quantileVal <- quantile(prDat$eggBomb, 0.1)
prDat[prDat$eggBomb < quantileVal, 1] 
```

```
## [1] 25 14  3 35
```

