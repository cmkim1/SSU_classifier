# SSU_classifier
23_1_current_biotech  
  
## 1. Data Collecting
* Data downloaded from [here](https://github.com/yphsieh/16S-ITGDB/tree/master/data)
* contains 486,640 16S sequqences and each taxa
* number of phyla: 317
   
|frequency more than| > 100 | > 1000 | > 10000 |
|-----|----|-----|------|
|number of phyla| 119 | 38 | 5 |

![스크린샷 2023-06-28 오후 5 47 02](https://github.com/cmkim1/SSU_classifier/assets/119988478/0935f846-b4f6-4315-be31-7dbcece006e1)  
  
## 2. Data Handling
* For example1, 20,000 sequences were randomly selected and each nucleotide mono-, di-, trimer frequencies were merged.
   
|firmicutes|others|
|-----|-----|
|4627|15373|  

* For example2, 20,000 sequences were randomly selected and each nucleotide pentamer frequencies were merged.
   
|firmicutes|others|
|-----|-----|
|4659|15340|  

* For example3, 20,000 sequences were randomly selected and each nucleotide pentamer frequencies were merged. Then highly correlated pentamer sequences were removed.
   
|firmicutes|others|       | data |collapsed data|
|-----|-----|-----|----|-----|
|4569|15431|      |1024|147|  

* For finding example, 1,000 seuqences of each phylum which has more than 1,000 sequences were randomly selected and other 1000 sequences were added.
  
## 3. Used Model - All process were worked in R environment
* Logistic regression
* LASSO
* Tree
* Random Forest
* Boosting
  
## 4. Results
#### Example1: mono-, di-, trimer
  * Summary
    |Model|time (s)|auc|bin_dev|
    |-----|-----|-----|-----|
    |LOGISTIC|0.846|0.9894948|806.6|
    |LASSO|9.456|0.9899495|784.4|
    |RF|24.025|0.9941032|1144.1|
    |BOOSTING|27.32|0.9929736|680.5|
      
  * top5 significant sequqences: 1. TA 2. AT 3. GG 4. TG 5. TT


            
#### Example2: pentamer
  * Summary
    |Model|time (s)|auc|bin_dev|
    |-----|-----|-----|-----|
    |LOGISTIC|219.47|0.9990533|238.8|
    |LASSO|18.317|0.9997188|113.7|
    |RF|373.869|0.9996336|576.5|
    |BOOSTING|1079.283|0.9995689|139.3|

  * top5 significant sequences: 1. TGATA 2. ACCGA 3. GTTCT 4. TTACG 5. CTGGT

#### Example3: collapsed pentamer
  Correlation coefficient threshold was determined by Example2 LASSO data.
  Appropriate number of parameters was around 125.
  
  ```
  length(which(coef(data2_cvfit, s ="lambda.min")>0))
  ```
  ```
  ## [1] 114
  ```

  The target number of parameters was 150.  
  
  * Summary - comparison with Example2
    |Model|time (s)|auc|bin_dev|time (s)|auc|bin_dev|
    |-----|-----|-----|-----|-----|-----|-----|
    |used|-|All pentamer|-|-|Collapsed pentamer|-|
    |LOGISTIC|219.47|0.9990533|238.8|2.017|0.9987319|297.7|
    |LASSO|18.317|0.9997188|113.7|9.597|0.9989139|278.0|
    |RF|373.869|0.9996336|576.5|36.716|0.9994395|737.9|
    |BOOSTING|1079.283|0.9995689|139.3|23.708|0.9985275|321.9|

  Time required were reduced especially in LOGISTIC, RF, BOOSTING model.
  It has additional time for correlation step (9.264 sec), but plenty of time has been saved
```
tic("correlation")
cor_matrix <- cor(data[,c(2:1025)])
cor_time <- toc()
```
```
## correlation: 9.264 sec elapsed
```

  
#### Example4: finding phylum
  sequence used:
  [sequence file]([sequence-2.txt](https://github.com/cmkim1/SSU_classifier/files/11893250/sequence-2.txt))
  blast result: *Escherichia Coli* (Proteobacteria)  
  * Result
    
  |rank|Phylum|lambda.min|
  |-----|-----|-----|
  |1|Proteobacteria|0.02066498|
  |2|Chloroflexi|0.64248237|
  |3|Bdellovibrionota|0.65915294|
  |4|Fusobacteriota|0.69814215|
  |5|Myxococcota|0.891009659|
  
  
```
tic("finding")
...
...
...
finding_time <- toc()
```
```
## finding: 64.383 sec elapsed
```

