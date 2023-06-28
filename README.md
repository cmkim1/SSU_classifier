# SSU_classifier
23_1_current_biotech  
  
## 1. Data Collecting
* Data download from [here](https://github.com/yphsieh/16S-ITGDB/tree/master/data)
* contains 486,640 16S sequqences and each taxa
* number of phylum: 317
   
|frequency more than| > 100 | > 1000 | > 10000 |
|-----|----|-----|------|
|number of phylum| 119 | 38 | 5 |

![스크린샷 2023-06-28 오후 5 47 02](https://github.com/cmkim1/SSU_classifier/assets/119988478/0935f846-b4f6-4315-be31-7dbcece006e1)  
  
## 2. Data Handling
* For example1, 20,000 sequences were randomly selected and each nucleotide mono-, di-, trimer frequencies were merged for distinguishing firmicutes sequence.
   
|firmicutes|others|
|-----|-----|
|4627|15373|  

* For example2, 20,000 sequences were randomly selected and each nucleotide pentamer frequencies were merged for distinguishing firmicutes sequence.
   
|firmicutes|others|
|-----|-----|
|4659|15340|  

* For example3, 20,000 sequences were randomly selected and each nucleotide pentamer frequencies were merged. Then highly correlated pentamer sequences were removed for data collapsing.
   
|firmicutes|others|       | data |collapsed data|
|-----|-----|-----|----|-----|
|4569|15431|      |1024|147|  

* For finding example, 1,000 seuqences of each phylum which has more than 1,000 sequences were randomly selected and other 1000 sequences were added for raining.
  
## 3. Used Model - All process were worked in R environment
* Logistic regression
* LASSO
* Tree
* Random Forest
* Boosting
  
## 4. Results
* Example1
  * Summary
    |method|time (s)|auc|bin_dev|
    |-----|-----|-----|-----|
    |LOGISTIC|0.846|0.9894948|806.6|
    |LASSO|9.456|0.9899495|784.4|
    |RF|24.025|0.9941032|1144.1|
    |BOOSTING|27.32|0.9929736|680.5|
      
  * top5 significant sequqences
    1. TA 2. AT 3. GG 4. TG 5. TT
* Example2
  * Summary



