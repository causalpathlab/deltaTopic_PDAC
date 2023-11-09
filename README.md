# Analysis of PDAC data using deltaTopic and Bayesian ETM

## Model output
### DeltaTopic
***Directory: models/BDeltaTopic_allgenes_ep2000_nlv16_bs1024_combinebyadd_lr0.01_train_size1.0***

Report (needs to be knitr to get updated):

- v1: scaled kl_beta by sample size, upweighted kl_beta by a factor of 10, kl_beta has more effects in training
https://causalpathlab.github.io/deltaTopic_PDAC/models/BDeltaTopic_allgenes_ep2000_nlv32_bs1024_combinebyadd_lr0.01_train_size1.0_pip0rho_0.1_pip0delta_0.1_klbeta_10.0v1/Report.html

- Archived:
https://causalpathlab.github.io/deltaTopic_PDAC/models/BDeltaTopic_allgenes_ep2000_nlv32_bs1024_combinebyadd_lr0.01_train_size1.0_v1/Report.html

### Bayesian ETM(Spliced)

***Directory: models/BETM_spliced_ep2000_nlv32_bs1024_lr0.01_train_size0.9_pip0rho_0.1_klbeta_10.0v1***

### deltaETM
***Directory: models/TotalDeltaETM_allgenes_ep1000_nlv32_bs512_combinebyadd_lr0.01_train_size1***

|  filenames            |  descriptions                          | type                     | 
|-----------------------|----------------------------------------|--------------------------|
|   rho_weights.csv.gz  | shared topics loadings                 | genes loadings for topics|
|   delta_weights.csv.gz|   delta topic loadings                 | genes loadings for topics|
|   genes.csv.gz        |   gene symbols                         | annotations              |
|   samples.csv.gz      |   cell_id, sample_id, cancer_type, sex | annotations              |
|  topics.csv.gz         |  topic proportions for each cell       | topics proportions      |
|  topics_untrain.csv.gz|  topics before Softmax for each cell   | topics proportions       |


## Data(training) are stored on *Astrocytes*
***Directory: /data/PDAC/deltaTopic/QC***

|  filenames            |  descriptions                          | 
|-----------------------|----------------------------------------|
|  final_qc_correct.cols| columns                                | 
|  final_qc_correct.rows| rows                                   |
|  final_qc_correct.mtx | MTX                                    | 
