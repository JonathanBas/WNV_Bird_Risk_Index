Data and code to reproduce:

**Mapping the Bird Risk Index for West Nile virus in Europe and its relationship with disease occurrence in humans**

J. Bastard, R. Metras, B. Durand

* *code_1.R*: prediction of WNV seroprevalence in 150 bird species (Model 1), 1000 bootstraps
* *code_2.py*: computes 1000 maps of the WNV Bird Risk Index (BRI) using predictions from *code_1.R* and bird spatial distribution maps from *Lumbierres et al, 2022*
* *code_3.R*: computes maps of BRI's mean, coefficient of variation, P50 and bounds of the confidence interval
* *code_4.R*: aggregates BRI maps at the NUTS region level
* *def_fun_predict_PIT.R*: defines useful functions
* *code_5.R*: relationship between the BRI and WNV cases (Model 2)
