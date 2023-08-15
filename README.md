# Biologically-Supervised-Graph-Embedding
In this work, we introduce a novel approach named Biologically Supervised Graph Embedding (BSE) to allow for selecting most relevant features to enhance the prediction accuracy of comorbid disease pairs. Our investigation into BSE's impact on both centered and uncentered embedding methods showcases its consistent superiority over the state-of-the-art techniques and its adeptness in selecting dimensions enriched with vital biological insights, thereby improving prediction performance significantly, up to 50% when measured by ROC for some variations. Further analysis indicates that BSE consistently and substantially improves the ratio of disease associations to gene connectivity, affirming its potential in uncovering latent biological factors affecting comorbidity. The statistically significant enhancements across diverse metrics underscore BSE's potential to introduce novel avenues for precise disease comorbidity predictions and other potential applications. 


# Examples
```pip install -r requirements.txt```

* run BSE with uncentered embedding method "vect" for RR0 dataset 

  ```python3 BSE.py RR0 vect```

* run BSE with ucentered embedding method Isomap for RR1 dataset 

  ```python3 BSE.py RR1 iso```

#  results
### AVERAGE METRIC SCORES FOR RR0
The tables below show the comparison of GEE and GEE_Sparse on the real datasets. 
| Metric	| iso_r	 | iso_s  | p_val	 | std	  | emb_r  | emb_s  | p_val    | std    | vect_r | vect_s |	p_val    |	std   |
| --------- |:------:|:------:|:--------:|:------:|:------:|:------:|:--------:|:------:|:------:|:------:|:--------:|:------:|
| precision	| 0.8644 | 0.9046 | 3.71E-09 | 0.0057 |	0.8260 | 0.9074 | 8.63E-12 | 0.0059 | 0.8352 | 0.9075 |	7.02E-11 | 0.0066 |
| recall	| 0.9846 | 0.9717 | 5.92E-05 | 0.0058 |	1.0000 | 0.9742	| 2.25E-08 | 0.0045 | 0.9977 | 0.9743 |	6.49E-07 | 0.0061 |
| f1	    | 0.9206 | 0.9369 | 1.05E-06 | 0.0045 |	0.9047 | 0.9396	| 2.15E-09 | 0.0047 | 0.9093 | 0.9397 |	1.81E-08 | 0.0052 |
| Accuracy	| 0.8596 | 0.8919 | 3.69E-07 | 0.0078 |	0.8260 | 0.8966	| 5.59E-10 | 0.0081 | 0.8355 | 0.8967 |	5.35E-09 | 0.0091 |
| roc_auc	| 0.6255 | 0.7424 | 4.43E-09 | 0.0170 |	0.5000 | 0.7511	| 5.13E-12 | 0.0172 | 0.5315 | 0.7512 |	5.59E-11 | 0.0196 |





