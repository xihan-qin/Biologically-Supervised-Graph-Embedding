# Biologically-Supervised-Graph-Embedding
In this work, we introduce a novel approach named Biologically Supervised Graph Embedding (BSE) to allow for selecting most relevant features to enhance the prediction accuracy of comorbid disease pairs. Our investigation into BSE's impact on both centered and uncentered embedding methods showcases its consistent superiority over the state-of-the-art techniques and its adeptness in selecting dimensions enriched with vital biological insights, thereby improving prediction performance significantly, up to 50% when measured by ROC for some variations. Further analysis indicates that BSE consistently and substantially improves the ratio of disease associations to gene connectivity, affirming its potential in uncovering latent biological factors affecting comorbidity. The statistically significant enhancements across diverse metrics underscore BSE's potential to introduce novel avenues for precise disease comorbidity predictions and other potential applications. 


# Run BSE
## Installation
```pip install -r requirements.txt```

## Examples
* run BSE with uncentered embedding method "emb" for RR0 dataset 

  ```python3 BSE.py RR0 emb```

* run BSE with ucentered embedding method Isomap for RR1 dataset 

  ```python3 BSE.py RR1 iso```


#  Results Showcase
For additional details, please consult the paper titled "Enhancing Disease Comorbidity Prediction Using Biologically Supervised Graph Embedding on the Human Interactome."
### AVERAGE METRIC SCORES FOR RR0
| Metric	| iso_r	 | iso_s  | p_val	 | std	  | emb_r  | emb_s  | p_val    | std    | vect_r | vect_s |	p_val    |	std   |
| --------- |:------:|:------:|:--------:|:------:|:------:|:------:|:--------:|:------:|:------:|:------:|:--------:|:------:|
| precision	| 0.8644 | 0.9046 | 3.71E-09 | 0.0057 |	0.8260 | 0.9074 | 8.63E-12 | 0.0059 | 0.8352 | 0.9075 |	7.02E-11 | 0.0066 |
| recall	| 0.9846 | 0.9717 | 5.92E-05 | 0.0058 |	1.0000 | 0.9742	| 2.25E-08 | 0.0045 | 0.9977 | 0.9743 |	6.49E-07 | 0.0061 |
| f1	    | 0.9206 | 0.9369 | 1.05E-06 | 0.0045 |	0.9047 | 0.9396	| 2.15E-09 | 0.0047 | 0.9093 | 0.9397 |	1.81E-08 | 0.0052 |
| Accuracy	| 0.8596 | 0.8919 | 3.69E-07 | 0.0078 |	0.8260 | 0.8966	| 5.59E-10 | 0.0081 | 0.8355 | 0.8967 |	5.35E-09 | 0.0091 |
| roc_auc	| 0.6255 | 0.7424 | 4.43E-09 | 0.0170 |	0.5000 | 0.7511	| 5.13E-12 | 0.0172 | 0.5315 | 0.7512 |	5.59E-11 | 0.0196 |

### AVERAGE METRIC SCORES FOR RR1
| Metric	| iso_r  | 	iso_s |	p_val    |	std	  | emb_r  | emb_s	| p_val    | std    | vect_r | vect_s |	p_val    | std    |
| --------- |:------:|:------:|:--------:|:------:|:------:|:------:|:--------:|:------:|:------:|:------:|:--------:|:------:|
| precision | 0.6975 | 0.7262 | 1.07E-05 | 0.0104 |	0.5859 | 0.7335 | 4.08E-11 | 0.0127 | 0.6456 | 0.7370 |	4.18E-09 | 0.0132 |
| recall    | 0.8250 | 0.8102 |	1.37E-02 | 0.0154 |	0.9900 | 0.8102 | 1.52E-10 | 0.0179 | 0.8649 | 0.8057 |	8.80E-07 | 0.0158 |
| f1        | 0.7558 | 0.7658 |	4.84E-03 | 0.0085 |	0.7361 | 0.7698 | 7.13E-06 | 0.0116 | 0.7393 | 0.7697 |	6.87E-06 | 0.0104 |
| Accuracy  | 0.6889 | 0.7109 |	1.21E-04 | 0.0108 |	0.5859 | 0.7173 | 3.37E-10 | 0.0143 | 0.6440 | 0.7187 |	5.29E-08 | 0.0144 |
| roc_auc   | 0.6616 | 0.6910 |	3.22E-05 | 0.0122 |	0.5048 | 0.6987 | 2.12E-11 | 0.0155 | 0.5997 | 0.7013 |	1.01E-08 | 0.0162 |

### Average metric scores along concatenated dimensions by BSE
![fig1](https://github.com/xihan-qin/Biologically-Supervised-Graph-Embedding/blob/main/plots/4_metrics_together_both_rrs.png)

### Top 20 genes with associated disease numbers for first 5 dimensions
![fig2](https://github.com/xihan-qin/Biologically-Supervised-Graph-Embedding/blob/main/plots/select_dims_top_20_genes_diseases_no_both_rrs.png)

### Top 20 genes and their degrees for the first 5 dimensions
![fig3](https://github.com/xihan-qin/Biologically-Supervised-Graph-Embedding/blob/main/plots/select_dims_top_20_genes_degrees_no_both_rrs.png)

### Disease count and degrees for top 20 genes and 5 dimensions
![fig4](https://github.com/xihan-qin/Biologically-Supervised-Graph-Embedding/blob/main/plots/disease_degree_all_methods_both_rr.png)