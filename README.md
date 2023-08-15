# Biologically-Supervised-Graph-Embedding
In this work, we introduce a novel approach named Biologically Supervised Graph Embedding (BSE) to allow for selecting most relevant features to enhance the prediction accuracy of comorbid disease pairs. Our investigation into BSE's impact on both centered and uncentered embedding methods showcases its consistent superiority over the state-of-the-art techniques and its adeptness in selecting dimensions enriched with vital biological insights, thereby improving prediction performance significantly, up to 50% when measured by ROC for some variations. Further analysis indicates that BSE consistently and substantially improves the ratio of disease associations to gene connectivity, affirming its potential in uncovering latent biological factors affecting comorbidity. The statistically significant enhancements across diverse metrics underscore BSE's potential to introduce novel avenues for precise disease comorbidity predictions and other potential applications. 


# Examples
```pip install -r requirements.txt```

* run BSE with uncentered embedding method "vect" for RR0 dataset 

  ```python3 BSE.py RR0 vect```

* run BSE with ucentered embedding method Isomap for RR1 dataset 

  ```python3 BSE.py RR1 iso```

#  results
## METRIC SCORES 

| Metric	| iso_r	 | iso_s  | p_val	 | std	  | emb_r  | emb_s  | p_val    | std    | vect_r | vect_s |	p_val    |	std   |
| --------- |:------:|:------:|:--------:|:------:|:------:|:------:|:--------:|:------:|:------:|:------:|:--------:|:------:|
| precision	| 0.8644 | 0.9046 | 3.71E-09 | 0.0057 |	0.8260 | 0.9074 | 8.63E-12 | 0.0059 | 0.8352 | 0.9075 |	7.02E-11 | 0.0066 |
| recall	| 0.9846 | 0.9717 | 5.92E-05 | 0.0058 |	1.0000 | 0.9742	| 2.25E-08 | 0.0045 | 0.9977 | 0.9743 |	6.49E-07 | 0.0061 |
| f1	    | 0.9206 | 0.9369 | 1.05E-06 | 0.0045 |	0.9047 | 0.9396	| 2.15E-09 | 0.0047 | 0.9093 | 0.9397 |	1.81E-08 | 0.0052 |
| Accuracy	| 0.8596 | 0.8919 | 3.69E-07 | 0.0078 |	0.8260 | 0.8966	| 5.59E-10 | 0.0081 | 0.8355 | 0.8967 |	5.35E-09 | 0.0091 |
| roc_auc	| 0.6255 | 0.7424 | 4.43E-09 | 0.0170 |	0.5000 | 0.7511	| 5.13E-12 | 0.0172 | 0.5315 | 0.7512 |	5.59E-11 | 0.0196 |


### Real Datasets
The tables below show the comparison of GEE and GEE_Sparse on the real datasets. 
<table>
<thead>
  <tr>
    <th rowspan="2">Data Set</th>
    <th colspan="2">Lap = T, Diag = T, Cor = T</th>
    <th colspan="2">Lap = T, Diag = T, Cor = F</th>
    <th colspan="2">Lap = T, Diag = F, Cor = T</th>
    <th colspan="2">Lap = T, Diag = F, Cor = F</th>
  </tr>
  <tr>
    <th>GEE(s)</th>
    <th>GEE_sparse(s)</th>
    <th>GEE(s)</th>
    <th>GEE_sparse(s)</th>
    <th>GEE(s)</th>
    <th>GEE_sparse(s)</th>
    <th>GEE(s)</th>
    <th>GEE_sparse(s)</th>
  </tr>
</thead>
<tbody>
  <tr>
    <td>CiteSeer</td>
    <td>0.096</td>
    <td><strong>0.040</strong></td>
    <td>0.097</td>
    <td><strong>0.032</strong></td>    
    <td>0.053</td>
    <td><strong>0.032</strong></td>    
    <td>0.051</td>
    <td><strong>0.030</strong></td>        
  </tr>
  <tr>
    <td>Cora</td>
    <td>0.083</td>
    <td><strong>0.028</strong></td>    
    <td>0.088</td>
    <td><strong>0.027</strong></td>        
    <td>0.068</td>
    <td><strong>0.033</strong></td>            
    <td>0.068</td>
    <td><strong>0.026</strong></td>              
  </tr>
  <tr>
    <td>Proteins-all</td>
    <td>2.268</td>
    <td><strong>0.438</strong></td>      
    <td>2.259</td>
    <td><strong>0.419</strong></td>          
    <td>1.866</td>
    <td><strong>0.391</strong></td>          
    <td>1.846</td>
    <td><strong>0.478</strong></td>      
  </tr>
  <tr>
    <td>PubMed</td>
    <td>0.776</td>
    <td><strong>0.194</strong></td>       
    <td>0.739</td>
    <td><strong>0.201</strong></td>       
    <td>0.673</td>
    <td><strong>0.208</strong></td>       
    <td>0.560</td>
    <td><strong>0.199</strong></td>   
  </tr>
  <tr>
    <td>CL-100K-1d8-L9</td>
    <td>4.926</td>
    <td><strong>1.091</strong></td>   
    <td>4.850</td>
    <td><strong>0.978</strong></td>   
    <td>4.166</td>
    <td><strong>0.992</strong></td>   
    <td>3.823</td>
    <td><strong>1.095</strong></td>   
  </tr>
  <tr>
    <td>CL-100K-1d8-L5</td>
    <td>604.018</td>
    <td><strong>174.552</strong></td>   
    <td>585.288</td>
    <td><strong>147.790</strong></td>   
    <td>633.746</td>
    <td><strong>118.705</strong></td>   
    <td>571.360</td>
    <td><strong>123.691</strong></td>   
  </tr>
</tbody>
</table>



<table>
<thead>
  <tr>
    <th rowspan="2">Data Set</th>
    <th colspan="2">Lap = F, Diag = T, Cor = T</th>
    <th colspan="2">Lap = F, Diag = T, Cor = F</th>
    <th colspan="2">Lap = F, Diag = F, Cor = T</th>
    <th colspan="2">Lap = F, Diag = F, Cor = F</th>
  </tr>
  <tr>
    <th>GEE(s)</th>
    <th>GEE_sparse(s)</th>
    <th>GEE(s)</th>
    <th>GEE_sparse(s)</th>
    <th>GEE(s)</th>
    <th>GEE_sparse(s)</th>
    <th>GEE(s)</th>
    <th>GEE_sparse(s)</th>
  </tr>
</thead>
<tbody>
  <tr>
    <td>CiteSeer</td>
    <td><strong>0.024</strong></td>
    <td>0.031</td>
    <td><strong>0.024</strong></td>
    <td>0.033</td>
    <td><strong>0.016</strong></td>    
    <td>0.034</td>
    <td><strong>0.014</strong></td>
    <td>0.031</td>
  </tr> 
  <tr>
    <td>Cora</td>
    <td><strong>0.024</strong></td>
    <td>0.025</td>
    <td><strong>0.026</strong></td>    
    <td><strong>0.026</strong></td>   
    <td><strong>0.023</strong></td>   
    <td>0.024</td>
    <td><strong>0.019</strong></td>      
    <td>0.025</td>
  </tr>
  <tr>
    <td>Proteins-all</td>
    <td>0.830</td>
    <td><strong>0.399</strong></td>      
    <td>0.623</td>
    <td><strong>0.411</strong></td>      
    <td>1.143</td>
    <td><strong>0.432</strong></td>      
    <td>0.518</td>
    <td><strong>0.462</strong></td>      
  </tr>
  <tr>
    <td>PubMed</td>
    <td>0.231</td>
    <td><strong>0.201</strong></td>        
    <td>0.201</td>
    <td>0.228</td>
    <td><strong>0.185</strong></td>    
    <td><strong>0.170</strong></td>    
    <td>0.188</td>
    <td><strong>0.177</strong></td>    
    <td>0.183</td>
  </tr>
  <tr>
    <td>CL-100K-1d8-L9</td>
    <td>1.330</td>
    <td><strong>0.909</strong></td>        
    <td>1.322</td>
    <td><strong>0.936</strong></td>    
    <td>1.114</td>
    <td><strong>0.924</strong></td>    
    <td><strong>1.058</strong></td>    
    <td>1.360</td>
  </tr>
  <tr>
    <td>CL-100K-1d8-L5</td>
    <td>203.860</td>
    <td><strong>108.977</strong></td>       
    <td>192.780</td>
    <td><strong>132.160</strong></td>       
    <td>171.838</td>
    <td><strong>125.935</strong></td>       
    <td>171.714</td>
    <td><strong>106.264</strong></td>   
  </tr>
</tbody>
</table>


