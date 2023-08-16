"""
    python3 BSE.py arg1 arg2 arg3
    arg1: RR0 or RR1
    arg2: iso, emb or vect
    example: python3 BSE.py RR0 vect
"""
################################################################################
import sys
import networkx as nx
import random
import numpy as np
from sklearn.metrics import roc_curve, auc
from sklearn.model_selection import StratifiedKFold
#Import svm model
from sklearn import svm
#Import scikit-learn metrics module for accuracy calculation
from sklearn import metrics
from utils import (file_to_array, file_to_matrix,
                    matrix_to_file, array_to_file,
                    save_list_to_file, get_graph_from_file, get_gene_idx_dict_from_file)
from sklearn.manifold import Isomap
import time
from datetime import timedelta
################################################################################
"""
    main body of the program
"""

def main():

    # from command line arguments
    rr_dataset = sys.argv[1]
    method = sys.argv[2]
    
    #------------------------Pre_setting_block_start---------------------------#
    input_folder = 'input'

    edge_list_file_path = f'{input_folder}/interactom_edges.txt' # stores the edges for the largest connected component in human Interactome
    node_file_path = f'{input_folder}/interactom_nodes.txt'   # stores the nodes for the largest connected component in human Interactome
    diseasepath = f'{input_folder}/Disease_genes' # stores the disease gene sets


    input_folder_2 = f'{input_folder}/{rr_dataset}/{method}'
    train_file_path = f'{input_folder_2}/train_set.tsv' 
    test_file_path = f'{input_folder_2}/test_set.tsv' 

    #output_folder
    output_folder = f'output/{rr_dataset}/{method}'  


    # set dynamic output paths in the run function below.
    #------------------------Pre_setting_block_end-----------------------------#
    #------------------------main_steps_start----------------------------------#
    # 1. get graph nodes
    node_idx_dict = get_gene_idx_dict_from_file(node_file_path)
    
    # 2. get selected disease pairs
    train_set_dict, train_disease_genes_dict = get_disease_sets(train_file_path)
    test_set_dict, test_disease_genes_dict = get_disease_sets(test_file_path)

    # 3. diease_gene dict for all
    disease_genes_dict = {**train_disease_genes_dict, **test_disease_genes_dict}
    
    # 4. select the way to generate feature vectors
    feature_opt = "concat"

    # 5. provid a classifier 
    kernel = 'rbf'
    c = 3.5
    gamma = None
    clf = svm.SVC(kernel=kernel, C = c)  

    # 6. get  vectors
    vec_opt = method  # "vect"  -- svd, U, "emb" -- svd, U*Lambda, "iso", -- isomapm, this does not have evals.

    # 7. set final dimensions for selection
    select_max_eigen_num = 20

    # 8. Accuracy before BSE 
    evecs, evals, embedding = get_vecs(vec_opt, edge_list_file_path, input_folder_2)
    acc = acc_before_BSE(clf, select_max_eigen_num, train_set_dict, test_set_dict, node_idx_dict, evecs, embedding, disease_genes_dict)

    print(f'vec opt: {vec_opt}')
    print(f'dim: {select_max_eigen_num}')
    print("Before BSE:")
    print("Accuracy:", acc)

    # 8. stratified folds. The folds are made by preserving the percentage of samples for each class. 
    skf5 = StratifiedKFold(n_splits=5, random_state= None, shuffle = False) # 5 folds for cross-valid selection

    # 10. run BSE and output metric score
    time_s = time.time()
    if embedding is not None:       # for opt "emb"
        max_eigen_emb, max_idx_list = BSE(clf, skf5, train_set_dict, test_set_dict, node_idx_dict, select_max_eigen_num, embedding, evals, disease_genes_dict)
        time_e = time.time()
        train_disease_pairs, train_feature_vecs, train_labels = get_feature_vec_concat(train_set_dict, node_idx_dict, max_eigen_emb, disease_genes_dict)
        test_disease_pairs, test_feature_vecs, test_labels = get_feature_vec_concat(test_set_dict, node_idx_dict, max_eigen_emb, disease_genes_dict)

    else:
        max_eigen_vecs, max_idx_list = BSE(clf, skf5, train_set_dict, test_set_dict, node_idx_dict, select_max_eigen_num, evecs, evals, disease_genes_dict)
        time_e = time.time()
        train_disease_pairs, train_feature_vecs, train_labels = get_feature_vec_concat(train_set_dict, node_idx_dict, max_eigen_vecs, disease_genes_dict)
        test_disease_pairs, test_feature_vecs, test_labels = get_feature_vec_concat(test_set_dict, node_idx_dict, max_eigen_vecs, disease_genes_dict)

    time_tlt = time_e - time_s

    print("After BSE:")
    acc = train_and_test_get_acc(clf, train_feature_vecs, train_labels, test_feature_vecs, test_labels)
    print("Accuracy:", acc)
    print(f"time: {str(timedelta(seconds=time_tlt))}")

    
    # 9. save the selected dimensions to files
    save_list_to_file(max_idx_list, f'{output_folder}/max_idxs.tsv')

    #------------------------main_steps_end------------------------------------#
################################################################################
#####--------------------------functions_start-----------------------------#####
def BSE(clf_ori, skf, train_set_dict, test_set_dict, node_idx_dict, select_max_eigen_num, dim_vecs_ori, dim_vals_ori, disease_genes_dict):
    max_idx_list = []   # list importance from largest to smallest

    dim_idxs_ori = [i for i in range(dim_vecs_ori.shape[1])]
    dim_idx_random = list(np.random.permutation(dim_idxs_ori))

    auc_avg_ori = 0
    for j in range(select_max_eigen_num):
        print(f"loop_{j}--started")
        auc_avg_delta_list = []
        auc_avg_new_list = []

        for i in range(len(dim_idx_random)):
            ran_idx = dim_idx_random[i]
            vec_ran = dim_vecs_ori[:,ran_idx].reshape(-1,1)
            if j == 0:
                dim_vecs_ran = vec_ran
            else:
                dim_vecs_ran = np.hstack((max_eigen_vecs, vec_ran)) # add column by column
            
            train_disease_pairs, train_feature_vecs, train_labels = get_feature_vec_concat(train_set_dict, node_idx_dict, dim_vecs_ran, disease_genes_dict)
            auc_avg_new, clf = auc_for_train_disease_pairs_cross_valid(clf_ori, skf, train_feature_vecs, train_labels)
            auc_avg_delta = auc_avg_new - auc_avg_ori 
            auc_avg_delta_list.append(auc_avg_delta)
            auc_avg_new_list.append(auc_avg_new)


        # in case there are more than one dimension that is the max, randomly choose one from them if there are more than one
        idxs_of_max_idx = [idx for idx, auc_avg_delta in enumerate(auc_avg_delta_list) if auc_avg_delta == max(auc_avg_delta_list)]
        idx_of_max_idx = random.choice(idxs_of_max_idx)
        max_idx = dim_idx_random[idx_of_max_idx]
        max_dimension_vec = dim_vecs_ori[:,max_idx].reshape(-1,1) # reshape to has one column, then can do vstack
        max_dimension_idx = dim_idxs_ori[max_idx]
        if dim_vals_ori is not None:
            max_dimension_val = dim_vals_ori[max_idx]
        if j == 0:
            max_eigen_vecs = max_dimension_vec
        else:
            max_eigen_vecs = np.hstack((max_eigen_vecs, max_dimension_vec)) # add column by column
        
        max_idx_list.append(max_dimension_idx)
        dim_idx_random.pop(idx_of_max_idx)

        auc_avg_ori = auc_avg_new_list[idx_of_max_idx]

        print(f"loop_{j}--finished")

    
    return max_eigen_vecs, max_idx_list    

#------------------------------------------------------------------------------# 
def acc_before_BSE(clf, dim, train_set_dict, test_set_dict, node_idx_dict, evecs, embedding, disease_genes_dict):
    if embedding is not None:       # for opt "emb"
        embedding = embedding[:, :dim]
        train_disease_pairs, train_feature_vecs, train_labels = get_feature_vec_concat(train_set_dict, node_idx_dict, embedding, disease_genes_dict)
        test_disease_pairs, test_feature_vecs, test_labels = get_feature_vec_concat(test_set_dict, node_idx_dict, embedding, disease_genes_dict)        
    else:
        evecs = evecs[:, :dim]
        train_disease_pairs, train_feature_vecs, train_labels = get_feature_vec_concat(train_set_dict, node_idx_dict, evecs, disease_genes_dict)
        test_disease_pairs, test_feature_vecs, test_labels = get_feature_vec_concat(test_set_dict, node_idx_dict, evecs, disease_genes_dict)
    
    acc = train_and_test_get_acc(clf, train_feature_vecs, train_labels, test_feature_vecs, test_labels)
    return acc

#------------------------------------------------------------------------------# 
def get_vecs(vec_opt, edge_list_file_path, input_folder):
    if vec_opt == "iso":
        G_sub = get_graph_from_file(edge_list_file_path)
        A = nx.adjacency_matrix(G_sub)  # sparse matrix, reordered the nodes so that is from small to large

        embedding = Isomap(n_components=100)
        vecs = embedding.fit_transform(A)

        return vecs, None, None
    
    elif vec_opt == "vect":
        selected_evals = file_to_array(f"{input_folder}/selected_eigval.tsv")       
        selected_evecs = file_to_matrix(f"{input_folder}/selected_eigvec.tsv")          # this file has dimension 100
        return selected_evecs, selected_evals, None
    
    elif vec_opt == "emb":
        selected_evals = file_to_array(f"{input_folder}/selected_eigval.tsv")
        selected_evecs = file_to_matrix(f"{input_folder}/selected_eigvec.tsv")          # this file has dimension 100
        diag_eigen_vals = np.zeros((100, 100), float)
        np.fill_diagonal(diag_eigen_vals, selected_evals)
        selected_embedding = np.matmul(selected_evecs, diag_eigen_vals)

        return selected_evecs, selected_evals, selected_embedding
#------------------------------------------------------------------------------#    
def get_disease_sets(file_path):
    set_dict = {}   #{dis_pair: rr}
    disease_genes_dict = {}     #{disease: [gene_1, gene_2, ...]}

    f = open(file_path, "r")
    head = True
    for line in f:
        if head:
            head = False
            continue
        
        row = line.strip().split("\t")
        dis_pair, disease_a_genes, disease_b_genes, all_genes, rr = row

        disease_a, disease_b = dis_pair.split("&")

        set_dict[dis_pair] = int(rr)
        disease_genes_dict[disease_a] = disease_a_genes.split(",")
        disease_genes_dict[disease_b] = disease_b_genes.split(",")


    f.close()

    return set_dict, disease_genes_dict

#------------------------------------------------------------------------------#
def get_feature_vec_concat(disease_pair_dict, node_idx_dict, evecs, disease_genes_dict, average = False):
    """
        paper:
        generate feature vectors based on the selected eigen vectors 
        one disease pair one feature vector -- disease_a & disease_b: sum over feature vector column vise for gene a, gene b seperately,
        and then concatenate them into one feature vector 
        fi is the eigen vector for each gene, before selected, the feature vector for each gene is the eigen vector with size == dim_ori
    """
    disease_pairs = []
    feature_vecs = []
    labels = []

    First = True
    for disease_pair, rr in disease_pair_dict.items():
        
        disease_a, disease_b = disease_pair.split("&")

        gene_a_list = disease_genes_dict[disease_a]
        gene_b_list = disease_genes_dict[disease_b]

        feature_vec_a = disease_mass_vec(gene_a_list, node_idx_dict, evecs)
        feature_vec_b = disease_mass_vec(gene_b_list, node_idx_dict, evecs)

        if average:
            feature_vec_a /= len(gene_a_list)
            feature_vec_b /= len(gene_b_list)
        
        feature_vec = np.concatenate((feature_vec_a, feature_vec_b), axis=0)

        disease_pairs.append(disease_pair)
        
        if First:
            feature_vecs = feature_vec
            First = False
        else:
            feature_vecs = np.vstack((feature_vecs, feature_vec))
        
        labels.append(rr)
    labels = np.array(labels)
    
    return disease_pairs, feature_vecs, labels

#------------------------------------------------------------------------------#
def disease_mass_vec(gene_list, node_idx_dict, evecs):
    node_idices = [node_idx_dict[gene] for gene in gene_list if gene in node_idx_dict]
    gene_evecs = evecs[node_idices, :]

    feature_vec = gene_evecs.sum(axis = 0) # sum along the columns
    return feature_vec

#------------------------------------------------------------------------------#
def get_sublist(list_ori, indices):
    return [list_ori[i] for i in indices]  

#------------------------------------------------------------------------------#
def auc_for_train_disease_pairs_cross_valid(clf, skf, train_feature_vecs, train_labels):
    auc_list = []
    folds = skf.split(train_feature_vecs, train_labels)
    for fold_no, (train_idx, test_idx) in enumerate(folds):

        #Train the model using the training sets
        clf.fit(train_feature_vecs[train_idx], train_labels[train_idx])

        #Predict the response for test dataset
        y_pred = clf.predict(train_feature_vecs[test_idx])

        # valid_acc = metrics.accuracy_score(valid_labels, y_pred)

        false_positive_rate, true_positive_rate, thresholds = roc_curve(train_labels[test_idx], y_pred)
        roc_auc = auc(false_positive_rate, true_positive_rate)

        auc_list.append(roc_auc)

    auc_avg = sum(auc_list)/len(auc_list)

    return auc_avg, clf

#------------------------------------------------------------------------------#
def get_acc(clf, test_set, true_labels):
    #Predict the response for test dataset
    y_pred = clf.predict(test_set)
    return metrics.accuracy_score(true_labels, y_pred)  # Model Accuracy: how often is the classifier correct?

#------------------------------------------------------------------------------#
def train_and_test_get_acc(clf, train_feature_vecs, train_labels, test_feature_vecs, test_labels):
    #Train the model using the training sets
    clf.fit(train_feature_vecs, train_labels)

    #Predict the response for test dataset
    y_pred = clf.predict(test_feature_vecs)

    return metrics.accuracy_score(test_labels, y_pred)  

#####--------------------------functions_end-------------------------------#####
################################################################################
if __name__ == '__main__':
    main()