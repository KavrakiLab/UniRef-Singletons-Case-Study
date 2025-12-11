"""
Protocol on how to call GOINTOSIM for calculating pairwise GO Term similarity.

"""
import multiprocessing as mp
import numpy as np
import pandas as pd
from goatools.base import get_godag
from tqdm import tqdm
from lib.go_into_sim.Similarity import Similarity_of_Two_GOTerms

def linclust_protocol_parallel(uniref_clusters_with_go_term_members, go_path, processes=120):
    with mp.Pool(processes=processes, initargs=(go_path,),initializer=init_worker) as pool:
        results = list(tqdm(
            pool.imap(process_row, [(row, key) for key, row in uniref_clusters_with_go_term_members.iterrows()]),
            total=len(uniref_clusters_with_go_term_members)
        ))
    data = {
        "cluster_name": [index for index, row in uniref_clusters_with_go_term_members.iterrows()],
        "consistency": results
    }
    return pd.DataFrame(data)


def init_worker(go_path):
    global GO_DAG
    GO_DAG = get_godag(go_path, optional_attrs={'relationship'})

def process_row(row_tuple):
    data = row_tuple[0]["go_terms"]
    hit = 0
    name = row_tuple[1].split("_")[-1]
    for index, member in enumerate(row_tuple[0]["members"]):
        if member == name:
            hit = index
            break
    go_term_list_A = data[hit]
    clust_similarities = [linclust_consistency(go_term_list_A, go_term_list_B)
                          for index, go_term_list_B in enumerate(data)]# if not index == hit]
    return np.array(clust_similarities)


def linclust_consistency(go_term_list_A, go_term_list_B):
    GO_DAG = get_godag("/home/Users/fmq1/PycharmProjects/UniRefAnalysis/go-basic.obo", optional_attrs={'relationship'})
    def valid_term(term):
        """Return True if term exists in GO_DAG, else False."""
        try:
            _ = GO_DAG[term]
            return True
        except KeyError:
            return False

    def max_similarity(term, other_terms):
        """Return the max similarity score of `term` against a list of terms."""
        best = 0.0
        for other in other_terms:
            if not valid_term(other):
                continue
            val = float(go_term_similarity(term, other, GO_DAG))
            if val > best:
                best = val
        return best

    # A → B
    valid_A = [t for t in go_term_list_A if valid_term(t)]
    sum_A = sum(max_similarity(t, go_term_list_B) for t in valid_A)
    # B → A
    valid_B = [t for t in go_term_list_B if valid_term(t)]
    sum_B = sum(max_similarity(t, go_term_list_A) for t in valid_B)

    denom = len(valid_A) + len(valid_B)
    if len(valid_A) == 0 or len(valid_B) == 0:
        set1, set2 = set(go_term_list_A), set(go_term_list_B)
        union = set1 | set2
        return len(set1 & set2) / len(union) if union else 0.0
    return (sum_A + sum_B) / denom

def go_term_similarity(go_term_1, go_term_2, go_dag):
    return Similarity_of_Two_GOTerms(go_term_1, go_term_2, go_dag, 'GOGO')
