#!/usr/bin/env python3

import numpy as np
import pandas as pd
import networkx as nx


def gene_network(df, cutoff, weight="weight", distance=True):
    """ Create a gene network based on a pandas::DataFrame

    Parameters
    ----------
    df : pandas::DataFrame
        Weighted edge list of the network. The basic header should be as
        follows::
            |gene_1|weight|gene_2|chrom_1|chrom_2|
        Where `chrom_1' and `chrom_2' are the chromosomes of gene_1 and gene_2,
        respectively. Column `weight' is the co-expression of `gene_1' and
        `gene_2' and can be specified in the keyword variable `weight'. 
        
        All columns have `str` as dype except for `weight` which is a float.
    
    
    cutoff: 2-tuple of ints/floats
        Interval `(weight_min, weight_max)' defining the edge set.
    
    weight : str, default 'weight'
        Column name of weight.
    
    distance : boolean, default True
        If True then `df` should include the columns `start_1' and `start_2'
        where `start_1' is the starting position of `gene_1' in `chrom_1`.
    """

    basic_columns = {"gene_1", "gene_2", "chrom_1", "chrom_2"}
    columns_error = ("'{}' must be columns of `df' "
                     "pandas.DataFrame".format(basic_columns))

    if not basic_columns.issubset(df.columns):
        raise Exception(columns_error)

    df = df[(df[weight] >= cutoff[0]) & (df[weight] < cutoff[1])]
    # Create edge attributes
    # Cis/Trans Flag
    df["type"] = np.where((df['chrom_1'] == df['chrom_2']), "cis", "trans")
    # Chromosomes separated with a hyphen
    df["chroms"] = df[["chrom_1", "chrom_2"]].apply(lambda x: "-".join(x),
                                                    axis=1).values
    edge_attributes = ["type", "chroms"]
    # Distance (only for cis)
    if distance:
        df["distance"] = 0
        df.loc[(df["type"] == "cis") & (df.start_1 < df.start_2), 'distance'] = (
            df.loc[(df["type"] == "cis") &
                   (df.start_1 < df.start_2), 'start_2'] -
            df.loc[(df["type"] == "cis") & (df.start_1 < df.start_2), 'start_1'])
        df.loc[(df["type"] == "cis") & (df.start_2 < df.start_1), 'distance'] = (
            df.loc[(df["type"] == "cis") &
                   (df.start_2 < df.start_1), 'start_1'] -
            df.loc[(df["type"] == "cis") & (df.start_2 < df.start_1), 'start_2'])
        edge_attributes.append("distance")

    G = nx.from_pandas_edgelist(df,
                                source="gene_1",
                                target="gene_2",
                                edge_attr=edge_attributes)

    return G

df = pd.DataFrame({"gene_1": ["ENSG00000167074", "ENSG00000159189"],
                   "gene_2": ["ENSG00000174738", "ENSG00000173369"],
                   "weight": [0.577364, 0.568734],
                   "chrom_1": ["22", "1"],
                   "chrom_2": ["3", "1"],
                   "start_1": [41367333, 22643630],
                   "start_2": [23945260, 22652762]})


G = gene_network(df, cutoff=(0,1), distance=True)
