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
    
    
    cutoff: 2-tuple
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
