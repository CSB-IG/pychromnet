#!/usr/bin/env python3

import numpy as np
import pandas as pd
import networkx as nx


def table_format(edgelist, biomart):
    """Edit edgelist to correct table format for `gene_network` function.
    
    Default edgelist |gene_1|weight|gene_2| needs additional attributes, namely
    `start_1`, `start_2`, `chrom_1`, and `chrom_2`. Those attributes are taken
    from a `BioMart` file.
    
    Parameters
    ----------
    edgelist : str
        Path to edgelist. Tab separated format, 3 columns.
        HEADER:  |gene_1|weight|gene_2|
    

    biomart : str
        Path to BioMart file. Tab separated format, 8 columns.
        HEADER: |Gene stable ID|Chromosome/scaffold name|Gene start (bp)|\
        |Gene end (bp)|Gene % GC content|Gene type|HGNC symbol|NCBI gene ID|
    """

    df_bio = pd.read_csv(biomart, sep="\t")
    df_bio = df_bio.drop_duplicates(subset="Gene stable ID")

    df = pd.read_csv(edgelist, sep="\t", names=["gene_1", "weight", "gene_2"])

    # First gene
    df = pd.merge(df,
                  df_bio[[
                      "Chromosome/scaffold name", "Gene stable ID",
                      'Gene start (bp)', 'Gene end (bp)'
                  ]],
                  how="left",
                  left_on="gene_1",
                  right_on="Gene stable ID")
    df.drop("Gene stable ID", axis=1, inplace=True)
    df.rename(columns={
        "Chromosome/scaffold name": "chrom_1",
        'Gene start (bp)': "start_1",
        "Gene end (bp)": "end_1"
    },
              inplace=True)
    # Second gene
    df = pd.merge(df,
                  df_bio[[
                      "Chromosome/scaffold name", "Gene stable ID",
                      'Gene start (bp)', 'Gene end (bp)'
                  ]],
                  how="left",
                  left_on="gene_2",
                  right_on="Gene stable ID")
    df.drop("Gene stable ID", axis=1, inplace=True)
    df.rename(columns={
        "Chromosome/scaffold name": "chrom_2",
        'Gene start (bp)': "start_2",
        "Gene end (bp)": "end_2"
    },
              inplace=True)

    return df


def gene_network(df, cutoff, weight="weight", chromosome_level=False):
    """Create a gene network based on a pandas::DataFrame

    Parameters
    ----------
    df : pandas::DataFrame
        Weighted edge list of the network. The basic header should be as
        follows::
            |gene_1|weight|gene_2|chrom_1|chrom_2|start_1|start_2|
        Where `chrom_1' and `chrom_2' are the chromosomes of gene_1 and gene_2,
        respectively. Column `weight' is the co-expression of `gene_1' and
        `gene_2' and can be specified in the keyword variable `weight'. 
    
        `start_1` and `start_2` define the starting position of `gene_1` and
        `gene_2`, respectively.
        
        All columns have `str` as dtype except for `weight` which is a float.
    
    
    cutoff: 2-tuple of ints/floats
        Interval `(weight_min, weight_max)' defining the edge set.
    
    weight : str, default 'weight'
        Column name of weight.
    
    chromosome_level: boolean, default False
        If True collapse nodes $n_1$ and $n_2$ both in $gene_1$, for all $n_1$
        and $n_2$.
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
    edge_attributes = ["weight", "type", "chroms"]
    # Distance (only for cis)
    df["distance"] = 0
    df.loc[(df["type"] == "cis") & (df.start_1 < df.start_2), 'distance'] = (
        df.loc[(df["type"] == "cis") & (df.start_1 < df.start_2), 'start_2'] -
        df.loc[(df["type"] == "cis") & (df.start_1 < df.start_2), 'start_1'])
    df.loc[(df["type"] == "cis") & (df.start_2 < df.start_1), 'distance'] = (
        df.loc[(df["type"] == "cis") & (df.start_2 < df.start_1), 'start_1'] -
        df.loc[(df["type"] == "cis") & (df.start_2 < df.start_1), 'start_2'])
    edge_attributes.append("distance")

    if not chromosome_level:
        G_genes = nx.from_pandas_edgelist(df,
                                          source="gene_1",
                                          target="gene_2",
                                          edge_attr=edge_attributes)

        return G_genes

    # Numer of genes per chromosome
    chrom_card = df["chrom_1"].value_counts() + df["chrom_2"].value_counts()

    aggregation = {"weight": ["count", "mean"], "distance": "mean"}
    df["edge_id"] = df["chrom_1"] + "-" + df["chrom_2"]
    df["edge_id"] = df["edge_id"].apply(lambda x: "-".join(sorted(x.split("-"))
                                                           ))
    df_chrom = df.groupby("edge_id").agg(aggregation).reset_index()

    df_chrom["chrom_1"] = df_chrom["edge_id"].apply(lambda x: x.split("-")[0])
    df_chrom["chrom_2"] = df_chrom["edge_id"].apply(lambda x: x.split("-")[1])

    df_final = pd.DataFrame()
    df_final["chrom_1"] = df_chrom["chrom_1"]
    df_final["chrom_2"] = df_chrom["chrom_2"]
    df_final["strength"] = df_chrom["weight"]["mean"].values
    df_final["n_pairs"] = df_chrom["weight"]["count"].values
    df_final["mean_distance"] = df_chrom["distance"]["mean"].values

    G_chrom = nx.from_pandas_edgelist(
        df_final,
        source="chrom_1",
        target="chrom_2",
        edge_attr=["mean_distance", "strength", "n_pairs"])

    # Add node attributes
    node_attrs = {}
    for chrom in chrom_card.index:
        node_attrs[chrom] = chrom_card.loc[chrom]
    nx.set_node_attributes(G_chrom, node_attrs, name="size")

    return G_chrom
