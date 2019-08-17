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


def gene_network(edgelist,
                 biomart,
                 cutoff=(0, 1),
                 block="chrom",
                 normalize=True,
                 self_loops=True,
                 start="start"):
    """Create a block network and its edgelist based on a gene regulatory network.

    Parameters
    ----------
    edgelist : str
        Path to gene-gene weigthed edgelist. IMPORTANT: format of
        <tab>-separated table must be:
            gene_1, weight, gene_2.
    
    biomart: str
        Path to biomart gene information <tab>-sperated table. Table should
        contain a column `Gene stable ID` to identify each gene. A `block'
        column to identify the granularity of the network (see related keyword
        variable `block`).
    
    cutoff: 2-tuple of floats, default (0, 1)
        Half-closed interval $[a, b) \subset (0, 1)$ defining the edge set.
    
    block : str, default 'chrom'
        Name of gene label to define the gene partition and the resulting
        quotient network, where each group of genes defined by `block` is
        contracted into a single node. This label needs to be the name of a
        column in `biomart` (e.g. 'chromosome', 'karyoband', etc.).
    
    normalize : boolean, default True *TO INCLUDE*
       If True, normalizes block sizes.
    
    self_loops : boolean, default True *TO INCLUDE*
       If True, considers self loops (intra-block loops) and computes average
       distance between adjacent genes within blocks.
    
    start : str, default "start"
        Name of start column.
    """

    # Global variables
    start_1 = start + "_1"
    start_2 = start + "_2"
    block_1 = block + "_1"
    block_2 = block + "_2"

    # Read edgelist
    edgelist_header = ["gene_1", "weight", "gene_2"]
    edgelist_dtypes = ["str", "float", "str"]
    df = pd.read_csv(edgelist,
                     delimiter="\t",
                     names=edgelist_header,
                     dtype=dict(zip(edgelist_header, edgelist_dtypes)))
    # Cutoff filter
    df = df[(df["weight"] >= cutoff[0]) & (df["weight"] < cutoff[1])]

    # Read BioMart
    df_bio = pd.read_csv(biomart, delimiter="\t")
    df_bio = df_bio.drop_duplicates(subset="Gene stable ID")

    # Merge blocks
    df_merge = pd.merge(df,
                        df_bio,
                        how="left",
                        left_on="gene_1",
                        right_on="Gene stable ID")
    df_merge = df_merge.drop("Gene stable ID", axis=1)
    df_merge = pd.merge(df_merge,
                        df_bio,
                        how="left",
                        left_on="gene_2",
                        right_on="Gene stable ID",
                        suffixes=("_1", "_2"))
    df_merge = df_merge.drop("Gene stable ID", axis=1)

    # Cis/Trans Flag
    df_merge["type"] = np.where((df_merge[block_1] == df_merge[block_2]), "cis", "trans")
    # Distance
    df_merge["distance"] = 0
    df_merge.loc[(df_merge["type"] == "cis"), "distance"] = np.abs(
        df_merge.loc[(df_merge["type"] == "cis"), start_1] -
        df_merge.loc[(df_merge["type"] == "cis"), start_2])

    # Add variables to main DF
    df["distance"] = df_merge["distance"]
    df["type"] = df_merge["type"]

    # Compute number of genes per block, for normalization.
    block_card = df_merge[block_1].value_counts() + df_merge[block_2].value_counts()

    # Contracted edges need an aggregation function. Based on `weight`, we
    # compute `strength` and `n_pairs`. Being mean weight shared by gene-pairs
    # between blocks and number of gene-pairs, respectively. We also
    # consider `distance` to be the mean distance of connected gene-pairs within
    # a block.
    aggregation = {"weight": ["count", "mean", "median", "max", "min", "std"],
                   "distance": ["mean", "median"]}
    df["edge_id"] = df_merge[block_1] + "-" + df_merge[block_2]
    # Some genes might not be listed in `biomart` file.
    df = df.dropna().reset_index()
    df["edge_id"] = df["edge_id"].apply(lambda x: "-".join(sorted(x.split("-"))
                                                           ))
    df_block = df.groupby("edge_id").agg(aggregation).reset_index()

    df_block[block_1] = df_block["edge_id"].apply(lambda x: x.split("-")[0])
    df_block[block_2] = df_block["edge_id"].apply(lambda x: x.split("-")[1])

    df_final = pd.DataFrame()
    df_final[block_1] = df_block[block_1]
    df_final[block_2] = df_block[block_2]
    # Final values
    df_final["strength"] = df_block["weight"]["mean"].values
    df_final["median_weight"] = df_block["weight"]["median"].values
    df_final["max_weight"] = df_block["weight"]["max"].values
    df_final["min_weight"] = df_block["weight"]["min"].values
    df_final["std_weight"] = df_block["weight"]["std"].values
    df_final["n_pairs"] = df_block["weight"]["count"].values
    df_final["mean_distance"] = df_block["distance"]["mean"].values
    df_final["median_distance"] = df_block["distance"]["median"].values

    G_block = nx.from_pandas_edgelist(
        df_final,
        source=block_1,
        target=block_2,
        edge_attr=["mean_distance", "strength", "n_pairs"])

    # Add node attributes
    node_attrs = {}
    for _block in block_card.index:
        node_attrs[_block] = block_card.loc[_block]
    nx.set_node_attributes(G_block, node_attrs, name="size")

    return G_block, df_final
