import argparse
import scanpy as sc
import pandas as pd
import numpy as np
import scipy.sparse as sp
from Bio import SeqIO

def load_pbmc_10x(input_h5):
    print(f"[INFO] Reading 10X h5 file from: {input_h5}")
    adata = sc.read_10x_h5(input_h5)
    adata.var["gene_id_clean"] = adata.var_names.str.split(".").str[0]
    return adata

def map_sequences(adata, fasta_path):
    print(f"[INFO] Mapping gene_id → sequence from FASTA: {fasta_path}")
    gene_seq_map = {}
    for record in SeqIO.parse(fasta_path, "fasta"):
        header_parts = record.description.split("|")
        for part in header_parts:
            if part.startswith("ENSG"):
                gene_id = part.split(".")[0]
                gene_seq_map[gene_id] = str(record.seq)
                break
    adata.var["sequence"] = adata.var["gene_id_clean"].map(gene_seq_map)
    mapped = adata.var["sequence"].notnull().sum()
    missing = adata.var["sequence"].isnull().sum()
    print(f"[INFO] Successfully mapped: {mapped} / {adata.var.shape[0]}")
    print(f"[INFO] Missing sequences: {missing}")
    return adata

def add_obs_metadata(adata):
    print("[INFO] Generating .obs metadata...")
    if adata.obs is None or adata.obs.empty:
        adata.obs = pd.DataFrame(index=adata.obs_names)
    adata.obs["source"] = "PBMC10k"
    adata.obs["batch"] = "PBMC10k"
    if "gene_ids" in adata.var.columns:
        mito_genes = adata.var_names.str.upper().str.startswith("MT-")
    else:
        mito_genes = adata.var.index.str.upper().str.startswith("MT-")
    adata.obs["mt_frac"] = adata[:, mito_genes].X.sum(axis=1).A1 / adata.X.sum(axis=1).A1
    adata.obs["n_counts"] = adata.X.sum(axis=1).A1
    adata.obs["n_genes"] = (adata.X > 0).sum(axis=1).A1
    adata.obs["cell_id"] = adata.obs_names
    print("[INFO] .obs columns:", adata.obs.columns.tolist())
    return adata

def inspect_matrix_summary(adata):
    print("\n>>> .X (expression matrix):")
    print("Shape:", adata.shape)
    print("Sparse Type:", type(adata.X))
    print("Non-zero Elements:", adata.X.nnz if sp.issparse(adata.X) else np.count_nonzero(adata.X))
    mean_expr = adata.X[:, :5].mean(axis=0).A1 if sp.issparse(adata.X) else adata.X[:, :5].mean(axis=0)
    print("Mean expression of first 5 genes:", mean_expr)

    print("\n>>> .obs preview:")
    print(adata.obs.head())

    print("\n>>> .var preview:")
    print(adata.var[["gene_id_clean", "sequence"]].head())
    if "sequence" in adata.var.columns:
        missing = adata.var["sequence"].isnull().sum() + (adata.var["sequence"] == "").sum()
        print(f"Missing sequences in `.var['sequence']`: {missing} / {adata.var.shape[0]}")

def main(input_h5, fasta_path, output_path):
    adata = load_pbmc_10x(input_h5)
    adata = map_sequences(adata, fasta_path)
    adata = add_obs_metadata(adata)
    inspect_matrix_summary(adata)
    print(f"[INFO] Writing final AnnData to: {output_path}")
    adata.write(output_path)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="PBMC preprocessing with sequence mapping and QC metrics.")
    parser.add_argument('--input_h5', required=True, help='Path to 10X raw_feature_bc_matrix.h5 file')
    parser.add_argument('--fasta', required=True, help='Path to GENCODE transcript FASTA file')
    parser.add_argument('--output', required=True, help='Output path for processed .h5ad file')
    args = parser.parse_args()
    main(args.input_h5, args.fasta, args.output)
