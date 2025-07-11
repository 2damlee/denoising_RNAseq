# pipeline/scripts/annotate_gene_sequences.py
import gzip, pandas as pd, anndata as ad
from Bio import SeqIO

def extract_gene_sequences_extended(fasta_path):
    rows = []
    with gzip.open(fasta_path, "rt") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            parts = record.description.split("|")
            if len(parts) >= 6:
                transcript_id = parts[0].split(".")[0]
                gene_id = parts[1].split(".")[0]
                gene_symbol = parts[5]
                rows.append({
                    "gene_id": gene_id,
                    "gene_symbol_from_fasta": gene_symbol,
                    "transcript_id": transcript_id,
                    "sequence": str(record.seq)
                })
    return pd.DataFrame(rows).drop_duplicates("gene_id")

def annotate_sequences_to_adata(h5ad_path, fasta_path, output_path):
    adata = ad.read_h5ad(h5ad_path)
    seq_df = extract_gene_sequences_extended(fasta_path)
    if "gene_id_clean" not in adata.var.columns:
        adata.var["gene_id_clean"] = adata.var.index.str.replace(r"\.\d+$", "", regex=True)
    adata.var = adata.var.merge(seq_df, left_on="gene_id_clean", right_on="gene_id", how="left")
    adata.var["sequence"] = adata.var["sequence"].fillna("")
    adata.var.set_index("gene_id_clean", inplace=True)
    adata.write(output_path)
