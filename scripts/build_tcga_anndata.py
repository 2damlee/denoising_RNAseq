# pipeline/scripts/build_tcga_anndata.py
import os, pandas as pd, numpy as np
from scipy import sparse
import anndata as ad
from lxml import etree
from glob import glob

def fetch_mapping_from_manifest(manifest_path, output_csv):
    with open(manifest_path) as f:
        lines = f.readlines()[1:]
        uuids = [line.split("\t")[0] for line in lines if line.strip()]
    url = "https://api.gdc.cancer.gov/files"
    all_records = []
    for i in range(0, len(uuids), 1000):
        batch = uuids[i:i + 1000]
        filters = {"op": "in", "content": {"field": "file_id", "value": batch}}
        params = {
            "filters": filters,
            "format": "JSON",
            "fields": "file_id,cases.submitter_id",
            "size": len(batch)
        }
        r = requests.post(url, json=params)
        r.raise_for_status()
        for h in r.json()["data"]["hits"]:
            if h["cases"]:
                all_records.append((h["file_id"], h["cases"][0]["submitter_id"]))
    df = pd.DataFrame(all_records, columns=["uuid", "barcode"]).drop_duplicates("uuid")
    df.to_csv(output_csv, index=False)
    return df.set_index("uuid")["barcode"].to_dict()

def merge_counts_with_mapping(count_dir, mapping, output_csv):
    dfs = []
    for subdir, _, files in os.walk(count_dir):
        uuid = os.path.basename(subdir)
        if uuid not in mapping:
            continue
        sample_id = "-".join(mapping[uuid].split("-")[:3])
        for f in files:
            if not f.endswith(".tsv"):
                continue
            df = pd.read_csv(os.path.join(subdir, f), sep="\t", comment="#")
            if "gene_id" in df.columns and "unstranded" in df.columns:
                df = df[["gene_id", "unstranded"]].set_index("gene_id")
                df.columns = [sample_id]
                dfs.append(df)
    merged = pd.concat(dfs, axis=1)
    merged.to_csv(output_csv)
    return merged

def parse_clinical_xml(xml_dir, output_csv):
    records = []
    for f in glob(f"{xml_dir}/**/*.xml", recursive=True):
        try:
            filename = os.path.basename(f)
            sample_id = "-".join(filename.split(".")[-2].split("-")[:3])
            tree = etree.parse(f)
            text_dict = {e.tag.split("}")[-1]: e.text for e in tree.iter()}
            gender = text_dict.get("gender")
            age = text_dict.get("age_at_diagnosis") or text_dict.get("age_at_initial_pathologic_diagnosis")
            vital_status = text_dict.get("vital_status")
            tumor_stage = text_dict.get("clinical_stage") or text_dict.get("pathologic_stage")
            batch = sample_id.split("-")[1]
            records.append({
                "sample_id": sample_id, "gender": gender, "age": age,
                "vital_status": vital_status, "tumor_stage": tumor_stage, "batch": batch
            })
        except Exception:
            continue
    df = pd.DataFrame(records).drop_duplicates("sample_id").set_index("sample_id")
    df.to_csv(output_csv)
    return df

def match_and_save(count_csv, clin_csv, gene_mapping_csv, output_h5ad):
    counts = pd.read_csv(count_csv, index_col=0)
    metadata = pd.read_csv(clin_csv, index_col=0)
    gene_map = pd.read_csv(gene_mapping_csv)
    counts.columns = ["-".join(c.split("-")[:3]) for c in counts.columns]
    metadata.index = ["-".join(i.split("-")[:3]) for i in metadata.index]
    common = sorted(set(counts.columns) & set(metadata.index))
    counts = counts[common]
    metadata = metadata.loc[common]
    gene_map = gene_map.drop_duplicates("ensembl_gene_id").set_index("ensembl_gene_id")
    var = pd.DataFrame(index=counts.index)
    var["gene_id_clean"] = counts.index.str.replace(r"\.\d+$", "", regex=True)
    var["gene_symbol"] = var["gene_id_clean"].map(gene_map["hgnc_symbol"])
    X = sparse.csr_matrix(counts.T.fillna(0).astype(np.float32))
    adata = ad.AnnData(X=X, obs=metadata, var=var)
    adata.write(output_h5ad)
    return adata
