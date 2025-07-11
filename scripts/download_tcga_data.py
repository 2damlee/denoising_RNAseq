# pipeline/scripts/download_tcga_data.py
import os, requests, subprocess, pandas as pd

def query_gdc_files(project, category, dtype, workflow_type=None):
    filters = {
        "op": "and",
        "content": [
            {"op": "in", "content": {"field": "cases.project.project_id", "value": [project]}},
            {"op": "in", "content": {"field": "files.data_category", "value": [category]}},
            {"op": "in", "content": {"field": "files.data_type", "value": [dtype]}},
            {"op": "in", "content": {"field": "files.access", "value": ["open"]}},
        ]
    }
    if workflow_type:
        filters["content"].append({
            "op": "in", "content": {"field": "analysis.workflow_type", "value": [workflow_type]}
        })

    params = {
        "filters": filters,
        "format": "JSON",
        "fields": "file_id",
        "size": 10000
    }
    r = requests.post("https://api.gdc.cancer.gov/files", json=params)
    r.raise_for_status()
    return [x["file_id"] for x in r.json()["data"]["hits"]]

def download_manifest(ids, filename):
    r = requests.post("https://api.gdc.cancer.gov/manifest", json={"ids": ids})
    with open(filename, "wb") as f:
        f.write(r.content)

def run_gdc_client(manifest, output_dir):
    subprocess.run(["gdc-client", "download", "-m", manifest, "-d", output_dir], check=True)
