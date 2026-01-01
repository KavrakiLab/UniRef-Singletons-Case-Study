import gzip
import logging
import os
import shutil
import tarfile
from pathlib import Path
from typing import Tuple, Any

import requests

from .consistency_scoring import linclust_protocol_parallel
from lib.const import CliOPTIONS, benchmark_root_dir
from funcs import parse_uniref_xmls
from src.workflow import Workflow
import dask.dataframe as dd


class UniProt2016Workflow(Workflow):

    @property
    def benchmark_year(self) -> str:
        return CliOPTIONS.Benchmark2016.value



    def download_sequence_db(self) -> tuple[Path, Path, Path]:
        base_url = "https://ftp.uniprot.org/pub/databases/uniprot/previous_major_releases/release-2016_03/uniref/uniref2016_03.tar.gz"
        work_dir = benchmark_root_dir / self.benchmark_working_dir
        output_path = work_dir / "uniref100_2016.tar.gz"
        work_dir.mkdir(parents=True, exist_ok=True)

        # Helper to extract tar files
        def extract_tar(tar_path: Path, dest_dir: Path):
            if not dest_dir.exists():
                dest_dir.mkdir(parents=True, exist_ok=True)
                with tarfile.open(tar_path, "r:*") as tar:
                    tar.extractall(path=dest_dir)

        # Helper to decompress gzip
        def decompress_gzip(gz_path: Path, out_path: Path):
            if not out_path.exists():
                with gzip.open(gz_path, "rb") as f_in, open(out_path, "wb") as f_out:
                    shutil.copyfileobj(f_in, f_out)

        # 1. Download if needed
        if not output_path.exists():
            with requests.get(base_url, stream=True) as r:
                r.raise_for_status()
                with open(output_path, "wb") as f:
                    for chunk in r.iter_content(chunk_size=8192):
                        if chunk:
                            f.write(chunk)

        # 2. Extract main tar.gz
        main_extract_dir = work_dir
        main_tar_path = work_dir / "uniref100.tar"
        if not main_tar_path.exists():
            with tarfile.open(output_path, "r:gz") as tar:
                tar.extractall(path=main_extract_dir)

        # 3. Process uniref100 / uniref90 / uniref50
        result_paths = {}
        for db in ["uniref100", "uniref90", "uniref50"]:
            tar_file = work_dir / f"{db}.tar"
            extract_dir = work_dir / db
            xml_gz = extract_dir / f"{db}.xml.gz"
            xml_out = extract_dir / f"{db}.xml"

            extract_tar(tar_file, extract_dir)
            decompress_gzip(xml_gz, xml_out)
            result_paths[db] = xml_out

        return result_paths["uniref100"], result_paths["uniref90"], result_paths["uniref50"]

    def download_gaf_db(self) -> Path:
        goa_uniprot_all = "https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/UNIPROT/goa_uniprot_all.gaf.156.gz"
        logging.info(f"Downloading GAF DB witht he following URL: {goa_uniprot_all}")
        output_path = benchmark_root_dir.joinpath(self.benchmark_working_dir).joinpath("../goa_file_2016.gz")
        if not output_path.exists():
            with requests.get(goa_uniprot_all, stream=True) as r:
                r.raise_for_status()
                with open(output_path, 'wb') as f:
                    for chunk in r.iter_content(chunk_size=8192):
                        if chunk:  # filter out keep-alive chunks
                            f.write(chunk)
            logging.debug("Done writing")
        return output_path

    def build_gaf_dbs(self, gaf_all, gaf_bf, gaf_mf, gaf_cc) -> Tuple[Any, Any, Any, Any]:
        logging.info(f"Building GAF databases for CD-HIT and LINCLUST.")
        if not gaf_all.exists() or not gaf_bf.exists() or not gaf_mf.exists() or not gaf_cc.exists():
            gaf_db = self.download_gaf_db()
            gaf_dbs = self.process_gaf_file(gaf_db,individual_components=True)
            cdhit_uniref90_df = parse_uniref_xmls(self.sequence_db_90.as_posix(),
                                                  self.sequence_db_90.parent.joinpath('uniref90_cd_hit.parquet').as_posix())
            cdhit_uniref50_df = parse_uniref_xmls(self.sequence_db_50.as_posix(),
                                                  self.sequence_db_50.parent.joinpath('uniref50_cd_hit.parquet').as_posix())
            cdhit_uniref90_df = self.post_process_cdhit_xml_dataframe(cdhit_uniref90_df)
            cdhit_uniref50_df = self.post_process_cdhit_xml_dataframe(cdhit_uniref50_df)

            mmseqs_cluster_90, mmseqs_cluster_50 = self.run_mmseqs_commands(self.sequence_db_100)
            mmseqs2_90_results: dd.DataFrame = self.post_process_mmseqs2_clustering(mmseqs_cluster_90)
            mmseqs2_50_results: dd.DataFrame = self.post_process_mmseqs2_clustering(mmseqs_cluster_50)
            linclust_cluster_90, linclust_cluster_50 = self.run_linclust_commands(self.sequence_db_100)
            linclust_90_results: dd.DataFrame = self.post_process_mmseqs2_clustering(linclust_cluster_90)
            linclust_50_results: dd.DataFrame = self.post_process_mmseqs2_clustering(linclust_cluster_50)

            logging.info("Merging clustering results into GAF ALL databases.")

            for db in gaf_dbs:
                db["cdhit_90"] = cdhit_uniref90_df["clusters"]
                db["cdhit_50"] = cdhit_uniref50_df["clusters"]
                db["mmseqs_90"] = mmseqs2_90_results["clusters"]
                db["mmseqs_50"] = mmseqs2_50_results["clusters"]
                db["linclust_90"] = linclust_90_results["clusters"]
                db["linclust_50"] = linclust_50_results["clusters"]
            gaf_dbs[0].to_parquet(gaf_all.as_posix())
            gaf_dbs[1].to_parquet(gaf_cc.as_posix())
            gaf_dbs[2].to_parquet(gaf_bf.as_posix())
            gaf_dbs[3].to_parquet(gaf_mf.as_posix())
            return gaf_dbs[0], gaf_dbs[1], gaf_dbs[2], gaf_dbs[3]
        else:
            return (dd.read_parquet(gaf_all.as_posix()).compute(), dd.read_parquet(gaf_cc.as_posix()).compute(),
                    dd.read_parquet(gaf_bf.as_posix()).compute(), dd.read_parquet(gaf_mf.as_posix()).compute())

    def run(self):
        gaf_all = benchmark_root_dir.joinpath(self.benchmark_working_dir).joinpath('gaf_2025_with_clustering_data_all')
        gaf_bf = benchmark_root_dir.joinpath(self.benchmark_working_dir).joinpath('gaf_2025_with_clustering_data_bf')
        gaf_mf = benchmark_root_dir.joinpath(self.benchmark_working_dir).joinpath('gaf_2025_with_clustering_data_mf')
        gaf_cc = benchmark_root_dir.joinpath(self.benchmark_working_dir).joinpath('gaf_2025_with_clustering_data_cc')
        gaf_all_db, gaf_cc_db, gaf_bf_db, gaf_mf_db = self.build_gaf_dbs(gaf_all, gaf_bf, gaf_mf, gaf_cc)
        gaf_dbs = {"all": gaf_all_db, "cc": gaf_cc_db, "bf": gaf_bf_db, "mf": gaf_mf_db}
        for db_type, db in gaf_dbs.items():
            logging.info(f"Running for {db_type}")
            cluster_dicts = {'cdhit_90': {},'cdhit_50': {},
                             'linclust_90': {}, 'linclust_50': {},
                             'mmseqs_90': {}, 'mmseqs_50': {}}
            datasets = []
            cluster_results_not_filtered = self.add_go_terms_to_dataframe(db, cluster_dicts, False)
            cluster_results = self.add_go_terms_to_dataframe(db, cluster_dicts, True)
            for cluster_algo in cluster_dicts.keys():
                datasets.append(cluster_results_not_filtered[cluster_algo], "non_filtered_" + cluster_algo + "_results.csv")
                datasets.append(cluster_results[cluster_algo], "filtered_" + cluster_algo + "_results.csv")
            for data, filename in datasets:
                self.process_and_save(data, filename)

    def process_and_save(self, data, filename):
        results = linclust_protocol_parallel(data, self.go_path)
        results = Workflow.post_process_consistency_data(results)
        output_path = benchmark_root_dir.joinpath(self.benchmark_working_dir, filename)
        results.to_csv(output_path)
