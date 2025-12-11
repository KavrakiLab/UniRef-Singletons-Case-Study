import logging
from typing import Any
import pandas as pd
import requests
import dask.dataframe as dd
from builtins import tuple

from consistency_scoring import linclust_protocol_parallel
from lib.const import CliOPTIONS, benchmark_root_dir
from lib.funcs import parse_grep_output
from src.workflow import Workflow
import dask
from pathlib import Path
from typing import Tuple
dask.config.set({"dataframe.convert-string": True})


class UniProt2025Workflow(Workflow):

    def rescue_gaf(self, gaf_file, mmseqs_90, mmseqs_50, cd_hit_90):
        """
        load uniref100 cluster and see if there is a match

        :param cd_hit_90:
        :param mmseqs_50:
        :param mmseqs_90:
        :param gaf_file:
        :return:
        """
        uniref100_cluster_data = dd.read_parquet(
            self.sequence_db_100.parent.joinpath("uniref100_clusters.parquet").as_posix())

        na_rows = gaf_file[gaf_file.isna().any(axis=1)]
        merged = na_rows.merge(
            uniref100_cluster_data,
            left_on='DB_Object_ID',
            right_on='member',
            how='left'
        )
        merged = merged.merge(mmseqs_90, left_on='rep', right_on='1', how='left')
        merged = merged.drop('mmseqs_90', axis=1)
        merged = merged.rename(columns={"clusters": "mmseqs_90"})
        merged = merged.merge(mmseqs_50, left_on='rep', right_on='1', how='left')
        merged = merged.drop('mmseqs_50', axis=1)
        merged = merged.rename(columns={"clusters": "mmseqs_50"})
        merged = merged.merge(cd_hit_90, left_on='rep', right_on='members', how='left')
        merged = merged.drop('cdhit_90', axis=1)
        merged = merged.rename(columns={"clusters": "cd_hit_90"})
        merged = merged.dropna()
        merged = merged.rename(columns={"member": "DB_Object_ID"})
        merged = merged.set_index("DB_Object_ID")
        merged = merged.drop('rep_x', axis=1)
        merged = merged.drop('rep_y', axis=1)

        repaired_gaf = dd.concat([gaf_file.dropna(), merged])
        return repaired_gaf

    def build_gaf_dbs(self, gaf_all: Path, gaf_bf: Path, gaf_mf: Path, gaf_cc: Path) -> Tuple[Any, Any, Any, Any]:
        """Build or load GAF databases for CD-HIT, MMseqs2, and Linclust."""
        logging.info("Building GAF databases for CD-HIT and LINCLUST.")

        # === Load from disk if already built ===
        if gaf_all.exists():
            logging.info(f"Loading precomputed GAF database from {gaf_all}")
            return dd.read_parquet(gaf_all.as_posix()).compute(), None, None, None

        # === Step 1: Download and process GAF database ===
        gaf_db = self.download_gaf_db()
        gaf_all_df = self.process_gaf_file(gaf_db)

        # === Step 2: Run CD-HIT ===
        cd_hit_90 = self.run_cd_hit_command(self.sequence_db_100)
        cd_hit_processed_path = cd_hit_90.parent / "cd_hit_processed.parquet"

        if not cd_hit_processed_path.exists():
            logging.info("Post-processing CD-HIT clusters...")
            cdhit_uniref90_df = Workflow.post_process_cd_hit_clustering(cd_hit_90)
            cdhit_uniref90_df.to_parquet(cd_hit_processed_path)
        else:
            logging.info("Loading cached CD-HIT clustering results.")
            cdhit_uniref90_df = dd.read_parquet(cd_hit_processed_path.as_posix())

        # === Step 3: Process Linclust outputs ===
        linclust_cluster_90 = parse_grep_output(
            self.sequence_db_90.parent / "parsed_xml_grep_output.txt",
            self.sequence_db_90.parent / "uniref90_linclust.parquet"
        )
        linclust_cluster_50 = parse_grep_output(
            self.sequence_db_50.parent / "parsed_xml_grep_output.txt",
            self.sequence_db_50.parent / "uniref50_linclust.parquet"
        )

        linclust_cluster_90_df = self.post_process_cdhit_xml_dataframe(linclust_cluster_90)
        linclust_cluster_50_df = self.post_process_cdhit_xml_dataframe(linclust_cluster_50)

        # === Step 4: Run MMseqs2 ===
        mmseqs2_cluster_90, mmseqs2_cluster_50 = self.run_mmseqs_commands(self.sequence_db_100)
        mmseqs2_90_df = self.post_process_mmseqs2_clustering(mmseqs2_cluster_90)
        mmseqs2_50_df = self.post_process_mmseqs2_clustering(mmseqs2_cluster_50)

        # === Step 5: Merge clustering results into each GAF db ===
        logging.info("Merging clustering results into GAF ALL databases.")

        gaf_all_df["cdhit_90"] = cdhit_uniref90_df["rep"] # need to fix dype to pyarrow string by changing index to members
        gaf_all_df["mmseqs_90"] = mmseqs2_90_df["clusters"]
        gaf_all_df["mmseqs_50"] = mmseqs2_50_df["clusters"]
        gaf_all_df["linclust_90"] = linclust_cluster_90_df["clusters"]
        gaf_all_df["linclust_50"] = linclust_cluster_50_df["clusters"]
        gaf_all_df.to_parquet("/home/Users/fmq1/intermediate_clustering.parquet")
        # === Step 6: Rescue and save ===
        logging.info("Rescuing and saving merged GAF database.")
        gaf_all_df = self.rescue_gaf(gaf_all_df, mmseqs2_90_df, mmseqs2_50_df, cdhit_uniref90_df)
        gaf_all_df.to_parquet(gaf_all.as_posix())

        return gaf_all_df, None, None, None

    def run(self):
        gaf_all = benchmark_root_dir.joinpath(self.benchmark_working_dir).joinpath('gaf_2025_with_clustering_data_all')
        gaf_bf = benchmark_root_dir.joinpath(self.benchmark_working_dir).joinpath('gaf_2025_with_clustering_data_bf')
        gaf_mf = benchmark_root_dir.joinpath(self.benchmark_working_dir).joinpath('gaf_2025_with_clustering_data_mf')
        gaf_cc = benchmark_root_dir.joinpath(self.benchmark_working_dir).joinpath('gaf_2025_with_clustering_data_cc')
        # gaf_all_db, gaf_cc_db, gaf_bf_db, gaf_mf_db = self.build_gaf_dbs(gaf_all, gaf_bf, gaf_mf, gaf_cc)
        gaf_all_db = self.build_gaf_dbs(gaf_all, gaf_bf, gaf_mf, gaf_cc)

        gaf_dbs = {"all": gaf_all_db[0]}  # , "cc": gaf_cc_db, "bf": gaf_bf_db, "mf": gaf_mf_db}
        for db_type, db in gaf_dbs.items():
            logging.info(f"Running for {db_type}")
            cluster_dicts = {'cd_hit_90': {},
                             'linclust_90': {}, 'linclust_50': {},
                             'mmseqs_90': {}, 'mmseqs_50': {}}
            datasets = []
            cluster_results_not_filtered = self.add_go_terms_to_dataframe(db, cluster_dicts, False)
            for cluster_algo in cluster_dicts.keys():
                datasets.append((cluster_results_not_filtered[cluster_algo], "non_filtered_" + cluster_algo + "_results_gogo.csv"))
            for data, filename in datasets:
                self.process_and_save(data, filename)
#            datasets=[]
#            cluster_results = self.add_go_terms_to_dataframe(db, cluster_dicts, True)
#            for cluster_algo in cluster_dicts.keys():
#                datasets.append((cluster_results[cluster_algo], "filtered_" + cluster_algo + "_results_gogo.csv"))
            for data, filename in datasets:
                self.process_and_save(data, filename)

    def process_and_save(self, data, filename):
        results = linclust_protocol_parallel(data, self.go_path)
        results = Workflow.post_process_consistency_data(results)
        output_path = benchmark_root_dir.joinpath(self.benchmark_working_dir, filename)
        results.to_csv(output_path)

    @property
    def benchmark_year(self) -> str:
        return CliOPTIONS.Benchmark2025.value

    def download_sequence_db(self) -> tuple:
        """Download UniRef100/90/50 FASTA files if they do not already exist."""
        base_url = "https://ftp.uniprot.org/pub/databases/uniprot/uniref"
        uniref_levels = ["uniref100", "uniref90", "uniref50"]

        paths = []

        for level in uniref_levels:
            # Define paths
            dir_path = benchmark_root_dir / self.benchmark_working_dir / level
            dir_path.mkdir(parents=True, exist_ok=True)
            fasta_path = dir_path / f"{level}.fasta"
            paths.append(fasta_path)

            # Download if not present
            if not fasta_path.exists():
                url = f"{base_url}/{level}/{level}.fasta.gz"
                logging.info(f"Downloading {level}...")
                response = requests.get(url, stream=True)
                response.raise_for_status()
                with open(fasta_path, "wb") as f:
                    for chunk in response.iter_content(chunk_size=8192):
                        if chunk:
                            f.write(chunk)
                logging.info(f"{level} downloaded to {fasta_path}")
        return tuple(paths)

    def download_gaf_db(self) -> Path:
        goa_uniprot_all = "http://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/goa_uniprot_all.gaf.gz"
        output_path = benchmark_root_dir.joinpath(self.benchmark_working_dir).joinpath("goa_file_2025.gz")
        if not output_path.exists():
            with requests.get(goa_uniprot_all, stream=True) as r:
                r.raise_for_status()
                with open(output_path, 'wb') as f:
                    for chunk in r.iter_content(chunk_size=8192):
                        if chunk:  # filter out keep-alive chunks
                            f.write(chunk)
        return output_path
