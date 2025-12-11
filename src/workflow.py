import abc
import logging
import os
import re
import shutil
from pathlib import Path
from typing import Optional, Tuple
import plotly.graph_objects as go
import plotly.express as px
import pandas as pd
import numpy as np
import subprocess
import tqdm
from goatools.base import get_godag

from lib.const import CliOPTIONS, benchmark_working_dir, benchmark_root_dir, EXP_EVIDENCE_CODES
from lib.funcs import xml_to_fasta, parse_uniref_xmls
import dask.dataframe as dd
import tarfile
import gzip


class Workflow(abc.ABC):

    def __init__(self, mmseqs2_executable_path: str, cd_hit_path: str, gene_ontology_path: str):
        self.go_path = gene_ontology_path
        self._mmseqs2_executable_path: Path = Path(mmseqs2_executable_path)
        self._cd_hit_path: Path = Path(cd_hit_path)
        self._mmseqs2_results: Optional[pd.DataFrame] = None
        self._cd_hit_results: Optional[pd.DataFrame] = None
        self._prepare_files()

    @property
    @abc.abstractmethod
    def benchmark_year(self) -> str:
        raise NotImplementedError

    def _prepare_files(self):
        """
        Load GAF annotation file, grab uniref100 for subsequent clustering. Build uniref100 fasta file if needed.
        :return:
        """
        logging.info("Preparing files before running...")
        sequence_db_100, sequence_db_90, sequence_db_50 = self.download_sequence_db()
        if self.benchmark_year == CliOPTIONS.Benchmark2016.value:
            logging.info("Running 2016 workflow")
            #sequence_db_100 = sequence_db_100.with_suffix(".fasta")
            if not sequence_db_100.exists():
                xml_to_fasta(sequence_db_100.as_posix(), sequence_db_100.with_suffix(".fasta").as_posix())
        gaf_file = self.download_gaf_db()
        gaf_parquet = self.process_gaf_file(gaf_file)
        self.gaf_parquet = gaf_parquet
        self.sequence_db_100: Path = sequence_db_100
        self.sequence_db_90: Path = sequence_db_90
        self.sequence_db_50: Path = sequence_db_50



    @property
    def benchmark_working_dir(self) -> str:
        return benchmark_working_dir[self.benchmark_year]

    @abc.abstractmethod
    def download_sequence_db(self) -> Path:
        raise NotImplementedError

    @abc.abstractmethod
    def download_gaf_db(self) -> Path:
        raise NotImplementedError

    @abc.abstractmethod
    def run(self):
        raise NotImplementedError

    def run_mmseqs_commands(self, sequence_db):
        logging.info("Running mmseqs commands")
        if not benchmark_root_dir.joinpath(self.benchmark_working_dir).joinpath('uniref100_db').exists():
            build_db: [str] = [self._mmseqs2_executable_path, 'createdb',
                               benchmark_root_dir.joinpath(self.benchmark_working_dir).joinpath(
                                   sequence_db).as_posix(),
                               benchmark_root_dir.joinpath(self.benchmark_working_dir).joinpath(
                                   'uniref100_db').as_posix()]
            logging.info(f"Running createdb command {build_db}")
            result_build_db = subprocess.run(build_db, capture_output=False, text=True, check=True)
        if not benchmark_root_dir.joinpath(self.benchmark_working_dir).joinpath('90_mmseqs').joinpath(
                '90_identity.index').exists():
            run_command_90 = [self._mmseqs2_executable_path, 'cluster',
                              benchmark_root_dir.joinpath(self.benchmark_working_dir).joinpath(
                                  'uniref100_db').as_posix(),
                              benchmark_root_dir.joinpath(self.benchmark_working_dir).joinpath('90_mmseqs').joinpath(
                                  '90_identity'),
                              "/home/Users/fmq1/tmp",
                              '--cov-mode', '1', '-c', '0.9', '--min-seq-id', '0.9', '--cluster-mode', '2']
            logging.info(f"Running mmseqs command {run_command_90}")
            result_command_90 = subprocess.run(run_command_90, capture_output=False, text=True, check=True)
        if not benchmark_root_dir.joinpath(self.benchmark_working_dir).joinpath('90_mmseqs').joinpath(
                'clusters.tsv').exists():
            create_tsv_command_90 = [self._mmseqs2_executable_path, 'createtsv',
                                     benchmark_root_dir.joinpath(self.benchmark_working_dir).joinpath(
                                         'uniref100_db').as_posix(),
                                     benchmark_root_dir.joinpath(self.benchmark_working_dir).joinpath(
                                         'uniref100_db').as_posix(),
                                     benchmark_root_dir.joinpath(self.benchmark_working_dir).joinpath(
                                         '90_mmseqs').joinpath(
                                         '90_identity'),
                                     benchmark_root_dir.joinpath(self.benchmark_working_dir).joinpath(
                                         '90_mmseqs').joinpath(
                                         'clusters.tsv'),
                                     ]
            logging.info(f"Running createtsv command {create_tsv_command_90}")
            create_tsv_command_90 = subprocess.run(create_tsv_command_90, capture_output=False, text=True, check=True)
        if not benchmark_root_dir.joinpath(self.benchmark_working_dir).joinpath('50_mmseqs').joinpath(
                '50_identity.index').exists():
            run_command_50 = [self._mmseqs2_executable_path, 'cluster',
                              benchmark_root_dir.joinpath(self.benchmark_working_dir).joinpath(
                                  'uniref100_db').as_posix(),
                              benchmark_root_dir.joinpath(self.benchmark_working_dir).joinpath('50_mmseqs').joinpath(
                                  '50_identity'),
                              "/home/Users/fmq1/tmp",
                              '--cov-mode', '1', '-c', '0.9', '--min-seq-id', '0.5', '--cluster-mode', '2']
            logging.info(f"Running mmseqs command {run_command_50}")
            result_command_50 = subprocess.run(run_command_50, capture_output=False, text=True, check=True)
        if not benchmark_root_dir.joinpath(self.benchmark_working_dir).joinpath('50_mmseqs').joinpath(
                'clusters.tsv').exists():
            create_tsv_command_50 = [self._mmseqs2_executable_path, 'createtsv',
                                     benchmark_root_dir.joinpath(self.benchmark_working_dir).joinpath(
                                         'uniref100_db').as_posix(),
                                     benchmark_root_dir.joinpath(self.benchmark_working_dir).joinpath(
                                         'uniref100_db').as_posix(),
                                     benchmark_root_dir.joinpath(self.benchmark_working_dir).joinpath(
                                         '50_mmseqs').joinpath(
                                         '50_identity_cluster'),
                                     benchmark_root_dir.joinpath(self.benchmark_working_dir).joinpath(
                                         '50_mmseqs').joinpath(
                                         'clusters.tsv'),
                                     ]
            logging.info(f"Running createtsv command {create_tsv_command_50}")
            create_tsv_command_50 = subprocess.run(create_tsv_command_50, capture_output=False, text=True, check=True)
        return (benchmark_root_dir.joinpath(self.benchmark_working_dir).joinpath('90_mmseqs').joinpath('clusters.tsv'),
                benchmark_root_dir.joinpath(self.benchmark_working_dir).joinpath('50_mmseqs').joinpath('clusters.tsv'))

    def run_linclust_commands(self, sequence_db):
        logging.info("Running linclust commands")
        if not benchmark_root_dir.joinpath(self.benchmark_working_dir).joinpath('uniref100_db').exists():
            build_db: [str] = [self._mmseqs2_executable_path, 'createdb',
                               benchmark_root_dir.joinpath(self.benchmark_working_dir).joinpath(
                                   sequence_db).as_posix(),
                               benchmark_root_dir.joinpath(self.benchmark_working_dir).joinpath(
                                   'uniref100_db').as_posix()]
            logging.info(f"Running createdb command {build_db}")
            result_build_db = subprocess.run(build_db, capture_output=False, text=True, check=True)
        if not benchmark_root_dir.joinpath(self.benchmark_working_dir).joinpath('90').joinpath(
                '90_identity_cluster.index').exists():
            run_command_90 = [self._mmseqs2_executable_path, 'linclust',
                              benchmark_root_dir.joinpath(self.benchmark_working_dir).joinpath(
                                  'uniref100_db').as_posix(),
                              benchmark_root_dir.joinpath(self.benchmark_working_dir).joinpath('90').joinpath(
                                  '90_identity_cluster'),
                              benchmark_root_dir.joinpath(self.benchmark_working_dir).joinpath('tmp'),
                              '--cov-mode', '1', '-c', '0.9', '--min-seq-id', '0.9', '--cluster-mode', '2']
            logging.info(f"Running linclust command {run_command_90}")
            result_command_90 = subprocess.run(run_command_90, capture_output=False, text=True, check=True)
        if not benchmark_root_dir.joinpath(self.benchmark_working_dir).joinpath('90').joinpath(
                'clusters.tsv').exists():
            create_tsv_command_90 = [self._mmseqs2_executable_path, 'createtsv',
                                     benchmark_root_dir.joinpath(self.benchmark_working_dir).joinpath(
                                         'uniref100_db').as_posix(),
                                     benchmark_root_dir.joinpath(self.benchmark_working_dir).joinpath(
                                         'uniref100_db').as_posix(),
                                     benchmark_root_dir.joinpath(self.benchmark_working_dir).joinpath('90').joinpath(
                                         '90_identity_cluster'),
                                     benchmark_root_dir.joinpath(self.benchmark_working_dir).joinpath('90').joinpath(
                                         'clusters.tsv'),
                                     ]
            logging.info(f"Running createtsv command {create_tsv_command_90}")
            create_tsv_command_90 = subprocess.run(create_tsv_command_90, capture_output=False, text=True, check=True)
        if not benchmark_root_dir.joinpath(self.benchmark_working_dir).joinpath('50').joinpath(
                '50_identity_cluster.index').exists():
            run_command_50 = [self._mmseqs2_executable_path, 'linclust',
                              benchmark_root_dir.joinpath(self.benchmark_working_dir).joinpath(
                                  'uniref100_db').as_posix(),
                              benchmark_root_dir.joinpath(self.benchmark_working_dir).joinpath('50').joinpath(
                                  '50_identity_cluster'),
                              "/home/Users/fmq1/tmp",
                              '--cov-mode', '1', '-c', '0.9', '--min-seq-id', '0.5', '--cluster-mode', '2']
            logging.info(f"Running linclust command {run_command_50}")
            result_command_50 = subprocess.run(run_command_50, capture_output=False, text=True, check=True)
        if not benchmark_root_dir.joinpath(self.benchmark_working_dir).joinpath('50').joinpath(
                'clusters.tsv').exists():
            create_tsv_command_50 = [self._mmseqs2_executable_path, '',
                                     benchmark_root_dir.joinpath(self.benchmark_working_dir).joinpath(
                                         'uniref100_db').as_posix(),
                                     benchmark_root_dir.joinpath(self.benchmark_working_dir).joinpath(
                                         'uniref100_db').as_posix(),
                                     benchmark_root_dir.joinpath(self.benchmark_working_dir).joinpath(
                                         '50').joinpath(
                                         '50_identity_cluster'),
                                     benchmark_root_dir.joinpath(self.benchmark_working_dir).joinpath(
                                         '50').joinpath(
                                         'clusters.tsv'),
                                     ]
            logging.info(f"Running createtsv command {create_tsv_command_50}")
            create_tsv_command_50 = subprocess.run(create_tsv_command_50, capture_output=False, text=True, check=True)
        return (benchmark_root_dir.joinpath(self.benchmark_working_dir).joinpath('90').joinpath('clusters.tsv'),
                benchmark_root_dir.joinpath(self.benchmark_working_dir).joinpath('50').joinpath('clusters.tsv'))

    def run_cd_hit_command(self, sequence_db) -> Path:
        logging.info("Running cd_hit clustering")
        if not benchmark_root_dir.joinpath(self.benchmark_working_dir).joinpath('90').joinpath(
                'cd_hit_90.clstr').exists():
            command_90: [str] = [self._cd_hit_path, '-n', '5', '-M', '0', '-c', '0.9', '-i',
                                 benchmark_root_dir.joinpath(self.benchmark_working_dir).joinpath(
                                     sequence_db).as_posix(),
                                 '-o', benchmark_root_dir.joinpath(self.benchmark_working_dir).joinpath('90').joinpath(
                    'cd_hit_90')]
            result_command_90 = subprocess.run(command_90, capture_output=True, text=True, check=True)
        return benchmark_root_dir.joinpath(self.benchmark_working_dir).joinpath('90').joinpath('cd_hit_90.clstr')

    @staticmethod
    def post_process_cd_hit_clustering(cd_hit_clustering_output_path: Path):
        clusters = []
        current_cluster = None
        with open(cd_hit_clustering_output_path, 'r') as file:
            for line in tqdm.tqdm(file):
                line = line.strip()
                if line.startswith(">Cluster"):
                    current_cluster = int(line.split()[1])
                else:
                    match = re.match(r"(\d+)\s+(\d+)aa, >(.*?)\.\.\. ?(?:at ([\d.]+)%|(\*))", line)
                    if match:
                        accession = match.group(3)
                        representative = match.group(5) == '*'
                        clusters.append({
                            "clusters": current_cluster,
                            "accession": accession,
                            "representative": representative})
        raw_df_translation = pd.DataFrame(clusters)
        rep_map = raw_df_translation.loc[raw_df_translation["representative"],
        ["clusters", "accession"]].set_index("clusters")["accession"]
        raw_df_translation["members"] = raw_df_translation["accession"]
        raw_df_translation["clusters"] = raw_df_translation["clusters"].map(rep_map)
        cleaned_df = raw_df_translation[["members", "clusters"]]
        cleaned_df.set_index("members")
        return cleaned_df

    @staticmethod
    def add_go_terms_to_dataframe(protein_database_with_clustering_data: pd.DataFrame,
                                  cluster_dicts: dict, filter_for_singletons: bool):
        logging.info(f"Adding GO terms to dataframe. Filtering for singletons {filter_for_singletons}")

        def add_to_dict(cluster_dict, cluster_id, entry_name, go_terms):
            """Helper to append members and GO terms to the cluster dictionary."""
            cluster_dict.setdefault(cluster_id, {"members": [], "go_terms": []})
            cluster_dict[cluster_id]["members"].append(entry_name)
            cluster_dict[cluster_id]["go_terms"].append(go_terms)

        for cluster_names in cluster_dicts.keys():
            protein_database_with_clustering_data[cluster_names] = (
                protein_database_with_clustering_data[cluster_names].str.removeprefix("UniRef90_"))
            protein_database_with_clustering_data[cluster_names] = (
                protein_database_with_clustering_data[cluster_names].str.removeprefix("UniRef50_")
            )
        grouped_exp_go_terms = protein_database_with_clustering_data.dropna()

        # Initialize dicts for each clustering method

        for index, row in tqdm.tqdm(grouped_exp_go_terms.iterrows()):
            entry_name = index
            go_terms = ["GO:" + term for term in row["GO_ID"].split("GO:")[1:]]
            go_terms = np.unique(go_terms).tolist()

            # Add to all cluster dictionaries
            for method in cluster_dicts:
                add_to_dict(cluster_dicts[method], row[method], entry_name, go_terms)

        # Convert dicts to DataFrames and optionally filter singletons
        cluster_dfs = {}
        for method, cdict in cluster_dicts.items():
            df = pd.DataFrame.from_dict(cdict, orient="index")
            if filter_for_singletons:
                df = df[df["members"].apply(len) > 1]
            cluster_dfs[method] = df
        return cluster_dfs

    @staticmethod
    def post_process_mmseqs2_clustering(mmseqs2_clustering_output_path: Path) -> pd.DataFrame:
        logging.info("Post processing MMSeqs2 Clustering results")
        if mmseqs2_clustering_output_path.with_suffix('.parquet').exists():
            logging.info(f"Parquet exists {mmseqs2_clustering_output_path.with_suffix('.parquet').as_posix()}"
                         f", so skipping the post processing.")
            clusters = dd.read_parquet(mmseqs2_clustering_output_path.with_suffix('.parquet').as_posix())
        else:
            clusters = dd.read_csv(mmseqs2_clustering_output_path, header=None, sep="\t")
            clusters = clusters.set_index(1)
            clusters = clusters.rename(columns={0: "clusters"})
            clusters.index = clusters.index.str.removeprefix("UniRef90_")
            clusters.index = clusters.index.str.removeprefix("UniRef100_")
            clusters = clusters[~clusters.index.duplicated(keep="first")]
            clusters.to_parquet(mmseqs2_clustering_output_path.with_suffix('.parquet'))
        return clusters

    @staticmethod
    def post_process_cdhit_xml_dataframe(cd_hit_df_path: Path) -> pd.DataFrame:
        logging.info("Post processing cd-hit xml dataframe.")
        if cd_hit_df_path.parent.joinpath(cd_hit_df_path.name.split('.')[0] + '_processed.parquet').exists():
            return dd.read_parquet(
                cd_hit_df_path.parent.joinpath(
                    cd_hit_df_path.name.split('.')[0] + '_processed.parquet').as_posix()).compute()
        cd_hit_df = dd.read_parquet(cd_hit_df_path.as_posix())
        cd_hit_df = cd_hit_df.set_index("member")
        cd_hit_df = cd_hit_df.rename(columns={"rep": "clusters"})
        cd_hit_df.index = cd_hit_df.index.str.removeprefix("UniRef90_")
        cd_hit_df.index = cd_hit_df.index.str.removeprefix("UniRef100_")
        #cd_hit_df.index     = cd_hit_df.index.str.split("_").str[0]
        cd_hit_df["clusters"] = cd_hit_df["clusters"].str.split("_").str[0]
        cd_hit_df = cd_hit_df.reset_index()
        cd_hit_df = cd_hit_df.drop_duplicates(subset="member", keep="first")
        cd_hit_df = cd_hit_df.set_index("member")
        cd_hit_df.to_parquet(
            cd_hit_df_path.parent.joinpath(cd_hit_df_path.name.split('.')[0] + '_processed.parquet').as_posix())
        return cd_hit_df

    @benchmark_year.setter
    def benchmark_year(self, value):
        self._benchmark_year = value

    @staticmethod
    def post_process_consistency_data(data):
        logging.info("Post processing consistency data")
        data["average"] = data["consistency"].apply(lambda x: np.mean(x))
        data["worst"] = data["consistency"].apply(lambda x: np.min(x))
        data["length"] = data["consistency"].apply(lambda x: len(x))
        return data

    def process_gaf_file(self, gaf_file: Path, individual_components: bool = False):
        """
        Process a GAF (Gene Annotation File) efficiently using Dask on a high-performance machine.

        This function:
        - Reads and filters GAF data using Dask.
        - Groups GO terms per protein (DB_Object_ID) across all aspects (C, F, P).
        - Optionally processes and saves per-aspect (CC, MF, BP) groupings.
        - Writes results to Parquet and returns a lazy Dask DataFrame.

        Parameters
        ----------
        gaf_file : Path
        Path to gzipped GAF/TSV file.
        individual_components : bool, default=False
            If True, also writes and returns per-aspect (C, F, P) Parquets.

        Returns
        -------
        dd.DataFrame
            Aggregated GO terms per DB_Object_ID (and optionally per aspect if saved).
        """
        # --- Setup paths ---
        workdir = Path(benchmark_root_dir) / self.benchmark_working_dir
        workdir.mkdir(parents=True, exist_ok=True)

        p_all = workdir / "goa_file_all.parquet"
        p_cc = workdir / "gaf_file_cc.parquet"
        p_bf = workdir / "gaf_file_bf.parquet"
        p_mf = workdir / "gaf_file_mf.parquet"

        # Determine rebuild condition
        rebuild = (
                not p_all.exists() or
                (individual_components and not all(p.exists() for p in [p_cc, p_bf, p_mf]))
        )

        if rebuild:
            # --- Load GAF file ---
            df = dd.read_csv(
                gaf_file,
                sep="\t",
                engine="pyarrow",
                names=["DB_Object_ID", "GO_ID", "Evidence_Code", "With_or_From", "Aspect"],
                usecols=[1, 4, 6, 7, 8],
                skiprows=9,
                compression="gzip",
                dtype="string[pyarrow]",
                assume_missing=True,

            )

            # Keep relevant columns
            df = df[["DB_Object_ID", "GO_ID", "Evidence_Code", "Aspect"]]

            # Optimize memory usage
            df["Aspect"] = df["Aspect"].astype("category")
            df["Evidence_Code"] = df["Evidence_Code"].astype("category")

            # Filter by evidence codes
            df = df[df["Evidence_Code"].isin(EXP_EVIDENCE_CODES)][["DB_Object_ID", "GO_ID", "Aspect"]]

            # Partition by DB_Object_ID for efficient groupby

            df = df.set_index("DB_Object_ID", shuffle="tasks")
            df["accension"] = df.index
            # --- Aggregation helper ---
            def aggregate_terms(dataframe: dd.DataFrame) -> dd.DataFrame:
                return dataframe.groupby(dataframe.index)["GO_ID"].apply(lambda s: list(set(s)), meta=("GO_terms", "object")).to_frame(name="GO_terms")

            # Aggregate all GO terms
            grouped_all = aggregate_terms(df)

            grouped_all.to_parquet(p_all, write_index=True, engine="pyarrow")

            # Aggregate by aspect if requested
            if individual_components:
                for aspect, path in zip(["C", "P", "F"], [p_cc, p_bf, p_mf]):
                    agg = aggregate_terms(df[df["Aspect"] == aspect])
                    agg.to_parquet(path, write_index=True, engine="pyarrow")

        # --- Load outputs lazily ---
        out_all = dd.read_parquet(p_all.as_posix(), engine="pyarrow")

        if individual_components:
            out_cc = dd.read_parquet(p_cc.as_posix(), engine="pyarrow")
            out_bf = dd.read_parquet(p_bf.as_posix(), engine="pyarrow")
            out_mf = dd.read_parquet(p_mf.as_posix(), engine="pyarrow")
            return out_all, out_cc, out_bf, out_mf

        return out_all
