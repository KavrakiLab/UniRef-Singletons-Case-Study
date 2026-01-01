import logging

import typer

from src.Dataset_workflow_2016 import UniProt2016Workflow
from src.Dataset_workflow_2025 import UniProt2025Workflow
from lib.const import CliOPTIONS

logging.basicConfig(
    level=logging.DEBUG,  # Set minimum level to capture
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)


def main(year_selection: CliOPTIONS, mmseqs2_executable_path, cd_hit_executable_path):
    if year_selection.value == CliOPTIONS.Benchmark2016.value:
        UniProt2016Workflow(mmseqs2_executable_path, cd_hit_executable_path,
                            '/home/Users/fmq1/PycharmProjects/UniRefAnalysis/go-basic.obo').run()
    elif year_selection.value == CliOPTIONS.Benchmark2025.value:
        UniProt2025Workflow(mmseqs2_executable_path, cd_hit_executable_path,
                            '/home/Users/fmq1/PycharmProjects/UniRefAnalysis/go-basic.obo').run()


typer.run(main)
