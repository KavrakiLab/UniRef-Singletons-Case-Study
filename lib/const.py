import enum
from pathlib import Path
import os
import typer
from goatools.obo_parser import GODag

root_dir: Path = Path(os.path.abspath(__file__)).parent.parent
benchmark_root_dir: Path = root_dir.joinpath("benchmarks")
EXP_EVIDENCE_CODES = ["EXP", "IDA", "IPI", "IMP", "IGI", "IEP", "HTP", "HDA", "HMP", "HGI", "HEP", "IC"]
COMPUTATIONAL_EVIDENCE_CODES = ["ISS", "ISO","ISA","ISM","IGC","RCA"]
GO_DAG = GODag('/home/Users/fmq1/PycharmProjects/UniRefAnalysis/go-basic.obo')


class CliOPTIONS(enum.Enum):
    Benchmark2016 = '2016'
    Benchmark2025 = '2025'


benchmark_working_dir = {
    CliOPTIONS.Benchmark2025.value: '2025UniRefResults',
    CliOPTIONS.Benchmark2016.value: '2016UniRefResults'
}
