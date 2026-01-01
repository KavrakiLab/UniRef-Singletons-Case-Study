import logging
import xml.etree.ElementTree as ET
import tqdm
from xml.etree import cElementTree
import pandas as pd
from pathlib import Path
from .const import benchmark_root_dir
import pyarrow as pa
import pyarrow.parquet as pq


def xml_to_fasta(xml_path: str, output_fasta_path: str) -> None:
    logging.info("Processing XML to FASTA")
    prev_seq = None
    # parser = etree.XMLParser(encoding="utf-8")
    if Path(output_fasta_path).exists():
        return None
    with open(output_fasta_path, 'w') as out:
        for event, elem in tqdm.tqdm(cElementTree.iterparse(xml_path, events=('end',))):
            # Adjust the parent element as needed based on your XML structure
            if elem.tag == '{http://uniprot.org/uniref}representativeMember':
                accession = elem[0][0].attrib["value"]
                for child in elem:
                    if child.tag == "{http://uniprot.org/uniref}sequence":
                        prev_seq = child.text.strip()
                        out.write(f">{accession}\n{prev_seq}\n")
            #                elem.clear()
            elif elem.tag == '{http://uniprot.org/uniref}member':
                if elem[0].attrib["type"] == "UniProtKB ID":
                    accession = elem[0][0].attrib["value"]
                elif elem[0].attrib["type"] == "UniParc ID":
                    accession = elem[0].attrib["id"]
                out.write(f">{accession}\n{prev_seq}\n")
            elif elem.tag == '{http://uniprot.org/uniref}entry':
                elem.clear()


def parse_rep_string(string_snippet):
    match = re.search(r'id="([^"]+)"', string_snippet)
    protein_id = (match.group(1))  #
    return protein_id


def parse_entry_string(string_snippet):
    match = re.search(r'value="([^"]+)"', string_snippet)
    protein_id = (match.group(1))  #
    return protein_id


def parse_grep_output(grep_file, parquet_file_output):
    dataframe_dict = {"member": [], "rep": []}
    if parquet_file_output.exists():
        return parquet_file_output
    with open(grep_file, 'r') as grep_output:
        for line in tqdm.tqdm(grep_output):
            if line == "<representativeMember>\n":
                db_reference = next(grep_output)
                rep_info = parse_entry_string(next(grep_output))
                dataframe_dict["rep"].append(rep_info)
                dataframe_dict["member"].append(rep_info)
            elif line == "<member>\n":
                db_reference = next(grep_output)
                member_info = parse_entry_string(next(grep_output))
                dataframe_dict["rep"].append(rep_info)
                dataframe_dict["member"].append(member_info)
    #            member_info = parse_entry_string(next(grep_output))
    #            dataframe_dict["member"][-1].append(member_info)
    pd.DataFrame(dataframe_dict).to_parquet(parquet_file_output)
    return parquet_file_output


def find_entry_offset(xml_path: str, entry_id: str) -> int:
    """Find the byte offset of a given UniRef entry ID in the XML file."""
    with open(xml_path, "rb") as f:
        for line in f:
            if entry_id.encode() in line:
                return f.tell() - len(line)
    return 0

#def run_grep_command(xml_file_path: str):
#    """Runs the command"""

#    "grep -A 2 -e '<representativeMember>' -e '<member>' tmp > tmp_grep"
def parse_uniref_xmls(xml_path: str, parquet_path: str, chunk_size: int = 500_000,
                      checkpoint_path: str = "checkpoint.json") -> Path:
    """
    Stream-parse UniRef XML into parquet with checkpointing + fast restart.
    Instead of re-parsing whole file, seeks to last written entry.
    """
    final_path = Path(parquet_path)
    checkpoint_file = Path(checkpoint_path)

    # Load checkpoint
    last_id = None
    if checkpoint_file.exists():
        with open(checkpoint_file) as f:
            last_id = json.load(f).get("last_id", None)
        logging.info(f"Resuming from last processed ID: {last_id}")

    # Open XML file
    xml_fh = open(xml_path, "rb")
    if last_id:
        offset = find_entry_offset(xml_path, last_id)
        logging.info(f"Seeking to byte offset {offset}")
        xml_fh.seek(offset)

    parser = ET.XMLPullParser("end")
    # context = cElementTree.iterparse(xml_fh, events=("end",))
    buffer = []
    writer = None
    last_seen = None

    for _, elem in tqdm.tqdm(context, desc="Parsing UniRef XML"):
        if elem.tag == "{http://uniprot.org/uniref}entry":
            rep = elem.attrib["id"]
            last_seen = rep

            # Skip the checkpoint entry itself
            if last_id and rep == last_id:
                elem.clear()
                last_id = None
                continue

            buffer.append({"rep": rep, "member": rep.split("_")[-1]})

            for child in elem:
                if child.tag == "{http://uniprot.org/uniref}member":
                    buffer.append({"rep": rep, "member": child[-1][0].attrib["value"]})
                child.clear()
            elem.clear()

            if len(buffer) >= chunk_size:
                df = pd.DataFrame(buffer)
                table = pa.Table.from_pandas(df, preserve_index=False)
                if writer is None:
                    writer = pq.ParquetWriter(final_path, table.schema,
                                              use_dictionary=True, compression="snappy")
                writer.write_table(table)
                buffer = []

                # Save checkpoint
                with open(checkpoint_file, "w") as f:
                    json.dump({"last_id": last_seen}, f)

    # Final flush
    if buffer:
        df = pd.DataFrame(buffer)
        table = pa.Table.from_pandas(df, preserve_index=False)
        if writer is None:
            writer = pq.ParquetWriter(final_path, table.schema,
                                      use_dictionary=True, compression="snappy")
        writer.write_table(table)
        with open(checkpoint_file, "w") as f:
            json.dump({"last_id": last_seen}, f)

    if writer is not None:
        writer.close()

    logging.info(f"Finished writing parquet: {final_path}")
    checkpoint_file.unlink(missing_ok=True)
    return final_path


from pathlib import Path
import requests
import tarfile
import gzip
import shutil


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


import re
import time
import json
import zlib
from xml.etree import ElementTree
from urllib.parse import urlparse, parse_qs, urlencode
import requests
from requests.adapters import HTTPAdapter, Retry

POLLING_INTERVAL = 3
API_URL = "https://rest.uniprot.org"

retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
session = requests.Session()
session.mount("https://", HTTPAdapter(max_retries=retries))


def check_response(response):
    try:
        response.raise_for_status()
    except requests.HTTPError:
        print(response.json())
        raise


def submit_id_mapping(from_db, to_db, ids):
    request = requests.post(
        f"{API_URL}/idmapping/run",
        data={"from": from_db, "to": to_db, "ids": ",".join(ids)},
    )
    check_response(request)
    job_id = request.json()["jobId"]
    if check_id_mapping_results_ready(job_id):
        link = get_id_mapping_results_link(job_id)
        results = get_id_mapping_results_search(link)
        return results


def get_next_link(headers):
    re_next_link = re.compile(r'<(.+)>; rel="next"')
    if "Link" in headers:
        match = re_next_link.match(headers["Link"])
        if match:
            return match.group(1)


def check_id_mapping_results_ready(job_id):
    while True:
        request = session.get(f"{API_URL}/idmapping/status/{job_id}")
        check_response(request)
        j = request.json()
        if "jobStatus" in j:
            if j["jobStatus"] in ("NEW", "RUNNING"):
                print(f"Retrying in {POLLING_INTERVAL}s")
                time.sleep(POLLING_INTERVAL)
            else:
                raise Exception(j["jobStatus"])
        else:
            return bool(j["results"] or j["failedIds"])


def get_batch(batch_response, file_format, compressed):
    batch_url = get_next_link(batch_response.headers)
    while batch_url:
        batch_response = session.get(batch_url)
        batch_response.raise_for_status()
        yield decode_results(batch_response, file_format, compressed)
        batch_url = get_next_link(batch_response.headers)


def combine_batches(all_results, batch_results, file_format):
    if file_format == "json":
        for key in ("results", "failedIds"):
            if key in batch_results and batch_results[key]:
                all_results[key] += batch_results[key]
    elif file_format == "tsv":
        return all_results + batch_results[1:]
    else:
        return all_results + batch_results
    return all_results


def get_id_mapping_results_link(job_id):
    url = f"{API_URL}/idmapping/details/{job_id}"
    request = session.get(url)
    check_response(request)
    return request.json()["redirectURL"]


def decode_results(response, file_format, compressed):
    if compressed:
        decompressed = zlib.decompress(response.content, 16 + zlib.MAX_WBITS)
        if file_format == "json":
            j = json.loads(decompressed.decode("utf-8"))
            return j
        elif file_format == "tsv":
            return [line for line in decompressed.decode("utf-8").split("\n") if line]
        elif file_format == "xlsx":
            return [decompressed]
        elif file_format == "xml":
            return [decompressed.decode("utf-8")]
        else:
            return decompressed.decode("utf-8")
    elif file_format == "json":
        return response.json()
    elif file_format == "tsv":
        return [line for line in response.text.split("\n") if line]
    elif file_format == "xlsx":
        return [response.content]
    elif file_format == "xml":
        return [response.text]
    return response.text


def get_xml_namespace(element):
    m = re.match(r"\{(.*)\}", element.tag)
    return m.groups()[0] if m else ""


def merge_xml_results(xml_results):
    merged_root = ElementTree.fromstring(xml_results[0])
    for result in xml_results[1:]:
        root = ElementTree.fromstring(result)
        for child in root.findall("{http://uniprot.org/uniprot}entry"):
            merged_root.insert(-1, child)
    ElementTree.register_namespace("", get_xml_namespace(merged_root[0]))
    return ElementTree.tostring(merged_root, encoding="utf-8", xml_declaration=True)


def print_progress_batches(batch_index, size, total):
    n_fetched = min((batch_index + 1) * size, total)
    print(f"Fetched: {n_fetched} / {total}")


def get_id_mapping_results_search(url):
    parsed = urlparse(url)
    query = parse_qs(parsed.query)
    file_format = query["format"][0] if "format" in query else "json"
    if "size" in query:
        size = int(query["size"][0])
    else:
        size = 500
        query["size"] = size
    compressed = (
        query["compressed"][0].lower() == "true" if "compressed" in query else False
    )
    parsed = parsed._replace(query=urlencode(query, doseq=True))
    url = parsed.geturl()
    request = session.get(url)
    check_response(request)
    results = decode_results(request, file_format, compressed)
    total = int(request.headers["x-total-results"])
    print_progress_batches(0, size, total)
    for i, batch in enumerate(get_batch(request, file_format, compressed), 1):
        results = combine_batches(results, batch, file_format)
        print_progress_batches(i, size, total)
    if file_format == "xml":
        return merge_xml_results(results)
    return results


def get_id_mapping_results_stream(url):
    if "/stream/" not in url:
        url = url.replace("/results/", "/results/stream/")
    request = session.get(url)
    check_response(request)
    parsed = urlparse(url)
    query = parse_qs(parsed.query)
    file_format = query["format"][0] if "format" in query else "json"
    compressed = (
        query["compressed"][0].lower() == "true" if "compressed" in query else False
    )
    return decode_results(request, file_format, compressed)


def parse_dat_file_for_fasta(dat_file_path, fasta_file_output):
    dat_file = open(dat_file_path, 'r')
    fasta_output = open(fasta_file_output, 'w')

    read_seq = False
    fasta = []
    for line in tqdm.tqdm(dat_file):
        if line.startswith("AC"):
            entry = line.split()[-1][:-1]
        elif line.startswith("SQ"):
            read_seq = True
        elif read_seq:
            if not "//" in line:
                fasta.append("".join(line.split()) + "\n")
            else:
                read_seq = False
                fasta_output.write(">" + entry + "\n")
                [fasta_output.write(seq) for seq in fasta]
                fasta = []
