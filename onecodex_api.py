import os
import sys
import time
import pprint
import itertools
import functools
import csv

from urllib.error import HTTPError
from enum import Enum

import requests
from Bio import Entrez

Entrez.email = "ecl@mail.med.upenn.edu"

OCX_API = "https://app.onecodex.com/api/v0/"

TaxRanks = [
    'superkingdom',
    'phylum',
    'class',
    'order',
    'family',
    'genus',
    'species'
] 

def _get_ocx_url(url, api_key):
    logging.debug(url)
    r = requests.get(url, auth=(api_key, ''))
    r.raise_for_status()
    return r


@functools.lru_cache()
def get_samples(api_key):
    """List samples on One Codex."""
    url = OCX_API + "samples"
    r = _get_ocx_url(url, api_key)
    return {s['filename']:s['id'] for s in r.json()}


@functools.lru_cache()
def get_analyses(api_key):
    """List completed analyses on One Codex."""
    url = OCX_API + "analyses"
    return _get_ocx_url(url, api_key).json()


def get_analyses_for_id(sample_id, api_key):
    """List the analyses for a sample id."""
    analyses = get_analyses(api_key)
    analysis_ids = [_['id'] for _ in analyses if _['sample_id'] == sample_id]
    for a_id in analysis_ids:
        url = OCX_API + "analyses/" + a_id
        yield _get_ocx_url(url, api_key).json()


def get_analyses_for_sample(sample_name, api_key):
    """List the analyses for a sample."""
    samples = get_samples(api_key)
    sample_id = samples[sample_name]
    return get_analyses_for_id(sample_id, api_key)


def get_raw_ocx_analysis_for_sample(sample_name, out_fp, api_key):
    """Download the raw analysis file for a sample."""
    analyses = get_analyses_for_sample(sample_name, api_key)
    for analysis in analyses:
        if analysis['reference_name'] == 'One Codex Database':
            url = OCX_API + "analyses/{}/raw".format(analysis['id'])
            out = _get_ocx_url(url, api_key)
            out_filename = os.path.join(out_fp, sample_name + ".ocx.raw.tsv.gz")
            with open(out_filename, 'wb') as outfile:
                logging.info("Downloading OCX analysis for {} to {}".format(
                    sample_name,
                    out_filename))
                for chunk in out.iter_content(1024):
                    outfile.write(chunk)
            return out_filename

        
def get_ocx_analysis_for_sample(sample_name, api_key):
    """Get the table-form results for a sample (warning: can be large)."""
    analyses = get_analyses_for_sample(sample_name, api_key)
    for analysis in analyses:
        if analysis['reference_name'] == 'One Codex Database':
            url = OCX_API + "analyses/{}/table".format(analysis['id'])
            table = _get_ocx_url(url, api_key).json()
            return table

        
def get_taxa_in_sample(sample_name, out_fp, api_key):
    """Write table with reads and full tax. info for a sample."""
    taxa = get_ocx_analysis_for_sample(sample_name, api_key)
    tax_ids = [_['tax_id'] for _ in taxa]
    # Post list to NCBI
    tax_info = list(itertools.chain.from_iterable(_ncbi_get_many_taxa(tax_ids)))
    tax_info = [{r['Rank']: r['ScientificName'] for r in t.get('LineageEx', ())} for t in tax_info]
    print(len(tax_info), len(taxa))
    assert(len(tax_info) == len(taxa))
    out_filename = os.path.join(out_fp, sample_name + ".taxa.tsv")
    with open(out_filename, 'w') as out:
        writer = csv.DictWriter(
            out,
            fieldnames=['tax_id', 'name', 'readcount'] + TaxRanks,
            delimiter='\t',
            quoting=csv.QUOTE_MINIMAL,
            extrasaction='ignore',
            restval='NA'
        )
        writer.writeheader()
        rows = ({**x, **y} for x, y in zip(taxa, tax_info))
        writer.writerows(rows)
    return out_filename


def _ncbi_get_many_taxa(ids, batch_size=5000):
    saved = Entrez.read(Entrez.epost('taxonomy', id=','.join(map(str, ids))))
    webenv = saved['WebEnv']
    query_key = saved['QueryKey']
    for start in range(0, len(ids), batch_size):
        end = min(len(ids), start + batch_size)
        attempt = 1
        while attempt <= 3:
            try:
                handle = Entrez.efetch(
                    db='taxonomy', retstart=start, retmax=batch_size,
                    webenv=webenv, query_key=query_key)
                break
            except HTTPError as err:
                if 500 <= err.code <= 599:
                    print("NCBI server error ({}); retrying ({}/3)".format(err.code, attempt))
                    attempt += 1
                    time.sleep(15)
                else:
                    raise
        yield(Entrez.read(handle))
    

def main():
    parser = argparse.ArgumentParser("Download and clean up analyses from OneCodex.")

    parser.add_argument(
        'files', nargs='+',
        help='Files to which to download corresponding analyses.')
    parser.add_argument(
        '--output_fp', help='Folder to write analysis tables. (default: %(default))',
        default='one_codex')
    parser.add_argument(
        '--api_key', help='One Codex API key (by default, read from $ONE_CODEX_API_KEY)',
        default=os.environ.get("ONE_CODEX_API_KEY"))
    parser.add_argument(
        '--verbose', help='Print lots of stuff', action='store_true')

    args = parser.parse_args()

    if not os.path.exists(args.output_fp):
        os.mkdir(args.output_fp)

    # Set up logging
    import logging
    log = logging.getLogger("ocx_analyses")
    log.setLevel(logging.DEBUG)

    log_format = logging.Formatter("%(asctime)s %(levelname)s: %(message)s")

    stream_log = logging.StreamHandler()
    if args.verbose:
        stream_log.setLevel(logging.DEBUG)
    else:
        stream_log.setLevel(logging.INFO)
    stream_log.setFormatter(log_format)
    log.addHandler(stream_log)
    
    file_log = logging.FileHandler(args.output_fp + "/log.txt", 'w')
    file_log.setLevel(logging.DEBUG)
    file_log.setFormatter(log_format)
    log.addHandler(file_log)

    log.info(vars(args))

    
if __name__ == "__main__":
    main()
