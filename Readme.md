## Download and parse One Codex analyses

Downloads analyses for a set of samples already uploaded to One Codex, then
expands the taxonomy of all the taxa from NCBI and writes summaries.

**Output**: 

- Two summary tables: a classic row-column matrix of taxa/reads and
  the complete taxonomy for each taxon ID.
- Individual analyses for all samples

### Requirements:
- Python 3.5+
    - Requests
    - BioPython
    - Snakemake
- One Codex API key stored in an environmental variable called `$ONE_CODEX_API_KEY` or in
`snakemake_config.yaml` under the key "api_key".

### Running:
```sh
$ snakemake
```

