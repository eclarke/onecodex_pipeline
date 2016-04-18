## Download and parse One Codex analyses

Downloads analyses for a set of samples already uploaded to One Codex, then
expands the taxonomy of all the taxa from NCBI.

**Output**: 
- Two summary tables: a classic row-column matrix of taxa/reads and
  the complete taxonomy for each taxon ID.
- Individual analyses for all samples

### Requirements:
- Python 3.5+
  - Requests
  - BioPython
  - Snakemake

### Running:
```sh
$ snakemake
```

