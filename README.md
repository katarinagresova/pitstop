<div align="center">

[![tests](https://github.com/katarinagresova/pitstop/actions/workflows/tests.yml/badge.svg)](https://github.com/katarinagresova/pitstop/actions/workflows/tests.yml)
[![python](https://img.shields.io/python/required-version-toml?tomlFilePath=https://raw.githubusercontent.com/katarinagresova/pitstop/main/pyproject.toml)](pyproject.toml)
[![docs](https://img.shields.io/badge/docs-online-E5412D)](https://katarinagresova.github.io/pitstop/)

<img alt="pitstop" src="docs/assets/pitstop-banner.svg" width="300">

**From Ribo-seq reads to translation: how much each gene is translated, which
upstream ORFs are translated, and where translation starts.**

</div>

pitstop is a Python toolkit for ribosome profiling (Ribo-seq). Give it
footprints aligned to the transcriptome, and it gives you plain tables to take
straight into your analysis.

## What it does

| Question | Command | You get |
|---|---|---|
| How much is each coding sequence translated? | `pitstop quant` | reads per CDS, with reads shared by isoforms or paralogs split between them |
| Which uORFs and other small ORFs are translated? | `pitstop orfs`, `pitstop score` | for each ORF in your table, how strongly its reads follow its reading frame, with a p-value |
| Where does translation start? | `pitstop starts` | start sites at ATG and near-cognate codons from harringtonine libraries, and a uORF catalogue from them |

## Why pitstop

- **Learns from your data.** The P-site offset of each read length is
  estimated from the reads, not assumed.
- **Keeps the reading frame.** For uORFs, every read keeps its frame, and the
  score allows for the off-frame noise of the coding sequence they overlap.
- **Reproducible.** No random numbers anywhere: the same input gives the same
  output, byte for byte.
- **Checked on real data**, not only on simulations. The docs show how.

## Install

```bash
pip install git+https://github.com/katarinagresova/pitstop.git
```

pitstop needs Python 3.10 or later.

## Quick start

pitstop reads footprints aligned to the transcript sequences (for example with
STAR or bowtie), and the GTF and genome FASTA the transcripts come from.

```bash
pitstop quant --bam sample.bam --gtf annotation.gtf --fasta genome.fa \
    --read-lengths 18-30 --out-prefix out/sample
```

This writes `out/sample.quant.tsv` (reads per CDS), `out/sample.offsets.tsv`
(the P-site offset of each read length) and a few more tables, among them one
for quality-control plots.

To test candidate uORFs, count reads on your ORF table in each replicate, then
score the ORFs over all of them:

```bash
pitstop orfs --bam rep1.bam --gtf annotation.gtf --fasta genome.fa \
    --read-lengths 18-30 --orfs candidates.tsv --out-prefix out/rep1
# the same for rep2, with --offsets out/rep1.offsets.tsv
pitstop score --orfs-prefix out/rep1 out/rep2 --out-prefix out/pooled
```

`out/pooled.scores.tsv` has a p- and q-value for each ORF. For start sites
from harringtonine libraries, see [Usage](https://katarinagresova.github.io/pitstop/usage/#start-sites-starts).

## Documentation

**[katarinagresova.github.io/pitstop](https://katarinagresova.github.io/pitstop/)**

- [Usage](https://katarinagresova.github.io/pitstop/usage/): every command, input format and output column
- [Method: quant](https://katarinagresova.github.io/pitstop/method/), [ORF counting](https://katarinagresova.github.io/pitstop/orfs/) and
  [start sites](https://katarinagresova.github.io/pitstop/starts/): how each step works, with interactive examples
  and validation on real data

## License

Apache License 2.0, see [LICENSE](LICENSE) and [NOTICE](NOTICE).
