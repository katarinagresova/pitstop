<div align="center">

[![tests](https://github.com/katarinagresova/pitstop/actions/workflows/tests.yml/badge.svg)](https://github.com/katarinagresova/pitstop/actions/workflows/tests.yml)
[![python](https://img.shields.io/python/required-version-toml?tomlFilePath=https://raw.githubusercontent.com/katarinagresova/pitstop/main/pyproject.toml)](pyproject.toml)
[![docs](https://img.shields.io/badge/docs-online-E5412D)](https://katarinagresova.github.io/pitstop/)

<p><picture><img alt="pitstop" src="docs/assets/pitstop-banner.svg" width="300" align="top"></picture></p>

**Ribo-seq reads in, translation out.**

</div>

<picture>
<img alt="pitstop at a glance: reading frames drawn as the three lanes of a race track. Ribosomes keep to their ORF's lane, so a uORF whose reads stay in its frame is translated (pitstop orfs + score); harringtonine holds ribosomes on start codons (pitstop starts); reads are counted per coding sequence (pitstop quant); each library is checked first, for the read lengths that keep their lane (pitstop inspect); and extra 5' nucleotides are clipped from its reads (pitstop repair)." src="docs/assets/pitstop-abstract-auto.svg" width="100%">
</picture>

pitstop is a Python toolkit for ribosome profiling (Ribo-seq). Give it
footprints aligned to the transcriptome, and it gives you plain tables to take
straight into your analysis.

## What it does

| Question | Command | You get |
|---|---|---|
| Is my library good, and which read lengths do I use? | `pitstop inspect` | an HTML report and tables: for each read length, how well its reads keep one reading frame, flags for common library problems, and the read lengths to use |
| Can my reads be repaired? | `pitstop repair` | a new BAM in which extra 5' nucleotides that do not match the transcript are clipped, and a table of the reads it changed |
| How much is each coding sequence translated? | `pitstop quant` | reads per CDS, with reads shared by isoforms or paralogs split between them |
| Which uORFs and other small ORFs are translated? | `pitstop orfs`, `pitstop score` | for each ORF in your table, how strongly its reads follow its reading frame, with a p-value |
| Where does translation start? | `pitstop starts` | start sites at ATG and near-cognate codons from harringtonine libraries, and a uORF catalogue from them |

<picture>
<img alt="One transcript with a uORF (CUG start) and a CDS (AUG start), and two read tracks above it. Harringtonine reads peak at both start codons, above what the Ribo-seq reads predict (pitstop starts). Ribo-seq P-sites are coloured by reading frame: on the uORF most keep the uORF's frame (pitstop orfs + score), on the CDS the CDS's frame (pitstop quant). A bar chart of the read lengths of a library marks the lengths whose reads keep one frame (pitstop inspect). Reads with an extra 5' nucleotide are shown before and after it is clipped (pitstop repair)." src="docs/assets/pitstop-transcript-auto.svg" width="100%">
</picture>

## Why pitstop

- **Learns from your data.** The P-site offset of each read length is
  estimated from the reads, not assumed.
- **Keeps the reading frame.** For uORFs, every read keeps its frame, and the
  score allows for the off-frame noise of the coding sequence they overlap.
- **Reproducible.** No random numbers anywhere: the same input gives the same
  output, byte for byte.
- **Checked on real data**, not only on simulations: see how the
  [uORF scores](https://katarinagresova.github.io/pitstop/orfs/#validation-on-real-data) and
  [start sites](https://katarinagresova.github.io/pitstop/starts/#validation-on-real-data) hold up.

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

Not sure which `--read-lengths` to use? `pitstop inspect` checks a library and
names the read lengths that keep their reading frame. If it finds extra
nucleotides at the 5' end of the reads, let the aligner clip them the next time
(STAR with `--quantTranscriptomeSAMoutput BanSingleEnd`, or bowtie2 `--local`).
For the BAM you have, `pitstop repair` clips them and writes a new BAM for the
other commands.

To test candidate uORFs, count reads on your ORF table in each replicate, then
score the ORFs over all of them:

```bash
pitstop orfs --bam rep1.bam --gtf annotation.gtf --fasta genome.fa \
    --read-lengths 18-30 --orfs candidates.tsv --out-prefix out/rep1
# the same for rep2, with --offsets out/rep1.offsets.tsv: libraries scored together share one offsets table
pitstop score --orfs-prefix out/rep1 out/rep2 --out-prefix out/pooled
```

`out/pooled.scores.tsv` has a p- and q-value for each ORF. For start sites
from harringtonine libraries, see [Usage](https://katarinagresova.github.io/pitstop/usage/#start-sites-starts).

## Documentation

**[katarinagresova.github.io/pitstop](https://katarinagresova.github.io/pitstop/)**

- [Usage](https://katarinagresova.github.io/pitstop/usage/): every command, input format and output column
- [Read QC flags](https://katarinagresova.github.io/pitstop/flags/): what each flag of `pitstop inspect` means,
  its cut and why, and what to do
- [Typical workflows](https://katarinagresova.github.io/pitstop/#typical-workflows): which commands to run, in which
  order, for elongation libraries, initiation and elongation libraries, and comparisons of conditions
- [Method: quant](https://katarinagresova.github.io/pitstop/method/): offsets, assignment and the EM, with interactive examples
- [Method: ORF counting](https://katarinagresova.github.io/pitstop/orfs/): the frame score, with an interactive example and
  validation on real data
- [Method: start sites](https://katarinagresova.github.io/pitstop/starts/): the harringtonine test, with validation on real data

## License

Apache License 2.0, see [LICENSE](LICENSE) and [NOTICE](NOTICE).
