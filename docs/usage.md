# Usage

Every command reads transcriptome alignments and writes tab-separated tables,
`<prefix>.<table>.tsv`. How each step works: [Method: quant](method.md),
[ORF counting](orfs.md) and [start sites](starts.md).

## Install

```bash
pip install git+https://github.com/katarinagresova/pitstop.git
```

pitstop needs Python 3.10 or later, with numpy, scipy, pandas and pysam. To
run the tests, or to use the exact versions they ran with, use the conda
environment instead:

```bash
git clone https://github.com/katarinagresova/pitstop.git
cd pitstop
mamba env create -p .env -f environment.yml
.env/bin/pip install --no-deps -e .
.env/bin/python -m pytest          # 46 tests on synthetic data with known truth
```

## Counts per CDS: `quant`

```bash
pitstop quant --bam sample.bam --gtf annotation.gtf --fasta genome.fa \
    --read-lengths 18-30 --out-prefix out/sample.human
# spike-in, reusing the human offsets (or leave --offsets out to estimate them from the spike-in reads)
pitstop quant --bam sample.spikein.bam --gtf yeast.gtf --fasta yeast.fa \
    --read-lengths 18-30 --offsets out/sample.human.offsets.tsv --out-prefix out/sample.yeast
```

Inputs:

- `--bam`: alignments to transcripts; reference names are transcript ids. Text
  after a `|` is ignored.
- `--gtf` / `--fasta`: the annotation and its genome. Transcript lengths must
  match the BAM header.
- `--read-lengths`: footprint lengths to use (required).

| Output | Content |
|---|---|
| `<prefix>.quant.tsv` | `Name Length EffectiveLength ritpm NumReads`, one row per CDS, sorted by Name. `Length` = CDS without the stop codon. `NumReads` = expected reads; they sum to the assigned reads. `ritpm` = reads per nt, scaled to sum to 1e6. `NA` where no read is compatible. |
| `<prefix>.offsets.tsv` | per (length, phase): `offset` (nt from the 5' end to the P-site codon), `support` (reads covering the first nt of the start codon or of the last sense codon; a read with n alignments counts 1/n), `z` (how firmly the data pin this phase's offset: score lead of the chosen window over the best window that gives this phase another offset, divided by the square root of the reads the two disagree on; below about 4 = weakly determined), `reads` |
| `<prefix>.ties.tsv` | `Name tie_group`: CDSs that no read tells apart (e.g. identical paralogs). Their split comes from the model, not from the data: equally long CDSs keep the EM's even start, otherwise the shortest gets almost all the reads. |
| `<prefix>.stats.tsv` | reads left after each filter, EM steps and log-likelihood |
| `<prefix>.psites.tsv` | `read Name psite length`, one row per alignment assigned to a CDS: `psite` is the 0-based transcript position of the P-site's first nt |
| `<prefix>.qc.tsv` | `length region phase reads`, for QC plots: every read once, at its footprint length, split evenly over its alignments to transcripts with a CDS (over all its alignments if it has none). `region` = where the P-site is: `leader`, `CDS` or `trailer` (the stop codon and after); `NA` without a CDS, outside `--read-lengths` or without an offset. `phase` = (5' end − CDS start) mod 3, as in `offsets.tsv`, if the 5' end is in the CDS; else `NA`. At every read length. |

## ORFs other than the CDS: `orfs` and `score`

`pitstop orfs` counts reads on a table of ORFs — upstream ORFs (uORFs),
ORFs overlapping the CDS out of frame (uoORFs), and others — alongside the
annotated CDS, and keeps every read's reading frame instead of snapping it to
the CDS. `pitstop score` turns those per-codon counts into frame evidence for
each ORF, pooled over one or more `orfs` runs (e.g. replicate libraries):

```bash
pitstop orfs --bam sample_rep1.bam --gtf annotation.gtf --fasta genome.fa \
    --read-lengths 18-30 --orfs candidates.tsv --out-prefix out/sample_rep1
# the same for sample_rep2, with --offsets out/sample_rep1.offsets.tsv
pitstop score --orfs-prefix out/sample_rep1 out/sample_rep2 --out-prefix out/pooled
```

`--orfs` is a TSV `ORF_id Name start end` (0-based transcript coordinates,
`start` = first nt of the start codon, `end` = one past the last sense codon,
same convention as the CDS). The annotated CDSs are added automatically.

To compare conditions, give every library the same GTF, ORF table and
offsets table, and score all libraries together:
[Typical workflow](orfs.md#typical-workflow). What the scores do and do
not say: [Reading the results](orfs.md#reading-the-results).

| Output | Content |
|---|---|
| `<prefix>.orfs.tsv` | `ORF_id Name type start end start_codon Length NumReads`, one row per ORF (from `--orfs`, plus the annotated CDSs) and per outside component (`<tx>:leader`, `<tx>:trailer`, `<tx>:transcript`: the positions in no ORF). `type` is `CDS`, `uORF` (ends at or before the CDS start), `uoORF` (starts upstream of the CDS, out of frame with it) or `other` for an ORF, and `leader`, `trailer` or `transcript` for an outside component. |
| `<prefix>.offsets.tsv` | the same table as `quant`'s (the same BAM and read lengths give the same offsets); either command's `--offsets` accepts it. `orfs` uses only each length's phase-0 offset. |
| `<prefix>.frames.tsv` | per read length: `offset` (the phase-0 offset, not a window), `reads`, and `frame0 frame1 frame2` — the share of P-sites in each frame of the CDS, from CDS interiors. |
| `<prefix>.codons.tsv` | `ORF_id codon length frame0 frame1 frame2`: P-sites per codon of every ORF but the annotated CDSs, by length and frame — `score`'s input. |
| `<prefix>.ties.tsv`, `<prefix>.stats.tsv`, `<prefix>.psites.tsv` | as for `quant`, but over every ORF and outside component: `ties.tsv` is `ORF_id tie_group`, and `psites.tsv` keeps every P-site, not only those assigned to a CDS. |
| `<prefix>.scores.tsv` | one row per scored ORF: `ORF_id Name type start end codons codons_with_reads reads in_frame_share leads expected z p min_p q` — `p`/`q` are the frame-evidence significance and its BH correction; `min_p` is the best `p` the ORF's length and reads could reach. |
| `<prefix>.decoys.tsv` | each ORF shifted by +1 and +2 nt and scored the same way, for null calibration. |
| `<prefix>.score_stats.tsv` | the read lengths used and how many ORFs/decoys were scored. |

## Start sites: `starts`

`pitstop starts` finds translation start sites in harringtonine (or
lactimidomycin) libraries. Each harringtonine library is compared with a
matched elongation library (same condition and replicate), so that a pause,
which both have, is not taken for a start. It reads the outputs of
`pitstop orfs` runs:

```bash
# one offset table per pool: quant on the merged harringtonine BAMs
pitstop quant --bam harr_merged.bam --gtf annotation.gtf --fasta genome.fa \
    --read-lengths 18-30 --out-prefix out/pool
# orfs on every library, harringtonine and elongation, with that table
pitstop orfs --bam harr_rep1.bam --gtf annotation.gtf --fasta genome.fa \
    --read-lengths 18-30 --offsets out/pool.offsets.tsv --out-prefix out/harr_rep1
# ... harr_rep2, elong_rep1, elong_rep2 the same way
pitstop starts --libraries libraries.tsv --gtf annotation.gtf --fasta genome.fa \
    --scan --out-prefix out/starts
```

`--libraries` is a TSV `pool harringtonine elongation`, one row per
replicate, each an out prefix of an `orfs` run. A pool is a set of
replicates of one condition. Several pools can be given, and their evidence
adds up. `--scan` tests every ATG and near-cognate start codon (CTG, GTG,
TTG, AAG, ACG, AGG, ATA, ATC, ATT) with harringtonine reads, and calls the
starts with BH q < 0.05; `--start-codons ATG,CTG` tests only these codons and
`--fdr 0.01` sets the cut. `--orfs table.tsv` tests the starts of an ORF
table instead (columns `ORF_id Name start`). The annotated starts are always
tested. How it works and how to read it: [Start sites](starts.md).

| Output | Content |
|---|---|
| `<prefix>.starts.tsv` | one row per tested start: `Name start codon region frame ORF_id harringtonine elongation expected enrichment p typical_p q multimapped shoulder`. `harringtonine` and `elongation` are the P-sites in the start's window (nt -1 to +2 of the start codon), summed over the replicates and pools; `expected` is the harringtonine P-sites that the elongation profile predicts there; `q` is BH per class (annotated starts; other starts per region: leader, CDS, trailer, transcript without a CDS); `typical_p` is the p the window would have with the reads of a typical used start; `multimapped` is the share of the window's reads with another alignment; `shoulder` is log2 harringtonine / elongation P-sites (+1 each) in the codons 1-5 after the start, minus the same in the 5 codons before its window (a used start raises its codons 1-5; not part of the test). With `--scan` also `called stopped_by end type`. |
| `<prefix>.start_orfs.tsv` | with `--scan`: the ORF of each called start, to its first in-frame stop, but the annotated CDSs: `ORF_id Name start end type start_codon` and the start evidence. `ORF_id` is `<transcript>_<1-based start>`. `type` is `uORF`, `uoORF`, `extension` (in frame upstream of the CDS), `truncation` (in frame inside it), `internal` (out of frame inside it), `dORF` (after it) or `other` (no annotated CDS). An ORF table for `pitstop orfs --orfs`. |
| `<prefix>.catalogue.tsv` | with `--scan`: the uORFs and uoORFs of the called starts, and every annotated CDS with its start evidence; the columns of `start_orfs.tsv`. An ORF table for `pitstop orfs --orfs`. |
| `<prefix>.kernel.tsv` | per pool, read length and nt from the annotated start (-15 to +32): harringtonine and elongation P-sites, and the kernel (harringtonine minus scaled elongation, as a share). |
| `<prefix>.factors.tsv` | per pool: the harringtonine / elongation factor of each bin of region and distance. |
| `<prefix>.starts_stats.tsv` | the start lengths, P-site totals, the prior and overdispersion per region, and per class the tests, the calls and the false calls expected from null windows. With `--scan`, the overdispersion is the second one, and `null_windows_near_called` is the number of null windows dropped near the starts of the first call. |
