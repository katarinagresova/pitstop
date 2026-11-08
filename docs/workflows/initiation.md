# Initiation and elongation libraries

The workflow for initiation libraries (harringtonine or lactimidomycin) with
matched elongation libraries: start sites, the ORFs they make, and the
translation of these ORFs.

## The experiment

- **Libraries.** For each replicate, two libraries from the same cells:
  - an **initiation library**: the cells get harringtonine or lactimidomycin
    for a few minutes before lysis. Ribosomes that start translation in that
    time stay on the start codon; the others run off the mRNA.
  - an **elongation library**: the same cells without the drug (or with
    cycloheximide only).
- **Make both with the same nuclease digestion and size selection.** `starts`
  compares the two libraries at each read length. If they were digested
  differently, the initiation reads have few elongation reads of the same
  length behind them, and the test loses its background ([step
  2](#2-qc-of-each-library-and-of-each-pair)).
- **Questions.** Where does translation start, at ATG and at near-cognate
  codons? Which uORFs and other ORFs do these starts make? Are they
  translated?
- **Not covered.** Initiation libraries without a matched elongation library,
  for example harringtonine with puromycin and no elongation library: `starts`
  needs an elongation library for each initiation library.

![The initiation and elongation workflow, one card per step. Step 2: pitstop quant on each library writes qc.tsv, for the QC of each library and each pair; it guides the choice of read lengths. Step 3: pitstop quant on the pool's elongation libraries, merged, writes offsets.tsv. Step 4: pitstop orfs on each library, with no ORF table and these offsets, writes psites.tsv. Step 5: pitstop starts --scan with libraries.tsv, then pitstop plot, writes catalogue.tsv and lengths.pdf. Step 6: pitstop orfs --orfs with the catalogue and the pool offsets on each elongation library, then pitstop score, writes orfs.tsv and scores.tsv.](../img/workflow-initiation.svg)

`orfs` runs two times. The first run, with no ORF table, gives the P-sites
that `starts` reads. The second run counts the reads of the ORFs that `starts`
found.

## Inputs

- A transcriptome BAM per library, all made the same way: [Preparing the
  reads](reads.md). Without UMIs, do not collapse identical reads: the start
  peaks of initiation libraries are many identical reads.
- The GTF and the genome FASTA of the alignment.
- A library sheet (step 5).

## Steps

The examples have one condition (pool `ctrl`) with two replicates:
`harr_1`, `harr_2` (initiation) and `elong_1`, `elong_2` (elongation).

### 1. Prepare the reads

[Preparing the reads](reads.md), the same for both kinds of library.

### 2. QC of each library and of each pair

```bash
for lib in harr_1 harr_2 elong_1 elong_2; do
    pitstop quant --bam $lib.bam --gtf annotation.gtf --fasta genome.fa \
        --read-lengths 18-34 --out-prefix out/$lib
done
```

Check the read lengths and the frame of each library as for [elongation
libraries](elongation.md#2-check-the-read-lengths). An initiation library
often has a weaker frame than its elongation library.

Then compare each initiation library with its elongation library. The share of
each read length, for the reads with the 5′ end in a CDS, and the distance
between the two distributions (0: the same; 1: no length in common):

```python
import pandas as pd


def lengths(path):
    qc = pd.read_csv(path, sep="\t").dropna(subset=["phase"])
    n = qc.groupby("length").reads.sum()
    return n / n.sum()


pair = pd.DataFrame({"initiation": lengths("out/harr_1.qc.tsv"),
                     "elongation": lengths("out/elong_1.qc.tsv")}).fillna(0)
print(pair.loc[20:36].round(3).to_string())
print("distance:", round(0.5 * (pair.initiation - pair.elongation).abs().sum(), 2))
```

The two distributions should peak at the same length and have the same
shape. If the initiation library peaks several nt longer (for example at
31-32 nt against 28 nt), it was digested less. Pair it with an elongation
library of the same digestion if you have one, or leave the pair out.

### 3. Pool offsets

A **pool** is the replicates of one condition. The libraries of a pool must
share their P-site offsets. Take them from the pool's elongation libraries,
merged:

```bash
samtools cat -o elong_ctrl.bam elong_1.bam elong_2.bam
pitstop quant --bam elong_ctrl.bam --gtf annotation.gtf --fasta genome.fa \
    --read-lengths 18-34 --out-prefix out/ctrl
```

Not from the initiation libraries: `quant` places the offsets with the reads
at both ends of the CDSs, and after the run-off an initiation library has few
reads at the stop codons. `starts` then moves each library's read lengths by
the 1-2 nt that its own digestion needs ([calibration](../starts.md#1-calibration-and-start-lengths)).

### 4. P-sites of each library

```bash
for lib in harr_1 harr_2 elong_1 elong_2; do
    pitstop orfs --bam $lib.bam --gtf annotation.gtf --fasta genome.fa \
        --read-lengths 18-34 --offsets out/ctrl.offsets.tsv --out-prefix out/$lib
done
```

No `--orfs`: this run is for the P-sites (`psites.tsv`). `--read-lengths`
must hold the lengths that `starts` can use: by default it chooses them in
24-34 nt.

### 5. Start sites

The library sheet, `libraries.tsv`: one row per replicate, each an out prefix
of step 4.

```text
pool    harringtonine    elongation
ctrl    out/harr_1       out/elong_1
ctrl    out/harr_2       out/elong_2
```

```bash
pitstop starts --libraries libraries.tsv --gtf annotation.gtf --fasta genome.fa \
    --scan --out-prefix out/starts
pitstop plot --prefix out/starts      # out/starts.lengths.pdf
```

`--scan` tests every ATG and near-cognate start codon with initiation reads,
and calls the starts with q < 0.05. `--start-codons ATG,CTG` tests fewer
codons, `--fdr` sets the cut, and `--orfs table.tsv` tests the starts of a
given table instead.

Before you use the calls, look at the read-length QC: `lengths.tsv` and the
PDF of `plot` ([Usage](../usage.md#qc-figures-plot)). A start length has most
of its P-sites in one frame and a peak at nt 0 of the annotated starts. Only
the start lengths count. If no length passes, `starts` stops; if the chosen
lengths look wrong, give them with `--start-lengths`.

### 6. Counts and frame evidence of the new ORFs

`catalogue.tsv` holds the uORFs and uoORFs of the called starts and every
annotated CDS: an ORF table for `orfs`. Count its reads in the elongation
libraries, and score them:

```bash
for lib in elong_1 elong_2; do
    pitstop orfs --bam $lib.bam --gtf annotation.gtf --fasta genome.fa \
        --read-lengths 18-34 --orfs out/starts.catalogue.tsv \
        --offsets out/ctrl.offsets.tsv --out-prefix out/$lib.catalogue
done
pitstop score --orfs-prefix out/elong_1.catalogue out/elong_2.catalogue --out-prefix out/catalogue
```

Use the elongation libraries here: frame evidence comes from ribosomes that
elongate, and in the initiation libraries most of them ran off. With several
pools, use one offsets table for all the libraries you count: [Comparing
conditions](conditions.md).

## What to use

| Question | Table |
|---|---|
| Which starts are used? | `starts.tsv` (every tested start, `q`; with `--scan`, `called`) |
| Which ORFs do the called starts make? | `start_orfs.tsv` (uORFs, uoORFs, extensions, truncations, internal ORFs, dORFs); `catalogue.tsv` (the uORFs and uoORFs, and the annotated CDSs) |
| Are the new ORFs translated in their frame? | `scores.tsv` of step 6 |
| Reads per ORF | `orfs.tsv` of step 6 |

## Notes

- **Start evidence is not translation evidence.** A call says that ribosomes
  stopped on the codon in the drug. Whether they go on to elongate in the
  ORF's frame is the question of `score`. An ORF with both has the strongest
  support.
- **Calls with care.** Leader starts with many expected reads, and starts in
  trailers, have more false calls than the rest: [Reading the
  results](../starts.md#reading-the-results).
- **Near starts.** Of two starts within 8 nt, or within 10 codons after a
  called start, one is called.
- **Multimappers.** Paralogs and pseudogenes share a start peak: `multimapped`
  in `starts.tsv` is the share of the window's reads with another alignment.
- **Drugs and treatment times.** Put libraries with another drug or another
  treatment time into their own pool. Each pool gets its own background and
  calibration, and their evidence adds up in one test.
- **Conditions.** Give each condition its own pool, in one `starts` run: the
  calls then do not depend on the condition. [Comparing
  conditions](conditions.md).
