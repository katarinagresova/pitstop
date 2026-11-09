# Elongation libraries

The workflow for Ribo-seq libraries of elongating ribosomes, with no
initiation libraries.

## The experiment

- **Libraries.** Standard Ribo-seq: footprints of elongating ribosomes, after
  nuclease digestion, with or without cycloheximide. Replicates of each
  condition. Often matched RNA-seq too, for translation efficiency; pitstop
  does not count RNA-seq reads.
- **Questions.** How much is each coding sequence (CDS) translated? Are given
  ORFs, such as upstream ORFs (uORFs), translated: do their reads keep the
  ORF's reading frame?
- **Not covered.** pitstop does not find new ORFs in elongation libraries
  alone. Give it an ORF table: from an ORF caller, from an annotation such as a
  list of Ribo-seq ORFs, or your own candidates. To find start sites with
  pitstop, add initiation libraries: [Initiation and elongation
  libraries](initiation.md).

![The elongation workflow, one card per step. Your BAMs go into step 1, pitstop inspect on each library, which writes inspect.html and summary.tsv. Step 2: check the read lengths in the reports and in MultiQC, choose --read-lengths, and repair a BAM if its report says so. Step 3: pitstop quant on each library writes quant.tsv. Step 4: pitstop quant on all libraries merged writes offsets.tsv. Step 5: pitstop orfs --orfs on each library, with the BAMs, your ORF table and offsets.tsv, writes orfs.tsv. Step 6: pitstop score on all libraries together writes scores.tsv.](../img/workflow-elongation.svg)

## Inputs

- A transcriptome BAM per library: [Preparing the reads](reads.md).
- The GTF and the genome FASTA of the alignment.
- For ORFs: a TSV with the columns `ORF_id Name start end`, in transcript
  coordinates of the same GTF.

```text
ORF_id                  Name                 start   end
ENST00000000001_uORF1   ENST00000000001      40      85
ENST00000000001_uORF2   ENST00000000001      112     130
```

`Name` is the transcript id. `start` is the first nt of the start codon
(0-based), and `end` is one past the last sense codon: the stop codon is not
part of the ORF. The annotated CDSs are added to the table. To make the table
from the output of an ORF caller, or from ORFs in genome coordinates: [ORF
tables from other tools](orf-tables.md).

## Steps

### 1. QC of each library

```bash
for lib in ctrl_1 ctrl_2 ctrl_3; do
    pitstop inspect --bam $lib.bam --gtf annotation.gtf --fasta genome.fa --out-prefix out/$lib
done
multiqc out/
```

`out/ctrl_1.inspect.html` is the report of one library: its flags, the read
lengths to use, and what to do next ([Usage](../usage.md#read-qc-inspect)).
MultiQC shows all libraries side by side: the read lengths, and the
top-phase share of each length.

If the report says so (the flag [untemplated 5′
nt](../flags/untemplated.md) fails), align the library again with soft clips
at the read ends ([Preparing the reads](reads.md#star-recommended)). For a
BAM that you cannot align again, repair it, and inspect the new BAM:

```bash
pitstop repair --bam ctrl_1.bam --gtf annotation.gtf --fasta genome.fa --out-prefix out/ctrl_1.repaired
pitstop inspect --bam out/ctrl_1.repaired.bam --gtf annotation.gtf --fasta genome.fa \
    --out-prefix out/ctrl_1.repaired
```

Then use `out/ctrl_1.repaired.bam` in all next steps, and take the read
lengths from its report.

### 2. Check the read lengths

"Next steps" in each report gives the read lengths to use, as
`--read-lengths` (`read_lengths` in `summary.tsv`): the clean lengths, those
with 0.7 or more of their reads in the CDS interior in one phase.

- **One range for all libraries** that you count together: the lengths that
  are clean in all of them. The plot of the top-phase share in MultiQC shows
  them. The examples below use 26-30 nt:

    ```bash
    LENGTHS=26-30
    ```

- **A top-phase share near 1/3 means no frame**: the P-site of these reads is
  not known. Leave such lengths out of `--read-lengths`.
- **No length with a clear frame** points to the library, not to pitstop.
  The flags of the report name the cause, if it is a known one: an
  untemplated 5′ nucleotide kept by the aligner, nucleotides left at a read
  end, or too little nuclease digestion ([Flags](../flags/index.md)).
- **Short footprints** (about 20-22 nt) come from ribosomes in another state
  than the 28-30 nt footprints ([short state](../flags/short-state.md)). They
  can have a frame of their own and an offset of their own: `offsets.tsv`
  gives one per length.

### 3. Counts per CDS

```bash
for lib in ctrl_1 ctrl_2 ctrl_3; do
    pitstop quant --bam $lib.bam --gtf annotation.gtf --fasta genome.fa \
        --read-lengths $LENGTHS --out-prefix out/$lib
done
```

`quant.tsv` has the reads per CDS (`NumReads`), with the reads of shared
sequence split between isoforms and paralogs.

### 4. One offsets table for all libraries

```bash
samtools cat -o all.bam ctrl_1.bam ctrl_2.bam ctrl_3.bam
pitstop quant --bam all.bam --gtf annotation.gtf --fasta genome.fa \
    --read-lengths $LENGTHS --out-prefix out/all
```

The offsets of short read lengths can differ between libraries by a codon or
more. For a short ORF, such a difference moves reads across its ends. The
counts are only comparable between libraries that share one offsets table:
here the table of all libraries merged (the table of the deepest library also
works).

### 5. Counts per ORF

```bash
for lib in ctrl_1 ctrl_2 ctrl_3; do
    pitstop orfs --bam $lib.bam --gtf annotation.gtf --fasta genome.fa \
        --read-lengths $LENGTHS --orfs orfs.tsv --offsets out/all.offsets.tsv \
        --out-prefix out/$lib
done
```

`orfs.tsv` has the reads of each ORF and of each annotated CDS, and of the
parts of each transcript in no ORF (`<transcript>:5UTR`, `:3UTR`). Every
P-site keeps its reading frame.

### 6. Frame evidence

```bash
pitstop score --orfs-prefix out/ctrl_1 out/ctrl_2 out/ctrl_3 --out-prefix out/all
```

`scores.tsv` has a p- and q-value per ORF: do its codons have more reads in
the ORF's frame than in the other two? `score` sums the libraries, so score
all replicates together. It uses the read lengths with a sharp frame in every
library (frame 0 share of 0.9 or more in `frames.tsv`). If no length passes
in every library, `score` stops: give `--frame-lengths` with the sharpest
lengths.

## What to use

| Question | Table | Column |
|---|---|---|
| CDS translation, differential translation | `quant.tsv` | `NumReads` |
| Reads per ORF, per 5′ UTR and 3′ UTR | `orfs.tsv` | `NumReads` |
| Is the ORF translated in its frame? | `scores.tsv` | `q`, and the group: called, not called, cannot be called |

## Notes

- **For CDSs, use `quant`.** `orfs.tsv` has the CDSs too, but the ORFs and the
  5′ UTRs and 3′ UTRs take some of their reads.
- **"Not called" is not "not translated".** Most short ORFs have too few
  codons with reads to pass at any depth. `scores.tsv` sorts every ORF into
  called, not called and cannot be called: [Reading the
  results](../orfs.md#reading-the-results).
- **A call is about the frame, not the start.** Two ORFs with the same stop
  share their codons, and the score cannot tell which start is used. Sum the
  counts of such nested ORFs: [Known limits](../orfs.md#known-limits).
- **Spike-ins.** Count spike-in reads with `quant` on the spike-in's GTF and
  FASTA, with the offsets of the main library (`--offsets`), for size factors:
  [Usage](../usage.md#counts-per-cds-quant).
- **Translation efficiency** needs RNA-seq counts from another tool, at the
  same level as pitstop's (transcripts or genes). Two or more conditions:
  [Comparing conditions](conditions.md).
- **Reruns.** pitstop uses no random numbers: the same BAM and options give the
  same outputs, byte for byte.
