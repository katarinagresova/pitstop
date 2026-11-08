# Untemplated 5′ nucleotide

**Many reads start with a nucleotide that is not in the transcript.** In
many library protocols, the reverse transcriptase adds a nucleotide that is
not in the RNA, and it becomes the first nucleotide of the read. If the
alignment keeps it, the footprint starts 1 nt too early, its P-site is 1 nt
off, and the frame of its length is blurred.

## The value

The **5′ mismatch** (`mismatch5` in `summary.tsv`): the share of the reads
whose first nucleotide does not match the transcript. The first nucleotide
is that of the read, a soft clip included. The background
(`mismatch5_background`) is the mismatch rate at positions 4-7 from the 5′
end. `clip5_share` is the share of the reads with a 5′ soft clip, and
`mismatch5_T` the share of T among the mismatched first nucleotides.

| State | Cut |
|---|---|
| fail | the 5′ mismatch is 5% or more and 5 times the background or more, and the aligner did not clip it |
| info | the same, but the aligner clipped it: the 5′ soft clips are half of the 5′ mismatch or more |
| pass | else |
| unknown | the BAM has no read sequences |

## Why the cut is there

- BAMs without the problem: 0.3-1.3%.
- Libraries whose alignments kept the nucleotide: 19-21% (fail). On one of
  them, `repair` took the clean share from 0.24 to 0.55.
- Libraries aligned with soft clips: the 5′ mismatch is high, but most of
  these reads have a 5′ soft clip (info).

## What to do

**Fail.** Two fixes, the first one is better:

1. **For the next alignment:** let the aligner clip the nucleotide. An
   aligner that soft-clips the read ends removes it, and every pitstop
   command leaves a 5′ soft clip out of the footprint.
   - STAR in local mode (its default), with `--quantTranscriptomeSAMoutput
     BanSingleEnd` for its transcriptome BAM. STAR's default for this option
     extends the soft clips again. In STAR 2.7.10b and older, the option is
     `--quantTranscriptomeBan Singleend`. See [Preparing the
     reads](../workflows/reads.md#star-recommended).
   - bowtie2 with `--local`.
   - bowtie (version 1) cannot soft-clip.
2. **For the BAM that you have** (from another pipeline, or from bowtie):
   `pitstop repair`. It makes the mismatched first 1-2 nucleotides of each
   read soft clips, and writes a new BAM. Then run `inspect` on the new BAM,
   and take the read lengths from that report. The report gives the clean
   share and the clean lengths that it expects after repair.

```bash
pitstop repair --bam sample.bam --gtf annotation.gtf --fasta genome.fa --out-prefix out/sample.repaired
pitstop inspect --bam out/sample.repaired.bam --gtf annotation.gtf --fasta genome.fa --out-prefix out/sample.repaired
```

No trim before the alignment can do this: the extra nucleotide is on some
reads only, and only the transcript shows which ones.

**What no rule can fix.** An added nucleotide that matches the transcript by
chance looks like a true footprint end. An added nucleotide matches the
transcript about one time in four, so these are about a third of the ones
that are seen. The report gives the estimate (`untemplated_hidden` in `summary.tsv`).

**Info.** Nothing to do: the aligner clipped the nucleotide, and pitstop
leaves a 5′ soft clip out of the footprint.
