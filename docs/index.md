---
hide:
  - toc
---

# pitstop

**Ribo-seq reads in, translation out.**

![pitstop at a glance: reading frames drawn as the three lanes of a race track. Ribosomes keep to their ORF's lane, so a uORF whose reads stay in its frame is translated (pitstop orfs + score); harringtonine holds ribosomes on start codons (pitstop starts); reads are counted per coding sequence (pitstop quant); each library is checked first, for the read lengths that keep their lane (pitstop inspect); and extra 5' nucleotides are clipped from its reads (pitstop repair).](assets/pitstop-abstract.svg)

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

![One transcript with a uORF (CUG start) and a CDS (AUG start), and two read tracks above it. Harringtonine reads peak at both start codons, above what the Ribo-seq reads predict (pitstop starts). Ribo-seq P-sites are coloured by reading frame: on the uORF most keep the uORF's frame (pitstop orfs + score), on the CDS the CDS's frame (pitstop quant). A bar chart of the read lengths of a library marks the lengths whose reads keep one frame (pitstop inspect). Reads with an extra 5' nucleotide are shown before and after it is clipped (pitstop repair).](assets/pitstop-transcript.svg)

## Typical workflows

Which commands you run, and in which order, depends on your libraries. Each
workflow page describes the experiment, the inputs, the commands with
examples, and what to watch for.

| Your libraries | Workflow | Commands |
|---|---|---|
| Ribo-seq of elongating ribosomes | [Elongation libraries](workflows/elongation.md) | `inspect`; `quant`; `orfs` and `score` on your ORF table |
| Harringtonine or lactimidomycin libraries, each with a matched elongation library | [Initiation and elongation libraries](workflows/initiation.md) | `inspect`; `quant`, `orfs`, `starts`; then `orfs` and `score` on the new ORFs |
| Either of these, in two or more conditions | [Comparing conditions](workflows/conditions.md) | the same, with one ORF table and one offsets table for all libraries; then your own tests on the counts |

Every workflow starts from transcriptome alignments: [Preparing the
reads](workflows/reads.md) covers the trimming, the contaminants, and the
alignment with STAR or bowtie. To count ORFs from an ORF caller or from genome
coordinates: [ORF tables from other tools](workflows/orf-tables.md).

Before the first command, `pitstop inspect` checks each library and names the
read lengths that keep their reading frame, for `--read-lengths`. Each of its
[flags](flags/index.md) has a page: what it means, its cut, and what to do. If
it finds extra nucleotides at the 5′ end of the reads, let the aligner clip
them the next time (STAR with `--quantTranscriptomeSAMoutput BanSingleEnd`, or
bowtie2 `--local`). For the BAM you have, `pitstop repair` clips them and
writes a new BAM for the other commands.

## Why

Ribo-seq quantification needs the isoform-length term accounted for, a fit
that actually converges from a fixed, deterministic start, and P-site offsets
estimated per read length from the data rather than assumed. Counting evidence for a
uORF additionally needs a reading frame that isn't thrown away before it's
checked, and a null that accounts for the host CDS's own off-frame noise.
pitstop does all of this, benchmarked against real runs rather than
simulation alone.

## Get started

- **[Usage](usage.md)**: install, the commands, their inputs, and every
  output column.
- **[Flags of the read QC](flags/index.md)**: what each flag of `inspect`
  means, its cut and why, and what to do.
- **[Typical workflows](#typical-workflows)**: which commands to run, in
  which order, from [preparing the reads](workflows/reads.md) to the tables
  for your tests.
- **[Method: quant](method.md)**: how offsets, P-site assignment, and
  quantification work, with interactive examples.
- **[Method: ORF counting](orfs.md)**: ORF types, outside components, the
  frame-weighted EM, and the codon lead score.
- **[Method: start sites](starts.md)**: harringtonine against elongation,
  the background, the test, the scan and the catalogue.
- **Using `orfs` and `score`**: a [typical workflow](orfs.md#typical-workflow)
  for two conditions, [how to read the scores](orfs.md#reading-the-results),
  and [what they give on real data](orfs.md#validation-on-real-data).

## How it works, in short

1. **Annotation.** One CDS per transcript, in transcript coordinates, without
   the stop codon. It is dropped if its length is not a multiple of 3 or no
   stop follows it. Non-ATG starts are kept.
2. **Reads.** Mapped, forward-strand, primary and secondary alignments. A 5'
   soft clip counts as untemplated nucleotides: the 5' end is the first
   aligned base and the length excludes the clip. Alignments with
   indels/splices are dropped. Each read name counts once.
3. **Offsets.** CDS inclusion: for each read length, the three consecutive
   offsets (one per phase) that put the most reads' P-sites inside CDSs.
4. **Assignment.** A read is compatible with a CDS if its P-site falls inside
   it.
5. **EM.** Expectation-maximisation over equivalence classes, with
   P(read | CDS) = 1/length, accelerated with SQUAREM. It starts
   deterministically and runs until an EM step moves no expected count by
   more than 1e-3 reads. It fails if it does not get there.

`orfs` and `score` add:

1. **ORF table.** uORFs, uoORFs and other ORFs besides the CDS, in transcript
   coordinates; dropped if not a multiple of 3 or no stop follows. Positions
   in no ORF form per-transcript leader/trailer/whole-transcript components.
   Alignments to BAM references that are not in the GTF are dropped.
2. **Frame-true P-sites.** One offset per read length (no phase snapping), so
   every P-site keeps its reading frame, plus a per-length frame profile from
   CDS interiors.
3. **EM with a frame term.** As above, but P(read | ORF) also depends on how
   often the read length lands in that ORF's frame, so reads shared with the
   CDS split by frame as well as by density.
4. **Codon lead score.** Per ORF, each codon (its start excluded) votes
   whether more of its reads are in the ORF's frame than in either other
   frame; the score is how unlikely that many leads are if the ORF were not
   translated, against a null that accounts for a host CDS's own frame bleed.

`starts` adds:

1. **Start lengths.** The read lengths whose harringtonine P-sites peak on
   the first nt of the annotated start codons in every library.
2. **Background.** The elongation P-sites times a harringtonine factor per
   region (leader, CDS, trailer) and distance from the CDS start, and a
   scale per transcript, fit on null codons: positions with no start codon
   within 8 nt.
3. **The test.** Per start, the harringtonine P-sites in nt -1 to +2 of its
   codon against the expected, negative binomial, with a gamma prior on the
   window's elongation rate and an overdispersion calibrated on the null
   codons. A scan calibrates it twice: the second time without the null
   codons in the window and codons 1-10 of the starts of a first call.
4. **The scan.** Start codons with the smallest p first; a start is called
   if it passes BH and no called start within 8 nt or 10 codons upstream
   stops it. Each called start gives an ORF to its first in-frame stop.

Full method, with worked examples: [Method: quant](method.md),
[ORF counting](orfs.md) and [start sites](starts.md).

## License

Apache License 2.0, see [LICENSE](https://github.com/katarinagresova/pitstop/blob/main/LICENSE)
and [NOTICE](https://github.com/katarinagresova/pitstop/blob/main/NOTICE).
