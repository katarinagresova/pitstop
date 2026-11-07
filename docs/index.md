# pitstop

A Python Ribo-seq toolkit, deterministic from end to end: P-site offset
estimation and CDS-level quantification from transcriptome alignments
(`pitstop quant`), reading-frame evidence for uORFs and other ORFs besides
the CDS (`pitstop orfs` / `pitstop score`), and translation start sites from
harringtonine libraries, with a uORF catalogue from them (`pitstop starts`).

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
