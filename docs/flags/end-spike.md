# End spike

**Mismatches pile up at fixed positions of the read ends.** Nucleotides that
the trimming left on the reads, such as UMIs, random nucleotides of the
library kit or adapter nucleotides, are not in the transcript. At their
positions about 3 of 4 reads have a mismatch. A nucleotide left at the 5′ end
moves the footprint's 5′ end; one left at the 3′ end makes the footprint
longer.

## The value

The mismatch rate at positions 2-3 from the 5′ end (`mismatch5_inner` in
`summary.tsv`, the higher of the two) and at positions 0-1 from the 3′ end
(`mismatch3`), each against the background of its end: the mismatch rate at
positions 4-7 from that end. Position 0 of the 5′ end is the
[untemplated nucleotide](untemplated.md), and has its own flag.

| State | Cut |
|---|---|
| warn | a rate is more than 3 times its background and more than 10% |
| pass | else |
| unknown | the BAM has no read sequences |

There is no fail: a spike points to the preprocessing, not to the library.

## Why the cut is there

- Left-over nucleotides give about 75% at their positions.
- The highest rate in all the libraries of the check was 5%. No library warns.
- A first cut of 2% warned on three libraries with 3-5% at position 2 of the
  5′ end. A mismatch at position 2 does not move the 5′ end, so the P-site and
  the frame do not change, and no action follows.

## What to do

Check the trimming against the read structure of the library kit ([Preparing
the reads](../workflows/reads.md#before-the-alignment)). Remove the fixed
nucleotides before the alignment, for example with cutadapt `-u`, and align
again. `repair` does not trim a fixed number of nucleotides.
