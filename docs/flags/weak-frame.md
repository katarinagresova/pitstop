# Weak frame

**Few reads are in read lengths that keep one reading frame.** The P-site of
a read comes from its 5′ end and the offset of its length. If the reads of a
length do not keep one frame, their P-sites are not known to the nucleotide.
The steps that use the frame (`orfs`, `score`, `starts`) are then weak on this
library. The counts per CDS (`quant`) still work: a whole CDS has many
codons, and a P-site that is off by 1 nt stays in it.

## The value

The **clean share**: the share of the reads in clean read lengths
(`clean_share` in `summary.tsv`). A read length is clean if 0.7 or more of
its reads in the CDS interior have the same 5′ phase (the top-phase share,
`top_frame` in `lengths.tsv`). The CDS interior: the read is at least 15 nt
inside the CDS at both ends, so the start and stop peaks do not count. A
length needs 100 reads or more in the CDS interior and 0.3% of the reads or
more to get a class. The 20-22 nt reads of the [short state](short-state.md)
are not clean lengths.

| State | Cut |
|---|---|
| pass | clean share 0.5 or more |
| warn | 0.3 or more |
| fail | below 0.3 |
| unknown | no read length has enough reads for a class |

## Why the cut is there

- Libraries without a known problem: 0.58-0.83. The 81 libraries of one
  human cell line: 0.60-0.94.
- A harringtonine library warned (0.44). Initiation libraries often have a
  weaker frame than their elongation libraries.
- A library with untemplated 5′ nucleotides kept in the alignment failed
  (0.24). After `repair` it had 0.55.
- The libraries digested with little nuclease failed (0).
- Yeast spike-ins mostly warn (75 of 81). Their counts per CDS do not need
  the frame.

## What to do

1. Read the other flags. Each gives a cause, and some causes have a fix:
   - [Untemplated 5′ nt](untemplated.md) fails: align again with soft clips
     at the read ends, or run `repair`. The report gives the clean share
     after repair.
   - [End spike](end-spike.md) warns: trim the reads again.
   - [Under-digestion](under-digestion.md) fails or warns: no fix for this
     library.
2. If no other flag explains it, check the preprocessing: the adapter, and
   every nucleotide that the library kit adds to the insert ([Preparing the
   reads](../workflows/reads.md#before-the-alignment)).
3. Use the clean lengths, if there are some: "Next steps" in the report gives
   them as `--read-lengths`.
