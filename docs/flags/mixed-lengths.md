# Mixed lengths

**At the start codons, the 5′ ends of a read length split between two
neighbouring nucleotides.** The reads of ribosomes on the start codon then
form two groups, 1 nt apart. A length has one offset, so it puts the smaller
group's P-sites 1 nt off. At a start codon this matters: a start is tested in
a window of a few nucleotides. In the CDS, the frame of the length shows how
many of its reads are off.

## The value

For each read length, the commonest 5′ end at the annotated start codons
(`start_peak` in `lengths.tsv`) and the larger of its two neighbours, as a
share of the peak (`start_beside`). A length is mixed if this share is 0.5 or
more and the peak has 20 reads or more (`mixed` in `lengths.tsv`; the list in
`mixed_lengths` of `summary.tsv`).

| State | Cut |
|---|---|
| info | one length or more is mixed |
| pass | no length is mixed |

There is no fail. A mixed length can still be a clean length: its class comes
from its frame in the CDS interior, as for every length, because `quant`,
`orfs` and `score` use that frame.

## Why the cut is there

The cut is the same as the cut of `starts` for its start lengths: one group
of 5′ ends with half of the other's reads or more. Which of the two
nucleotides holds the peak can change between libraries of one experiment, so
the offset of a mixed length is not stable at the starts.

## What to do

- **Start sites.** Leave the mixed lengths out of the start lengths of
  `pitstop starts`. `starts` checks this itself in each harringtonine library
  and drops a length that is mixed in one of them. `--start-lengths`
  overrides its choice.
- **Counts and frame evidence.** Use the clean lengths, mixed or not. Their
  frame in the CDS decides.
- If the flag [untemplated 5′ nt](untemplated.md) fails, the added nucleotide
  can make the second group. Inspect the repaired BAM, or align again with
  soft clips.

The heat map "5′ ends around the start codon" of the report shows the two
groups, and marks each mixed length.
