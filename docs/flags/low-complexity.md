# Low complexity

**Few of the reads are distinct footprints.** A footprint is a transcript, a
5′ end and a length. Many reads of the same footprint come from PCR copies of
one molecule, or from many ribosomes on the same codon: at pauses, at start
codons, and on highly expressed genes. Without UMIs, the two look the same.

## The value

The **distinct share at 1 million reads** (`distinct_share_at_depth` in
`summary.tsv`): the expected number of distinct footprints in a subsample of
1 million reads, divided by 1 million. It is exact, from the reads per
footprint: the mean over all subsamples of that size, not one random
subsample.

The distinct share of all reads (`distinct_share`) falls as the depth goes
up, so it cannot compare libraries of different depths. The flag compares
them at one depth.

| State | Cut |
|---|---|
| pass | 0.4 or more |
| warn | below 0.4 |
| unknown | fewer than 1 million reads |

There is no fail.

## Why the cut is there

- The lowest value of the libraries of the check with 1 million reads or more
  was 0.47 (a harringtonine library). No library warns.
- The distinct share of all reads warned on three libraries with 12-24
  million reads, but at the same depth they were in the range of the
  libraries that passed.
- Libraries whose PCR copies were removed with UMIs had the same distinct
  share as the other libraries. Thus the value does not see PCR copies well,
  and a low value can be biology. The flag only warns.

## What to do

- Check the duplicate removal of the preprocessing. With UMIs, remove PCR
  copies with them. Without UMIs, do not collapse identical reads: many are
  different ribosomes ([Preparing the
  reads](../workflows/reads.md#before-the-alignment)).
- Look at the curve "Distinct footprints against the reads" in the report. If
  it still rises at full depth, more sequencing adds new footprints. If it is
  flat at a low distinct share, the library is exhausted: a new library, with
  more input, gives more.
