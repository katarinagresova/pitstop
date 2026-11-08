# Short state

**Are the 20-22 nt reads ribosomes?** Ribosomes with an open A site protect a
shorter footprint, about 21 nt, than the 28-30 nt of the other states (Lareau
et al. 2014; Wu et al. 2019). Their 5′ end is that of a 28 nt read, and the
3′ end is short. They are a state of the ribosome, not damaged footprints:
do not repair them, and do not count them as 28 nt reads. But RNA fragments
of other origin are also often this short.

## The value

The share of the P-sites of the 20-22 nt reads that fall in the CDS
(`short_cds_share` in `summary.tsv`). Ribosomes are mostly in the CDS. RNA
fragments are not.

| State | Cut |
|---|---|
| info | half of their P-sites or more are in the CDS: the short state |
| warn | less than half: RNA fragments |
| unknown | fewer than 10,000 reads of 20-22 nt (`short_reads`) |

A read length of 20-22 nt gets the class "short state" in `lengths.tsv` if
its 5′ ends peak 12 nt before the start codons (20 reads or more there), and
half of its P-sites or more are in the CDS.

## Why the cut is there

- With 10,000 reads, the standard error of the CDS share is 0.005 or less. A
  first cut of 50,000 reads made a fifth of the libraries of one human cell
  line unknown.
- The CDS share was 0.64-0.97 in all libraries of the check. No library
  warns.
- A limit: with a GTF without leaders and trailers (as for many yeast
  spike-ins), every P-site is in the CDS. The share is then 1, and the flag is
  info.

## What to do

- **Info.** The 20-22 nt reads are ribosomes. They have their own offset
  (`offsets.tsv` of `quant` gives one per length) and often a frame of their
  own. They are not clean lengths, so they are not in the suggested read
  lengths. To count them, widen `--read-lengths` (one range: it then also
  takes the lengths between), and check their frame in the report.
- **Warn.** Leave the 20-22 nt reads out. Look at the regions per length in
  the report: if many P-sites of all lengths are outside the CDS, the library
  has contamination.

## References

- Lareau LF, Hite DH, Hogan GJ, Brown PO (2014). Distinct stages of the
  translation elongation cycle revealed by sequencing ribosome-protected mRNA
  fragments. eLife 3:e01257.
  [doi:10.7554/eLife.01257](https://doi.org/10.7554/eLife.01257)
- Wu CC-C, Zinshteyn B, Wehner KA, Green R (2019). High-resolution ribosome
  profiling defines discrete ribosome elongation states and translational
  regulation during cellular stress. Molecular Cell 73:959-970.
  [doi:10.1016/j.molcel.2018.12.009](https://doi.org/10.1016/j.molcel.2018.12.009)
