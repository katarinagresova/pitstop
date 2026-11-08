# Under-digestion

**The nuclease did not trim the footprints to the ribosome.** A ribosome
protects about 28-30 nt. With enough nuclease, the reads of ribosomes on the
same codon have the same ends. With too little, the ends of the reads are at
many places, the reads are longer, and their frame is weak.

## The value

The **3′ end spread** (`end3_spread` in `summary.tsv`, in nt). At the
annotated start codons, `inspect` finds the commonest 5′ end of the reads of
23-49 nt (`end3_peak`, nt -20 to -4 from the start codon). For the reads with
their 5′ end there (`end3_reads`), it takes the standard deviation of their
lengths. The 5′ end is fixed, so the length gives the position of the 3′ end.

The reads of 22 nt or less are not in the spread. They are the
[short state](short-state.md): the 5′ end of a 28 nt read and a short 3′ end,
not a trim problem.

| State | Cut |
|---|---|
| pass | below 1.4 nt |
| warn | below 1.7 nt |
| fail | 1.7 nt or more |
| unknown | fewer than 1,000 reads at the peak |

A GTF without leaders gives unknown: no read starts before the start codon.

## Why the cut is there

| Libraries | 3′ end spread (nt) |
|---|---|
| without a known problem | 0.89-1.30 |
| 81 libraries of one human cell line | 0.97-1.17 |
| digested with a medium amount of nuclease | 1.56 (warn) |
| digested with little nuclease (elongation, harringtonine, lactimidomycin) | 1.82-2.46 (fail) |

The [long reads](long-reads.md) do not find the library of medium digestion:
it has 3% long reads, and its read lengths are almost those of a clean
library. Its frame is weak, and its 3′ ends spread.

The spread of the 5′ ends at the stop codons does not separate the libraries:
2.15 nt in the library of medium digestion, 2.02 nt in a library without a
problem.

## What to do

- **`repair` does not help.** No rule for each read finds the true ends of a
  footprint that the nuclease did not trim.
- Use the clean lengths, if there are some. Often there are none.
- **Pairs.** An initiation library and its elongation library must have the
  same digestion: `starts` compares them at each read length. Compare the
  read lengths of the two ([Initiation and elongation
  libraries](../workflows/initiation.md#2-qc-of-each-library-and-of-each-pair)),
  and pair the library only with a library of the same digestion.
- For the next libraries: digest more.
