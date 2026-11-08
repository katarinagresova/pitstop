# Long reads

**Many reads are 32 nt or longer.** A ribosome protects about 28-30 nt of the
mRNA from the nuclease. Many longer reads mean that the nuclease did not trim
the footprints to the ribosome, or that the size selection kept longer
fragments.

## The value

The share of the reads of 32 nt or more (`long_share` in `summary.tsv`). The
length is the footprint length: the read length without a 5′ soft clip.

| State | Cut |
|---|---|
| pass | below 10% |
| warn | below 25% |
| fail | 25% or more |

## Why the cut is there

- Libraries without a known problem: 9% or less. The 81 libraries of one
  human cell line: 1.5% or less.
- Libraries digested with little nuclease: 32-55% (fail).
- Some libraries have long reads with a 3′ end at few places: 15-21% (warn)
  and 38% (fail). Their 3′ ends do not spread, so they are not
  [under-digested](under-digestion.md). Look at their frame: the long reads
  of these libraries had no frame.

The share of long reads alone does not find every under-digested library. A
library digested with a medium amount of nuclease had 3% long reads, and a
length distribution like that of a clean library. The flag
[under-digestion](under-digestion.md) finds it.

## What to do

- Read the flag [under-digestion](under-digestion.md). If it fails or warns,
  the nuclease is the cause.
- The long reads with a frame are clean lengths, and you can use them.
  "Next steps" in the report gives the clean lengths.
- For the next libraries: more nuclease, or a narrower size selection.
- Pair an initiation library only with an elongation library of the same
  digestion ([Initiation and elongation
  libraries](../workflows/initiation.md#2-qc-of-each-library-and-of-each-pair)).
