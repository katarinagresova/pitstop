# Flags of the read QC

`pitstop inspect` checks one library and gives each check a flag. The flags
are in the report (`<prefix>.inspect.html`, the table under "Verdict and next
steps"), in `summary.tsv` (the columns `flag_<name>`), and in MultiQC (the
number of flags that fail or warn per library). Each flag of the report links
to its page here.

## The states

| State | Meaning |
|---|---|
| pass | The value is on the good side of the cut. |
| warn | Look at the library before you use it. The page of the flag says what to do. |
| fail | The library has the problem. The page of the flag says what fixes it, if anything does. |
| unknown | Too few reads (or no read sequences in the BAM) to measure the value. Unknown is not pass. |
| info | A finding with no cut: it does not mean that the library is bad, but it can change what you do next. |

A fail does not change the exit code of `inspect`: a flag does not stop a
workflow.

## The flags

| Flag | It checks | Can be |
|---|---|---|
| [Weak frame](weak-frame.md) | the share of the reads in read lengths that keep one reading frame | pass, warn, fail, unknown |
| [Long reads](long-reads.md) | the share of the reads of 32 nt or more | pass, warn, fail |
| [Under-digestion](under-digestion.md) | how far the 3′ ends of the reads spread, for reads with the same 5′ end | pass, warn, fail, unknown |
| [Untemplated 5′ nt](untemplated.md) | an extra nucleotide at the 5′ end of the reads that is not in the transcript | pass, info, fail, unknown |
| [End spike](end-spike.md) | mismatches at fixed positions of the read ends: left-over UMI or adapter nucleotides | pass, warn, unknown |
| [Mixed lengths](mixed-lengths.md) | read lengths whose 5′ ends split between two nucleotides at the start codons | pass, info |
| [Short state](short-state.md) | whether the 20-22 nt reads are ribosomes or RNA fragments | info, warn, unknown |
| [Low complexity](low-complexity.md) | the share of distinct footprints at a fixed depth | pass, warn, unknown |

## How the cuts were set

The cuts are in the code (`src/pitstop/inspect.py`), not in a file of
options. They were set on two groups of libraries:

- 20 libraries of human and mouse cells, with known problems and without:
  libraries digested with little nuclease, libraries with untemplated 5′
  nucleotides kept in the alignment, harringtonine and lactimidomycin
  libraries, and well-made libraries of several protocols;
- 81 libraries of one human cell line, made with one protocol and without a
  known problem, and their yeast spike-ins.

Each page gives the values of the two groups, and why the cut is where it is.
A cut was set so that each known problem fails or warns, and so that the
libraries without a problem pass.
