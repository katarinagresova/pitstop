# Start sites from harringtonine

`pitstop starts` finds translation start sites from harringtonine (or
lactimidomycin) Ribo-seq libraries. It compares each harringtonine library
with an elongation library of the same condition and replicate, tests given
starts or every start codon of the transcriptome, calls starts, and writes
their ORFs and a catalogue of uORFs and uoORFs.

## Background

**Harringtonine stops ribosomes at the start codon.** Harringtonine binds the
large ribosomal subunit and blocks the first steps of elongation, but not
initiation ([Fresno et al. 1977](#references)). Ribosomes that start after the
drug is added stay on the start codon; ribosomes that were already elongating
run off the mRNA. After a few minutes, footprints pile up on start codons
([Ingolia et al. 2011](#references)). Lactimidomycin acts in the same way
([Lee et al. 2012](#references)).

**A peak in a harringtonine library is not yet a start.** Three things also
make harringtonine footprints:

- **Run-off is not complete.** Ribosomes far into a long CDS have not run off
  yet, so the harringtonine / elongation ratio grows along the CDS. Leaders
  have a higher ratio than CDSs.
- **Codons 1-10 after a start are raised.** Ribosomes that made a few peptide
  bonds before the drug acted stay there: 4-8x at codons 1-5, about 3x at
  codons 6-10, against codons 11-20.
- **Pauses.** A codon where elongating ribosomes stall has a peak in every
  library, with or without the drug.

Thus `starts` does not look for peaks in the harringtonine library alone. It
asks if a start codon has more harringtonine P-sites than the elongation
P-sites at the same position predict.

## Inputs

`starts` reads the outputs of `pitstop orfs` runs
([ORF counting](orfs.md)), one per library, and a library sheet:

```text
pool        harringtonine          elongation
control     out/harr_control_rep1  out/elong_control_rep1
control     out/harr_control_rep2  out/elong_control_rep2
treated     out/harr_treated_rep1  out/elong_treated_rep1
treated     out/harr_treated_rep2  out/elong_treated_rep2
```

Each row is one replicate: a harringtonine library and its matched elongation
library, given as out prefixes of `orfs`. A **pool** is the replicates of one
condition. Several pools can be given; each keeps its own background, and
their evidence adds up.

The libraries of a pool must have the same P-site offsets, or their P-sites
are not comparable. Give each `orfs` run of a pool the offsets of
`pitstop quant` on the pool's elongation BAMs, merged into one
(`--offsets`). Not the harringtonine BAMs: `quant` places the offsets with
the reads at both ends of the CDSs, and after the run-off a harringtonine
library has few reads at the stop codons. Its offsets can then be a codon
too long. `starts` then moves each library's read lengths by the 1-2 nt
that its own digestion needs ([calibration](#1-calibration-and-start-lengths)).
The libraries also need the same transcripts: give them the same GTF.
Transcripts that are not in every `orfs` run are left out.

```bash
samtools cat -o elong_control.bam elong_control_rep1.bam elong_control_rep2.bam
pitstop quant --bam elong_control.bam --gtf annotation.gtf --fasta genome.fa \
    --read-lengths 18-30 --out-prefix out/control
for lib in harr_control_rep1 harr_control_rep2 elong_control_rep1 elong_control_rep2; do
    pitstop orfs --bam $lib.bam --gtf annotation.gtf --fasta genome.fa \
        --read-lengths 18-30 --offsets out/control.offsets.tsv --out-prefix out/$lib
done
# the same for the treated pool, with its own offsets
pitstop starts --libraries libraries.tsv --gtf annotation.gtf --fasta genome.fa \
    --scan --out-prefix out/starts
```

## Method

### 1. Calibration and start lengths

One offset table serves every library of a pool, but the libraries are not
digested alike. In a harringtonine library, most 5' ends of a read length can
be 1 nt longer than in the elongation library. Its P-sites then sit 1 nt
before the codon: at the start codon and at every pause. `starts` therefore
calibrates each library, per read length, with two anchors:

- the **start peak**: the nt (-15 to +32) where its P-sites peak around the
  annotated start codons;
- the **frame**: the frame of the CDS that holds most of its P-sites well
  inside the CDSs (15 nt or more from both ends).

A length is **calibrated** when its start peak is 2 nt or less from nt 0
(with 20 P-sites or more at the peak), and the frame does not clearly put
the P-sites elsewhere. "Clearly" means at least 100 P-sites inside the CDSs,
and a dominant frame with 1.5 times the P-sites of the frame that the move
puts onto the codon. Its P-sites then move by minus the start peak. The two
anchors need each other: the frame cannot tell a codon (a peak at -1 and a
peak at +2 are both frame 2), and a start peak alone could be a property of
the start only, not of every P-site of the length. A length that is not
calibrated moves by its frame if that is clear, else it stays.

A length is **mixed** when a nt beside its start peak has half the peak's
P-sites or more (with 20 P-sites or more at the peak). Its 5' ends then split
between two offsets, often a long read with 1 extra nt at the 5' end or
without it. A length has one offset, so the P-sites of the other half sit
1 nt off. At a start this costs little, because a start's window is 4 nt
(step 2). But it puts the other half 1 nt off at every other codon too, and
whether the peak or its neighbour wins can change from library to library.

The **start lengths** come from a QC of each read length in each library. A
start length:

- is in the **window**, 24-34 nt by default (`--length-window`). An 80S
  ribosome protects about 28 nt; shorter reads are mostly RNA fragments and
  other ribosome states, with a weak frame and harringtonine peaks that the
  elongation libraries do not have;
- is calibrated, and not mixed, in every harringtonine library;
- has a **clear frame** in every library: 100 P-sites or more well inside the
  CDSs, and 40% or more of them in one frame (a third is no frame).

Of the runs of consecutive lengths that pass, the run with the most P-sites
wins. `--start-lengths` sets the start lengths instead, mixed lengths too.
Only their P-sites count, one per alignment, in both kinds of library.

`lengths.tsv` has the QC, per library and length: the P-sites (all, their
share of the library's, and per region), the P-sites per frame inside the
CDSs, the start peak and the larger nt beside it, whether the length is mixed
and calibrated, its move, its offset (the pool's offset plus the move), and
the verdict (`pass`, or why not). `profiles.tsv` has the profiles behind it:
per library and length, the P-sites around the annotated starts and stop
codons and along the CDSs, from the pool's offsets, before the calibration.
`pitstop plot` draws both ([Usage](usage.md#qc-figures-plot)). Look at it
before you trust the start lengths: a length that the calibration moves peaks
at -1 or +1 in its library's profile, and its row gives its offset and the
move; a mixed length has two bars of about the same height at its peak.

The **kernel** is the shape of a start: per length, the harringtonine P-sites
minus the elongation P-sites (scaled on codons 11-20), as a share of nt -15
to +32. `kernel.tsv` has it per pool.

### 2. The window

A start's evidence is the P-sites in its **window**: nt -1 to +2 of its codon
(the codon and the nt before it). The kernel puts about half of a start's
P-sites on nt 0, a few % on nt -1, and the rest on codons 1-10. A wider
window, such as the codon +-1 codon, holds a few more of them, but it also
holds nt 0 of the codon before or after. Then two start codons one codon
apart get the same reads, and the scan cannot tell which one is used. With
the narrow window, the start codon wins against its neighbours in 90-95% of
the annotated starts with 30 or more harringtonine P-sites in the three
windows.

### 3. The background

Without a start, the harringtonine P-sites of a window follow the elongation
P-sites there, times a factor that depends on where the window is:

$$
\text{expected } h \;=\; s_t \, f_b \, \lambda .
$$

- $\lambda$ is the window's elongation rate (step 4).
- $f_b$ is a factor per **bin** of region and distance: leaders by nt
  before the CDS start (0, 30, 100, 300 and more), the CDS by codon from its
  start (0, 20, 50, 100, 200, 400, 700, 1,000 and more), trailers by nt after
  the stop codon, and transcripts without an annotated CDS.
- $s_t$ is a scale per transcript: harringtonine against elongation reads
  differ between transcripts.

The factors are fit on **null codons**, by a Poisson fit with alternating
updates. A null codon is a position with no start codon (ATG or one of its nine
near-cognates) within 8 nt, outside the window and codons 1-10 of an annotated
start, and with mostly unique reads. (Processed pseudogenes share their
parent's reads, start peak included.) The scale of a transcript comes from all
its positions but the windows and codons 1-10 of annotated starts.
`factors.tsv` has the factors.

### 4. The test

Write $R = s_t f_b$, $h$ and $e$ for the harringtonine and elongation P-sites
in the window. The window's elongation rate has a gamma prior. Its mean
$\lambda_0$ is the elongation density of the window's nt in their transcript:
per zone (leader, CDS frame 0, 1 and 2, trailer), as a window of 4 nt holds
one or two frame-0 nt of a CDS. Its shape $a$ is fit per region on the null
windows. Then

$$
h \sim \text{NB}\Big(\text{mean } m = R \max\!\big(e,\; \tfrac{a + e}{a/\lambda_0 + 1}\big),\;
\text{shape } a + e\Big),
$$

with an extra overdispersion $\phi$ (variance $+\,\phi m^2$). The posterior
mean lifts a window whose elongation reads are few by chance; the maximum keeps
an elongation peak (a pause) from being shrunk, so a pause is not called.

- **Pools** add their means and variances.
- **$\phi$** is calibrated per region and bin of expected reads on the null
  windows (one per 4 nt): the smallest $\phi$ that puts at most 1% of them at
  $p \le 0.01$. A scan calibrates it twice ([section 5](#5-the-scan)).
- $p = P(H \ge h)$; `enrichment` is $h / m$.
- $q$ is Benjamini-Hochberg per **class**: the annotated starts, and the other
  starts per region (leader, CDS, trailer, transcript without a CDS). With
  `--scan`, BH counts every start codon of the class: the ones without
  harringtonine reads are tests with $p = 1$.
- `typical_p` is the $p$ the window would have with the reads of a typical
  used start there: the median enrichment of the annotated starts with 3 or
  more expected reads. If it is large, the start cannot be called.
- `shoulder` is not part of the test. It is
  $\log_2 \frac{h + 1}{e + 1}$ over the codons 1-5 after the start (nt 3 to
  17), minus the same over the 5 codons before its window (nt -18 to -4),
  with the reads of all pools. Ribosomes that made a few peptide bonds before
  harringtonine acted raise the codons 1-5 of a used start. A pause, or a peak
  at a codon that cannot start, mostly does not raise them.
- `starts_stats.tsv` gives, per class, the rate of null windows at the $p$ of
  the last call, and the false calls that rate gives over the class's tests.

### 5. The scan

With `--scan`, the candidates are every ATG, CTG, GTG, TTG, AAG, ACG, AGG,
ATA, ATC and ATT with harringtonine reads in its window, and the annotated
starts. `--start-codons` gives other start codons, for example
`--start-codons ATG,CTG,GTG,TTG`. Then BH counts only these codons, but the
null codons stay clear of all ten. Starts near each other share reads, so
they are called in order:

1. The codons 1-10 of an annotated start are not called: harringtonine raises
   them when the start is used, and the test cannot tell them from a start.
2. The candidates are taken by $p$, smallest first. A candidate with
   $q < 0.05$ (`--fdr`) is **called** unless a called start stops it.
3. A called start stops the candidates within 8 nt of it and in its codons
   1-10. `stopped_by` names the start that stops a candidate.

Thus of two starts one or two codons apart, only the stronger is called.

**The null windows of a scan.** A used start raises its codons 1-10, an
annotated start or not: ribosomes that made a few peptide bonds before the
drug acted stay there. The null codons stay clear of the annotated starts
only, so null windows in the codons 1-10 of a uORF start look like starts and
make $\phi$ larger. Thus a scan calls twice. The first call gives the called
starts. The null windows in their window and codons 1-10 (nt -3 to +32) are
dropped, $\phi$ is calibrated again on the rest, and the second call is the
result. The factors and the prior keep all null windows. `starts_stats.tsv`
gives the second $\phi$ and the number of null windows dropped
(`null_windows_near_called`).

### 6. ORFs and the catalogue

Each called start runs to its first in-frame stop codon. Its type comes from
the coordinates against the annotated CDS: `uORF` (ends at or before the CDS
start), `uoORF` (starts before the CDS, out of frame, ends in it), `extension`
(starts before the CDS, in frame), `CDS`, `truncation` (inside the CDS, in
frame), `internal` (inside the CDS, out of frame), `dORF` (after the CDS) or
`other` (no annotated CDS). `start_orfs.tsv` has them, but the annotated CDSs.
`catalogue.tsv` has the uORFs and uoORFs, and every annotated CDS with its
start evidence. Both are ORF tables for `pitstop orfs --orfs`: count the
catalogue's ORFs in other libraries, and score their frames with
`pitstop score`.

## Reading the results

- **A call is start evidence, not translation evidence.** It says that
  pitstop found ribosomes stopped on that codon by harringtonine. Whether
  ribosomes go on to elongate is a frame question: `orfs` and `score` on the
  ORF table.
- **Power depends on the elongation reads.** A start with few elongation reads
  in its window has a wide posterior for its elongation rate, so it needs many
  harringtonine reads. Look at `expected` and `typical_p`.
- **Leaders need more.** Leaders have a higher harringtonine / elongation
  ratio than CDSs near their start, so a leader start needs more
  harringtonine reads than a CDS start with the same elongation reads.
- **Deep leader starts.** The expected false calls in `starts_stats.tsv`
  come from the null windows of all depths, and most null windows are
  shallow. Leader windows with many expected reads have more peaks at codons
  that cannot start. In the five pools of the validation, an estimated 11% of
  the leader calls with fewer than 30 expected reads are false, but 26% (95%
  interval 7-64%) of those with 30 or more. Most deep calls (64%) also raise
  their codons 1-5, as a used start does. Thus the deep calls stay in the
  output. Find them by `expected`, and look at their `shoulder` before you
  use one alone. For one start the shoulder is noisy: 28% of the annotated
  starts with 3 or more expected reads have a shoulder of 0 or less (median
  1.0), and so do 43-56% of the null-codon peaks that reach the $p$ of a
  call (median 0.0-0.2).
- **Trailers.** Harringtonine libraries can hold RNA fragments that the
  elongation libraries do not have, mostly in trailers and of short reads.
  `starts_stats.tsv` gives many expected false calls in trailers; use those
  calls with care.
- **Multimappers.** Each alignment counts once, so paralogs and pseudogenes
  share a start peak. `multimapped` gives the share of the window's reads
  with another alignment.

## Validation on real data

A human cell line, four experiments: five pools of three replicates (one
experiment in two conditions), 15 harringtonine libraries and their 15
elongation libraries, on a transcriptome with one transcript per gene (12,519
with a CDS).

- **Start lengths**: 26-29 nt. 30 nt is mixed in 14 of the 15 harringtonine
  libraries (the nt beside the start peak has 45-93% of its P-sites) and 25
  nt in 13; 26-29 nt in none (37% or less). The kernel holds 37-54% of its
  P-sites at nt 0 and 1-5% at nt -1.
- **Background**: leader factors are 1.3-3.3 near the CDS start and 0.7-1.6
  300 nt or more before it; CDS factors are 0.16-0.49 at codons 20-200 and
  0.8-0.9 at codons 700-1,000 (all against codons 1,000 and more).
- **Null windows** (held-out transcripts: $\phi$ from the null windows of
  half the transcripts, $p$ of the other half): 0.14-0.15x alpha at
  $p \le 0.01$ and 0.8-0.9x alpha at $p \le 0.001$ in leaders; 0.30-0.33x
  and 0.9-1.1x in CDSs.
- **Annotated starts** (scan): 60% called (7,526 of 12,519); 22% of those
  with less than 0.3 expected reads in the window, 26% with 0.3-1, 54% with
  1-3, 93% with 3 or more. Of the annotated starts with a call within 8 nt,
  the call is on the annotated codon in 99.1%.
- **Swaps**: with an elongation library as the "harringtonine" library of
  another replicate (replicate 3 on 2, and replicate 1 on 2, in each pool),
  2 and 1 of the 12,519 annotated starts are called, and 127 and 134 leader
  starts at q < 0.05 (6,105 with the harringtonine libraries).
- **Start context**: called leader starts have a purine at -3 in 66-69% and a
  G at +4 in 45-52%; start codons in leaders that are not called, at the
  same depth, 54-58% and 30-35% ([Kozak 1986](#references)). The sequence
  plays no part in the test.
- **Leader calls**: 6,105 leader starts have q < 0.05, and the scan calls
  4,692 of them (the near-start rule); the catalogue has 3,383 uORFs and 908
  uoORFs. False calls among the leader calls, estimated per bin of expected
  reads (the null windows' rate at the cut times the leader start codons of
  the bin): 11% below 30 expected reads, 14% in all. At 30 or more expected
  reads, the null codons of that depth give 26% (95% interval 7-64%). The
  harringtonine reads in codons 1-5 after the call do not point to more: 36%
  of the deep calls raise these codons more than 2-fold against the codons
  before the call, as do 36-41% of the calls with 1-30 expected reads. Too
  few null codons of that depth reach the cut (4) to give their share.

## References

- Benjamini Y, Hochberg Y (1995). Controlling the false discovery rate: a
  practical and powerful approach to multiple testing. *J R Stat Soc B* 57:289-300.
- Fresno M, Jiménez A, Vázquez D (1977). Inhibition of translation in
  eukaryotic systems by harringtonine. *Eur J Biochem* 72:323-330.
- Ingolia NT, Lareau LF, Weissman JS (2011). Ribosome profiling of mouse
  embryonic stem cells reveals the complexity and dynamics of mammalian
  proteomes. *Cell* 147:789-802.
- Kozak M (1986). Point mutations define a sequence flanking the AUG initiator
  codon that modulates translation by eukaryotic ribosomes. *Cell* 44:283-292.
- Lee S, Liu B, Lee S, Huang SX, Shen B, Qian SB (2012). Global mapping of
  translation initiation sites in mammalian cells at single-nucleotide
  resolution. *PNAS* 109:E2424-E2432.
