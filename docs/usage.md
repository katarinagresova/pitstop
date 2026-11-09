# Usage

Every command reads transcriptome alignments and writes tab-separated tables,
`<prefix>.<table>.tsv`; `inspect` also writes an HTML report, and `repair` a
new BAM. How each step works: [Method: quant](method.md),
[ORF counting](orfs.md) and [start sites](starts.md). What the flags of
`inspect` mean: [Flags](flags/index.md).

## Install

```bash
pip install git+https://github.com/katarinagresova/pitstop.git
```

pitstop needs Python 3.10 or later. pip also installs numpy, scipy, pandas,
pysam and matplotlib.

To run the tests, or to use the exact versions they ran with, use the conda
environment instead:

```bash
git clone https://github.com/katarinagresova/pitstop.git
cd pitstop
mamba env create -p .env -f environment.yml
.env/bin/pip install --no-deps -e .
.env/bin/python -m pytest          # tests on synthetic data with known truth
```

## Read QC: `inspect`

`pitstop inspect` checks one library: which read lengths keep one reading
frame, where their reads fall, what is at the read ends, and how many distinct
footprints the library has. It gives flags for common library problems, the
read lengths to use, and an HTML report. Run it on each library before the
other commands.

```bash
pitstop inspect --bam sample.bam --gtf annotation.gtf --fasta genome.fa --out-prefix out/sample
```

The inputs are those of `quant` (no `--read-lengths`: `inspect` looks at
every length). It reads the alignments that the other commands read: sense
alignments of matches with soft clips at the ends, on transcripts of the GTF.
A 5′ soft clip is not part of the footprint, and the length of a read is its
footprint length. A read with several alignments is split over them (over its
alignments to transcripts with a CDS, if it has some), so each read counts
once. The mismatch tables need the read sequences in the BAM.

Words of the tables:

- **CDS interior**: the read is at least 15 nt inside the CDS at both ends.
- **Phase** of a read: (5′ end − CDS start) mod 3. The **top-phase share**
  (`top_frame`) of a length: the share of its top phase in the CDS interior.
- **Start peak** of a length: its commonest 5′ end at the annotated start
  codons, in nt -20 to -4 from the first nt of the start codon. **Beside**:
  the larger of the two nt beside the peak, as a share of the peak.
- **Offset** of a length: minus its start peak, or 12 nt with fewer than 20
  reads at the peak. The P-site of a read is its 5′ end plus the offset.
- **Class** of a length, the first rule that holds:

    | Class | Rule |
    |---|---|
    | `few reads` | fewer than 100 reads in the CDS interior, or less than 0.3% of the reads |
    | `short state` | 20-22 nt, the start peak at nt -12 with 20 reads or more, and half of the P-sites or more in the CDS ([short state](flags/short-state.md)) |
    | `clean` | top-phase share 0.7 or more |
    | `no frame` | else |

- **Suggested read lengths**: of the runs of consecutive clean lengths, the
  one with the most reads.
- **After repair**: the values that the reads would have after `pitstop
  repair` (the same rule, on each read).

| Output | Content |
|---|---|
| `<prefix>.inspect.html` | the report (below) |
| `<prefix>.summary.tsv` | one row. The inputs (`library bam gtf fasta`: file names, no paths) and `pitstop` (the version). The reads: `reads_forward alignments_dropped_cigar reads_cigar_ok refs_not_in_gtf alignments_dropped_ref_not_in_gtf` (as `quant`'s `stats.tsv`), `reads` (the reads that count), `reads_on_cds_transcript`, `reads_nh1 reads_nh2 reads_nh3 reads_nh4plus` (by number of alignments). Complexity: `distinct_share` (distinct footprints per read), `new_per_read` (the new footprints per read in the last 10% of the reads), `distinct_share_at_depth` (in a subsample of 1 million reads; `NA` with fewer). Lengths: `common_length max_length`, `long_share` (32 nt or more), `short_share` (22 nt or less), `disome_share` (50-70 nt). The 3′ end spread: `end3_peak end3_reads end3_spread` ([under-digestion](flags/under-digestion.md)). Regions: `cds_share` (of the P-sites), `short_reads short_cds_share` (20-22 nt), `length_distance` (total-variation distance of the read lengths on CDSs and on transcripts without a CDS). Frame: `classed_lengths` (lengths with a class other than `few reads`), `clean_share clean_lengths read_lengths`, the same with `_repaired`, `best_length best_top_frame` (the length with the highest top-phase share, short state and few reads left out), `mixed_lengths`. Read ends: `mismatch5` (first nt), `mismatch5_inner` (the higher of nt 2 and 3), `mismatch5_background` (nt 4-7), `mismatch5_T` (the share of T among the mismatched first nt), `mismatch3` (the higher of the last two nt), `mismatch3_background`, `clip5_share` (reads with a 5′ soft clip), `untemplated_hidden` (an estimate: added nt that match the transcript), `bias_max bias_max_at` (the strongest nt preference at a cut, and where). The flags: `flag_weak_frame flag_long_reads flag_under_digestion flag_untemplated flag_end_spike flag_mixed_lengths flag_short_state flag_low_complexity` ([Flags](flags/index.md)). `repair`: `yes` if `repair` helps this BAM. |
| `<prefix>.lengths.tsv` | per footprint length: `reads share`; `interior` (reads in the CDS interior), `phase0 phase1 phase2` (their shares), `top_phase top_frame`; at the annotated starts, `start_reads` (reads with the 5′ end in the window), `start_peak start_peak_share start_beside`; at the annotated stops the same for the 3′ ends (`stop_reads stop_peak stop_peak_share stop_beside`, nt from the first nt of the stop codon); `offset`; the share of the P-sites per region, `5UTR CDS 3UTR transcript` (`3UTR`: the stop codon and after; `transcript`: on transcripts without a CDS); `class`; `mixed` (the start peak has 20 reads or more and beside it 0.5 of them or more: [mixed lengths](flags/mixed-lengths.md)); after repair, `reads_repaired top_frame_repaired class_repaired`. |
| `<prefix>.anchors.tsv` | `length anchor end x reads share`: per length, the reads by the position of the 5′ end at nt -30 to 0 from the annotated starts (`anchor` `start`), and by the position of the 3′ end at nt -5 to 30 from the first nt of the stop codons (`stop`). `share`: of the length's reads in that window. |
| `<prefix>.ends.tsv` | per length and end (`end` 5 or 3), and for all lengths (`length` `all`): `reads`; `clip1 clip2 clip3plus` (the share with a soft clip of 1, 2, 3 or more nt); `mm0`-`mm7` (the mismatch rate at positions 0-7 from that end; a soft clip included, so position 0 is the first nt of the read); `mm0_A`-`mm0_T` (the read nt at a mismatched position 0); `cut_A`-`cut_T` (the transcript nt at the cut: the first or the last nt of the footprint); `bias_A`-`bias_T` (the share of a nt at the cut over its share one codon inside the footprint). |
| `<prefix>.complexity.tsv` | `fraction reads distinct`: per subsample of 10%, 20%, ... 100% of the reads (fixed by the read name, no random numbers), the reads and the distinct footprints (transcript, 5′ end, length). |
| `<prefix>.inspect_general_mqc.json`, `<prefix>.inspect_lengths_mqc.json`, `<prefix>.inspect_frame_mqc.json` | for MultiQC (below). |

### The report

`<prefix>.inspect.html` is one file: it needs no other file and no network,
so you can send it or keep it with the data. It follows the light or dark
theme of the browser. Each number in it is also in a table, and a rerun gives
a byte-identical report.

- **Verdict and next steps**: the read lengths to use, as `--read-lengths`;
  the fix for the alignment and `repair`, if they help; advice for each
  problem that they cannot fix; and the table of the flags. Each flag links
  to its section of the report and to its [page](flags/index.md).
- **The reads that pitstop reads**: the reads after each filter, and by
  number of alignments.
- **Depth and complexity**: the distinct footprints against the reads.
- **Read lengths**: the reads per length, by class.
- **Regions per length**: where the P-sites of each length fall.
- **Frame per length**: the phases in the CDS interior, against the cut of a
  clean length, and after repair.
- **Ends at the annotated starts and stops**: the 5′ ends at the start codons
  and the 3′ ends at the stop codons, per length.
- **Read ends against the transcript**: the mismatches by position from each
  end, the soft clips, and the nucleotide preference at the cuts.
- **Not in this report**: what other tools check (raw reads, rRNA).

### Several libraries: MultiQC

`inspect` writes MultiQC files (names that end in `_mqc.json`). MultiQC
finds them, with no configuration, and shows all libraries in one report: a
row of the general statistics per library (the clean share, now and after
repair, the 5′ mismatch, the share of 32 nt or more, the distinct share at 1
million reads, the flags that fail and warn), the read lengths of all
libraries in one plot, and the top-phase share per length of all libraries.

```bash
multiqc out/
```

## Read repair: `repair`

`pitstop repair` clips the untemplated 5′ nucleotide: an extra nucleotide at
the 5′ end of a read that is not in the transcript. `inspect` says when it
helps (the flag [untemplated 5′ nt](flags/untemplated.md), and `repair` in
`summary.tsv`). Its report gives the fix at the alignment first: an aligner
that soft-clips the read ends removes the nucleotide, and a new alignment is
better than a repair. `repair` is for the BAM that you have: from another
pipeline, or from an aligner that cannot soft-clip.

```bash
pitstop repair --bam sample.bam --gtf annotation.gtf --fasta genome.fa --out-prefix out/sample.repaired
pitstop inspect --bam out/sample.repaired.bam --gtf annotation.gtf --fasta genome.fa --out-prefix out/sample.repaired
```

The rule looks at the first 2 nt of each read. From the first aligned nt on,
a nt that does not match the transcript becomes a soft clip. It stops at the
first match or at an N, and at least one nt stays aligned. A 5′ soft clip of
the aligner counts towards the 2 nt, so a second run changes nothing. The
rule applies to each alignment of a read (each isoform), and only to the
alignments that the other commands read; it copies all other records. It
changes the position, the CIGAR, and the tags NM and MD where the alignment
has them; not the sequence, the qualities or other tags. Every pitstop
command leaves a 5′ soft clip out of the footprint, so all of them read the
new BAM as it is. `inspect`'s values after repair use the same rule: they
are those of `inspect` on the new BAM.

| Output | Content |
|---|---|
| `<prefix>.bam` | the BAM, in the order of the input, with a `@PG` line in the header. A BAM sorted by coordinate stays sorted, and gets an index, `<prefix>.bam.bai`. |
| `<prefix>.repair.tsv` | `length clipped alignments reads`: per footprint length before the repair (the read length without a 5′ soft clip) and number of clipped nt (0-2), the alignments that `repair` looked at, and the reads (by their primary alignment). |

What `repair` does not do: it does not trim a fixed number of nucleotides
(trim them before the alignment: [end spike](flags/end-spike.md)), and no
rule for each read fixes [under-digestion](flags/under-digestion.md),
[mixed lengths](flags/mixed-lengths.md), or an added nucleotide that matches
the transcript by chance. It does not touch the 20-22 nt reads of the [short
state](flags/short-state.md): they are ribosomes. It does not filter by
length: each command has `--read-lengths`.

## Counts per CDS: `quant`

```bash
pitstop quant --bam sample.bam --gtf annotation.gtf --fasta genome.fa \
    --read-lengths 18-30 --out-prefix out/sample.human
# spike-in, reusing the human offsets (or leave --offsets out to estimate them from the spike-in reads)
pitstop quant --bam sample.spikein.bam --gtf yeast.gtf --fasta yeast.fa \
    --read-lengths 18-30 --offsets out/sample.human.offsets.tsv --out-prefix out/sample.yeast
```

Inputs:

- `--bam`: alignments to transcripts; reference names are transcript ids. Text
  after a `|` is ignored.
- `--gtf` / `--fasta`: the annotation and its genome. Transcript lengths must
  match the BAM header.
- `--read-lengths`: footprint lengths to use (required).

| Output | Content |
|---|---|
| `<prefix>.quant.tsv` | `Name Length EffectiveLength ritpm NumReads`, one row per CDS, sorted by Name. `Length` = CDS without the stop codon. `NumReads` = expected reads; they sum to the assigned reads. `ritpm` = reads per nt, scaled to sum to 1e6. `NA` where no read is compatible. |
| `<prefix>.offsets.tsv` | per (length, phase): `offset` (nt from the 5' end to the P-site codon), `support` (reads covering the first nt of the start codon or of the last sense codon; a read with n alignments counts 1/n), `z` (how firmly the data pin this phase's offset: score lead of the chosen window over the best window that gives this phase another offset, divided by the square root of the reads the two disagree on; below about 4 = weakly determined), `reads` |
| `<prefix>.ties.tsv` | `Name tie_group`: CDSs that no read tells apart (e.g. identical paralogs). Their split comes from the model, not from the data: equally long CDSs keep the EM's even start, otherwise the shortest gets almost all the reads. |
| `<prefix>.stats.tsv` | reads left after each filter, EM steps and log-likelihood |
| `<prefix>.psites.tsv` | `read Name psite length`, one row per alignment assigned to a CDS: `psite` is the 0-based transcript position of the P-site's first nt |
| `<prefix>.qc.tsv` | `length region phase reads`, for QC plots: every read once, at its footprint length, split evenly over its alignments to transcripts with a CDS (over all its alignments if it has none). `region` = where the P-site is: `5UTR`, `CDS` or `3UTR` (the stop codon and after); `NA` without a CDS, outside `--read-lengths` or without an offset. `phase` = (5' end − CDS start) mod 3, as in `offsets.tsv`, if the 5' end is in the CDS; else `NA`. At every read length. |

## ORFs other than the CDS: `orfs` and `score`

`pitstop orfs` counts reads on a table of ORFs — upstream ORFs (uORFs),
ORFs overlapping the CDS out of frame (uoORFs), and others — alongside the
annotated CDS, and keeps every read's reading frame instead of snapping it to
the CDS. `pitstop score` turns those per-codon counts into frame evidence for
each ORF, pooled over one or more `orfs` runs (e.g. replicate libraries):

```bash
pitstop orfs --bam sample_rep1.bam --gtf annotation.gtf --fasta genome.fa \
    --read-lengths 18-30 --orfs candidates.tsv --out-prefix out/sample_rep1
# the same for sample_rep2, with --offsets out/sample_rep1.offsets.tsv
pitstop score --orfs-prefix out/sample_rep1 out/sample_rep2 --out-prefix out/pooled
```

`--orfs` is a TSV `ORF_id Name start end` (0-based transcript coordinates,
`start` = first nt of the start codon, `end` = one past the last sense codon,
same convention as the CDS). The annotated CDSs are added automatically.
From the output of an ORF caller, or from genome coordinates: [ORF tables from
other tools](workflows/orf-tables.md).

To compare conditions, give every library the same GTF, ORF table and
offsets table, and score all libraries together:
[Typical workflow](orfs.md#typical-workflow). What the scores do and do
not say: [Reading the results](orfs.md#reading-the-results).

| Output | Content |
|---|---|
| `<prefix>.orfs.tsv` | `ORF_id Name type start end start_codon Length NumReads`, one row per ORF (from `--orfs`, plus the annotated CDSs) and per outside component (`<tx>:5UTR`, `<tx>:3UTR`, `<tx>:transcript`: the positions in no ORF). `type` is `CDS`, `uORF` (ends at or before the CDS start), `uoORF` (starts upstream of the CDS, out of frame with it) or `other` for an ORF, and `5UTR`, `3UTR` or `transcript` for an outside component. |
| `<prefix>.offsets.tsv` | the same table as `quant`'s (the same BAM and read lengths give the same offsets); either command's `--offsets` accepts it. `orfs` uses only each length's phase-0 offset. |
| `<prefix>.frames.tsv` | per read length: `offset` (the phase-0 offset, not a window), `reads`, and `frame0 frame1 frame2` — the share of P-sites in each frame of the CDS, from CDS interiors. |
| `<prefix>.codons.tsv` | `ORF_id codon length frame0 frame1 frame2`: P-sites per codon of every ORF but the annotated CDSs, by length and frame — `score`'s input. |
| `<prefix>.ties.tsv`, `<prefix>.stats.tsv`, `<prefix>.psites.tsv` | as for `quant`, but over every ORF and outside component: `ties.tsv` is `ORF_id tie_group`, and `psites.tsv` keeps every P-site, not only those assigned to a CDS. |
| `<prefix>.scores.tsv` | one row per scored ORF: `ORF_id Name type start end codons codons_with_reads reads in_frame_share leads expected z p min_p q` — `p`/`q` are the frame-evidence significance and its BH correction; `min_p` is the best `p` the ORF's length and reads could reach. |
| `<prefix>.decoys.tsv` | each ORF shifted by +1 and +2 nt and scored the same way, for null calibration. |
| `<prefix>.score_stats.tsv` | the read lengths used and how many ORFs/decoys were scored. |

## Start sites: `starts`

`pitstop starts` finds translation start sites in harringtonine (or
lactimidomycin) libraries. Each harringtonine library is compared with a
matched elongation library (same condition and replicate), so that a pause,
which both have, is not taken for a start. It reads the outputs of
`pitstop orfs` runs:

```bash
# one offset table per pool: quant on the merged elongation BAMs
pitstop quant --bam elong_merged.bam --gtf annotation.gtf --fasta genome.fa \
    --read-lengths 18-30 --out-prefix out/pool
# orfs on every library, harringtonine and elongation, with that table
pitstop orfs --bam harr_rep1.bam --gtf annotation.gtf --fasta genome.fa \
    --read-lengths 18-30 --offsets out/pool.offsets.tsv --out-prefix out/harr_rep1
# ... harr_rep2, elong_rep1, elong_rep2 the same way
pitstop starts --libraries libraries.tsv --gtf annotation.gtf --fasta genome.fa \
    --scan --out-prefix out/starts
```

`--libraries` is a TSV `pool harringtonine elongation`, one row per
replicate, each an out prefix of an `orfs` run. A pool is a set of
replicates of one condition. Several pools can be given, and their evidence
adds up. `--scan` tests every ATG and near-cognate start codon (CTG, GTG,
TTG, AAG, ACG, AGG, ATA, ATC, ATT) with harringtonine reads, and calls the
starts with BH q < 0.05; `--start-codons ATG,CTG` tests only these codons and
`--fdr 0.01` sets the cut. `--orfs table.tsv` tests the starts of an ORF
table instead (columns `ORF_id Name start`). The annotated starts are always
tested. Only the start lengths count: by default the read lengths that pass a
QC in 24-34 nt (`--length-window`), or those of `--start-lengths LO-HI`. How
it works and how to read it: [Start sites](starts.md). Read the report
([below](#the-report-of-starts)) before you use the calls.

| Output | Content |
|---|---|
| `<prefix>.starts.html` | the report ([below](#the-report-of-starts)) |
| `<prefix>.starts.tsv` | one row per tested start: `Name start codon region frame ORF_id harringtonine elongation expected enrichment p typical_p q multimapped shoulder`. `harringtonine` and `elongation` are the P-sites in the start's window (nt -1 to +2 of the start codon), summed over the replicates and pools; `expected` is the harringtonine P-sites that the elongation profile predicts there; `q` is BH per class (annotated starts; other starts per region: 5′ UTR, CDS, 3′ UTR, transcript without a CDS); `typical_p` is the p the window would have with the reads of a typical used start; `multimapped` is the share of the window's reads with another alignment; `shoulder` is log2 harringtonine / elongation P-sites (+1 each) in the codons 1-5 after the start, minus the same in the 5 codons before its window (a used start raises its codons 1-5; not part of the test). With `--scan` also `called stopped_by end type`. |
| `<prefix>.start_orfs.tsv` | with `--scan`: the ORF of each called start, to its first in-frame stop, but the annotated CDSs: `ORF_id Name start end type start_codon` and the start evidence. `ORF_id` is `<transcript>_<1-based start>`. `type` is `uORF`, `uoORF`, `extension` (in frame upstream of the CDS), `truncation` (in frame inside it), `internal` (out of frame inside it), `dORF` (after it) or `other` (no annotated CDS). An ORF table for `pitstop orfs --orfs`. |
| `<prefix>.catalogue.tsv` | with `--scan`: the uORFs and uoORFs of the called starts, and every annotated CDS with its start evidence; the columns of `start_orfs.tsv`. An ORF table for `pitstop orfs --orfs`. |
| `<prefix>.lengths.tsv` | the QC of the read lengths, per library and read length: `library pool kind length psites share 5UTR CDS 3UTR transcript frame0 frame1 frame2 interior_share top_frame periodic start_peak at_peak beside_peak beside_nt beside_share mixed pair calibrated shift offset qc start_length`. `psites` are its P-sites, `share` their share of the library's, `5UTR`-`transcript` the same per region; `frame0-2` are its P-sites per CDS frame, 15 nt or more inside the CDSs, `interior_share` their share of the library's such P-sites, and `top_frame` the share of the top frame; `periodic`: 100 of these P-sites or more and `top_frame` 0.4 or more. `start_peak` is the nt (-15 to +32 from the annotated starts) where its P-sites peak, `at_peak` their number there, `beside_peak` the number at the larger of the two nt beside it, `beside_nt` that nt and `beside_share` `beside_peak / at_peak`; `mixed`: `beside_share` is 0.5 or more, with 20 P-sites or more at the peak (two offsets); `pair`: the length is mixed in half of the harringtonine libraries or more, and its anchor is the pair of `start_peak` and `beside_nt`, moved onto nt 0 and 1. `shift` is the nt its P-sites move, and `offset` the offset they end at: nt from the 5' end to the P-site, the `orfs` run's offset plus `shift`. `qc` is the verdict on the length, `pass` or why not ([start lengths](starts.md#1-calibration-and-start-lengths)). It is written before the test, also when no length passes. |
| `<prefix>.profiles.tsv` | the profiles behind `lengths.tsv`, per library and read length, from the pool's offsets, before the calibration: `library length anchor x frame psites`. `anchor` `start`: P-sites at nt `x` (-15 to +62) from the annotated starts; `stop`: at nt `x` (-45 to +14) from the first nt of the stop codons; `cds`: in bin `x` of 20 along the CDSs. `frame` is the frame of the CDS. |
| `<prefix>.pairs.tsv` | one row per row of `--libraries` (a pair): `pool harringtonine elongation distance`. `distance` is the read-length distance of the two libraries: half the summed absolute difference of their `interior_share` per length (0: the same lengths; 1: none in common). Pairs of one digestion are below 0.3. Written also when no length passes. |
| `<prefix>.libraries.tsv` | per library, the P-sites of the start lengths after the calibration: `library pool kind psites start_psites start_share 5UTR CDS 3UTR transcript interior frame0 frame1 frame2 at_start start_enrichment`. `start_share` is `start_psites / psites`; `5UTR`-`transcript` their shares per region; `interior` those 15 nt or more inside the CDSs, `frame0-2` their shares per frame; `at_start` the P-sites at the first nt of the annotated start codons, and `start_enrichment` the start peak: `at_start` per annotated start over the P-sites per in-frame codon inside the CDSs. |
| `<prefix>.start_codons.tsv` | per class (`annotated`, `5UTR`, `CDS`, `3UTR`, `transcript`) and start codon: `class codon tests test_share called call_share`. `tests` are the start codons the class tests (with `--scan`, all of the class, also those without harringtonine reads), `called` the calls (with `--scan`, `called` of `starts.tsv`; else q below `--fdr`), and the shares are of the class's. |
| `<prefix>.kernel.tsv` | per pool, read length and nt from the annotated start (-15 to +32): harringtonine and elongation P-sites (after the calibration), and the kernel (harringtonine minus scaled elongation, as a share). |
| `<prefix>.factors.tsv` | per pool: the harringtonine / elongation factor of each bin of region and distance. |
| `<prefix>.starts_stats.tsv` | the inputs (`pitstop`, the version; the file names of `libraries_sheet gtf fasta`; `candidates_from`, `scan` or the ORF table; `fdr`; `start_lengths_given`), the length window, the start lengths, P-site totals, the prior and overdispersion per region, and per class the tests, the calls and the false calls expected from null windows. With `--scan`, the overdispersion is the second one, and `null_windows_near_called` is the number of null windows dropped near the starts of the first call. |

### The report of `starts`

`<prefix>.starts.html` is one file, as the report of `inspect`: no other
file, no network, the light or dark theme of the browser, each number also in
a table, and a byte-identical rerun. `starts` writes it also when no read
length passes the QC, before it stops, so the report says why.

The sections before **Pair detail** are for the whole run; the less
important parts of a section are folded until you open them.

- **Verdict and next steps**: the checks of the run. Beside them, **Used in
  this run** lists what `starts` did (the start lengths, the mixed and pair
  lengths, the offsets that the calibration moved: nothing to pass on), and **Next
  steps** lists what to do about the checks.
- **Calls**: per class the tests, the calls and the false calls that the null
  windows expect; with `--scan`, the ORFs of the called starts by type; the
  start codons of the calls against those of the tests (real starts are rich
  in ATG, then CTG and GTG; calls whose codons look like the tests are mostly
  noise: restrict `--start-codons`, or lower `--fdr`).
- **Pairs**: a pair is a row of `--libraries`. Its two libraries must have
  the same digestion: the read-length distance of their P-sites inside the
  CDS is below 0.3 (initiation libraries digested to longer reads were at
  0.59-0.74 from their elongation libraries). And the drug must work: the initiation library has the higher
  start peak. A pair that fails either check gives a wrong background to its
  initiation library: pair it with an elongation library of the same
  digestion, or leave it out.
- **Start lengths**: the QC verdict on each read length in the window. A
  start length has most of its P-sites in frame 0 and a peak at nt 0 of the
  annotated starts. Folded, per pool: the P-sites per length of each library
  (the start lengths dark, mixed lengths hatched, the offset above each bar),
  with the regions and frames of the start-length P-sites; and the
  harringtonine minus the elongation P-sites around the annotated starts,
  per length (where the drug works and the offset is right, most of a row is
  at nt 0).
- **Pair detail**: one pair at a time, chosen in the list at its top (the
  list stays in view as you scroll). The read lengths of the two libraries
  and their frames per read length; the P-sites of each start length at the
  annotated starts, initiation above and elongation below, each with its own
  scale. Folded: the other read lengths, the numbers per read length, the
  stop codons, and the P-sites along the CDS (in an initiation library they
  climb from the 5' end, where the ribosomes ran off).
