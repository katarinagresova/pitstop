# Comparing conditions

The workflow for two or more conditions, for example a control and a
knockdown, a drug or a stress, with initiation and elongation libraries. For
elongation libraries only, leave out the start sites (step 3) and give your
own ORF table.

## The experiment

- **Libraries.** Replicates of each condition:
  - elongation libraries;
  - RNA-seq libraries, if you want translation efficiency (TE);
  - initiation libraries with their matched elongation libraries, if you want
    start sites and new ORFs ([Initiation and elongation
    libraries](initiation.md)).
- **One protocol for all.** The same digestion, size selection and library
  kit in every condition. Otherwise the read lengths differ between the
  conditions, and so do the frame and the offsets.
- **Questions.** Which CDSs change their translation, or their translation
  efficiency? Which uORFs change, relative to the CDS of their transcript?
- **What pitstop does.** It gives the counts and the labels: reads per CDS
  (`quant`), reads per ORF, leader and trailer (`orfs`), frame evidence
  (`score`), and start sites with a catalogue of ORFs (`starts`).
- **What it does not do.** The tests of differences. Run them on pitstop's
  counts with a model for count data, for example DESeq2
  ([Love et al. 2014](#references)), as deltaTE does for translation
  efficiency ([Chothani et al. 2019](#references)).

![The workflow for comparing conditions, one card per step. Step 2: pitstop quant on every library writes quant.tsv and qc.tsv. Step 3: one pitstop starts --scan run with a pool per condition (or your own ORF table) writes catalogue.tsv. Step 4: pitstop quant on all elongation libraries merged writes offsets.tsv. Step 5: pitstop orfs --orfs on each elongation library, with the catalogue and these offsets, writes orfs.tsv. Step 6: pitstop score on all libraries of both conditions writes scores.tsv. Step 7, outside pitstop: your tests, for example with DESeq2 or deltaTE, take quant.tsv, orfs.tsv and scores.tsv.](../img/workflow-conditions.svg)

## Three rules

1. **One GTF, one ORF table and one offsets table** for every library that you
   count for the comparison. The offsets of short read lengths can differ
   between libraries by a codon or more. For a short ORF, different offsets
   move reads across its ends differently in each library: a difference
   between the conditions that is not biological.
2. **Find the ORFs and score them without the condition labels.** Make the
   catalogue from all conditions in one `starts` run, with one pool per
   condition, and score all libraries of all conditions together. ORFs called
   in one condition and then tested between conditions are biased towards the
   ORFs with more reads in that condition.
3. **Labels describe; they do not select.** Many uORFs are too short, or have
   too few reads, for `score` to call them ([Reading the
   results](../orfs.md#reading-the-results)). Test all ORFs of the table, and
   use the groups called, not called and cannot be called to read the results.

## Steps

The examples have two conditions, `ctrl` and `treat`, with replicates 1-3:
elongation libraries `elong_ctrl_1` ... `elong_treat_3`, initiation libraries
`harr_ctrl_1` ... `harr_treat_3`.

### 1. Prepare the reads

[Preparing the reads](reads.md), the same for every library.

### 2. QC and CDS counts

```bash
for lib in elong_{ctrl,treat}_{1,2,3} harr_{ctrl,treat}_{1,2,3}; do
    pitstop quant --bam $lib.bam --gtf annotation.gtf --fasta genome.fa \
        --read-lengths 18-34 --out-prefix out/$lib
done
```

Check each library's read lengths and frame ([Elongation
libraries](elongation.md#2-check-the-read-lengths)), each initiation library
against its elongation library ([Initiation and elongation
libraries](initiation.md#2-qc-of-each-library-and-of-each-pair)), and the
read lengths of the conditions against each other. Choose `--read-lengths`
for the next steps from these checks.

`quant.tsv` of the elongation libraries gives the CDS counts for the tests.
With each library's own offsets, the CDS counts are fine: a whole CDS has
many codons, and an offset changes only the reads at its two ends.

### 3. Start sites and the catalogue

As in [Initiation and elongation libraries](initiation.md), steps 3-5, with
one pool per condition: each pool gets offsets from its own elongation
libraries, and one `starts` run takes both pools.

```text
pool    harringtonine         elongation
ctrl    out/harr_ctrl_1       out/elong_ctrl_1
ctrl    out/harr_ctrl_2       out/elong_ctrl_2
ctrl    out/harr_ctrl_3       out/elong_ctrl_3
treat   out/harr_treat_1      out/elong_treat_1
treat   out/harr_treat_2      out/elong_treat_2
treat   out/harr_treat_3      out/elong_treat_3
```

Each pool keeps its own background, and the evidence of the pools adds up.
The calls thus come from both conditions. Without initiation libraries, use
an ORF table from another source ([ORF tables from other
tools](orf-tables.md)).

### 4. One offsets table for the counts

```bash
samtools cat -o elong_all.bam elong_{ctrl,treat}_{1,2,3}.bam
pitstop quant --bam elong_all.bam --gtf annotation.gtf --fasta genome.fa \
    --read-lengths 18-34 --out-prefix out/elong_all
```

### 5. Counts per ORF

```bash
for lib in elong_{ctrl,treat}_{1,2,3}; do
    pitstop orfs --bam $lib.bam --gtf annotation.gtf --fasta genome.fa \
        --read-lengths 18-34 --orfs out/starts.catalogue.tsv \
        --offsets out/elong_all.offsets.tsv --out-prefix out/$lib.catalogue
done
```

### 6. Frame evidence, over both conditions

```bash
pitstop score --orfs-prefix out/elong_{ctrl,treat}_{1,2,3}.catalogue --out-prefix out/catalogue
```

### 7. Your tests

pitstop's counts go into your own model. Two common tests:

- **Translation efficiency.** The CDS counts of the elongation libraries
  (`NumReads` in `quant.tsv`) with the RNA-seq counts of the same transcripts
  or genes (sum pitstop's transcripts per gene if the RNA-seq counts are per
  gene). Test the interaction of assay (Ribo-seq, RNA-seq) and condition, as
  deltaTE does.
- **A uORF relative to its CDS.** The uORF's `NumReads` and the `NumReads` of
  the CDS of the same transcript, from each library's `orfs.tsv` (same `Name`).
  Test the interaction of ORF type (uORF, CDS) and condition. Both counts come
  from one model, in which reads that the uORF and the CDS share are split by
  their frame, so the ratio is consistent.

`NumReads` are expected counts, not whole numbers. Round them for a model
that needs integers.

## Notes

- **Read a uORF change against the global change.** A condition can move
  ribosomes from the CDSs to the leaders as a whole. Then most uORFs go up
  relative to their CDS, and a test against no change finds many. Compare
  with the change of all leaders (the `<transcript>:leader` rows of
  `orfs.tsv`) or of all uORFs.
- **Per gene, the uORF counts are small.** A test of the interaction over many
  uORFs at once (all uORFs against all CDSs) has more power than one test per
  uORF.
- **Nested ORFs.** ORFs with the same stop share their reads in a way the data
  cannot resolve: sum their counts before the test ([Known
  limits](../orfs.md#known-limits)).
- **Spike-ins.** If the conditions change the total amount of translation,
  normalise with a spike-in: count its reads with `quant` and the offsets of
  the main library ([Usage](../usage.md#counts-per-cds-quant)).
- **Start usage between conditions.** `starts` calls starts per pool. It does
  not give start counts per library, so it does not test start usage between
  conditions. An initiation peak also counts ribosomes held on a start during
  the drug treatment, and a ribosome held on a uORF start can block scanning
  ribosomes behind it: a change in a peak is not the same as a change in how
  often the start is used.

## References

- Chothani S, Adami E, Ouyang JF, Viswanathan S, Hubner N, Cook SA, Schafer S,
  Rackham OJL (2019). deltaTE: Detection of translationally regulated genes by
  integrative analysis of Ribo-seq and RNA-seq data. *Curr Protoc Mol Biol*
  129:e108. [doi:10.1002/cpmb.108](https://doi.org/10.1002/cpmb.108)
- Love MI, Huber W, Anders S (2014). Moderated estimation of fold change and
  dispersion for RNA-seq data with DESeq2. *Genome Biol* 15:550.
  [doi:10.1186/s13059-014-0550-8](https://doi.org/10.1186/s13059-014-0550-8)
