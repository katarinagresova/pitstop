# Preparing the reads

Every pitstop command reads footprints aligned to the transcript sequences (a
transcriptome BAM), with the GTF and the genome FASTA that the transcripts come
from. This page says what that BAM must hold and how to make it. Prepare every
library of an experiment the same way: elongation and initiation libraries,
every condition, every replicate.

## What pitstop needs

| Input | Requirement |
|---|---|
| BAM | Alignments to transcripts. Each reference is one transcript, named by its transcript id. The BAM does not need to be sorted or indexed. |
| GTF | The annotation that the transcripts come from. The length of each BAM reference must equal the length of that transcript in the GTF. A GTF with fewer transcripts works: alignments to transcripts that are not in it are ignored. |
| Genome FASTA | The genome of the GTF. pitstop builds the transcript sequences from the two. |

In a reference name, text after a `|` is ignored, so the headers of
GENCODE's transcript FASTA work as they are.

From the BAM, pitstop uses ([Method: reads](../method.md#2-reads)):

- **Sense alignments.** Ribo-seq libraries are stranded, with the read in the
  sense orientation of the RNA. If your protocol sequences the other strand,
  reverse-complement the reads before the alignment. `stats.tsv` shows the
  problem: few `reads_forward`.
- **Every alignment of a read.** pitstop splits a read over all its
  alignments: the isoforms of a gene and its paralogs. The aligner must report
  all of them, not one picked at random.
- **Matches with soft clips at the ends.** A 5′ soft clip is taken as
  untemplated nucleotides, and the footprint starts at the first aligned base.
  Alignments with insertions, deletions or splices are dropped.
- **One read per read name.** Read names must be unique. This matters when you
  merge BAMs: reads straight from the sequencer or from SRA have unique names.

## Before the alignment

1. **Adapters and read structure.** Remove the 3′ adapter and every extra
   nucleotide your library kit adds to the insert: UMIs and random
   nucleotides, at either end. Follow the kit's read structure exactly. A
   nucleotide left at the 5′ end moves the footprint's 5′ end; one left at the
   3′ end makes the footprint longer. Both blur the reading frame.
2. **Insert lengths.** Keep inserts of footprint size, for example 18-40 nt,
   and drop the reads in which no adapter was found. The range only needs to
   hold the lengths you will analyse: every command chooses its lengths again
   with `--read-lengths`.
3. **Contaminants.** Remove the reads of rRNA, tRNA, snoRNA and other abundant
   non-coding RNAs: align to their sequences first and keep the reads that do
   not align. Fragments of these RNAs that resemble parts of annotated
   transcripts otherwise make peaks that are not ribosomes.
4. **Duplicates.** If the library has UMIs, remove PCR duplicates with them
   (for example with UMI-tools). Without UMIs, do not collapse identical reads.
   Many ribosomes stop on the same codon, at pauses and above all at the start
   codons of initiation libraries, so identical reads are often different
   footprints.

For a kit without random nucleotides, one cutadapt call does steps 1 and 2:

```bash
cutadapt -a ADAPTER -m 18 -M 40 -o sample.trimmed.fastq.gz sample.fastq.gz
```

A read without the adapter keeps its full length, which is longer than 40 nt,
so `-M 40` drops it (this needs reads longer than the longest insert you
keep). For a kit that adds, say, 4 random nucleotides at each end of the
insert, cut the adapter first and the random nucleotides after, in a second
call (cutadapt removes `-u` bases before it looks for the adapter):

```bash
cutadapt -a ADAPTER --discard-untrimmed -o sample.noadapter.fastq.gz sample.fastq.gz
cutadapt -u 4 -u -4 -m 18 -M 40 -o sample.trimmed.fastq.gz sample.noadapter.fastq.gz
```

## Alignment

pitstop works with any aligner whose BAM meets the requirements above. We
recommend STAR. bowtie is an alternative.

### STAR (recommended)

STAR aligns the reads to the genome with the annotated splice junctions, and
writes the same alignments in transcript coordinates
(`--quantMode TranscriptomeSAM`). One run gives a genome BAM, for genome
browsers and genome-based tools, and a transcriptome BAM for pitstop. In its
default local mode STAR soft-clips an untemplated 5′ nucleotide, and pitstop
leaves it out of the footprint.

The index, once per genome and annotation:

```bash
STAR --runMode genomeGenerate --runThreadN 8 --genomeDir star_index \
    --genomeFastaFiles genome.fa --sjdbGTFfile annotation.gtf
```

The alignment, per library:

```bash
STAR --runThreadN 8 --genomeDir star_index \
    --readFilesIn sample.trimmed.fastq.gz --readFilesCommand zcat \
    --outSAMtype BAM Unsorted --outSAMattributes NH HI AS NM --outSAMmode NoQS \
    --quantMode TranscriptomeSAM --quantTranscriptomeSAMoutput BanSingleEnd \
    --seedSearchLmax 10 --outFilterMultimapNmax 255 --outFilterMismatchNmax 2 \
    --outFilterMultimapScoreRange 0 --outFilterIntronMotifs RemoveNoncanonical \
    --outFileNamePrefix out/sample.
```

`out/sample.Aligned.toTranscriptome.out.bam` is the `--bam` of pitstop. Give
pitstop the GTF and the genome FASTA of the index.

| Option | Why |
|---|---|
| `--quantMode TranscriptomeSAM` | writes each alignment in the coordinates of every annotated transcript it fits |
| `--quantTranscriptomeSAMoutput BanSingleEnd` | keeps the soft clips in the transcriptome BAM. The default of this option (`BanSingleEnd_BanIndels_ExtendSoftclip`, made for RSEM) extends every soft clip back over the clipped bases, so an untemplated 5′ nucleotide becomes part of the footprint again. In STAR 2.7.10b and older, the option is `--quantTranscriptomeBan Singleend` |
| local alignment (STAR's default: do not set `--alignEndsType EndToEnd`) | soft-clips read ends that do not match, such as an untemplated 5′ nucleotide |
| `--outFilterMultimapNmax 255`, `--outFilterMultimapScoreRange 0` | keeps reads with up to 255 loci, with their best alignments only. pitstop splits such reads over their alignments. With STAR's default of 10 loci, reads from large gene families are lost |
| `--outFilterMismatchNmax 2` | at most 2 mismatches, for reads this short |
| `--seedSearchLmax 10` | shorter seeds, so that short reads find their loci |
| `--outFilterIntronMotifs RemoveNoncanonical` | drops alignments over non-canonical junctions |
| `--outSAMattributes NH HI AS NM`, `--outSAMmode NoQS` | tags for other tools, and no base qualities (smaller files). pitstop does not need them |

Sort and index the genome BAM for a genome browser. Leave the transcriptome
BAM as STAR writes it.

### bowtie (an alternative)

bowtie aligns straight to the transcript sequences. It is simple and fast, but
it gives no genome BAM, and it aligns end to end: an untemplated 5′ nucleotide
stays in the alignment as a mismatch, and that footprint then starts 1 nt too
early. Use bowtie only if your libraries have few untemplated nucleotides
([check it](#check-the-alignment)).

The transcript sequences, from the same GTF and genome, and their index:

```bash
gffread -w transcripts.fa -g genome.fa annotation.gtf
bowtie-build --threads 8 transcripts.fa tx_index
```

The alignment, per library:

```bash
bowtie -p 8 -v 2 --norc -k 256 --best --strata --reorder --sam -x tx_index \
    sample.trimmed.fastq.gz | samtools view -b -o sample.bam -
```

| Option | Why |
|---|---|
| `--norc` | sense alignments only |
| `-k 256 --best --strata` | all alignments of a read with its fewest mismatches, up to 256. bowtie's default reports one alignment, so a read of a gene with several isoforms would count for one of them only |
| `-v 2` | at most 2 mismatches, end to end |
| `--reorder` | writes the reads in input order with several threads, so that a rerun gives the same BAM and the same pitstop outputs |

## Check the alignment

Run `pitstop quant` on each library (see the workflows). Then:

- **`stats.tsv`** counts the reads after each filter: the sense reads, the
  reads with a usable alignment, the reads assigned to a CDS. A large loss at
  one step points to its cause: antisense reads, alignments with indels, or a
  reference that does not match the GTF.
- **`qc.tsv`** gives the reads per length and frame. Footprints of good
  libraries have most reads at a few lengths, and at those lengths most reads
  in one frame. How to read it: [Elongation
  libraries](elongation.md#2-check-the-read-lengths).
- **Untemplated nucleotides (bowtie).** bowtie writes the MD tag, which starts
  with `0` and the reference base when the first base of the read does not
  match. The share of reads with such a first base:

    ```bash
    samtools view sample.bam | awk '$1 != r {r = $1; n++; if (/MD:Z:0[ACGT]/) m++} END {print m / n}'
    ```

    A share of a few percent or less is normal. Tens of percent means that many
    reads carry an untemplated nucleotide: align these libraries with STAR.
