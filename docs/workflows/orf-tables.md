# ORF tables from other tools

`pitstop orfs --orfs` counts the reads of the ORFs in a table that you give
it. This page shows how to make that table from the output of an ORF caller,
or from ORFs in genome coordinates.

## The format

A TSV file with a header line and these columns:

| Column | Value |
|---|---|
| `ORF_id` | a name for the ORF, unique in the table |
| `Name` | the transcript id, as in the GTF |
| `start` | the first nt of the start codon, 0-based |
| `end` | one past the last sense codon: the 0-based position of the first nt of the stop codon |

The positions are in the transcript sequence, counted from its 5′ end.
`end - start` is the length of the ORF without its stop codon, a multiple of
3. pitstop adds the annotated CDSs to the table.

```text
ORF_id                  Name                 start   end
ENST00000000001_uORF1   ENST00000000001      40      85
```

Here nt 40-84 of the transcript are the codons of the uORF, and nt 85-87 are
its stop codon.

## Use the same annotation

Run the ORF caller with the GTF and the genome that you give pitstop.
Transcript coordinates are only valid for the transcripts that the tool used:
in another annotation release, a transcript with the same id can have other
exons, and the same positions then point to other nucleotides. The transcript
ids must be the same as in the GTF, with the same version suffix.

## From an ORF caller

Many ORF callers give each ORF in the coordinates of its transcript. Convert
the columns:

| Tool | Table | `Name` | `start` | `end` |
|---|---|---|---|---|
| [RiboCode](https://github.com/xryanglab/RiboCode) | the ORF table (`.txt` or `_collapsed.txt`) | `transcript_id` | `ORF_tstart - 1` | `ORF_tstop - 3` |
| [Ribo-TISH](https://github.com/zhpn1024/ribotish) | the output of `ribotish predict` | `Tid` | `Start` | `Stop - 3` |
| [RiboTIE](https://github.com/TRISTAN-ORF/RiboTIE) | the table of predictions (`.csv`) | `transcript_id` | `TIS_pos - 1` | `TTS_pos - 1` |

- **RiboCode** (version 1.2.15): `ORF_tstart` and `ORF_tstop` are 1-based,
  and `ORF_tstop` is the last nt of the stop codon.
- **Ribo-TISH** (version 0.2.8): `Start` and `Stop` are 0-based, and `Stop`
  is one past the stop codon. An ORF that runs to the end of its transcript
  has no stop codon, and pitstop drops it.
- **RiboTIE**: `TIS_pos` and `TTS_pos` are 1-based, and `TTS_pos` is the first
  nt of the stop codon.

For example, for RiboCode:

```python
import pandas as pd

r = pd.read_csv("ribocode.txt", sep="\t")
pd.DataFrame({"ORF_id": r["ORF_ID"], "Name": r["transcript_id"],
              "start": r["ORF_tstart"] - 1, "end": r["ORF_tstop"] - 3}
             ).to_csv("orfs.tsv", sep="\t", index=False)
```

For Ribo-TISH, which has no ORF id column, make one, for example from `Tid`
and `Start`. For a tool that is not in the table, find out from its
documentation whether its positions are 0-based or 1-based and whether its
ORFs include the stop codon. Then [check the table](#check-the-table).

## From genome coordinates

Other tools and ORF catalogues give the ORFs in genome coordinates, often as a
BED file. To convert an ORF, map two of its positions to the transcript: the
first nt of its start codon and its last nt. Both lie in exons, so this also
works for ORFs that span introns. Use the exons of the GTF that you give
pitstop.

### With a transcript for each ORF: Python

Your table gives the transcript of each ORF and two genome positions, 1-based:
`first`, the first nt of the start codon, and `last`, the last nt of the ORF.
On the minus strand, `last` is smaller than `first`.

```text
ORF_id   transcript_id     first     last
orf1     ENST00000000001   1000141   1003260
```

```python
import re

import pandas as pd

STOP_IN_ORF = True  # does each ORF of your table end with its stop codon?

exons = {}  # transcript id -> (strand, its exons from 5' to 3', 1-based)
for line in open("annotation.gtf"):
    f = line.split("\t")
    if len(f) > 8 and f[2] == "exon":
        tid = re.search(r'transcript_id "([^"]+)"', f[8]).group(1)
        exons.setdefault(tid, (f[6], []))[1].append((int(f[3]), int(f[4])))
for strand, ex in exons.values():
    ex.sort(reverse=strand == "-")


def to_transcript(tid, pos):
    """The 0-based transcript position of the 1-based genome position pos."""
    strand, ex = exons[tid]
    done = 0
    for s, e in ex:
        if s <= pos <= e:
            return done + (pos - s if strand == "+" else e - pos)
        done += e - s + 1
    raise ValueError(f"{pos} is not in an exon of {tid}")


orfs = pd.read_csv("orfs_genome.tsv", sep="\t")  # ORF_id transcript_id first last
start = [to_transcript(t, p) for t, p in zip(orfs["transcript_id"], orfs["first"])]
end = [to_transcript(t, p) + 1 - 3 * STOP_IN_ORF for t, p in zip(orfs["transcript_id"], orfs["last"])]
pd.DataFrame({"ORF_id": orfs["ORF_id"], "Name": orfs["transcript_id"], "start": start, "end": end}
             ).to_csv("orfs.tsv", sep="\t", index=False)
```

From a BED file: on the plus strand, `first` is `chromStart + 1` and `last`
is `chromEnd`; on the minus strand, `first` is `chromEnd` and `last` is
`chromStart + 1`. If the BED lines are whole transcripts with the ORF as the
thick part, use `thickStart` and `thickEnd` in their place.

### Without a transcript for each ORF: R

The function `mapToTranscripts` of the Bioconductor package
[GenomicFeatures](https://bioconductor.org/packages/GenomicFeatures/) finds
every transcript that holds a position. The script below keeps the
transcripts that hold all of the ORF: those on which the ORF has the length
of its BED blocks. An ORF can then be in the table more than once, once for
each of its transcripts. It also needs the packages
[txdbmaker](https://bioconductor.org/packages/txdbmaker/) and
[rtracklayer](https://bioconductor.org/packages/rtracklayer/).

```r
library(GenomicFeatures)  # Bioconductor, with txdbmaker and rtracklayer
stop_in_orf <- TRUE       # does each ORF of your BED end with its stop codon?

txdb <- txdbmaker::makeTxDbFromGFF("annotation.gtf")
exons <- exonsBy(txdb, by = "tx", use.names = TRUE)
orfs <- rtracklayer::import("orfs.bed")  # BED12: one ORF per line, its exons as blocks

first <- mapToTranscripts(resize(orfs, 1, fix = "start"), exons)  # the first nt of the start codon
last <- mapToTranscripts(resize(orfs, 1, fix = "end"), exons)     # the last nt of the ORF
m <- merge(data.frame(i = first$xHits, Name = as.character(seqnames(first)), start = start(first) - 1),
           data.frame(i = last$xHits, Name = as.character(seqnames(last)), end = end(last) - 3 * stop_in_orf))
# keep the transcripts whose exons hold all of the ORF
m <- m[m$end - m$start == sum(width(orfs$blocks))[m$i] - 3 * stop_in_orf, ]
write.table(data.frame(ORF_id = paste0(orfs$name[m$i], "_", m$Name), Name = m$Name, start = m$start, end = m$end),
            "orfs.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
```

We tested both scripts on the CDSs of 600 transcripts, given as genome
positions and as a BED file: on both strands, with and without introns, each
CDS came back at its coordinates in pitstop's annotation.

## Check the table

- **Before the run.** `end - start` is a multiple of 3, and equals the length
  of the ORF that your tool gives, without the stop codon.
- **After `pitstop orfs`.** `stats.tsv` counts the rows of your table
  (`orf_table_rows`) and the rows dropped, by reason:
  `orfs_dropped_transcript_not_in_gtf`, `orfs_dropped_not_multiple_of_3`,
  `orfs_dropped_off_transcript` and `orfs_dropped_no_stop`. With a correct
  conversion and the same annotation, few rows are dropped. Many
  `no_stop` drops mean that the positions are off: by 3 nt if the stop codon
  was counted into the ORF, by 1 nt if 1-based positions were used as
  0-based.
- **The start codons.** `orfs.tsv` has the `start_codon` of each ORF. If few
  of them are ATG or a near-cognate codon, `start` is off.
