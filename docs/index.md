# pitstop

A Python Ribo-seq toolkit, deterministic from end to end: P-site offset
estimation and CDS-level quantification from transcriptome alignments
(`pitstop quant`), reading-frame evidence for uORFs and other ORFs besides
the CDS (`pitstop orfs` / `pitstop score`), and translation start sites from
harringtonine libraries, with a uORF catalogue from them (`pitstop starts`).

## Why

Ribo-seq quantification needs the isoform-length term accounted for, a fit
that actually converges from a fixed, deterministic start, and P-site offsets
estimated per read length from the data rather than assumed. Counting evidence for a
uORF additionally needs a reading frame that isn't thrown away before it's
checked, and a null that accounts for the host CDS's own off-frame noise.
pitstop does all of this, benchmarked against real runs rather than
simulation alone.

## Get started

- **Install and run**: see the [README](https://github.com/katarinagresova/pitstop#readme)
  for setup, the CLI, and output file formats.
- **[Method: quant](method.md)**: how offsets, P-site assignment, and
  quantification work, with interactive examples.
- **[Method: ORF counting](orfs.md)**: ORF types, outside components, the
  frame-weighted EM, and the codon lead score.
- **[Method: start sites](starts.md)**: harringtonine against elongation,
  the background, the test, the scan and the catalogue.
- **Using `orfs` and `score`**: a [typical workflow](orfs.md#typical-workflow)
  for two conditions, [how to read the scores](orfs.md#reading-the-results),
  and [what they give on real data](orfs.md#validation-on-real-data).

## License

Apache License 2.0, see [LICENSE](https://github.com/katarinagresova/pitstop/blob/main/LICENSE)
and [NOTICE](https://github.com/katarinagresova/pitstop/blob/main/NOTICE).
