"""The output tables: tab-separated, a header line, no index."""
import pandas as pd


def write_tsv(df, path, na_rep=""):
    """df as a TSV at path. na_rep: the text of a missing value ("" is pandas' default)."""
    df.to_csv(path, sep="\t", index=False, na_rep=na_rep)


def write_stats(stats, path):
    """A dict of stats as a TSV with the columns stat and value, in the order of the dict."""
    write_tsv(pd.DataFrame(list(stats.items()), columns=["stat", "value"]), path)
