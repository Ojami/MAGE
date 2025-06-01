"""
bFileReaderDep.py

A Dask‐based re‐implementation of the Java class bFileReaderDep.
All function names and signatures are kept identical to the original,
so your MATLAB toolbox can continue to use them via the py module.

This module always uses dd.read_csv to read tabular files (whether compressed or not).
The CSV parser uses the highly optimized C engine and, when possible,
splits the file into partitions (using the blocksize parameter) so that operations are vectorized and parallelized.
If the delimiter parameter is not provided (set to "NA", None, or empty),
the code auto-detects it from the header.

Note: For gzipped files, blocksize is ignored since gzip is not splittable.
All functions (including filtering, extracting rows, etc.) now use dd.read_csv exclusively.
"""

import sys
import multiprocessing
# Force new processes to use pythonw.exe on Windows to avoid extra console windows.
if sys.executable.lower().endswith("python.exe"):
    new_executable = sys.executable[:-len("python.exe")] + "pythonw.exe"
    multiprocessing.set_executable(new_executable)

import dask.dataframe as dd
import gzip
import re
import dask
import operator
import functools
from dask.distributed import Client, LocalCluster

# Global client variable
client = None

def get_client():
    """Return the global Dask client, creating it if necessary."""
    global client
    if client is None:
        client = Client(LocalCluster(n_workers=32, threads_per_worker=1, local_directory="C:\\Temp"))
    return client

def shutdown_client():
    """Shut down the global Dask client if it exists and set it to None."""
    global client
    if client is not None:
        client.shutdown()
        client.close()
        client = None
        return "Cluster shut down."
    else:
        return "No cluster running."

# -------------
# Helper functions for file reading and processing
# -------------

def _get_compression(filename):
    """Return 'gzip' if the file is gzipped, else None."""
    return 'gzip' if filename.endswith('.gz') else None

def detect_delimiter(filename, skip=0, blocksize=None):
    """
    Detect the delimiter by reading the first non-empty line after skipping.
    Considers common candidates: tab, comma, semicolon, and pipe.
    """
    comp = _get_compression(filename)
    if comp == 'gzip':
        f = gzip.open(filename, 'rt')
    else:
        f = open(filename, 'r')
    for i in range(skip):
        f.readline()
    header = ""
    while not header:
        header = f.readline().strip()
        if header == "":
            break
    f.close()
    if not header:
        return ','
    candidates = ['\t', ',', ';', '|']
    counts = {d: header.count(d) for d in candidates}
    best = max(counts, key=counts.get)
    return best if counts[best] > 0 else ','

def get_delim(filename, delim, skip=0, blocksize=None):
    """
    If delim is "NA", None, or empty, auto-detect the delimiter from the file header.
    """
    if delim == "NA" or delim is None or delim == "":
        return detect_delimiter(filename, skip, blocksize)
    return delim

def _dd_read_csv(filename, delim, skip, blocksize=None):
    """
    Read a file (compressed or uncompressed) as a Dask DataFrame.
    The blocksize parameter is used only for uncompressed files.
    """
    return dd.read_csv(filename, delimiter=delim, engine='c', dtype=str, skiprows=skip, blocksize=blocksize)

def _dd_combine_rows(df, delim):
    """
    Combine all columns of a DataFrame into a single string for each row using apply().
    """
    combined = df.astype(str).apply(lambda row: delim.join(row.tolist()), axis=1, meta=('x', 'object'))
    return combined

def _combine_dd_rows(df, delim, selected_cols=None):
    """
    Combine rows of a DataFrame into strings.
    If selected_cols is provided, use only those columns.
    """
    if selected_cols is not None and selected_cols[0] is not None:
        df = df.iloc[:, selected_cols]
    return _dd_combine_rows(df, delim)

def _combine_regex(terms):
    """
    Combine a list of regex patterns into one pattern using non-capturing groups.
    """
    return "(?:" + ")|(?:".join(terms) + ")"

def getCustomCols(s, delim, cols):
    """
    Split the string s by delim, select the columns in 'cols', and rejoin them.
    """
    parts = s.split(delim)
    selected = [parts[i] for i in cols if i < len(parts)]
    return delim.join(selected)

# -------------
# The main class
# -------------

class bFileReaderDep:
    """
    A general purpose file-reading tool that supports tabular files (gzipped or uncompressed).
    All functions use dd.read_csv so that the C engine and vectorized operations are used
    to speed up reading and filtering. The interface remains similar to the original Java version.
    """

    @staticmethod
    def getClient():
        """Return the Dask client (creating it if needed)."""
        return get_client()

    @staticmethod
    def shutdownCluster():
        """Shut down the Dask cluster and return a message."""
        return shutdown_client()

    @staticmethod
    def readAll(filename, delim="NA", skip=0, bench="sp", blocksize=None):
        delim = get_delim(filename, delim, skip, blocksize)
        df = _dd_read_csv(filename, delim, skip, blocksize)
        series = _dd_combine_rows(df, delim)
        lines = series.compute().tolist()
        return delim.join(lines)

    @staticmethod
    def readAll2(filename, skip=0, bench="sp", blocksize=None):
        delim = get_delim(filename, "NA", skip, blocksize)
        df = _dd_read_csv(filename, delim, skip, blocksize)
        series = _dd_combine_rows(df, delim)
        return series.compute().tolist()

    @staticmethod
    def readHeader(filename, skip=0, blocksize=None):
        delim = get_delim(filename, "NA", skip, blocksize)
        df = dd.read_csv(filename, delimiter=delim, engine='c', dtype=str,
                         skiprows=skip, nrows=1, blocksize=blocksize)
        head = df.head(1)
        if len(head) == 0:
            return ""
        row = delim.join(head.iloc[0].astype(str).tolist())
        return row

    @staticmethod
    def lineCount2(filename, bench="sp", blocksize=None):
        delim = get_delim(filename, "NA", 0, blocksize)
        df = dd.read_csv(filename, delimiter=delim, engine='c', dtype=str, blocksize=blocksize)
        return int(df.shape[0].compute())

    @staticmethod
    def lineCount(filename):
        if filename.endswith('.gz'):
            f = gzip.open(filename, 'rt')
        else:
            f = open(filename, 'r')
        count = 0
        chunk_size = 1024
        while True:
            chunk = f.read(chunk_size)
            if not chunk:
                break
            count += chunk.count('\n')
        f.close()
        return count if count > 0 else 1

    @staticmethod
    def compare(filename, terms, delim="NA", skip=0, bench="sp", retrunCols=None, blocksize=None):
        delim = get_delim(filename, delim, skip, blocksize)
        pattern = _combine_regex(terms)
        df = _dd_read_csv(filename, delim, skip, blocksize)
        series = _dd_combine_rows(df, delim)
        filtered = series[series.str.contains(pattern, regex=True, na=False)]
        lines = filtered.compute().tolist()
        return delim.join(lines)

    @staticmethod
    def compareToCols(filename, terms, delim="NA", col=None, skip=0, bench="sp", retrunCols=None, blocksize=None):
        delim = get_delim(filename, delim, skip, blocksize)
        pattern = _combine_regex(terms)
        df = _dd_read_csv(filename, delim, skip, blocksize)
        conditions = [df.iloc[:, i].str.contains(pattern, regex=True, na=False) for i in col]
        overall = functools.reduce(operator.or_, conditions)
        filtered = df[overall]
        if retrunCols is not None and retrunCols[0] is not None:
            filtered = filtered.iloc[:, retrunCols]
        series = _dd_combine_rows(filtered, delim)
        lines = series.compute().tolist()
        return delim.join(lines)

    @staticmethod
    def compareToJoinedCols(filename, terms, delim="NA", col=None, skip=0, bench="sp", retrunCols=None, chunk=None, blocksize=None):
        delim = get_delim(filename, delim, skip, blocksize)
        pattern = _combine_regex(terms)
        df = _dd_read_csv(filename, delim, skip, blocksize)
        df['joined'] = df.iloc[:, col].astype(str).apply(lambda row: "_".join(row.tolist()), axis=1, meta=('x','object'))
        filtered = df[df['joined'].str.contains(pattern, regex=True, na=False)]
        if retrunCols is not None and retrunCols[0] is not None:
            filtered = filtered.iloc[:, retrunCols]
        series = _dd_combine_rows(filtered, delim)
        lines = series.compute().tolist()
        return delim.join(lines)

    @staticmethod
    def filterCol(filename, delim="NA", col=None, term=None, cond=None, skip=0, bench="sp", retrunCols=None, blocksize=None):
        delim = get_delim(filename, delim, skip, blocksize)
        ops = {"gt": operator.gt, "lt": operator.lt, "eq": operator.eq, "ge": operator.ge, "le": operator.le}
        df = _dd_read_csv(filename, delim, skip, blocksize)
        conds = []
        for i, op_key in enumerate(term):
            col_series = df.iloc[:, col[i]].astype(float)
            conds.append(ops[op_key](col_series, cond[i]))
        overall = functools.reduce(operator.and_, conds)
        filtered = df[overall]
        if retrunCols is not None and retrunCols[0] is not None:
            filtered = filtered.iloc[:, retrunCols]
        series = _dd_combine_rows(filtered, delim)
        lines = series.compute().tolist()
        return delim.join(lines)

    @staticmethod
    def compareFilterCol(filename, skip=0, delim="NA", patt=None, colT=None, operation=None, colN=None, cond=None, bench="sp", retrunCols=None, blocksize=None):
        delim = get_delim(filename, delim, skip, blocksize)
        pattern = _combine_regex(patt)
        ops = {"gt": operator.gt, "lt": operator.lt, "eq": operator.eq, "ge": operator.ge, "le": operator.le}
        df = _dd_read_csv(filename, delim, skip, blocksize)
        num_conds = []
        for i, op_key in enumerate(operation):
            num_conds.append(ops[op_key](df.iloc[:, colN[i]].astype(float), cond[i]))
        overall_num = functools.reduce(operator.and_, num_conds)
        text_conds = [df.iloc[:, j].str.contains(pattern, regex=True, na=False) for j in colT]
        overall_text = functools.reduce(operator.or_, text_conds)
        filtered = df[overall_num & overall_text]
        if retrunCols is not None and retrunCols[0] is not None:
            filtered = filtered.iloc[:, retrunCols]
        series = _dd_combine_rows(filtered, delim)
        lines = series.compute().tolist()
        return delim.join(lines)

    @staticmethod
    def multiCompareFilterCol(filename, skip=0, delim="NA", patt=None, colT=None, operation=None, colN=None, cond=None, bench="sp", retrunCols=None, blocksize=None):
        delim = get_delim(filename, delim, skip, blocksize)
        ops = {"gt": operator.gt, "lt": operator.lt, "eq": operator.eq, "ge": operator.ge, "le": operator.le}
        df = _dd_read_csv(filename, delim, skip, blocksize)
        num_conds = []
        for i, op_key in enumerate(operation):
            num_conds.append(ops[op_key](df.iloc[:, colN[i]].astype(float), cond[i]))
        overall_num = functools.reduce(operator.and_, num_conds)
        text_conds = [df.iloc[:, colT[i]].str.contains(patt[i], regex=True, na=False) for i in range(len(patt))]
        overall_text = functools.reduce(operator.and_, text_conds)
        filtered = df[overall_num & overall_text]
        if retrunCols is not None and retrunCols[0] is not None:
            filtered = filtered.iloc[:, retrunCols]
        series = _dd_combine_rows(filtered, delim)
        lines = series.compute().tolist()
        return delim.join(lines)

    @staticmethod
    def multiCompareToCols(filename, skip=0, delim="NA", patt=None, colT=None, bench="sp", retrunCols=None, blocksize=None):
        delim = get_delim(filename, delim, skip, blocksize)
        df = _dd_read_csv(filename, delim, skip, blocksize)
        text_conds = [df.iloc[:, colT[i]].str.contains(patt[i], regex=True, na=False) for i in range(len(patt))]
        overall = functools.reduce(operator.and_, text_conds)
        filtered = df[overall]
        if retrunCols is not None and retrunCols[0] is not None:
            filtered = filtered.iloc[:, retrunCols]
        series = _dd_combine_rows(filtered, delim)
        lines = series.compute().tolist()
        return delim.join(lines)

    @staticmethod
    def getColumn(filename, delim="NA", cols=None, skip=0, bench="sp", rows=None, chunk=None, blocksize=None):
        delim = get_delim(filename, delim, skip, blocksize)
        df = _dd_read_csv(filename, delim, skip, blocksize)
        if rows is not None and len(rows) > 0 and rows[0] is not None:
            df = df.iloc[[r-1 for r in rows if r-1 < int(df.shape[0].compute())]]
        if cols is not None:
            df = df.iloc[:, cols]
        series = _dd_combine_rows(df, delim)
        return delim.join(series.compute().tolist())

    @staticmethod
    def getColumnBuffer(filename, delim="NA", col=None, skip=0, bench="sp"):
        delim = get_delim(filename, delim, skip)
        result_lines = []
        if filename.endswith('.gz'):
            f = gzip.open(filename, 'rt')
        else:
            f = open(filename, 'r')
        lineCtr = 0;
        for line in f:
            if lineCtr < skip:
                lineCtr += 1;
                continue
            line = line.rstrip('\n')
            parts = line.split(delim)
            selected = [parts[i] for i in col if i < len(parts)]
            result_lines.append(delim.join(selected))
            lineCtr += 1
        f.close()
        return delim.join(result_lines)

    @staticmethod
    def getRowColumn(filename, delim="NA", cols=None, rows=None, skip=0, bench="sp", chunk=None, blocksize=None):
        delim = get_delim(filename, delim, skip, blocksize)
        df = _dd_read_csv(filename, delim, skip, blocksize)
        if rows is not None and len(rows) > 0:
            df = df.iloc[[r-1 for r in rows if r-1 < int(df.shape[0].compute())]]
        if cols is not None:
            df = df.iloc[:, cols]
        series = _dd_combine_rows(df, delim)
        return delim.join(series.compute().tolist())

# Expose getCustomCols at module level
import sys
sys.modules[__name__].getCustomCols = getCustomCols

# New methods for cluster management
@staticmethod
def getClient():
    return get_client()

@staticmethod
def shutdownCluster():
    return shutdown_client()

bFileReaderDep.getClient = getClient
bFileReaderDep.shutdownCluster = shutdownCluster
