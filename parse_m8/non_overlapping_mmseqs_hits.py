"""
Parse mmseqs hits and report those that don't overlap the hit with the lowest e-value
"""

import os
import sys
import argparse
from hackathon_lib import stream_blast_results

__author__ = 'Rob Edwards'


def non_overlapping_mmseqs(m8file, mineval, verbose=False):
    """
    Find non-overlapping hits
    """

    hits = []
    last_query =  None
    positions = set()
    for r in stream_blast_results(m8file, verbose):
        if r.evalue > mineval:
            continue
        if r.query != last_query:
            if not hits:
                print(f"There were no hits for {last_query}", file=sys.stderr)
            for h in hits:
                print(str(h))
            hits = []
            positions = set()
            last_query = r.query
        keep = True
        qsta = r.query_start
        qend = r.query_end
        if qsta > qend:
            (qend, qsta) = (qsta, qend)
        for i in range(qsta, qend + 1):
            if i in positions:
                keep = False
        if keep:
            for i in range(qsta, qend + 1):
                positions.add(i)
            hits.append(r)
    print(list(map(str, hits)))







if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Parse m8 format mmseqs results')
    parser.add_argument('-f', help='mmseqs m8 format output file', required=True)
    parser.add_argument('-e', help=f'minimum e value [Default: %(default)f]', default=1e-5, type=float)
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    print(f"E valu: {args.e}", file=sys.stderr)

    non_overlapping_mmseqs(args.f, args.e, args.v)
