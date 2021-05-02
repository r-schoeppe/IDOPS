#!/usr/bin/env python3

import os
import sys
from logging import getLogger

from idops.cli import cli
from idops.genomic_analysis.analyze_synteny import analyze_synteny
from idops.hmm.scan_models import scan_models
from idops.tree.make_trees import make_trees


def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]
    args = cli(argv)
    # TODO: configure logging
    logger = getLogger("")

    # Scan HMM DB
    df_scan = scan_models(args.sequence_file, args)
    logger.info(f"IDOPS: Identified {len(df_scan)} pesticidal sequences.\n")
    # Debug save
    # scan_result.to_pickle("scan_result.pickle")
    df_scan.to_csv(os.path.join(args.output, "idops_scan.tsv"), sep="\t", header=True, index=False)

    # Future continue?
    # df_scan = pd.read_csv(os.path.join(args.output, "idops_scan.tsv"), sep="\t",)

    # Genomic Environment
    if args.synteny:
        analyze_synteny(df_scan, args)

    if args.tree:
        make_trees(df_scan, args)
        df_scan.to_csv(os.path.join(args.output, "idops_scan.tsv"), sep="\t", header=True, index=False)

    return


if __name__ == "__main__":
    main(sys.argv[1:])
