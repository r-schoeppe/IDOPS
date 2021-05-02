import subprocess
import sys
import tempfile
from logging import getLogger

import pandas as pd
from Bio import SearchIO

from idops.hmm import idops_hmm_db
from idops.io.sequences import get_input_sequences


def scan_models(files, args):
    logger = getLogger(__name__)
    results = []
    for input_file in files:
        logger.info(f"IDOPS: Scanning {input_file}")
        with tempfile.NamedTemporaryFile(mode="r+") as tblfile:
            try:
                pr = subprocess.run(["hmmscan", "--cut_ga"*(not args.no_tc),  # use ga cutoffs
                                     "--tblout", tblfile.name,  # tbl out
                                     idops_hmm_db,  # hmm
                                     "-"  # Stdinstream
                                     ],
                                    input="".join(get_input_sequences(input_file)), text=True,
                                    check=True, capture_output=True)
            except subprocess.CalledProcessError as err:
                logger = getLogger(__name__)
                logger.error(err.stderr)
                logger.error("hmmscan failed - make sure your environment is properly configured")
                sys.exit(1)

            try:
                hmmscan_results = [r for r in SearchIO.parse(tblfile, format="hmmer3-tab")]
            except ValueError as e:
                pass

        best_hits = [max(query_result, key=lambda hit: hit.bitscore) for query_result in hmmscan_results]
        # Not necessary anymore:  if hit.id != "CryM11" else hit.bitscore*-1

        results.append(
            pd.DataFrame([[input_file, h.query_id, h.id, h.evalue, h.bitscore, h.description] for h in best_hits],
                      columns=["inputfile", "seq_id", "model", "evalue", "bitscore", "model_description"]))
    df = pd.concat(results)
    return df
