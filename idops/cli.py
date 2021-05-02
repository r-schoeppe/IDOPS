import os
import argparse
# import shutil
import time
import logging
import sys


def cli(argv):
    parser = argparse.ArgumentParser(prog="idops",
                                     description="IDentification Of Pesticidal Sequences")

    # Mandatory
    parser.add_argument("sequence_file", nargs="+",
                        help="File(s) containing input sequences. "
                             "Supported formats: *.faa [Protein], *.fasta [Protein], *.gbk [DNA with CDS Features]")
    parser.add_argument("-v", "--verbosity", action="count", help="increase output verbosity")
    parser.add_argument("-o", "--output", type=str, default="IDOPS_%DATE_%TIME",
                        help="Output directory, default 'IDOPS_{DATE_TIME}'")
    # parser.add_argument("-f", "--force", dest="force", action="store_true",
    #                     help="Force overwrite output directory."
    #                          "Careful (!) the previous output folder will be deleted and lost!")

    # Optional

    parser.add_argument("-d", "--disable_tc", dest="no_tc", action="store_true",
                        help="No cutoff for hmmscan")
    parser.add_argument("-i", "--identify-conserved-env", dest="synteny", action="store_true",
                        help="Analyze genomic environment of hits (requires properly configured prokka) "
                             "and plot Easyfig")
    parser.add_argument("-c", "--cluster-cutoff", dest="cutoff", type=float, metavar='value', default=0.6,
                        help="Uses the value as distance cutoff for annotation distance single linkage clustering:"
                             "1 (one cluster) >= value >= 0 (no clusters) , default 0.6")
    parser.add_argument("-k", "--keep-annotations", action="store_true", dest="keep_annotations",
                        help="switches off the reannotation with Prokka. "
                             "BEWARE: for proper genomic analysis it is crucial "
                             "that the annotations are done by the same tool!")
    parser.add_argument("-t", "--tree", action="store_true",
                        help="Create phylogenetic tree for each hit with the 10 closest sequences of the corresponding protein group. "
                             "Adds the column 'NearestNeighbor' in summary table.")

    # Parse
    args = parser.parse_args(argv)
    if args.output == "IDOPS_%DATE_%TIME":
        s_time = time.strftime("%y%m%d_%H%M%S")
        args.output = f"IDOPS_{s_time}"
    # if os.path.exists(args.output):
    #     if args.force:
    #         shutil.rmtree(args.output)
    os.mkdir(args.output)

    # init logging
    _init_logging(args)

    if args.synteny and not all([s[-4:] == ".gbk" for s in args.sequence_file]):
        raise ValueError("Analysis of genomic environment requires a .gbk file!")

    return args


def _init_logging(args):
    # verbosity
    lvl = logging.WARNING
    if isinstance(args.verbosity, int):
        if args.verbosity == 1:
            lvl = logging.INFO
        if args.verbosity >= 2:
            lvl = logging.DEBUG

    # BasicConfig has to be DEBUG to catch everything
    fmt = '%(levelname)-11s:%(message)s'
    logging.basicConfig(filename='/dev/null', filemode='w', level=logging.DEBUG, format=fmt)
    logger = logging.getLogger('')

    fileHandlerINFO = logging.FileHandler(filename=os.path.join(f'{args.output}', 'idops.log'), mode="w")
    fileHandlerINFO.setLevel(logging.INFO)
    fileHandlerINFO.formatter = logger.handlers[0].formatter
    logger.addHandler(fileHandlerINFO)

    streamHandler = logging.StreamHandler(stream=sys.stdout)
    streamHandler.setLevel(lvl)
    streamHandler.formatter = logger.handlers[0].formatter
    logger.addHandler(streamHandler)
    return
