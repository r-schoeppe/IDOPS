import os
import subprocess
import tempfile
from logging import getLogger

from Bio import AlignIO
from Bio import Phylo
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor

import idops
from idops.io.sequences import get_selected_sequences


def make_trees(df_scan, args):
    logger = getLogger(__name__)
    logger.info("Make trees from profile alignments")

    os.makedirs(os.path.join(args.output, "trees"), exist_ok=True)
    filename = os.path.join(args.output, "trees" , "{id}_{model}.tree")

    for model, df_group in df_scan.groupby("model"):
        logger.info("Do alignments with HMM " + model)
        aln = hmmalign_hits(model, df_group)
        idops_hits = set(df_group["seq_id"])

        # Fix Alphabet
        for rec in aln:
            rec.seq = rec.seq.upper()

        calculator = DistanceCalculator('blosum45')

        for rec in aln:
            if rec.id in idops_hits:
                logger.info("")
                neighbors = get_alignment_with_closest_sequences(rec, aln, filter=idops_hits)
                df_scan.loc[df_scan["seq_id"] == rec.id, "NearestNeighbor"] = neighbors[0].id
                constructor = DistanceTreeConstructor(calculator, 'nj')
                tree = constructor.build_tree(neighbors)

                with open(filename.format(id=rec.id, model=model), "w") as f:
                    f.write(tree.format("newick"))
                with open(filename.format(id=rec.id, model=model) + "_ascii", "w") as f:
                    Phylo.draw_ascii(tree, file=f)

        # #dm = calculator.get_distance(aln)
        # #logger.info("Did distances of group " + model)
        # constructor = DistanceTreeConstructor(calculator, 'nj')
        # tree = constructor.build_tree(aln)
        # logger.info("Did (distances?) + tree of group " + model)
        #
        # with open(os.path.join(args.output, f"{model}.tree"), "w") as f:
        #     f.write(tree.format("newick"))
        # with open(os.path.join(args.output, f"{model}.tree_ascii"), "w") as f:
        #     Phylo.draw_ascii(tree, file=f)

        # dm.matrix
        # np.array(dm.matrix)
        # import numpy as np
        # np.array(dm.matrix)
        # type(dm.matrix)
        # np.zeros(shape=(len(dm.matrix), len(dm.matrix)))
        # m_dm = np.zeros(shape=(len(dm.matrix), len(dm.matrix)))
        # for i, r in enumerate(dm.matrix):
        #     for j, c in enumerate(r):
        #         m_dm[i, j] = m_dm[j, i] = c
        # m_
        # m_dm

    return


def hmmalign_hits(model, df_group):
    dict_files_seqids = {}
    for input_file, _df in df_group.groupby("inputfile"):
        dict_files_seqids[input_file] = list(_df["seq_id"])

    with tempfile.NamedTemporaryFile(mode="r+") as alnfile:
        try:
            subprocess.run(["hmmalign",
                             "-o", alnfile.name,
                            os.path.join(os.path.dirname(idops.__file__), "ressources", "hmms", model + ".hmm"),
                            "-"  # Stdinstream
                            ],
                           input=("".join(get_selected_sequences(dict_files_seqids)) +
                                  open(os.path.join(os.path.dirname(idops.__file__), "ressources", "seqdb",
                                                    model.rsplit("M")[0] + ".faa")).read()
                                  ),
                           text=True,
                           check=True, capture_output=True)
        except subprocess.CalledProcessError as err:
            logger = getLogger(__name__)
            logger.error(err.stderr)
            logger.error("hmmalign failed")
            raise
        aln = AlignIO.read(alnfile.name, "stockholm")
    return aln


def get_alignment_with_closest_sequences(rec, aln, filter=None):
    calculator = DistanceCalculator('blosum45')

    distances = []
    for i, seq in enumerate(aln):
        if isinstance(filter, set) and seq.id in filter:
            continue
        distances.append((calculator._pairwise(rec, seq), i))

    return MultipleSeqAlignment([aln[i] for _, i in sorted(distances)[:10]] + [rec])
