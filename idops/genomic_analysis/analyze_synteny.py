import difflib
import os
import shutil
import subprocess
import sys
from glob import glob
from logging import getLogger

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import SeqIO
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
from scipy.spatial.distance import squareform

from idops.genomic_analysis.colorify_gbk import colorify_gbks
from idops.hmm import idops_hmm_db

surround_distance = 5000


def analyze_synteny(df_scan, args):
    # proximity annotation
    surround_snippet_records = get_hit_environment(df_scan)
    estimate_syntheny(surround_snippet_records, args)
    return


def get_hit_environment(df_scan):
    surround_snippet_records = {}
    for inputfile in df_scan.loc[:, "inputfile"].unique():
        # Fix due to multiple gb records in one gbk file
        # missing_ids = []
        for gb_record in SeqIO.parse(inputfile, "genbank"):
            feature_indices = index_genbank_features(gb_record, "CDS", "locus_tag")

            for _, hit in df_scan.loc[df_scan["inputfile"] == inputfile].iterrows():
                hit_id = hit["seq_id"]
                #
                # # sanity fixture
                if hit_id not in feature_indices.keys():
                    # missing_ids.append(hit_id)
                    continue
                # elif hit_id in missing_ids:
                #     missing_ids.remove(hit_id)

                feature_index = feature_indices[hit_id]
                feature = gb_record.features[feature_index]

                # Overwrite annotation?
                # feature.qualifiers["product"] = [hit["model"]]

                snippet_end = min([max(feature.location) + surround_distance, len(gb_record)])
                snippet_start = max([0, min(feature.location) - surround_distance])

                gb_snippet = gb_record[snippet_start: snippet_end]

                # unify toxin direction in figure
                if feature.location.strand == -1:
                    gb_snippet = gb_snippet.reverse_complement()
                surround_snippet_records[f"{inputfile}:{hit_id}"] = gb_snippet

        # if len(missing_ids):
        #     getLogger(__name__).error("Missed some ids: %r", missing_ids)
    return surround_snippet_records


def estimate_syntheny(surround_snippet_records, args):
    annotation_dir = os.path.join(args.output, "annotation_tmp/")

    # annotate snippets and get snippet label with its annotations
    labels, l_products = annotate_snippets(surround_snippet_records, annotation_dir, args)

    # Create Distance Matrix (based on "annotation distance")
    dmat_snippets_features = pd.DataFrame(_calc_distance_matrix(l_products, method=pairwise_annotation_distance),
                                          index=labels, columns=labels)

    # Single Linkages (NJ)
    linkages = linkage(squareform(dmat_snippets_features.values), method="single", )
    # Cut off on accumulated distance
    cluster_dist_cutoff = args.cutoff
    # Cluster ids of labels
    cluster_ids = fcluster(linkages, t=cluster_dist_cutoff, criterion="distance")

    # Create Dendogram figure
    result_dendo = create_synteny_dendogram(labels, linkages, cluster_dist_cutoff, args)

    # Easyfigs
    create_easyfig_figures_on_synteny_clusters(labels, cluster_ids, result_dendo, annotation_dir, args)
    return


def annotate_snippets(surround_snippet_records, annotation_dir, args):

    l_products = []
    labels = []
    os.makedirs(annotation_dir, exist_ok=False)  # True?
    for i, (desc, surround_snippet_record) in enumerate(surround_snippet_records.items()):
        label = desc.rsplit(":")[1]
        basename = os.path.join(f"{annotation_dir}", f"{i}_{label}")

        # Direct genbank output
        if args.keep_annotations:
            with open(basename + ".gbk", "w") as f:
                f.write(surround_snippet_record.format("genbank"))
        # Annotate with prokka
        else:
            with open(basename + ".fa", "w") as f:
                f.write(surround_snippet_record.format("fasta"))
            _prokka_annotation(prefix=f"{i}_{label}", basename=basename)

        l_products.append(get_gbk_cds_products(basename + ".gbk"))
        labels.append(label)
    return labels, l_products


def _prokka_annotation(prefix, basename):
    logger = getLogger(__name__)

    cmd = ["prokka",
           "--hmms", idops_hmm_db,
           "--outdir", basename,
           "--prefix", prefix,
           basename + ".fa"]
    logger.info("> " + " ".join(cmd))
    try:
        subprocess.run(cmd, check=True, capture_output=True)
    except subprocess.CalledProcessError as err:
        logger.error("> " + " ".join(cmd))
        logger.error(err.stderr)
        logger.error("prokka failed - make sure your environment is properly configured")
        sys.exit(1)

    shutil.move(glob(os.path.join(basename, "*.gbk"))[0], basename + ".gbk")
    shutil.rmtree(basename)
    return


def _calc_distance_matrix(v, method):
    n = len(v)
    dmat = np.zeros((n, n), )
    for i in range(n):
        for j in range(n):
            if j >= i:
                dmat[i, j] = dmat[j, i] = method(v[i], v[j])
    return dmat


def create_synteny_dendogram(labels, linkages, cluster_dist_cutoff, args):
    plt.figure(figsize=(6, 9))
    result_dendo = dendrogram(linkages, get_leaves=True, labels=np.array(labels), color_threshold=cluster_dist_cutoff,
                              orientation="left", leaf_rotation=0, distance_sort="ascending")
    plt.axvline(x=cluster_dist_cutoff, c='grey', lw=1, linestyle='dashed')
    plt.title("Annotation Distance of hit environment")
    plt.tight_layout()
    plt.savefig(os.path.join(args.output, "dendogram_annotation_distance.svg"))
    return result_dendo


def create_easyfig_figures_on_synteny_clusters(labels, cluster_ids, result_dendo, annotation_dir, args):
    logger = getLogger(__name__)

    df_clusters = pd.DataFrame(np.array([labels, cluster_ids]).T,
                               index=range(len(labels)), columns=["label", "cluster"]).sort_index()

    for i_cl, df in df_clusters.groupby("cluster"):
        if len(df.index) <= 1:
            continue

        # sort by dendogram order/joining
        index_sorted = sorted(df.index, key=lambda i: result_dendo["ivl"].index(df.loc[i, "label"]))

        selected_gbks = [glob(os.path.join(f"{annotation_dir}", f"{j}_*.gbk"))[0] for j in index_sorted]

        # include colors for features
        colorify_gbks(selected_gbks, add_gene_to_product=True)
        # generate figure - run Easyfic
        cmd = ["Easyfig_idops.py",
               # Svg output
               "-o", os.path.join(args.output, f"{i_cl}.svg"), "-svg",
               # Alignment of genomes
               "-aln", "best",
               # Blast legend + red
               "-f1", "T",  # "-blast_col", "red",
               # legend settings
               "-legend", "double", "-leg_name", "product",
               # input files
               *selected_gbks]
        logger.info("> " + " ".join(cmd))
        try:
            subprocess.run(cmd, check=True, capture_output=True)
        except subprocess.CalledProcessError as err:
            logger.error("> " + " ".join(cmd))
            logger.error(err.stderr)
            logger.error("Easyfig failed!")
            sys.exit(1)
    return


def index_genbank_features(gb_record, feature_type, qualifier):
    answer = dict()
    for (index, feature) in enumerate(gb_record.features):
        if feature.type == feature_type:
            if qualifier in feature.qualifiers:
                # There should only be one locus_tag per feature, but there
                # are usually several db_xref entries
                for value in feature.qualifiers[qualifier]:
                    if value in answer:
                        getLogger(__name__).warning("Duplicate key %s for %s features %i and %i" %
                                                    (value, feature_type, answer[value], index))
                    else:
                        answer[value] = index
    return answer


def get_gbk_cds_products(gbk_file):
    return [f.qualifiers["product"][0]
            for f in SeqIO.read(gbk_file, "genbank").features
            if f.type == "CDS"]


def pairwise_annotation_distance(a, b):
    seq_matcher = difflib.SequenceMatcher(a=a,
                                          b=b,
                                          # isjunk=lambda s: s == "hypothetical protein",
                                          )
    n_matches = sum(triple[-1] for triple in seq_matcher.get_matching_blocks())
    seq_ratio = n_matches / max([len(a), len(b), 1])

    dist = 1 - seq_ratio
    return dist
