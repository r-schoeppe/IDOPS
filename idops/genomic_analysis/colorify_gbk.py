import sys
import numpy as np
from Bio import SeqIO
from matplotlib import cm


def colorify_gbks(filenames, add_gene_to_product=False):
    features = {}
    for gbk_file in filenames:
        for gb_record in SeqIO.parse(gbk_file, "genbank"):
            for feature in _iter_cds_features(gb_record):
                if feature.type != "CDS":
                    continue
                if "product" in feature.qualifiers and not feature.qualifiers["product"][0] == "hypothetical protein":
                    feat_name = feature.qualifiers["product"][0]
                elif "gene" in feature.qualifiers:
                    feat_name = feature.qualifiers["gene"][0]
                else:
                    continue
                features[feat_name] = features.get(feat_name, 0) + 1 + 10 * ("IDOPS" in feat_name)  # last part for idops prio

    features = sorted(features.items(), key=lambda i: i[1], reverse=True)
    color_map = dict(zip([i[0] for i in features[:20]],
                         np.array(cm.get_cmap("tab20").colors)))
    #print(features)
    # else: if n <= 20:
    #     color_map = dict(zip(features, cm.get_cmap("gist_rainbow")(np.linspace(0.0, 1.0, n))[:, :3]))

    for gbk_file in filenames:
        records = []
        for gb_record in SeqIO.parse(gbk_file, "genbank"):
            records.append(gb_record)
            for feature in _iter_cds_features(gb_record):
                if "product" in feature.qualifiers and not feature.qualifiers["product"][0] == "hypothetical protein":
                    color = color_map.get(feature.qualifiers["product"][0], None)
                elif "gene" in feature.qualifiers:
                    color = color_map.get(feature.qualifiers["gene"][0], None)
                    if add_gene_to_product and "product" in feature.qualifiers:
                        feature.qualifiers["product"][0] += f" GN={feature.qualifiers['gene'][0]}"
                else:
                    continue

                # skip ignored features
                if color is None:
                    continue
                feature.qualifiers["color"] = ["%d %d %d" % tuple(np.rint(color*255))]
        with open(gbk_file, "w") as f_gbk:
            for gb_record in records:
                f_gbk.write(gb_record.format("genbank"))
    return


def _iter_cds_features(gb_record):
    for feature in gb_record.features:
        if feature.type == "CDS":
            yield feature


if __name__ == "__main__":
    colorify_gbks(sys.argv[1:])