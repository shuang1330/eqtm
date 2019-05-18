import os
import pandas as pd


def flip_zscores(to_be_flipped_path, save_filepath, sep="\t"):
    to_be_flipped = pd.read_csv(to_be_flipped_path, sep=sep, index_col=0)

    def flip(zscore, allele):
        if allele == "T":
            return -zscore
        else:
            return zscore

    to_be_flipped["flippedZscore"] = [flip(item[0], item[1]) for item in
                                      zip(to_be_flipped["OverallZScore"],
                                          to_be_flipped["AlleleAssessed"])]
    to_be_flipped.to_csv(save_filepath, sep=sep)
    print("Flipped file saved in %s" % save_filepath)
    return None


if __name__ == "__main__":
    project_rootdir = "/groups/umcg-gcc/tmp03/umcg-sli/"
    to_be_flipped_path = "/groups/umcg-gcc/tmp03/umcg-sli/replication_output/concat_replicate_results/liver_metaTCGA_replicates_in_bios.txt"
    save_filepath = "/groups/umcg-gcc/tmp03/umcg-sli/replication_output/concat_replicate_results/liver_metaTCGA_replicates_in_bios_flipped.txt"
    flip_zscores(to_be_flipped_path, save_filepath, sep="\t")

    raise NotImplementedError
    # flip bios part1 and part2
    to_be_flipped_name = "BIOS_Part2_all_features_added"
    to_be_flipped_path = os.path.join(project_rootdir, "data", "eqtmZscores",
                                      "withExpressionTSSMethyCpgOverlapPromoter",
                                      "BIOS_Part2_all_features_added.txt")
    save_filepath = os.path.join(project_rootdir, "data", "eqtmZscores",
                                 "withExpressionTSSMethyCpgOverlapPromoter",
                                 to_be_flipped_name+"_flipped.txt")
    flip_zscores(to_be_flipped_path, save_filepath, sep="\t")