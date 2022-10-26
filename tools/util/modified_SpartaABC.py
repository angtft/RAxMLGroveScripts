#!/usr/bin/env python3
import collections
import copy
import math
import os
import statistics
import subprocess
import sys

import tools.util.msa_parser as msa_parser
#import msa_parser

"""
The functions in this file contain modifications (or extensions) of distances used in the improved ("differentiating 
insertions from deletions") SpartaABC tool by Loewenthal et al. (https://doi.org/10.1093/molbev/msab266).

We converted parts of their C++ code to Python and made some modifications.
"""


key_list = [
        "avg_indel_len",
        "alignment_len",
        "msa_max_len",
        "msa_min_len",
        "total_num_of_indels",
        "num_of_indels_of_len_one",
        "num_of_indels_of_len_two",
        "num_of_indels_of_len_three",
        "num_of_indels_of_len_at_least_four",
        "avg_unique_indel_len",
        "total_num_of_unique_indels",
        "num_of_indels_of_len_one_in_one_pos",
        "num_of_indels_of_len_one_in_two_pos",
        "num_of_indels_of_len_one_in_n_minus_1_pos",
        "num_of_indels_of_len_two_in_one_pos",
        "num_of_indels_of_len_two_in_two_pos",
        "num_of_indels_of_len_two_in_n_minus_1_pos",
        "num_of_indels_of_len_three_in_one_pos",
        "num_of_indels_of_len_three_in_two_pos",
        "num_of_indels_of_len_three_in_n_minus_1_pos",
        "num_of_indels_of_len_at_least_four_in_one_pos",
        "num_of_indels_of_len_at_least_four_in_two_pos",
        "num_of_indels_of_len_at_least_four_in_n_minus_1_pos",
        "num_of_msa_pos_with_0_gaps",
        "num_of_msa_pos_with_1_gaps",
        "num_of_msa_pos_with_2_gaps",
        "num_of_msa_pos_with_n_minus_1_gaps"
    ]
extended_key_list = copy.deepcopy(key_list)
extended_key_list.extend([
    "num_patterns",
    "max_pattern_weight",
    "avg_pattern_weight"
])

key_list_norm_groups = {
    "num_sites": [
        # "avg_indel_len",
        "alignment_len",
        "msa_max_len",
        "msa_min_len",
        # "total_num_of_indels",
        "num_of_indels_of_len_one",
        "num_of_indels_of_len_two",
        "num_of_indels_of_len_three",
        "num_of_indels_of_len_at_least_four",
        # "avg_unique_indel_len",
        # "total_num_of_unique_indels",
        "num_of_indels_of_len_one_in_one_pos",
        "num_of_indels_of_len_one_in_two_pos",
        "num_of_indels_of_len_one_in_n_minus_1_pos",
        "num_of_indels_of_len_two_in_one_pos",
        "num_of_indels_of_len_two_in_two_pos",
        "num_of_indels_of_len_two_in_n_minus_1_pos",
        "num_of_indels_of_len_three_in_one_pos",
        "num_of_indels_of_len_three_in_two_pos",
        "num_of_indels_of_len_three_in_n_minus_1_pos",
        "num_of_indels_of_len_at_least_four_in_one_pos",
        "num_of_indels_of_len_at_least_four_in_two_pos",
        "num_of_indels_of_len_at_least_four_in_n_minus_1_pos",
        "num_of_msa_pos_with_0_gaps",
        "num_of_msa_pos_with_1_gaps",
        "num_of_msa_pos_with_2_gaps",
        "num_of_msa_pos_with_n_minus_1_gaps",

        "max_pattern_weight",
    ],
    "num_patterns": [
        # extended features
        "num_patterns",
        "avg_pattern_weight",
    ],
    "sites_times_taxa": [
        "total_num_of_indels",
        "total_num_of_unique_indels",
    ],
    "normalized": [
        "avg_indel_len",
        "avg_unique_indel_len",
    ]
}


def _create_unique_indel_map(sequences):
    num_sequences = len(sequences)
    msa_len = len(sequences[0].sequence)

    unique_indel_map = {}
    for j in range(num_sequences):
        seq = sequences[j].sequence

        previous_is_indel = 0
        curr_start_indel_point = -1
        curr_end_indel_point = -1

        for i in range(msa_len):
            if seq[i] == "-" and previous_is_indel == 0:
                previous_is_indel = 1
                curr_start_indel_point = i
                curr_end_indel_point = i
            elif seq[i] == "-" and previous_is_indel == 1:
                curr_end_indel_point += 1
            else:
                if curr_start_indel_point == -1:
                    previous_is_indel = 0
                    continue

                curr_pair = (curr_start_indel_point, curr_end_indel_point)
                curr_len = curr_end_indel_point - curr_start_indel_point + 1
                if curr_pair not in unique_indel_map:
                    curr_value = [curr_len, 1]
                    unique_indel_map[curr_pair] = curr_value
                else:
                    unique_indel_map[curr_pair][1] += 1
                previous_is_indel = 0
                curr_start_indel_point = -1
                curr_end_indel_point = -1
        if curr_start_indel_point != -1:
            curr_pair = (curr_start_indel_point, curr_end_indel_point)
            curr_len = curr_end_indel_point - curr_start_indel_point + 1
            if curr_pair not in unique_indel_map:
                curr_value = [curr_len, 1]
                unique_indel_map[curr_pair] = curr_value
            else:
                unique_indel_map[curr_pair][1] += 1
    return unique_indel_map


def _calc_gap_features(sequences_):
    sequences = copy.deepcopy(sequences_)
    """
    max_len = 0
    delete_sites = []
    for i in range(len(sequences[0].sequence)):
        num_gaps = 0

        for j in range(len(sequences)):
            print(len(sequences[j].sequence))
            c = sequences[j].sequence[i]
            if c == "-":
                num_gaps += 1
        if num_gaps == len(sequences):
            for j in range(len(sequences)):
                sequences[j].sequence = sequences[j].sequence[:i] + sequences[j].sequence[i+1:]"""

    num_seqs = len(sequences)
    msa_len = len(sequences[0].sequence)
    features = {}
    for key in key_list:
        features[key] = 0




    indel_counter = [0 for _ in range(msa_len)]
    for i in range(num_seqs):
        for j in range(msa_len):
            if sequences[i].sequence[j] == "-": indel_counter[j] += 1

    aligned_seqs = sequences
    for i in range(num_seqs):
        j = 0
        seq_end = len(aligned_seqs[i].sequence)
        while j < seq_end:
            place_in_indel_counter = j
            try:
                if indel_counter[place_in_indel_counter] == num_seqs:
                    aligned_seqs[i].sequence = aligned_seqs[i].sequence[:j] + aligned_seqs[i].sequence[j+ 1:]
                if j != seq_end: j += 1
            except Exception as e:
                print(indel_counter)
                print(len(indel_counter))
                print(len(aligned_seqs[i].sequence))
                print(msa_len)
                print(place_in_indel_counter)
                raise e

    for i in range(len(indel_counter)):
        if indel_counter[i] == 0:
            features["num_of_msa_pos_with_0_gaps"] += 1
        elif indel_counter[i] == 1:
            features["num_of_msa_pos_with_1_gaps"] += 1
        elif indel_counter[i] == 2:
            features["num_of_msa_pos_with_2_gaps"] += 1
        elif indel_counter[i] == num_seqs - 1:
            features["num_of_msa_pos_with_n_minus_1_gaps"] += 1




    unique_indel_map = _create_unique_indel_map(sequences)

    num_seqs = len(sequences)
    total_num_of_gap_chars = 0
    total_num_of_unique_gap_chars = 0

    msa_len = len(sequences[0].sequence)
    features["alignment_len"] = msa_len

    for gap_indices in unique_indel_map:
        length_and_count = unique_indel_map[gap_indices]

        total_num_of_gap_chars += length_and_count[0] * length_and_count[1]
        features["total_num_of_indels"] += length_and_count[1]

        features["total_num_of_unique_indels"] += 1
        total_num_of_unique_gap_chars += length_and_count[0]

        if length_and_count[0] == 1:
            features["num_of_indels_of_len_one"] += length_and_count[1]
            if length_and_count[1] == 1: features["num_of_indels_of_len_one_in_one_pos"] += 1
            if length_and_count[1] == 2: features["num_of_indels_of_len_one_in_two_pos"] += 1
            if length_and_count[1] == num_seqs - 1: features["num_of_indels_of_len_one_in_n_minus_1_pos"] += 1

        if length_and_count[0] == 2:
            features["num_of_indels_of_len_two"] += length_and_count[1]
            if length_and_count[1] == 1: features["num_of_indels_of_len_two_in_one_pos"] += 1
            if length_and_count[1] == 2: features["num_of_indels_of_len_two_in_two_pos"] += 1
            if length_and_count[1] == num_seqs - 1: features["num_of_indels_of_len_two_in_n_minus_1_pos"] += 1

        if length_and_count[0] == 3:
            features["num_of_indels_of_len_three"] += length_and_count[1]
            if length_and_count[1] == 1: features["num_of_indels_of_len_three_in_one_pos"] += 1
            if length_and_count[1] == 2: features["num_of_indels_of_len_three_in_two_pos"] += 1
            if length_and_count[1] == num_seqs - 1: features["num_of_indels_of_len_three_in_n_minus_1_pos"] += 1

        if length_and_count[0] > 3:
            features["num_of_indels_of_len_at_least_four"] += length_and_count[1]
            if length_and_count[1] == 1: features["num_of_indels_of_len_at_least_four_in_one_pos"] += 1
            if length_and_count[1] == 2: features["num_of_indels_of_len_at_least_four_in_two_pos"] += 1
            if length_and_count[1] == num_seqs - 1: features["num_of_indels_of_len_at_least_four_in_n_minus_1_pos"] += 1

    if features["total_num_of_indels"] > 0:
        features["avg_indel_len"] = total_num_of_gap_chars / features["total_num_of_indels"]
        features["avg_unique_indel_len"] = total_num_of_unique_gap_chars / features["total_num_of_unique_indels"]

    max_len = 0
    min_len = msa_len
    for i in range(num_seqs):
        current_len = 0
        for j in range(msa_len):
            if sequences[i].sequence[j] != "-":
                current_len += 1
        if current_len > max_len:
            max_len = current_len
        if current_len < min_len:
            min_len = current_len

    features["msa_max_len"] = max_len
    features["msa_min_len"] = min_len




    return features


def count_patterns(sequences):
    patterns = get_patterns(sequences)
    return len(patterns.keys())


def get_patterns(sequences):
    patterns = collections.defaultdict(lambda: 0)
    for j in range(len(sequences[0].sequence)):
        temp = ""
        for i in range(len(sequences)):
            temp += sequences[i].sequence[j]
        patterns[temp] += 1
    return patterns


def count_gap_proportion(sequences):
    msa_len = len(sequences[0].sequence)
    num_gaps = 0
    for sequence in sequences:
        seq = sequence.sequence
        for c in range(len(seq)):
            if c == "-":
                num_gaps += 1
    return num_gaps / (msa_len * len(sequences))


def get_features_from_sparta(msa_path1, msa_path2):
    command = [
        os.path.join("tools", "SpartaABC", "cpp_code", "SpartaABC"),
        "stat-dist", msa_path1, msa_path2
    ]

    out = b"0"
    try:
        out = subprocess.check_output(command, timeout=60)
    except subprocess.TimeoutExpired:
        return 0
    except Exception as e:
        print(e)
        return 0

    out_lines = out.decode("utf8").split("\n")

    features1 = [float(x) for x in out_lines[0].split()][:-1]
    features2 = [float(x) for x in out_lines[1].split()][:-1]

    return features1, features2


def get_gap_features(sequences):
    feature_dict = _calc_gap_features(sequences)

    feature_vec = [feature_dict[key] for key in key_list]
    return feature_vec


def get_extended_gap_features(sequences):
    feature_dict = _calc_gap_features(sequences)

    patterns = get_patterns(sequences)
    msa_len = len(sequences[0].sequence)
    num_taxa = len(sequences)
    num_patterns = len(patterns.keys())
    patterns_by_sites = num_patterns/msa_len
    pattern_weights = [patterns[key] for key in patterns.keys()]
    max_pattern_weight = max(pattern_weights)
    avg_pattern_weight = statistics.mean(pattern_weights)

    feature_dict["num_patterns"] = num_patterns
    feature_dict["max_pattern_weight"] = max_pattern_weight
    feature_dict["avg_pattern_weight"] = avg_pattern_weight
    feature_dict["num_taxa"] = num_taxa

    feature_vec = [feature_dict[key] for key in key_list]
    feature_vec.append(patterns_by_sites)
    feature_vec.append(max_pattern_weight)
    feature_vec.append(avg_pattern_weight)

    return feature_dict


def get_simple_features(sequences):
    msa_len = len(sequences[0].sequence)
    num_patterns = count_patterns(sequences)
    num_gaps = count_gap_proportion(sequences)

    feature_vec = [
        msa_len,
        num_patterns,
        num_gaps
    ]
    return feature_vec


def msa_blind_dist1(features1, features2, weights=[]):
    # use features1 as reference and normalize accordingly
    msa_len = features1[0]
    features1 = [features1[0] / msa_len, features1[1] / msa_len, features1[2]]
    features2 = [features2[0] / msa_len, features2[1] / msa_len, features2[2]]

    if weights:
        diffs = [weights[i] * abs(features1[i] - features2[i]) for i in range(len(weights))]
    else:
        diffs = [abs(features1[i] - features2[i]) for i in range(len(features1))]
    return sum(diffs)


def msa_blind_dist2(features1, features2, weights=[]):
    # use features1 as reference and normalize accordingly
    msa_len = features1[0]
    features1 = [features1[0] / msa_len, features1[1] / msa_len, features1[2]]
    features2 = [features2[0] / msa_len, features2[1] / msa_len, features2[2]]

    if weights:
        squared_diffs = [weights[i] * ((features1[i] - features2[i]) ** 2) for i in range(len(weights))]
    else:
        squared_diffs = [(features1[i] - features2[i]) ** 2 for i in range(len(features1))]
    return math.sqrt(sum(squared_diffs))


def sparta_dist(features1, features2, weights=[]):
    if weights:
        squared_diffs = [(weights[i] ** 2) * ((features1[i] - features2[i]) ** 2) for i in range(len(weights))]
    else:
        squared_diffs = [((features1[i] - features2[i]) ** 2) for i in range(len(features1))]
    return math.sqrt(sum(squared_diffs))


def sparta_extended_dist(feature_dict1, feature_dict2, weights=[]):
    features1 = []
    features2 = []

    norm_dict = []
    for g in key_list_norm_groups:
        norm_dict.extend([(x, g) for x in key_list_norm_groups[g]])
    norm_dict = dict(norm_dict)

    num_sites = feature_dict1["alignment_len"]
    num_patterns = feature_dict1["num_patterns"]
    sites_times_taxa = num_sites * feature_dict1["num_taxa"]

    for key in extended_key_list:
        if norm_dict[key] == "num_sites":
            features1.append(feature_dict1[key] / num_sites)
            features2.append(feature_dict2[key] / num_sites)
        elif norm_dict[key] == "num_patterns":
            features1.append(feature_dict1[key] / num_patterns)
            features2.append(feature_dict2[key] / num_patterns)
        elif norm_dict[key] == "sites_times_taxa":
            features1.append(feature_dict1[key] / sites_times_taxa)
            features2.append(feature_dict2[key] / sites_times_taxa)
        elif norm_dict[key] == "normalized":
            features1.append(feature_dict1[key])
            features2.append(feature_dict2[key])
        else:
            raise ValueError(f"unknown feature: {key}")

    if not weights:
        # 24 is the number of the gap metrics. we weight site and pattern lengths to make them as important
        # as the gap statistics -> weight sites ~ weight patterns ~ weight of 24 gap stats
        features1[-3] *= 24
        features2[-3] *= 24
        features1[1] *= 24
        features2[1] *= 24
        squared_diffs = [(features1[i] - features2[i]) ** 2 for i in range(len(features1))]
    else:
        squared_diffs = [(weights[i] ** 2) * ((features1[i] - features2[i]) ** 2) for i in range(len(weights))]
    return math.sqrt(sum(squared_diffs))


class DistContainer:
    def __init__(self):
        self.feature_function = lambda x: [0.0]
        self.distance_function = lambda x1, x2, w: 1

    def features(self, sequences: list[msa_parser.Sequence]) -> list[float]:
        return self.feature_function(sequences)

    def distance(self, reference: list[float], other: list[float], weights: list[float] = []) -> float:
        return self.distance_function(reference, other, weights)

    def features_from_part_dict(self, part_dict: dict) -> list[float]:
        return [0.0]


class SpartaDist(DistContainer):
    def __init__(self):
        super().__init__()
        self.feature_function = get_gap_features
        self.distance_function = sparta_dist


class ExtendedSpartaDist(DistContainer):
    def __init__(self):
        super().__init__()
        self.feature_function = get_extended_gap_features
        self.distance_function = sparta_extended_dist

    def features_from_part_dict(self, part_dict: dict) -> list[float]:
        raise NotImplementedError(f"{self.__class__.__name__}.features_from_part_dict()")


class BlindDist(DistContainer):
    def __init__(self):
        super().__init__()
        self.feature_function = get_simple_features
        self.distance_function = msa_blind_dist1

    def features_from_part_dict(self, part_dict: dict) -> list[float]:
        return [part_dict["NUM_ALIGNMENT_SITES"], part_dict["NUM_PATTERNS"], part_dict["GAPS"]]


def main():
    msa_path1 = sys.argv[1]

    sequences = msa_parser.parse_and_fix_msa_somehow(msa_path1, "DNA")
    features = _calc_gap_features(sequences)
    print(len(features.keys()))

    feature_vec = [features[key] for key in key_list]
    print(feature_vec)

    f1, f2 = get_features_from_sparta(msa_path1, msa_path1)

    for i in range(len(f1)):
        if feature_vec[i] != f1[i]:
            print(f"!= {key_list[i]}: {feature_vec[i]} - {f1[i]}")


if __name__ == "__main__":
    main()
