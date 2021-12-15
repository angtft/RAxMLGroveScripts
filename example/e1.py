#!/usr/bin/env python3

import collections
import concurrent.futures
import json
import os
import random
import subprocess
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

"""
(Execute from RAxMLGroveScripts/examples/ directory)

In this script we present a possible toy study which uses RAxMLGrove (RG) data to find out about possible relations 
between tree parameters and tree reconstructibility with RAxML-NG. We define tree reconstructibility (TR) as the ability 
to reconstruct a given true tree t after simulating MSAs based on t and its inferred parameters as determined during the 
phylogenetic tree inferece process of t. We assess the quality of TR using the Robinson-Foulds distance (RF-distance).

The experiments go as follows:
1. Pick non-outlier data sets from RG with the help of RAxMLGroveScripts (RGS). We also use the functionality of RGS
   to simulate MSAs for every picked data set.
2. We infer n trees t_i with RAxML-NG using the simulated MSA
3. We compute the RF-distance between every t_i and t and compute an average
4. For a column c (from a subset of tree parameter columns), we map for every data set d the average rf-distance to the 
   value of c from d... and plot the graphs
"""

RGS_PATH = "../"
RAXML_NG_PATH = os.path.abspath("../../../raxml-ng_v1.0.1_source/bin/raxml-ng")  # "PATH_TO_RAXML-NG"
OUT_PATH = "./out"
NUM_THREADS = 4
NUM_REPEATS = 10
NUM_START_TREES = 10


def concat_tree_files(root_path, file_names):
    """
    Concatenates multiple trees from different files into one file
    @param root_path: data set directory
    @param file_names: tree file names of the trees to be concatenated
    @return: 
    """
    out_file_name = "all_trees.newick"
    out_path = os.path.join(root_path, out_file_name)
    with open(out_path, 'w') as outfile:
        for fname in file_names:
            with open(fname) as infile:
                for line in infile:
                    outfile.write(line)
    return out_file_name


def run_raxml_inference(folder_path, prefix):
    """
    Run RAxML-NG with standard parameters
    @param folder_path: data set directory
    @param prefix: prefix for the output files (to keep them apart)
    @return:
    """
    msa_path = "assembled_sequences.fasta"
    command = f"{RAXML_NG_PATH} --msa {msa_path} --model GTR+G --prefix {prefix} --seed {random.randint(1,10000)} " \
              "--tree pars{" + str(NUM_START_TREES) + "},rand{" + str(NUM_START_TREES) + "} " \
              + f"--threads {NUM_THREADS}"
    try:
        subprocess.check_output(command.split(), cwd=folder_path)
    except Exception as e:
        print(e)


def run_raxml_rfdist(folder_path, all_tree_path):
    """
    Uses RAxML-NG to compute the RF-distance
    @param folder_path: data set directory
    @param all_tree_path: path to file with concatenated trees (for which we want to know the pairwise RF-distances)
    @return:
    """
    command = f"{RAXML_NG_PATH} --rfdist {all_tree_path}"

    try:
        out = subprocess.check_output(command.split(), cwd=folder_path)
    except Exception as e:
        print(e)


def get_average_rf_from_log(log_path):
    """
    Reads the RAxML-NG -rfdist log and returns the average RF-distance noted there
    @param log_path:
    @return:
    """
    ret = -1
    with open(log_path) as file:
        for line in file:
            line = line.rstrip()
            if line.startswith("Average relative RF distance in this tree set: "):
                ret = float(line.split()[-1])
    return ret


def get_tree_info(tree_dict_path, info_dict):
    """
    Takes the tree dict from path and fills the passed info_dict with values from 'interesting' columns. The info_dict
    shall contain for every interesting column a list of values found in the selected data sets
    @param tree_dict_path: path to tree dict file
    @param info_dict: dict to add data into
    @return:
    """
    def cast_float(value):
        try:
            ret = float(value)
            return ret
        except:
            return -1

    with open(tree_dict_path) as file:
        dct = json.load(file)

    keys = [
        "NUM_TAXA", "TREE_LENGTH", "TREE_DIAMETER", "TREE_HEIGHT", "BRANCH_LENGTH_MEAN", "BRANCH_LENGTH_VARIANCE",
        "OVERALL_NUM_ALIGNMENT_SITES", "OVERALL_NUM_PATTERNS", "INVARIANT_SITES", "OVERALL_NUM_PARTITIONS",
        "MISSING_DATA_RATE"
    ]
    for key in keys:
        info_dict[key].append(cast_float(dct[0][key]))

    info_dict["NUM_SITES_by_NUM_TAXA"].append(dct[0]["OVERALL_NUM_ALIGNMENT_SITES"] / dct[0]["NUM_TAXA"])
    info_dict["NUM_PATTERNS_by_NUM_TAXA"].append(dct[0]["OVERALL_NUM_PATTERNS"] / dct[0]["NUM_TAXA"])


def toy_study(num_points, num_repeats=NUM_REPEATS):
    """
    Conducts a toy study as described above
    @param num_points: number of data sets to include
    @param num_repeats: number of repeats per data set (since results might be noisy)
    @return:
    """

    # Use RGS to pick non-outlier trees and simulate MSAs
    generate_command_str = f"{RGS_PATH}/org_script.py generate --filter-outliers --num-msas {num_points} -o {OUT_PATH}"
    generate_command = generate_command_str.split()
    # generate_command.extend(["-q", "NUM_TAXA < 20"])  # uncomment for a faster test
    subprocess.check_output(generate_command)
    
    folders = os.listdir(OUT_PATH)
    info_dict = collections.defaultdict(lambda: [])
    y_list = []

    # Iterate over data set folders
    for folder in folders:
        folder_path = os.path.join(OUT_PATH, folder)
        # Run RAxML instances in paralles for different tree inferences for the same data set
        with concurrent.futures.ThreadPoolExecutor() as executor:
            func_calls = [(run_raxml_inference, folder_path, prefix) for prefix in range(0, num_repeats)]
            futures = [executor.submit(f, a, b) for f, a, b in func_calls]
            print(f"folder: {folder}")
            print([f.result() for f in futures])

            concat_tree_path = concat_tree_files(folder_path, [os.path.join(folder_path, f"{i}.raxml.bestTree") for i in range(0, num_repeats)])
            run_raxml_rfdist(folder_path, concat_tree_path)
            # Read the log file and get the average RF-distance
            rf = get_average_rf_from_log(os.path.join(folder_path, "all_trees.newick.raxml.log"))
            # Read the tree dict as generated by RGS for easier access of the tree parameters,
            # and add the parameters to info_dict
            get_tree_info(os.path.join(folder_path, "tree_dict.json"), info_dict)

            y_list.append(rf)

    print("rf dists:")
    print(y_list)

    # Draw plots
    for key in info_dict:
        print(f"key: {key}")
        x_list = info_dict[key]

        print(x_list)

        fig = plt.figure()
        plt.plot(x_list, y_list, "o")
        plt.xlabel(key)
        plt.ylabel("RF dist")
        plt.savefig(f"out_fig_{key}.png")
        plt.close(fig)


def main():
    toy_study(int(sys.argv[1]), int(sys.argv[2]))


if __name__ == "__main__":
    main()
