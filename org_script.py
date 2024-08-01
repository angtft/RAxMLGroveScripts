#!/usr/bin/env python3

import argparse
import collections
import copy
import itertools
import json
import math
import os
import random
import re
import shutil
import sqlite3
import statistics
import subprocess
import sys
import traceback
from io import StringIO
from typing import Type
from urllib.request import urlopen

from Bio import Phylo, SeqIO, Seq, SeqRecord
from skopt import Optimizer, space, gp_minimize
from skopt.utils import use_named_args

import tools.util.msa_parser as msa_parser
import tools.util.modified_SpartaABC as sparta_util

global_max_tree_file_len = 0

global_exception_counter = 0
global_num_of_too_big_trees = 0
global_num_of_checked_jobs = 0
global_tree_dict_list = []
global_part_list_dict = {}
global_db_object = None

global_node_counter = 0
global_tree_name_dict = {}
BASE_TREE_FORMAT = "newick"
BASE_TREE_NAME = "tree_{}.{}"
BASE_TREE_DICT_NAME = "tree_dict.json"
BASE_NODE_NAME = "taxon{}"

# Columns to be present in the SQLite database (as lists of tuples of name and data type).
META_COLUMNS = [
    ("RG_VERSION", "CHAR(100)"), ("COMMIT_HASH", "CHAR(100)"), ("URL", "CHAR(100)")
]

COLUMNS = [
    ("TREE_ID", "CHAR(255)"), ("NUM_TAXA", "INT"), ("TREE_LENGTH", "FLOAT"), ("TREE_DIAMETER", "FLOAT"),
    ("TREE_HEIGHT", "FLOAT"), ("BRANCH_LENGTH_MEAN", "FLOAT"), ("BRANCH_LENGTH_VARIANCE", "FLOAT"),
    ("IS_INDELIBLE_COMPATIBLE", "INT"), ("OVERALL_NUM_ALIGNMENT_SITES", "INT"), ("OVERALL_NUM_PATTERNS", "INT"),
    ("OVERALL_GAPS", "FLOAT"), ("RAXML_NG", "INT"),
    ("OVERALL_NUM_PARTITIONS", "INT"), ("MISSING_DATA_RATE", "FLOAT")
]

# TODO: fix STATIONARY_FREQ_STR (database had two different entries THAT and FREQ_STR)

PARTITION_COLUMNS = [
    ("MODEL", "CHAR(50)"), ("ALPHA", "FLOAT"), ("RATE_AC", "FLOAT"), ("RATE_AG", "FLOAT"), ("RATE_AT", "FLOAT"),
    ("RATE_CG", "FLOAT"), ("RATE_CT", "FLOAT"), ("RATE_GT", "FLOAT"), ("FREQ_A", "FLOAT"), ("FREQ_C", "FLOAT"),
    ("FREQ_G", "FLOAT"), ("FREQ_T", "FLOAT"), ("NUM_ALIGNMENT_SITES", "INT"), ("NUM_PATTERNS", "INT"),
    ("GAPS", "FLOAT"), ("INVARIANT_SITES", "FLOAT"), ("DATA_TYPE", "CHAR(50)"),
    ("RATE_STR", "CHAR(5000)"), ("PARTITION_NUM", "INT"),
    ("STATIONARY_FREQ_STR", "CHAR(100)"), ("PROPORTION_INVARIANT_SITES_STR", "CHAR(100)"),
    ("AMONG_SITE_RATE_HETEROGENEITY_STR", "CHAR(100)"), ("ASCERTAINMENT_BIAS_CORRECTION_STR", "CHAR(100)"),
    ("CUSTOM_CHAR_TO_STATE_MAPPING", "CHAR(100)"), ("PARENT_ID", "CHAR(255)")
]

SUBSTITUTION_MODELS = {
    "DNA": ['JC', 'K80', 'F81', 'HKY', 'TN93ef', 'TN93', 'K81', 'K81uf', 'TPM2', 'TPM2uf', 'TPM3', 'TPM3uf', 'TIM1',
            'TIM1uf', 'TIM2', 'TIM2uf', 'TIM3', 'TIM3uf', 'TVMef', 'TVM', 'SYM', 'GTR'],
    "AA": ['Blosum62', 'cpREV', 'Dayhoff', 'DCMut', 'DEN', 'FLU', 'HIVb', 'HIVw', 'JTT', 'JTT-DCMut', 'LG', 'mtART',
           'mtMAM', 'mtREV', 'mtZOA', 'PMB', 'rtREV', 'stmtREV', 'VT', 'WAG', 'LG4M', 'LG4X', 'PROTGTR'],
    "BIN": ["BIN"],
    "UNPHASED_DIPLOID_GENOTYPES": ['GTJC', 'GTHKY4', 'GTGTR4', 'GTGTR'],
    "MULTISTATE": ["MULTI"]
    # "USERDEFINED": ""     # TODO
}

BASE_GITHUB_LINK = "https://github.com/{}/{}/raw/{}/trees/{}/{}"  # c6ec6f73eedc42b20a08707060a2782d0b515599 hash of RG v0.2 commit
BASE_GITHUB_REPO_NAME = "RAxMLGrove"
BASE_GITHUB_REPO_OWNER = "angtft"

BASE_DB_FILE_NAME = "latest.db"
BASE_FILE_DIR = os.path.dirname(os.path.abspath(__file__))
BASE_OUT_DIR = os.path.join(BASE_FILE_DIR, "out")
BASE_STAT_OUT_FILE = os.path.join(BASE_FILE_DIR, "statistics.csv")
BASE_SEQ_FILE_FORMAT = "Phylip"  # TODO: maybe think some more about formats (since Phylip only allows taxon names up to 10 characters, and SeqGen doesn't output other formats?)
BASE_DAWG_SEQ_FILE_FORMAT = "Fasta"
BASE_SQL_FIND_COMMAND = "SELECT * FROM TREE t INNER JOIN PARTITION p ON t.TREE_ID = p.PARENT_ID"
BLANK_SYMBOL = "-"

DEFAULT_DB_FILE_PATH = os.path.join(BASE_FILE_DIR, BASE_DB_FILE_NAME)
DAWG_PATH = os.path.join(BASE_FILE_DIR, "tools", "dawg-1.2")
SEQGEN_PATH = os.path.join(BASE_FILE_DIR, "tools", "Seq-Gen-1.3.4")
ALISIM_PATH = os.path.join(BASE_FILE_DIR, "tools", "iqtree-2.2.0-beta-Linux")
GENESIS_PATH = os.path.join(BASE_FILE_DIR, "tools", "genesis-0.24.0")
RAXML_NG_PATH = os.path.join(BASE_FILE_DIR, "tools", "raxml-ng_v1.1.0_linux_x86_64", "raxml-ng")
SPARTAABC_PATH = os.path.join(BASE_FILE_DIR, "tools", "SpartaABC", "cpp_code", "SpartaABC")

SIMULATION_TIMEOUT = 20         # after this number of seconds, the MSA simulating process will be cancelled
SIMULATION_OPT_STOP_THRESH = 0.01
SPARTA_BURNIN_NUM = 1000       # default values 10k/100k (burn-in/sim)
SPARTA_SIM_NUM = 10000
EXTENDED_SPARTA_SIM_NUM = 100   # default 100
BONK_SIM_NUM = 100

BIG_NUMBER = sys.maxsize


def get_tukeys_fences(lst, k=1.5):
    """
    Computes the Tukey's fences for values in a given list and returns the low and high fences. Values below the lower
    fence and above the higher fence are considered to be outlier.
    @param lst: list with numerical values
    @param k: a nonnegative constant, usually 1.5
    @return: low fence, high fence
    """

    def is_float(value):
        try:
            c = float(value)
        except Exception as e:
            return False
        return True

    try:
        filtered_lst = sorted(float(v) for v in list(filter(lambda x: x != "None" and is_float(x), lst)))

        midpoint = int(round(len(filtered_lst) / 2.0))
        low_half_of_list = filtered_lst[:midpoint]
        high_half_of_list = filtered_lst[midpoint:]

        q1 = statistics.median(low_half_of_list)
        q3 = statistics.median(high_half_of_list)
        iqr = q3 - q1

        low_fence = q1 - k * iqr
        high_fence = q3 + k * iqr
    except Exception as e:
        raise e

    return low_fence, high_fence


def traverse_and_rename_nodes(clade):
    """
    Currently not used (?). Function is supposed to substitute taxon names in the tree, recursively.
    @param clade: root of the current subtree
    @return:
    """
    global global_node_counter
    global global_tree_name_dict
    if clade.name:
        if clade.name not in global_tree_name_dict:
            global_node_counter += 1
            global_tree_name_dict[clade.name] = BASE_NODE_NAME.format(global_node_counter)
        clade.name = global_tree_name_dict[clade.name]
    for c in clade.clades:
        traverse_and_rename_nodes(c)
    return 0


def copy_tree_file(src_path, dest_path):
    """
    Copies trees from source to destination paths. Used to anonymize trees as well at some point
    (currently RAxMLGrove trees are already anonymized).
    @param src_path: source tree file path
    @param dest_path: destination tree file path
    @return:
    """
    global global_node_counter
    global global_tree_name_dict
    global_tree_name_dict = {}
    try:
        tree_strings = []
        trees = []
        with open(src_path) as file:
            for line in file:
                tree_strings.append(line)

        for ts in tree_strings:
            global_node_counter = 0
            handle = StringIO(ts)
            tree = Phylo.read(handle, BASE_TREE_FORMAT)
            # traverse_and_rename_nodes(tree.root)  # TODO: remove?
            trees.append(tree)

        Phylo.write(trees, dest_path, BASE_TREE_FORMAT, plain=False, format_branch_length="%s")
    except Exception as e:
        print("Exception in copy_tree_file {}: {}".format(src_path, e))
        return False
    return True


def read_tree(path):
    """
    Reads a one-line file (the author specifically expects newick files here) and returns the string
    @param path: path to tree file
    @return: newick tree in string format
    """
    with open(path) as file:
        tree_string = file.read().rstrip()
    return tree_string


def create_dir_if_needed(path):
    """
    Creates a directory at path if that directory does not exist yet
    @param path: path to directory
    @return:
    """
    try:
        os.makedirs(path)
    except OSError as e:
        pass


def clear_directory(path):
    for file_name in os.listdir(path):
        os.remove(os.path.join(path, file_name))


def find_unused_tree_folder_name(root_dir_path, name):
    """
    This function finds the smallest number such that name + _number is an unused directory name in the specified
    root directory
    @param root_dir_path: root directory
    @param name: desired base name
    @return:
    """
    counter = 0
    while True:  # TODO: JPL ICS would not approve!
        if counter == 0:
            suggested_name = name
        else:
            suggested_name = f"{name}_{counter}"
        if not os.path.isdir(os.path.join(root_dir_path, suggested_name)):
            return suggested_name
        else:
            counter += 1


class Simulator:
    def __init__(self, path):
        self.path = path
        self.seed_line = ""

    def set_seed(self, seed):
        self.seed_line = seed

    def execute(self, tree_path, out_path, tree_params, new_seq_len=0, indel_rates=(), indel_distr=(), indel_distr_params=(), timeout=0):
        raise NotImplementedError()


class Dawg(Simulator):
    """
    We use Dawg here (https://github.com/reedacartwright/dawg)
    """

    def __init__(self, path):
        """
        @param path: path to Dawg root folder
        """
        super().__init__(path)

        self.template_path = os.path.join(self.path, "examples", "template.dawg")
        self.config_path = os.path.join(self.path, "examples", "template_modified.dawg")
        self.execute_path = os.path.join(self.path, "build", "src", "dawg")
        if os.path.isdir(self.path) and not os.path.isfile(self.execute_path):
            self.__compile()

    def set_seed(self, seed):
        """
        Sets the seed for the pseudo random number generator
        @param seed: seed
        @return:
        """
        self.seed_line = f"Seed = {seed}"

    def execute(self, tree_path, out_path, tree_params, new_seq_len=0, indel_rates=(), indel_distr=(), indel_distr_params=(), timeout=SIMULATION_TIMEOUT):
        """
        Executes Dawg to generate MSAs
        @param tree_path: tree file path
        @param out_path: output directory path
        @param tree_params: tree dict with tree information (such as model, substitution rates)
        @return:
        """
        out_dir = os.path.dirname(os.path.abspath(out_path))
        self.config_path = os.path.join(out_dir, "template_modified.dawg")

        self.__configure(tree_path, tree_params, new_seq_len=new_seq_len, indel_rates=indel_rates, out_path=out_path)
        try:
            call = [self.execute_path, self.config_path]
            if timeout:
                subprocess.run(call, cwd=out_dir, stdout=subprocess.DEVNULL, timeout=timeout)
            else:
                subprocess.run(call, cwd=out_dir, stdout=subprocess.DEVNULL)
        except subprocess.TimeoutExpired:
            print(f"Dawg.execute() timeout (after {SIMULATION_TIMEOUT}s)!")
            return 1
        except Exception as e:
            print(e)
            print(traceback.print_exc())
            return 1
        return 0

    def __compile(self):
        """
        Tries to compile Dawg
        @return:
        """
        print("Compiling Dawg...")

        build_path = os.path.join(self.path, "build")
        create_dir_if_needed(build_path)
        calls = [
            ["cmake", ".."],
            ["make"]
        ]

        try:
            for call in calls:
                subprocess.check_output(call, cwd=build_path)
        except Exception as e:
            print(e)
        print("Done!")

    def __configure(self, tree_path, tree_params, new_seq_len=0, indel_rates=(), out_path=""):
        """
        Modifies a Dawg configuration file template with the information found in a tree dict and writes the modified
        version to config_path
        @param tree_path: path to tree file
        @param tree_params: tree dict with tree information (such as model, substitution rates)
        @return:
        """
        tree_string = read_tree(tree_path)
        if not new_seq_len:
            seq_len = tree_params["NUM_ALIGNMENT_SITES"]
        else:
            seq_len = new_seq_len
        out_lines = []

        with open(self.template_path) as example_file:
            for line in example_file:
                if line.startswith("Tree"):
                    out_lines.append("Tree = " + tree_string + "\n")
                elif line.startswith("Length"):
                    out_lines.append(f"Length = {seq_len}\n")
                elif line.startswith("Params"):
                    out_lines.append("Params = {"
                                     + f'{tree_params["RATE_AC"]}, {tree_params["RATE_AG"]}, {tree_params["RATE_AT"]}, {tree_params["RATE_CG"]}, {tree_params["RATE_CT"]}, {tree_params["RATE_GT"]}'
                                     + "}\n")
                elif line.startswith("Freqs"):
                    out_lines.append("Freqs = {"
                                     + f'{tree_params["FREQ_A"]}, {tree_params["FREQ_C"]}, {tree_params["FREQ_G"]}, {tree_params["FREQ_T"]}'
                                     + "}\n")
                elif line.startswith("Model"):
                    out_lines.append('Model = "' + tree_params["MODEL"].split("+")[0] + '"\n')
                elif out_path and line.startswith("File"):
                    out_lines.append(f'File = "{os.path.basename(out_path)}"\n')
                else:
                    out_lines.append(line)
            out_lines.append(f"\n{self.seed_line}\n")
        if tree_params["ALPHA"] != "None":
            out_lines.append(f"Gamma = {tree_params['ALPHA']}\n")
        if indel_rates:
            out_lines.append("Lambda = {"
                             f"{indel_rates[0]}, {indel_rates[1]}"
                             "}\n")

        with open(self.config_path, "w+") as config_file:
            config_file.write("".join(out_lines))


class SeqGen(Simulator):
    """
    We also use Seq-Gen here (https://github.com/rambaut/Seq-Gen)
    """

    def __init__(self, path):  # TODO: fix out paths (for Dawg as well)
        """
        @param path: path to Seq-Gen base directory
        """
        super().__init__(path)
        self.path = path
        self.seed = -1
        self.execute_path = os.path.join(self.path, "source", "seq-gen")
        if not os.path.isfile(self.execute_path):
            self.__compile()

    def set_seed(self, seed):  # TODO: error handling (if seed is not an integer > 0)
        """
        Sets the seed for the pseudo random number generator
        @param seed: seed
        @return:
        """
        self.seed = seed

    def execute(self, tree_path, out_path, tree_params, new_seq_len=0, indel_rates=(), indel_distr=(), indel_distr_params=()):
        """
        @deprecated
        Executes Seq-Gen to generate MSA files
        @param tree_path: path to tree file
        @param out_path: output directory path
        @param tree_params: tree dict with the relevant model/rate/frequency information
        @param num_repeats: number of sequences to be generated (deprecated)
        @param num_of_sequence: sequence number (deprecated)
        @return:
        """
        out_path_abs = out_path if os.path.isabs(out_path) else os.path.join(BASE_FILE_DIR, out_path)
        try:
            call = [self.execute_path, tree_path, f"-m{tree_params['MODEL']}",
                    f"-f{tree_params['FREQ_A']},{tree_params['FREQ_C']},{tree_params['FREQ_G']},{tree_params['FREQ_T']}",
                    f"-r{tree_params['RATE_AC']},{tree_params['RATE_AG']},{tree_params['RATE_AT']},{tree_params['RATE_CG']},{tree_params['RATE_CT']},{tree_params['RATE_GT']}",
                    f"-l{tree_params['NUM_ALIGNMENT_SITES']}", "-or"]
            if self.seed > 0:
                call.append("-z")
                call.append(f"{self.seed}")
            if num_of_sequence == 0:
                for i in range(0, num_repeats):
                    out = subprocess.run(call, cwd=BASE_FILE_DIR, stderr=subprocess.DEVNULL, timeout=SIMULATION_TIMEOUT)
                    with open(os.path.join(out_path_abs, f"seq_{i}.{BASE_SEQ_FILE_FORMAT}"), "w+") as file:
                        file.write(out.decode("UTF-8"))
            else:
                out = subprocess.run(call, cwd=BASE_FILE_DIR, stderr=subprocess.DEVNULL, timeout=SIMULATION_TIMEOUT)
                with open(os.path.join(out_path_abs, f"seq_{num_of_sequence}.{BASE_SEQ_FILE_FORMAT}"), "w+") as file:
                    file.write(out.decode("UTF-8"))
        except Exception as e:
            print("Exception in SeqGen.execute(): {}".format(e))
            exit(0)

    def __compile(self):
        """
        Tries to compile Seq-Gen
        @return:
        """
        print("Compiling Seq-Gen...")

        build_path = os.path.join(self.path, "source")
        calls = [
            ["make"]
        ]

        try:
            for call in calls:
                subprocess.check_output(call, cwd=build_path)
        except Exception as e:
            print(e)
        print("Done!")


class AliSim(Simulator):
    """
    We use the AliSim here, which is now part of IQ-TREE2 (see https://github.com/iqtree/iqtree2/wiki/AliSim)
    """

    def __init__(self, path):
        """
        @param path: path to IQ-TREE2 root folder
        """
        super().__init__(path)
        self.execute_path = os.path.join(self.path, "bin", "iqtree2")

    def set_seed(self, seed):
        self.seed_line = f"{seed}"

    def replace_model_name(self, model):
        models = [
            "JC", "HKY", "GTR", "JTT", "LG", "BLOSUM62", "DAYHOFF"
        ]
        subs = {
            "BLOSUM62": "Blosum62",
            "PROTGTR": "GTR20",
            "PROTGTRGAMMA": "GTR20"
        }
        for key in models:
            if model.startswith("PROT"):
                if model.endswith(key):
                    if key in subs:
                        return subs[key]
                    return key
            else:
                if model.startswith(key):
                    if key in subs:
                        return subs[key]
                    return key
        return model

    def execute(self, tree_path, out_path, tree_params, new_seq_len=0, indel_rates=(), indel_distr=(), indel_distr_params=(), timeout=SIMULATION_TIMEOUT):
        out_dir = os.path.dirname(os.path.abspath(out_path))
        file_name = os.path.basename(out_path)

        try:
            if new_seq_len:
                seq_len = new_seq_len
            else:
                seq_len = tree_params["NUM_ALIGNMENT_SITES"]

            model = self.replace_model_name(tree_params['MODEL'])
            model_string = f"{model}"

            if tree_params["RATE_STR"] != "None" and "{" in tree_params["RATE_STR"]:
                rate_str = tree_params["RATE_STR"].split("{")[1]
                rate_str = rate_str.replace("{", "").replace("}", "")

                if tree_params["DATA_TYPE"] == "AA":
                    if model == "GTR":  # TODO: maybe handle these exceptions more gracefully...
                        model = "GTR20"
                        model_string = f"{model}"
                        # remove the base 1.0 in the end since AliSim does not understand it for some reason (?)
                        rate_str = "/".join(rate_str.split("/")[:-1])

                model_string += f"{{{rate_str}}}"
            if tree_params["STATIONARY_FREQ_STR"] != "None":
                freq_str = tree_params["STATIONARY_FREQ_STR"].split("{")[1]
                freq_str = freq_str.replace("{", "").replace("}", "")
                model_string += f"+F{{{freq_str}}}"

            if tree_params["ALPHA"] != "None":
                model_string += '+G{' + f'{tree_params["ALPHA"]}' + "}"
            if tree_params["INVARIANT_SITES"] != "None" and tree_params["RAXML_NG"] == 0:
                model_string += "+I{" + f"{float(tree_params['INVARIANT_SITES']) / 100}" + "}"
            elif tree_params["PROPORTION_INVARIANT_SITES_STR"] != "None":
                invar_str = tree_params['PROPORTION_INVARIANT_SITES_STR']
                invar_str = invar_str.split("{")[1].replace("}", "")
                model_string += f"+I{{{invar_str}}}"

            call = [
                self.execute_path,
                "--alisim",
                out_path,
                "-m", model_string,
                "-t", tree_path,
                "--length", f"{seq_len}",
                "-af", "fasta",
                "--no-export-sequence-wo-gaps"
            ]
            if indel_rates:
                call.extend([
                    "--indel", f'{indel_rates[0]},{indel_rates[1]}'
                ])
            if indel_distr and indel_distr_params:
                call.extend([
                    "--indel-size", f"{indel_distr[0]}" + "{" +
                                    f"{indel_distr_params[0][0]}/{indel_distr_params[0][1]}" + "}," +
                                    f"{indel_distr[1]}" + "{" +
                                    f"{indel_distr_params[1][0]}/{indel_distr_params[1][1]}" + "}"
                ])
            if self.seed_line:
                call.extend([
                    "--seed", self.seed_line
                ])

            subprocess.run(call, cwd=BASE_FILE_DIR, stdout=subprocess.DEVNULL, timeout=timeout)
        except subprocess.TimeoutExpired:
            print(f"AliSim.execute() timeout (after {timeout}s)!")
            return 1
        except Exception as e:
            print(e)
            print(traceback.print_exc())
            return 1

        if os.path.isfile(out_path + ".fa"):
            # deal with the annoying AliSim naming...
            os.rename(out_path + ".fa", out_path)
            return 0
        else:
            return 1


class BetterTreeDataBase(object):
    """
    The main SQLite database (db) object which is used to read and store tree information
    """

    def __init__(self, path, force_rewrite=False):
        """
        @param path: path to database file
        @param force_rewrite: if True, deletes old database file and creates a new one
        """
        if os.path.isabs(path):
            self.db_path = path
        else:
            self.db_path = os.path.join(BASE_FILE_DIR, path)

        setup_required = False
        if not os.path.isfile(self.db_path) or force_rewrite:
            setup_required = True

        self.conn = sqlite3.connect(self.db_path)
        self.conn.row_factory = sqlite3.Row
        self.cursor = self.conn.cursor()

        if setup_required:
            self.__prepare_empty_table()

    def __prepare_empty_table(self):
        self.cursor.execute("DROP TABLE IF EXISTS TREE")
        self.cursor.execute("DROP TABLE IF EXISTS PARTITION")
        self.cursor.execute("DROP TABLE IF EXISTS META_DATA")

        command = \
            f"""
                        CREATE TABLE META_DATA(
                            {", ".join([f"{entry} {type}" for entry, type in META_COLUMNS])}
                        );
                    """
        self.cursor.execute(command)

        command = \
            f"""
                CREATE TABLE TREE(
                    {", ".join([f"{entry} {type}" for entry, type in COLUMNS])}
                );
            """
        self.cursor.execute(command)

        command = \
            f"""
                CREATE TABLE PARTITION(
                    {", ".join([f"{entry} {type}" for entry, type in PARTITION_COLUMNS])},
                    FOREIGN KEY(PARENT_ID) REFERENCES TREE(TREE_ID)
                );
            """
        self.cursor.execute(command)
        self.conn.commit()

    def close(self):
        self.conn.close()

    def fill_database(self, tree_dict_list, partition_list_dict, meta_info_dict):
        """
        Fills the db file with information from the tree and partition dicts
        @param tree_dict_list: list of tree dicts
        @param partition_list_dict: list of partition dicts
        @param meta_info_dict: dict with meta information (such as the relevant RG commit hash)
        @return:
        """
        for entry, _ in META_COLUMNS:
            if entry not in meta_info_dict:
                meta_info_dict[entry] = None

        try:
            command = \
                f"""
                    INSERT INTO META_DATA({", ".join([f"{entry}" for entry, _ in META_COLUMNS])})
                    VALUES ({", ".join([f"'{meta_info_dict[entry]}'" for entry, _ in META_COLUMNS])});
                """
            self.cursor.execute(command)
        except Exception as e:
            print(f"Exception in fill_database during writing of meta information: {e}")

        for key in partition_list_dict:
            for part in partition_list_dict[key]:
                for entry, _ in PARTITION_COLUMNS:
                    if entry not in part:
                        part[entry] = None

        for dct in tree_dict_list:
            for entry, _ in COLUMNS:
                if entry not in dct:
                    dct[entry] = None

            try:
                tree_id = dct["TREE_ID"]
                for part_dct in partition_list_dict[tree_id]:
                    part_command = \
                        f"""
                        INSERT INTO PARTITION({
                        ", ".join([f"{entry}" for entry, _ in PARTITION_COLUMNS])
                        })
                        VALUES ({
                        ", ".join([f"'{part_dct[entry]}'" for entry, _ in PARTITION_COLUMNS])
                        });
                        """
                    self.cursor.execute(part_command)

                command = \
                    f"""
                        INSERT INTO TREE({", ".join([f"{entry}" for entry, _ in COLUMNS])}) 
                        VALUES ({", ".join([f"'{dct[entry]}'" for entry, _ in COLUMNS])});
                    """
                self.cursor.execute(command)
            except Exception as e:
                print(f"Exception in fill_database: {e}")
                continue
        self.conn.commit()

    def database_entry_exists(self, id):
        """
        Checks if a given tree id is already present in the db
        @param id: tree id
        @return: True if tree is present in the db
                 False otherwise
        """
        command = \
            f"""
                SELECT {id} FROM TREE;
            """
        self.cursor.execute(command)
        result = self.cursor.fetchall()
        if len(result):
            return True
        else:
            return False

    def execute_command(self, command):
        """
        Executes a given SQL command on the db
        @param command: command to execute
        @return:
        """
        # Allow only reading access here!
        self.cursor.execute(command)

    def find(self, command):
        """
        Expects a command to perform a query on the db and to return a list of dicts with found entries
        @param command: "SELECT [...]" command to execute on the database
        @return: list of results (results being dicts of column entries)
        """
        self.cursor.execute(command)
        result = self.cursor.fetchall()
        return [dict(row) for row in result]

    def get_meta_info(self):
        """
        Reads the meta info of the current db from the META_DATA table. If something goes wrong,
        returns a default dict with the commit hash of RAxMLGrove v0.2.
        @return: meta info dict
        """
        meta_info_dict = {"COMMIT_HASH": "c6ec6f73eedc42b20a08707060a2782d0b515599"}  # commit hash of RG v0.2
        for entry, _ in META_COLUMNS:
            if entry not in meta_info_dict:
                meta_info_dict[entry] = None

        try:
            query = f"SELECT * FROM META_DATA"
            self.cursor.execute(query)
            results = self.cursor.fetchall()
            meta_dict = dict(results[0])
            return meta_dict
        except Exception as e:
            print(f"Exception in get_meta_info: {e}")
            return meta_info_dict
        return meta_info_dict


class RaxmlNGLogReader(object):
    """
    Used to read the log and model files produced by RAxML-NG
    and extract the information needed for the columns of the database ("create")
    """

    def __init__(self, path, model_path):
        """
        @param path: path to log file
        @param model_path: path to model output file
        """
        self.path = path
        self.partitions_dict = {}
        self.partition_names = []
        try:
            self.__read()
            self.__read_model(model_path)
            self.__fill_general_info()
        except Exception as e:
            print(f"Exception in raxml-ng: {self.path}\n{e}")
            global global_exception_counter
            global_exception_counter += 1
            raise e

    def __read(self):
        """
        Extracts information from log file, puts it into partition dict
        @return:
        """
        with open(self.path) as file:
            part_info_dict = {}
            for line in file:
                line = line.strip()
                if line.startswith("Partition"):
                    if "name" in part_info_dict:
                        self.partitions_dict[part_info_dict["name"]].update(part_info_dict)
                        self.partition_names.append(part_info_dict["name"])
                        part_info_dict = {}

                    part_info_dict["name"] = line.split()[-1]
                    if part_info_dict["name"] not in self.partitions_dict:
                        self.partitions_dict[part_info_dict["name"]] = {}

                if line.startswith("Alignment sites / patterns:"):
                    part_info_dict["NUM_ALIGNMENT_SITES"] = int(line.split()[-3])
                    part_info_dict["NUM_PATTERNS"] = int(line.split()[-1])
                if line.startswith("Gaps:"):
                    part_info_dict["GAPS"] = float(line.split()[-2])
                if line.startswith("Invariant sites:"):
                    part_info_dict["INVARIANT_SITES"] = float(line.split()[-2])
                if "Loaded alignment with" in line:     # this line is currently not be available in anonymized data!
                    tres = re.findall("taxa and (.*?) sites", line)
                    sl = int(tres[0])
                    part_info_dict["NUM_ALIGNMENT_SITES"] = sl

            self.partition_names.append(part_info_dict["name"])
            self.partitions_dict[part_info_dict["name"]].update(part_info_dict)

    def __get_modifier(self, modifier_str):
        modifier = modifier_str.split("{")[0]
        return modifier

    def __get_values_from_modifiers(self, modifier_str):
        temp_res = re.search("{(.*)}", modifier_str)
        if temp_res:
            inner_values = temp_res.group(1)
            return inner_values.split("/")
        else:
            return []

    def __read_model_fix(self, path):
        # TODO: implement (for cases where log file contains duplicated partition descriptions, for some reason...)
        modifier_dict = {
            "STATIONARY_FREQ_STR": ["F", "FC", "FO", "FE", "FU"],  # stationary frequencies
            "PROPORTION_INVARIANT_SITES_STR": ["I", "IO", "IC", "IU"],  # proportion of invariant sites
            "AMONG_SITE_RATE_HETEROGENEITY_STR": ["G", "G4m", "R"],  # among-site rate heterogeneity model
            "ASCERTAINMENT_BIAS_CORRECTION_STR": ["ASC_LEWIS", "ASC_FELS", "ASC_STAM"],  # ascertainment bias correction
            "CUSTOM_CHAR_TO_STATE_MAPPING": ["M", "Mi"]
        }
        with open(path) as file:
            pass

    def __read_model(self, path):
        """
        Reads the model output file and fills the partition dict with that info as well
        @param path: path to model file
        @return:
        """
        modifier_dict = {
            "STATIONARY_FREQ_STR": ["F", "FC", "FO", "FE", "FU"],  # stationary frequencies
            "PROPORTION_INVARIANT_SITES_STR": ["I", "IO", "IC", "IU"],  # proportion of invariant sites
            "AMONG_SITE_RATE_HETEROGENEITY_STR": ["G", "G4m", "R"],  # among-site rate heterogeneity model
            "ASCERTAINMENT_BIAS_CORRECTION_STR": ["ASC_LEWIS", "ASC_FELS", "ASC_STAM"],  # ascertainment bias correction
            "CUSTOM_CHAR_TO_STATE_MAPPING": ["M", "Mi"]
        }

        with open(path) as file:
            lines = file.readlines()
            if len(lines) - 1 > len(self.partition_names):
                raise ValueError(f"lines {len(lines)} > part_names {len(self.partition_names)}")

            for i in range(len(lines)):
                try:
                    line = lines[i].strip()
                    if line:
                        part_key = self.partition_names[i]
                    else:
                        continue
                except Exception as e:
                    print(f"Exception in raxml-ng __read_model: {self.path}\n{e}")
                    print(f"lines {len(lines)}  part_names {len(self.partition_names)}")
                    global global_exception_counter
                    global_exception_counter += 1
                    continue

                modifiers_info = line.rstrip().split(",")[0].split("+")

                model = self.__get_modifier(modifiers_info[0])
                self.partitions_dict[part_key]["RATE_STR"] = modifiers_info[0]
                self.partitions_dict[part_key]["MODEL"] = model
                rates = self.__get_values_from_modifiers(self.partitions_dict[part_key]["RATE_STR"])
                if model == "GTR":  # TODO: maybe remove these fields completely
                    self.partitions_dict[part_key]["RATE_AC"] = rates[0]
                    self.partitions_dict[part_key]["RATE_AG"] = rates[1]
                    self.partitions_dict[part_key]["RATE_AT"] = rates[2]
                    self.partitions_dict[part_key]["RATE_CG"] = rates[3]
                    self.partitions_dict[part_key]["RATE_CT"] = rates[4]
                    self.partitions_dict[part_key]["RATE_GT"] = rates[5]

                for key in SUBSTITUTION_MODELS:
                    if model in SUBSTITUTION_MODELS[key]:
                        self.partitions_dict[part_key]["DATA_TYPE"] = key

                for mi in modifiers_info:
                    modifiers = re.findall(r"(.*?)\{[\d|\.|\/]*\}\+*", mi)
                    if modifiers:
                        modifier = modifiers[0]
                    else:
                        modifier = None

                    for category in modifier_dict:
                        if modifier in modifier_dict[category]:
                            values = re.findall(r"\{(.*?)\}", mi)
                            self.partitions_dict[part_key][category] = mi
                            if category == "AMONG_SITE_RATE_HETEROGENEITY_STR":
                                alpha = values[0]
                                self.partitions_dict[part_key]["ALPHA"] = float(alpha)
                                if "/" in alpha:
                                    print(self.path)
                                    raise ValueError(alpha)

    def __fill_general_info(self):
        """
        use after __read() and __read_model() were executed
        """
        for part_key in self.partitions_dict:
            if "ASCERTAINMENT_BIAS_CORRECTION_STR" in self.partitions_dict[part_key] and \
                    "G" in self.partitions_dict[part_key]["ASCERTAINMENT_BIAS_CORRECTION_STR"]:
                self.partitions_dict[part_key]["ALPHA"] = \
                    self.__get_values_from_modifiers(
                        self.partitions_dict[part_key]["ASCERTAINMENT_BIAS_CORRECTION_STR"])[0]

            if "STATIONARY_FREQ_STR" in self.partitions_dict[part_key] and \
                    "DATA_TYPE" in self.partitions_dict[part_key] and \
                    self.partitions_dict[part_key]["DATA_TYPE"] == "DNA":
                values = self.__get_values_from_modifiers(self.partitions_dict[part_key]["STATIONARY_FREQ_STR"])
                if values:
                    self.partitions_dict[part_key]["FREQ_A"] = values[0]
                    self.partitions_dict[part_key]["FREQ_C"] = values[1]
                    self.partitions_dict[part_key]["FREQ_G"] = values[2]
                    self.partitions_dict[part_key]["FREQ_T"] = values[3]

    def get_partition_info(self):
        return self.partitions_dict


class OldRaxmlReader(object):
    """
    Used to read the log of a RAxML file and extract the information needed for the columns of the database ("create")
    """

    def __init__(self, path):
        """
        @param path: path to the RAxML log file
        """
        self.path = path
        self.partitions_dict = {}
        self.__read()
        self.__fill_model_info()

    def __parse_known_args(self, command_line):
        parser = argparse.ArgumentParser()
        parser.add_argument("-m")
        parser.add_argument("-A")
        parser.add_argument("-b")
        parser.add_argument("-B")
        parser.add_argument("-c")
        parser.add_argument("-C", action='store_true')
        parser.add_argument("-d", action='store_true')
        parser.add_argument("-D", action='store_true')
        parser.add_argument("-e")
        parser.add_argument("-f")
        parser.add_argument("-F", action='store_true')
        parser.add_argument("-G")
        parser.add_argument("-H", action='store_true')
        parser.add_argument("-i")
        parser.add_argument("-I")
        parser.add_argument("-j", action='store_true')
        parser.add_argument("-J")
        parser.add_argument("-k", action='store_true')
        parser.add_argument("-K", action='store_true')
        parser.add_argument("-L")
        parser.add_argument("-M", action='store_true')
        parser.add_argument("-O", action='store_true')
        parser.add_argument("-p")
        parser.add_argument("-P")
        parser.add_argument("-T")
        parser.add_argument("-u", action='store_true')
        parser.add_argument("-U", action='store_true')
        parser.add_argument("-v", action='store_true')
        parser.add_argument("-V", action='store_true')
        parser.add_argument("-W")
        parser.add_argument("-x")
        parser.add_argument("-X", action='store_true')
        parser.add_argument("-y", action='store_true')
        parser.add_argument("-#")
        parser.add_argument("-N")

        args = None
        try:
            args, _ = parser.parse_known_args(command_line.split()[1:])
        except Exception as e:
            print(f"Exception in __parse_known_args: {e}")
            print(traceback.print_exc())
        return args

    def __read(self):
        """
        Very big function to read all the different parameters listed in a log file and write them to the object's
        partitions dict
        @return:
        """

        temp_part_dict = {
            "NUM_ALIGNMENT_SITES": [],
            "NUM_PATTERNS": [],
            "DATA_TYPE": [],
            "MODEL": [],
            "BASE_FREQUENCIES": [],
            "ALPHA": [],
            "RATES": [],
            "TREE_LENGTH": [],
            "GAPS": []
        }

        overall_num_sites = 0
        proportion_of_gaps = 0

        def get_value(line):
            return line.split(":")[-1].strip()

        with open(self.path) as file:
            current_rates = []
            current_freqs = []
            alphas = []  # TODO: currently we only take one assignment of alphas and rates, even if
            rates = []  # multiple searches were performed
            command_line = ""

            for line in file:
                line = line.rstrip()
                if line.startswith("Alignment Patterns:"):
                    value = int(get_value(line))
                    temp_part_dict["NUM_PATTERNS"].append(value)
                elif line.startswith("Proportion of gaps and completely undetermined characters in this alignment:"):
                    value_pct = get_value(line)
                    value = float(value_pct.split("%")[0])
                    proportion_of_gaps = value
                elif line.startswith("Alignment sites:"):
                    value = line.split(":")
                    if len(value) > 1:
                        rside = " ".join(value[1].strip().split())
                        rside = rside.split()
                        if len(rside) > 1:
                            overall_num_sites = int(rside[1])
                        elif len(rside) == 1:
                            overall_num_sites = int(rside[0])
                elif line.startswith("sites partition_"):
                    value = line.split("=")
                    if not "None" in value[1]:
                        try:
                            intervals = value[1].strip().split(",")
                            part_size = 0
                            for interval in intervals:
                                temp_split = interval.split("\\")
                                summands = temp_split[0].replace(" ", "").split("-")
                                divisor = int(temp_split[1]) if len(temp_split) > 1 else 1

                                s1 = int(summands[0])
                                if len(summands) == 2:
                                    s2 = int(summands[1])
                                elif len(summands) == 1:
                                    s2 = s1
                                else:
                                    raise ValueError(f"len(summands) = {len(summands)}")
                                part_size += int((s2 - s1 + 1) / divisor)
                                if part_size <= 0:
                                    raise ValueError(f"Partition size <= 0: {part_size}")
                            temp_part_dict["NUM_ALIGNMENT_SITES"].append(part_size)
                        except Exception as e:
                            print(f"Exception in old_raxml __read partition sites: {self.path}\n{e}")
                            print(traceback.print_exc())
                            global global_exception_counter
                            global_exception_counter += 1
                            temp_part_dict["NUM_ALIGNMENT_SITES"].append(None)
                elif line.startswith("DataType:"):
                    value = get_value(line)
                    temp_part_dict["DATA_TYPE"].append(value)
                elif line.startswith("Substitution Matrix:"):
                    value = get_value(line)
                    temp_part_dict["MODEL"].append(value)
                elif line.startswith("Base frequencies:"):
                    value = get_value(line)
                    try:
                        value_list = [float(x) for x in value.split()]
                        temp_part_dict["BASE_FREQUENCIES"].append(value_list)
                    except Exception:
                        pass
                elif line.startswith("Tree-Length:"):
                    value = float(get_value(line))
                    temp_part_dict["TREE_LENGTH"].append(value)
                elif line.startswith("alpha:"):
                    value = float(get_value(line))
                    temp_part_dict["ALPHA"].append(value)

                if line.startswith("raxml"):
                    command_line = line

                if line.startswith("rate "):
                    value = float(get_value(line))
                    current_rates.append(value)
                else:
                    if current_rates:
                        temp_part_dict["RATES"].append(copy.deepcopy(current_rates))
                        current_rates = []

                if line.startswith("freq "):
                    value = float(get_value(line))
                    current_freqs.append(value)
                else:
                    if current_freqs:
                        temp_part_dict["BASE_FREQUENCIES"].append(copy.deepcopy(current_freqs))
                        current_freqs = []

                if "[" in line:
                    alphas = re.findall(r"alpha\[(.*?)\]: (.*?) ", f"{line} ")
                    rates = re.findall(r"rates\[(.*?)\] ac ag at cg ct gt: (.*?) (.*?) (.*?) (.*?) (.*?) (.*?) ",
                                       f"{line} ")

            if current_rates:
                temp_part_dict["RATES"].append(copy.deepcopy(current_rates))
            if current_freqs:
                temp_part_dict["BASE_FREQUENCIES"].append(copy.deepcopy(current_freqs))

            num_partitions = len(temp_part_dict["NUM_PATTERNS"])  # TODO: this should hopefully be representative

            if not len(temp_part_dict["MODEL"]) == len(temp_part_dict["NUM_PATTERNS"]):
                print(f"{len(temp_part_dict['MODEL'])} {len(temp_part_dict['NUM_PATTERNS'])}")
                raise ValueError(f"len(models) len(patterns) mismatch in __read")

            alpha_idx = 0
            rate_idx = 0
            for i in range(num_partitions):
                if alpha_idx < len(alphas):
                    ta = alphas[alpha_idx]
                    alpha_num = int(ta[0])
                    if alpha_num == i:
                        temp_part_dict["ALPHA"].append(float(ta[1]))
                        alpha_idx += 1
                    else:
                        temp_part_dict["ALPHA"].append(None)

                if rate_idx < len(rates):
                    tr = rates[rate_idx]
                    rate_num = int(tr[0])
                    if rate_num == i:
                        temp_part_dict["RATES"].append(tr[1:])
                        rate_idx += 1
                    else:
                        temp_part_dict["RATES"].append(None)

                temp_part_dict["GAPS"].append(proportion_of_gaps)

            if num_partitions == 1:
                temp_part_dict["NUM_ALIGNMENT_SITES"] = [overall_num_sites]
                parsed_args = self.__parse_known_args(command_line)
                temp_part_dict["MODEL"] = [parsed_args.m]

            for i in range(num_partitions):
                new_part = {}
                for key in temp_part_dict:
                    if temp_part_dict[key] and len(temp_part_dict[key]) == num_partitions:
                        new_part[key] = temp_part_dict[key][i]
                self.partitions_dict[str(i)] = copy.deepcopy(new_part)

    def __fill_model_info(self):
        """
        For the DNA data sets inferred under the GTR model (which are overall most of the data sets), fills the
        partition dict with subsitution rates and frequencies found in the log
        @return:
        """
        for part_key in self.partitions_dict:
            part = self.partitions_dict[part_key]
            if "MODEL" in part and "RATES" in part:
                part["RATE_STR"] = f"{part['MODEL']}{{{'/'.join([str(x) for x in part['RATES']])}}}"
                if part["MODEL"] == "GTR" and "DATA_TYPE" in part and part["DATA_TYPE"] == "DNA":
                    part["RATE_AC"] = part["RATES"][0]
                    part["RATE_AG"] = part["RATES"][1]
                    part["RATE_AT"] = part["RATES"][2]
                    part["RATE_CG"] = part["RATES"][3]
                    part["RATE_CT"] = part["RATES"][4]
                    part["RATE_GT"] = part["RATES"][5]
                    if "BASE_FREQUENCIES" in part:
                        part["FREQ_A"] = part["BASE_FREQUENCIES"][0]
                        part["FREQ_C"] = part["BASE_FREQUENCIES"][1]
                        part["FREQ_G"] = part["BASE_FREQUENCIES"][2]
                        part["FREQ_T"] = part["BASE_FREQUENCIES"][3]

            if "BASE_FREQUENCIES" in part:
                part["STATIONARY_FREQ_STR"] = f"{{{'/'.join([str(x) for x in part['BASE_FREQUENCIES']])}}}"

    def get_partition_info(self):
        """
        Returns the partitions dict, which maps partition numbers to dicts with information about that partition
        @return: partition dict
        """
        return self.partitions_dict


class GenesisTreeDiameter(object):  # TODO: Add genesis to ./tools/
    """
    We use Genesis here (https://github.com/lczech/genesis) to make our life easier when collecting tree parameters.
    """

    def __init__(self, path):
        """
        @param path: path to Genesis base directory
        """
        self.path = path
        self.executable_path = os.path.join(self.path, "bin", "apps", "tree_diameter")
        if not os.path.isfile(self.executable_path):
            self.__compile()

    def __compile(self):
        """
        Compiles (or at least tries to) Genesis
        @return:
        """
        print("Compiling genesis...")

        build_path = self.path
        calls = [
            ["make"]
        ]

        try:
            for call in calls:
                subprocess.check_output(call, cwd=build_path)
        except Exception as e:
            print(e)
            print(traceback.print_exc())
        print("Done!")

    def get_len_and_diam_and_height(self, tree_path1):
        """
        Uses Genesis to compute some tree parameters
        @param tree_path1: path to tree file
        @return: length, diameter and height of the tree, if successful
                 -1, -1, -1 otherwise
        """
        call = [self.executable_path, tree_path1]

        try:
            out = subprocess.check_output(call).decode()
            values = out.split()
            return float(values[-3]), float(values[-2]), float(values[-1])
        except Exception as e:
            print(e)
            print(traceback.print_exc())
            return -1, -1, -1


def init_args(arguments):
    """
    Parses command line arguments
    @param arguments: command line arguments
    @return: object with set arguments
    """

    parser = argparse.ArgumentParser()
    parser.add_argument("operation",
                        choices=["create", "add", "execute", "find", "generate", "stats", "justgimmeatree",
                                 "download"],
                        # TODO: rework 'find' command?
                        help="'create' iterates over a RAxML out files archive parametrized with '-a' and "
                             "writes a database (db) file (the name can be parametrized with '-n'). Default db name "
                             f"{BASE_DB_FILE_NAME}.\n"
                             "'add' adds trees from archive path to the db.\n"
                             "'execute' executes a command parametrized with '-c' on the db.\n"
                             "'find' tries to find trees in the db satisfying the conditions set with '-q'.\n"
                             "'download' downloads all the trees returned by the '-q' query to the directory "
                                 "specified with '-o' without asking questions.\n"
                             "'generate' generates simulated MSAs using randomly drawn trees from the db.\n"
                             "'stats' prints some statistical information about the database entries.")
    parser.add_argument("-n", "--db-name", help=f"Sets the name of the database file. If not set, this tool will "
                                                f"try to create/access the '{BASE_DB_FILE_NAME}' file.")
    parser.add_argument("-a", "--archive-path", help="Sets the raxml out files archive path. (create)")
    parser.add_argument("--rg-commit-hash", help="Sets a specific commit hash of RAxMLGrove in the db metadata. "
                                                 "(create)")
    parser.add_argument("--download-link", default=BASE_GITHUB_LINK,
                        help="Sets the source link for the data sets. This argument expects a format string of a "
                             "GitHub repository link with replacement fields for (1) commit hash, (2) data set id, "
                             f"and (3) file name. Default: '{BASE_GITHUB_LINK}'. (create)")
    parser.add_argument("-c", "--command", help="The command to execute on the database. (execute)")
    parser.add_argument("-q", "--query", help="Part of the statement after the 'WHERE' clause "
                                              "to find trees in the db. CAUTION: We do not sanitize, "
                                              "don't break your own database...")
    parser.add_argument("--list", action='store_true', help="Lists all found trees, "
                                                            "without asking for download. (find)")
    parser.add_argument("--num-msas", default=1, help="CURRENTLY NOT WORKING! "  # TODO: make this work again!
                                                      "Number of MSAs to be generated. To generate "
                                                      "a MSA we randomly draw a tree from the db (can be "
                                                      "used with -q) and run a sequence generator with that tree. "
                                                      "(generate)")
    parser.add_argument("--local", default=True, help="FOR TESTING PURPOSES (don't use it).")
    parser.add_argument("--force-rewrite", action='store_true', help="Forces to rewrite the default db or the db "
                                                                     "specified using '-n'. (create)")
    parser.add_argument("--seq-len", default=0, help="Sets the sequence length of the generated MSA (starting "
                                                     "sequence length if indels are used). (generate)")
    parser.add_argument("--filter-outliers", action='store_true', help="Filters trees with uncommon characteristics "
                                                                       "(using Tukey's fences). (generate)")
    parser.add_argument("--insert-matrix-gaps", action="store_true", help="Uses the presence/absence matrices to "
                                                                          "insert gaps into the simulated sequences. "
                                                                          "(generate)")
    parser.add_argument("--use-all-trees", action="store_true", help="Forces the usage of all found trees instead of "
                                                                     "pulling trees randomly from the set of found "
                                                                     "trees. (generate)")
    parser.add_argument("--indel",
                        help=argparse.SUPPRESS)
                        #help="CURRENTLY NOT WORKING! Comma separated insertion/deletion rates for the MSA generation. "  # TODO: make this work
                        #                "(generate)")
    parser.add_argument("--use-bonk", action="store_true",    # TODO: rework this
                        help="Use an optimizer to experiment with different gap indel rates during the MSA simulation "
                             "to get the sequence length, number of patterns and gap proportions of the simulated "
                             "MSA as close as possible to the values in the according RG database entry. (generate)")
    parser.add_argument("-g", "--generator", choices=["dawg", "seq-gen", "alisim"], default="alisim",
                        help="Selects the sequence generator used. (generate)")
    parser.add_argument("-o", "--out-dir", default=BASE_OUT_DIR,
                        help=f"Output directory (default: {BASE_OUT_DIR}). (find, generate)")
    parser.add_argument("--seed", help="Sets the seed for the random number generator of this script, "
                                       "as well as the sequence generators. (generate)")
    parser.add_argument("--weights",
                        help=argparse.SUPPRESS)
                        #help="DON'T USE IT! format: mw1,mw2,mw3,range_factor")         # TODO: this should be removed at this point
    parser.add_argument("--use-sparta", action="store_true",
                        help=argparse.SUPPRESS)
                        #help="DON'T (maybe) USE IT! Uses the slightly modified SpartaABC tool "     # TODO: this
                        #     "to infer indel rates and sizes from either the available original "
                        #     "msa.fasta or the according feature vector in the db. (generate)")
    parser.add_argument("--use-modified-sparta", action="store_true",
                        help=argparse.SUPPRESS)
                        #help="DON'T USE IT! Same idea as with '--use-sparta'. Here, however, we give "
                        #     "a modified distance function to the optimizer. (generate)")
    parser.add_argument("--avoid-empty-sequences", action="store_true",
                        help="Makes sure the simulated MSA contains no empty sequences by resimulating the "
                             "MSA or not inserting indels into one of the partitions. (generate)")
    parser.add_argument("--use-local-db",
                        help=argparse.SUPPRESS)
                        #help="DON'T USE THIS. Sets the specified local path as the source "
                        #                       "of data sets for MSA simulations. (generate)")
    parser.add_argument("--no-simulation", action="store_true", help=argparse.SUPPRESS)     # TODO: workaround for our eval pipeline...

    args = parser.parse_args(arguments)
    if args.operation == "create" and not args.archive_path:
        parser.error("'create' requires '-a'")
    if args.operation == "add" and not args.archive_path:
        parser.error("'add' requires '-a'")
    if args.operation == "execute" and not args.command:
        parser.error("'execute' requires '-c'")
    if args.operation == "find" and not args.query:
        parser.error("'find' requires '-q'")
    if args.operation in ["find", "generate"] and args.out_dir == BASE_OUT_DIR:
        print(f"No output directory specified. Using default '{BASE_OUT_DIR}'.")
        create_dir_if_needed(BASE_OUT_DIR)
    if args.list and not args.operation == "find":
        parser.error("'--list' may only be used with 'find'")
    if args.seed:
        try:
            s = int(args.seed)
            if s < 1:
                parser.error("The seed must be a non-zero, positive integer.")
        except Exception as e:
            parser.error("The seed must be a non-zero, positive integer.")
        random.seed(int(args.seed))
    elif args.indel and args.use_bonk:
        parser.error("Mutually exclusive arguments used: --indel and --use-bonk")
    if args.use_local_db:
        args.use_local_db = os.path.abspath(args.use_local_db)
    if args.use_sparta and args.use_modified_sparta:
        parser.error("Mutually exclusive arguments used: --use-sparta and --use-modified-sparta")
    if args.generator != "alisim" and (args.use_bonk or args.use_sparta or args.use_modified_sparta or args.indel):
        parser.error("Indels are currently only usable with AliSim!")
    if args.seq_len:
        try:
            args.seq_len = int(args.seq_len)
            if args.seq_len <= 0:
                parser.error("Sequence length may not be <= 0!")
        except:
            parser.error(f"Error converting sequence length {args.seq_len} to int!")

    return args


def yes_no_prompt(message):
    """
    Asks the user for a 'yes' or 'no'
    @param message: displayed question
    @return: True if user picks 'yes', False if user picks 'no' otherwise
    """
    print(message + " y/n")
    while True:
        user_in = input('>>> ')
        if user_in in ('y', 'yes'):
            return True
        elif user_in in ('n', 'no'):
            return False
        else:
            print('Please answer with "y" (yes) or "n" (no)!')


def save_tree_dict(path, result):
    """
    Saves the tree dict to file
    @param path: output file path
    @param result: tree dict
    @return:
    """
    try:
        json_obj = json.dumps(result, indent=4)
        out_path = os.path.join(path, BASE_TREE_DICT_NAME)
        with open(out_path, "w+") as file:
            file.write(json_obj)
    except Exception as e:
        print(e)
        print(f"Error: could not export tree dict to file: {path}")
        print(traceback.print_exc())


def write_partitions_file(dir_path, partition_results, file_name=""):
    sorted_parts = copy.deepcopy(partition_results)
    print(sorted_parts)
    try:
        sorted_parts.sort(key=lambda x: int(x["PARTITION_NUM"]))
    except ValueError:
        if sorted_parts[0]["OVERALL_NUM_PARTITIONS"] != 1:
            raise ValueError(f"{sorted_parts['TREE_ID']} contains partition(s) with invalid PARTITION_NUM")

    if not file_name:
        file_path = os.path.join(dir_path, "partitions.txt")
    else:
        file_path = os.path.join(dir_path, file_name)
    with open(file_path, "w+") as file:
        current_site_num = 1
        for part in sorted_parts:
            model_str = part["MODEL"]
            if part["ALPHA"] != "None":
                model_str += "+G"
            if part["INVARIANT_SITES"] != "None" and part["RAXML_NG"] == 0:
                model_str += "+I"
            elif part["PROPORTION_INVARIANT_SITES_STR"] != "None":
                invar_str = part['PROPORTION_INVARIANT_SITES_STR']
                invar_str = invar_str.split("{")[0]
                model_str += f"+{invar_str}"
            line = f"{part['DATA_TYPE']}, {model_str}, partition_{part['PARTITION_NUM']} = " \
                   f"{current_site_num}-{current_site_num + int(part['NUM_ALIGNMENT_SITES']) - 1}\n"
            current_site_num += int(part["NUM_ALIGNMENT_SITES"])
            file.write(line)


def get_repo_info_from_url(url):
    repo_dir = os.path.dirname(url)
    repo_name = ".".join(os.path.basename(url).split(".")[:-1])
    repo_owner = os.path.basename(repo_dir)

    return repo_owner, repo_name


def download_trees(dest_path, grouped_result, meta_info_dict, keys=[], amount=0, forced_out_dir="", source_dir=""):
    """
    Downloads a data set from GitHub, also can save the tree dict (result) in the destination directory
    @param dest_path: destination directory path
    @param grouped_result: may be passed if the tree dict is to be saved as well
    @param meta_info_dict: meta information about the repository to download the datasets from
    @param keys: if set, will use the tree ids in that list to download trees (instead of iterating over keys in
                 grouped_results)
    @param amount: amount of data sets to download (if result dict > amount)
    @param forced_out_dir: will use that directory name if set, otherwise uses the tree id as the directory name
    @param source_dir: database directory, in case a local source is used
    @return: list of paths to the downloaded data sets
    """

    commit_hash = meta_info_dict["COMMIT_HASH"]

    returned_paths = []
    tree_keys = keys if keys else list(grouped_result.keys())
    try:
        for i in range(amount):
            current_index = i % len(tree_keys)
            tree_id = tree_keys[current_index]
            dct = grouped_result[tree_id]

            print(f"downloading {tree_id}")
            dir_path = os.path.join(dest_path, tree_id) if not forced_out_dir else os.path.join(dest_path,
                                                                                                forced_out_dir)
            possible_files = [
                "tree_best.newick", "tree_part.newick", "log_0.txt", "model_0.txt", "iqt.pr_ab_matrix", "msa.fasta",
                f"{tree_id}.tar.gz"
            ]
            create_dir_if_needed(dir_path)
            if grouped_result:
                save_tree_dict(dir_path, grouped_result[tree_id])

                if len(grouped_result[tree_id]) > 1:
                    write_partitions_file(dir_path, grouped_result[tree_id])

            if not source_dir:
                for file_name in possible_files:
                    try:
                        if "URL" in meta_info_dict:
                            repo_owner, repo_name = get_repo_info_from_url(meta_info_dict["URL"])
                            link = BASE_GITHUB_LINK.format(repo_owner, repo_name, commit_hash, tree_id, file_name)
                        else:
                            link = BASE_GITHUB_LINK.format(BASE_GITHUB_REPO_OWNER, BASE_GITHUB_REPO_NAME,
                                                           commit_hash, tree_id, file_name)
                        with urlopen(link) as webpage:
                            content = webpage.read()
                        with open(os.path.join(dir_path, file_name), "wb+") as output:
                            output.write(content)
                    except Exception as e:
                        pass
            else:
                for file_name in possible_files:
                    try:
                        file_path = os.path.join(source_dir, tree_id, file_name)
                        dest_file_path = os.path.join(dir_path, file_name)
                        shutil.copy(file_path, dest_file_path)
                    except Exception as e:
                        pass

            returned_paths.append(dir_path)
    except Exception as e:
        print("Error while downloading: {}".format(e))
        print(traceback.print_exc())

    return returned_paths


def count_tree_leaves(clade):
    """
    Recursively determines the numbers of leaves in a tree and also carries branch lengths
    @param clade: root of the current subtree
    @return: number of leaves, list of branch lengths
    """
    leaf_counter = 0
    branch_length_list = [clade.branch_length] if clade.branch_length else []
    if clade.name:
        leaf_counter = 1
    for c in clade.clades:
        nl, bl = count_tree_leaves(c)
        leaf_counter += nl
        branch_length_list.extend(bl)
    return leaf_counter, branch_length_list


def get_tree_info(src_path, tree_id):
    """
    Creates a dict with statistical information about a tree with the help of Genesis
    @param src_path: path to tree file (newick format)
    @param tree_id: unique id of the tree
    @return: tree dict
    """
    global global_num_of_too_big_trees

    num_leaves = 0
    ret_dct = {}
    try:
        tree = Phylo.read(src_path, "newick")
        num_leaves, branch_length_list = count_tree_leaves(tree.root)
        try:
            diamcalc = GenesisTreeDiameter(GENESIS_PATH)  # TODO: maybe don't do it with genesis
            tree_len, tree_diam, tree_height = diamcalc.get_len_and_diam_and_height(src_path)
        except Exception as e:
            print("Genesis exception: {}".format(e))
            diamcalc = -1
            tree_len, tree_diam, tree_height = -1

        ret_dct["TREE_ID"] = tree_id  # TODO: remember to change!
        ret_dct["NUM_TAXA"] = num_leaves
        ret_dct["TREE_LENGTH"] = tree_len
        ret_dct["TREE_DIAMETER"] = tree_diam
        ret_dct["TREE_HEIGHT"] = tree_height
        ret_dct["BRANCH_LENGTH_MEAN"] = statistics.mean(branch_length_list)
        ret_dct["BRANCH_LENGTH_VARIANCE"] = statistics.variance(branch_length_list)
        ret_dct["IS_INDELIBLE_COMPATIBLE"] = 1
        for bl in branch_length_list:
            if "e" in str(bl):
                ret_dct["IS_INDELIBLE_COMPATIBLE"] = 0
                break

    except Exception as e:
        print("Exception in get_tree_info: {}".format(e))
        raise e
    return ret_dct


def file_exists(path, substr):
    """
    Checks if a file with a specific substring in its name can be found
    @param path: path to file
    @param substr: substring that the desired file must contain
    @return: True if found, otherwise False
    """
    files = os.listdir(path)
    for file_path in files:
        if substr in file_path:
            return True
    return False


def read_pr_ab_matrix(path):
    """
    Reads the presence/absence matrix and returns the number of 0s and 1s
    as well as the matrix itself (currently as list of lists...)
    @param path: path to matrix file
    @return: number of 0s, number of 1s in the matrix, and the matrix itself
    """
    with open(path) as file:
        first_line = True
        num_bits = 0
        num_0 = 0
        num_1 = 0
        matrix = []

        for line in file:
            line_spl = line.rstrip().split()
            if first_line:
                first_line = False
                num_bits = int(line_spl[1])
                continue

            temp_list = []

            for i in range(1, num_bits + 1):
                val = int(line_spl[-i])
                if val == 0:
                    num_0 += 1
                elif val == 1:
                    num_1 += 1

                temp_list.append(val)
            temp_list.reverse()
            matrix.append(temp_list)

    return num_0, num_1, matrix


def read_pr_ab_matrix2(path):
    entries = collections.defaultdict(lambda: collections.defaultdict(lambda: 1))
    with open(path) as file:
        first_line = True
        num_bits = 0

        for line in file:
            split_line = line.strip().split()
            if first_line:
                first_line = False
                num_bits = int(split_line[1])
                continue

            taxon = split_line[0]
            for i in range(1, num_bits + 1):
                entries[i - 1][taxon] = int(split_line[i])
    return entries


def group_partitions_in_result_dicts(results):
    """
    Groups partitions of the same tree
    @param results: as returned by db.find()
    @return: dict which maps tree ids to lists of the partitions of that tree
    """
    grouped_results = {}
    for result in results:
        tree_id = result["TREE_ID"]
        if tree_id in grouped_results:
            grouped_results[tree_id].append(result)
        else:
            grouped_results[tree_id] = [result]
    return grouped_results


def filter_incomplete_groups(grouped_results):
    """
    Filters data sets with incomplete partition sets (incomplete due to some partitions being excluded with the
    SQL query).
    @param grouped_results: dict mapping tree id to list of tree partition dicts
    @return:
    """
    filtered_grouped_results = {}
    for id in grouped_results:
        num_parts = grouped_results[id][0]["OVERALL_NUM_PARTITIONS"]
        if len(grouped_results[id]) == num_parts:
            filtered_grouped_results[id] = grouped_results[id]
    return filtered_grouped_results


def get_msa_gap_rate(msa_path):
    records = SeqIO.parse(msa_path, BASE_DAWG_SEQ_FILE_FORMAT.lower())  # TODO: debug with other file formats than fasta
    num_gaps = 0
    seq_len = 0
    num_taxa = 0
    for record in records:
        num_taxa += 1
        seq = record.seq
        if not seq_len:
            seq_len = len(seq)
        num_gaps += seq.count("_") + seq.count("-")
    return num_gaps / (seq_len * num_taxa)


def assemble_sequences(path_list, out_path, matrix_path="", in_format=BASE_DAWG_SEQ_FILE_FORMAT.lower(),
                       out_format="fasta"):
    """
    Concatenates multiple MSAs into a single file. We use this to create a big MSA file out of separately generated
    per-partition-MSAs.
    @param path_list: list of MSA file paths
    @param out_path: output path of the assembled MSA file
    @param matrix_path: path to the presence/absence matrix to use
                        if set: sequence i in partition p will be filled with blank symbols if matrix[i, p] == 0
                        if not set: no missing data will be introduced
    @param in_format: format of the sequence files to be assembled
    @param out_format: format of the MSA file
    @return: path of assembled msa file
    """
    # out_path = os.path.join(out_dir, f"assembled_sequences.fasta")
    with open(out_path, "w+") as file:
        file.write("")

    records_list = []
    if os.path.isfile(matrix_path):
        _, _, matrix = read_pr_ab_matrix(matrix_path)
    else:
        matrix = []

    for path in path_list:
        records = SeqIO.parse(path, in_format)
        records_list.append(list(records))

    for i in range(len(records_list[0])):
        sequence = ""
        for j in range(len(path_list)):
            bit = matrix[i][j] if matrix else 1

            if bit:
                sequence += str(records_list[j][i].seq)
            else:
                sequence += BLANK_SYMBOL * len(str(records_list[j][i].seq))
        new_rec = SeqRecord.SeqRecord(Seq.Seq(sequence), id=records_list[0][i].id, description="")

        with open(out_path, "a+") as file:
            SeqIO.write(new_rec, file, out_format)

    return out_path


def fix_msa_for_iqt(msa_path):
    """
    DEPRECATED and probably shouldn't be used
    Iterates over all sequences of the MSA. If a sequence is empty, sets the first character of that sequence to 'A'.
    This is done due to the fact that IQTREE2 would exit with an error if empty sequences are present in the MSA.
    @param msa_path: path to the msa file to be "fixed"
    @return:
    """

    records = SeqIO.parse(msa_path, BASE_DAWG_SEQ_FILE_FORMAT.lower())
    new_records = []
    rewrite = False

    for record in records:
        sequence = list(str(record.seq))
        if sequence.count(BLANK_SYMBOL) == len(sequence):
            sequence[0] = "A"
            rewrite = True
        new_rec = SeqRecord.SeqRecord(Seq.Seq("".join(sequence)), id=record.id, description="")
        new_records.append(new_rec)

    if rewrite:
        with open(msa_path, "w+") as file:
            for new_rec in new_records:
                SeqIO.write(new_rec, file, "fasta")


def fix_msa_for_sparta(source_path, dest_path, data_type="DNA"):
    sequences = msa_parser.parse_and_fix_msa_somehow(source_path, data_type)
    msa_parser.save_msa(sequences, dest_path, msa_format="fasta")


def check_for_empty_matrix_sequences(matrix_path):
    with open(matrix_path) as file:
        first_line = True
        seq_len = 0
        for line in file:
            line = line.strip()
            if first_line:
                seq_len = int(line.split()[-1])
                first_line = False
                continue
            if not line:
                continue
            sequence = "".join(line.split()[1:])
            if sequence.count("0") == seq_len:
                return True
    return False


def check_for_empty_sequences(msa_path):
    counter = 0
    records = SeqIO.parse(msa_path, BASE_DAWG_SEQ_FILE_FORMAT.lower())
    for record in records:
        if record.seq.count("-") == len(record.seq):
            counter += 1
    return counter


def blank_sequences_in_partition(part_seq_path, part_num, pr_ab_matrix_path):
    """
    Fills sequence s[i][p], for partition p and sequence index i,
    with blank symbols if presence/absence matrix m[i][p] == 0
    @param part_seq_path: path to partition MSA
    @param part_num: partition index
    @param pr_ab_matrix_path: path to presence/absence matrix
    @return:
    """
    if pr_ab_matrix_path:
        _, _, pr_ab_matrix = read_pr_ab_matrix(pr_ab_matrix_path)
        records = list(SeqIO.parse(part_seq_path, BASE_DAWG_SEQ_FILE_FORMAT.lower()))
        new_records = []
        for rec_num in range(len(records)):
            record = records[rec_num]
            if pr_ab_matrix[rec_num][part_num] == 1:
                new_records.append(record)
            else:
                sequence = BLANK_SYMBOL * len(record.seq)
                new_record = SeqRecord.SeqRecord(Seq.Seq(sequence), id=record.id, description="")
                new_records.append(new_record)
        with open(part_seq_path, "w+") as file:
            for new_rec in new_records:
                SeqIO.write(new_rec, file, BASE_DAWG_SEQ_FILE_FORMAT.lower())


def split_msa(msa_path, part_path="", msa_format="fasta"):
    msa_dir = os.path.dirname(msa_path)
    sequences = msa_parser.parse_msa_somehow(msa_path)
    seq_len = len(sequences[0].sequence)
    names = [s.id for s in sequences]

    part_dict = {}
    part_sites = collections.defaultdict(lambda: {})
    part_seqs = collections.defaultdict(lambda: collections.defaultdict(lambda: []))
    if not os.path.isfile(part_path):
        temp_msa_path = os.path.join(msa_dir, f"temp_msa_part_0.{msa_format}")
        part_dict[0] = temp_msa_path
        msa_parser.save_msa(sequences, temp_msa_path, msa_format=msa_format)
    else:
        with open(part_path) as file:
            for line in file:
                line = line.strip()
                model_string, all_site_string = line.split("=")
                part_id = int(model_string.split(",")[-1].split("_")[-1])
                site_strings = [s.strip() for s in all_site_string.split(",")]
                for site_string in site_strings:
                    if "-" in site_string:
                        start_n, end_n = map(int, site_string.split("-"))
                        for i in range(start_n, end_n + 1):
                            part_sites[part_id][i] = 1
                    elif "/" in site_string:
                        offset, step = map(int, site_string.split("/"))
                        for i in range(offset, seq_len, step=step):
                            part_sites[part_id][i] = 1
                    else:
                        site_n = int(site_string)
                        part_sites[part_id][site_n] = 1
        for sequence in sequences:
            seq = sequence.sequence
            name = sequence.id

            for part_id in part_sites:
                part_seq = "".join([seq[i-1] if i in part_sites[part_id] else "" for i in range(1, len(seq)+1)])
                part_seqs[part_id][name] = part_seq

        part_len_sum = sum([len(part_seqs[part_id][names[0]]) for part_id in part_seqs])
        if part_len_sum != seq_len:
            raise ValueError(f"something went wrong during MSA split (lengths mismatch: {part_len_sum} vs {seq_len})!")

        for part_id in part_seqs:
            temp_seqs = [msa_parser.Sequence(name, part_seqs[part_id][name], comment="") for name in names]
            temp_msa_path = os.path.join(msa_dir, f"temp_msa_part_{part_id}.{msa_format}")

            part_dict[part_id] = temp_msa_path
            msa_parser.save_msa(temp_seqs, temp_msa_path, msa_format=msa_format)

    return part_dict

def fix_tree_for_indelible(source_path, dest_path):
    """
    INDELible has strange issues reading newick tree files. We try to avoid some of these.
    @param source_path: source tree path
    @param dest_path: dest path of "fixed" tree
    @return:
    """
    def traverse_and_fix_branches(clade):
        if clade.branch_length:
            clade.branch_length = float(clade.branch_length)
        for c in clade.clades:
            traverse_and_fix_branches(c)

    with open(source_path) as file:
        for line in file:
            handle = StringIO(line)
    tree = Phylo.read(handle, "newick")
    traverse_and_fix_branches(tree.root)

    Phylo.write(tree, dest_path, BASE_TREE_FORMAT, plain=False, format_branch_length="%.16f")
    with open(dest_path) as file:
        new_tree_str = file.read()
    new_tree_str = new_tree_str.replace(":0.0000000000000000;", ";")
    with open(dest_path, "w+") as file:
        file.write(new_tree_str)


class SimMSAEval:
    def __init__(self, part_dict, tree_path, msa_out_path, reference_features, distance_class, generator=None, weights=[],
                 data_type="DNA", matrix_path=""):

        if generator:
            self.generator = generator
        else:
            self.generator = AliSim(ALISIM_PATH)
        self.part_dict = part_dict
        self.tree_path = tree_path
        self.msa_out_path = msa_out_path
        self.msa_dir = os.path.dirname(msa_out_path)
        self.reference_features = reference_features
        self.distance_obj = distance_class()
        self.weights = weights
        self.data_type = data_type
        self.matrix_entries = None
        if matrix_path:
            self.matrix_entries = read_pr_ab_matrix2(matrix_path)

        self.min_dist = BIG_NUMBER
        self.step_counter = 0

    @classmethod
    def using_sequences(cls, part_dict, tree_path, msa_out_path, sequences, distance_class, **kwargs):
        reference_features = distance_class().features(sequences)
        return cls(part_dict, tree_path, msa_out_path, reference_features, distance_class, **kwargs)


    def evaluate(self, sim_params, num_repeats=5):
        temp_msa_path = os.path.join(self.msa_dir, "temp_sim_msa_eval.fasta")
        temp_dists = []

        print(f"step {self.step_counter}: {sim_params}")
        self.step_counter += 1
        for _ in range(num_repeats):
            sim_error = self.generator.execute(
                self.tree_path,
                temp_msa_path,
                self.part_dict,
                new_seq_len=sim_params[0],
                indel_rates=(sim_params[1], sim_params[2]),
                indel_distr=("POW", "POW"),
                indel_distr_params=((sim_params[3], 500), (sim_params[4], 500))
            )
            dist = "N/A"
            try:
                if sim_error:
                    temp_dists.append(BIG_NUMBER)
                else:
                    temp_sequences = msa_parser.parse_and_fix_msa_somehow(temp_msa_path, self.data_type)
                    if self.matrix_entries:
                        for sequence in temp_sequences:
                            if self.matrix_entries[self.part_dict["PARTITION_NUM"]][sequence.id] == 0:
                                sequence.sequence = "-" * len(sequence.sequence)

                    temp_features = self.distance_obj.features(temp_sequences)
                    dist = self.distance_obj.distance(self.reference_features, temp_features, self.weights)

                    if dist < self.min_dist:
                        self.min_dist = dist
                        shutil.copy(temp_msa_path, self.msa_out_path)
                        print("    new min!")
                    temp_dists.append(dist)

                print(f"    {dist}")
            except Exception as e:
                temp_dists.append(BIG_NUMBER)
                print(f"    {dist}")

        return statistics.mean(temp_dists)


def simulate_msa_with_sparta(part_dict, msa_path, tree_path, msa_out_path, generator):
    fixed_tree_path = tree_path + "_mod"
    fix_tree_for_indelible(tree_path, fixed_tree_path)
    fixed_msa_path = msa_path + "_mod"
    fix_msa_for_sparta(msa_path, fixed_msa_path)

    msa_dir = os.path.dirname(msa_path)
    command = [
        "--path", msa_dir,
        "--msaf", fixed_msa_path,
        "--trf", fixed_tree_path,
        "--mode", "nuc",     # TODO: expand this eventually
        "--freq", str((part_dict["FREQ_A"], part_dict["FREQ_C"], part_dict["FREQ_G"], part_dict["FREQ_T"])),
        "--submodel", "GTR",
        "--rates", str((part_dict["RATE_AC"], part_dict["RATE_AG"], part_dict["RATE_AT"], part_dict["RATE_CG"],
                   part_dict["RATE_CT"])),
        "--nsim", f"{SPARTA_SIM_NUM}",     # TODO: adjust this
        "--nburnin", f"{SPARTA_BURNIN_NUM}",
         "--only-indelible",
    ]
    try:
        sparta.run(command)
    except Exception as e:
        print(e)

    sparta_out_path = os.path.join(msa_dir, "SpartaABC_data_name_iddif.posterior_params")   # TODO: we use the RIM stuff
    rl = r_i = r_d = a_i = a_d = 0
    weights = []
    with open(sparta_out_path) as file:
        for line in file:
            line = line.strip()
            if line.startswith("The posterior expectation for the root length (RL) is:"):
                rl = int(float(line.split(":")[-1]))
            if line.startswith("The posterior expectation for the A Insertion parameter is:"):
                a_i = float(line.split(":")[-1])
            if line.startswith("The posterior expectation for the A Deletion parameter is:"):
                a_d = float(line.split(":")[-1])
            if line.startswith("The posterior expectation for the indel rate (InsR) is:"):
                r_i = float(line.split(":")[-1])
            if line.startswith("The posterior expectation for the indel rate (DelR) is:"):
                r_d = float(line.split(":")[-1])
            if line.startswith("weights"):
                weights = [float(s) for s in line.split()[1:]]

    sequences = msa_parser.parse_and_fix_msa_somehow(msa_path, part_dict["DATA_TYPE"])
    eval = SimMSAEval.using_sequences(part_dict, tree_path, msa_out_path, sequences, weights=weights,
                      distance_class=sparta_util.SpartaDist, generator=generator)
    eval.evaluate([rl, r_i, r_d, a_i, a_d])


def simulate_msa_with_extended_sparta(part_dict, msa_path, tree_path, msa_out_path, generator,
                           dist_function=sparta_util.sparta_extended_dist):

    rl = part_dict["NUM_ALIGNMENT_SITES"]
    rl_s = 0.5
    num_rounds = EXTENDED_SPARTA_SIM_NUM

    search_space = [
        space.Integer(int(rl * (1 - rl_s)), rl, name="rl"),
        space.Real(0.0, 0.05, name="r_i"),    # TODO: try log-uniform distributions
        space.Real(0.0, 0.05, name="r_d"),
        space.Real(1.001, 2, name="a_i"),
        space.Real(1.001, 2, name="a_d")
    ]

    # TODO: fix that mess with weights somehow
    weights = []
    sequences = msa_parser.parse_and_fix_msa_somehow(msa_path, part_dict["DATA_TYPE"])
    evaluator = SimMSAEval.using_sequences(part_dict, tree_path, msa_out_path, sequences, weights=weights,
                                      distance_class=sparta_util.ExtendedSpartaDist, generator=generator)
    gp_res = gp_minimize(evaluator.evaluate, search_space, n_calls=num_rounds)


def simulate_msa_with_bonk(part_dict, tree_path, msa_out_path, generator, matrix_path=""):
    rl = part_dict["NUM_ALIGNMENT_SITES"]
    rl_s = 0.7
    num_rounds = BONK_SIM_NUM

    search_space = [
        space.Integer(int(rl * (1 - rl_s)), rl, name="rl"),
        space.Real(0.0, 0.05, name="r_i"),    # TODO: try log-uniform distributions
        space.Real(0.0, 0.05, name="r_d"),
        space.Real(1.001, 2, name="a_i"),
        space.Real(1.001, 2, name="a_d")
    ]

    # TODO: fix that mess with weights somehow
    weights = [9, 10, 1]
    dist_class = sparta_util.BlindDist
    reference_features = dist_class().features_from_part_dict(part_dict)
    evaluator = SimMSAEval(part_dict, tree_path, msa_out_path, reference_features, dist_class,
                           generator=generator, weights=weights, data_type=part_dict["DATA_TYPE"],
                           matrix_path=matrix_path)
    gp_res = gp_minimize(evaluator.evaluate, search_space, n_calls=num_rounds)


def create_generator(args):
    final_alignment_path = "{}"
    if args.generator == "seq-gen":
        raise ValueError("seq-gen currently not working!")
        generator = SeqGen(SEQGEN_PATH)
        alignment_file_ext = f".{BASE_DAWG_SEQ_FILE_FORMAT}"
    elif args.generator == "dawg":
        generator = Dawg(DAWG_PATH)
        alignment_file_ext = f".{BASE_DAWG_SEQ_FILE_FORMAT}"
    else:
        generator = AliSim(ALISIM_PATH)
        alignment_file_ext = f""

    if args.seed:
        generator.set_seed(args.seed)
        print(f"Seed: {args.seed}")
    return generator, alignment_file_ext, final_alignment_path


def generate_sequences(grouped_results, args, meta_info_dict, forced_out_dir=""):
    generator, alignment_file_ext, final_alignment_path = create_generator(args)
    print(f"Using {args.generator}.")

    out_dir = os.path.join(os.path.abspath(args.out_dir))
    create_dir_if_needed(out_dir)

    returned_paths = []
    returned_results = {}

    rand_key = random.choice(list(grouped_results.keys()))
    rand_tree_data = grouped_results[rand_key][0]

    print(rand_tree_data)

    # Download tree
    tree_dir_name = forced_out_dir
    if not forced_out_dir:
        if not args.force_rewrite:
            tree_dir_name = find_unused_tree_folder_name(out_dir, rand_tree_data["TREE_ID"])
        else:
            tree_dir_name = os.path.join(out_dir, rand_tree_data["TREE_ID"])
        dl_tree_path = os.path.join(out_dir, tree_dir_name)
    else:
        dl_tree_path = os.path.join(out_dir, forced_out_dir)

    # TODO: currently only one data set is downloaded, maybe change it again in future
    returned_paths.extend(download_trees(out_dir, grouped_results, meta_info_dict,
                                         keys=[rand_key],
                                         amount=1, forced_out_dir=tree_dir_name,
                                         source_dir=args.use_local_db))
    if args.no_simulation:
        write_partitions_file(dl_tree_path, grouped_results[rand_key], "sim_partitions.txt")
        returned_results[rand_key] = grouped_results[rand_key]
        return returned_results, returned_paths

    num_retries = 1
    pr_ab_matrix_path = os.path.join(dl_tree_path, "iqt.pr_ab_matrix")
    # TODO: maybe also use pr/absence matrices with the sparta stuff
    if not os.path.isfile(pr_ab_matrix_path) or args.use_sparta or args.use_modified_sparta:
        pr_ab_matrix_path = ""

    assembled_msa_path = os.path.join(dl_tree_path, "assembled_sequences.fasta")
    backup_seq_part_path_tuples = []
    sorted_seq_paths = []

    if args.avoid_empty_sequences:
        num_retries = 3
    for i in range(num_retries):
        no_indel_partition_saved = False
        no_indel_partition_file = "no_indel_partition_{}.fasta"
        backup_seq_part_path_tuples = []

        returned_results[rand_key] = grouped_results[rand_key]
        tree_path = os.path.join(dl_tree_path, BASE_TREE_NAME.format("best", BASE_TREE_FORMAT))

        seq_part_path_tuples = []
        partitions = grouped_results[rand_key]
        part_msa_paths = {}

        if args.use_sparta or args.use_modified_sparta:
            # TODO: rework. if msa.fasta is not available, the (currently not existing) feature vector
            #  from db should be used
            msa_path = os.path.join(dl_tree_path, "msa.fasta")
            model_path = os.path.join(dl_tree_path, "model_0.txt")
            if not os.path.isfile(msa_path) or not os.path.isfile(model_path):
                print(f"cannot use sparta since original msa or model file not found for tree {rand_key}!")
                return

            part_msa_paths = split_msa(msa_path, part_path=model_path, msa_format="fasta")

        for part in partitions:
            if args.seq_len:
                tmp_nsites = max(args.seq_len // part["OVERALL_NUM_PARTITIONS"], 1)
                part["NUM_ALIGNMENT_SITES"] = tmp_nsites
                if part["PARTITION_NUM"] in [0, "0"]:
                    part["NUM_ALIGNMENT_SITES"] = tmp_nsites + (args.seq_len % part["OVERALL_NUM_PARTITIONS"])

            if part["PARTITION_NUM"] == "None":
                part["PARTITION_NUM"] = 0
            current_part_num = part["PARTITION_NUM"]

            # TODO: do something with the formats...
            seq_part_path = os.path.join(dl_tree_path, f"seq_{i}.part{current_part_num}.fasta")
            formatted_seq_path = final_alignment_path.format(seq_part_path)

            # generate an MSA in any case without weights in a "dry run", even if we reoptimize later
            generator.execute(tree_path, seq_part_path, part, timeout=SIMULATION_TIMEOUT)
            seq_part_path_tuples.append((formatted_seq_path, current_part_num))

            # TODO: this should be (maybe) reworked eventually!
            # we save the first "dry run" partition. in case --avoid-empty-sequences is used and we fail to
            # simulate a MSA without completely empty sequences after multiple attempts, we use this
            # first partition without indels instead of the normally simulated partition MSA with indels.
            if not no_indel_partition_saved and args.avoid_empty_sequences:
                no_indel_partition_file = no_indel_partition_file.format(current_part_num)
                no_indel_partition_file = os.path.join(dl_tree_path, no_indel_partition_file)
                shutil.copy(formatted_seq_path, no_indel_partition_file)
                no_indel_partition_saved = True

                backup_seq_part_path_tuples.append((no_indel_partition_file, current_part_num))
            else:
                backup_seq_part_path_tuples.append((formatted_seq_path, current_part_num))

            if args.use_sparta:
                simulate_msa_with_sparta(part, part_msa_paths[current_part_num], tree_path, seq_part_path, generator)
            elif args.use_modified_sparta:
                simulate_msa_with_extended_sparta(part, part_msa_paths[current_part_num], tree_path, seq_part_path,
                                                  generator, dist_function=sparta_util.sparta_extended_dist)
            elif args.use_bonk or args.weights:
                simulate_msa_with_bonk(part, tree_path, seq_part_path, generator, matrix_path=pr_ab_matrix_path)
            else:
                generator.execute(tree_path, seq_part_path, part, timeout=SIMULATION_TIMEOUT)

        # Assemble MSAs
        seq_part_path_tuples.sort(key=lambda x: x[1])
        sorted_seq_paths = [x[0] for x in seq_part_path_tuples]
        assemble_sequences(sorted_seq_paths, assembled_msa_path, matrix_path=pr_ab_matrix_path)

        if args.avoid_empty_sequences:
            if not check_for_empty_sequences(assembled_msa_path):
                break
        else:
            break

    if args.avoid_empty_sequences:
        if check_for_empty_sequences(assembled_msa_path):
            backup_seq_part_path_tuples.sort(key=lambda x: x[1])
            sorted_seq_paths = [x[0] for x in backup_seq_part_path_tuples]
            assemble_sequences(sorted_seq_paths, assembled_msa_path, matrix_path=pr_ab_matrix_path)

            with open(os.path.join(dl_tree_path, "failed_to_insert_gaps"), "w+") as file:
                file.write("1\n")

    final_part_list = []
    part_id_map = {}
    for part in grouped_results[rand_key]:
        part_id_map[int(part["PARTITION_NUM"])] = part

    for i, seq_path in enumerate(sorted_seq_paths):
        seqs = msa_parser.parse_msa_somehow(seq_path)
        msa_len = len(seqs[0].sequence)
        part_dct = copy.deepcopy(part_id_map[i])
        part_dct["NUM_ALIGNMENT_SITES"] = msa_len
        final_part_list.append(part_dct)

    write_partitions_file(dl_tree_path, final_part_list, "sim_partitions.txt")

    return returned_results, returned_paths


def get_archive_meta_data(archive_path):
    """
    Collects relevant meta data about the RAxMLGrove (RG) repository to put into the SQLite db.
    Currently we require the GitPython package for this function to work. If GitPython is not installed,
    a warning will be printed. We use this function to keep track of the different commits to RG, which might
    change the mapping of TREE_ID in SQLite db to the directories in RG, by saving the commit hash and using it
    to access directories in the RG repository on GitHub (when files are being downloaded).
    @param archive_path: path to RG root directory
    @return: dict with meta data
    """
    meta_dict = {}

    try:
        import git  # GitPython import here, as we use it for this functionality only
        repo = git.Repo(archive_path, search_parent_directories=True)
        commit_hash = repo.head.object.hexsha
        url = repo.remotes[0].config_reader.get("url")

        meta_dict["COMMIT_HASH"] = commit_hash
        meta_dict["URL"] = url
        print(commit_hash)
    except Exception as e:
        print("WARNING: the script was unable to get the current commit hash of the RG repository (supplied with -a).\n"
              "The commit hash will not be noted in the SQLite db, which might mess up the mapping of TREE_IDs in the "
              "SQLite db to the data set directories of the RG repository on GitHub at some point in future "
              "(or even now).")

    return meta_dict


def hopefully_somewhat_better_directory_crawl(root_path, db_object, add_new_files_only=False, local=False):
    """
    Crawls the RAxML Grove directory and creates a dict with tree information for every job (for "create" command).
    Since this function is recursive, we just store the dicts globally in global_tree_dict_list and in
        global_part_list_dict (maybe not the most beautiful way, but one of the simplest)
    @param root_path: archive path as passed by -a
    @param db_object: our standard db_object
    @param add_new_files_only: if True:
                                    only adds new entries (trees with not present ids) to the db
                                if False:
                                    deletes the old db and creates a new one (if --force-rewrite is used) or just
                                    creates a db if db not present
    @param local: deprecated flag which was once used to create local databases
                  (which did not download the trees from git)
    @return: Nothing...
    """

    global global_tree_dict_list
    global global_part_list_dict

    global global_exception_counter
    global global_num_of_checked_jobs
    global global_max_tree_file_len

    files = os.listdir(root_path)

    tree_dicts = []
    file_dict = {}
    is_rax_ng = False

    dir_counter = 0
    for file_path in files:
        current_path = os.path.join(root_path, file_path)

        if os.path.isdir(current_path):
            """dir_counter += 1
            if dir_counter > 100:
                return"""
            hopefully_somewhat_better_directory_crawl(current_path, db_object, add_new_files_only, local)
        else:
            if "RAxML_bestTree" in file_path:
                file_dict["BEST_TREE"] = current_path
            elif file_path.endswith(".raxml.bestTree") or "tree_best.newick" in file_path:
                file_dict["BEST_TREE"] = current_path
            elif "RAxML_info" in file_path:
                file_dict["INFO"] = current_path
            elif ".raxml.bestPartitionTrees" in file_path or "tree_part.newick" in file_path:
                file_dict["PART_TREES"] = current_path
            elif ".raxml.bestModel" in file_path or "model_0" in file_path:
                file_dict["BEST_MODEL"] = current_path
                is_rax_ng = True
            elif ".raxml.log" in file_path or "log_0" in file_path:
                file_dict["INFO"] = current_path
            elif "iqt.pr_ab_matrix" == file_path:
                file_dict["PR_AB_MATRIX"] = current_path

    try:
        if is_rax_ng:
            if "BEST_TREE" not in file_dict or \
                    "INFO" not in file_dict or \
                    "BEST_MODEL" not in file_dict:
                return

            if local:
                tree_id = file_dict["BEST_TREE"]
            else:
                tree_id = os.path.basename(os.path.dirname(file_dict["BEST_TREE"]))

            log_reader = RaxmlNGLogReader(file_dict["INFO"], file_dict["BEST_MODEL"])
            tree_info = get_tree_info(file_dict["BEST_TREE"], tree_id)
            partitions_info = log_reader.get_partition_info()

            tree_info["RAXML_NG"] = 1
            tree_info["OVERALL_NUM_ALIGNMENT_SITES"] = 0
            tree_info["OVERALL_NUM_PATTERNS"] = 0
            tree_info["OVERALL_GAPS"] = 0

            i = 0
            for key in partitions_info:
                temp_dict = {
                    "PARENT_ID": tree_info["TREE_ID"],
                    "PARTITION_NUM": i if len(partitions_info) > 1 else None
                }
                i += 1
                temp_dict.update(partitions_info[key])
                tree_dicts.append(temp_dict)

                tree_info["OVERALL_NUM_ALIGNMENT_SITES"] += partitions_info[key]["NUM_ALIGNMENT_SITES"] \
                    if "NUM_ALIGNMENT_SITES" in partitions_info[key] and partitions_info[key][
                    "NUM_ALIGNMENT_SITES"] else 0
                tree_info["OVERALL_NUM_PATTERNS"] += partitions_info[key]["NUM_PATTERNS"] \
                    if "NUM_PATTERNS" in partitions_info[key] and partitions_info[key]["NUM_PATTERNS"] else 0
                tree_info["OVERALL_GAPS"] += partitions_info[key]["GAPS"] \
                    if "GAPS" in partitions_info[key] and partitions_info[key]["GAPS"] else 0

            tree_info["OVERALL_GAPS"] = tree_info["OVERALL_GAPS"] / len(partitions_info)
            tree_info["OVERALL_NUM_PARTITIONS"] = len(partitions_info)

            global_num_of_checked_jobs += 1

        else:
            if "BEST_TREE" not in file_dict or \
                    "INFO" not in file_dict:
                return

            if local:
                tree_id = file_dict["BEST_TREE"]
            else:
                tree_id = os.path.basename(os.path.dirname(file_dict["BEST_TREE"]))

            log_reader = OldRaxmlReader(file_dict["INFO"])
            tree_info = get_tree_info(file_dict["BEST_TREE"], tree_id)
            partitions_info = log_reader.get_partition_info()

            tree_info["RAXML_NG"] = 0
            tree_info["OVERALL_NUM_ALIGNMENT_SITES"] = 0
            tree_info["OVERALL_NUM_PATTERNS"] = 0
            tree_info["OVERALL_GAPS"] = 0

            i = 0
            for key in partitions_info:
                temp_dict = {
                    "PARENT_ID": tree_info["TREE_ID"],
                    "PARTITION_NUM": i if len(partitions_info) > 1 else None
                }
                i += 1
                temp_dict.update(partitions_info[key])
                tree_dicts.append(temp_dict)

                tree_info["OVERALL_NUM_ALIGNMENT_SITES"] += partitions_info[key]["NUM_ALIGNMENT_SITES"] \
                    if "NUM_ALIGNMENT_SITES" in partitions_info[key] and partitions_info[key][
                    "NUM_ALIGNMENT_SITES"] else 0
                tree_info["OVERALL_NUM_PATTERNS"] += partitions_info[key]["NUM_PATTERNS"] \
                    if "NUM_PATTERNS" in partitions_info[key] and partitions_info[key]["NUM_PATTERNS"] else 0
                tree_info["OVERALL_GAPS"] = partitions_info[key]["GAPS"] \
                    if "GAPS" in partitions_info[key] and partitions_info[key]["GAPS"] else 0

            tree_info["OVERALL_NUM_PARTITIONS"] = len(partitions_info)

            global_num_of_checked_jobs += 1

        if "PR_AB_MATRIX" in file_dict:
            num_0, num_1, _ = read_pr_ab_matrix(file_dict["PR_AB_MATRIX"])
            tree_info["MISSING_DATA_RATE"] = num_0 / (num_0 + num_1)

    except Exception as e:
        print(f"Exception in directory crawl at {root_path}:\n    {e}")
        print(traceback.print_exc())
        global_exception_counter += 1
        return

    if global_num_of_checked_jobs % 1000 == 0:
        print("________________________\n{}\n________________________".format(global_num_of_checked_jobs))
        print(f"exceptions: {global_exception_counter}")

    global_tree_dict_list.append(tree_info)
    global_part_list_dict[tree_info["TREE_ID"]] = tree_dicts


def print_statistics(db_object, query):
    """
    Prints statistical information about the entries of the columns in the db (stats operation). If -q is used,
    stats of the filtered results are printed.
    @param db_object: our standard database object
    @param query: query, as expected by -q argument
    """

    results = db_object.find(f"SELECT * FROM TREE t INNER JOIN PARTITION p ON t.TREE_ID = p.PARENT_ID WHERE {query};")
    grouped_results = group_partitions_in_result_dicts(results)
    grouped_results = filter_incomplete_groups(grouped_results)

    cat_float_values = {}
    cat_str_values = {}

    tree_columns = []
    for col, _ in COLUMNS:
        tree_columns.append(col)

    for tree_id in grouped_results:
        first = True
        for tree_dict in grouped_results[tree_id]:
            for cat in tree_dict:
                value = tree_dict[cat]
                if value == "None" or cat in ["TREE_ID", "PARENT_ID", "IS_INDELIBLE_COMPATIBLE", "RAXML_NG",
                                              "PROPORTION_INVARIANT_SITES_STR", "PARTITION_NUM"]:
                    continue

                if cat in tree_columns and not first:
                    continue

                if cat in ["MODEL", "DATA_TYPE"]:
                    if cat in cat_str_values:
                        cat_str_values[cat].append(value)
                    else:
                        cat_str_values[cat] = [value]
                else:
                    try:
                        cvalue = float(value)
                        if cat in cat_float_values:
                            cat_float_values[cat].append(cvalue)
                        else:
                            cat_float_values[cat] = [cvalue]
                    except Exception as e:
                        pass
            first = False

    print(f"Number of trees: {len(list(grouped_results.keys()))}")
    print()

    for cat in cat_float_values:
        lower_fence, upper_fence = get_tukeys_fences(cat_float_values[cat])
        filtered_values = list(filter(lambda x: lower_fence <= x <= upper_fence, cat_float_values[cat]))

        perc = 0.95
        lp = 0
        hp = perc
        sorted_values = sorted(cat_float_values[cat])

        print(f"{cat}:\n"
              f"    min {min(cat_float_values[cat])} max {max(cat_float_values[cat])}\n"
              f"    mean {statistics.mean(cat_float_values[cat])} "
              f"median {statistics.median(cat_float_values[cat])}")

        if len(filtered_values) > 0:
            print(f"  without outliers (k=1.5):\n"
                  f"    lower/upper fence {(lower_fence, upper_fence)}\n"
                  f"    min {min(filtered_values)} max {max(filtered_values)}\n"
                  f"    mean {statistics.mean(filtered_values)} "
                  f"median {statistics.median(filtered_values)}"
                  )
            hp_index = int(len(sorted_values) * hp) + 1
            if hp_index == len(sorted_values) and len(sorted_values) > 0:
                hp_index = hp_index - 1

            print(f"  {perc}:\n"
                  f"    low {sorted_values[int(len(sorted_values) * lp)]}  high {sorted_values[hp_index]}\n"
                  )

    for cat in cat_str_values:
        str_buckets = {}
        str_buckets_list = []
        for entry in cat_str_values[cat]:
            if entry in str_buckets:
                str_buckets[entry] += 1
            else:
                str_buckets[entry] = 1
        for key in str_buckets:
            str_buckets_list.append(
                (key, str_buckets[key], "{0:.2f}%".format(str_buckets[key] * 100 / len(cat_str_values[cat]))))
        str_buckets_list.sort(key=lambda x: x[1], reverse=True)

        print(f"{cat}: {str_buckets_list}\n")


def main(args_list, is_imported=True):
    """
    Main function. We do stuff depending on the selected operation.
    Throughout most of the code, db_object will contain our SQLite database and answer queries. The returned
    values will be called "results". results are to contain a list of dicts of trees with their respective
    attributes. We use two tables:
        TREE for the overall tree information
        PARTITION to store information of different partitions for one TREE entry
    So we have a many-to-one relation between PARTITIONs and TREEs. When we return joined results, we will have a
    dict for every partition in the database, which is why we usually group results
    (with group_partitions_in_result_dicts()) by tree ids before usage.

    The operations
        "create", "add", "execute" were made for database creation/manipulation, and should probably not
    be very interesting to most users.
        "find" will try to find data sets according to the passed query "-q query".
        "generate" will find data sets which were inferred under the GTR substitution model, download the data sets,
    and generate sequences for these sets (according to the other arguments/flags)
        "stats" prints statistical information about the data sets in the database.
    """

    args = init_args(args_list)

    db_path = args.db_name if args.db_name else DEFAULT_DB_FILE_PATH

    force_rewrite = True if args.operation == "create" and args.force_rewrite else False
    db_object = BetterTreeDataBase(db_path, force_rewrite=force_rewrite)
    meta_info_dict = {
        "COMMIT_HASH": ""
    }

    # We return found data sets / data sets which were used for the generation of MSAs, in case main() is called
    # directly. returned_paths contains the directory paths the data sets were downloaded to.
    returned_results = {}
    returned_paths = []

    if args.operation == "create" or args.operation == "add":
        archive_path = args.archive_path
        if not args.rg_commit_hash:
            meta_info_dict = get_archive_meta_data(archive_path)
        else:
            meta_info_dict = {"COMMIT_HASH": args.rg_commit_hash}

        # Crawl the RAxMLGrove archive, parse RAxML output files, fill global_tree_dict_list and global_part_list_dict,
        # write them into a SQLite database afterwards.
        hopefully_somewhat_better_directory_crawl(archive_path, db_object, add_new_files_only=(args.operation == "add"),
                                                  local=False)
        db_object.fill_database(global_tree_dict_list, global_part_list_dict, meta_info_dict)

        print("\nExceptions: {}".format(global_exception_counter))
        print("Num too big trees: {}".format(global_num_of_too_big_trees))
    elif args.operation == "execute":
        db_object.execute_command(args.command)
    elif args.operation == "find":
        if not args.use_local_db:
            meta_info_dict = db_object.get_meta_info()
        if args.rg_commit_hash:
            meta_info_dict["COMMIT_HASH"] = args.rg_commit_hash

        result = db_object.find(
            f"SELECT * FROM TREE t INNER JOIN PARTITION p ON t.TREE_ID = p.PARENT_ID WHERE {args.query};")
        grouped_result = group_partitions_in_result_dicts(result)
        grouped_result = filter_incomplete_groups(grouped_result)

        num_results = len(list(grouped_result.keys()))

        printed_results = 0
        printed_dcts = {}

        if not args.list:
            for dct in result:
                current_tree_id = dct["TREE_ID"]
                if current_tree_id not in printed_dcts:
                    print(dct)
                    printed_dcts[dct["TREE_ID"]] = 1
                    printed_results += 1

                if printed_results >= 10:
                    print(f"...({num_results - printed_results} more)...")
                    break

            print("\nNumber of results: {}".format(num_results))

            if num_results > 0:
                ask_download = yes_no_prompt("Would you like to download these trees?")
                if ask_download:
                    tree_dest_dir = args.out_dir or input('Destination (default cwd): ')

                    amount_to_download = input('How many trees to download (default all): ')
                    returned_paths = download_trees(tree_dest_dir, grouped_result, meta_info_dict,
                                                    keys=[],
                                                    amount=(
                                                        int(amount_to_download) if amount_to_download else num_results),
                                                    source_dir=args.use_local_db)
        else:
            if not is_imported:
                for dct in result:
                    if dct["TREE_ID"] not in printed_dcts:
                        print(dct)
                        printed_dcts[dct["TREE_ID"]] = 1

        returned_results = grouped_result

    elif args.operation == "download":
        if not args.use_local_db:
            meta_info_dict = db_object.get_meta_info()
        if args.rg_commit_hash:
            meta_info_dict["COMMIT_HASH"] = args.rg_commit_hash

        result = db_object.find(
            f"SELECT * FROM TREE t INNER JOIN PARTITION p ON t.TREE_ID = p.PARENT_ID WHERE {args.query};")
        grouped_result = group_partitions_in_result_dicts(result)
        grouped_result = filter_incomplete_groups(grouped_result)

        returned_paths = download_trees(args.out_dir, grouped_result, meta_info_dict,
                                        amount=len(grouped_result), source_dir=args.use_local_db)

    elif args.operation == "generate" or args.operation == "justgimmeatree":
        if not args.use_local_db:
            meta_info_dict = db_object.get_meta_info()
        if args.rg_commit_hash:
            meta_info_dict["COMMIT_HASH"] = args.rg_commit_hash  # TODO: warn if hash was present and is overwritten?

        if not args.filter_outliers:
            query = args.query if args.query else "OVERALL_NUM_ALIGNMENT_SITES > 0"  # TODO: expand possible models
            results = db_object.find(
                f"{BASE_SQL_FIND_COMMAND} WHERE OVERALL_NUM_ALIGNMENT_SITES > 0 AND {query};")
        else:
            # Categories we currently filter outliers for
            categories = {  # TODO: make this work for AA (and all the other stuff)
                "NUM_TAXA": 0,
                "OVERALL_NUM_ALIGNMENT_SITES": 0,
                #"TREE_DIAMETER": 0,
                #"BRANCH_LENGTH_VARIANCE": 0,
                #"RATE_AC": 0,
                #"RATE_AG": 0,
                #"RATE_AT": 0,
                #"RATE_CG": 0,
                #"RATE_CT": 0,
                #"RATE_GT": 0,
                #"FREQ_A": 0,
                #"FREQ_C": 0,
                #"FREQ_G": 0,
                #"FREQ_T": 0,
                #"ALPHA": 0,
                #"NUM_PATTERNS": 0,
                #"GAPS": 0
                # "TREE_LENGTH": 0     # TODO: check if available for all
            }
            all_tree_data = db_object.find(f"{BASE_SQL_FIND_COMMAND};")
            filter_list = []
            for cat in categories:
                column_data = [entry[cat] for entry in all_tree_data]
                categories[cat] = get_tukeys_fences(column_data, k=1.5)

                filter_list.append(f"{cat} >= {categories[cat][0]} AND {cat} <= {categories[cat][1]}")
            if args.query:
                query = args.query + " AND " + " AND ".join(filter_list)
            else:
                query = " AND ".join(filter_list)
            results = db_object.find(
                f"{BASE_SQL_FIND_COMMAND} WHERE OVERALL_NUM_ALIGNMENT_SITES > 0 AND {query};")

        grouped_results = group_partitions_in_result_dicts(results)
        grouped_results = filter_incomplete_groups(grouped_results)

        print("Found {} trees".format(len(list(grouped_results.keys()))))

        if len(results) > 0:
            if args.operation == "generate":
                returned_results, returned_paths = generate_sequences(grouped_results, args, meta_info_dict)
            else:
                returned_results, returned_paths = generate_sequences(grouped_results, args, meta_info_dict,
                                                                      forced_out_dir="default")
    elif args.operation == "stats":
        print_statistics(db_object, args.query if args.query else "1")

    db_object.close()

    return returned_results, returned_paths


if __name__ == "__main__":
    if len(sys.argv) > 1:
        main(sys.argv[1:], is_imported=False)
    else:
        main(["-h"])
