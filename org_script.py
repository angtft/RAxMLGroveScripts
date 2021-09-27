#!/usr/bin/env python3

import sys
import os
import sqlite3
import statistics
import subprocess
import argparse
import random
import re
import copy
from io import StringIO
from urllib.request import urlopen
from Bio import Phylo, SeqIO, Seq, SeqRecord

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
BASE_NODE_NAME = "taxon{}"

"""COLUMNS = [
    ("TREE_ID", "CHAR(255)"), ("NUM_TAXA", "INT"), ("TREE_LENGTH", "FLOAT"), ("TREE_DIAMETER", "FLOAT"),
    ("MODEL", "CHAR(50)"), ("ALPHA", "FLOAT"), ("RATE_AC", "FLOAT"), ("RATE_AG", "FLOAT"), ("RATE_AT", "FLOAT"),
    ("RATE_CG", "FLOAT"), ("RATE_CT", "FLOAT"), ("RATE_GT", "FLOAT"), ("FREQ_A", "FLOAT"), ("FREQ_C", "FLOAT"),
    ("FREQ_G", "FLOAT"), ("FREQ_T", "FLOAT"), ("BRANCH_LENGTH_MEAN", "FLOAT"), ("BRANCH_LENGTH_VARIANCE", "FLOAT"),
    ("IS_INDELIBLE_COMPATIBLE", "INT"), ("NUM_ALIGNMENT_SITES", "INT"), ("NUM_PATTERNS", "INT"),
    ("GAPS", "FLOAT"), ("INVARIANT_SITES", "FLOAT"), ("RAXML_NG", "INT"), ("DATA_TYPE", "CHAR(50)"),
    ("RATE_STR", "cHAR(5000)"), ("FREQ_STR", "CHAR(2000)"), ("PARTITION_NUM", "INT"),
    ("STATIONARY_FREQ_STR", "CHAR(100)"), ("PROPORTION_INVARIANT_SITES_STR", "CHAR(100)"),
    ("AMONG_SITE_RATE_HETEROGENEITY_STR", "CHAR(100)"), ("ASCERTAINMENT_BIAS_CORRECTION_STR", "CHAR(100)"),
    ("CUSTOM_CHAR_TO_STATE_MAPPING", "CHAR(100)")  # , ("IS_ON_TERRACE", "INT")
    # TODO: get that is_on_terrace bool (newer versions of RAxML-NG?)
    # TODO: test AA trees!
]"""
COLUMNS = [
    ("TREE_ID", "CHAR(255)"), ("NUM_TAXA", "INT"), ("TREE_LENGTH", "FLOAT"), ("TREE_DIAMETER", "FLOAT"),
    ("BRANCH_LENGTH_MEAN", "FLOAT"), ("BRANCH_LENGTH_VARIANCE", "FLOAT"),
    ("IS_INDELIBLE_COMPATIBLE", "INT"), ("OVERALL_NUM_ALIGNMENT_SITES", "INT"), ("OVERALL_NUM_PATTERNS", "INT"),
    ("OVERALL_GAPS", "FLOAT"), ("INVARIANT_SITES", "FLOAT"), ("RAXML_NG", "INT"),
    ("OVERALL_NUM_PARTITIONS", "INT"), ("MISSING_DATA_RATE", "FLOAT")   # TODO: missing data using presence/absence matrices
]

PARTITION_COLUMNS = [
    ("MODEL", "CHAR(50)"), ("ALPHA", "FLOAT"), ("RATE_AC", "FLOAT"), ("RATE_AG", "FLOAT"), ("RATE_AT", "FLOAT"),
    ("RATE_CG", "FLOAT"), ("RATE_CT", "FLOAT"), ("RATE_GT", "FLOAT"), ("FREQ_A", "FLOAT"), ("FREQ_C", "FLOAT"),
    ("FREQ_G", "FLOAT"), ("FREQ_T", "FLOAT"), ("NUM_ALIGNMENT_SITES", "INT"), ("NUM_PATTERNS", "INT"),
    ("GAPS", "FLOAT"), ("INVARIANT_SITES", "FLOAT"), ("DATA_TYPE", "CHAR(50)"),
    ("RATE_STR", "CHAR(5000)"), ("FREQ_STR", "CHAR(2000)"), ("PARTITION_NUM", "INT"),
    ("STATIONARY_FREQ_STR", "CHAR(100)"), ("PROPORTION_INVARIANT_SITES_STR", "CHAR(100)"),
    ("AMONG_SITE_RATE_HETEROGENEITY_STR", "CHAR(100)"), ("ASCERTAINMENT_BIAS_CORRECTION_STR", "CHAR(100)"),
    ("CUSTOM_CHAR_TO_STATE_MAPPING", "CHAR(100)"), ("PARENT_ID", "INT")
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

BASE_GITHUB_LINK = "https://raw.githubusercontent.com/angtft/RAxMLGrove/main/trees/{}/{}"
BASE_DB_FILE_NAME = "latest.db"
BASE_FILE_DIR = os.path.dirname(__file__)
BASE_OUT_DIR = os.path.join(BASE_FILE_DIR, "out")
BASE_SEQ_FILE_FORMAT = "Phylip"  # TODO: maybe think some more about formats (since Phylip only allows taxon names up to 10 characters, and SeqGen doesn't output other formats?)
BASE_DAWG_SEQ_FILE_FORMAT = "Fasta"
BASE_SQL_FIND_COMMAND = "SELECT * FROM TREE t INNER JOIN PARTITION p ON t.TREE_ID = p.PARENT_ID"

global_num_partitioned_trees = 0  # TODO: remove this!


def get_tukeys_fences(lst):
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

        k = 1.5

        low_fence = q1 - k * iqr
        high_fence = q3 + k * iqr
    except Exception as e:
        raise e

    return low_fence, high_fence


def traverse_and_rename_nodes(clade):
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
    with open(path) as file:
        tree_string = file.read().rstrip()
    return tree_string


def create_dir_if_needed(path):
    try:
        os.makedirs(path)
    except OSError as e:
        pass


def local_tree_copy(dest_dir, results, amount=1):
    if type(results) is list:
        tree_pile = results
    else:
        tree_pile = [results]
    try:
        for i in range(amount):
            dct = tree_pile[i]
            tree_path = dct["TREE_ID"]  # TREE_ID should contain the absolute path
            copy_tree_file(tree_path, os.path.join(dest_dir, BASE_TREE_NAME.format(i, BASE_TREE_FORMAT)))
    except Exception as e:
        print("Error in local_tree_copy: {}".format(e))


class Dawg(object):
    """
    We use Dawg here (https://github.com/reedacartwright/dawg)
    """

    def __init__(self, path):
        self.path = path
        self.template_path = os.path.join(self.path, "examples", "template.dawg")
        self.config_path = os.path.join(self.path, "examples", "template_modified.dawg")
        self.execute_path = os.path.join(self.path, "build", "src", "dawg")
        if not os.path.isfile(self.execute_path):
            self.__compile()

    def execute(self, tree_path, out_path, tree_params, num_repeats=1, num_of_sequence=0):
        self.__configure(tree_path, tree_params)
        try:
            if num_repeats > 1:
                for i in range(0, num_repeats):
                    call = [self.execute_path, self.config_path, "-o",
                            # os.path.join(out_path, f"seq_{i}.{BASE_DAWG_SEQ_FILE_FORMAT}")]
                            out_path]
                    subprocess.check_output(call, cwd=BASE_FILE_DIR)
            else:
                call = [self.execute_path, self.config_path, "-o",
                        # os.path.join(out_path, f"seq_{num_of_sequence}.{BASE_DAWG_SEQ_FILE_FORMAT}")]
                        out_path]
                subprocess.check_output(call, cwd=BASE_FILE_DIR)
        except Exception as e:
            print(e)
            exit(0)

    def __compile(self):
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

    def __configure(self, tree_path, tree_params):
        tree_string = read_tree(tree_path)
        seq_len = tree_params["NUM_ALIGNMENT_SITES"]
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
                else:
                    out_lines.append(line)
        with open(self.config_path, "w+") as config_file:
            config_file.write("".join(out_lines))


class SeqGen(object):
    def __init__(self, path):
        self.path = path
        self.execute_path = os.path.join(self.path, "source", "seq-gen")
        if not os.path.isfile(self.execute_path):
            self.__compile()

    def execute(self, tree_path, out_path, tree_params, num_repeats=1, num_of_sequence=0):
        out_path_abs = out_path if os.path.isabs(out_path) else os.path.join(BASE_FILE_DIR, out_path)
        try:
            call = [self.execute_path, tree_path, f"-m{tree_params['MODEL']}",
                    f"-f{tree_params['FREQ_A']},{tree_params['FREQ_C']},{tree_params['FREQ_G']},{tree_params['FREQ_T']}",
                    f"-r{tree_params['RATE_AC']},{tree_params['RATE_AG']},{tree_params['RATE_AT']},{tree_params['RATE_CG']},{tree_params['RATE_CT']},{tree_params['RATE_GT']}",
                    f"-l{tree_params['NUM_ALIGNMENT_SITES']}", "-or"]
            if num_of_sequence == 0:
                for i in range(0, num_repeats):
                    out = subprocess.check_output(call, cwd=BASE_FILE_DIR, stderr=subprocess.DEVNULL)
                    with open(os.path.join(out_path_abs, f"seq_{i}.{BASE_SEQ_FILE_FORMAT}"), "w+") as file:
                        file.write(out.decode("UTF-8"))
            else:
                out = subprocess.check_output(call, cwd=BASE_FILE_DIR, stderr=subprocess.DEVNULL)
                with open(os.path.join(out_path_abs, f"seq_{num_of_sequence}.{BASE_SEQ_FILE_FORMAT}"), "w+") as file:
                    file.write(out.decode("UTF-8"))
        except Exception as e:
            print("Exception in SeqGen.execute(): {}".format(e))
            exit(0)

    def __compile(self):
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


class BetterTreeDataBase(object):
    def __init__(self, path, force_rewrite=False):
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

    def fill_database(self, tree_dict_list, partition_list_dict):
        for key in partition_list_dict:
            for part in partition_list_dict[key]:
                for entry, _ in PARTITION_COLUMNS:
                    if entry not in part:
                        part[entry] = None

        print(tree_dict_list)
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
                print("Exception in fill_database: {}".format(e))
                continue
        self.conn.commit()

    def database_entry_exists(self, id):
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
        # Allow only reading access here!
        self.cursor.execute(command)

    def find(self, command):
        self.cursor.execute(command)
        result = self.cursor.fetchall()
        return [dict(row) for row in result]


class RaxmlNGLogReader(object):
    def __init__(self, path, model_path):
        self.path = path
        self.partitions_dict = {}
        self.__read()
        self.__read_model(model_path)
        self.__fill_general_info()

    def _read_fix(self):
        # TODO: implement
        with open(self.path) as file:
            pass

    def __read(self):
        with open(self.path) as file:
            part_info_dict = {}
            for line in file:
                line = line.strip()
                if line.startswith("Partition"):
                    if "name" in part_info_dict:
                        self.partitions_dict[part_info_dict["name"]].update(part_info_dict)
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
                if line.startswith("Rate heterogeneity:") and "NONE" not in line:  # TODO: add the other params
                    part_info_dict["ALPHA"] = float(line.split()[-7])
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
        modifier_dict = {
            "STATIONARY_FREQ_STR": ["F", "FC", "FO", "FE", "FU"],  # stationary frequencies
            "PROPORTION_INVARIANT_SITES_STR": ["I", "IO", "IC", "IU"],  # proportion of invariant sites
            "AMONG_SITE_RATE_HETEROGENEITY_STR": ["G", "G4m", "R"],  # among-site rate heterogeneity model
            "ASCERTAINMENT_BIAS_CORRECTION_STR": ["ASC_LEWIS", "ASC_FELS", "ASC_STAM"],  # ascertainment bias correction
            "CUSTOM_CHAR_TO_STATE_MAPPING": ["M", "Mi"]
        }

        with open(path) as file:
            lines = file.readlines()
            # assert len(lines) == len(self.partitions_dict)    # TODO: maybe do this check...

            i = 0
            for part_key in self.partitions_dict:
                try:
                    line = lines[i]
                except Exception as e:
                    print(f"Exception in raxml-ng __read_model: {self.path}\n{e}")
                    continue
                modifiers_info = line.rstrip().split("+")

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
                            self.partitions_dict[part_key][category] = mi
                i += 1

    def __fill_general_info(self):
        """
        use after __read() and __read_model() were executed
        """
        for part_key in self.partitions_dict:
            if "ASCERTAINMENT_BIAS_CORRECTION_STR" in self.partitions_dict[part_key] and \
                    "G" in self.partitions_dict[part_key]["ASCERTAINMENT_BIAS_CORRECTION_STR"]:
                self.partitions_dict[part_key]["ALPHA"] = \
                self.__get_values_from_modifiers(self.partitions_dict[part_key]["ASCERTAINMENT_BIAS_CORRECTION_STR"])[0]

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
    def __init__(self, path):
        self.path = path
        self.partitions_dict = {}
        self.__read()
        self.__fill_model_info()

    def __read(self):
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
                    try:
                        part_line = value[1].strip() + ","
                        intervals = re.findall(r"(.*?),[ ]*", part_line)

                        summands = re.findall(r"(\d+?)[-\\]+", intervals[0] + "-")
                        s1 = int(summands[0])
                        s2 = int(summands[1])
                        part_size = s2 - s1 + 1
                        temp_part_dict["NUM_ALIGNMENT_SITES"].append(part_size)
                    except Exception as e:
                        print(f"Exception in old_raxml __read: {self.path}\n{e}")
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

            for i in range(num_partitions):
                new_part = {}
                for key in temp_part_dict:
                    if temp_part_dict[key] and len(temp_part_dict[key]) == num_partitions:
                        new_part[key] = temp_part_dict[key][i]
                self.partitions_dict[str(i)] = copy.deepcopy(new_part)

    def __fill_model_info(self):
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
                part["FREQ_STR"] = f"{{{'/'.join([str(x) for x in part['BASE_FREQUENCIES']])}}}"

    def get_partition_info(self):
        return self.partitions_dict


class GenesisTreeDiameter(object):  # TODO: Add genesis to ./tools/
    """
    We use Genesis here (https://github.com/lczech/genesis)
    """

    def __init__(self, path):
        self.path = path
        self.executable_path = os.path.join(self.path, "bin", "apps", "tree_diameter")
        if not os.path.isfile(self.executable_path):
            self.__compile()

    def __compile(self):
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
        print("Done!")

    def get_diam(self, tree_path1):
        call = [self.executable_path, tree_path1]

        try:
            out = subprocess.check_output(call).decode()
            return float(out.split()[-1])
        except Exception as e:
            print(e)
            return 1

    def get_len_and_diam(self, tree_path1):
        call = [self.executable_path, tree_path1]

        try:
            out = subprocess.check_output(call).decode()
            return float(out.split()[-2]), float(out.split()[-1])
        except Exception as e:
            print(e)
            return -1, -1


def init_args(arguments):
    parser = argparse.ArgumentParser()
    parser.add_argument("operation", choices=["create", "add", "execute", "find", "generate"],
                        # TODO: rework 'find' command?
                        help="'create' iterates over a raxml out files archive parametrized with '-a' and "
                             "writes a database (db) file (the name can be parametrized with '-n'). Default db name "
                             f"{BASE_DB_FILE_NAME}.\n"
                             "'add' adds trees from archive path to the db.\n"
                             "'execute' executes a command parametrized with '-c' on the db.\n"
                             "'find' tries to find trees in the db satisfying the conditions set with '-q'.\n"
                             "'generate' generates sequences with Dawg using randomly drawn trees from the db.")
    parser.add_argument("-n", "--db-name", help=f"Sets the name of the database file. If not set, this tool will "
                                                f"try to create/access the '{BASE_DB_FILE_NAME}' file.")
    parser.add_argument("-a", "--archive-path", help="Sets the raxml out files archive path. (create)")
    parser.add_argument("-c", "--command", help="The command to execute on the database. (execute)")
    parser.add_argument("-q", "--query", help="Part of the statement after the 'WHERE' clause "
                                              "to find trees in the db. CAUTION: We do not sanitize, "
                                              "don't break your own database...")
    parser.add_argument("--num-sequences", default=1, help="Amount of sequences to be generated. To generate "
                                                           "a sequence we randomly draw a tree from the db (can be "
                                                           "used with -q) and run a sequence generator with that tree. "
                                                           "(generate)")
    parser.add_argument("--local", default=True, help="FOR TESTING PURPOSES (don't use it).")
    parser.add_argument("--force-rewrite", action='store_true', help="Forces to rewrite the default db or the db "
                                                                     "specified using '-n'. (create)")
    parser.add_argument("--seq-len", default=8000, help="Default generated sequence length if it is not specified in "
                                                        "the tree log data. (generate)")
    parser.add_argument("--set-seq-len", help="Sets the generated sequence length IGNORING the sequence length "
                                              "specified in the tree log data. (generate)")
    parser.add_argument("--filter-outliers", action='store_true', help="Filters trees with uncommon characteristics "
                                                                       "(using Tukey's fences). (generate)")
    parser.add_argument("--insert-matrix-gaps", action="store_true", help="Uses the presence/absence matrices to "
                                                                          "insert gaps into the simulated sequences. "
                                                                          "(generate)")
    parser.add_argument("-g", "--generator", choices=["dawg", "seq-gen"], default="dawg",
                        help="Selects the sequence generator used. (generate)")
    parser.add_argument("-o", "--out-dir", default=BASE_OUT_DIR,
                        help=f"Output directory (default: {BASE_OUT_DIR}).")

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

    return args


def yes_no_prompt(message):
    print(message + " y/n")
    while True:
        user_in = input('>>> ')
        if user_in in ('y', 'yes'):
            return True
        elif user_in in ('n', 'no'):
            return False
        else:
            print('Please answer with "y" (yes) or "n" (no)!')


def download_trees(dest_path, result, amount=0):
    try:
        for i in range(amount):
            dct = result[i]
            tree_id = dct["TREE_ID"]
            dir_path = os.path.join(dest_path, tree_id)
            possible_files = [
                "tree_best.newick", "tree_part.newick", "log_0.txt", "model_0.txt", "iqt.pr_ab_matrix"
            ]
            create_dir_if_needed(dir_path)

            for file_name in possible_files:
                try:
                    with urlopen(BASE_GITHUB_LINK.format(tree_id, file_name)) as webpage:
                        content = webpage.read().decode()
                    with open(os.path.join(dir_path, file_name), "w+") as output:
                        output.write(content)
                except Exception as e:
                    pass
    except Exception as e:
        print("Error while downloading: {}".format(e))


def count_tree_leaves(clade):
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
    global global_num_of_too_big_trees

    num_leaves = 0
    ret_dct = {}
    try:
        tree = Phylo.read(src_path, "newick")
        num_leaves, branch_length_list = count_tree_leaves(tree.root)
        try:
            diamcalc = GenesisTreeDiameter("./tools/genesis-0.24.0/")  # TODO: don't do it with genesis
        except Exception as e:
            print("Genesis exception: {}".format(e))
            diamcalc = -1

        if num_leaves < 20000:  # TODO: fix(?). trees above this size make genesis eat too much ram :'(
            tree_len, tree_diam = diamcalc.get_len_and_diam(src_path)
        else:
            print("tree too big! num taxa: {}".format(num_leaves))
            global_num_of_too_big_trees += 1
            tree_len = tree_diam = -1

        ret_dct["TREE_ID"] = tree_id  # TODO: remember to change!
        ret_dct["NUM_TAXA"] = num_leaves
        ret_dct["TREE_LENGTH"] = tree_len
        ret_dct["TREE_DIAMETER"] = tree_diam
        ret_dct["BRANCH_LENGTH_MEAN"] = statistics.mean(branch_length_list)
        ret_dct["BRANCH_LENGTH_VARIANCE"] = statistics.variance(branch_length_list)
        ret_dct["IS_INDELIBLE_COMPATIBLE"] = 1
        for bl in branch_length_list:
            if "e" in str(bl):
                ret_dct["IS_INDELIBLE_COMPATIBLE"] = 0
                break

    except Exception as e:
        print("Exception in get_tree_info: {}".format(e))
        return -1
    return ret_dct


def file_exists(path, substr):
    files = os.listdir(path)
    for file_path in files:
        if substr in file_path:
            return True
    return False


def read_pr_ab_matrix(path):
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


def assemble_sequences(path_list, out_dir, matrix_path="", tree_id=0):
    out_path = os.path.join(out_dir, f"test_{tree_id}.fasta")
    with open(out_path, "w+") as file:
        file.write("")

    records_list = []
    if os.path.isfile(matrix_path):
        _, _, matrix = read_pr_ab_matrix(matrix_path)
    else:
        matrix = []

    for path in path_list:
        records = SeqIO.parse(path, BASE_DAWG_SEQ_FILE_FORMAT.lower())
        records_list.append(list(records))

    for i in range(len(records_list[0])):
        sequence = ""
        for j in range(len(path_list)):
            bit = matrix[i][j] if matrix else 1

            if bit:
                sequence += str(records_list[j][i].seq)
            else:
                sequence += "_" * len(str(records_list[j][i].seq))
        new_rec = SeqRecord.SeqRecord(Seq.Seq(sequence), id=records_list[0][i].id, description="")

        with open(out_path, "a+") as file:
            SeqIO.write(new_rec, file, "fasta")


def generate_sequences(results, args):
    if args.generator == "seq-gen":
        generator = SeqGen(os.path.join(BASE_FILE_DIR, "tools", "Seq-Gen-1.3.4"))
    else:
        generator = Dawg(os.path.join(BASE_FILE_DIR, "tools", "dawg-1.2"))
    print(f"Using {args.generator}.")

    temp_tree_dir = os.path.join(BASE_FILE_DIR, "temp")
    create_dir_if_needed(temp_tree_dir)

    grouped_results = {}
    for result in results:
        tree_id = result["TREE_ID"]
        if tree_id in grouped_results:
            grouped_results[tree_id].append(result)
        else:
            grouped_results[tree_id] = [result]

    for key in list(grouped_results.keys()):
        try:
            num_parts = grouped_results[key][0]["OVERALL_NUM_PARTITIONS"]
            if num_parts != len(grouped_results[key]):
                del grouped_results[key]
        except Exception as e:
            print(e)
            del grouped_results[key]

    # TODO: think of ways to fix the following
    for i in range(int(args.num_sequences)):
        rand_key = random.choice(list(grouped_results.keys()))   # TODO: currently we just randomly pick results out of the whole set. however, it could be better to
                                                                 #       iterate through all keys if num_sequences > len(grouped_results)
        rand_tree_data = grouped_results[rand_key][0]
        print(f"\nRun {i}. Selected tree:\n {rand_tree_data}")

        if rand_tree_data["NUM_ALIGNMENT_SITES"] == "None":
            rand_tree_data["NUM_ALIGNMENT_SITES"] = int(args.seq_len)
        if args.set_seq_len:
            rand_tree_data["NUM_ALIGNMENT_SITES"] = int(args.set_seq_len)

        download_trees(temp_tree_dir, [rand_tree_data], amount=1)
        dl_tree_path = os.path.join(temp_tree_dir, rand_tree_data["TREE_ID"])
        tree_path = os.path.join(dl_tree_path, BASE_TREE_NAME.format("best", BASE_TREE_FORMAT))
        if args.insert_matrix_gaps:
            pr_ab_matrix_path = os.path.join(dl_tree_path, "iqt.pr_ab_matrix")
        else:
            pr_ab_matrix_path = ""

        seq_part_paths = []
        partitions = grouped_results[rand_key]

        for p_num in range(len(partitions)):
            seq_part_path = os.path.join(dl_tree_path, f"seq_{i}.part{p_num}.{BASE_DAWG_SEQ_FILE_FORMAT}")  # TODO: do something with the formats...

            part = partitions[p_num]
            generator.execute(tree_path, seq_part_path, part)
            seq_part_paths.append(seq_part_path)

        assemble_sequences(seq_part_paths, out_dir=args.out_dir, matrix_path=pr_ab_matrix_path, tree_id=rand_tree_data["TREE_ID"])


def hopefully_somewhat_better_directory_crawl(root_path, db_object, add_new_files_only=False, local=False):
    global global_tree_dict_list
    global global_part_list_dict

    global global_exception_counter
    global global_num_of_checked_jobs
    global global_max_tree_file_len

    global global_num_partitioned_trees  # TODO: remove this

    files = os.listdir(root_path)

    tree_dicts = []
    file_dict = {}
    is_rax_ng = False

    for file_path in files:
        current_path = os.path.join(root_path, file_path)

        if os.path.isdir(current_path):
            hopefully_somewhat_better_directory_crawl(current_path, db_object, add_new_files_only, local)
        else:
            if "RAxML_bestTree" in file_path:
                file_dict["BEST_TREE"] = current_path
            elif ".raxml.bestTree" in file_path or "tree_best.newick" in file_path:
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

            if len(partitions_info) > 1:
                global_num_partitioned_trees += 1

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
                    if "NUM_ALIGNMENT_SITES" in partitions_info[key] and partitions_info[key]["NUM_ALIGNMENT_SITES"] else 0
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

            if len(partitions_info) > 1:
                global_num_partitioned_trees += 1

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
                    if "NUM_ALIGNMENT_SITES" in partitions_info[key] and partitions_info[key]["NUM_ALIGNMENT_SITES"] else 0
                tree_info["OVERALL_NUM_PATTERNS"] += partitions_info[key]["NUM_PATTERNS"] \
                    if "NUM_PATTERNS" in partitions_info[key] and partitions_info[key]["NUM_PATTERNS"] else 0
                tree_info["OVERALL_GAPS"] = partitions_info[key]["GAPS"] \
                    if "GAPS" in partitions_info[key] and partitions_info[key]["GAPS"] else 0

            tree_info["OVERALL_NUM_PARTITIONS"] = len(partitions_info)

            global_num_of_checked_jobs += 1

        if "PR_AB_MATRIX" in file_path:
            num_0, num_1, _ = read_pr_ab_matrix(file_dict["PR_AB_MATRIX"])
            tree_info["MISSING_DATA_RATE"] = num_0 / (num_0 + num_1)

    except Exception as e:
        print(f"Exception in directory crawl: {e}")
        return

    if global_num_of_checked_jobs % 1000 == 0:
        print("________________________\n{}\n________________________".format(global_num_of_checked_jobs))

    global_tree_dict_list.append(tree_info)
    global_part_list_dict[tree_info["TREE_ID"]] = tree_dicts


def main(args_list):
    args = init_args(args_list)

    db_path = args.db_name if args.db_name else BASE_DB_FILE_NAME
    db_object = BetterTreeDataBase(db_path, force_rewrite=args.force_rewrite)

    if args.operation == "create" or args.operation == "add":
        archive_path = args.archive_path

        hopefully_somewhat_better_directory_crawl(archive_path, db_object, add_new_files_only=(args.operation == "add"),
                                                  local=False)  # TODO: change local

        db_object.fill_database(global_tree_dict_list, global_part_list_dict)

        print("\nExceptions: {}".format(global_exception_counter))
        print("Num too big trees: {}".format(global_num_of_too_big_trees))
    elif args.operation == "execute":
        db_object.execute_command(args.command)
    elif args.operation == "find":
        result = db_object.find(f"SELECT * FROM TREE t INNER JOIN PARTITION p ON t.TREE_ID = p.PARENT_ID WHERE {args.query};")
        num_results = len(result)

        printed_results = 0
        for dct in result:
            printed_results += 1
            print(dct)

            if printed_results >= 100:
                print(f"...({num_results - printed_results} more)...")
                break

        print("\nNumber of results: {}".format(num_results))

        if num_results > 0:
            ask_download = yes_no_prompt("Would you like to download these trees?")
            if ask_download:
                tree_dest_dir = input('Destination (default cwd): ')
                amount_to_download = input('How many trees to download (default all): ')
                download_trees(tree_dest_dir, result,
                               amount=(int(amount_to_download) if amount_to_download else num_results))
                """local_tree_copy(tree_dest_dir, result,
                                amount=(int(amount_to_download) if amount_to_download else num_results))"""
    elif args.operation == "generate":
        print(f"Using {args.generator}.")

        if not args.filter_outliers:
            query = args.query if args.query else "MODEL LIKE 'GTR%' AND RATE_AC AND FREQ_A AND OVERALL_NUM_ALIGNMENT_SITES > 0"      # TODO: expand possible models
            results = db_object.find(f"{BASE_SQL_FIND_COMMAND} WHERE {query};")
        else:
            categories = {  # TODO: make this work for AA (and all the other stuff)
                "NUM_TAXA": 0,
                "TREE_DIAMETER": 0,
                "BRANCH_LENGTH_VARIANCE": 0,
                "RATE_AC": 0,
                "RATE_AG": 0,
                "RATE_AT": 0,
                "RATE_CG": 0,
                "RATE_CT": 0,
                "RATE_GT": 0,
                "FREQ_A": 0,
                "FREQ_C": 0,
                "FREQ_G": 0,
                "FREQ_T": 0,
                "ALPHA": 0,
                "NUM_PATTERNS": 0,
                "GAPS": 0
                # "TREE_LENGTH": 0     # TODO: check if available for all
            }
            all_tree_data = db_object.find(f"{BASE_SQL_FIND_COMMAND};")
            filter_list = []
            for cat in categories:
                column_data = [entry[cat] for entry in all_tree_data]
                categories[cat] = get_tukeys_fences(column_data)

                filter_list.append(f"{cat} >= {categories[cat][0]} AND {cat} <= {categories[cat][1]}")
            if args.query:
                query = args.query + " AND " + " AND ".join(filter_list)
            else:
                query = " AND ".join(filter_list)
            results = db_object.find(f"{BASE_SQL_FIND_COMMAND} WHERE MODEL LIKE 'GTR%' AND OVERALL_NUM_ALIGNMENT_SITES > 0 AND {query};")

        print("Found {} trees".format(len(set([r["TREE_ID"] for r in results]))))
        generate_sequences(results, args)

    db_object.close()


if __name__ == "__main__":
    if len(sys.argv) > 1:
        main(sys.argv[1:])

        print("Num part trees: {}".format(global_num_partitioned_trees))
    else:
        main(["-h"])
