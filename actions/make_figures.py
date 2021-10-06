#!/usr/bin/env python3
import statistics
import sys
import os
import sqlite3
import pandas as pd
import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# sys.path.append("../org_script")      # TODO: fix this
# from org_script import get_tukeys_fences, BetterTreeDataBase

BASE_FILE_DIR = os.path.dirname(__file__)
BASE_README_IMG_FORMAT = '<img src="{}" width="{}%"></img>'
CSV_PATH = "./latest.csv"


class BetterTreeDataBase(object):
    def __init__(self, path, force_rewrite=False):
        """if os.path.isabs(path):
            self.db_path = path
        else:
            self.db_path = os.path.join(BASE_FILE_DIR, path)"""
        self.db_path = path

        self.conn = sqlite3.connect(self.db_path)
        self.conn.row_factory = sqlite3.Row
        self.cursor = self.conn.cursor()

    def close(self):
        self.conn.close()

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


def rewrite_readme(figure_path_list):
    out_lines = []
    readme_path = os.path.join(BASE_FILE_DIR, "..", "README.md")
    with open(readme_path) as file:
        for line in file:
            line = line.rstrip()
            if "./figures/" not in line and "Some Figures" not in line:
                out_lines.append(line)

    with open(readme_path, "w+") as file:
        for line in out_lines:
            file.write(line + "\n")
        file.write("## Some Figures\n\n")
        for path in figure_path_list:
            file.write(BASE_README_IMG_FORMAT.format(path, 45))


def get_values_for_cat(lst, cat):
    ret = []
    for dct in lst:
        ret.append(dct[cat])
    return ret


def is_float(value):
    try:
        c = float(value)
    except Exception as e:
        return False
    return True


def do_test(name):
    db_path = os.path.join(BASE_FILE_DIR, "..", name)

    if not os.path.isfile(CSV_PATH):
        export_to_csv(db_path)

    df = pd.read_csv(CSV_PATH)
    db_object = BetterTreeDataBase(db_path)
    plt.ioff()

    columns = [
        "NUM_TAXA",
        "TREE_DIAMETER",
        "TREE_LENGTH",
        "TREE_HEIGHT",
        "BRANCH_LENGTH_VARIANCE",
        "ALPHA",
        "NUM_ALIGNMENT_SITES",
        "GAPS",
        "RATE_AC",
        "RATE_AG",
        "RATE_AT",
        "RATE_CG",
        "RATE_CT",
        "RATE_GT",
        "FREQ_A",
        "FREQ_C",
        "FREQ_G",
        "FREQ_T"
    ]
    num_buckets = 100
    scale = 1
    fig_paths = []

    for cat in columns:
        raw_values = df.loc[:, cat]
        raw_values = [float(v) for v in list(filter(lambda x: x != "None" and is_float(x), raw_values))]

        low_fence, high_fence = get_tukeys_fences(raw_values)
        results = db_object.find(
            f"SELECT * FROM TREE t INNER JOIN PARTITION p ON t.TREE_ID = p.PARENT_ID WHERE {cat} >= {low_fence} AND {cat} <= {high_fence} "
            f"AND (PARTITION_NUM == 'None' OR PARTITION_NUM == 0);")
        values = get_values_for_cat(results, cat)

        vmin = min(values)
        vmax = max(values)
        bucket_size = (vmax - vmin) / num_buckets

        fig = plt.figure()
        bins = [bucket_size * x for x in range(num_buckets)]
        # plt.hist(np.clip(raw_values, bins[0], bins[-1]), bins=bins)
        plt.hist(values, bins=bins)

        plt.xlabel(cat)
        plt.ylabel("Number of trees")
        fig_path = f"./figures/test_{cat}.png"
        plt.savefig(fig_path)
        plt.close(fig)

        fig_paths.append(fig_path)

    rewrite_readme(fig_paths)


def export_to_csv(db_path):
    db_path = os.path.join(BASE_FILE_DIR, "..", db_path)
    print(db_path)
    try:
        conn = sqlite3.connect(db_path, isolation_level=None, detect_types=sqlite3.PARSE_COLNAMES)
        db_df = pd.read_sql("SELECT * FROM TREE t INNER JOIN PARTITION p ON t.TREE_ID = p.PARENT_ID "
                            "WHERE (PARTITION_NUM == 'None' OR PARTITION_NUM == 0)", conn)
    except Exception as e:
        print(e)
        print("-------t2--------")
        try:
            conn = sqlite3.connect("test", isolation_level=None, detect_types=sqlite3.PARSE_COLNAMES)
            db_df = pd.read_sql("SELECT * FROM TREE t INNER JOIN PARTITION p ON t.TREE_ID = p.PARENT_ID "
                                "WHERE (PARTITION_NUM == 'None' OR PARTITION_NUM == 0)", conn)
        except Exception as e:
            print(e)
            print("-------t3-------")
            conn = sqlite3.connect("test.db", isolation_level=None, detect_types=sqlite3.PARSE_COLNAMES)
            db_df = pd.read_sql("SELECT * FROM TREE t INNER JOIN PARTITION p ON t.TREE_ID = p.PARENT_ID "
                                "WHERE (PARTITION_NUM == 'None' OR PARTITION_NUM == 0)", conn)

    db_df.to_csv(CSV_PATH, index=False)


def main():
    if len(sys.argv) != 2:
        print("Illegal number of arguments. Run as 'script.py path_to_db'")
        exit(1)

    db_path = sys.argv[1]

    if not os.path.isfile(CSV_PATH):
        export_to_csv(db_path)
    do_test(db_path)


if __name__ == "__main__":
    main()
