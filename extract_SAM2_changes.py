import json
import numpy as np
import matplotlib.pyplot as plt 
import os
import pandas as pd
import csv
import math
from collections import defaultdict
import argparse
from pathlib import Path
import glob

json_path = "/Users/ellag/Desktop/PhD/academic_projects/3d_reconstruction/data/error/json_files"

sam2_correct = []
ii = 0
for filename in os.listdir(json_path):

    ii = ii + 1

    filename_path = os.path.join(json_path, filename)
    print(filename_path)

    if not filename.lower().endswith(".json"):
        continue

    #for each json file
    with open(filename_path, "r") as file:
        d =json.load(file)
        file.close()
        print(filename_path)

        data = []
        i = 0
        for items in d['annotations']:
            i = i + 1
            #print(items)
            image_id = items['image_id']
            iscrowd = items['iscrowd']
            row = [image_id, iscrowd]
            data.append(row)

        df = pd.DataFrame(data)
        df.columns = ["image_ID", "iscrowd"]

        prop_correct = df['iscrowd'].sum()/i

    row1 = [filename, prop_correct]
    sam2_correct.append(row1)

sam2_correct_df = pd.DataFrame(sam2_correct)
sam2_correct_df.columns = ["filename", "prop_correct"]
    


sam2_correct_df.to_csv('sam2_correct.csv')

