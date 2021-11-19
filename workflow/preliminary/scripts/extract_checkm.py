import sys
import os
import pandas as pd
import re
import json

def extract_checkm_stats(file):
    
    with open(file, "r") as handle:
        f = handle.read()
        
    f = re.sub(".*{", "{", f)
    f = re.sub("'", "\"", f)
    f = json.loads(f)
    
    completeness = f["Completeness"]
    contamination = f["Contamination"]
    gc = f["GC"]
    genome_size = f["Genome size"]
    coding_density = f["Coding density"]
    
    value = [gc, genome_size, coding_density, completeness, contamination]
    return value

metadata = pd.read_csv(sys.argv[1], sep = "\t")

keys = metadata["id"].tolist()

df = []

for key in keys:
    checkm_file = "results/checkm/Bacteria_{}/storage/bin_stats_ext.tsv".format(key)
    
    answer = [key] + \
        extract_checkm_stats(checkm_file)
    
    df.append(answer)

df = pd.DataFrame(df)
df.columns = ["id", "GC", "Genome_Size", "Coding_Density", "Bacteria_Completeness", "Bacteria_Contamination"]
# df.set_index("Accession", inplace = True)

df = metadata.merge(df, left_on = "id", right_on = "id")

df.sort_values(["Bacteria_Completeness", "Bacteria_Contamination"],
    ascending = [False, True],
    inplace = True)
df.to_csv(sys.stdout, sep = "\t", index = False)
