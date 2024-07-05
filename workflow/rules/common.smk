import glob
import pandas as pd
from snakemake.utils import validate


validate(config, schema="../schemas/config.schema.yaml")

windows = (
    pd.read_csv(config["windows"], sep=",", dtype={"chr": str, "start": int, "end": int, "name": str})
    .set_index("name", drop=False)
    .sort_index()
)
validate(windows,schema="../schemas/windows.schema.yaml")

wildcard_constraints:
    window="|".join(windows["name"])

def get_coords(wildcards):
    """Get coordinates of a given window."""
    coords = windows.loc[(wildcards.window), ["chr","start","end"]].dropna()
    thin = config["thin"]
    return f"--chr {coords['chr']} --from-bp {coords['start']} --to-bp {coords['end']} --thin {thin} --recode-INFO-all"

WINDOWS=windows['name'].tolist()

def get_trees():
    final_output = expand("results/{window}.raxml.bestTree",window=WINDOWS)
    return final_output


fileType = ["markedTree.pdf","pca.pdf"]
def get_R_outs():
    final_output = expand("results/{window}.{figure}",window=WINDOWS,figure=fileType)
    final_output.append("results/metrics.tsv")
    return final_output

# def get_R_outs():
#    final_output = expand("results/{window}.markedTree.pdf",window=WINDOWS)
#    final_output.extend(expand("results/{window}.pca.pdf",window=WINDOWS))
#    final_output.append("results/metrics.tsv")
#    return final_output