# Usage: python main.py aml_morita
 
import argparse
import os
import time
import pandas as pd
import copy
import json
import sys
from scipy.spatial import distance
from anytree import RenderTree, PreOrderIter
from anytree.importer import JsonImporter

from utils import *
from utils_aml_morita import getTreesAMLMorita, getTreesAMLCompass
from utils_noble2022 import getTreesNoble2022
from utils_scatrex import getScatrexTrees
from utils_tracerx_format import getTreesInTracerxFormat
from utils_conipher import getConipherTrees
from utils_umap import computeUMAP
from utils_oncotreevis import createOncotreeVISInput

parser = argparse.ArgumentParser()
parser.add_argument("cancer_type", choices=[
    "tupro_aml", 
    "aml_morita", "aml_compass", "morita_side_by_side", 
    "noble_2022", 
    "brca_razavi", 
    "tracerx_lung", 
    "tracerx421"])
parser.add_argument("--metadata", default='', type=str)
parser.add_argument("--embeddings", default='', type=str)
parser.add_argument("--distance_threshold", type=float, default='1', required=False) # distance 0 means perfect match
parser.add_argument("--out", default='out', type=str)
args = parser.parse_args()

## Read json trees.
def read_json_trees(path, samples_to_remove=[]):
  with open(path) as json_file:
    sample_map = json.load(json_file)
    for key in samples_to_remove:
      sample_map.pop(key, None)
    return sample_map

def isNaN(num):
  return num != num

_CANCER_TYPE = args.cancer_type

timestamp = str(int(time.time()))
filaname = "_".join([timestamp, _CANCER_TYPE + ".json"])
create_dir(args.out)
output_path = os.path.join(args.out, filaname)

metadata = {}
if os.path.isfile(args.metadata):
    metadata = pd.read_csv(args.metadata).to_dict("index")
    keys = copy.deepcopy(list(metadata.keys()))
    for key in keys:
      sample_name_key = ""
      if "sample_name" in metadata[key]:
        sample_name_key = "sample_name"
      elif "Patient_ID" in metadata[key]:
        sample_name_key = "Patient_ID"
      else:
        break
      sample_name = metadata[key][sample_name_key]
      metadata[sample_name] = metadata[key]
      for kkey in metadata[key]:
        if isNaN(metadata[key][kkey]):
          metadata[sample_name][kkey] = ""
      del metadata[sample_name][sample_name_key]
      del metadata[key]

df_distances = None
if os.path.isfile(args.embeddings):
  df_embeddings = pd.read_csv(args.embeddings, index_col=0)
  samples_to_remove = df_embeddings.filter(regex='#2$', axis=0).index
  df_embeddings = df_embeddings.drop(index=samples_to_remove)
  samples = df_embeddings.index
  df_distances = pd.DataFrame(0, columns=samples, index=samples).astype(float)
  for sample_1 in samples:
    if "#" in sample_1:
      continue
    for sample_2 in samples:  
      if "#" in sample_2:
        continue
      embedding_1 = list(df_embeddings.loc[sample_1, :])
      embedding_2 = list(df_embeddings.loc[sample_2, :])
      df_distances.loc[sample_1, sample_2] = distance.cosine(embedding_1, embedding_2)

importer = JsonImporter()
anytrees = {}
clusters = None
if _CANCER_TYPE == "tupro_aml":
  anytrees = read_json_trees("data/tupro/samples_aml_v1.15_prioriry_genes.js")
  for sample_name in anytrees:
    anytrees[sample_name] = importer.import_(json.dumps(anytrees[sample_name]["event_tree"]))

elif _CANCER_TYPE == "aml_morita":
  raw_data_dir = "data/aml_morita/trees" 
  anytrees = getTreesAMLMorita(raw_data_dir)
  clusters = [['AML-07-002', 'AML-18-002', 'AML-48-001', 'AML-104-001', 'AML-56-001', 'AML-64-001', 'AML-115-001', 'AML-08-001', 'AML-60-001', 'AML-70-001', 'AML-108-001', 'AML-109-001', 'AML-13-001', 'AML-40-001', 'AML-52-001', 'AML-63-005'], ['AML-79-001', 'AML-101-001', 'AML-110-001', 'AML-106-001', 'AML-78-001', 'AML-107-002', 'AML-117-001'], ['AML-91-001', 'AML-95-001', 'AML-93-001', 'AML-94-001', 'AML-112-001', 'AML-116-001', 'AML-118-001', 'AML-85-001', 'AML-80-001', 'AML-90-001', 'AML-103-001'], ['AML-47-001', 'AML-43-001', 'AML-23-001', 'AML-21-002', 'AML-10-001'], ['AML-87-001', 'AML-54-001', 'AML-03-001', 'AML-19-001', 'AML-01-002'], ['AML-28-001', 'AML-62-001', 'AML-15-001', 'AML-59-001', 'AML-16-001'],  ['AML-55-001', 'AML-33-001', 'AML-57-001', 'AML-11-001'],  ['AML-77-001'], ['AML-67-001'], ['AML-29-001'], ['AML-84-001'], ['AML-66-003'], ['AML-45-001'], ['AML-92-001'], ['AML-44-001'], ['AML-105-001'], ['AML-25-001'], ['AML-42-001'], ['AML-98-001'], ['AML-111-001'], ['AML-89-001'], ['AML-86-001'], ['AML-120-001', 'AML-88-002', 'AML-39-002', 'AML-73-001', 'AML-122-001', 'AML-119-001'],  ['AML-38-003', 'AML-61-001', 'AML-41-001', 'AML-74-001', 'AML-14-001'], ['AML-72-001', 'AML-46-001'], ['AML-50-001', 'AML-102-001'], ['AML-09-002', 'AML-12-001', 'AML-113-001', 'AML-75-001'], ['AML-24-001', 'AML-123-001', 'AML-02-001', 'AML-32-001', 'AML-82-001', 'AML-121-001', 'AML-26-001', 'AML-114-001', 'AML-30-001', 'AML-20-001', 'AML-99-005', 'AML-58-001', 'AML-83-002', 'AML-97-006', 'AML-37-001', 'AML-65-001', 'AML-27-001', 'AML-68-001'],  ['AML-34-001', 'AML-22-001', 'AML-06-001', 'AML-71-001', 'AML-69-001', 'AML-17-001', 'AML-49-001', 'AML-35-001', 'AML-53-001'], ['AML-36-001', 'AML-05-001', 'AML-81-001', 'AML-96-001'],  ['AML-04-003'], ['AML-31-001'], ['AML-100-001'], ['AML-76-001', 'AML-51-001']]

elif _CANCER_TYPE == "aml_compass":
  raw_data_dir = "data/aml_compass/trees"
  anytrees = getTreesAMLCompass(raw_data_dir)

elif _CANCER_TYPE == "morita_side_by_side":
  anytrees_morita = getTreesAMLMorita("data/aml_morita/trees")
  anytrees_compass = getTreesAMLCompass("data/aml_compass/trees")
  keys = set(anytrees_morita.keys()).union(set(anytrees_compass.keys()))
  anytrees = {}
  for key in keys:
    if key in anytrees_morita:
      anytrees[key + "_morita"] = anytrees_morita[key]
    if key in anytrees_compass:
      anytrees[key + "_compass"] = anytrees_compass[key]

elif _CANCER_TYPE == "noble_2022":
  raw_data_dir = "data/noble_2022/trees"  
  anytrees = getTreesNoble2022(raw_data_dir)
  clusters = [["AML-16-001", "AML-33-001", "AML-35-001", "AML-73-001", "AML-05-001", "MED012", "MED001", "AML-77-001", "MED034", "AML-02-001", "AML-55-001"], ["UMM069", "UMM059", "UMM063", "UMM061", "UMM064", "UMM065", "UMM062", "CRUK0029", "UMM066", "CRUK0065", "PD9852", "PD9849"], ["TN5", "K255", "TN7", "TN2", "TN1", "TN4", "TN8"], ["TN6", "TN3", "K136", "K252", "K153", "MED024", "MED023", "K448", "PD9694", "CRUK0096", "MED027", "CRUK0062", "CRUK0071"]]

elif _CANCER_TYPE == "brca_razavi":
  raw_data_dir = "data/brca_razavi/trees_brca_razavi.txt"
  anytrees = getTreesInTracerxFormat(raw_data_dir)

  selected_genes = ["PIK3CA", "NF1", "ESR1", "CDH1", "GATA3", "TP53", "PTEN", "MAP3K1", "KMT2D", "KMT2C", "FOXA1",
       "RB1", "EPHA7", "TSC2", "RHOA", "PIK3R1", "PRDM1",  "PBRM1", "CD79A"]

  def getNodeGene(node):
    if node.parent:
      return list(node.gene_events.keys())[0]
    else:
      return "Root"

  unique_trees = {}
  for key, tree in anytrees.items():
    include_tree = False
    nodes_to_delete = []
    initial_num_nodes = len(list(PreOrderIter(tree)))
    for node in PreOrderIter(tree):
      if node.parent:
        assert len(node.gene_events.keys()) == 1
        gene = getNodeGene(node)
        if gene not in selected_genes:
          nodes_to_delete.append(node)
    for node in nodes_to_delete:
      if node.parent:
        parent = node.parent
        for child in node.children:
          child.parent = parent
        node.parent = None

    if len(list(PreOrderIter(tree))) == 1:
      continue

    links = []
    for node in PreOrderIter(tree):
      if node.parent:
        links.append(getNodeGene(node.parent) + "-" + getNodeGene(node))
    tree_string = ";".join(links)
    patient_id = key.split("_")[0]
    unique_trees[patient_id] = {}
    unique_trees[patient_id][tree_string] = {"key":key, "tree":tree}

  anytrees = {}
  for patient_id, tree in unique_trees.items():
    for tree_string, data in tree.items():
      anytrees[data["key"]] = data["tree"]

elif _CANCER_TYPE == "tracerx_lung":
  raw_data_dir = "data/tracerx_lung/tracerx_raw_lung_trees.txt"  
  anytrees = getTreesInTracerxFormat(raw_data_dir)
  clusters = [["CRUK0001_7", "CRUK0001_6", "CRUK0001_2", "CRUK0001_10", "CRUK0001_5", "CRUK0001_9", "CRUK0001_11", "CRUK0001_8", "CRUK0001_1", "CRUK0001_3", "CRUK0001_4"], ["CRUK0011_1", "CRUK0038_1", "CRUK0042_1", "CRUK0044_1", "CRUK0059_1", "CRUK0069_1", "CRUK0017_1"], ["CRUK0054_1", "CRUK0012_1", "CRUK0019_1"], ["CRUK0013_1", "CRUK0013_2", "CRUK0061_1", "CRUK0013_3", "CRUK0013_4"], ["CRUK0016_7", "CRUK0016_4", "CRUK0016_6", "CRUK0016_14", "CRUK0016_11", "CRUK0016_13", "CRUK0016_10", "CRUK0016_12", "CRUK0016_8", "CRUK0016_9", "CRUK0016_2", "CRUK0016_5", "CRUK0016_1", "CRUK0016_3"], ["CRUK0051_1", "CRUK0063_1", "CRUK0063_4", "CRUK0063_2", "CRUK0063_5", "CRUK0063_3", "CRUK0063_6"], ["CRUK0097_1", "CRUK0098_1", "CRUK0074_1", "CRUK0071_1", "CRUK0052_1"], ["CRUK0014_1", "CRUK0027_1"], ["CRUK0073_1", "CRUK0062_1", "CRUK0062_2"], ["CRUK0022_1", "CRUK0004_1", "CRUK0058_1"], ["CRUK0070_1", "CRUK0067_1", "CRUK0075_1"], ["CRUK0031_1", "CRUK0002_1", "CRUK0002_2"], ["CRUK0034_1"], ["CRUK0068_1"], ["CRUK0009_1"], ["CRUK0087_1"], ["CRUK0085_1"], ["CRUK0036_1"], ["CRUK0081_1"], ["CRUK0055_1"], ["CRUK0065_1"], ["CRUK0005_1"], ["CRUK0078_1"], ["CRUK0099_1"], ["CRUK0035_1"], ["CRUK0032_1"], ["CRUK0082_1"], ["CRUK0010_1"], ["CRUK0093_1"], ["CRUK0095_1"], ["CRUK0077_1"], ["CRUK0046_1"], ["CRUK0049_1"], ["CRUK0096_1"], ["CRUK0007_1"], ["CRUK0072_1"], ["CRUK0043_1"], ["CRUK0089_1"], ["CRUK0088_1"], ["CRUK0039_1"], ["CRUK0086_1"], ["CRUK0048_1"], ["CRUK0033_1"], ["CRUK0047_1"], ["CRUK0030_1"], ["CRUK0037_1"], ["CRUK0008_1"], ["CRUK0018_1"], ["CRUK0041_1"], ["CRUK0080_1"], ["CRUK0084_1"], ["CRUK0028_1"], ["CRUK0025_1"], ["CRUK0060_1"], ["CRUK0091_1"], ["CRUK0050_1"], ["CRUK0083_1"], ["CRUK0040_1"], ["CRUK0092_1"], ["CRUK0021_1"], ["CRUK0094_1"], ["CRUK0079_1"], ["CRUK0090_1"], ["CRUK0029_1"], ["CRUK0057_1"], ["CRUK0003_1"], ["CRUK0100_1"], ["CRUK0026_1"], ["CRUK0045_1", "CRUK0020_1", "CRUK0020_2"], ["CRUK0024_1", "CRUK0024_2"], ["CRUK0066_1", "CRUK0076_1"], ["CRUK0006_2", "CRUK0006_1", "CRUK0006_3"], ["CRUK0023_1", "CRUK0023_2", "CRUK0015_1", "CRUK0064_1", "CRUK0056_1"]]

elif _CANCER_TYPE == "tracerx421":
  rds_file = "data/tracerx421_conipher/conipher_tracerx421.RDS"
  metadata = "data/tracerx421_conipher/20221109_TRACERx421_all_patient_df.rds"
  anytrees, clusters = getConipherTrees(rds_file, metadata)

# OncotreeVIS.  
trees_oncotreevis = createOncotreeVISInput(anytrees, metadata, _CANCER_TYPE)
json_oncotreevis = {}
json_oncotreevis["trees"] = trees_oncotreevis
if clusters:
  json_oncotreevis["clusters"] = clusters

if df_distances is not None:
  pairwise_distances = []
  for sample_1 in df_distances.index:
    for sample_2 in df_distances.index:
      pairwise_distances.append({"sample_1":sample_1, "sample_2":sample_2, "distance":df_distances[sample_1][sample_2]})
  json_oncotreevis["pairwise_tree_distances"] = pairwise_distances

  umap_filaname_prefix = "_".join([timestamp, _CANCER_TYPE])
  umap_path_prefix = os.path.join(args.out, umap_filaname_prefix)
  computeUMAP(df_embeddings, df_distances, clusters, umap_path_prefix)

file = open(output_path, "w")
file.write(json.dumps(json_oncotreevis))
file.close()
