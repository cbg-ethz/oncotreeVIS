import copy
import json
import pandas as pd
import sys

from anytree.importer import JsonImporter
from anytree.exporter import JsonExporter
from anytree import Node, RenderTree, PreOrderIter

def getScatrexTrees(dir_path, metadata):
  anytrees = {}
  root_label = "Root"
  gene_expressions = []  
  gene_expression_threshold = 0

  for filename in os.listdir(dir_path):
    sample_name = filename.split('.')[0].split("-")[0]
    with open(os.path.join(dir_path, filename)) as json_file:
      nodes = json.load(json_file)

    cn_neutral_state = 2
    if metadata[sample_name]["has_wgd"]:
      cn_neutral_state = 4
    genes = nodes["genes"]

    anytree_nodes = {}
    # Add empty root.
    root = Node(root_label, node_id=root_label)
    anytree_nodes[root_label] = root

    for node_id in nodes:
      if node_id == "genes": 
        continue

      node = nodes[node_id]
      if node["parent"] == "NULL":
        parent_id = root_label
      else:
        parent_id = node["parent"]

      cnv_values = dict(zip(genes, node["cnv"]))
      exp_values = dict(zip(genes, node["total_exp"]))
      #cnvs = [gene + "_amp" for gene in cnv_values if cnv_values[gene] > cn_neutral_state]
      #cnvs.extend([gene + "_del" for gene in cnv_values if cnv_values[gene] < cn_neutral_state])
      #gene_expressions.extend(list(exp_values.values()))
      #gene_exp = [gene + "_up" for gene in exp_values if exp_values[gene] > 2]
      #gene_exp.extend([gene + "_down" for gene in exp_values if exp_values[gene] < -2])

      anytree_node = Node(node_id,
          node_id=node_id,
          parent_id=parent_id,
          gene_events=""
      )
      anytree_nodes[node_id] = anytree_node

    # Populate parent.
    for node_id, node in anytree_nodes.items():
      if node_id == root_label: # if root
        continue
      node.parent = anytree_nodes[node.parent_id]
      del node.parent_id

  '''
  # Print distribution of gene differential expressions.
  heatmap_dir = os.path.join(out_dir, "heatmaps")
  if not os.path.exists(heatmap_dir):
    os.makedirs(heatmap_dir)

  plt.figure()
  sns.set(font_scale=1)
  histogram = sns.histplot(data = gene_expressions)
  plt.yscale('log')
  histogram.set(ylabel='Log counts', xlabel='Expression levels')
  histogram_fig = histogram.get_figure()
  histogram_filename = os.path.join(heatmap_dir, "histogram_expression_levels_" + cancer_type + ".png")
  histogram_fig.savefig(histogram_filename, format='png', dpi=300)
  '''
 
  return anytrees_event
