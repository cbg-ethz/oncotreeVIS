import json
from anytree.exporter import JsonExporter
from anytree import PreOrderIter, RenderTree
#from utils_anytree import *
import sys

def createOncotreeVISInput(anytrees, metadata, dataset_name):
  tree_map = {}
  exporter = JsonExporter(indent=2, sort_keys=False)
  for sample_name in anytrees:
    tree = anytrees[sample_name]
    for node in PreOrderIter(tree):
      if node.parent: #not root
        gene_events = {}
        if "tupro" in dataset_name:
          if hasattr(node, 'region_cn_state'):
            del node.region_cn_state
          del node.gene_state
          del node.node_label
          del node.node_depth
          del node.num_cells
          for gene in node.gene_cn_events:
            gene_events[gene] = {}
            gene_events[gene]["CNA"] = node.gene_cn_events[gene]
          del node.gene_cn_events
          node.matching_label = 0

        elif "morita" in dataset_name:
          gene = node.label[0]
          gene_events[gene] = {}
          variant = ""
          if hasattr(node, 'variant'): 
            variant = node.variant
          gene_events[gene]["SNV"] = variant
          del node.label
          del node.variant
  
        elif "compass" in dataset_name: 
          for affected_gene in node.label:
            gene_split = affected_gene.split("+")
            gene = gene_split[0]
            if gene not in gene_events:
              gene_events[gene] = {}
            if len(gene_split) == 1:
              gene_events[gene]["SNV"] = ""
            else:
              if gene_split[1] == "del":
                gene_events[gene]["CNA"] = "-"
              elif gene_split[1] == "gain":
                gene_events[gene]["CNA"] = "+"
              elif gene_split[1] == "loh":
                gene_events[gene]["CNA"] = "loh"
              else:
                gene_events[gene]["SNV"] = gene_split[1]
          del node.label
        del node.name
        if gene_events:
          node.gene_events = gene_events
        
    json_tree = json.loads(exporter.export(tree))
    tree_map[sample_name] = {}
    tree_map[sample_name]["tree"] = json_tree
    sample_basename = sample_name.split("_")[0]
    if "morita" in dataset_name or "compass" in dataset_name:
      sample_basename = sample_name.rsplit("-", 1)[0]
    if "tracerx" in dataset_name:
      sample_basename = sample_name.rsplit("_", 1)[0]
    if sample_basename in metadata:
      tree_map[sample_name]["metadata"] = metadata[sample_basename]
    else:
      tree_map[sample_name]["metadata"] = {}
  return tree_map
