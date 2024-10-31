import numpy as np
import pandas as pd
import rpy2.robjects.numpy2ri
from rpy2.robjects.packages import importr
rpy2.robjects.numpy2ri.activate()
import rpy2.robjects as robjects

import anytree
from anytree import Node
from anytree import RenderTree, PreOrderIter
from anytree.exporter import JsonExporter
import json
import sys

def create_anytree_node(id, parent=None):
  node = Node(id, node_id=id, matching_label=id)
  if parent:
    node.parent = parent
  return node

def tree_matrix_to_anytree(matrix, branch_ends_list):
  anytree_node_map = {}
  num_edges = int(matrix.nrow)
  for idx in range(num_edges):
    parent_id = matrix[idx]
    child_id = matrix[idx + num_edges]

    if parent_id in anytree_node_map:
      parent = anytree_node_map[parent_id]
    else:
      parent = create_anytree_node(parent_id)
      anytree_node_map[parent_id] = parent

    if child_id in anytree_node_map:
      child = anytree_node_map[child_id]
    else:
      child = create_anytree_node(child_id, parent)
      if child_id in branch_ends_list:
        child.branch_end = True
      anytree_node_map[child_id] = child
    child.parent = parent

  for id, node in anytree_node_map.items():
    if not node.parent: # root
      return node

def getConipherTrees(rds_file, metadata_file):
  readRDS = robjects.r['readRDS']
  rds = readRDS(rds_file)
  exporter = JsonExporter(indent=2, sort_keys=False)

  #metadata = readRDS(metadata_file)
  #print(type(metadata))

  anytrees = {}
  clusters = []
  for sample_name, data in dict(zip(rds.names, rds)).items():
    tree_data = dict(zip(data.names, data))["graph_pyclone"]
    conserved_branch_ends = [pair[1].split(":")[1] for pair in tree_data.rx('consensus_branches')[0].items()]
    cluster = []

    original_tree = tree_data.rx('Corrected_tree')
    original_tree_matrix = dict(zip(original_tree.names, original_tree))["Corrected_tree"]
    anytree = tree_matrix_to_anytree(original_tree_matrix, conserved_branch_ends)
    original_tree_key = sample_name + "_0"
    anytrees[original_tree_key] = anytree
    cluster.append(original_tree_key)

    alt_trees = tree_data.rx('alt_trees')
    tree_list = dict(zip(alt_trees.names, alt_trees))["alt_trees"]
    for idx, tree_matrix in enumerate(tree_list):
      anytree = tree_matrix_to_anytree(tree_matrix, conserved_branch_ends)
      key = sample_name + "_" + str(idx+1)
      cluster.append(key)
      anytrees[key] = anytree
      if idx > 100:
        break
    clusters.append(cluster)

  return anytrees, clusters

