import pandas as pd
import os

from anytree import Node, RenderTree

def add_to_node_map(node_map, node_id):
  if node_id not in node_map:
    node_map[node_id] = Node(node_id, node_id=node_id)

def is_float(string):
    if str(string).replace(".", "").isnumeric():
        return True
    else:
        return False

def getTreesNoble2022 (dir):
  anytrees = {}
  for filename in os.listdir(dir):
    sample_name = filename.split('.')[0]
    tree_df = pd.read_csv(os.path.join(dir, filename), header=0)

    anytree_nodes = {}
    node_id_map = {}

    # Root has id 0, so add it in the first place.
    node_id_map[0] = 0
    add_to_node_map(anytree_nodes, 0)
    root = anytree_nodes[0]

    for index, row in tree_df.iterrows():
      node_id = row['Identity']
      if is_float(node_id):
        node_id = int(node_id)
      if node_id not in node_id_map:
        node_id_map[node_id] = len(node_id_map)
      add_to_node_map(anytree_nodes, node_id_map[node_id])#, str(row['Population']))

    for index, row in tree_df.iterrows():
      node_id = row['Identity']
      if is_float(node_id):
        node_id = int(node_id)
      if node_id != 0: # Not the root.
        parent_id = row['Parent']
        if is_float(parent_id):
          parent_id = int(parent_id)
        anytree_nodes[node_id_map[node_id]].parent = anytree_nodes[node_id_map[parent_id]]
    anytrees[sample_name] = root
  return anytrees     

