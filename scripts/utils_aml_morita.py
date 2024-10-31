import os
import re
import sys
import json

from anytree import Node, RenderTree

def getTreesAMLMorita(dir):
  anytrees = {}
  gene_ids = {}
  for file in os.listdir(dir):
    sample_name = re.match('.*all_(.+)\.gv', file).groups()[0]
    with open(dir + "/" + file) as f:
      lines = list(f)
      anytree_nodes = {}
      root = None
      for line in lines:
        # Match the node label lines.
        match = re.match('^(\d+)\[.+label="(.+)",', line)
        if match:
          match = match.groups()
          assert len(match) == 2
          node_id = match[0]
          parse_mutation = match[1].split("_")
          gene = parse_mutation[0]
          if gene not in gene_ids:          
            gene_ids[gene] = len(gene_ids) + 1

          variant = ""
          if len(parse_mutation) > 1:
            variant = parse_mutation[1]
          anytree_nodes[node_id] = Node(node_id, node_id=node_id, label=[gene], variant=variant)
          if gene == "Root":
            root = anytree_nodes[node_id]
          continue

        # Match the edge lines.
        match = re.match('^(\d+) -> (\d+)', line)
        if match:
          match = match.groups()
          assert len(match) == 2
          parent_id = match[0]
          node_id = match[1]
          anytree_nodes[node_id].parent = anytree_nodes[parent_id]

        # Match the subclone sizes.
        match = re.match('^s_(\d+)\[label=\"(.+)%', line)
        if match:
          match = match.groups()
          assert len(match) == 2
          node_id = match[0]
          size_percent = round(float(match[1]) / 100, 4)
          anytree_nodes[node_id].size_percent = size_percent

    for node_id,node in anytree_nodes.items():
      node.matching_label = gene_ids[node.label[0]]
       
    anytrees[sample_name] = root
  return anytrees

def getTreesAMLCompass(dir_path):
  anytrees = {}
  for filename in os.listdir(dir_path):
    sample_name = filename.split('.')[0].split("_")[0]
    with open(os.path.join(dir_path, filename)) as json_file:
      nodes = json.load(json_file)["nodes"]

    node_id_object_map = {}
    # Add empty root.
    root = Node(name=0, node_id=0, label="Root")
    node_id_object_map[0] = root
    for node in nodes:
      mutations = []
      for mutation in node["SNV"]:
        snv_split = mutation.split(" ")
        if snv_split[0] == "SNP":
          snv_split = snv_split[1:]
        if len(snv_split) == 1:
          mutations.append(snv_split[0])
        else:
          mutations.append(snv_split[0] + "+" + snv_split[1])
      cnas = []
      for cn in node["CNA"]:
        cn_split = cn.split(" ")
        if cn_split[0] == "Gain":
          cnas.append(cn_split[1] + "+gain")
        elif cn_split[0] == "Loss":
          cnas.append(cn_split[1] + "+del")
        elif cn_split[0] == "CNLOH":
          cnas.append(cn_split[1] + "+loh")
      label = list(set(mutations + cnas))

      node_id = int(node["name"].split(" ")[-1]) + 1 # Reserve 0 for the root.
      if node["parent"] ==  '-':
        parent_id = 0
      else:
        parent_id = int(node["parent"].split(" ")[-1]) + 1

      anytree_node = Node(name=node_id,
          node_id=node_id,
          parent_id=parent_id,
          label=label)
      node_id_object_map[node_id] = anytree_node

    # Remove nodes with empty labels. 
    for node_id, node in node_id_object_map.items():
      if not node.label: # empty node -> remove empty nodes
        node_id_object_map[node_id] = node_id_object_map[node.parent_id] # replace the node by its parent

    # Populate parent and node_depth.
    for node_id, node in node_id_object_map.items():
      if node.node_id == 0: # if root
        continue
      node.parent = node_id_object_map[node.parent_id]
      del node.parent_id
    anytrees[sample_name] = root

  return anytrees
