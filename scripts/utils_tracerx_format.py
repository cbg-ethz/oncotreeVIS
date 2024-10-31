import re
from anytree import Node, RenderTree

ROOT_LABEL = "GL"

def createNode(mutation_id, gene_string):
  node = Node(mutation_id, node_id=mutation_id)
  if mutation_id: # not root
    node.gene_events = getGeneEvents(gene_string)
    node.matching_label = mutation_id 
  return node

def getMutationId(gene_string, mutation_ids):
  if gene_string == ROOT_LABEL:
    return 0
  if gene_string not in mutation_ids:
    mutation_ids[gene_string] = len(mutation_ids) + 1
  return mutation_ids[gene_string]

def getGeneEvents(gene_string):
  gene_set = gene_string.split(";")
  gene_events = {}
  for gene in gene_set:
    gene_events[gene] = {"SNV":""}
  return gene_events

def getTreesInTracerxFormat(path):
  with open(path) as f:
    tree_file = list(f)

  num_patients = int(re.match('(\d+) #patients', tree_file[0]).groups()[0])
  line_idx = 1
  patient_idx = 1
  mutation_ids = {}
  anytrees = {}
  while True:
    match = re.match('(\d+) #trees for.* (.+)\n', tree_file[line_idx]).groups()
    num_trees = int(match[0])
    patient_id = match[1]
    line_idx += 1
    for tree_idx in range(1, num_trees+1):
      match = re.match('(\d+) #edges', tree_file[line_idx]).groups()
      num_edges = int(match[0])
      line_idx += 1 
      tree_nodes = {}
      for edge_idx in range(0, num_edges):
        match = re.match('([^ ]+) ([^ ]+)\n', tree_file[line_idx]).groups()
        gene_set_1 = match[0]
        gene_set_2 = match[1]
        mutation_id_1 = getMutationId(gene_set_1, mutation_ids)
        mutation_id_2 = getMutationId(gene_set_2, mutation_ids)
        if mutation_id_1 not in tree_nodes:
           tree_nodes[mutation_id_1] = createNode(mutation_id_1, gene_set_1)
        if mutation_id_2 not in tree_nodes:
           tree_nodes[mutation_id_2] = createNode(mutation_id_2, gene_set_2)
        tree_nodes[mutation_id_2].parent = tree_nodes[mutation_id_1]
        line_idx += 1
      anytrees[patient_id + "_" + str(tree_idx)] = tree_nodes[0]
    patient_idx += 1
    if patient_idx == num_patients+1:
      break

  return anytrees


