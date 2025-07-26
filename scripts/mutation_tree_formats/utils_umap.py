import pandas as pd
import matplotlib.pyplot as plt
import umap 
#import umap.plot
from sklearn.preprocessing import StandardScaler
import plotly.express as px
from sklearn.manifold import TSNE
from sklearn.manifold import MDS
import sys

def savePlot(projections, color_idx, out_file):
  fig = plt.scatter(
    projections[:, 0],
    projections[:, 1],
    c=color_idx
  )
  plt.savefig(out_file)
  plt.close()

def computeUMAP(df_embeddings, df_distances, clusters, out_file_prefix):

  samples_to_remove = df_embeddings.filter(regex='#2$', axis=0).index
  df_embeddings = df_embeddings.drop(index=samples_to_remove)

  color_idx = None
  if clusters:
    sample_color_map = {}
    for idx, x in enumerate(clusters):
      for sample in clusters[idx]:
        sample_color_map[sample] = idx

    color_idx = []
    for sample in df_embeddings.index:
      color_idx.append(sample_color_map[sample])

  embeddings = df_embeddings.to_numpy()
  reducer = umap.UMAP()
  projections = reducer.fit_transform(embeddings)
  savePlot(projections, color_idx, out_file_prefix + "_umap_embeddings.png")

  distances = df_distances.to_numpy()
  reducer = umap.UMAP(metric='precomputed')
  projections = reducer.fit_transform(distances)
  savePlot(projections, color_idx, out_file_prefix + "_umap_distances.png")

  tsne = TSNE(n_components=2, verbose=1, perplexity=40, n_iter=300)
  projections = tsne.fit_transform(embeddings)
  savePlot(projections, color_idx, out_file_prefix + "_tsne.png")

  mds = MDS(random_state=0)
  points = mds.fit_transform(distances)  
  savePlot(points, color_idx, out_file_prefix + "_mds.png")

  '''
  umap_colors_map = {} 
    for idx, cluster in enumerate(tree_clusters):  
      if len(cluster) == 1:
        umap_colors_map[cluster[0]] = len(tree_clusters)
      else:
        for sample in cluster:
          umap_colors_map[sample] = idx
  umap_colors = [str(umap_colors_map[sample]) for sample in df_embeddings.index]
  
  fig = px.scatter(
      projections, x=0, y=1
      #color=umap_colors, labels={'color': 'tree cluster'}, hover_data={"sample": df_embeddings.index},
      #opacity=0.75
  )
  '''
