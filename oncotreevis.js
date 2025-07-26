///////////////////////
//// Main function ////
///////////////////////
function async_func(input, callback) {
  setTimeout(function () {
    callback(input);
    $("body").removeClass("wait");
  }, 0);
  $("body").addClass("wait");
}     
    
function oncotreeVIS(data, container_div_id) {
  var args = {
    "data": data,
    "container_div_id": container_div_id,
  }
  async_func(args, oncotreeVIS_slow)
}

function oncotreeVIS_slow(args) {
  data = args.data 
  container_div_id = args.container_div_id

  var num_max_node_colors = 10
  var tree_cohort_div_id = "tree_cohort"
  var tree_info_div_id = "tree_info"

  // Compute additional information from the data.
  trees = data["trees"]

  // Convert future key values (`matching_label`) to strings. 
  for (var sample_name in trees) {
    var tree_json = trees[sample_name]["tree"]
    var node_list = getTreeNodes(tree_json)
    for (var node of node_list) {  
      if (!node.parent && !node.data.matching_label) { // root
        node.matching_label = -1
      } 
      if (node.data.matching_label) {
        node.data.matching_label = node.data.matching_label.toString()
      }
    }
    // Propagate gene states.
    populateGeneStates(tree_json)
  }

  // Add data to all trees.
  matching_label_tree_map = new Map()
  unique_matching_labels_cohort = {} // map {matching_label: bool}
  for (sample_name in trees) {
    tree = trees[sample_name]["tree"]
    matching_labels = new Set(Array.from(getTreeMalignantMatchingLabels(tree).keys()))
    for (let matching_label of matching_labels) {
      key = matching_label
      if (key in unique_matching_labels_cohort) {
        unique_matching_labels_cohort[key] = false
      } else {
        unique_matching_labels_cohort[key] = true
      }
      // Populate aux variable.
      if (matching_label_tree_map.has(key)) {
        matching_label_tree_map.get(key).push(sample_name)
      } else {
        matching_label_tree_map.set(key, [sample_name])
      }
    }
  }

  display_matching_labels = false
  num_non_unique_matching_labels = Object.values(unique_matching_labels_cohort).filter(x => x == false).length
  if(num_non_unique_matching_labels > num_max_node_colors) {
    display_matching_labels = true
  }
  cohort_node_color_map = {}
  cohort_node_color_generator = new ColorGenerator(node_colors, shuffle=true, transparent_flavour=0.3)
  for (key in unique_matching_labels_cohort) {
    if (!unique_matching_labels_cohort[key]) {
      cohort_node_color_map[key] = cohort_node_color_generator.next()
    } 
  }

  // Populate the node colors.
  for (sample_name in trees) {
    tree = trees[sample_name]["tree"]
    var node_list = getTreeNodes(trees[sample_name]["tree"])
    for(node of node_list) {
      if (node.parent && node.data.matching_label) { // not root.
        key = node.data.matching_label
        if (key in cohort_node_color_map) {
          node.data.color_in_cohort = cohort_node_color_map[key]
          node.data.display_node_label_cohort = display_matching_labels
        }
      }
    }
  } 

  // Clusters.
  if ("clusters" in data && data["clusters"].length != 0) {
    clusters = data["clusters"]
  } else {
    clusters = [Object.keys(trees)] // One big cluster.
  }

  // Colors at cluster level.
  used_colors = {}
  metadata_color_map = {}
  cluster_color_generator = new ColorGenerator(background_colors, shuffle=false, transparent_flavour=0.7)
  cluster_node_color_generator = new ColorGenerator(node_colors, shuffle=true, transparent_flavour=0.3)
  for (const [i, cluster] of clusters.entries()) {
    conserved_branches_cluster = getBranches(trees[cluster[0]]["tree"])
    unique_matching_labels_cluster = {}
    for (sample_name of cluster) {
      // Cluster colors.
      trees[sample_name]["cluster_color"] = "white"
      trees[sample_name]["2d_color"] = "white"
     
      // Nodes and branches.
      tree = trees[sample_name]["tree"]
      conserved_branches_cluster = conserved_branches_cluster.intersection(getBranches(tree))           
      matching_labels = new Set(Array.from(getTreeMalignantMatchingLabels(tree).keys()))
      for (matching_label of matching_labels) {
        key = matching_label
        if (key in unique_matching_labels_cluster) {
          unique_matching_labels_cluster[key] = false
        } else {
          unique_matching_labels_cluster[key] = true
        }
      } 
    }

    display_matching_labels = false
    num_non_unique_matching_labels = Object.values(unique_matching_labels_cluster).filter(x => x == false).length
    if(num_non_unique_matching_labels > num_max_node_colors) {
      display_matching_labels = true
    }

    var intersecting_clones = new Map()
    if (cluster.length > 1 && clusters.length > 1) {
      cluster_color = cluster_color_generator.next()
      var cluster_node_color_map = new Map() // matching_label:color
      for (key in unique_matching_labels_cluster) {
        if (!unique_matching_labels_cluster[key]) {
          if (key in used_colors) {
            cluster_node_color_map.set(key, used_colors[key])
          } else {
            color = cluster_node_color_generator.next()
            cluster_node_color_map.set(key, color)
            used_colors[key] = color
          }
        } 
      }

      for (sample_name of cluster) {
        // Update cluster colors.
        trees[sample_name]["cluster_color"] = cluster_color
        trees[sample_name]["2d_color"] = tinycolor(cluster_color).darken(50).toHexString()

        // For each cluster tree populate node colors and conserved branches
        // and compute common cluster events.
        var node_list = getTreeNodes(trees[sample_name]["tree"])
        for (var node of node_list) {
          // Common cluster events
          if (node.data.matching_label && node.data.gene_events) {
            var matching_label_key = node.data.matching_label
            if (cluster_node_color_map.has(matching_label_key)) {
              var color = cluster_node_color_map.get(matching_label_key)
              if (!(intersecting_clones.has(color))) {
                var empty_array = new Array() 
                intersecting_clones.set(color, empty_array)
              }
              var arr = intersecting_clones.get(color)
              arr.push(node.data.gene_events)
              intersecting_clones.set(color, arr)
            }
          }
        }
        for (var node of node_list) {
          if(node.parent) { // not root.
            // Conserved branches.
            if (conserved_branches_cluster.has(node.parent.data.matching_label + "_" + node.data.matching_label)) {
              node.data.conserved_parent_node_edge_cluster = true
            } 
            // Matching nodes.
            if (node.parent && node.data.matching_label) {
              key = node.data.matching_label
              if (cluster_node_color_map.has(key)) { 
                node.data.color_in_cluster = cluster_node_color_map.get(key)
                node.data.display_node_label_cluster = display_matching_labels
              }
            } 
          }
        }
      }
      // Compute gene event intersection.
      for (let [color, list] of intersecting_clones.entries()) {
        intersection = computeGeneEventIntersection(list)
        intersecting_clones.set(color, intersection)
      }
    }

    if (cluster.length > 1) {
      // Cluster metadata.
      var sample_metadata_map = {}
      for (sample_name of cluster) {
        if ("metadata" in trees[sample_name]) {
          sample_metadata = trees[sample_name]["metadata"]
          sample_metadata_map[sample_name] = sample_metadata
        }
      }
      [sample_metadata_colors, metadata_color_map] = getMetadataColorMap(sample_metadata_map, metadata_color_map)
      table_color_codes = getColorCodesTable(metadata_color_map)
      trees[cluster[0]]["sample_metadata_colors"] = sample_metadata_colors
      trees[cluster[0]]["table_color_codes"] = table_color_codes
      trees[cluster[0]]["matching_nodes_cluster_details"] = intersecting_clones
    }
  }

  // Metadata all samples.
  sample_list = []
  for (const [i, cluster] of clusters.entries()) {
    sample_list.push(...cluster)
  }
  var sample_metadata_map = {}
  for (sample_name of sample_list) {
    if ("metadata" in trees[sample_name]) {
      sample_metadata = trees[sample_name]["metadata"]
      sample_metadata_map[sample_name] = sample_metadata
    }
  }
  [sample_metadata_colors, metadata_color_map] = getMetadataColorMap(sample_metadata_map, metadata_color_map)
  table_color_codes = getColorCodesTable(metadata_color_map)
  trees[clusters[0][0]]["cohort_metadata"] = sample_metadata_colors
  trees[clusters[0][0]]["cohort_table_color_codes"] = table_color_codes

  // kNN
  knn_json = {}
  if ("pairwise_subclone_distances" in data) {
    // TODO
  }
  else if (mapSize(trees) < 1500 && mapSize(cohort_node_color_map) < 50) {
    for (let sample_1 in trees) {
      tree_1 = trees[sample_1]["tree"]
      //var matching_tree_ids = Object.keys(trees)
      var matching_labels = new Set(Array.from(getTreeMalignantMatchingLabels(tree_1).keys()))
      matching_tree_ids = new Set()
      for (let key of matching_labels) {
        for (let sample of matching_label_tree_map.get(key)) {
          matching_tree_ids.add(sample)
        }
      }
      for (let sample_2 of matching_tree_ids) {
        if (sample_1 != sample_2) {
          var tree_2 = trees[sample_2]["tree"]
          var matching_node_pairs = getMatchingNodesWithMatchingLabels(getTreeNodes(tree_1), getTreeNodes(tree_2))
          if (!matching_node_pairs[0]) {
            continue
          }
          const [nodes_1, links_1, max_depth_1] = getJSONNodes(sample_1, tree_1)
          const [nodes_2, links_2, max_depth_2] = getJSONNodes(sample_2, tree_2)
          const max_depth = Math.max(max_depth_1, max_depth_2)
          const max_similarity = Math.max(nodes_1.length, nodes_2.length)
          nodes = nodes_1.concat(nodes_2)
          links = links_1.concat(links_2)
          for (let pair of matching_node_pairs) {
            links.push({
              "source": sample_1 + "_" + pair[0].data.node_id,
              "target": sample_2 + "_" + pair[1].data.node_id,
              "similarity": 1,
            })
          }
          knn_json[sample_1 + "_" + sample_2] = {
            "sample_1": sample_1,
            "sample_2": sample_2,
            "nodes": nodes,
            "links": links,
            "similarity": matching_node_pairs.length / max_similarity,
            "max_depth": max_depth
          }
        } 
      }
    }
  }
  data["matching_trees"] = knn_json
  // END Populate additional data.

  var useful_variables = {
    "tree_cohort_div_id": tree_cohort_div_id, 
    "tree_info_div_id": tree_info_div_id,
    "data": data,
    "sorting": 0,
    "matching": 0,
    "rectangle": 0,
  }
  
  // DGIdb data.
  var gene_drug_map = parseGeneDrugInteractions(gene_drug_interaction)

  // Create the HTML elements. 
  addHTMLElements(container_div_id, useful_variables)
  
  // Populate with tree data. 
  populateTreeView_slow(useful_variables)
}

function getJSONNodes(sample_name, tree_json) {
  json_node_list = new Array()
  json_links = new Array()
  max_depth = 0
  for (let node of getTreeNodes(tree_json)) {
    is_root = false
    if (!node.parent) {
      is_root = true
    }
    json_node_list.push({
        "id": sample_name + "_" + node.data.node_id,
        "depth": node.depth,
        "matching_label": node.data.matching_label,
        "color": "white",
        "dx": node.x,
        "is_root": is_root,
        "is_neutral": node.data.is_neutral
    })
    if (node.parent) {
      json_links.push({
        "source": sample_name + "_" + node.data.node_id,
        "target": sample_name + "_" + node.parent.data.node_id,
        "similarity": -1,
      })
    }
    if (max_depth < node.depth) {
      max_depth = node.depth
    }
  }
  return [json_node_list, json_links, max_depth]
}

///////////////////////
//// HTML elements ////
///////////////////////
function createBlueBorder(){
  var div = document.createElement('span')
  div.style.border = "4px solid #c2dffe"
  div.style.borderRadius = "8px"
  div.style.paddingLeft = "0.5px"
  div.style.paddingTop = "5px"
  div.style.paddingBottom = "6px"
  return div
}

function createActionIcon(icon_class, id=null) {
  var button = document.createElement('span')
  if (id) {
    button.setAttribute("id", id)
  }
  button.innerHTML = '<a class="action_icon"><i class="' + icon_class +
      '" style="font-size:19px; cursor: pointer;"></i></a>&nbsp;'
  return button
}

function downloadDiv(original_element) {
  $("body").addClass("wait")

  const element = original_element.cloneNode(true)
  document.body.appendChild(element)

  var svgElements = element.getElementsByTagName('svg');
  Array.from(svgElements).forEach((svg) => {
    const svgData = new XMLSerializer().serializeToString(svg);
    const svgBase64 = 'data:image/svg+xml;base64,' + btoa(unescape(encodeURIComponent(svgData)));
      
    // Create a new <img> element
    const img = document.createElement('img');
    img.src = svgBase64;
    img.width = svg.clientWidth;
    img.height = svg.clientHeight;
    // Replace <svg> with <img>
    svg.parentNode.replaceChild(img, svg);
  });

  html2pdf()
    .set({
      margin: 0.5,
      filename: 'figure.pdf',
      image: {type: 'pdf', quality: 0.7},
      html2canvas: { 
          scale: 3, 
          useCORS: true,
          svgRendering: true,
      },
      jsPDF: {unit: 'in', format: 'letter', orientation: 'portrait'}
    })
    .from(element)
    .save()
    .then(() => {
      $("body").removeClass("wait");
      //document.body.removeChild(element)
    })
}

function addHTMLElements(container_div_id, args) {
  var div_container = document.getElementById(container_div_id)
  div_container.innerHTML = ""
  var tree_cohort_div_id = args.tree_cohort_div_id
  var outer_div = document.createElement('div')

  // Button zoom in.
  var button_zoomIn = createActionIcon("fa fa-search-plus")
  button_zoomIn.addEventListener('click', (event) => {
    var target_div = document.getElementById(event.currentTarget.target_div_id)
    target_div.style.zoom = parseFloat(target_div.style.zoom) + 0.1
  })
  button_zoomIn.target_div_id = tree_cohort_div_id
  outer_div.appendChild(button_zoomIn)

  // Button zoom out.
  var button_zoomOut = createActionIcon("fa fa-search-minus")
  button_zoomOut.addEventListener('click', (event) => {
    var target_div = document.getElementById(event.currentTarget.target_div_id)
    target_div.style.zoom = parseFloat(target_div.style.zoom) - 0.1
  })
  button_zoomOut.target_div_id = tree_cohort_div_id
  outer_div.appendChild(button_zoomOut)

  // Button zoom reset. 
  var button_zoomReset = createActionIcon("fa fa-expand-arrows-alt")
  button_zoomReset.addEventListener('click', (event) => {
    var target_div = document.getElementById(event.currentTarget.target_div_id)
    target_div.style.zoom = DEFAULT_ZOOM
  })
  button_zoomReset.target_div_id = tree_cohort_div_id
  outer_div.appendChild(button_zoomReset)

  // Button download.
  var button_dwl = createActionIcon("fa fa-download")
  button_dwl.id = "dwl_button"
  addInfoBoxToElement(button_dwl, "Export each tree cohort view as a camera-ready PDF figure.",
      bg_color="#0868d2", width=115, margin_left="", position="top", line_height="17px")
  button_dwl.addEventListener('click', (event) => {
    var target_div = document.getElementById(event.currentTarget.target_div_id)
    downloadDiv(target_div)
  })
  button_dwl.target_div_id = tree_cohort_div_id
  outer_div.appendChild(button_dwl)

  // Trev view buttons.
  var tree_view_div = createBlueBorder() 
  var tree_view_button = document.createElement('button')
  tree_view_button.id = "tree_view_button"
  tree_view_button.className = "button-15"
  tree_view_button.addEventListener('click', ()=>{ populateTreeView(args); })
  tree_view_button.innerHTML = '<i class="fa fa-tree" style="font-size:19px"></i>&nbsp;TREE VIEW'
  addInfoBoxToElement(tree_view_button, "Show mutation trees side by side, grouped by a given clustering (if provided). " + 
      "Nodes correspond to subclones of different sizes. Each node is labeled with the provided" +
      "set of gene mutations (SNVs, CNAs, etc) acquired by the subclone - labels are also displayed on the incoming edges. "+
      "Matching subclones and conserved edges are highlighted. Neutral clones (if labeled) are colored in " +
      "<font color=lightyellow><b>lightyellow</b></font>.", bg_color="#0868d2", width=365, margin_left=135, position="top", line_height="13px")
  tree_view_div.appendChild(tree_view_button)

  // Button sort trees.
  var button_sortTrees_id = "sort_trees"
  var button_sortTrees = createActionIcon("fa fa-sort-alpha-down", button_sortTrees_id)
  addInfoBoxToElement(button_sortTrees, "Change tree order:<br/><b>clustered</b> (default),<br/><b>alphabetical</b>," +
      "<br/>or <b>random</b>.", bg_color="#0868d2", width=115, margin_left="", position="top", line_height="17px")
  button_sortTrees.addEventListener('click', (event) => {
    var this_button = document.getElementById(event.currentTarget.id)
    var state = event.currentTarget.state
    icon_class_state_0 = "fa fa-sort-alpha-down"
    icon_class_state_1 = "fa fa-random"
    icon_class_state_2 = "fa fa-refresh"
    if (state == 0) {
      this_button.innerHTML = this_button.innerHTML.replace(icon_class_state_0, icon_class_state_1)  
    } else if (state == 1) {
      this_button.innerHTML = this_button.innerHTML.replace(icon_class_state_1, icon_class_state_2)    
    } else if (state == 2) {
      this_button.innerHTML = this_button.innerHTML.replace(icon_class_state_2, icon_class_state_0)    
    }
    var next_state = (state + 1) % 3
    this_button.state = next_state
    args = event.currentTarget.args
    args.sorting = next_state 
    populateTreeView(args)
  })
  button_sortTrees.id = button_sortTrees_id
  button_sortTrees.state = 0  
  button_sortTrees.args = args
  appendSpace(tree_view_div)
  tree_view_div.appendChild(button_sortTrees)

  // Button matching.
  var button_matching_id = "matching"
  var button_matching = createActionIcon("fa fa-times-circle", button_matching_id)
  addInfoBoxToElement(button_matching, "Highlight <b>matching subclones and conserved branches</b> (default), " +
      "<b>conserved branches only</b>, or <b>matching subclones only</b>. Matching nodes have the same " +
      "<i>matching_label</i>.", bg_color="#0868d2", width=190, margin_left="", position="top", line_height="17px")
  button_matching.addEventListener('click', (event) => {
    var this_button = document.getElementById(event.currentTarget.id)
    var state = event.currentTarget.state
    icon_class_state_0 = "fa fa-times-circle"
    icon_class_state_1 = "fa fa-circle"
    icon_class_state_2 = "fas fa-refresh"
    if (state == 0) {
      this_button.innerHTML = this_button.innerHTML.replace(icon_class_state_0, icon_class_state_1)
    } else if (state == 1) {
      this_button.innerHTML = this_button.innerHTML.replace(icon_class_state_1, icon_class_state_2)
    } else if (state == 2) {
      this_button.innerHTML = this_button.innerHTML.replace(icon_class_state_2, icon_class_state_0)
    }
    var next_state = (state + 1) % 3
    this_button.state = next_state
    args = event.currentTarget.args
    args.matching = next_state
    populateTreeView(args)
  })
  button_matching.id = button_matching_id
  button_matching.state = 0
  button_matching.args = args
  tree_view_div.appendChild(button_matching) 
  outer_div.appendChild(tree_view_div)
   
  // Heatmap view buttons.
  if ("pairwise_tree_distances" in args.data) { 
    var heatmap_view_div = createBlueBorder()
    var heatmap_view_button = document.createElement('button')
    heatmap_view_button.className = "button-15"
    heatmap_view_button.id = "heatmap_view_button"
    heatmap_view_button.addEventListener('click', ()=>{ populateHeatmapView(args) })
    heatmap_view_button.innerHTML = '<i style="font-size:19px" class="fa fa-th-large"></i> HEATMAP VIEW'
    addInfoBoxToElement(heatmap_view_button, "Show a heatmap visualization of the given pairwise distances between the mutation trees.",
        bg_color="#0868d2", width=200, margin_left="", position="top", line_height="17px")
    heatmap_view_div.appendChild(heatmap_view_button)

    // Button remove rectangles.
    var button_removeRect_id = "remove_rect"
    var button_removeRect = createActionIcon("fa fa-square-o", button_removeRect_id)
    addInfoBoxToElement(button_removeRect, "Highlight / Remove clusters.", 
         bg_color="#0868d2", width=157, margin_left="", position="top", line_height="17px")
    button_removeRect.addEventListener('click', (event) => {
      var this_button = document.getElementById(event.currentTarget.id)
      var state = event.currentTarget.state
      icon_class_state_0 = "fa fa-square-o"
      icon_class_state_1 = "fa fa-square"
      if (state == 0) {
        this_button.innerHTML = this_button.innerHTML.replace(icon_class_state_0, icon_class_state_1)
      } else if (state == 1) {
        this_button.innerHTML = this_button.innerHTML.replace(icon_class_state_1, icon_class_state_0)
      }
      var next_state = (state + 1) % 2
      this_button.state = next_state
      args = event.currentTarget.args
      args.rectangle = next_state
      populateHeatmapView(args)
    })
    button_removeRect.id = button_removeRect_id
    button_removeRect.state = 0 
    button_removeRect.args = args
    appendSpace(heatmap_view_div)
    heatmap_view_div.appendChild(button_removeRect)
    appendSpace(outer_div)
    appendSpace(outer_div) 
    outer_div.appendChild(heatmap_view_div)    
    div_container.appendChild(outer_div) 

    var umap_view_button = document.createElement('button')
    umap_view_button.className = "button-15"
    umap_view_button.id = "umap_view_button"
    umap_view_button.addEventListener('click', ()=>{ populate2DView(args); })
    umap_view_button.innerHTML = '<i style="font-size:19px" class="fa fa-dot-circle-o"></i> 2D VIEW'
    addInfoBoxToElement(umap_view_button, "Show a 2D projection of the tree points" +
        "based on a given tree pairwise distances, using Multidimensional scaling (MDS).",
        bg_color="#0868d2", width=190, margin_left="", position="top", line_height="17px")
    appendSpace(outer_div)
    outer_div.appendChild(umap_view_button)    
  }
  div_container.appendChild(outer_div)

  var outer_div = document.createElement('div')
  var width_outer_div = outer_div.clientWidth
  // Tree cohort div.
  tree_cohort_div = createDivContainer(tree_cohort_div_id)
  tree_cohort_div.style.zoom = 1
  tree_cohort_div.style.float = "left"
  tree_cohort_div.style.position = "relative"
  tree_cohort_div.style.padding = "3px"    
  tree_cohort_div.style.paddingTop = "10px" 
  outer_div.appendChild(tree_cohort_div)

  // Tree info div.
  var tree_info_div_id = args.tree_info_div_id
  tree_info_div = createDivContainer(tree_info_div_id)
  //tree_info_div.style.border = "1px solid darkgray"
  tree_info_div.style.borderRadius = "8px"
  tree_info_div.style.fontSize = "13px"
  tree_info_div.style.position = "fixed"
  tree_info_div.style.width = "max(24%, 300px)" 
  tree_info_div.style.marginTop = "10px"

  tree_info_div.style.padding = "15px"
  tree_info_div.style.float = "right"
  //tree_info_div.style.left = "calc(73% + 2px)"
  tree_info_div.style.right = "15px"
  tree_info_div.style.bottom = "10px"
  tree_info_div.style.overflowX = "scroll"
  tree_info_div.style.overflowY = "scroll"
  tree_info_div.style.overflow = "auto"
  tree_info_div.style.direction = "rtl"
  tree_info_div.style.resize = "horizontal"
  tree_info_div.innerHTML = ""
  outer_div.appendChild(tree_info_div)
  div_container.appendChild(outer_div) 
  var div_tree_cohort_top_offset = tree_cohort_div.getBoundingClientRect().top
  tree_info_div.style.top = div_tree_cohort_top_offset + "px"

  tree_cohort_div_width = window.innerWidth - tree_info_div.clientWidth - 30
  tree_cohort_div_percent = tree_cohort_div_width * 100 / window.innerWidth
  tree_cohort_div.style.width = tree_cohort_div_percent.toString() + "%"

  // Adjust initial zoom.
  canvas_width = tree_cohort_div.clientWidth - 6
  remaining_margin_right = canvas_width % (TREE_CANVAS_WIDTH + 40)
  if (remaining_margin_right > 30) {
    num_tree_on_row = Math.floor(canvas_width / (TREE_CANVAS_WIDTH + 2 * TREE_CANVAS_PADDING))
    // Make it fit one more tree.
    new_tree_canvas_size = Math.floor(canvas_width / (num_tree_on_row + 1))
    DEFAULT_ZOOM = new_tree_canvas_size / (TREE_CANVAS_WIDTH + 2 * TREE_CANVAS_PADDING) - 0.01
    tree_cohort_div.style.zoom = DEFAULT_ZOOM
  }
}

///////////////////
//// Tree view ////
///////////////////
function populateTreeView(args) {
  async_func(args, populateTreeView_slow)
}

function populateTreeView_slow(args){
  // Arguments:
  // trees: tree information; 
  // clusters: list of lists with sample names which match the tree keys;
  // args: dictionary with values for:
  //   - data (with trees, clusters, knn),
  //   - tree_cohort_div_id and tree_info_div_id,
  //   - sorting: 0=clusters, 1=alphabetic order, 2=random,
  //   - matching: 0=matching nodes, 1=matching edges, 2=both.

  // Data.
  data = args.data
  trees = data["trees"]
  clusters = [Object.keys(trees)] // One big cluster.
  if (args.sorting == 1) {
    clusters = [clusters[0].sort()]
  } else if (args.sorting == 2) { 
    shuffleArray(clusters[0])
  } else if ("clusters" in data && data["clusters"].length != 0) {
    clusters = data["clusters"]
  }

  // Prepare tree cohort container.
  var tree_cohort_div_id = args.tree_cohort_div_id
  div_tree_cohort = document.getElementById(tree_cohort_div_id)
  div_tree_cohort.innerHTML = ""
  var tree_info_div_id = args.tree_info_div_id
  tree_info_div = document.getElementById(tree_info_div_id)
  tree_info_div.innerHTML = '<div style="text-align:center;height:100%;display:flex;flex-direction:row;' +
      'align-items:center;justify-content:center;"><i>The tree clusters are indicated by<br/>' +
      ' different background <font color=#A87676>c</font><font color=#E493B3>o</font><font color=#B784B7>l</font><font color=#8E7AB5>o</font>'+
      '<font color=#F6995C>r</font><font color=#88AB8E>s</font>.&lrm;<br/>' +
      ' Zoom out <i class="fa fa-search-minus"></i> to get the<br/> full overview of the tree clusters.&lrm;<br/><br/>***<br/><br/>' +
      ' Click on the trees to visualize<br/>details of each subclone.&lrm;<br/><br/>' +
      ' Click on the " <i class="fa fa-desktop fa-sm"></i> show details of cluster X" icons<br/>' +
      ' to visualize additional cluster information.&lrm;'+
      '</i></div>'
  tree_info_div.style.backgroundColor = "white"

  for (const [i, cluster] of clusters.entries()) {
    for (sample_name of cluster) {
      var node_list = getTreeNodes(trees[sample_name]["tree"]) 
      for(node of node_list) {
        node.data.color = "white"
        node.data.conserved_parent_node_branch = false
        node.data.display_node_label = false
        if (args.sorting == 0 && clusters.length > 1){ // clusters
          if (args.matching == 0) { // display node colors and conserved branches
            node.data.color = node.data.color_in_cluster
            node.data.display_node_label = node.data.display_node_label_cluster
            if (node.data.conserved_parent_node_edge_cluster) {
              node.data.conserved_parent_node_branch = true
            }
          } 
          else if (args.matching == 1) { // conserved edges only
            if (node.data.conserved_parent_node_edge_cluster) {
              node.data.conserved_parent_node_branch = true
            }
          }
          else if (args.matching == 2) { // display node colors only
            node.data.color = node.data.color_in_cluster
            node.data.display_node_label = node.data.display_node_label_cluster
          } 
        } else { // alternative sorting
          if (args.matching == 0 || args.matching == 2) {
            node.data.display_node_label = node.data.display_node_label_cohort
            node.data.color = node.data.color_in_cohort
          }
        }
      }
    }

    // Populate click events.
    cluster.forEach(function(sample_name, j) {
      tree_data = trees[sample_name]
      tree_json = tree_data["tree"]
      if (args.sorting == 0 && clusters.length > 1) {
        cluster_color = tree_data["cluster_color"]
      } else {
        cluster_color = "white"
      }
      tree_metadata = {}
      if ("metadata" in tree_data) {
        tree_metadata = tree_data["metadata"]
      }

      // Create div structure: <div><cohort_div><outer_div><tree_div></div></div></div>
      var outer_div = document.createElement("div")
      outer_div.style.display = "inline-block"
      outer_div.style.backgroundColor = cluster_color
      div_tree_cohort.appendChild(outer_div)

      if (clusters.length > 1 && j == 0 && cluster.length > 1) {
        // Show cluster details button.
        var click_cluster_details_div = document.createElement("div")
        click_cluster_details_div.id = "cluster_details_" + i
        click_cluster_details_div.style.cursor = "pointer"
        var text_color = tinycolor(cluster_color).darken(30).desaturate(40).toHexString()
        click_cluster_details_div.innerHTML = '&nbsp;<b><i class="fa fa-desktop fa-sm" style="color:' + text_color +
          '"></i> <i><font color="' + text_color  + '"> &thinsp; show details of cluster ' + i + ' </i></b>'

        args_cluster = {}
        args_cluster["matching_nodes_details"] = trees[cluster[0]]["matching_nodes_cluster_details"]
        args_cluster["tree_info_div_id"] = tree_info_div_id
        args_cluster["cluster_bg_color"] = cluster_color
        args_cluster["cluster_metadata"] = trees[cluster[0]]["sample_metadata_colors"]
        args_cluster["table_color_codes"] = trees[cluster[0]]["table_color_codes"]
        args_cluster["highlighted_genes"] = data["highlighted_genes"];
        (function(args_cluster) {
          click_cluster_details_div.addEventListener('click', () => { showClusterInfo(args_cluster) })
        })(args_cluster);
        outer_div.appendChild(click_cluster_details_div)
      }  else {
        outer_div.innerHTML = "<br/>"
      } 
      var timestamp = Date.now()
      var tree_div_id = sample_name + "_" + timestamp + "_tree"
      var tree_div = createDivContainer(tree_div_id)
      var border_color = tinycolor(cluster_color).darken(5).desaturate(20).toHexString()
      tree_div.setAttribute("style", "border: 1px solid " + border_color + 
          ";cursor: pointer; border-radius: 8px; display: inline-block; background-color:" + cluster_color)
      tree_div.addEventListener('click', ()=>{ showTreeInfo(sample_name, args); }) 
      $('#' + tree_div_id).lazyload()

      outer_div.appendChild(tree_div)
      displayTree(tree_div_id, sample_name, tree_json, "", "", null, tree_info_view=false)
    })
  }

  $(function(){
    var container = $('#' + tree_info_div_id);
    var initial_position = container.position().top
  
    $(document).scroll(function() {
      scroll_offset = $(document).scrollTop()
      container.css('top', Math.max(0, initial_position - scroll_offset));
    });
  });
}

//////////////
//// MISC ////
//////////////
background_colors = ["#F9F6EE", "#FFFDD0", "#FFD7BA", // yellow
                     "#e1f8fc", "#C3E2E8", "#DEFDEF", "#c2e6c1", // blue and green
                     "#F8EAEC", "#FFD3DD", "#f8ebfc", "#E0D9E4"] // red and violet

node_colors = ["#A87676", "#E493B3",
        "#B784B7", "#8E7AB5", "#F6995C", "#EEC759",
        "#E9B384", "#88AB8E", "#4f6f52", "#65647C",
        "#8B7E74", "#FF8Dc7", "#A7D2CB", "#554994"]

const sleep = (milliseconds) => {
  return new Promise(resolve => setTimeout(resolve, milliseconds))
}

function deepCopy(oldValue) {
  var newValue
  strValue = JSON.stringify(oldValue)
  return newValue = JSON.parse(strValue)
}

function mapSize(map) {
  return Object.keys(map).length
} 

function shuffleArray(array) {
  for (let i = array.length - 1; i > 0; i--) {
    const j = Math.floor(Math.random() * (i + 1));
    [array[i], array[j]] = [array[j], array[i]];
  }
}

function linspace(startValue, stopValue, cardinality) {
  var arr = [];
  var step = (stopValue - startValue) / (cardinality - 1);
  for (var i = 0; i < cardinality; i++) {
    arr.push(startValue + (step * i));
  }
  return arr;
}

function capitalizeFirstLetter(string) {
    return string.charAt(0).toUpperCase() + string.slice(1);
}

function round(value, precision) {
    var multiplier = Math.pow(10, precision || 0);
    return Math.round(value * multiplier) / multiplier;
}

/////////////
//// kNN ////
/////////////

// Descendents without including the start node.
function getDescendents(node) { 
  var treemap = d3.tree()
  var descendents = treemap(node).descendants()
  return descendents.slice(1,descendents.length+1)
}

// Ancestors without including the start node.
function getAncestors(start_node) {
  var treemap = d3.tree()

  var nodes_to_remove = treemap(start_node).descendants()
  var node_ids_to_remove = nodes_to_remove.map((node) => node.data.node_id)

  var root = null
  for (let node of treemap(start_node).ancestors()) {
    if (!node.parent) {
      root = node
      break
    }
  }

  node_list = new Array()
  var root_copy = root.copy() 
  for (let node of treemap(root_copy).descendants()) {
    // Exclude the start node subtree and the root.
    if (!node.parent) { // root
      continue
    }
    if (node.data.node_id == start_node.data.node_id) {
      if (node.parent) {
        node.parent.children = null
      }
    }
    if (node_ids_to_remove.includes(node.data.node_id)) {
      continue
    }
    if (node.parent.data.node_id == root_copy.data.node_id){ // child of the root becomes new root
      node.parent = null
    }
    node_list.push(node)
  }
  return node_list
}

function getMatchingNodesWithMatchingLabels(node_list_1, node_list_2) {
  if (node_list_1.length == 0 || node_list_2.length == 0) {
    return [null, null]
  }

  node_map_1 = getMalignantMatchingLabels(node_list_1)
  node_map_2 = getMalignantMatchingLabels(node_list_2)
  matching_labels = new Set(Array.from(node_map_1.keys())).intersection(new Set(Array.from(node_map_2.keys()))) 
  if (matching_labels.size == 0) {
    return [null, null]
  }

  match_depth = new Map()
  for (const label of matching_labels) {
    const depth = node_map_1.get(label).depth + node_map_2.get(label).depth
    if (!match_depth.has(depth)) {
      match_depth.set(depth, [label])
    } else {
      match_depth.get(depth).push(label)
    }
  }

  const min_depth = Math.min(...Array.from(match_depth.keys()))
  const label_min_depth = match_depth.get(min_depth)[0]
  const node_1 = node_map_1.get(label_min_depth)
  const node_2 = node_map_2.get(label_min_depth)
  const matching_pair = new Array(node_1, node_2)
  
  // Match nodes in children and ancestor subtrees.
  const node_1_wo_parent = node_1.copy() // deep copy and node is root
  const node_2_wo_parent = node_2.copy() // deep copy and node is root
  const matching_pair_descendents = getMatchingNodesWithMatchingLabels(getDescendents(node_1_wo_parent), getDescendents(node_2_wo_parent)) 
  const matching_pair_ancestors = getMatchingNodesWithMatchingLabels(getAncestors(node_1), getAncestors(node_2))

  var matching_pairs = new Array()
  matching_pairs.push(matching_pair)
  if (matching_pair_descendents[0]) {
    for (let pair of matching_pair_descendents) {
      matching_pairs.push(pair)
    }
  } 
  if (matching_pair_ancestors[0]) {
    for (let pair of matching_pair_ancestors) {
      matching_pairs.push(pair)
    }
  }  
  return matching_pairs 
}



//////////////////////////
//// UTILS D3.JS TREE ////
//////////////////////////
function getTreeNodes(tree_json) {
  var node_hierarchy = d3.hierarchy(tree_json)
  var treemap = d3.tree()
  return treemap(node_hierarchy).descendants()
}

function getMalignantMatchingLabels(node_list){
  var node_map = new Map()
  for (node of node_list) {
    if ((!node.parent && mapSize(node.data.gene_states) && node.data.matching_label) ||
        (node.parent && node.data.matching_label && !node.data.is_neutral)) {
      node_map.set(node.data.matching_label, node)
    }
  }
  return node_map
}

function getTreeMalignantMatchingLabels(tree_json){
  return getMalignantMatchingLabels(getTreeNodes(tree_json))
}

// Returns a list of node-parent strings.
function getBranches(tree_json) {
  branches = new Set()
  node_list = getTreeNodes(tree_json)
  for (node of node_list) {
    if(node.parent && node.parent.data.matching_label && node.data.matching_label) {
      branches.add(node.parent.data.matching_label + "_" + node.data.matching_label)
    }
  }
  return branches
}

function CNValueToInt(value) {
  if (value == "+") {
    return 1
  } else if (value == "-") {
    return -1
  } else {
    return parseInt(value)
  }
}

function updateGeneStates(gene_events_map, parent_gene_states) {
  var new_gene_states = deepCopy(parent_gene_states)
  // Update the new_gene_states with each current gene event.
  for (const [gene, events_map] of Object.entries(gene_events_map)) {
    if (!(gene in new_gene_states)) {
      new_gene_states[gene] = deepCopy(events_map)
    } else {
      for (const [old_event_type, old_value] of Object.entries(new_gene_states[gene])) {
        for (const [new_event_type, new_value] of Object.entries(events_map)) {
          // Apply each event from gene_events
          if(old_event_type == new_event_type) {
            if (old_event_type == "CNA") {
              var new_cn_value = CNValueToInt(old_value) + CNValueToInt(new_value)
              if (new_cn_value == 0) {
                delete new_gene_states[gene]["CNA"]
                if (mapSize(new_gene_states[gene]) == 0) {
                  delete new_gene_states[gene]
                }
              } else {
                events_map["CNA"] = new_cn_value
              }
            } else {
              new_gene_states[gene][new_event_type] =  old_value + ";" + new_value
            }
          } else {
            if (!(gene in new_gene_states)) {
              new_gene_states[gene] = {}
            }
            new_gene_states[gene][new_event_type] = new_value
          }
        }
      }
    }
  }
  return new_gene_states
}

function populateGeneStates(tree){
  var node_list = getTreeNodes(tree) // Nodes are in breadth first order.
  for(node of node_list) {
    // Populate node with the summaries for the events (to be displayed on the tree branch) and
    // propagate the events from parents to the children (gene states).
    if (node.data.gene_events){
      const [gene_events, gene_events_with_details] = getGeneCategoriesInNode(node)
      node.data.gene_event_categories = gene_events
      node.data.gene_event_categories_details = gene_events_with_details
      node_gene_events = deepCopy(node.data.gene_events)
      if (!node.parent) { // root
        node.data.gene_states = node_gene_events
      } else {
        node.data.gene_states = updateGeneStates(node_gene_events, node.parent.data.gene_states)
      }
    } else {
      node.data.gene_states = new Map()
    }
  }
}

//////////////////////
//// HTML helpers ////
//////////////////////
function appendLineBreak(div, num=1) {
  div.appendChild(document.createElement("p"))
  if (num > 1) {
    var br_div = document.createElement("div")
    br_div.innerHTML = "<br/>"
    for (let i= 1; i<num; i++) {
      div.appendChild(br_div)
    }
  }
}

function appendHalfLineBreak(div) {
  div.innerHTML += "<p style='margin:6px;'></p>"
}

function appendSpace(div) {
  new_div = document.createElement("div")
  new_div.style.display = "inline-block"
  new_div.innerHTML = "&nbsp;"
  div.appendChild(new_div)
} 

function createDivContainer(id) {
  var div_container = document.createElement('div')
  div_container.setAttribute("id", id)
  return div_container
}

function addInfoBoxToElement(elem, text, bg_color="#0868d2", width=250, margin_left="", position="top", line_height="17px") {
  var span = document.createElement('span')
  span.style.backgroundColor = bg_color
  span.style.width = width + "px"
  span.style.lineHeight = line_height
  if (position == "top" || position == "bottom") {
    if (margin_left=="") {
      span.style.marginLeft = "-" + width/2 + "px"
    } else {
      span.style.marginLeft = "-" + margin_left + "px"
    }
  }

  span.style.whiteSpace = "normal"
  span.classList.add("tooltiptext-" + position)
  span.innerHTML = text
  elem.classList.add("info-tooltip")
  elem.appendChild(span)
}

function createInfoTooltip() {
  var div = document.createElement('span')
  div.classList.add("info-tooltip")
  div.innerHTML = "&nbsp;<i class='fa fa-question-circle' style='color:onyx; font-size:16px;'></i>"
  return div
}

////////////////////////////
//// Populate tree info ////
////////////////////////////
function createInfoHeader(html_text, color_motif="#d2cae6"){
  div_container = document.createElement('div')
  div_container.style.backgroundColor = color_motif
  div_container.innerHTML = "&nbsp;" + html_text
  return div_container
}

function createExpandBox(div_id, reverse=false) {
  var button_class = "fa-caret-down"
  var show_div = true
  if(reverse) {
    button_class = "fa-caret-up"
    show_div = false
  }
  var button_expand = document.createElement('span')
  var button_id = "button_" + div_id
  button_expand.setAttribute("id", button_id)
  button_expand.innerHTML += "&nbsp;<a><i class='fa " + button_class + " fa-lg' style='color:onyx; cursor:pointer'></i></a>"
  button_expand.addEventListener('click', (event) => {
    var div = document.getElementById(event.currentTarget.div_id)
    var button = document.getElementById(event.currentTarget.button_id)
    var show_div = event.currentTarget.show_div
    if (show_div) {
      div.style.display = 'block'
      button.innerHTML = button.innerHTML.replace("down", "up")
      info_text = "&nbsp;&nbsp;Click to hide."
    } else {
      div.style.display = 'none'
      button.innerHTML = button.innerHTML.replace("up", "down")
      info_text = "Click to expand."
    }
    button.show_div = !show_div
    addInfoBoxToElement(button, info_text,
          bg_color="#353935", width=95, margin_left="", position="top", line_height="17px")
  })
  button_expand.div_id = div_id
  button_expand.button_id = button_id
  button_expand.show_div = show_div

  info_text = "Click to expand."
  if(reverse) {
    info_text = "&nbsp;&nbsp;Click to hide."
  }
  addInfoBoxToElement(button_expand, info_text,
        bg_color="#353935", width=95, margin_left="", position="top", line_height="17px")

  return button_expand
}

function showTreeInfo(sample_name, args) {
  // Input: sample_name and a dictionary with the following attached information:
  // sample name, tree hierarchy, tree metadata, matching trees map,
  // div_id (where the tree info is displayed).

  // Tree information
  var tree_object = data["trees"][sample_name]["tree"]
  var metadata = data["trees"][sample_name]["metadata"]
  var knn = data["matching_trees"]

  // Div container.
  var tree_info_div_id = args.tree_info_div_id
  var tree_info_div = document.getElementById(tree_info_div_id)
  tree_info_div.innerHTML = ""
  tree_info_div.style.border = "1px solid lightgray"
  tree_info_div.style.backgroundColor = "white"
  tree_info_div.style.width = "24%"
  tree_info_div.style.minWidth = "300px"

  // Compute gene list and variables related to the selected target gene. 
  const [gene_categories, gene_categories_with_details] = getGeneCategoriesInTree(tree_object)
  var drug_gene_map = {}
  var drug_name_map = {}
  var drug_list = []
  if (gene_categories.size) {
    var gene_list = Array.from(Array.from(gene_categories.values())[0])
    var gene_list_with_details = Array.from(Array.from(gene_categories_with_details.values())[0])
    var drug_gene_map = getDrugToGeneMap(gene_list, gene_drug_map)
    var drug_list = findDrugsAffectingGeneList(gene_list, drug_gene_map, tree_object)
  }

  // Display sample name.
  header_sample_name = createInfoHeader("<strong><b>" + sample_name + "</b></strong>")
  header_sample_name.style.direction = "ltr"
  tree_info_div.appendChild(header_sample_name)
  appendLineBreak(tree_info_div)

  // Gene selection.
  if (gene_categories.size) {
    gene_selection_dropdown = document.createElement("select")
    gene_selection_dropdown.setAttribute("id", "gene_selection")
    gene_selection_dropdown.classList.add("gene_selection")
    gene_selection_dropdown.style.width = "200px"
    option = document.createElement('option');
    option.textContent = "-- Select target gene event --"
    gene_selection_dropdown.appendChild(option);
    keys = Array.from(gene_categories.keys()).sort()
    for (event of keys) {
      array = gene_categories.get(event)
      sorted_genes = Array.from(array).sort()
      addGeneToGeneList(gene_selection_dropdown, sorted_genes, event)
    }
    div_container = document.createElement('div')
    div_container.style.direction = "ltr"
    div_container.style.width = "fit-content"
    div_container.style.margin = "auto"
    div_container.appendChild(gene_selection_dropdown)

    var info_icon = createInfoTooltip("fa fa-question-circle")
    info_text = "Subclones affected by <b>selected target gene</b> are highlighted with colors in the mutation tree below: " +
      "<b><font color=tomato>red</font></b>  for CN amplification, <b><font color=lightsteelblue>blue</font></b> for CN deletion " +
      "and <b><font color=#b4a7d6>violet</font></b> for any other mutation event. An event in one subclone of the tree " +
      "affects all the subsequent nodes in the child subtree." 
    addInfoBoxToElement(info_icon, info_text, bg_color="#353935", width=200, margin_left="", position="left", line_height="17px")
    div_container.appendChild(info_icon)

    tree_info_div.appendChild(div_container)
    appendHalfLineBreak(tree_info_div)

    // Drug selection.
    drug_selection_dropdown = document.createElement("select")
    drug_selection_dropdown.setAttribute("id", "drug_selection")
    drug_selection_dropdown.classList.add("drug_selection")
    drug_selection_dropdown.style.width = "200px"
    option = document.createElement('option');
    option.textContent = "-- Select target drug --"
    drug_selection_dropdown.appendChild(option)
    drug_name_map = {}
    for (idx = 0; idx < drug_list.length; idx++) {
      drug_short_name = drug_list[idx].substring(0,20)
      if (drug_short_name.length < drug_list[idx].length) {
        drug_short_name += "..."
      }
      drug_name_map[drug_short_name] = drug_list[idx]
      drug_selection_dropdown.options[drug_selection_dropdown.options.length] = new Option(drug_short_name, drug_short_name);
    }
    div_container = document.createElement('div')
    div_container.style.direction = "ltr"
    div_container.style.width = "fit-content"
    div_container.style.margin = "auto"
    div_container.appendChild(drug_selection_dropdown)

    var info_icon = createInfoTooltip("fa fa-question-circle")
    info_text = "Subclones with affected genes that have a theoretical " +
      "interaction with the <b> selected target drug</b> (according to DGIdb) are indicated with green squares in the mutation tree below. " + 
       "Drugs are listed in descending order by the number of cells they affect. Only drug interactions with at least one citation are considered."
    addInfoBoxToElement(info_icon, info_text, bg_color="#353935", width=230, margin_left="", position="left", line_height="17px")
    div_container.appendChild(info_icon)

    tree_info_div.appendChild(div_container)
    appendLineBreak(tree_info_div)

    target_gene = getSeletedItem(gene_selection_dropdown)
    target_drug = drug_name_map[getSeletedItem(drug_selection_dropdown)]
  }

  // Display tree.
  var tree_box_id = "tree_box"
  var tree_box_div = createDivContainer(tree_box_id)
  tree_box_div.style.direction = "ltr"
  tree_box_div.style.display = 'block'
  tree_info_div.appendChild(tree_box_div)
  appendLineBreak(tree_info_div)
  tree_div_height = Math.min(400, 2*screen.height/3)
  displayTree(tree_box_id, sample_name, tree_object, target_gene, target_drug, drug_gene_map, height=tree_div_height)

  // Displayed metadata.
  if (metadata && mapSize(metadata)) {
    var header_metadata = createInfoHeader("<b>Metadata</b>")
    header_metadata.style.direction = "ltr"
    tree_info_div.appendChild(header_metadata)

    var metadata_box_id = "metadata"
    var metadata_box_div = createDivContainer(metadata_box_id)
    metadata_box_div.style.direction = "ltr"
    metadata_box_div.style.display = 'none'
    for (const [key, value] of Object.entries(metadata)) {
      metadata_box_div.innerHTML += '<i><b>' + key + '</b>: ' + value + '</i></br>'
    }
    tree_info_div.appendChild(metadata_box_div)
    header_metadata.appendChild(createExpandBox(metadata_box_id))
    appendLineBreak(tree_info_div)
  }

  // Top DGIdb drugs.
  if (drug_list.length) {
    var header_top_drugs = createInfoHeader("<b>Top DGIdb drugs associated with target gene</b>")
    header_top_drugs.style.direction = "ltr"
    tree_info_div.appendChild(header_top_drugs)

    var top_drugs_box_id = "top_drugs"
    var top_drugs_box_div = createDivContainer(top_drugs_box_id)
    top_drugs_box_div.style.direction = "ltr"
    top_drugs_box_div.style.display = 'none'
    top_drugs_box_div.innerHTML = populateGeneDrugInfoHTML(target_gene, gene_drug_map)
    tree_info_div.appendChild(top_drugs_box_div)
    header_top_drugs.appendChild(createExpandBox(top_drugs_box_id))
    appendLineBreak(tree_info_div)
  }

  // kNN matching trees.
  var header_knn = createInfoHeader("<b>K-nearest tree neighbors</b>")
  header_knn.style.direction = "ltr"
  tree_info_div.appendChild(header_knn)

  var knn_box_id = "knn"
  var knn_box_div = createDivContainer(knn_box_id)
  knn_box_div.style.direction = "ltr"
  knn_box_div.style.display = 'none'
  tree_info_div.appendChild(knn_box_div)
  if (knn && mapSize(knn)) {
    knn_keys = Object.keys(knn)
    sample_matches = knn_keys.filter((key) => key.startsWith(sample_name))
    if (sample_matches.length) {
      sample_matches.sort((key_1, key_2) => knn[key_2]["similarity"] - knn[key_1]["similarity"])
      for(const [i, key] of sample_matches.entries()) {
        async_displayTreeMatching(knn_box_id, knn[key], idx=i+1)
      }
    } else {
      knn_box_div.innerHTML = "<i>No matching trees found.</i>"
    }
  } else {
    knn_box_div.innerHTML = "<i>kNN not computed either because node `matching_label` attributes are missing, or the number of trees " +
        "is too large.</i><br/>"
  }
  header_knn.appendChild(createExpandBox(knn_box_id))
  var info_icon = createInfoTooltip("fa fa-question-circle")
  info_text = "K-nearest matching trees to the selected tree based on the provided <i>matching_labels</i>, " + 
      "computed using a greedy approximation algorithm for the maximum matching problem with ordering constraints. " +
      "The panel shows pairs of matching trees (where the left tree always corresponds to the target sample) and the " +
      "corresponding matching subclones."
  addInfoBoxToElement(info_icon, info_text, bg_color="#353935", width=230, margin_left="", position="top", line_height="20px")
  header_knn.appendChild(info_icon)
  appendLineBreak(tree_info_div)

  // Add event listeners after the DOM is complete. 
  $('.gene_selection, .drug_selection').on('change', 
    {tree_div_id: tree_box_id,
     tree_object: tree_object,
     gene_dropdown_id: 'gene_selection',
     drug_dropdown_id: 'drug_selection',
     dgi_div_id: top_drugs_box_id,
     gene_drug_map: gene_drug_map,
     drug_gene_map: drug_gene_map,
     drug_name_map: drug_name_map,
     tree_div_height: tree_div_height}, 
    updateTree)
}

function getEventLabel(event, value) {
  if (event == "CNA") {
    if (parseInt(value) > 0 || value == "+") {
      return "amplified"
    }
    else if (parseInt(value) < 0 || value == "-") {
      //if (parseInt(value) < 0) {
      //  return "deleted " + parseInt(value).toString()
      //}
      return "deleted"
    }
    else {
      return value
    }
  }
  return event
}

function geneEventMapToEventGeneMap(gene_events_map) {
  map = new Map()
  for (const [gene, events_map] of Object.entries(gene_events_map)) {
    for (const [event, value] of Object.entries(events_map)) {
      var event_type = getEventLabel(event, value)
      if (!(map.has(event_type))) {
        map.set(event_type, new Set())
      }
      map_set = map.get(event_type)
      map_set.add(gene)   
      map.set(event_type, map_set)
    }
  }
  return map
}

function getGeneCategoriesInNode(node) {
  var gene_categories = new Map()
  var gene_categories_with_details = new Map()
  if (node.data.gene_events){
    for (const [gene, events] of Object.entries(node.data.gene_events)) {
      for (const [event, value] of Object.entries(events)) {
        var label = getEventLabel(event, value)
        if (!(gene_categories.has(label))) {
          gene_categories.set(label, new Set())
          gene_categories_with_details.set(label, new Set())
        }
        gene_categories.get(label).add(gene)
        if (event == "CNA") {
          gene_categories_with_details.get(label).add(gene)
        } else {
          gene_string = gene
          if (value) {
            gene_string += "_" + value
          }
          gene_categories_with_details.get(label).add(gene_string)
        }
      }
    }
  }
  return [gene_categories, gene_categories_with_details]
}

function getGeneCategoriesInTree(tree_object) {
  var node_list = getTreeNodes(tree_object)
  gene_categories = new Map()
  gene_categories_with_details = new Map()
  for (node of node_list) {
    const [node_gene_events, node_gene_events_with_details] = getGeneCategoriesInNode(node)
    for (let [label, set] of node_gene_events.entries()) {
      if (!(gene_categories.has(label))) {
        gene_categories.set(label, set)
        gene_categories_with_details.set(label, node_gene_events_with_details.get(label))
      } else {
        gene_categories.set(label, gene_categories.get(label).union(set))
        gene_categories_with_details.set(label, gene_categories.get(label).union(node_gene_events.get(label)))
      }
    }
  }
  return [gene_categories, gene_categories_with_details]
}

function computeGeneEventIntersection(list_gene_events){
  event_categories = {}
  set_0 = geneEventMapToEventGeneMap(list_gene_events[0])
  for (var gene_events of list_gene_events.slice(1,list_gene_events.length+1)){
    set = geneEventMapToEventGeneMap(gene_events)
    for (let [event, gene_set] of set.entries()) {
      if (set_0.has(event)) {
        intersection = set_0.get(event).intersection(gene_set)
        set_0.set(event,intersection)
      }
    }
  }
  return set_0
}

function getSeletedItem(select_element) {
  if(select_element.selectedIndex == 0) {
    return ""
  } else {
    return select_element.options[select_element.selectedIndex].text
  }
}

function addGeneToGeneList(html_element, gene_list, category){
  optgroup = document.createElement('optgroup')
  optgroup.label = category
  for (var gene of gene_list){
    option = document.createElement('option')
    option.textContent = gene
    option.value = gene
    optgroup.appendChild(option)
  }
  html_element.appendChild(optgroup);
}

function getCellCountForGene(tree_object, target_gene) {
  var node_list = getTreeNodes(tree_object)
  var cell_percentage = 0
  for (var node of node_list) {
    if (node.data.gene_states) {
      if (target_gene in node.data.gene_states) {
        if (node.data.size_percent) {
          cell_percentage += node.data.size_percent 
        } else {
          cell_percentage --
        }
      }
    }
  }
  return cell_percentage
}

function getCellCountForDrugInteraction(tree_object, target_drug, drug_gene_map) {
  if(!(drug_gene_map.has(target_drug))) {
    return 0
  }
  var node_list = getTreeNodes(tree_object)
  var gene_list = drug_gene_map.get(target_drug)
  var cell_percentage = 0
  for (node of node_list) {
    if (node.data.gene_events){
      node_genes = Array.from(Object.keys(node.data.gene_states))
      for (var gene_1 of gene_list){
        if (node_genes.includes(gene_1)) {
          cell_percentage += node.data.size_percent
          break
        }
      }
    }
  }
  return cell_percentage
}

function updateTree(event) {
  var tree_div_id = event.data.tree_div_id
  var tree_object = event.data.tree_object
  var gene_dropdown_id = event.data.gene_dropdown_id
  var drug_dropdown_id = event.data.drug_dropdown_id
  var dgi_div_id = event.data.dgi_div_id
  var gene_drug_map = event.data.gene_drug_map
  var drug_gene_map = event.data.drug_gene_map
  var drug_name_map = event.data.drug_name_map
  var tree_div_height = event.data.tree_div_height

  var gene_dropdown = document.getElementById(gene_dropdown_id)
  var target_gene = gene_dropdown.options[gene_dropdown.selectedIndex].text
  target_gene_cell_percent = getCellCountForGene(tree_object, target_gene)

  var drug_dropdown = document.getElementById(drug_dropdown_id)
  var target_drug = drug_name_map[drug_dropdown.options[drug_dropdown.selectedIndex].text]
  target_drug_cell_percent = getCellCountForDrugInteraction(tree_object, target_drug, drug_gene_map) 

  document.getElementById(tree_div_id).innerHTML = ""
  displayTree(tree_div_id, "", tree_object, target_gene, target_drug, drug_gene_map, height=tree_div_height) 
  
  cell_count_div_id = "cell_count_div"
  remove_element = document.getElementById(cell_count_div_id)
  if (typeof(remove_element) != 'undefined' && remove_element != null) {
    remove_element.remove();
  }
  var cell_count_div = document.createElement("div")
  cell_count_div.setAttribute("id", cell_count_div_id)
  cell_count_div.style.textAlign = "center"
  if (!isNaN(target_gene_cell_percent) && target_gene_cell_percent > 0) {
    cell_count_div.innerHTML += "<i>Target gene affected in " + round(target_gene_cell_percent*100, 1) + "%</i> cells.<br/>"
  }
  if (!isNaN(target_drug_cell_percent) && target_drug_cell_percent != 0) {
    cell_count_div.innerHTML += "<i>" + round(target_drug_cell_percent*100, 1) + "%</i> cells affected by target drug.<br/>"
  }
  document.getElementById(tree_div_id).append(cell_count_div)

  html_string = populateGeneDrugInfoHTML(target_gene, gene_drug_map)
  document.getElementById(dgi_div_id).innerHTML = html_string
} 

function populateGeneDrugInfoHTML(target_gene, gene_drug_map) {
  if (target_gene == "") {
    return "<i>No target gene selected.</i><br/>"
  }
  var drugs = gene_drug_map.get(target_gene)
  if (!(gene_drug_map.has(target_gene))) {
    gene_without_delimiter = target_gene.split(/_|-/)
    if (gene_without_delimiter.length > 1 && gene_drug_map.has(gene_without_delimiter[0])) {
      drugs = gene_drug_map.get(gene_without_delimiter[0])
    }
    else {
      return "No drugs associated with target gene.<br/>"
    }
  }
  drugs.sort(function(a, b){return b["drug_score"] - a["drug_score"]});
  activators = []
  inhibitors = []
  others = []
  for (var drug of drugs) {
    //if (drug.num_citations < 3) {
    //  continue
    //}
    if(isInhibitor(drug.interaction_types)) { 
      inhibitors.push(drug)
    }
    else if (isActivator(drug.interaction_types)) {
      activators.push(drug)
    }
    else {
      others.push(drug)
    }
  }

  if (activators.length + inhibitors.length + others.length == 0) {
    return "No drugs associated with target gene in at least one citation.<br/>"
  }

  html_string = ""
  if(activators.length) {
    html_string += "<b><u>Activators:</u></b><br/>"
    for (var drug of activators) {
      html_string += "<a href='https://www.dgidb.org/results?searchType=drug&searchTerms=" + drug.drug_name + 
        "' target='dgidb'>" + drug.drug_name + "</a><br/>"
    }
  }

  if(inhibitors.length) {
    html_string += "<b><u>Inhibitors:</u></b><br/>"
    for (var drug of inhibitors) {
      html_string += "<a href='https://www.dgidb.org/results?searchType=drug&searchTerms=" + drug.drug_name +
        "' target='dgidb'>" + drug.drug_name + "</a><br/>"
    }
  }

  if (others.length  && activators.length + inhibitors.length !=0) {
    html_string += "<b><u>Others:</u></b><br/>"
  }
  if(others.length) {
    for (var drug of others) {
      html_string += "<a href='https://www.dgidb.org/results?searchType=drug&searchTerms=" + drug.drug_name +
        "' target='dgidb'>" + drug.drug_name + "</a><br/>"
    }
  }
  return html_string
}

/////////////////////
//// DGIdb utils ////
/////////////////////
function parseGeneDrugInteractions(gene_drug_interaction){
  // Returns a map with gene key and drug info value.
  gene_drug_map = new Map()
  for (var item of gene_drug_interaction["data"]["genes"]["nodes"]) {
    gene = item["name"]
    gene_drug_map.set(gene, [])
    for (var drug of item["interactions"])  {
      drug_info = {}
      drug_info["drug_name"] = drug["drug"]["name"]
      drug_info["interaction_types"] = drug.interactionTypes
      drug_info["num_citations"] = drug["publications"].length
      drug_info["drug_score"] = drug["interactionScore"]
      gene_drug_map.get(gene).push(drug_info)
    }
  }
  return gene_drug_map
}

function getDrugToGeneMap(gene_list, gene_drug_map) {
  var drug_gene_map = new Map()
  for (var gene of gene_list){
    if(gene_drug_map.has(gene)) {
      for (var drug of gene_drug_map.get(gene)){
        //if(drug["num_citations"] < 3) {
        //  continue
        //}
        drug_name = drug.drug_name.substring(0,30)
        if (drug_gene_map.has(drug_name)) {
          drug_gene_map.get(drug_name).push(gene)
        }
        else {
          drug_gene_map.set(drug_name, [gene])
        }
      }
    }
  }
  return drug_gene_map
}

function custom_compare (drug_1, drug_2) {
  return drug_1.cnt - drug_2.cnt
}

function findDrugsAffectingGeneList(tumor_gene_list, drug_gene_map, tree_object) {
  var drug_cnt_list = []
  for (let [drug, drug_gene_list] of drug_gene_map.entries()) {
    var cells_affecte_by_drug = 0
    for (var tumor_gene of tumor_gene_list) {
      if (drug_gene_list.includes(tumor_gene)) {
        count = getCellCountForGene(tree_object, tumor_gene)
        if (count < 0){ // i.e., node.size_percent not set 
          cells_affecte_by_drug += Math.abs(count)
        } else {
          cells_affecte_by_drug += count
        }
        break
      }
    }
    if (isNaN(cells_affecte_by_drug)) {
      drug_cnt_list.push({"drug": drug, "cnt": 1})
    } else {
      drug_cnt_list.push({"drug": drug, "cnt": cells_affecte_by_drug})
    }
  }
  drug_cnt_list.sort(custom_compare).reverse()
  drug_list = []
  for (item of drug_cnt_list) {
    drug_list.push(item.drug)
  }
  return drug_list
}

function isInhibitor(interaction_types) {
  for (var interaction of interaction_types){
    if (interaction["directionality"]) {
      string = interaction["directionality"].toUpperCase()
      if(string.includes("INHIBITORY")) {
        return true
      }
    }
  }
  return false
}

function isActivator(interaction_types) {
  for (var interaction of interaction_types){
    if(interaction["directionality"]) {
      string = interaction["directionality"].toUpperCase()
      if(string.includes("ACTIVATING")) {
        return true
      }
    }
  }
  return false
}

///////////////////////////////
//// Populate cluster info ////
///////////////////////////////
function applyTextStyle(gene, highlighted_genes) {
  if (highlighted_genes) {
    if (gene in highlighted_genes) {
      style = highlighted_genes[gene]
      if (gene in gene_chr_map) {
        gene = "<span style='cursor:default;' title='chr" + gene_chr_map[gene] + "'>" + gene + "</span>"
      } 
      if (style == "bold") {
        return "<b>" + gene + "</b>"
      } else if (style == "italic") {
        return "<i>" + gene + "</i>"
      } else if (style == "underline") {
        return "<u>" + gene + "</u>"
      } else {
        return "<font color=" + style + ">" + gene + "</font>"
      }
    } else {
      return "<font color=lightgray>" + gene + "</font>"
    }
  } else {
    return gene
  }
}

function showClusterInfo(args) {
  async_func(args, showClusterInfo_slow)
} 

function showClusterInfo_slow(args) {
  var cluster_bg_color = args.cluster_bg_color
  var cluster_metadata = args.cluster_metadata
  var table_color_codes = args.table_color_codes
  var node_intersections = args.node_intersections 
  var matching_nodes_details = args.matching_nodes_details
  var highlighted_genes = args.highlighted_genes

  var tree_info_div_id = args.tree_info_div_id
  var tree_info_div = document.getElementById(tree_info_div_id)
  tree_info_div.innerHTML = ""
  tree_info_div.style.backgroundColor = cluster_bg_color
  tree_info_div.style.width = "24%"
  tree_info_div.style.minWidth = "300px"

  metadata_found = false
  for (sample_name in cluster_metadata){
    if (mapSize(cluster_metadata[sample_name])) {
      metadata_found = true
      break
    }
  }
  if (!metadata_found) {
    tree_info_div.innerHTML = "<div style='text-align:center;height:100%;display:flex;flex-direction:row;" +
      "align-items:center;justify-content:center;'><i>No input clinical data provided.</i></div>"
    return; 
  }

  background_color = tinycolor(cluster_bg_color).darken(30).desaturate(40).toHexString()

  // Matching details.
  if (matching_nodes_details.size) {
    appendLineBreak(tree_info_div)
    var header_genes = createInfoHeader("<b>Matching details</b>", color_motif=background_color)
    header_genes.style.direction = "ltr"
    tree_info_div.appendChild(header_genes)

    var genes_box_id = "cluster_genes"
    var genes_box_div = createDivContainer(genes_box_id)
    genes_box_div.style.direction = "ltr"
    genes_box_div.style.display = 'none'

    // Add text.
    appendLineBreak(genes_box_div)
    for (let [color, events] of matching_nodes_details.entries()) {
      // Matching genes.
      genes_box_div.innerHTML += '<i class="fa fa-circle" style="font-size:18px;color:' + color + '"></i> &nbsp;'
      for (let [event, gene_set] of events.entries()) {
        genes = Array.from(gene_set).sort()
        genes = genes.filter(e => !e.includes('.'))
        genes_box_div.innerHTML += event + ": " + applyTextStyle(genes[0], highlighted_genes)
        for (var gene of genes.slice(1,genes.length+1)) {
          genes_box_div.innerHTML += ", " + applyTextStyle(gene, highlighted_genes)
        }
        genes_box_div.innerHTML += "; "
      }
      appendHalfLineBreak(genes_box_div)

      // Affected chromosomes.
      affected_chromosomes = new Set()
      for (let [event, gene_set] of events.entries()) {
        for (let gene of Array.from(gene_set)) {
          if (gene in gene_chr_map) {
            affected_chromosomes.add(gene_chr_map[gene])
          }
        }
      }
      if (affected_chromosomes.size) {
        genes_box_div.innerHTML += "<i><font size=1> Affected chromosomes:</font></i>"
        appendHalfLineBreak(genes_box_div)
        genes_box_div.append(createChromosomeTable(Array.from(affected_chromosomes)))
        genes_box_div.innerHTML += "<hr>"
        appendLineBreak(genes_box_div)
      }
    }

    tree_info_div.appendChild(genes_box_div)
    header_genes.appendChild(createExpandBox(genes_box_id))
    var info_icon = createInfoTooltip("fa fa-question-circle")
    info_text = "List of gene mutations shared between matching subclones. " +
        "The matching subclones are indicated by the corresponding node color."
    addInfoBoxToElement(info_icon, info_text, bg_color="#353935", width=170, margin_left="", position="bottom", line_height="20px")
    header_genes.appendChild(info_icon)
    appendLineBreak(tree_info_div)
  }

  // Metadata.
  var header_meta = createInfoHeader("<b>Cluster metadata</b>", color_motif=background_color)
  header_meta.style.direction = "ltr"
  tree_info_div.appendChild(header_meta)

  var meta_box_id = "cluster_meta"
  var meta_box_div = createDivContainer(meta_box_id)
  meta_box_div.style.direction = "ltr"

  div_table_1 = createDivContainer("table_1")
  div_table_1.style.display = "inline-block"
  div_table_1.appendChild(table_color_codes)
  meta_box_div.appendChild(div_table_1)
  appendLineBreak(meta_box_div)
  appendLineBreak(meta_box_div)

  var metadata_table = getMetadataTable(cluster_metadata)
  metadata_table.style.display = "inline-table"
  meta_box_div.appendChild(metadata_table)
  appendLineBreak(meta_box_div)

  tree_info_div.appendChild(meta_box_div)
  header_meta.appendChild(createExpandBox(meta_box_id, reverse=true))
  var info_icon = createInfoTooltip("fa fa-question-circle")
  info_text = "Summary of the clinical data and the corresponding color legend for the samples in the selected cluster." 
  addInfoBoxToElement(info_icon, info_text, bg_color="#353935", width=170, margin_left="", position="bottom", line_height="20px")
  header_meta.appendChild(info_icon)

  appendLineBreak(tree_info_div)
}

function  getRandomColor(alpha=0) {
  var letters = '0123456789ABCDEF';
  var color = '#';
  for (var i = 0; i < 6; i++) {
    color += letters[Math.floor(Math.random() * 16)];
  }

  const fg = tinycolor(color).setAlpha(alpha)
  const bg = tinycolor("rgb (255, 255, 255)")
  const a = alpha
  const final_color = tinycolor({
      r: (1-a)*fg._r+a*bg._r,
      g: (1-a)*fg._g+a*bg._g,
      b: (1-a)*fg._b+a*bg._b 
  })
  return final_color.toHexString()
} 

class ColorGenerator {
  constructor(default_colors, shuffle=false, transparent_flavour=0) {
    this.color_idx = -1
    this.default_colors = [...default_colors]
    this.transparent_flavour = transparent_flavour
    if (shuffle) {
      shuffleArray(this.default_colors)
    }
  }

  next() {
    this.color_idx++
    if (this.color_idx < this.default_colors.length) { 
      return this.default_colors[this.color_idx]
    } else {
      return getRandomColor(this.transparent_flavour)
    }
  }
}

function getMetadataColorMap(sample_metadata_map, existing_metadata_color_map) {
  metadata_color_map = deepCopy(existing_metadata_color_map)
  for (const [sample, metadata] of Object.entries(sample_metadata_map)) {
    for (const [key, value] of Object.entries(metadata)) {
      if (!(key in metadata_color_map)) {
        metadata_color_map[key] = {}
      }
      if(!(value in metadata_color_map[key])) {
        metadata_color_map[key][value] = getRandomColor(0.2)
      }
    }
  } 
  let sample_metadata_colors = deepCopy(sample_metadata_map)
  for (const [sample, metadata] of Object.entries(sample_metadata_map)) {
    for (const [key, value] of Object.entries(metadata)) {
      sample_metadata_colors[sample][key] = metadata_color_map[key][value]
    }
  }
  return [sample_metadata_colors, metadata_color_map]
}

function createChromosomeTable(affected_chromosomes) {
  chr_list = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16',
      '17', '18', '19', '20', '21', '22', 'X', 'Y']
  var table = document.createElement("table")
  table.style.border = "1px solid darkgray"
  var tr = document.createElement('TR')
  table.appendChild(tr)
  for (chr of chr_list) {
    var td = document.createElement('TD')
    td.style.border = "1px solid lightgray"
    //td.fontWeight = "bold"
    td.style.fontSize = "11px"
    td.style.paddingRight = "2px"
    td.style.paddingLeft = "2px"
    td.appendChild(document.createTextNode(chr))
    if (affected_chromosomes.includes(chr)) {
      td.style.backgroundColor = "tomato"
    }
    tr.appendChild(td)
  }
  return table
}

function getMetadataTable(metadata_samples) {
  var table = document.createElement("table")
  table.style.cssFloat = "left"
  table.style.margin = "20px"
  table.style.marginTop = "5px"

  first_row = true
  for (var row in metadata_samples){
    if(first_row) {
      first_row = false
      var tr = document.createElement('TR');
      table.appendChild(tr);
      //table.style.tableLayout="fixed"
      var td = document.createElement('TD');
      td.style.fontSize = "12px"
      tr.appendChild(td)
      td.appendChild(document.createTextNode(""))
      for (var col in metadata_samples[row]){
        var td = document.createElement('TD');
        td.style.fontSize = "13px"
        td.style.padding = "5px"
        td.style.transform = "rotate(180deg)"
        /*td.style.msTransform = "rotateY(180deg)"
        td.style.webkitTransform = "rotateY(180deg)"
        td.style.mozTransform = "rotateY(180deg)"
        td.style.oTransform = "rotateY(180deg)"*/
        td.style.writingMode = "vertical-rl"
        td.style.whiteSpace = "nowrap" 
        td.appendChild(document.createTextNode(col));
        tr.appendChild(td)
      }
    }
    var tr = document.createElement('TR');
    table.appendChild(tr);
    var td = document.createElement('TD');
    td.style.fontSize = "13px"
    td.style.paddingRight = "2px"
    td.style.whiteSpace = "nowrap"
    tr.appendChild(td)
    td.appendChild(document.createTextNode(row));
    for (var col in metadata_samples[row]){
      var td = document.createElement('TD');
      td.style.backgroundColor=metadata_samples[row][col]
      td.appendChild(document.createTextNode(""));
      tr.appendChild(td)
    }
  }
  // Add 2 empty rows.
  for (var idx of [1,2]) {
    var tr = document.createElement('TR');
    table.appendChild(tr);
    var td = document.createElement('TD');
    td.style.backgroundColor="rgba(0, 0, 0, 0)"
    td.appendChild(document.createTextNode("."));
    td.style.color="rgba(0, 0, 0, 0)"
    tr.appendChild(td)
  }
  return table
}

function getColorCodesTable(map) {
  var table = document.createElement("table")
  table.style.cssFloat = "left"
  table.style.margin = "20px"
  table.style.marginBottom = "0px"

  max_num_columns = 0
  for (var row in map){
    max_num_columns = Math.max(max_num_columns, (map[row].size))
  }

  for (var row in map){
    var tr = document.createElement('TR');
    table.appendChild(tr);
    var td = document.createElement('TD');
    tr.appendChild(td)
    td.appendChild(document.createTextNode(row));
    td.style.fontSize = "13px"
    td.style.paddingRight = "2px"
    for (var col in map[row]){
      var td = document.createElement('TD');
      td.style.backgroundColor=map[row][col]
      td.style.textAlign="center"
      td.style.fontSize = "12px"
      td.style.paddingRight = "2px"
      td.style.paddingLeft = "2px"

      td.style.color="white"
      td.appendChild(document.createTextNode(col));
      tr.appendChild(td)
    }
  }
  return table
}

//////////////////////
//// Heatmap view ////
//////////////////////
function populateHeatmapView(args) {
  async_func(args, populateHeatmapView_slow)
}

function populateHeatmapView_slow(args) {
  // Canvas.
  var tree_cohort_div_id = args.tree_cohort_div_id
  div_tree_cohort = document.getElementById(tree_cohort_div_id)
  div_tree_cohort.innerHTML = ""
  div_tree_cohort.style.zoom = 1

  var tree_info_div_id = args.tree_info_div_id
  div_tree_info = document.getElementById(tree_info_div_id)
  div_tree_info.style.border = "0px"
  div_tree_info.style.backgroundColor = "transparent"
  div_tree_info.innerHTML = "<div style='text-align:center;height:100%;display:flex;flex-direction:row;" +
      "align-items:center;justify-content:center;'><i>Click on the sample names <br/>and" +
      " on the colored frames<br/>(corresponding to the given tree clusters)<br/> to visualize additional information.&lrm;</i></div>"

  // Data.
  data = args.data
  distances = data["pairwise_tree_distances"]
  sample_list = Object.keys(data["trees"])
  clusters = [sample_list]
  if ("clusters" in data) {
    clusters = args.data["clusters"]
    sample_list = []
    for (cluster of clusters) {
      for (sample_name of cluster) {
        sample_list.push(sample_name)
      }
    }
  }
  
  // Visualization.
  var tooltip = d3.select("div#" + tree_cohort_div_id)
    .append("div")
    .style("opacity", 0)
    .attr("class", "tooltip")
    .style("background-color", "white")
    .style("border", "solid")
    .style("border-width", "2px")
    .style("border-radius", "5px")
    .style("padding", "3px")
    .style("visibility", "hidden")
    .style("display", "none")

  var margin = {top: 30, right: 80, bottom: 80, left: 50},
    width = 600 - margin.left - margin.right,
    height = 600 - margin.top - margin.bottom;

  var svg = d3.select("#" + tree_cohort_div_id)
    .append("svg")
      .attr("width", width + margin.left + margin.right)
      .attr("height", height + margin.top + margin.bottom)
    .append("g")
      .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

  // Build X scales and axis.
  var x = d3.scaleBand()
    .range([0, width])
    .domain(sample_list)
    .padding(0)
  svg.append("g")
    .attr("transform", "translate(0," + height + ")")
    .call(d3.axisBottom(x))
    .selectAll("text")
      .attr("transform", "translate(-10,10)rotate(-90)")
      .style("text-anchor", "end")
      .style("font-size", 6)
      .style("cursor", "pointer")
      .on("click", function(d){
        showTreeInfo(d, args)
      })

  // Build Y scales and axis/
  var y = d3.scaleBand()
    .range([0, height])
    .domain(sample_list)
    .padding(0);
  svg.append("g")
    .call(d3.axisLeft(y))
    .selectAll("text")
      .style("font-size", 6)
      .style("cursor", "pointer")
      .on("click", function(d){
        showTreeInfo(d, args)
      })


  vmin = 1
  vmax = 0
  distances.forEach(function(item) {
    if (parseFloat(item.distance) > vmax) {
      vmax = parseFloat(item.distance)
    }
    if (parseFloat(item.distance) < vmin) {
      vmin = parseFloat(item.distance)
    }
  })
   
  // Build color scale
  legend_colors = ['#073B6F', '#0b559f', '#2b7bba', '#539ecd', '#89bedc', '#bad6eb', '#dbe9f6'] // Seaborn Blues
  var myColor = d3.scaleLinear()
    .range(legend_colors)
    .domain(linspace(vmin,vmax, 7))

  // Render heatmap.
  svg.selectAll()
    .data(distances) 
    .enter()
    .append("rect")
      .attr("x", function(d) { return x(d.sample_1) })
      .attr("y", function(d) { return y(d.sample_2) })
      .attr("width", x.bandwidth() )
      .attr("height", y.bandwidth() )
      .style("fill", function(d) { return myColor(d.distance)} )
    .on('mousemove', function(d) {
      tooltip
        .html("&nbsp;Distance b/w " + d.sample_1 + " & " + d.sample_2 + ": " + Math.round(d.distance * 100) / 100 + "&nbsp;") 
        .style("left", (d3.mouse(this)[0]) + "px")
        .style("top", (d3.mouse(this)[1]) + "px")
        .style("visibility", "visible")
        .style("display", "block")
        .style("position", "absolute")
        .style("z-index" ,10)
        .style("font-size", "10px") 
    })
    .on('mouseover', function(d) {
      tooltip.style("opacity", 1)
    })
    .on('mouseleave', function() {
      tooltip.style("opacity", 0)
        .style("visibility", "hidden")
        .style("display", "none") 
    })

  cluster_starts = {}
  clusters.forEach(function (cluster, idx) {
    cluster_starts[cluster[0]] = {
        "cluster_idx": idx,
        "size": cluster.length, 
        "color": data["trees"][cluster[0]]["2d_color"],
        "cluster_bg_color": data["trees"][cluster[0]]["cluster_color"],
        "cluster_metadata": trees[cluster[0]]["sample_metadata_colors"],
        "table_color_codes": trees[cluster[0]]["table_color_codes"],
        "matching_nodes_cluster_details": trees[cluster[0]]["matching_nodes_cluster_details"]
    }
  })

  if (args.rectangle == 0) {
    // Add rectangle for all the samples.
    num_samples = sample_list.length
    svg.selectAll()
      .data(distances)
      .enter()
      .append("rect")
        .attr("x", function(d) {
          if (d.sample_1 == clusters[0][0]  && d.sample_1 == d.sample_2) {
            return x(d.sample_1)
          }
        })
        .attr("y", function(d) {
          if (d.sample_1 == clusters[0][0] && d.sample_1 == d.sample_2) {
            return y(d.sample_1)
          }
        })
        .attr("width", function(d) {
          if (d.sample_1 == clusters[0][0] && d.sample_1 == d.sample_2) {
            return  num_samples * x.bandwidth()
          }
        })
        .attr("height", function(d) {
          if (d.sample_1 == clusters[0][0] && d.sample_1 == d.sample_2) {
            return  num_samples * y.bandwidth()
          }
        })
        .style("fill", "none")
        .style("stroke", function(d) { //"#0b559f"
          if (d.sample_1 == clusters[0][0]){
            return "LightGray"
          }
        })
        .style("stroke-width", "5")
        .style("cursor", "pointer")
        .on("click", function(d){
            args_cluster = {}
            args_cluster["matching_nodes_details"] = new Map()
            args_cluster["tree_info_div_id"] = tree_info_div_id
            args_cluster["highlighted_genes"] = []
            args_cluster["cluster_bg_color"] = "white"
            args_cluster["cluster_metadata"] = trees[clusters[0][0]]["cohort_metadata"]
            args_cluster["table_color_codes"] = trees[clusters[0][0]]["cohort_table_color_codes"]
            showClusterInfo(args_cluster)
        })
        .on('mousemove', function(d) {
          tooltip
            .html("<div align=center>Click to show metadata<br/>for all the samples.</font>")
            .style("left", (d3.mouse(this)[0]) + "px")
            .style("top", (d3.mouse(this)[1]) + "px")
            .style("visibility", "visible")
            .style("display", "block")
            .style("position", "absolute")
            .style("z-index" ,10)
            .style("font-size", "10px")
        })
        .on('mouseover', function(d) {
          tooltip.style("opacity", 1)
        })
        .on('mouseleave', function() {
          tooltip.style("opacity", 0)
            .style("visibility", "hidden")
            .style("display", "none")
        })

    // Add cluster rectangles.
    svg.selectAll()
      .data(distances)
      .enter()
      .append("rect")
        .attr("x", function(d) {
          if (d.sample_1 in cluster_starts && d.sample_1 == d.sample_2 && cluster_starts[d.sample_1]["size"] > 1) {
            return x(d.sample_1)
          }
        })
        .attr("y", function(d) {
          if (d.sample_1 in cluster_starts && d.sample_1 == d.sample_2 && cluster_starts[d.sample_1]["size"] > 1) {
            return y(d.sample_1)
          }
        })
        .attr("width", function(d) {
          if (d.sample_1 in cluster_starts && d.sample_1 == d.sample_2 && cluster_starts[d.sample_1]["size"] > 1) {
            return  cluster_starts[d.sample_1]["size"] * x.bandwidth()
          }
        })
        .attr("height", function(d) {
          if (d.sample_1 in cluster_starts && d.sample_1 == d.sample_2 && cluster_starts[d.sample_1]["size"] > 1) {
            return  cluster_starts[d.sample_1]["size"] * y.bandwidth()
          }
        })
        .style("fill", "none")
        .style("stroke", function(d) { //"#0b559f"
          if (d.sample_1 in cluster_starts){
            return cluster_starts[d.sample_1]["color"]
          }
        })
        .style("stroke-width", "3")
        .style("cursor", "pointer")
        .on("click", function(d){
            args_cluster = {}
            args_cluster["matching_nodes_details"] = cluster_starts[d.sample_1]["matching_nodes_cluster_details"] 
            args_cluster["tree_info_div_id"] = tree_info_div_id
            args_cluster["highlighted_genes"] = data["highlighted_genes"] 
            args_cluster["cluster_bg_color"] = cluster_starts[d.sample_1]["cluster_bg_color"]
            args_cluster["cluster_metadata"] = cluster_starts[d.sample_1]["cluster_metadata"]
            args_cluster["table_color_codes"] = cluster_starts[d.sample_1]["table_color_codes"];
            showClusterInfo(args_cluster)
        }) 
        .on('mousemove', function(d) {
          tooltip
            .html("&nbsp;Click to show details of cluster " + cluster_starts[d.sample_1]["cluster_idx"] + ".")
            .style("left", (d3.mouse(this)[0]) + "px")
            .style("top", (d3.mouse(this)[1]) + "px")
            .style("visibility", "visible")
            .style("display", "block")
            .style("position", "absolute")
            .style("z-index" ,10)
            .style("font-size", "10px")
        })
        .on('mouseover', function(d) {
          tooltip.style("opacity", 1)
        })
        .on('mouseleave', function() {
          tooltip.style("opacity", 0)
            .style("visibility", "hidden")
            .style("display", "none")
        })
  }

  // Color legend
  var grad = svg.append('defs')
    .append('linearGradient')
    .attr('id', 'grad')
    .attr('x1', '0%')
    .attr('x2', '0%')
    .attr('y1', '0%')
    .attr('y2', '100%');

  grad.selectAll('stop')
    .data(legend_colors)
    .enter()
    .append('stop')
    .style('stop-color', function(d){ return d; })
    .attr('offset', function(d,i){
      return 100 * (i / (legend_colors.length - 1)) + '%';
    })

  svg.append('rect')
    .attr('x', width + 10)
    .attr('y', 0)
    .attr('width', 25)
    .attr('height', 150)
    .style('fill', 'url(#grad)');
  
  svg.append("text")
    .attr("id", "rectangleText")
    .attr("class", "visible")
    .attr("x", width + 12)
    .attr("y", -5) 
    .attr("width",100)
    .style("font-size", "10px")
    .text("max");

  svg.append("text")
    .attr("id", "rectangleText")
    .attr("class", "visible")
    .attr("x", width + 13)
    .attr("y", 163) 
    .attr("width",100)
    .style("font-size", "10px")
    .text("min");

  svg.append("text")
    .attr("id", "rectangleText")
    .attr("class", "visible")
    .attr("x", -37)
    .attr("y", width + 36)
    .attr("transform", "translate(-10, 10) rotate(-90)")
    .style("text-anchor", "end")
    .attr("width",100)
    .style("font-size", "15px")
    .attr("fill", "white")
    .text("similarity");
}

/////////////////
//// 2D view ////
/////////////////
function MDS(distances, dimensions) {
  // Code taken from https://github.com/benfred/mds.js
  dimensions = dimensions || 2;

  // square distances
  var M = numeric.mul(-0.5, numeric.pow(distances, 2));

  // double centre the rows/columns
  function mean(A) { return numeric.div(numeric.add.apply(null, A), A.length); }
    var rowMeans = mean(M),
    colMeans = mean(numeric.transpose(M)),
    totalMean = mean(rowMeans);

    for (var i = 0; i < M.length; ++i) {
      for (var j =0; j < M[0].length; ++j) {
        M[i][j] += totalMean - rowMeans[i] - colMeans[j];
      }
    }

    // take the SVD of the double centred matrix, and return the
    // points from it
    var ret = numeric.svd(M)
    var eigenValues = numeric.sqrt(ret.S)
    return ret.U.map(function(row) {
      return numeric.mul(row, eigenValues).splice(0, dimensions);
    });
}

function populate2DView(args) {
  async_func(args, populate2DView_slow)
} 

function populate2DView_slow(args) {
  // Canvas.
  var tree_cohort_div_id = args.tree_cohort_div_id
  div_tree_cohort = document.getElementById(tree_cohort_div_id)
  div_tree_cohort.innerHTML = ""
  div_tree_cohort.style.zoom = 1

  outer_div = createDivContainer("outer_div_2d")
  outer_div.style.border = "solid"
  outer_div.style.borderWidth = "2px"
  outer_div.style.border = "solid"
  outer_div.style.borderColor = "darkgray"
  outer_div.style.borderRadius = "5px"
  outer_div.style.padding = "3px"
  outer_div.style.backgroundColor = "white"
  outer_div.style.display = "none"
  outer_div.style.zIndex = 10
  div_tree_cohort.appendChild(outer_div)

  var tree_info_div_id = args.tree_info_div_id
  div_tree_info = document.getElementById(tree_info_div_id)
  div_tree_info.style.border = "0px"
  div_tree_info.style.backgroundColor = "transparent"
  div_tree_info.innerHTML = "<div style='text-align:center;height:100%;display:flex;flex-direction:row;" +
      "align-items:center;justify-content:center;'><i>Click on the 2D points" +
      "<br/>to visualize additional information<br/>about the corresponding mutation tree.&lrm;</i></div>"

  // Data.
  distances = args.data["pairwise_tree_distances"]
  var sample_id_map = {}
  distances.forEach(function(item) {
    sample = item["sample_1"]
    if (!(sample in sample_id_map)) {
      sample_id_map[sample] = mapSize(sample_id_map)
    }
  })

  var samples = Object.keys(sample_id_map)
  var num_samples = samples.length
  var id_sample_map = {}
  for (var sample in sample_id_map) {
    id_sample_map[sample_id_map[sample]] = sample
  }

  var distance_matrix = Array(num_samples).fill(null).map(() => Array(num_samples).fill(0));
  distances.forEach(function(item) {
    i = sample_id_map[item["sample_1"]]
    j = sample_id_map[item["sample_2"]]
    distance_matrix[i][j] = parseFloat(item.distance)
  })

  var points = MDS(distance_matrix)
  var point_map = points.map(function (item, idx) {
    return { 
      "x": parseFloat(item[0]), 
      "y": parseFloat(item[1]), 
      "color": data["trees"][samples[idx]]["2d_color"],
      "sample_name": samples[idx],
      "tree_json":args.data["trees"][samples[idx]]["tree"] 
    } 
  });

  // Render the plot.
  var margin = {top: 30, right: 80, bottom: 80, left: 50},
    width = 600 - margin.left - margin.right,
    height = 600 - margin.top - margin.bottom;

  var svg = d3.select("#" + tree_cohort_div_id)
  .append("svg")
    .attr("width", width + margin.left + margin.right)
    .attr("height", height + margin.top + margin.bottom)
  .append("g")
    .attr("transform",
          "translate(" + margin.left + "," + margin.top + ")");

  var minX = 10000000; maxX = -10000000;
  var minY = 10000000; maxY = -10000000;
  point_map.forEach(function(item) {
    maxX = Math.max(maxX, item.x)
    maxY = Math.max(maxY, item.y)
    minX = Math.min(minX, item.x)
    minY = Math.min(minY, item.y)
  })
  minX = minX - (maxX - minX) / 10
  maxX = maxX + (maxX - minX) / 10
  minY = minY - (maxY - minY) / 10
  maxY = maxY + (maxY - minY) / 10

  // Add X axis
  var x = d3.scaleLinear()
    .domain([minX, maxX])
    .range([ 0, width ]);
  svg.append("g")
    .attr("transform", "translate(0," + height + ")")
    .call(d3.axisBottom(x));

  // Add Y axis
  var y = d3.scaleLinear()
    .domain([minY, maxY])
    .range([ height, 0]);
  svg.append("g")
    .call(d3.axisLeft(y));

  // Add dots
  svg.append('g')
    .selectAll("dot")
    .data(point_map)
    .enter()
    .append("circle")
      .attr("cx", function (d) { return x(d.x); } )
      .attr("cy", function (d) { return y(d.y); } )
      .attr("r", 6)
      .attr('fill-opacity', 0.5)
      .style("fill", function (d) { return d.color; } )
      .style("stroke", "lightgray")
      .style("cursor", "pointer")
      .on("click", function(d){
        showTreeInfo(d.sample_name, args)
      })
    .on('mouseover', function(d) {
      outer_div = document.getElementById("outer_div_2d")
      outer_div.innerHTML = ""
      outer_div.style.position = "absolute"
      outer_div.style.display = "inline-block"

      text_div = createDivContainer("text_div")
      outer_div.appendChild(text_div)
      text_div.innerHTML = d.sample_name + "&nbsp;&nbsp;"
      text_div.style.fontSize = "10px"
      text_div.style.transform = "rotate(180deg)"
      text_div.style.writingMode = "vertical-lr"
      text_div.style.verticalAlign = "top"
      text_div.style.textAlign = "right"
      text_div.style.position = "relative"
      text_div.style.display = "inline-block"

      tree_div_id = "tree_div_2d"
      tree_div = createDivContainer(tree_div_id)
      outer_div.appendChild(tree_div)
      displayTree(tree_div_id, "", d.tree_json, "", "", null, tree_info_view=true)
      tree_div.style.zoom = 0.2
      tree_div.style.display = "inline-block"
      tree_div.style.position = "relative"

      outer_div.style.top = (d3.mouse(this)[1]) + "px"
      outer_div.style.left = (d3.mouse(this)[0] + outer_div.offsetWidth) + "px"
    })
    .on('mouseleave', function() {
      div = document.getElementById("outer_div_2d")
      div.style.display = "none"
    })
}

////////// Help functions display tree ///////////
function circle_ray(percentage) {
  if(percentage == null) {
    return 10
  }
  if(percentage == 0) {
    return 5
  }
  return Math.max(8, Math.min(30, percentage * 60))
}

function get_cn_event_for_gene(node, target_gene) {
  if (target_gene in node.data.gene_events) {
    if ("CNA" in node.data.gene_events[target_gene]) {
      return node.data.gene_events[target_gene]["CNA"]
    }
  }
  return 0;
}

function get_color_for_cn_event(cn_event) {
  if (cn_event > 0 || cn_event == "+") {
    return "tomato";
  } else if (cn_event < 0 || cn_event == "-" || cn_event == "loh") {
    return "lightsteelblue";
  }
  return "darkgray";
}

var violet_color = "#b4a7d6"
function get_node_color(node, target_gene) {
  if (target_gene in node.data.gene_states) {
    gene_events = node.data.gene_states[target_gene]
    if ("CNA" in gene_events) {
      gene_state = gene_events["CNA"]
      return get_color_for_cn_event(gene_state)
    } else {
      return violet_color
    }
  }
  return "";
}

////////// Display tree ///////////
var TREE_CANVAS_WIDTH = 200
var TREE_CANVAS_PADDING = 20
var TREE_CANVAS_HEIGHT = 400
var DEFAULT_ZOOM = 1
function displayTree(div_id, tree_label, tree_data, target_gene, 
    target_drug, drug_gene_map, tree_info_view=true, 
    width=TREE_CANVAS_WIDTH, height=TREE_CANVAS_HEIGHT) {

  // Node info box.
  var tooltip = d3.select("div#" + div_id)
    .append("div")
    .attr("id", "tooltip")
    .style("opacity", 1)
    .attr("class", "tooltip")
    .style("background-color", "white")
    .style("border", "solid")
    .style("border-color", "darkgray")
    .style("border-width", "2px")
    .style("border-radius", "5px")
    .style("padding", "5px")
    .style("white-space", "nowrap")
    .style("overflowY", "scroll")
    .style("visibility", "hidden")
    .style("display", "none")

  // Set the dimensions and margins of the diagram.
  var margin = { top: 30, right: 20, bottom: 35, left: 20 };
  height = height - margin.top - margin.bottom;

  var treemap = d3.tree().size([width, height]);
  var nodes = d3.hierarchy(tree_data);
  nodes = treemap(nodes);
  var svg = d3
    .select("div#" + div_id)
    .append("svg")

  svg.attr("width", width + margin.left + margin.right)
    .attr("height", height + margin.top + margin.bottom)
    .style("display", "block")
    .style("margin", "auto")

  var g = svg.append("g")
          .attr("transform",
            "translate(" + margin.left + "," + margin.top + ")");

  // Add sample name.
  if(!tree_info_view) {
    g.append("text")
      .attr("id", "rectangleText")
      .attr("class", "visible")
      .attr("x", -10)
      .attr("y", -10)
      .attr("width",10)
      .style("font-size", "13px")
      .text(tree_label);
  }

  // Adds the links between the nodes.
  var link = g
    .selectAll(".link")
    .data(nodes.descendants().slice(1))
    .enter()
    .append("g")
    .attr("class", "link");
 
  // Link color.
  link
    .append("path")
    .attr("id", function(d) {
      return tree_label + "_" + d.data.node_id
    })
    .attr("d", function(d) {
      return (
        "M" + d.x + "," + d.y + 
        "C" + d.x + "," + (d.y + d.parent.y) / 2 + " " + d.parent.x + "," + 
        (d.y + d.parent.y) / 2 + " " + d.parent.x + "," + d.parent.y);
    })
    .attr("orient", "auto")
    .style("stroke", function(d) {
      if(d.data.conserved_parent_node_branch && !tree_info_view) { 
        return "#800000"//"dimgray"
      } 
      if (d.parent != null && d.data.gene_events) {
        cn_event = get_cn_event_for_gene(d, target_gene);
        return get_color_for_cn_event(cn_event);
      }
      return "darkgray"
    })

  // Add nodes.
  var node = g
    .selectAll(".node")
    .data(nodes.descendants())
    .enter()
    .append("g")
    .attr("class", function(d) {
      return "node" + (d.children ? " node--internal" : " node--leaf");
    })
    .attr("transform", function(d) {
      return "translate(" + d.x + "," + d.y + ")";
    })

  map_y_coords = {}
  node.each(function(d) {
    if (d.parent) { // not root
      coord_y = Math.round(d.y)
      if (coord_y in map_y_coords) {
        map_y_coords[Math.round(d.y)].push(d)
      } else {
        map_y_coords[Math.round(d.y)] = [d]
      }
    }    
  })
  for (coord in map_y_coords) {
    map_y_coords[coord].sort(function(d_1, d_2) {
      return d_1.x - d_2.x
    });
  }
  map_closest_nodes = {} 
  for (y_coord in map_y_coords) {
    x_coords = map_y_coords[y_coord]
    for (let i=0; i<x_coords.length; i++) {
      n = x_coords[i]
      node_id = n.data.node_id
      map_closest_nodes[node_id] = {"left":null, "right":null, "position":"right"}
      if (x_coords.length > 1) {
        if (i==0) {
          map_closest_nodes[node_id]["position"] = "left"
        } else if (i == x_coords.length - 1) {
          map_closest_nodes[node_id]["position"] = "right"
        } else {
          map_closest_nodes[node_id]["right"] = x_coords[i+1]
          map_closest_nodes[node_id]["left"] = x_coords[i-1]
          dist_left = n.x - x_coords[i-1].x
          dist_right = x_coords[i+1].x - n.x
          if (dist_left <= dist_right) {
            map_closest_nodes[node_id]["position"] = "right"
          } else {
            map_closest_nodes[node_id]["position"] = "left"
          }
        }
      } 
    }
  }
  
  // Link text labels.
  if (tree_info_view) {
    link
      .append("text")
      //.attr("x", -45)
      .attr("text-anchor", "middle")
      .attr("stroke-width", "1px")
      .attr("stroke", "gray")
      .attr("fill", "gray")
      .attr("font-size", "16px")
      .attr("transform", function(d) {
        if (d.parent) {
          var path_element = d3.select('#' + tree_label + "_" + d.data.node_id).node()
          var branch_length = path_element.getTotalLength()
          var coord = path_element.getPointAtLength(branch_length/2)
          var coord_x = coord.x// - d.parent.x
          var coord_y = coord.y// - d.parent.y 
          if(d.x < d.parent.x) {
            coord_x -= 15
          } else {
            coord_x += 15
          }
          return ("translate(" + coord_x + "," + coord_y + ")")
        }
      })
      .text(function(d) {
        if (d.parent && d.data.gene_events) {
          cn_event = get_cn_event_for_gene(d, target_gene);
          if (cn_event) {
            string = "";
            if (cn_event > 0) {
              string = "+";
            }
            string += cn_event;
            return string;
          }
        }
      })
  } else {
    link
      .append("text")
      .attr("stroke-width", "0.25px")
      .attr("font-size", "11px")
      .attr("stroke", "silver")
      .attr("fill", "silver") 
      .attr("transform", function(d) {
        if (d.parent) {
          x_coord = d.parent.x 
          circle_ray_parent = circle_ray(d.parent.data.size_percent)
          y_coord = d.parent.y 
          return ("translate(" + x_coord + "," + y_coord + ")")
        }
      })
      .each(function(d) {
        if (d.parent && d.data.gene_event_categories) {
          node_id = d.data.node_id
          var text_element = d3.select(this)
          var path_element = d3.select('#' + tree_label + "_" + d.data.node_id).node()
          branch_length = path_element.getTotalLength()
          ray_parent = circle_ray(d.parent.data.size_percent)
          ray_node = circle_ray(d.data.size_percent)
          line_height = 13
          bottom_padding = 3
          num_max_lines = (branch_length - ray_node - ray_parent - bottom_padding - 6) / line_height

          // Get path coordinates for neighbouring node.
          node_neighbor = map_closest_nodes[node_id]["right"]
          if (map_closest_nodes[node_id]["position"] == "left") {
            node_neighbor = map_closest_nodes[node_id]["left"]
          }
          node_neighbor_path = null 
          if(node_neighbor) { 
            node_neighbor_path = d3.select('#' + tree_label + "_" + node_neighbor.data.node_id).node()
            ray_node_neighbor = circle_ray(node_neighbor.data.size_percent)
          }

          var idx = 0
          for (let [key, gene_set] of d.data.gene_event_categories.entries()) {
            straight_height = ray_node + idx * line_height + 3
            coord = path_element.getPointAtLength(straight_height)
            coord_x = coord.x - d.parent.x
            coord_y = coord.y - d.parent.y - bottom_padding

            // Display text if the spance is at least 25px.
            if(node_neighbor_path) {
              coord_node_neighbor = node_neighbor_path.getPointAtLength(Math.max(circle_ray(node_neighbor.data.size_percent), straight_height))
              if (Math.abs(coord_node_neighbor.x - coord.x) - ray_node_neighbor < 25) {
                continue
              }
            }

            text = ""
            short_key = key.toString().substring(0,3).toLowerCase()
            num_genes = gene_set.size
            if (num_genes > 1) {
              text = text + num_genes + " " + short_key
            } else {
              text = [...gene_set][0]
            }

            if (map_closest_nodes[node_id]["position"] == "left") {
              coord_x -= text.length * 8
            } else { // right
              coord_x += 7 
            }

            text_element.append("tspan")
              .attr('x', coord_x)
              .attr('y', coord_y)
              .text(text)

            idx++
            if (idx == num_max_lines) {
              break
            }         
          }
        }
      })
  }

  if (tree_info_view) {
    max_genes_in_list = 4
    node.on('mousemove', function(d) {
      // Compute html info.
      if (d.parent || d.data.gene_event_categories) {
        html_info = ""
        if (d.data.matching_label) {
          html_info += "<b>Matching label:</b> " + d.data.matching_label + "<br/>"
        }
        html_info += "<b>Node id:</b> " + d.data.node_id + "<br/>"
        if (d.data.size_percent != undefined) {
          html_info += "<b>Clone size:</b> " + round(d.data.size_percent * 100, 1) + "%<br/>"
        }
        if (d.data.is_neutral) {
          html_info += "<b>Is neutral:</b> " + d.data.is_neutral + "<br/>"
        }
        var gene_events = d.data.gene_event_categories
        var gene_events_with_details = d.data.gene_event_categories_details
        if (gene_events && gene_events.size) {
          html_info += "<p style='margin:6px;'></p><b>Gene events:</b><br/>"
          var keys = Array.from(gene_events.keys()).sort().reverse()
          for (event of keys) {
            var gene_list = Array.from(gene_events_with_details.get(event)).sort()
            html_info += capitalizeFirstLetter(event) + ": " + gene_list.slice(0,max_genes_in_list).join(', ')
            if (gene_list.length > max_genes_in_list) {
              html_info += "..."
            }
            html_info += "<br/>"
          }
        }
      } else if (!(d.parent)){
        html_info = "root"
      }

      // Populate tooltip.
      tooltip 
        .html(html_info)
        .style('width', 99+'%')
        .style("left", (d3.mouse(this)[0]) + "px")
        .style("top", (d3.mouse(this)[1] + 50) + "px")
        .style("visibility", "visible")
        .style("display", "block")
        .style("position", "absolute")
        .style("z-index" ,10)
    })
    .on('mouseover', function(d) {
      tooltip.style("opacity", 1)
    })
    .on('mouseleave', function() {
      tooltip.style("opacity", 0)
        .style("visibility", "hidden")
        .style("display", "none")
    })
    .style("cursor", "pointer")
  } else {
    node.on("mouseover", function(d) {
      d3.select(this)
        .append("text")
        .classed("info", true)
        .attr("x", 0)
        .attr("y", 3)
        .style("font-size", "10px")
        .attr("text-anchor", "middle")
        .attr("fill", "black")
        .attr("stroke", "white")
        .attr("stroke-width", 0.1)
        .text(function(d) {
          if (d.parent && d.data.matching_label && d.data.display_node_label && !tree_info_view) {
            label = (d.data.matching_label).toString()
            // Count 5.5px per digit.
            if (2 * circle_ray(d.data.size_percent) < label.length * 5.5) {
              return label
            }
          }
        })
    })
    .on("mouseout", function() {
      // Remove the info text on mouse out.
      d3.select(this)
        .select("text.info")
        .remove();
    })
  }

  // adds circle of certain size to the node
  node
    .append("circle")
    .attr("r", function(d) {
      if (d.parent == null) { // root
        return 3
      }
      return circle_ray(d.data.size_percent);
    })
    .style("fill", function(d) {
      if (d.data.is_neutral) {
        return "lightyellow"
      }
      if(tree_info_view) {
        if (d.parent != null && d.data.gene_events) {
          return get_node_color(d, target_gene);
        }
      } else {
        return d.data.color
      }
      return "white"
    })
    .style("stroke", function(d) {
      if (d.data.size_percent == 0) {
        return "lightgray"
      }
      //return "gray"
    })
    .on("mouseover", function(d) {
      d3.select(this)
        .append("text")
        .attr("x", 0)
        .attr("y", 3)
        .style("font-size", "10px")
        .attr("text-anchor", "middle")
        .attr("fill", "black")
        .attr("stroke", "white")
        .attr("stroke-width", 0.1)
        .text(function(d) {
          if (d.parent && d.data.matching_label && d.data.display_node_label && !tree_info_view) {
            label = (d.data.matching_label*1000).toString()
            // Count 5.5px per digit.
            if (2 * circle_ray(d.data.size_percent) < label.length * 5.5) {
              return label
            }
          }
        })
    })

  node.append("text")
    .attr("x", 0)
    .attr("y", 3)
    .style("font-size", "10px")
    .attr("text-anchor", "middle")
    .attr("fill", "black")
    .attr("stroke", "white")
    .attr("stroke-width", 0.1)
    .text(function(d) {
      if (d.parent && d.data.matching_label && d.data.display_node_label && !tree_info_view) {
        label = (d.data.matching_label).toString()
        // Count 5.5px per digit.
        if (2 * circle_ray(d.data.size_percent) >= label.length * 5.5) {
          return label
        }
      }
      return ""
    })

  // Drug interaction rectangle. 
  if (tree_info_view) {
    node.append("rect")
    .attr("x", function(d) {
      if (d.parent != null) {
        ray = circle_ray(d.data.size_percent)
        rect_size = ray 
        return -ray/2 - (rect_size - ray) / 2
      }
    })
    .attr("y", function(d) {
      if (d.parent != null) {
        ray = circle_ray(d.data.size_percent)
        rect_size = ray 
        return -ray/2 - (rect_size - ray) / 2
      }
    })
    .attr('width', function(d) {
      if (d.parent != null) {
        ray = circle_ray(d.data.size_percent)
        return ray
      }
    })
    .attr('height', function(d) {
      if (d.parent != null) {
        ray = circle_ray(d.data.size_percent)
        return ray 
      }
    })
    .style("fill", function(d) {
      if (d.parent != null) {
        if(d.data.is_neutral) {
          return "lightyellow"
        }
        if(!d.data.gene_events) {
          return "white"
        }
        if (drug_gene_map && drug_gene_map.has(target_drug)) {
          gene_list = drug_gene_map.get(target_drug)
          for (var gene_1 of gene_list){
            for (var gene_2 of Object.keys(d.data.gene_states)){
              if (gene_2.includes(gene_1)) {
                return "green"
              }
            }
          }
        }
        var color = get_node_color(d, target_gene)
        if(!color) { 
          color = "white"
        } 
        return color
      }
    });
  }    
}

////////// Display matching tree pairs ///////////
function getYShift(depth) {
  return depth * 40 + 60
}

function getXShift(dx, group) {
  if (group == 2) {
    return dx * 100 + 150
  }
  return dx * 100 
}

async function async_displayTreeMatching(div_id, data, idx="") {
  displayTreeMatching(div_id, data)
  await sleep(2000)
}

function displayTreeMatching(div_id, data) {
  var colors = ['#a6cee3', '#1f78b4', '#b2df8a', '#33a02c', '#fb9a99', '#e31a1c', 
                '#fdbf6f', '#ff7f00', '#cab2d6', '#6a3d9a', '#ffff99', '#b15928'] //d3.schemePaired
  var color_index = 0
  var node_colors = {}

  sample_1 = data["sample_1"] 
  sample_2 = data["sample_2"]
  pair_id = sample_1 + "-" + sample_2
  nodes = data["nodes"]
  links = data["links"]
  score = data["similarity"]  
  max_depth = data["max_depth"]
  
  var margin = { top: 0, right: 0, bottom: 0, left: 0};
  var width = 500
  var height = max_depth * 40 + 100

  var svg = d3.select("div#" + div_id)
    .append("svg")
    .attr("id", pair_id)
    .attr("width", width + margin.left + margin.right)
    .attr("height", height + margin.top + margin.bottom)
    .attr("transform", "translate(" + margin.left + "," + margin.top + ")")

  svg.append("text")
    .attr("id", "rectangleText")
    .attr("class", "visible")
    .attr("x", 20)
    .attr("y", 30)
    .attr("width",10)
    .style("font-size", "13px")
    .text("(" + idx + ") " + sample_1 + " - " + sample_2);

  for (var node of nodes) {
    group = 2
    if (node.id.startsWith(sample_1)) {
      group = 1
    }
    node.x = getXShift(node.dx, group)
    node.y = getYShift(node.depth)
  }

  // Fetch the links between nodes.
  var link = svg.selectAll(".link")
    .data(links)
    .enter()
    .append("path")
    .attr("class", "link")
    .style("stroke", function(d){
      if(d.similarity > 0) {
        link_color = colors[color_index]
        node_colors[d.source] = link_color
        node_colors[d.target] = link_color
        color_index += 1
        return link_color
      }
      return "darkgray"
    })
    .style("stroke-width", "5px")
    .attr('marker-end','url(#arrowhead)')
 
  var edgepaths = svg.selectAll(".edgepath")
  .data(links)
  .enter()
  .append('path')
  .attrs({
    'class': 'edgepath',
    'fill-opacity': 0,
    'id': function (d, i) {
      return 'edgepath_' + pair_id + "_" + i}
    })
  .style("pointer-events", "none");

  var edgelabels = svg.append('g').selectAll(".edgelabel")
    .data(links)
    .enter()
    .append('text')
    .style("pointer-events", "none")
    .attrs({
      'class': 'edgelabel',
      'id': function (d, i) {return 'edgelabel_' + pair_id + "_" + i},
      'font-size': 10,
      'fill': '#aaa'
    });

  edgelabels.append('textPath')
    .attr('xlink:href', function (d, i) {return '#edgepath_' + pair_id + "_" + i})
    .style("text-anchor", "middle")
    .style("pointer-events", "none")
    .attr("startOffset", "50%")
    //.attr("shapePadding", "5px")
    .attr("font-size", 16)
    .text(function (d) {
      if(d.similarity >= 0) {
        //return "sim = " + d.similarity   
      }
    })

  var node = svg.selectAll(".node")
    .data(nodes)
    .enter()
    .append("g")
    .attr("class", "node " + "node_" + pair_id)

  node.append("circle")
    .attr("r", function(d) {
      if (d.is_root) {
        return 3;
      } 
      return 8
  })
  .style("fill", function(d) {
     return node_colors[d.id]
  })

  node.append("text")
    .attr("dx", -4)
    .attr("dy", 4)
    .text(function (d) {
      return ""
    })

  var simulation = d3.forceSimulation()
    .force("link", d3.forceLink().id(function (d) {return d.id;}).strength(0))

  // add the nodes to the simulation
  simulation
    .nodes(nodes)
    .on("tick",function(){
        node
          .attr("transform", function (d) {
            return "translate(" + d.x + ", " + d.y  + ")"
          })
        link.attr("d", getPath)
        edgepaths.attr("d", getPath)
    }) 

  // add the links to the simulation
  simulation
    .force("link").links(links)
  
}

// links are drawn as curved paths between nodes,
// through the intermediate nodes
function curvedPath(d) {
  source_dx = getXShift(d.source.dx, 1)
  target_dx = getXShift(d.target.dx, 2)
  source_dy = getYShift(d.source.depth)
  target_dy = getYShift(d.target.depth)

  var offset = 30;
  var midpoint_x = (source_dx + target_dx) / 2;
  var midpoint_y = (source_dy + target_dy) / 2;

  var dx = (target_dx - source_dx);
  var dy = (target_dy - source_dy);

  var normalise = 70 //Math.sqrt((dx * dx) + (dy * dy));

  var offSetX = midpoint_x + offset*(dy/normalise);
  var offSetY = midpoint_y + offset*(dx/normalise);

  return "M" + d.source.x + "," + d.source.y +
      "S" + offSetX + "," + offSetY +
      " " + d.target.x + "," + d.target.y;
}

function getPath(d) {
  if (d.similarity >= 0){
    return curvedPath(d)
  }
  else {
    return 'M ' + d.source.x + ' ' + d.source.y + ' L ' + d.target.x + ' ' + d.target.y
  }
}

