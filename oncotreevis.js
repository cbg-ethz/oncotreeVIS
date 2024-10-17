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
    target_div.style.zoom = 1
  })
  button_zoomReset.target_div_id = tree_cohort_div_id
  outer_div.appendChild(button_zoomReset)

  // Trev view buttons.
  var tree_view_div = createBlueBorder() 

  var tree_view_button = document.createElement('button')
  tree_view_button.className = "button-15"
  tree_view_button.addEventListener('click', ()=>{ populateTreeView(args); })
  tree_view_button.innerHTML = '<i class="fa fa-tree" style="font-size:19px"></i>&nbsp;TREE VIEW'
  addInfoBoxToElement(tree_view_button, "Show mutation trees side by side, grouped by a given clustering (by default). " + 
      "Nodes correspond to clones. Each node (also shown on the incoming edges) is labeled with the set of provided gene mutations (SNVs, CNAs, etc) acquired by the subclone. "+
      "Matching clones and conserved edges are highlighted. Neutral clones (if specified) are colored in " +
      "<font color=lightyellow><b>lightyellow</b></font>.", bg_color="#0868d2", width=360, margin_left=135, position="top", line_height="13px")
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
  addInfoBoxToElement(button_matching, "Highlight <b>matching clones and conserved branches</b> (default), " +
      "<b>conserved edges only</b>, or <b>matching clones only</b>. Matching nodes have the same " +
      "<i>matching_label</i>.", bg_color="#0868d2", width=165, margin_left="", position="top", line_height="17px")
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
  if ("pairwise_distances" in args.data) { 
    var heatmap_view_button = document.createElement('button')
    heatmap_view_button.className = "button-15"
    heatmap_view_button.addEventListener('click', ()=>{ populateHeatmapView(args) })
    heatmap_view_button.innerHTML = '<i style="font-size:19px" class="fa fa-th-large"></i> HEATMAP VIEW'
    addInfoBoxToElement(heatmap_view_button, "Show a heatmap visualization of the given pairwise distances between the mutation trees.",
        bg_color="#0868d2", width=200, margin_left="", position="top", line_height="17px")
    appendSpace(outer_div) 
    outer_div.appendChild(heatmap_view_button)    
    div_container.appendChild(outer_div) 

    var umap_view_button = document.createElement('button')
    umap_view_button.className = "button-15"
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
  // Tree cohort div.
  tree_cohort_div = createDivContainer(tree_cohort_div_id)
  tree_cohort_div.style.zoom = 1
  tree_cohort_div.style.float = "left"
  tree_cohort_div.style.position = "relative"
  tree_cohort_div.style.width = "73%"
  tree_cohort_div.style.padding = "3px"     
  outer_div.appendChild(tree_cohort_div)
  // Tree info div.
  var tree_info_div_id = args.tree_info_div_id
  tree_info_div = createDivContainer(tree_info_div_id)
  tree_info_div.style.borderRadius = "8px"
  tree_info_div.style.fontSize = "13px"
  tree_info_div.style.position = "fixed"
  tree_info_div.style.width = "24%"
  tree_info_div.style.maxWidth = "430px"
  tree_info_div.style.minWidth = "300px"
  tree_info_div.style.padding = "7px"
  tree_info_div.style.float = "right"
  tree_info_div.style.left = "calc(73% + 2px)"
  tree_info_div.style.bottom = "3px"
  tree_info_div.style.overflowX = "scroll"
  tree_info_div.style.overflowY = "scroll"
  tree_info_div.innerHTML = ""
  outer_div.appendChild(tree_info_div)
  div_container.appendChild(outer_div) 
  var div_tree_cohort_top_offset = tree_cohort_div.getBoundingClientRect().top
  tree_info_div.style.top = div_tree_cohort_top_offset + "px"
}

///////////////////
//// Tree view ////
///////////////////
function populateTreeView(args){
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

  knn = data["matching_trees"] 
  display_text_label = false
  if ("display_text_label" in data && data["display_text_label"]){
    display_text_label = true
  }

  // Prepare tree cohort container.
  var tree_cohort_div_id = args.tree_cohort_div_id
  div_tree_cohort = document.getElementById(tree_cohort_div_id)
  div_tree_cohort.innerHTML = ""
  var tree_info_div_id = args.tree_info_div_id
  tree_info_div = document.getElementById(tree_info_div_id)
  tree_info_div.innerHTML = '<div style="text-align:center;height:100%;display:flex;flex-direction:row;' +
      'align-items:center;justify-content:center;"><i>Click on the trees' +
      '<br/>or on the "show cluster details" icons<br/>to visualize additional information.</i></div>'
  tree_info_div.style.backgroundColor = "white"

  for (const [i, cluster] of clusters.entries()) {
    for (sample_name of cluster) {
      var node_list = getTreeNodes(trees[sample_name]["tree"]) 
      for(node of node_list) {
        node.data.color = "white"
        node.data.conserved_parent_node_branch = false
        node.data.display_node_label = false
        if (args.sorting == 0){ // clusters
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
    //for (const [j, sample_name] of cluster.entries()) {
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
        click_cluster_details_div.style.cursor = "pointer"
        var text_color = tinycolor(cluster_color).darken(30).desaturate(40).toHexString()
        click_cluster_details_div.innerHTML = '&nbsp;<i class="fa fa-desktop fa-sm" style="color:' + text_color +
          '"></i> <i><font color="' + text_color  + '"> &thinsp; show cluster details </i>'

        args_cluster = {}
        args_cluster["matching_nodes_details"] = {} // TODO matching_clones_color_map
        args_cluster["tree_info_div_id"] = tree_info_div_id
        args_cluster["cluster_bg_color"] = cluster_color
        args_cluster["cluster_metadata"] = trees[cluster[0]]["sample_metadata_colors"]
        args_cluster["table_color_codes"] = trees[cluster[0]]["table_color_codes"];
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
      outer_div.appendChild(tree_div)
      displayTree(tree_div_id, sample_name, tree_json, "", "", null, tree_info_view=false)
    })
  }
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

// Read sample_map
function loadFileVariable(filename){
  var script = document.createElement("script");
  script.src = filename;
  document.head.appendChild(script);
  sleep(2000).then(() => {})
}

const objectToMap = obj => {
   const keys = Object.keys(obj);
   const map = new Map();
   for (let i = 0; i < keys.length; i++){
      //inserting new key value pair inside map
      map.set(keys[i], obj[keys[i]]);
   };
   return map;
};

function deepCopy(oldValue) {
  var newValue
  strValue = JSON.stringify(oldValue)
  return newValue = JSON.parse(strValue)
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

function mapSize(map) {
  return Object.keys(map).length
}

function round(value, precision) {
    var multiplier = Math.pow(10, precision || 0);
    return Math.round(value * multiplier) / multiplier;
}

//////////////////////////
//// UTILS D3.JS TREE ////
//////////////////////////
function getTreeNodes(tree_json) {
  var node_hierarchy = d3.hierarchy(tree_json)
  var treemap = d3.tree()
  return treemap(node_hierarchy).descendants()
}

function getNodeMatchingLabels(tree_json){
  label_node_map = {}
  node_list = getTreeNodes(tree_json)
  for (node of node_list) {
    if (node.parent != null && node.data.matching_label && !node.data.is_neutral){ // discard root and neutral nodes.
      label = node.data.matching_label
      if(!(label in label_node_map)) {
        label_node_map[label] = [] 
      }
     label_node_map[label].push(node)
    }
  }
  return label_node_map
}

function getMatchingLabels(tree_json){
  var matching_labels = new Set()
  var node_list = getTreeNodes(tree_json)
  for (node of node_list) {
    if(node.data.matching_label) {
      matching_labels.add(node.data.matching_label)
    }
  }
  return matching_labels
}

// Returns a list of node-parent strings.
function getBranches(tree_json) {
  branches = new Set()
  node_list = getTreeNodes(tree_json)
  for (node of node_list) {
    if(node.parent) {
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

function updateGeneStates(gene_events, parent_gene_states) {
  gene_events_map = objectToMap(gene_events)
  var new_gene_states = new Map(parent_gene_states)
  // Update the new_gene_states with each current gene event.
  for (let [gene, events] of gene_events_map.entries()) {
    events_map = objectToMap(events)
    if (!(new_gene_states.has(gene))) {
      new_gene_states.set(gene, events_map)
    } else {
      for (let [old_event_type, old_value] of new_gene_states.get(gene).entries()) {
        for (let [new_event_type, new_value] of events_map.entries()) {
          // Apply each event from gene_events
          if(old_event_type == new_event_type) {
            if (old_event_type == "CNA") {
              var new_cn_value = CNValueToInt(old_value) + CNValueToInt(new_value)
              if (new_cn_value == 0) {
                events_map.delete("CNA")
                if (events_map.size == 0) {
                  new_gene_states.delete(gene)
                }
              } else {
                events_map.set("CNA", new_cn_value)
              }
            } else {
              new_gene_states.get(gene).set(new_event_type,  old_value + ";" + new_value)
            }
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
      if (!node.parent) { // root
        node.data.gene_states = objectToMap(node.data.gene_events)
      } else {
        node.data.gene_states = updateGeneStates(node.data.gene_events, node.parent.data.gene_states)
      }
    } else {
      node.data.gene_states = new Map()
    }
  }
}

//////////////////////
//// HTML helpers ////
//////////////////////
function appendLineBreak(div) {
  div.appendChild(document.createElement("p"))
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
  if (position == "top") {
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
  div.innerHTML = "&nbsp;<i class='fa fa-question-circle' style='color:onyx'></i>"
  return div
}

////////////////////////////
//// Populate tree info ////
////////////////////////////
function createInfoHeader(html_text){
  var color_motif = "#d2cae6"
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
    } else {
      div.style.display = 'none'
      button.innerHTML = button.innerHTML.replace("up", "down")
    }
    button.show_div = !show_div
  })
  button_expand.div_id = div_id
  button_expand.button_id = button_id
  button_expand.show_div = show_div

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

  // Compute gene list and variables related to the selected target gene. 
  const [gene_categories, gene_categories_with_details] = getGeneCategoriesInTree(tree_object)
  var drug_gene_map = {}
  var drug_name_map = {}
  var drug_list = []
  var target_gene = ""
  var target_drug = ""
  if (mapSize(gene_categories)) {
    var gene_list = Array.from(Object.values(gene_categories)[0])
    var gene_list_with_details = Array.from(Object.values(gene_categories_with_details)[0])
    drug_gene_map = getDrugToGeneMap(gene_list, gene_drug_map)
    drug_list = findDrugsAffectingGeneList(gene_list, drug_gene_map, tree_object)
  }

  // Display sample name.
  header_sample_name = createInfoHeader("<strong><b>" + sample_name + "</b></strong>")
  tree_info_div.appendChild(header_sample_name)
  appendLineBreak(tree_info_div)

  // Gene selection.
  if (mapSize(gene_categories)) {
    gene_selection_dropdown = document.createElement("select")
    gene_selection_dropdown.setAttribute("id", "gene_selection")
    gene_selection_dropdown.classList.add("gene_selection")
    gene_selection_dropdown.style.width = "200px"
    option = document.createElement('option');
    option.textContent = "-- Select target gene event --"
    gene_selection_dropdown.appendChild(option);
    keys = Object.keys(gene_categories).sort()
    for (event of keys) {
      array = gene_categories[event]
      sorted_genes = Array.from(array).sort()
      addGeneToGeneList(gene_selection_dropdown, sorted_genes, event)
    }
    div_container = document.createElement('div')
    div_container.style.width = "fit-content"
    div_container.style.margin = "auto"
    div_container.appendChild(gene_selection_dropdown)

    var info_icon = createInfoTooltip("fa fa-question-circle")
    info_text = "Clones affected by <b>selected target gene</b> are highlighted with colors in the mutation tree below: " +
      "<b><font color=tomato>red</font></b>  for CN amplification, <b><font color=lightsteelblue>blue</font></b> for CN deletion " +
      "and <b><font color=#b4a7d6>violet</font></b> for any other mutation event. An event in one clone of the tree affects all the subsequent nodes in the child subtree." 
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
      drug_selection_dropdown.options[drug_selection_dropdown.options.length] = new Option(drug_short_name, idx);
    }
    div_container = document.createElement('div')
    div_container.style.width = "fit-content"
    div_container.style.margin = "auto"
    div_container.appendChild(drug_selection_dropdown)

    var info_icon = createInfoTooltip("fa fa-question-circle")
    info_text = "Clones with affected genes that have a theoretical " +
      "interaction with the <b> selected target drug</b> (according to DGIdb) are indicated with green squares in the mutation tree below. " + 
       "Drugs are listed in descending order by the number of cells they affect. Only drug interactions with at least 3 citations are considered."
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
  tree_box_div.style.display = 'block'
  tree_info_div.appendChild(tree_box_div)
  appendLineBreak(tree_info_div)
  tree_div_height = Math.min(400, 2*screen.height/3)
  displayTree(tree_box_id, "", tree_object, target_gene, target_drug, drug_gene_map, height=tree_div_height)

  // Displayed metadata.
  if (metadata && mapSize(metadata)) {
    var header_metadata = createInfoHeader("<b>Clinical data</b>")
    tree_info_div.appendChild(header_metadata)

    var metadata_box_id = "metadata"
    var metadata_box_div = createDivContainer(metadata_box_id)
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
    tree_info_div.appendChild(header_top_drugs)

    var top_drugs_box_id = "top_drugs"
    var top_drugs_box_div = createDivContainer(top_drugs_box_id)
    top_drugs_box_div.style.display = 'none'
    top_drugs_box_div.innerHTML = populateGeneDrugInfoHTML(target_gene, gene_drug_map)
    tree_info_div.appendChild(top_drugs_box_div)
    header_top_drugs.appendChild(createExpandBox(top_drugs_box_id))
    appendLineBreak(tree_info_div)
  }

  // kNN matching trees.
  if (knn) { 
    var header_knn = createInfoHeader("<b>K-nearest tree neighbors</b>")
    tree_info_div.appendChild(header_knn)

    var knn_box_id = "knn"
    var knn_box_div = createDivContainer(knn_box_id)
    knn_box_div.style.display = 'none'
    tree_info_div.appendChild(knn_box_div)
    for (const [id, data] of Object.entries(knn)) {
      if (id.startsWith(sample_name)) {
        async_display_tree_matching(knn_box_id, data)
      }
    }
    header_knn.appendChild(createExpandBox(knn_box_id))
    appendLineBreak(tree_info_div)
  }

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

function getGeneCategoriesInNode(node) {
  var gene_categories = {}
  var gene_categories_with_details = {}
  if (node.data.gene_events){
    gene_map = objectToMap(node.data.gene_events)
    gene_map.forEach((events, gene) => {
      events = objectToMap(events)
      events.forEach((value, event) => {
        label = event
        if (event == "CNA") {
          if (parseInt(value) > 0 || value == "+") {
            label = "amplified"
          }
          else if (parseInt(value) < 0 || value == "-") {
            label = "deleted"
          }
          else {
            label = value
          }
        }
        if (!(label in gene_categories)) {
          gene_categories[label] = new Set()
          gene_categories_with_details[label] = new Set()
        }
        gene_categories[label].add(gene)
        if (event == "CNA") {
          gene_categories_with_details[label].add(gene)
        } else {
          gene_string = gene
          if (value) {
            gene_string += "_" + value
          }
          gene_categories_with_details[label].add(gene_string)
        }
      })
    })
  }
  return [gene_categories, gene_categories_with_details]
}

function getGeneCategoriesInTree(tree_object) {
  var node_list = getTreeNodes(tree_object)
  gene_categories = {}
  gene_categories_with_details = {}
  for (node of node_list) {
    const [node_gene_events, node_gene_events_with_details] = getGeneCategoriesInNode(node)
    for (label in node_gene_events) {
      if (!(label in gene_categories)) {
        gene_categories[label] = node_gene_events[label]
        gene_categories_with_details[label] = node_gene_events_with_details[label]
      } else {
        gene_categories[label] = gene_categories[label].union(node_gene_events[label])
        gene_categories_with_details[label] = gene_categories[label].union(node_gene_events[label])
      }
    }
  }
  return [gene_categories, gene_categories_with_details]
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
    option = document.createElement('option');
    option.textContent = gene;
    optgroup.appendChild(option);
  }
  html_element.appendChild(optgroup);
}

function getCellCountForGene(tree_object, target_gene) {
  var node_list = getTreeNodes(tree_object)
  var cell_percentage = 0
  for (var node of node_list) {
    if (node.data.gene_states) {
      if (node.data.gene_states.has(target_gene)) {
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
  if(!(target_drug in drug_gene_map)) {
    return 0
  }
  var node_list = getTreeNodes(tree_object)
  var gene_list = objectToMap(drug_gene_map).get(target_drug)
  var cell_percentage = 0
  for (node of node_list) {
    if (node.data.gene_events){
      node_genes = Array.from(node.data.gene_states.keys())
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
    return "No target gene selected.<br/>"
  }
  var drugs = gene_drug_map[target_gene]
  if (!(target_gene in gene_drug_map)) {
    gene_without_delimiter = target_gene.split(/_|-/)
    if (gene_without_delimiter.length > 1 && gene_without_delimiter[0] in gene_drug_map) {
      drugs = gene_drug_map[gene_without_delimiter[0]]
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
    if (drug.drug_score < 3) {
      continue
    }
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
    return "No drugs associated with target gene in at least 3 citations.<br/>"
  }

  html_string = ""
  if(activators.length) {
    html_string += "<b><u>Activators:</u></b><br/>"
    for (var drug of activators) {
      html_string += "<a href='https://www.dgidb.org/results?searchType=drug&searchTerms=" + drug.drug_name + 
        "' target='dgidb'>" + drug.drug_name + "</a> (" + drug.drug_score + " citations)<br/>"
    }
  }

  if(inhibitors.length) {
    html_string += "<b><u>Inhibitors:</u></b><br/>"
    for (var drug of inhibitors) {
      html_string += "<a href='https://www.dgidb.org/results?searchType=drug&searchTerms=" + drug.drug_name +
        "' target='dgidb'>" + drug.drug_name + "</a> (" + drug.drug_score + " citations)<br/>"
    }
  }

  if (others.length  && activators.length + inhibitors.length !=0) {
    html_string += "<b><u>Others:</u></b><br/>"
  }
  if(others.length) {
    for (var drug of others) {
      html_string += "<a href='https://www.dgidb.org/results?searchType=drug&searchTerms=" + drug.drug_name +
        "' target='dgidb'>" + drug.drug_name + "</a> (" + drug.drug_score + " citations)<br/>"
    }
  }
  return html_string
}

/////////////////////
//// DGIdb utils ////
/////////////////////
function parseGeneDrugInteractions(gene_drug_interaction){
  // Returns a map with gene key and drug info value.
  gene_drug_map = {}
  for (var item of gene_drug_interaction["matchedTerms"]) {
    gene = item["searchTerm"]
    gene_drug_map[gene] = []
    for (var drug of item["interactions"])  {
      drug_info = {}
      drug_info["drug_name"] = drug.drugName
      drug_info["interaction_types"] = drug.interactionTypes
      drug_info["drug_sources"] = drug.sources
      drug_info["drug_score"] = drug.score
      gene_drug_map[gene].push(drug_info)
    }
  }
  return gene_drug_map
}

function getDrugToGeneMap(gene_list, gene_drug_map) {
  var drug_gene_map = {}
  for (var gene of gene_list){
    if(gene in gene_drug_map) {
      for (var drug of gene_drug_map[gene]){
        if(drug["drug_score"] < 3) {
          continue
        }
        drug_name = drug.drug_name.substring(0,30)
        if (drug_name in drug_gene_map) {
          drug_gene_map[drug_name].push(gene)
        }
        else {
          drug_gene_map[drug_name] = [gene]
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
  for (drug in drug_gene_map) {
    drug_gene_list = drug_gene_map[drug]
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
    string = interaction.toLowerCase()
    if(string.includes("inhibitor") || string.includes("modulator") || string.includes("antagonist")) {
      return true
    }
  }
  return false
}

function isActivator(interaction_types) {
  for (var interaction of interaction_types){
    string = interaction.toLowerCase()
    if(string.includes("activator") || string.includes("inducer") || string.includes("agonist")) {
      return true
    }
  }
  return false
}

function isAntibody(interaction_types) {
  for (var interaction of interaction_types){
    string = interaction.toLowerCase()
    if(string.includes("antibody") || string.includes("binder")) {
      return true
    }
  }
  return false
}

///////////////////////////////
//// Populate cluster info ////
///////////////////////////////
function showClusterInfo(args) {
  var cluster_bg_color = args.cluster_bg_color
  var cluster_metadata = args.cluster_metadata
  var table_color_codes = args.table_color_codes
  var node_intersections = args.node_intersections 
  var matching_nodes_details = args.matching_nodes_details

  var tree_info_div_id = args.tree_info_div_id
  var tree_info_div = document.getElementById(tree_info_div_id)
  tree_info_div.innerHTML = ""
  tree_info_div.style.backgroundColor = cluster_bg_color

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

  /*tree_info_div.innerHTML = "<strong style='background-color:#d2cae6; display:block;'>&nbsp;Genes in matching nodes<br/></strong>"
  appendLineBreak(tree_info_div)
  for (match_color in matching_nodes_details) {
    tree_info_div.innerHTML += '<i class="fa fa-circle" style="font-size:18px;color:' + match_color + '"></i> &nbsp;'
    for (event in matching_nodes_details[match_color]) {
      genes = Array.from(matching_nodes_details[match_color][event]).sort()
      tree_info_div.innerHTML += event + ": " + genes.join(', ') + "; &nbsp;"
    }
    appendLineBreak(tree_info_div)
  }
  appendLineBreak(tree_info_div)*/

  bg_color = tinycolor(cluster_bg_color).darken(30).desaturate(40).toHexString()
  tree_info_div.innerHTML += "<strong style='background-color:" + bg_color + "; display:block;'>&nbsp;Cluster metadata<br/></strong>"
  tree_info_div.appendChild(table_color_codes)
  appendLineBreak(tree_info_div)
  var metadata_table = getMetadataTable(cluster_metadata)
  tree_info_div.appendChild(metadata_table)
  appendLineBreak(tree_info_div)
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

function getMetadataColorMap(sample_metadata_map) {
  metadata_color_map = {}
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
      " on the colored frames<br/>(corresponding to the given tree clusters)<br/>to visualize additional information.</i></div>"

  // Data.
  data = args.data
  distances = data["pairwise_distances"]
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
  for (var cluster of clusters) {
    cluster_starts[cluster[0]] = {
        "size": cluster.length, 
        "color": data["trees"][cluster[0]]["2d_color"],
        "cluster_bg_color": data["trees"][cluster[0]]["cluster_color"],
        "cluster_metadata": trees[cluster[0]]["sample_metadata_colors"],
        "table_color_codes": trees[cluster[0]]["table_color_codes"]
    }
  }

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
          args_cluster["matching_nodes_details"] = {} // TODO matching_clones_color_map
          args_cluster["tree_info_div_id"] = tree_info_div_id
          args_cluster["cluster_bg_color"] = cluster_starts[d.sample_1]["cluster_bg_color"]
          args_cluster["cluster_metadata"] = cluster_starts[d.sample_1]["cluster_metadata"]
          args_cluster["table_color_codes"] = cluster_starts[d.sample_1]["table_color_codes"];
          showClusterInfo(args_cluster)
      })

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
  //tree_div = createDivContainer("tree_div_2d")
  //outer_div.appendChild(tree_div)
  div_tree_cohort.appendChild(outer_div)

  var tree_info_div_id = args.tree_info_div_id
  div_tree_info = document.getElementById(tree_info_div_id)
  div_tree_info.style.border = "0px"
  div_tree_info.style.backgroundColor = "transparent"
  div_tree_info.innerHTML = "<div style='text-align:center;height:100%;display:flex;flex-direction:row;" +
      "align-items:center;justify-content:center;'><i>Click on the 2D points" +
      "<br/>to visualize additional information<br/>about the corresponding mutation tree.</i></div>"

  // Data.
  distances = args.data["pairwise_distances"]
  var opt = {epsilon: 10}; // epsilon is learning rate (10 = default)
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

