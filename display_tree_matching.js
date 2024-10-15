
function getYShift(depth) {
  return depth * 60 + 60
}

function getXShift(dx, group) {
  if (group == 2) {
    return dx * 4 + 170
  }
  return dx * 4 + 30
}

async function async_display_tree_matching(div_id, data) {
  display_tree_matching(div_id, data)
  await sleep(2000)
}

function display_tree_matching(div_id, data) {
  var colors = ['#a6cee3', '#1f78b4', '#b2df8a', '#33a02c', '#fb9a99', '#e31a1c', '#fdbf6f', '#ff7f00', '#cab2d6', '#6a3d9a', '#ffff99', '#b15928'] //d3.schemePaired
  var color_index = 0
  //var colors = d3.scaleOrdinal(d3.schemeCategory10);
  var node_colors = {}

  sample_1 = data["sample_1"]
  sample_2 = data["sample_2"]
  tree_type = data["tree_type"]

  nodes = data["nodes"]
  links = data["links"]
  score = data["score"]  
  max_score = data["max_score"]
  percentage = data["percentage"]
  num_significant_matches = data["num_significant_matches"]
  
  var margin = { top: 0, right: 0, bottom: 0, left: 0};
  var width = 500 //- margin.left - margin.right;
  var height = Math.min(400, 2*screen.height/3) - margin.top - margin.bottom;
 
  pair_id = sample_1 + "-" + sample_2 + "-" + tree_type  

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
    .text(sample_1 + " - " + sample_2);

  for (var i = 0; i < nodes.length; i++){
    node = nodes[i]
    group = 2
    if (node.id.includes(sample_1)) {
      group = 1
    }
    nodes[i].x = getXShift(node.dx, group)
    nodes[i].y = getYShift(node.depth)
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
    })
    .style("stroke-width", "4px")
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
        return "sim = " + d.similarity   
        return ((d.target).plit('_')[1])//.toString()
      }
    })

  var node = svg.selectAll(".node")
    .data(nodes)
    .enter()
    .append("g")
    .attr("class", "node " + "node_" + pair_id)

  node.append("circle")
    .attr("r", function(d) {
      if (d.matching_label == 0) {
        return 3;
      } 
      return 10
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

