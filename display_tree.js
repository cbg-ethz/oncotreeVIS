////////// Help functions ///////////
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
function displayTree(div_id, tree_label, tree_data, target_gene, target_drug, drug_gene_map, tree_info_view=true, width=200, height=400) {

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
        return "dimgray"
      } 
      if (d.parent != null && d.data.gene_events) {
        cn_event = get_cn_event_for_gene(d, target_gene);
        return get_color_for_cn_event(cn_event);
      }
      return "darkgray"
    });

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
      if (d.parent) {
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
        if (gene_events) {
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
      } else {
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
        return label
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

