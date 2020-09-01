var outerRadius = 1000 / 2,
    innerRadius = outerRadius - 240;
var formatNum = d3.format(".3");
var arc = d3.arc()
    .innerRadius(innerRadius)
    .outerRadius(innerRadius+24);

var arc2 = d3.arc()
    .innerRadius(innerRadius+24)
    .outerRadius(innerRadius+50);


function tooltip(g, dom, element, content) {
  g.selectAll(element)
    .on("mouseover", function(d) {
      d3.select("text#"+"dc"+d.data.sorter).remove(); // boundary case
      d3.select(dom)
        .append("text")
        .style("background", "black")
        .style("color", "white")
        .attr("dy", "0.31em")
        .attr("id", "dc"+d.data.sorter)
        .attr("transform", function() { return "translate(" + (innerRadius+180) + "," + (innerRadius + 240) + ")"+"translate("+arc.centroid(d)+")"; })
        .text(function() {
          if (content == "rp")
            return d.data.name + " " + formatNum(d.data["rp"]);
          else if (content == "p") {
            return d.data.name + " " + formatNum(d.data["p"]);
          } else {
            return d.data.name;
          }});
    });
  g.selectAll(element)
    .on("mouseout", function(d) {
      d3.select("text#"+"dc"+d.data.sorter)
        .transition()
        .duration(600)
        .style("opacity", 0)
        .remove();
        //.attr("transform", "translate(" + (innerRadius+120) + "," + (innerRadius+120) + ")"+"translate(" + arc.centroid(d) +")")
    });
}

function draw(error, life, mytype, dom) {

  if (error) throw error;

  var cluster = d3.cluster()
      .size([360, innerRadius])
      .separation(function(a, b) { return 1.8; });

  var svg = d3.select(dom)
      .attr("width", outerRadius * 2)
      .attr("height", outerRadius * 2);

  var chart = svg.append("g")
      .attr("transform", "translate(" + (outerRadius) + "," + (outerRadius) + ")");

  var root = d3.hierarchy(life)
      .sum(function(d) { return d.children ? 0 : 1; })
      .sort(function(a, b) { return (a.value - b.value) || d3.ascending(a.data.length, b.data.length); });

  if (mytype=="margerp") {
  function getFamily(d) {
    if (d.children) {
      family.push(d.data.name);
      d.children.forEach(getFamily);
    }
  }

  //var family = [];
  //getFamily(root);
  //family = new Set(family.slice(1));
  //family = [...family];

  //var color = d3.scaleOrdinal(d3.schemeCategory10)
  //    .domain(family);

  //var legend = svg.append("g")
  //    .attr("class", "legend")
  //    .selectAll("g")
  //    .data(color.domain())
  //    .enter().append("g")
  //    .attr("transform", function(d, i) { return "translate(" + (outerRadius * 2 - 240) + "," + (i * 20 + 10) + ")"; });

  //legend.append("rect")
  //  .attr("x", -18)
  //  .attr("width", 18)
  //  .attr("height", 18)
  //  .attr("fill", color);

  //legend.append("text")
  //  .attr("x", 1)
  //  .attr("y", 9)
  //  .attr("dy", ".35em")
  //  .attr("text-anchor", "start")
  //    .text(function(d) { return d; });
  }

  cluster(root);

  // var input = d3.select("#show-length input").on("change", changed),
  // timeout = setTimeout(function() { input.property("checked", true).each(changed); }, 1000);

  setRadius(root, root.data.length = 0, innerRadius / maxLength(root));

  // // Set the color of each node by recursively inheriting.
  // function setColor(d) {
  //   var name = d.data.name;
  //   d.color = color.domain().indexOf(name) >= 0 ? color(name) : d.parent ? d.parent.color : null;
  //   if (d.children) d.children.forEach(setColor);
  // }

  // if (mytype=='margerp')
  //   setColor(root);

  var linkExtension = chart.append("g")
      .attr("class", "link-extensions")
      .selectAll("path")
      .data(root.links().filter(function(d) { return !d.target.children; }))
    .enter().append("path")
    .each(function(d) { d.target.linkExtensionNode = this; })
    .attr("d", linkExtensionConstant);

  var link = chart.append("g")
      .attr("class", "links")
      .selectAll("path")
      .data(root.links())
    .enter().append("path")
    .each(function(d) { d.target.linkNode = this; })
    .attr("d", linkConstant)
      .attr("stroke", function(d) { return d.target.color; });

  var g = chart.append("g")
      .attr("class", "labels");
  var leaves = root.leaves();
  g.selectAll("text")
    .data(root.leaves())
    .enter().append("text")
    .attr("dy", "0.31em")
    .attr("transform", function(d) { return "rotate(" + (d.x - 90) + ")translate(" + (innerRadius + 60) + ",0)" + (d.x < 180 ? "" : "rotate(180)"); })
    .attr("text-anchor", function(d) { return d.x < 180 ? "start" : "end"; })
    .text(function(d) { return d.data.name.replace(/::/g, " "); })
    .on("mouseover", mouseovered(true))
    .on("mouseout", mouseovered(false))
    .style("opacity", 0)
    .transition()
    .duration(2000)
    .style("opacity", 1);

  var x = d3.scaleLinear()
      .range([0, 2 * Math.PI]);

  var y = d3.scaleSqrt()
      .range([0, outerRadius]);

  var pd = new Array();


  for (var i = 0; i < leaves.length; i++) {
    if (mytype == 'margerp')
      v = leaves[i].data[mytype];
    else
      v = -Math.log10(leaves[i].data[mytype]);
    color_values = {"sorter": i, "pie": 1, "rp": v, "p": -Math.log10(leaves[i].data.p), "name": leaves[i].data.name};
    pd.push(color_values);
  }

  var pext = d3.extent(pd, function(d) {return d.p;});

  var rext = d3.extent(pd, function(d) {return d.rp;});

  var pie = d3.pie().value(function(d) { return d.pie; })
      .sort(function(a,b) { return d3.ascending(a.sorter, b.sorter);});

  var rainbow = d3.scaleLinear()
    .domain(pext)
    .range(["white", "blue"])
    .interpolate(d3.interpolateRgb);

  var rainbow2 = d3.scaleLinear()
    .domain(rext)
    .range(["white", "red"])
    .interpolate(d3.interpolateRgb);

  g.selectAll("path.arc")
    .data(pie(pd)).enter()
    .append("path")
    .attr("class", "arc")
    .attr("fill", function(d) {
      if (d.data.p == "nan")
        return "white";
      else
        return rainbow(d.data.p);
    }).attr("d", arc)
    .style("opacity", 0)
    .transition()
    .duration(900)
    .style("opacity", 1);

  tooltip(g, dom, "path.arc", "p");

  g.selectAll("path.arc2")
    .data(pie(pd)).enter()
    .append("path")
    .attr("class", "arc2")
    .attr("fill", function(d) {
      if (d.data.rp == "nan")
        return "white";
      else
        return rainbow2(d.data.rp);
    }).attr("d", arc2)
    .style("opacity", 0)
    .transition()
    .duration(1500)
    .style("opacity", 1);

  tooltip(g, dom, "path.arc2", "rp");

  function changed() {
    clearTimeout(timeout);
    var t = d3.transition().duration(750);
    linkExtension.transition(t).attr("d", this.checked ? linkExtensionVariable : linkExtensionConstant);
    link.transition(t).attr("d", this.checked ? linkVariable : linkConstant);
  }

  function mouseovered(active) {
    return function(d) {
      d3.select(this).classed("label--active", active);
      d3.select(d.linkExtensionNode).classed("link-extension--active", active).each(moveToFront);
      do d3.select(d.linkNode).classed("link--active", active).each(moveToFront); while (d = d.parent);
    };
  }

  function moveToFront() {
    this.parentNode.appendChild(this);
  }
};

// Compute the maximum cumulative length of any node in the tree.
function maxLength(d) {
  return d.data.length + (d.children ? d3.max(d.children, maxLength) : 0);
}

// Set the radius of each node by recursively summing and scaling the distance from the root.
function setRadius(d, y0, k) {
  d.radius = (y0 += d.data.length) * k;
  if (d.children) d.children.forEach(function(d) { setRadius(d, y0, k); });
}


function linkVariable(d) {
  return linkStep(d.source.x, d.source.radius, d.target.x, d.target.radius);
}

function linkConstant(d) {
  return linkStep(d.source.x, d.source.y, d.target.x, d.target.y);
}

function linkExtensionVariable(d) {
  return linkStep(d.target.x, d.target.radius, d.target.x, innerRadius);
}

function linkExtensionConstant(d) {
  return linkStep(d.target.x, d.target.y, d.target.x, innerRadius);
}

// Like d3.svg.diagonal.radial, but with square corners.
function linkStep(startAngle, startRadius, endAngle, endRadius) {
  var c0 = Math.cos(startAngle = (startAngle - 90) / 180 * Math.PI),
  s0 = Math.sin(startAngle),
  c1 = Math.cos(endAngle = (endAngle - 90) / 180 * Math.PI),
  s1 = Math.sin(endAngle);
  return "M" + startRadius * c0 + "," + startRadius * s0
    + (endAngle === startAngle ? "" : "A" + startRadius + "," + startRadius + " 0 0 " + (endAngle > startAngle ? 1 : 0) + " " + startRadius * c1 + "," + startRadius * s1)
    + "L" + endRadius * c1 + "," + endRadius * s1;
}

var margin = {top: 22, right: 30, bottom: 22, left: 30},
    width = 450 - margin.left - margin.right,
    height = 420 - margin.top - margin.bottom;

var x = d3.scaleLinear().range([0, width]),
    y = d3.scaleBand().rangeRound([height, 0]),
    y2 = d3.scaleLinear().range([height, 0]);

function bar(data, dom) {
  var xAxis = d3.axisBottom(x).ticks(8);
  var yAxis = d3.axisLeft(y);

  var svg = d3.select(dom)
      .attr("width", width + margin.left + margin.right)
      .attr("height", height + margin.top + margin.bottom)
      .append("g")
      .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

  data.sort(function(a,b) { return d3.ascending(parseFloat(a.value),parseFloat(b.value));});
  // var extent = d3.extent(data, function(d) { return d.value; });
  // extent = [parseFloat(extent[0]) - 0.5, parseFloat(extent[1]) + 0.2];
  var ye = data.map(function(d) { return d.name; });
  x.domain([-1, 1]);
  y.domain(ye);

  svg.selectAll(".bar")
    .data(data)
    .enter()
    .append("rect")
    .attr("class", function(d) { return "bar bar--" + (d.value<0 ? "negative": "positive");})
    .attr("x", function(d) {return x(Math.min(0, d.value));})
    .attr("y", function(d) {return y(d.name);})
    .attr("width", function(d) { return Math.abs(x(d.value)-x(0));})
    .attr("height", function(d) {return y.bandwidth(); });

  svg.append("g")
    .attr("class", "axis")
    .attr("transform", "translate(0," + height +")")
    .call(xAxis).selectAll("text")
    .style("text-anchor", "end")
    .attr("dx", "-.8em")
    .attr("dy", ".15em")
    .attr("transform", function(d) {
      return "rotate(-65)";
    });

  var ya = svg.append("g")
    .attr("class", "yaxis")
    .attr("transform", "translate(" + x(0) +",0)")
    .call(yAxis).selectAll("text")
    .attr("x", function(d, i) {
      if (data[i].value < 0)
        return y(0);
      else
        return y(0);
    })
    .attr("text-anchor", function(d, i){
      if (data[i].value <= 0)
        return "start";
      else
        return "end";
    });
}

function auc_curve(root, dom) {
  var x = d3.scaleLinear().range([0, width]),
  y2 = d3.scaleLinear().range([height, 0]);
  x.domain([0,1]);
  y2.domain([0,1]);

  var data = new Array();
  for (var i = 0; i < root.fpr.length; i++) {
    var obj = {"fpr": root.fpr[i], "tpr": root.tpr[i]};
    data.push(obj);
  }

  var svg = d3.select(dom)
    .attr("width", width + margin.left + margin.right)
    .attr("height", height + margin.top + margin.bottom)
    .append("g")
    .attr("transform", "translate(" + margin.left + "," + margin.top + ")");
  var valueline = d3.line()
    .curve(d3.curveBasis)
    .x(function(d, i) { return x(d.fpr); })
    .y(function(d, i) { return y2(d.tpr); });


  var xAxis = d3.axisBottom(x).ticks(5);
  var yAxis = d3.axisLeft(y2).ticks(5);

  svg.append("path")
    .style("stroke", "steelblue")
    .style("stroke-width", "2.5px")
    .style("fill", "none")
    .attr("d", valueline(data));

  svg.append("g")
    .attr("class", "x axis")
    .attr("transform", "translate(0," + height + ")")
    .call(xAxis);
  svg.append("g")
    .attr("class", "x axis")
    .call(yAxis);

  var legend = svg.append("g")
    .attr("transform", "translate(" + width/10 + "," + height/10 + ")");
  var formatNum = d3.format(".3");

  legend.append("rect")
    .attr("x", x(0)-margin.left+5)
    .attr("y", y2(1)-margin.top+5)
    .attr("fill", "steelblue")
    .attr("width", width/25)
    .attr("height", height/25);

  legend.append("text")
    .attr("x", x(0)- margin.left+25)
    .attr("y", y2(1)-margin.top+23)
    .attr("dy", "-0.51em")
    .attr("dx", "0.58em")
    .text("AUC:" + formatNum(root.performance[0]));
}


function run_all(gl) {
    d3.tsv("H3K27ac_" + gl + ".coef", function(error, data) {
    bar(data, "svg.bar");
    });
    d3.json("H3K27ac_" + gl + "_performance.json", function(error, data) {
    auc_curve(data, "svg.auc");
    });
    var mytype = 'margerp';
    var dom = 'svg.tree';
    d3.json(gl + "_delta_batch_H3K27ac_withpromoter_dnase_99.json",
    function(error, data) {
    draw(error, data, mytype, dom);
    });

    d3.json(gl + "_delta_batch_H3K27ac_withpromoter_dnase_cistrome_dc.json", function(error, data) {
    draw(error, data, "beta", "svg.dc_tree");
    });

    $('.nav-tabs a[href="#DNase"]').click(function(){
    $(this).tab('show');

    d3.json("DNase_" + gl + "_performance.json", function(error, data) {
    auc_curve(data, "svg.auc2");
    });

    d3.tsv("DNase_" + gl + ".coef", function(error, data) {
    bar(data, "svg.bar2");
    });

    var mytype = 'margerp';
    var dom = 'svg.tree2';
    d3.json(gl + "_delta_batch_DNase_withpromoter_dnase_99.json",
    function(error, data) {
    draw(error, data, mytype, dom);
    });

    d3.json(gl + "_delta_batch_DNase_withpromoter_dnase_cistrome_dc.json", function(error, data) {
    draw(error, data, "beta", "svg.dc_tree2");
    });
    })

    $('.nav-tabs a[href="#H3K4me1"]').click(function(){
    d3.json("H3K4me1_" + gl + "_performance.json", function(error, data) {
    auc_curve(data, "svg.auc3");
    });
    d3.tsv("H3K4me1_" + gl +".coef", function(error, data) {
    bar(data, "svg.bar3");
    })
    var mytype = 'margerp';
    var dom = 'svg.tree3';
    d3.json(gl + "_delta_batch_H3K4me1_withpromoter_dnase_99.json",
    function(error, data) {
    draw(error, data, mytype, dom);
    });

    d3.json(gl + "_delta_batch_H3K4me1_withpromoter_dnase_cistrome_dc.json", function(error, data) {
    draw(error, data, "beta", "svg.dc_tree3");
    });
    })

    $('.nav-tabs a[href="#H3K4me3"]').click(function(){
    d3.json("H3K4me3_" + gl + "_performance.json", function(error, data) {
    auc_curve(data, "svg.auc4");
    });
    d3.tsv("H3K4me3_" + gl + ".coef", function(error, data) {
    bar(data, "svg.bar4");
    })
    var mytype = 'margerp';
    var dom = 'svg.tree4';
    d3.json(gl + "_delta_batch_H3K4me3_withpromoter_dnase_99.json",
    function(error, data) {
    draw(error, data, mytype, dom);
    });

    d3.json(gl + "_delta_batch_H3K4me3_withpromoter_dnase_cistrome_dc.json", function(error, data) {
    draw(error, data, "beta", "svg.dc_tree4");
    });
    })

    $('.nav-tabs a[href="#H3K4me2"]').click(function(){
    d3.json("H3K4me2_" + gl + "_performance.json", function(error, data) {
    auc_curve(data, "svg.auc5");
    });
    d3.tsv("H3K4me2_" + gl +".coef", function(error, data) {
    bar(data, "svg.bar5");
    })
    var mytype = 'margerp';
    var dom = 'svg.tree5';
    d3.json(gl + "_delta_batch_H3K4me2_withpromoter_dnase_99.json",
    function(error, data) {
    draw(error, data, mytype, dom);
    });
    d3.json(gl + "_delta_batch_H3K4me2_withpromoter_dnase_cistrome_dc.json", function(error, data) {
    draw(error, data, "beta", "svg.dc_tree5");
    });
    });

    $('.nav-tabs a[href="#ATAC-seq"]').click(function(){
    d3.json("ATAC-seq_" + gl + "_performance.json", function(error, data) {
    auc_curve(data, "svg.auc6");
    });
    d3.tsv("ATAC-seq_" + gl + ".coef", function(error, data) {
    bar(data, "svg.bar6");
    })
    var mytype = 'margerp';
    var dom = 'svg.tree6';
    d3.json(gl + "_delta_batch_ATAC-seq_withpromoter_dnase_99.json",
    function(error, data) {
    draw(error, data, mytype, dom);
    });
    d3.json(gl + "_delta_batch_ATAC-seq_withpromoter_dnase_cistrome_dc.json", function(error, data) {
    draw(error, data, "beta", "svg.dc_tree6");
    });
    })

    $('.nav-tabs a[href="#H3K9ac"]').click(function(){
    d3.json("H3K9ac_" + gl + "_performance.json", function(error, data) {
    auc_curve(data, "svg.auc7");
    });
    d3.tsv("H3K9ac_" + gl + ".coef", function(error, data) {
    bar(data, "svg.bar7");
    })
    var mytype = 'margerp';
    var dom = 'svg.tree7';
    d3.json(gl + "_delta_batch_H3K9ac_withpromoter_dnase_99.json",
    function(error, data) {
    draw(error, data, mytype, dom);
    });
    d3.json(gl + "_delta_batch_H3K9ac_withpromoter_dnase_cistrome_dc.json", function(error, data) {
    draw(error, data, "beta", "svg.dc_tree7");
    });
    });


    $('.nav-tabs a[href="#H3K27me3"]').click(function(){
    d3.json("H3K27me3_" + gl + "_performance.json", function(error, data) {
    auc_curve(data, "svg.auc8");
    });
    d3.tsv("H3K27me3_" + gl + ".coef", function(error, data) {
    bar(data, "svg.bar8");
    })
    var mytype = 'margerp';
    var dom = 'svg.tree8';
    d3.json(gl + "_delta_batch_H3K27me3_withpromoter_dnase_99.json",
    function(error, data) {
    draw(error, data, mytype, dom);
    });
    d3.json(gl + "_delta_batch_H3K27me3_withpromoter_dnase_cistrome_dc.json", function(error, data) {
    draw(error, data, "beta", "svg.dc_tree8");
    });
    });

    $('.nav-tabs a[href="#H3K36me3"]').click(function(){
    d3.json("H3K36me3_" + gl + "_performance.json", function(error, data) {
    auc_curve(data, "svg.auc9");
    });
    d3.tsv("H3K36me3_" + gl + ".coef", function(error, data) {
    bar(data, "svg.bar9");
    })
    var mytype = 'margerp';
    var dom = 'svg.tree9';
    d3.json(gl + "_delta_batch_H3K36me3_withpromoter_dnase_99.json",
    function(error, data) {
    draw(error, data, mytype, dom);
    });
    d3.json(gl + "_delta_batch_H3K36me3_withpromoter_dnase_cistrome_dc.json", function(error, data) {
    draw(error, data, "beta", "svg.dc_tree9");
    });
    });


    $('.nav-tabs a[href="#H3K9me3"]').click(function(){
    d3.json("H3K9me3_" + gl + "_performance.json", function(error, data) {
    auc_curve(data, "svg.auc10");
    });
    d3.tsv("H3K9me3_" + gl + ".coef", function(error, data) {
    bar(data, "svg.bar10");
    })
    var mytype = 'margerp';
    var dom = 'svg.tree10';
    d3.json(gl + "_delta_batch_H3K9me3_withpromoter_dnase_1010.json",
    function(error, data) {
    draw(error, data, mytype, dom);
    });
    d3.json(gl + "_delta_batch_H3K9me3_withpromoter_dnase_cistrome_dc.json", function(error, data) {
    draw(error, data, "beta", "svg.dc_tree10");
    });
    });



}
