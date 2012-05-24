/*
There is a matrix = grid of temperatures. The temperatures on the edges are fixed, while the inner
vertices should react physically plausible to their surrounding vertices over time. The algorithm
can either iteratively model the heat distribution, or try to find the end result directly.
*/

var M = 20, // row count of heat grid
    N = 20, // col count of heat grid
    heat = []; // array of temperature matrices

// visualization parameters
var w = 960,    // visualization width in pixels
    h = 580,    // visualization height in pixels
    colormap = d3.interpolateHsl("blue", "red");
    
var timer_id = null
   ,step_count = 0
   ,err_margin = 1e-3
   ,methods = [new MethodIterative(), new MethodItSolver(), new MethodItSolverCG()]
   //,methods = [new MethodItSolverCG()]
   ,use_sparse = true
   ,iteration_counter = {n: 0}
   ,timestep = 100;

function MethodIterative() {
  this.name = 'Local Updates';
  this.err = Infinity;
  this.eps = 0;
  this.init = function(heat, eps) { this.A = heat; this.eps = eps; };
  this.iteration = 0;
  this.step = function() {
    var dt = 1;
    var err = 0;
    for (var i=1; i<M-1; i++) {
      for (var j=1; j<N-1; j++) {
        var delta = dt*(this.A[i-1][j] + this.A[i+1][j] + this.A[i][j-1] + this.A[i][j+1] - 4*this.A[i][j])*0.25;
        err += delta*delta;
        this.A[i][j] += delta;
      }
    }
    this.err = Math.sqrt(err);
    if (this.err < this.eps) this.finished = true;
    this.iteration++;
  };
}

function MethodItSolver() {
  this.name = 'Iterative Solver';
  this.err = Infinity;
  this.eps = 0;
  this.iteration = 0;
  this.init = function(heat, eps) {
    this.A = heat;
    this.solver = init_solver(heat, false);
    this.eps = eps;
  };
  this.step = function() {
    this.err = this.solver.step();
    setHeatFromX(this.A, this.solver.x);
    if (this.err < this.eps) this.finished = true;
    this.iteration++;
  };
}

function MethodItSolverCG() {
  this.name = 'Iterative Solver CG';
  this.err = Infinity;
  this.finished = false;
  this.iteration = 0;
  this.init = function(heat, eps) {
    this.A = heat;
    this.solver = init_solver(heat, true);  
    this.eps = eps;
  };
  this.step = function() {
     this.err = this.solver.step();
     setHeatFromX(this.A, this.solver.x);
     if (this.err < this.eps) this.finished = true;
     this.iteration++;
  };
}

function init() {
  // initialize temperature matrices
  for (var m=0; m<methods.length; m++) {
    heat[m] = Matrix.construct(M, N, 0);
    // set the right edge and half of left edge to hot value
    for (var i=0; i<M; i++) heat[m][i][0] = 1;
    for (var i=Math.round(M/2); i<M; i++) heat[m][i][N-1] = 1;
  }
  // show the heat maps
  init_vis(heat);

  // initialize methods
  for (var m=0; m<methods.length; m++) methods[m].init(heat[m], err_margin);

  // start animation
  console.log('Using', use_sparse ? 'sparse' : 'dense','matrices.');
  console.log('Variables:', M*N - 2*(M+N) + 4);
  
  var finished = 0;
  timer_id = setInterval(function() {
    iteration_counter.n++;
    for (var m=0; m<methods.length; m++) {
      var meth = methods[m];
      if (!meth.finished) {
        meth.step();
        if (meth.finished) {
          console.log('Method', meth.name, 'finished after', step_count, 'iterations.');
          finished++;
        }
      }
    }
    if (finished >= methods.length) clearInterval(timer_id);
    update_vis();
  }, timestep);  
}

function step(dt) {
  var err = 0;
  for (var i=1; i<M-1; i++) {
    for (var j=1; j<N-1; j++) {
      var delta = dt*(heat[i-1][j] + heat[i+1][j] + heat[i][j-1] + heat[i][j+1] - 4*heat[i][j])*0.25;
      err += delta*delta;
      heat[i][j] += delta;
    }
  }
  return Math.sqrt(err);
}

/// Uses the iterative solver to find the equilibriuum state of the heat grid.
function init_solver(heat, use_cg) {
  var len = M*N - 2*(M+N) + 4;
  var A = use_sparse ? new SparseMatrix(len, len) : Matrix.construct(len, len, 0);
  var b = Vector.construct(len, 0);
  var to_col = function(idx) { return 1 + Math.floor(idx/(M-2)); }
  var to_row = function(idx) { return 1 + idx % (M-2); }
  var to_idx = function(row, col) { return (row-1)+(col-1)*(M-2); }
  var neighbors = [[-1,0],[1,0],[0,-1],[0,1]];
  for (var idx=0; idx<len; idx++) {
    for (var n in neighbors) {
      var i = to_row(idx) + neighbors[n][0];
      var j = to_col(idx) + neighbors[n][1];
      if (i==0 || i==M-1 || j==0 || j==N-1) b[idx] += heat[i][j];
      else use_sparse ? A.Set(idx, to_idx(i,j), -1) : A[idx][to_idx(i,j)] = -1;
    }
    use_sparse ? A.Set(idx, idx, 4) : A[idx][idx] = 4;;
  }
  return use_cg ? new IterativeSolverCG(A,b) : new IterativeSolver(A, b);
}

function setHeatFromX(heat, x) {
  var to_col = function(idx) { return 1 + Math.floor(idx/(M-2)); }
  var to_row = function(idx) { return 1 + idx % (M-2); }
  for (var i=0; i<x.length; i++) {
    heat[to_row(i)][to_col(i)] = x[i];
  }  
}

function init_vis(heat) {
  // construct a table for each heat matrix
  for (var h=0; h<heat.length; h++) {
    // each element's data is {A: heat_matrix, i: row, j: col]
    var div = d3.select("body").append("div").classed("content", true);
    div.selectAll("p").data([methods[h]]).enter().append("p").classed("name", true);
    var tr = div.append("table").selectAll("tr")
      .data(function() { var h=[]; for (var i=0; i<M; i++) h.push(i); return h;})
      .enter().append("tr");

    var td = tr.selectAll("td")
      .data(function(d) {var x = []; for (var j=0; j<N; j++) x.push({A: heat[h], i:d, j:j}); return x;})
      .enter().append("td")
      .style('width', function(d) { return (d.j == 0 || d.j == N-1) ? '3px' : '10px'})
      .style('height', function(d) { return (d.i == 0 || d.i == M-1) ? '3px' : '10px'});
   
    div.selectAll("p.counter").data([methods[h]]).enter().append("p").classed("counter", true);
  }
  update_vis();
}

function update_vis() {
  d3.selectAll("table").selectAll("td")
    .style("background-color", function(d) { return colormap(d.A[d.i][d.j]); });
  d3.selectAll("p.name").text(function(d) { return d.name + (d.finished ? ' (converged)' : ''); });
  d3.selectAll("p.counter").text(function(d) { return 'Iteration: ' + d.iteration; });
}
