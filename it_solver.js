/// Can iteratively solve the equation system Ax = b for x.
IterativeSolver = function(A, b) {
  if (A.N != b.length) throw "dimension mismatch, A must have as many rows as b";
  this.A = A;
  this.b = b;
  this.x = Vector.construct(b.length, 0); // current best solution
  this.r = b.copy(); // current residuum
  this.alpha_scaling = 1; // alpha scaling: set to value < 1 (e.g. 0.35) to avoid zig-zag course
}

/// Performs one iteration step towards the solution improving the vector x. Returns the norm of the
/// new residuum.
IterativeSolver.prototype.step = function() {
  var Ar = this.A.mul(this.r); // this takes the most time --> use sparse matrices if possible!
  var alpha = this.alpha_scaling*this.r.mul(this.r) / this.r.mul(Ar); // step width: alpha = r(t)*r(t) / r(t)*A*r(t)
  this.x.Add(this.r.scale(alpha)); // update x: x(t+1) = x(t) + alpha*r(t)
  this.r.Sub(Ar.scale(alpha)); // update residuum: r(t+1) = r(t) - alpha*A*r(t)
  var err = this.r.len()*alpha;
  return err;
}

/// Conjugate Gradients Iterative Solver.
/// Can iteratively solve the equation system Ax = b for x.
IterativeSolverCG = function(A, b) {
  if (A.N != b.length) throw "dimension mismatch, A must have as many rows as b";
  this.A = A;
  this.b = b;
  this.x = Vector.construct(b.length, 0); // current best solution
  this.r = b.copy(); // current residuum
  this.p = b.copy(); // current direction
}

/// Performs one iteration step towards the solution improving the vector x. Returns the norm of the
/// new residuum.
IterativeSolverCG.prototype.step = function() {
  var Ap = this.A.mul(this.p); // this takes the most time --> use sparse matrices if possible!
  var r_sq = this.r.mul(this.r);
  var alpha = r_sq / this.p.mul(Ap); // step width: alpha = r(t)*r(t) / p(t)*A*p(t)
  this.x.Add(this.p.scale(alpha)); // update x: x(t+1) = x(t) + alpha*p(t)
  this.r.Sub(Ap.scale(alpha)); // update residuum: r(t+1) = r(t) - alpha*A*p(t)
  var beta = this.r.mul(this.r) / r_sq; // step width: beta = r(t+1)*r(t+1) / r(t)*r(t)
  this.p.Scale(beta); this.p.Add(this.r); // update direction: p(t+1) = r(t+1) + beta*p(t)
  var err = this.p.len()*alpha;
  return err;
}
