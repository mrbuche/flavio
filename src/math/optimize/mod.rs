mod newton;

// Newton doesn't need the objective function,
// so whether it's actually an optimization or root finding problem,
// the i/o is the same, so it can live in optimize safely.
