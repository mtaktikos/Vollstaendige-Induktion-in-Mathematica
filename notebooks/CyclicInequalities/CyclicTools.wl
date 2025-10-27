(* CyclicTools.wl - Utilities for cyclic sums and an automatic proof of a 3-point cyclic inequality *)

(* Copyright (c) 2025 *)

(* Core utilities for cyclic sums and simplification *)
ClearAll[VectorCyclicSum, CyclicSum, CyclicAverage, CyclicCanonicalize, FindEdgeSOSAnsatz, ProveExampleInequality];

(* Vectorized builder for local terms of the form f(u_i, u_{i-1}, u_{i+1}) *)
VectorCyclicSum[f_, u_List] := Module[
  {left = RotateLeft[u], right = RotateRight[u]},
  Total@MapThread[f, {u, right, left}]
];

(* Symbolic cyclic sum over indices with u[0]=u[n], u[n+1]=u[1] *)
CyclicSum[term_, n_Integer?Positive] := Module[{ii},
  Total@Table[
    term /. {u[0] -> u[n], u[n + 1] -> u[1], i -> ii}
  , {ii, 1, n}]
];

(* Cyclic average and canonicalization under rotations *)
CyclicAverage[expr_, vars_List] := Mean@Table[
  expr /. Thread[vars -> RotateLeft[vars, k]],
  {k, 0, Length[vars] - 1}
];

CyclicCanonicalize[expr_, vars_List] := Module[
  {rots = Table[expr /. Thread[vars -> RotateLeft[vars, k]], {k, 0, Length[vars] - 1}]},
  First@MinimalBy[rots, ToString]
];

(* Discover edge-wise SOS representation for a two-variable cubic template
   a x^3 + b y^3 + p x^2 y + q x y^2 == (x - y)^2 (c x + d y) *)
FindEdgeSOSAnsatz[p_, q_] := Module[{x, y, a, b, c, d, poly, eqs},
  poly = a x^3 + b y^3 + p x^2 y + q x y^2 - (x - y)^2 (c x + d y);
  eqs = Thread[Flatten@CoefficientList[Expand@poly, {x, y}] == 0];
  Solve[eqs, {a, b, c, d}, Reals]
];

(* Prove: Sum[u_i (u_{i-1}^2 - 4 u_i^2 + 3 u_{i+1}^2)] <= 0 for u_i >= 0 *)
ProveExampleInequality[N_Integer?Positive] := Module[
  {u, left, right, S, sol, a, b, c, d, P, idOK, nonnegOK, concl},
  u = Array[u, N];
  left = RotateLeft[u];
  right = RotateRight[u];
  S = Total[u*(right^2 - 4 u^2 + 3 left^2)] // Expand;
  (* Find coefficients for (u_i - u_{i+1})^2 (c u_i + d u_{i+1}) >= 0 *)
  sol = First@FindEdgeSOSAnsatz[-1, -3]; (* matches -u_i^2 u_{i+1} - 3 u_i u_{i+1}^2 *)
  {a, b, c, d} = {a, b, c, d} /. sol;  (* Expect 5/3, 7/3, 5/3, 7/3 *)
  
  (* Edge-sum nonnegative polynomial *)
  P = Total[(u - left)^2*(c u + d left)] // Expand;
  
  (* Identity: P == (a + b) Sum u_i^3 - Sum u_i u_{i-1}^2 - 3 Sum u_i u_{i+1}^2 *)
  idOK = Simplify[
    Expand[P - ((a + b) Total[u^3] - Total[u*right^2] - 3 Total[u*left^2])]
  ] === 0;
  
  (* Nonnegativity for u_i >= 0 *)
  nonnegOK = Simplify[P >= 0, Assumptions -> And @@ Thread[u >= 0]];
  
  (* Desired conclusion *)
  concl = Simplify[S <= 0, Assumptions -> And @@ Thread[u >= 0]];
  
  <|
    "IdentityHolds" -> idOK,
    "EdgeWeights" -> <|"a" -> a, "b" -> b, "c" -> c, "d" -> d|>,
    "NonnegativityCertified" -> nonnegOK,
    "Conclusion" -> concl
  |>
];

(* End of CyclicTools.wl *)