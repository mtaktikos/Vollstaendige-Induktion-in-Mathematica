(* ::Package:: *)

(* CyclicInequalities Package *)
(* Utilities for cyclic sums and automatic proof of cyclic inequalities *)

BeginPackage["CyclicInequalities`"]

(* Usage messages *)
VectorCyclicSum::usage = "VectorCyclicSum[f, u] builds a cyclic sum of local terms of the form f(u_i, u_{i-1}, u_{i+1}) where u is a list.";

CyclicSum::usage = "CyclicSum[term, n] computes a symbolic cyclic sum over indices with cyclic boundary conditions u[0]=u[n], u[n+1]=u[1].";

CyclicAverage::usage = "CyclicAverage[expr, vars] computes the average of expr under all cyclic rotations of vars.";

CyclicCanonicalize::usage = "CyclicCanonicalize[expr, vars] returns the lexicographically minimal rotation of expr with respect to vars.";

FindEdgeSOSAnsatz::usage = "FindEdgeSOSAnsatz[p, q] finds the edge-wise SOS representation for a x^3 + b y^3 + p x^2 y + q x y^2 == (x - y)^2 (c x + d y).";

ProveExampleInequality::usage = "ProveExampleInequality[N] proves the inequality Sum[u_i (u_{i-1}^2 - 4 u_i^2 + 3 u_{i+1}^2)] <= 0 for u_i >= 0 with N variables.";

Begin["`Private`"]

(* Vectorized builder for local terms of the form f(u_i, u_{i-1}, u_{i+1}) *)
VectorCyclicSum[f_, u_List] := Module[{left = RotateLeft[u], right = RotateRight[u]},
  Total @ MapThread[f, {u, right, left}]
];

(* Symbolic cyclic sum over indices with u[0]=u[n], u[n+1]=u[1] *)
CyclicSum[term_, n_Integer?Positive] := Module[{ii},
  Total @ Table[term /. {u[0] -> u[n], u[n + 1] -> u[1], i -> ii}, {ii, 1, n}]
];

(* Cyclic average and canonicalization under rotations *)
CyclicAverage[expr_, vars_List] := Mean @ Table[expr /. Thread[vars -> RotateLeft[vars, k]], {k, 0, Length[vars] - 1}];

CyclicCanonicalize[expr_, vars_List] := Module[
  {rots = Table[expr /. Thread[vars -> RotateLeft[vars, k]], {k, 0, Length[vars] - 1}]},
  First @ MinimalBy[rots, ToString]
];

(* Discover edge-wise SOS representation for a two-variable cubic template
   a x^3 + b y^3 + p x^2 y + q x y^2 == (x - y)^2 (c x + d y) *)
FindEdgeSOSAnsatz[p_, q_] := Module[{x, y, a, b, c, d, poly, eqs},
  poly = a x^3 + b y^3 + p x^2 y + q x y^2 - (x - y)^2 (c x + d y);
  eqs = Thread[Flatten @ CoefficientList[Expand @ poly, {x, y}] == 0];
  Solve[eqs, {a, b, c, d}, Reals]
];

(* Prove: Sum[u_i (u_{i-1}^2 - 4 u_i^2 + 3 u_{i+1}^2)] <= 0 for u_i >= 0 *)
ProveExampleInequality[N_Integer?Positive] := Module[
  {u, left, right, S, sol, a, b, c, d, P, idOK, nonnegOK, concl},
  u = Array[u, N];
  left = RotateLeft[u];
  right = RotateRight[u];

  (* Target cyclic sum *)
  S = Total[u*(right^2 - 4 u^2 + 3 left^2)] // Expand;

  (* Edge-wise SOS coefficients for (u_i - u_{i+1})^2 (c u_i + d u_{i+1}) *)
  sol = First @ FindEdgeSOSAnsatz[-1, -3]; (* matches -u_i^2 u_{i+1} - 3 u_i u_{i+1}^2 *)
  {a, b, c, d} = {a, b, c, d} /. sol;     (* Expected: a=5/3, b=7/3, c=5/3, d=7/3 *)

  (* Nonnegative certificate *)
  P = Total[(u - left)^2*(c u + d left)] // Expand;

  (* Identity: P == (a + b) Sum u_i^3 - Sum u_i u_{i-1}^2 - 3 Sum u_i u_{i+1}^2 *)
  idOK = Simplify[
    Expand[P - ((a + b) Total[u^3] - Total[u*right^2] - 3 Total[u*left^2])]
  ] === 0;

  (* Nonnegativity and conclusion under u_i >= 0 *)
  nonnegOK = Simplify[P >= 0, Assumptions -> And @@ Thread[u >= 0]];
  concl = Simplify[S <= 0, Assumptions -> And @@ Thread[u >= 0]];

  <|
    "IdentityHolds" -> idOK,
    "EdgeWeights" -> <|"a" -> a, "b" -> b, "c" -> c, "d" -> d|>,
    "Certificate" -> P,
    "NonnegativityCertified" -> nonnegOK,
    "Conclusion" -> concl
  |>
];

End[]

EndPackage[]
