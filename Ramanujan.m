(* ::Package:: *)

(* Ramanujan Machine Package for Mathematica *)
(* Translation of functions from https://github.com/RamanujanMachine/RamanujanMachine *)

BeginPackage["Ramanujan`"]

(* Public function declarations *)
MobiusTransform::usage = "MobiusTransform[{{a,b},{c,d}}] represents a Mobius transformation (ax+b)/(cx+d).";
ApplyMobiusTransform::usage = "ApplyMobiusTransform[transform, x] applies the Mobius transformation to x.";
ComposeMobiusTransforms::usage = "ComposeMobiusTransforms[t1, t2] composes two Mobius transformations.";
InverseMobiusTransform::usage = "InverseMobiusTransform[transform] returns the inverse Mobius transformation.";
ReciprocalMobiusTransform::usage = "ReciprocalMobiusTransform[transform] returns the reciprocal transformation.";

GeneralizedContinuedFraction::usage = "GeneralizedContinuedFraction[an, bn, depth] computes a generalized continued fraction with coefficients an and bn.";
EvaluateGCF::usage = "EvaluateGCF[an, bn, depth] evaluates a generalized continued fraction to the specified depth with default precision of 50 digits.";
SimpleContinuedFraction::usage = "SimpleContinuedFraction[constant, depth] computes simple continued fraction expansion.";

FindMobiusTransform::usage = "FindMobiusTransform[x, y, limit, threshold] finds integer Mobius transform T such that T(x) = y. The constraint a*x + b >= 1 avoids trivial solutions and removes redundancy.";

CreateLHSTable::usage = "CreateLHSTable[constants, searchRange] creates a lookup table of Mobius transforms of constants. Options: Threshold (default 10^-10), Precision (default 50).";

GetRamanujanConstants::usage = "GetRamanujanConstants[] returns a list of mathematical constants used in Ramanujan Machine.";

PolynomialSeries::usage = "PolynomialSeries[coefs, n] evaluates a polynomial series with given coefficients at index n.";
IterPolynomialSeries::usage = "IterPolynomialSeries[coefs, maxN] generates polynomial series values from 1 to maxN.";

Begin["`Private`"]

(* ============================================================================ *)
(* Constants *)
(* ============================================================================ *)

GetRamanujanConstants[] := Association[
    "e" -> E,
    "pi" -> Pi,
    "pi_squared" -> Pi^2,
    "catalan" -> Catalan,
    "golden_ratio" -> GoldenRatio,
    "euler_gamma" -> EulerGamma,
    "zeta" -> Zeta
];

(* ============================================================================ *)
(* Mobius Transform Functions *)
(* ============================================================================ *)

(* Normalize a Mobius transform by dividing by GCD of all elements *)
normalizeMobiusTransform[mat_List] := Module[{a, b, c, d, g},
    {{a, b}, {c, d}} = mat;
    g = GCD[a, b, c, d];
    If[g == 0 || g == 1, mat, mat / g]
];

(* Apply Mobius transformation to a value *)
ApplyMobiusTransform[{{a_, b_}, {c_, d_}}, x_] := (a*x + b)/(c*x + d);
ApplyMobiusTransform[{{a_, b_}, {c_, d_}}] := b/d; (* x = None case, returns constant term *)

(* Compose two Mobius transformations (matrix multiplication) *)
ComposeMobiusTransforms[t1_List, t2_List] := Module[{result},
    result = t1 . t2;
    normalizeMobiusTransform[result]
];

(* Inverse of a Mobius transform *)
InverseMobiusTransform[{{a_, b_}, {c_, d_}}] := Module[{det, result},
    det = a*d - b*c;
    result = det * {{d, -b}, {-c, a}};
    normalizeMobiusTransform[result]
];

(* Reciprocal transform: transforms T(x) to 1/T(x) *)
ReciprocalMobiusTransform[{{a_, b_}, {c_, d_}}] := {{c, d}, {a, b}};

(* Symbolic expression of Mobius transform *)
MobiusTransformExpression[{{a_, b_}, {c_, d_}}, x_] := (a*x + b)/(c*x + d);

(* ============================================================================ *)
(* Generalized Continued Fraction Functions *)
(* ============================================================================ *)

(* Efficient evaluation of generalized continued fraction using convergents *)
(* Uses recursive formula for convergents: https://en.wikipedia.org/wiki/Generalized_continued_fraction *)
EvaluateGCF[an_List, bn_List, depth_Integer, precision_Integer:50] := Module[
    {prevA, A, prevB, B, i, tmpA, tmpB, len},
    
    len = Min[depth, Length[an], Length[bn]];
    
    If[len == 0, Return[0]];
    
    (* Initialize convergents *)
    prevA = 0;
    A = 1;
    prevB = 1;
    B = an[[1]];
    
    (* Iterate through the continued fraction using convergent recurrence *)
    Do[
        tmpA = A;
        tmpB = B;
        A = an[[i]] * A + bn[[i]] * prevA;
        B = an[[i]] * B + bn[[i]] * prevB;
        prevA = tmpA;
        prevB = tmpB,
        {i, 2, len}
    ];
    
    (* Return the final convergent with specified precision (default 50 digits) *)
    If[A == 0, 0, N[B/A, precision]]
];

(* Generalized continued fraction with Mobius transform representation *)
GeneralizedContinuedFraction[an_List, bn_List, depth_Integer] := Module[
    {mobius, i, mat, a0},
    
    If[Length[an] == 0 || Length[bn] == 0, Return[0]];
    
    a0 = an[[1]];
    mobius = IdentityMatrix[2];
    
    (* Build the continued fraction using Mobius transforms *)
    Do[
        mat = {{0, bn[[i]]}, {1, an[[i + 1]]}};
        mobius = ComposeMobiusTransforms[mobius, mat],
        {i, 1, Min[depth - 1, Length[bn], Length[an] - 1]}
    ];
    
    (* Return a0 + b0/(continued fraction) *)
    a0 + ApplyMobiusTransform[mobius]
];

(* Simple continued fraction from a constant *)
SimpleContinuedFraction[const_, depth_Integer] := Module[
    {x, result, i, ai},
    
    x = N[const, 100];
    result = {};
    
    Do[
        ai = Floor[x];
        AppendTo[result, ai];
        If[x - ai == 0, Break[]];
        x = 1/(x - ai),
        {i, depth}
    ];
    
    result
];

(* ============================================================================ *)
(* Polynomial Series Functions *)
(* ============================================================================ *)

(* Evaluate polynomial series in compact form: n(n(...(a[0]*n + a[1]) + a[2]) + ...) + a[k] *)
PolynomialSeries[coefs_List, n_Integer] := Module[{tmp},
    tmp = 0;
    Do[tmp = tmp * n + coefs[[i]], {i, 1, Length[coefs]}];
    tmp
];

(* Generate polynomial series values *)
IterPolynomialSeries[coefs_List, maxN_Integer] := 
    Table[PolynomialSeries[coefs, n], {n, 1, maxN}];

(* ============================================================================ *)
(* LHS Hash Table Functions *)
(* ============================================================================ *)

(* Options for CreateLHSTable *)
Options[CreateLHSTable] = {
    "Threshold" -> 10^-10,  (* Decimal threshold for comparison *)
    "Precision" -> 50        (* Working precision for constant evaluation *)
};

(* Create lookup table of Mobius transforms of constants *)
CreateLHSTable[constants_List, searchRange_Integer, opts:OptionsPattern[]] := Module[
    {threshold, precision, nConstants, coefRange, numeratorCoefs, denominatorCoefs, 
     lhsTable, constVals, num, denom, value, key, keyFactor},
    
    threshold = OptionValue["Threshold"];
    precision = OptionValue["Precision"];
    
    (* Key factor determines hash precision: 1/threshold gives number of significant digits *)
    keyFactor = 1 / threshold;
    
    (* Evaluate constants to numerical values with specified precision *)
    constVals = {1.0} ~ Join ~ N[constants, precision];
    nConstants = Length[constVals];
    
    coefRange = Range[-searchRange, searchRange];
    lhsTable = <||>;
    
    (* Enumerate all coefficient combinations *)
    numeratorCoefs = Tuples[coefRange, nConstants];
    denominatorCoefs = Tuples[coefRange, nConstants];
    
    Monitor[
        Do[
            num = numCoef . constVals;
            
            (* Only consider positive numerators to avoid duplication (T and -T give same ratio) *)
            If[num > 0,
                Do[
                    denom = denomCoef . constVals;
                    
                    (* Skip if denominator is zero *)
                    If[denom != 0,
                        value = N[num/denom, 30];
                        
                        (* Create hash key from first significant digits based on threshold *)
                        key = Floor[value * keyFactor];
                        
                        (* Store the coefficients that produce this value *)
                        If[!KeyExistsQ[lhsTable, key],
                            lhsTable[key] = {}
                        ];
                        AppendTo[lhsTable[key], {numCoef, denomCoef}]
                    ],
                    {denomCoef, denominatorCoefs}
                ]
            ],
            {numCoef, numeratorCoefs}
        ],
        Row[{"Building LHS table: ", Length[Keys[lhsTable]], " entries"}]
    ];
    
    lhsTable
];

(* ============================================================================ *)
(* Find Mobius Transform Functions *)
(* ============================================================================ *)

(* Find integer Mobius transform such that T(x) = y within given coefficient limits *)
FindMobiusTransform[x_?NumericQ, y_?NumericQ, limit_Integer, threshold_:10^-7] := Module[
    {a, b, c, d, eq, sol, minError, bestSol},
    
    (* The equation is: ax + b - cxy - dy = 0 *)
    (* We want to find integer solutions a, b, c, d within [-limit, limit] *)
    
    minError = Infinity;
    bestSol = None;
    
    (* Brute force search for small limits *)
    Do[
        Do[
            Do[
                Do[
                    (* Check the equation error and avoid trivial/redundant solutions *)
                    (* The constraint a*x + b >= 1 prevents trivial solutions and reduces redundancy *)
                    If[Abs[a*x + b - c*x*y - d*y] < threshold && a*x + b >= 1,
                        If[Abs[a*x + b - c*x*y - d*y] < minError,
                            minError = Abs[a*x + b - c*x*y - d*y];
                            bestSol = {{a, b}, {c, d}};
                        ]
                    ],
                    {d, -limit, limit}
                ],
                {c, -limit, limit}
            ],
            {b, -limit, limit}
        ],
        {a, -limit, limit}
    ];
    
    If[bestSol =!= None,
        normalizeMobiusTransform[bestSol],
        None
    ]
];

(* ============================================================================ *)
(* Utility Functions *)
(* ============================================================================ *)

(* Pretty print a Mobius transform *)
PrintMobiusTransform[{{a_, b_}, {c_, d_}}, x_:x] := Module[{},
    Print[TraditionalForm[(a*x + b)/(c*x + d)]]
];

(* Print generalized continued fraction in readable form *)
PrintGCF[an_List, bn_List, depth_Integer:3] := Module[{expr, x},
    x = Symbol["..."];
    expr = x;
    Do[
        expr = an[[i]] + bn[[i]]/expr,
        {i, depth, 1, -1}
    ];
    Print[TraditionalForm[expr]]
];

(* Compute convergence rate of a GCF *)
ConvergenceRate[an_List, bn_List, maxDepth_Integer:100] := Module[
    {values, i},
    
    values = Table[
        EvaluateGCF[an[[1;;i]], bn[[1;;i]], i],
        {i, 10, Min[maxDepth, Length[an], Length[bn]]}
    ];
    
    ListPlot[values, 
        PlotLabel -> "GCF Convergence",
        AxesLabel -> {"Depth", "Value"},
        PlotStyle -> PointSize[Medium],
        Joined -> True
    ]
];

End[]

EndPackage[]

(* Usage example:
   
   (* Load the package *)
   Get["Ramanujan.m"]
   
   (* Get constants *)
   consts = Ramanujan`GetRamanujanConstants[]
   
   (* Compute simple continued fraction of Pi *)
   cfPi = Ramanujan`SimpleContinuedFraction[Pi, 10]
   
   (* Evaluate a generalized continued fraction *)
   an = {1, 1, 1, 1, 1, 1, 1, 1};
   bn = {1, 1, 1, 1, 1, 1, 1, 1};
   Ramanujan`EvaluateGCF[an, bn, 8]
   
   (* Create Mobius transform *)
   t1 = {{1, 2}, {3, 4}};
   Ramanujan`ApplyMobiusTransform[t1, 5]
*)
