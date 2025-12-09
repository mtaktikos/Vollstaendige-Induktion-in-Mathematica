#!/usr/bin/env wolframscript
(* Quick Start Guide for Ramanujan Machine Package *)

(* Load the package *)
Get["Ramanujan.m"]

Print["Ramanujan Machine - Quick Start Examples"]
Print["==========================================\n"]

(* Example 1: Simple Continued Fraction of Pi *)
Print["1. Simple Continued Fraction of Pi (first 10 terms):"]
cfPi = Ramanujan`SimpleContinuedFraction[Pi, 10];
Print[cfPi]
Print[]

(* Example 2: Golden Ratio from GCF with all 1s *)
Print["2. Golden Ratio from Generalized CF with all 1s:"]
an = Table[1, {10}];
bn = Table[1, {10}];
gcfGolden = Ramanujan`EvaluateGCF[an, bn, 10];
Print["GCF result: ", N[gcfGolden, 20]]
Print["Golden Ratio: ", N[GoldenRatio, 20]]
Print[]

(* Example 3: Approximate E-1 using GCF *)
Print["3. Approximate E-1 using GCF where an = bn = n:"]
anE = Table[n, {n, 1, 15}];
bnE = Table[n, {n, 1, 15}];
gcfE = Ramanujan`EvaluateGCF[anE, bnE, 15];
Print["GCF result: ", N[gcfE, 20]]
Print["E-1: ", N[E-1, 20]]
Print[]

(* Example 4: Mobius Transform *)
Print["4. Mobius Transform T(x) = (x+2)/(3x+4):"]
t = {{1, 2}, {3, 4}};
result = Ramanujan`ApplyMobiusTransform[t, 5];
Print["T(5) = ", result, " = ", N[result]]
Print[]

(* Example 5: Polynomial Series *)
Print["5. Polynomial Series n^2 + 2n + 3:"]
coefs = {1, 2, 3};
series = Ramanujan`IterPolynomialSeries[coefs, 10];
Print[series]
Print[]

(* Example 6: Mathematical Constants *)
Print["6. Available Mathematical Constants:"]
consts = Ramanujan`GetRamanujanConstants[];
Print[Keys[consts]]
Print[]

Print["=========================================="]
Print["For more examples, see RamanujanExamples.nb"]
Print["For documentation, see Ramanujan_Usage.txt"]
