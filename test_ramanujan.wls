#!/usr/bin/env wolframscript

(* Load the package *)
Get["Ramanujan.m"]

Print["Testing Ramanujan Machine Package"]
Print["===================================="]

(* Test 1: Get Constants *)
Print["\nTest 1: Get Constants"]
consts = Ramanujan`GetRamanujanConstants[];
Print["Available constants: ", Keys[consts]]

(* Test 2: Simple Continued Fraction *)
Print["\nTest 2: Simple Continued Fraction of Pi"]
cfPi = Ramanujan`SimpleContinuedFraction[Pi, 10];
Print["CF(Pi): ", cfPi]

(* Test 3: Simple Continued Fraction of E *)
Print["\nTest 3: Simple Continued Fraction of E"]
cfE = Ramanujan`SimpleContinuedFraction[E, 10];
Print["CF(E): ", cfE]

(* Test 4: Golden Ratio from GCF *)
Print["\nTest 4: Golden Ratio from GCF with all 1s"]
an = Table[1, {10}];
bn = Table[1, {10}];
gcfGolden = Ramanujan`EvaluateGCF[an, bn, 10];
Print["GCF result: ", gcfGolden]
Print["Golden Ratio: ", N[GoldenRatio, 20]]
Print["Difference: ", Abs[gcfGolden - N[GoldenRatio, 20]]]

(* Test 5: Mobius Transform *)
Print["\nTest 5: Mobius Transform"]
t1 = {{1, 2}, {3, 4}};
result = Ramanujan`ApplyMobiusTransform[t1, 5];
Print["T(5) where T(x) = (x+2)/(3x+4): ", result]
Print["Expected: ", N[(5+2)/(3*5+4)]]

(* Test 6: Compose Transforms *)
Print["\nTest 6: Compose Transforms"]
t2 = {{2, 1}, {1, 3}};
t3 = Ramanujan`ComposeMobiusTransforms[t1, t2];
r1 = Ramanujan`ApplyMobiusTransform[t3, 5];
r2 = Ramanujan`ApplyMobiusTransform[t1, Ramanujan`ApplyMobiusTransform[t2, 5]];
Print["T1(T2(5)): ", r2]
Print["T3(5): ", r1]
Print["Equal? ", r1 == r2]

(* Test 7: Inverse Transform *)
Print["\nTest 7: Inverse Transform"]
tInv = Ramanujan`InverseMobiusTransform[t1];
x = 7;
forward = Ramanujan`ApplyMobiusTransform[t1, x];
backward = Ramanujan`ApplyMobiusTransform[tInv, forward];
Print["x: ", x]
Print["T(x): ", forward]
Print["T_inv(T(x)): ", backward]
Print["Recovered x? ", Abs[x - backward] < 0.0001]

(* Test 8: Polynomial Series *)
Print["\nTest 8: Polynomial Series"]
coefs = {1, 2, 3};
result = Ramanujan`PolynomialSeries[coefs, 5];
expected = 5^2 + 2*5 + 3;
Print["PolynomialSeries[{1,2,3}, 5]: ", result]
Print["Expected (5^2 + 2*5 + 3): ", expected]

(* Test 9: E from GCF *)
Print["\nTest 9: Approximate E-1 using GCF"]
anE = Table[n, {n, 1, 20}];
bnE = Table[n, {n, 1, 20}];
gcfE = Ramanujan`EvaluateGCF[anE, bnE, 20];
Print["GCF result: ", gcfE]
Print["E-1: ", N[E-1, 20]]
Print["Difference: ", Abs[gcfE - N[E-1, 20]]]

Print["\n===================================="]
Print["All tests completed!"]
