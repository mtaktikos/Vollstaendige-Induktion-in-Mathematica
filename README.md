# Vollstaendige-Induktion-in-Mathematica
Versuch der Automatisierung der Vollständigen Induktion mit Mathematica

## Packages

### Ramanujan.m
Translation of functions from the [Ramanujan Machine](https://github.com/RamanujanMachine/RamanujanMachine) to Mathematica. The Ramanujan Machine is an algorithmic approach to discover new mathematical conjectures, particularly focused on finding formulas relating fundamental constants like Pi, E, and the Riemann zeta function to various continued fractions.

**Features:**
- Mobius transformations
- Generalized continued fractions (GCF)
- Simple continued fraction expansions
- Polynomial series generation
- LHS hash tables for finding constant relationships
- Mathematical constants (Pi, E, GoldenRatio, Catalan, etc.)

**Quick Start:**
```mathematica
(* Load the package *)
<<Ramanujan`

(* Simple continued fraction of Pi *)
cfPi = Ramanujan`SimpleContinuedFraction[Pi, 10]

(* Evaluate generalized continued fraction *)
an = Table[1, {10}]; bn = Table[1, {10}];
Ramanujan`EvaluateGCF[an, bn, 10]  (* Returns GoldenRatio *)
```

See `RamanujanExamples.nb` for comprehensive examples and `Ramanujan_Usage.txt` for full documentation.

### CyclicInequalities.m
Tools for working with cyclic inequalities.
