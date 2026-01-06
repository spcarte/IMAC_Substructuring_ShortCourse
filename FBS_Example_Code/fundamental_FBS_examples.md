# Fundamental FBS Examples
The section will go over fundamental substructuring examples that show different ways and programmatic methods for performing FBS coupling. It will go over three examples with the two DOF systems:
1. [Fully Manual One-way and Two-way FBS Coupling](manual_FBS)
2. [LM-FBS Coupling with Bespoke Functions](bespoke_LM-FBS_function)
3. [LM-FBS Coupling with SDynPy Functions](sdynpy_LM-FBS_function)

```{note}
Some of these examples will operate on the raw data from the associated `TransferFunctionArray` to demonstrate how the computations work. However, it is suggested that a software package, which automatically manages the bookkeeping for the FRFs be used to mitigate common sources of error in the substructuring process. 
```