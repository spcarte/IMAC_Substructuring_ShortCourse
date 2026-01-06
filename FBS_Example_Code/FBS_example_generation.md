# Creating the Example Systems
This section will develop the SDynPy `System` and `TransferFunctionArray` objects that will be used in the example problems. These objects define the fundamental information for the different example systems:

- SDynPy `System` object - this object type stores the mass, stiffness, and damping matrices for the different systems with an accompanying DOF list (referred to as a `CoordinateArray` in SDynPy), which simplifies the bookkeeping 
- SDynPy `TransferFunctionArray` object - this object stores the FRF ordinate and abscissa quantities for the different systems with the accompanying reference and response DOF lists 

The examples leverage two different systems:

1. [Two DOF system](generate_2DOF_example_system) - The two DOF example system is used for relatively simple examples, where it is a convenient to be able to look at all the quantities in the FRF matrix.
2. [Beam system](generate_beam_example_system) - the example beam system is used when it is convenient to have a more complex system to demonstrate some of the practical problems that are experienced in substructuring. The models for the beam system have a significant number of nodes, which makes it more difficult to analyze the substructuring results compared to the two DOF system. 