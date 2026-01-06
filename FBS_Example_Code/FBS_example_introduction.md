# IMAC Substructuring Shortcourse FBS Examples
This documentation will go through example problems to develop a deeper understanding of some of the concepts that were discussed in the IMAC XLIV substructuring short course. It will show examples that fit into the following themes:

1. Fundamental examples of how FBS can be performed
2. Examples demonstrating the significance of common errors

```{note}
While these examples leverage FBS for the substructure coupling, many of the lessons learned for the common errors can be applied to the other substructuring domains.
```

## Navigating this Book
Each example was built off an accompanying Jupyter notebook. The code and data in all these notebooks are self-contained (where the data is imported from the root directory of the file structure), so the code can be independently ran and modified. Note that some of the narrative content will build in a logical flow from one chapter to another. 

The package dependencies for running the code in the notebooks includes:
- Python3
- Jupyter Lab or Jupyter Notebook
- NumPy
- SciPy
- Matplotlib
- SDynPy (v0.21.0 or greater)

Note that these examples assume that the reader has a basic familiarity with Python and will not explain the basic commands. Some of the SDynPy specific commands will be explained to describe their use. Refer to the [SDynPy documentation](https://github.com/sandialabs/sdynpy) for more details on how to use the SDynPy functions. 

## A Note on Data Format
The data for these examples will be saved and imported as SDynPy `System` and `TransferFunctionArray` objects for simplicity of bookkeeping. However, the fundamental FBS examples will operate directly on the raw data from the objects to explicitly demonstrate how the code works. It is unnecessary to operate directly on the raw data in most practical cases and predefined functions in SDynPy will be used to simplify the coding in the examples for common substructuring errors.  

## Acknowledgements
Sandia National Laboratories is a multi-mission laboratory managed and operated by National Technology & Engineering Solutions of Sandia, LLC (NTESS), a wholly owned subsidiary of Honeywell International Inc., for the U.S. Department of Energy’s National Nuclear Security Administration (DOE/NNSA) under contract DE-NA0003525. This written work is authored by an employee of NTESS. The employee, not NTESS, owns the right, title and interest in and to the written work and is responsible for its contents. Any subjective views or opinions that might be expressed in the written work do not necessarily represent the views of the U.S. Government. The publisher acknowledges that the U.S. Government retains a non-exclusive, paid-up, irrevocable, world-wide license to publish or reproduce the published form of this written work or allow others to do so, for U.S. Government purposes. The DOE will provide public access to results of federally sponsored research in accordance with the DOE Public Access Plan. 
