# Developability of Heightfields
Public code release for ["Developability of Heightfields via Rank Minimization"](http://dgp.toronto.edu/~sgsellan/pdf/compressed-developables.pdf), presented at SIGGRAPH 2020 and authored by [Silvia Sellán](http://dgp.toronto.edu/~sgsellan/), [Noam Aigerman](https://noamaig.github.io) and [Alec Jacobson](http://www.cs.toronto.edu/~jacobson/). Please note that while this Matlab implementation of our method is hereby released under MIT License, the method itself is pending a US patent filed in 2020 by Adobe Inc.



## Installation
To install this library, please start by cloning the repository recursively
```
git clone --recursive https://github.com/sgsellan/developability-of-heightfields.git
```
After this, we will build the mex functions in the `gptoolbox` directory:
```
cd developability-of-heightfields/gptoolbox/mex
mkdir build
cd build
cmake ..
make
```
Then, build our own mex files by entering Matlab and, in Matlab, adding this repository in its entirety to your Matlab path (for instance, by running `addpath(genpath(path/to/developability-of-heightfields))`), navigating to `developability-of-heightfields/mex` and running `setup_mex` in the Matlab console.

## Use
To use our code in the exact same way we did for most our paper results, run in
Matlab
```
gui_developables(path/to/some/mesh)
```
and follow the instructions in the console. For example, you can start by trying
```
gui_developables('data/bunny.obj')
```
 
We also allow you to replicate the results from our developable interpolation
results (Figs. 18 & 19). Simply look at `interpolation_grayscale.m`, follow the
instructions to use your own constraints and link to the image in Line 14. You
can also run it as is to recover our output from Fig. 18.

## Known Issues
Please do not hesitate to contact
[sgsellan@cs.toronto.edu](mailto:sgsellan@cs.toronto.edu) if you find any issues
or bugs in this code, or you struggle to run it in any way.



