# Developability of Heightfields
Public code release for ["Developability of Heightfields via Rank Minimization"](http://dgp.toronto.edu/~sgsellan/pdf/compressed-developables.pdf), presented at SIGGRAPH 2020 and authored by [Silvia Sellán](http://dgp.toronto.edu/~sgsellan/), [Noam Aigerman](https://noamaig.github.io) and [Alec Jacobson](http://www.cs.toronto.edu/~jacobson/). Please note that while this Matlab implementation of our method is hereby released under MIT License, the method itself is pending a US patent filed in 2020 by Adobe Inc.



## Installation
We have only succesfully validated this installation on MacOS 10.15 with Matlab 2020a, but are confident it should work on most Unix-based machines and versions of Matlab from the past five years. 

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
 
We also allow you to replicate the results from our paper exactly, by runnning the scripts in the `figures/` directory. Instructions to run and understand the output of each can be found in each scripts' first commented lines. In most cases, the scripts will generate input and output `.obj` files which we then rendered using Blender for the paper figures. You can look at our Blender setup in `render/` and substitute the existing meshes with the input and output from our scripts to exactly replicate our paper figures up to very minor lighting direction and orientation choices.

## C++ implementation

This Matlab implementation is the one we used to generate all the examples in the paper, and it has been thouroughly tested. I strongly encourage you to use it. However, if you do not have access to Maltab for whichever reason and running our code would be useful for you, I also provide a limited, guarantee-free C++ implementation of our algorithm [here](https://github.com/sgsellan/developability-of-heightfields-cpp).  

## Known Issues
Please do not hesitate to contact
[sgsellan@cs.toronto.edu](mailto:sgsellan@cs.toronto.edu) if you find any issues
or bugs in this code, or you struggle to run it in any way.



