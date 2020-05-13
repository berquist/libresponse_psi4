# libresponse_psi4

A Psi4 plugin that interfaces to [libresponse](https://github.com/LambrechtLab/libresponse), a generalized molecular property solver.

## Build

1. Get a working Psi4 install. The easiest way is through Anaconda:
```bash
conda create -n p4env python=3.7 psi4 psi4-dev -c psi4 -c psi4/label/dev
conda activate p4env
```
2. Install additional libresponse dependencies. Currently, this is only [Armadillo](http://arma.sourceforge.net/). There is no hard restriction on the version, as long as everything compiles later on:
```bash
conda install armadillo
```
3. Get a copy of this plugin and libresponse, always using the latest version of both (from their `master` branches):
```bash
git clone https://github.com/berquist/libresponse_psi4.git
cd libresponse_psi4
git clone -b psi4 https://github.com/LambrechtLab/libresponse.git
```
4. Prepare the plugin for compilation, then compile.
```bash
`psi4 --plugin-compile` -DCMAKE_EXPORT_COMPILE_COMMANDS=1
make
```
5. Run the example input:
```bash
psi4 input.dat
```
