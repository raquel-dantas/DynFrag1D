# DynFrag1D

A Python FEM code fora a dynamic fragmentation of a 1D bar by two approaches to modeling cracks:
- CZM by dynamic insertion of cohesive elements
- Diffuse damage by the Lip-field approach to fracture, follwing the methodology proposed in Moës, Lé, and Stershic  (2022).


The repository contains two src as follows:
- `src`: contains the code for using cohesive elements an the Lip-field.
- `src_akantu`: contains the code for the same problems but only using cohesive elements by benefiting from the FE open-source library Akantu (https://akantu.ch/) developed at LSMS.

The user has to source the setup file to install the required libraries and use the virtual environment to run the python files:
```bash
source setup.sh
```
Note: To use Akantu, the user needs to download, build and install Akantu. For this, follow the step-by-step in this link https://gitlab.com/akantu/akantu. Then, update the path to the python environment of Akantu in the `setup.sh`

## Usage

To run a simulation the user has to:

Settle the input values in a Python file and import it inside the `DFMesh.py` as `inputdata`. You can use the files inside the `input_files` directory to create a mesh, and random limit stress field.

Exemple:
```python
import my_input_file as inputdata
```
See the `input_data.py` file to see the required inputs and the format for it.

After settle the inputs just run the `main.py`.

The outputs files for each time-step are going to be save in a pickle file. You can use `verify_results.py` to read it.
