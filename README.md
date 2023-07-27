ConTesse: Accurate Occluding Contours for Subdivision Surfaces
==============================================================

<strong>Chenxi Liu<sup>1</sup>, Pierre Bénard<sup>2</sup>, Aaron Hertzmann<sup>3</sup>, Shayan Hoshyari<sup>4</sup></strong>

<small><sup>1</sup>University of British Columbia, <sup>2</sup>Univ.Bordeaux, CNRS, Bordeaux INP, INRIA, LaBRI, UMR 5800, <sup>3</sup>Adobe Research, <sup>4</sup>Adobe</small>

<img src="https://dgp.toronto.edu/~hertzman/contesse/representative.jpg"/>

This repository contains the source code for the paper
[ConTesse: Accurate Occluding Contours for Subdivision Surfaces](https://dgp.toronto.edu/~hertzman/contesse/).

 - `src` contains the source code of ConTesse.
 - `py` contains Python helper scripts to run on given models.
 - `cmake` contains cmake files to fetch external libraries.
 - `external` contains external codes.

The data (configuration files, camera files, obj models) is included in `data`; test models from [Computing Smooth Surface Contours with Accurate Topology](https://inria.hal.science/hal-00924273) can be downloaded from
[here](https://dgp.toronto.edu/~hertzman/contesse/contesse-data.zip).

Building
--------

1. Setup contess development env by running
    ```
    conda env update --file conda_contesse.yaml
    conda activate contess
    ```
2. Build the code
    ```bash
    mkdir build && cd build
    cmake .. -DCMAKE_BUILD_TYPE=Release
    # Build using the selected generator
    ```

Running our method
--------

### Run individual test cases

An individual test case is defined by a given `obj` model and a camera setting in a `json` file. For instance,
```
{
  "camera_info": {
    "index": 1,
    "value": {
      "camera_x": 0.0,
      "camera_y": 0.0,
      "camera_z": 0.0,
      "full_info": false,
      "lookat": [
        0.0,
        0.0,
        0.0
      ]
    }
  },
  "description": "Opensubdiv torus quad view 1: subdivided",
  "model_name": "torus_quad_v1",
  "model_path": "torus_quad"
  "subdiv_level": 1
}
```
Here `model_path` is relative to `data/meshes` folder in this repo. Configuration `json` files of input models used in our paper can be found in `data/json` (`pipeline` contains static models; `bunny-p1`, `bunny-p2`, `walking`, `Angela` contains test sequences from [Bénard et al. 2014]).

To run an individual test case, call, for instance,
```bash
# From build directory
python3 ../py/scripts/run_model.py ../data/json/pipeline/torus_quad_v1.json ./TestsContourTessellation --out [output_folder]
```
The resulting files `torus_quad_v1_out.svg` (ConTesse output), `torus_quad_v1_trivial.svg` (comparison output by trivial ray casting visibility test) are generated into `[output_folder]`.

### Run test cases in batch

Multiple test cases defined by individual `json` files can be run in batch mode. To define a batch of inputs, use a `json` file, for instance,
```
{
  "multicamera_test_cases": [
    "pipeline/split.json",
    "pipeline/torus_quad_v1.json"
  ],
  "pipeline_test_cases": [
    "pipeline/ribbon_open.json"
  ]
}
```
Here `multicamera_test_cases` automatically test 26 camera positions centered to the model in a grid; `pipeline_test_cases` uses the camera setting given in the individual `json` files.

To run multiple test cases in batch mode, call, for instance,
```bash
# From build directory
python3 ../py/scripts/run_model.py [batch_json] ./TestsContourTessellation --max-procs [max_number_of_processes] --out [output_folder]
```
The resulting files are generated with maximum `[max_number_of_processes]` processes into `[output_folder]`.

License
-------

The source code (everything under `src` and `py`) is licensed under [Version 2.0 of the Apache License](LICENSE.md).
The input models (which need to be downloaded separately) are licensed under separate licenses. Please refer to `license.yml` for license information for each model.


BibTeX
------

```
@article{ConTesse,
  title = {ConTesse: Accurate Occluding Contours for Subdivision Surfaces},
  author = {Liu, Chenxi and B\'{e}nard, Pierre and Hertzmann, Aaron and Hoshyari, Shayan},
  year = {2023},
  journal = {ACM Trans. Graph.},
  month = feb,
  volume = {42},
  number = {1},
  articleno = {5},
  url = {https://doi.org/10.1145/3544778},
  doi = {10.1145/3544778}
}
```