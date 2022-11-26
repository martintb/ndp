# Data reduction for neutron depth profiling measurements

## What is in this git:
1. src/ndp/reduce.py - Reduces neutron depth profiling data from a single sample according to a schema.json file
2. src/ndp/schema.py - Creates the schema.json file used to control the flow of data processing in Reduce.py
3. src/ndp/example_files/config.json - instrument configuration file in JSON format, read by Reduce.py
4. src/ndp/jupyter/ndpReduce.ipynb - Jupyter notebook interface for reducing ndp data with Reduce.py
5. src/ndp/jupyter/schema.ipynb - Jupyter notebook interface for creating schema files

## To install:
1. Install Python>=3.9 via anaconda.com (recommended) or python.org
2. (optional) Create an environment for ndp within Anaconda or using virtualenv
3. Open the Anaconda prompt or your favorite command line emulator and run:
```bash
$ pip install ndp
```

## To run this code:
1. Navigate in the terminal to your working directory and run
```bash
$ python -m ndp.get_notebooks
```
2. Start a jupyter server from this directory
```bash
$ jupyter-lab
```
3. Modify schema.ipynb to have all of the info specific to your data sets
4. Run schema.ipynb to create your schema file
3. Modify ndpReduce.ipynb to load your new schema file
4. Modify ndpReduce.ipynb to reflect the directory where you are keeping Reduce.py
5. Run all the cells
6. Adjust the plots to your preferences

