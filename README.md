# Data reduction for neutron depth profiling measurements

## What is in this git:
1. Reduce.py - Reduces neutron depth profiling data from a single sample according to a schema.json file
2. Schema.py - Creates the schema.json file used to control the flow of data processing in Reduce.py
3. config.json - instrument configuration file in JSON format, read by Reduce.py
4. Jupyter/ndpReduce.ipynb - Jupyter notebook interface for reducing ndp data with Reduce.py
5. Jupyter/schema.ipynb - Jupyter notebook interface for creating schema files

## To install:
1. Download Anaconda (there are other options, including conda and pip, but hard to beat Anaconda for ease)
2. Create an environment for ndp within Anaconda
3. Within that environment, install Python 3.x (im using 3.9, but earlier versions probably will work)
4. Install latest numpy and matplotlib libraries
5. Clone this repository to your own machine (using the green "code" menu on the main page of this git)

## To run this code:
1. Copy both Jupyter notebooks into the lowest level directory that contains all of the data to be reduced. 
2. Run Jupyter (anaconda is great for this, run it in the ndp environment)
3. Modify the schema cells to have all of the info specific to your data sets.
3. Modify ndpReduce.ipynb to load your new schema file
4. Modify ndpReduce.ipynb to reflect the directory where you are keeping Reduce.py
5. Run all the cells
6. Adjust the plots to your preferences

