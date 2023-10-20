# Cell annotation with decoupleR tutorial


This is a tutorial that can be used automatically annotate single-cell RNA-seq data by cell-type using decoupleR's python implementation.

To use, follow the steps below:

## Python requirements

- Python version >= 3.9

- Create virtual environment
   - for windows:
       ``````
       python -m venv venv

       venv\Scripts\activate

       pip install -r requirements.txt
       ``````

    - for linux:

        ``````
        python -m venv venv

        venv/bin/activate

        pip install -r requirements.txt
        ``````

    - for conda:

        ``````
        conda create -n <your_env_name>

        conda activate <your_env_name>

        conda install --file requirements.txt
        ``````

## How to run

- Open console.

- Activate virtual environment.

- Run the following command:

    ``````
    jupyter notebook
    ``````
- select the "decoupler_cell_annotation" notebook and run the example

### Acknowledgement

- This tutorial has been inspired by decoupleR's automatic cell annotation [vignette](https://decoupler-py.readthedocs.io/en/latest/notebooks/cell_annotation.html)