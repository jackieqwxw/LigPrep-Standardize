# LigPrep-Standardize
This tool assists in standardizing compounds during the ligand preparation stage before virtual screening.

## Workflow
Install RDKit toolkit on Linux using conda:

```conda create --name rdkit-tools python=3.8```

Activate the RDKit enviroment:

```conda activate rdkit-tools```

Run the script:
```python3 standardize.py inputfilename.smi outputfilename.smi```
