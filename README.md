# README

A simplistic python version of the [scHCL R package ](https://github.com/ggjlab/scHCL),
which maps query cells onto the closest representative in the HCL dataset.
The main use case is cell type identification!

This is currently in a pretty rough state, here's some usage:

```python
from scHCLpy.main import load_reference, scHCL

ref=load_reference()
scHCL(query_df)
```
