# Contributing to MNiShed

## Reporting bugs and requesting features

Please open an issue on [GitHub](https://github.com/MNiMORPH/MNiShed/issues).
Include a minimal reproducible example and your Python and package versions where relevant.

## Contributing code

1. Fork the repository and create a branch from `master`.
2. Install the package in editable mode with documentation dependencies:
   ```
   pip install -e ".[docs]"
   ```
3. Make your changes. Run the Cannon River example to confirm the model still produces a sensible result:
   ```
   cd examples/cannon_forward
   python run_forward.py
   ```
4. If you changed the documentation, build it locally to check for errors:
   ```
   cd docs
   make html
   ```
5. Open a pull request against `master` with a clear description of what changed and why.

## Code style

This project follows [PEP 8](https://peps.python.org/pep-0008/) and the
[CSDMS developer guidelines](https://csdms.colorado.edu/wiki/Development_best_practices),
with the intentional exceptions listed below.

You can check for unintentional style deviations with:
```
pycodestyle --max-line-length=99 mnished/
```

Write [NumPy-style docstrings](https://numpydoc.readthedocs.io/en/latest/format.html)
for all public methods.

### Intentional PEP 8 exceptions

| Code | Reason |
|------|--------|
| E741 | Single-letter scientific variables are permitted where conventional (e.g., `I` for the Thornthwaite thermal index). |
| E221 | Extra spaces before `=` are allowed for vertical alignment of related assignments. |
| E225 | Whitespace around `+` may be omitted inside multi-line string literals for readability. |
| E251 | Spaces around `=` in keyword arguments are allowed when aligning a multi-line constructor call. |

### CSDMS BMI interface

The `Buckets` class implements the
[CSDMS Basic Model Interface](https://bmi.readthedocs.io/).
Changes to `initialize`, `update`, `run`, or `finalize` must preserve
this interface contract.

## Contact

Andrew D. Wickert — [awickert@umn.edu](mailto:awickert@umn.edu) — [ORCID 0000-0002-9545-3365](https://orcid.org/0000-0002-9545-3365)
