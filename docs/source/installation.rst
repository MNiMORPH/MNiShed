Installation
=============

Requirements
~~~~~~~~~~~~

* Python 3.9 or later
* NumPy
* Pandas
* PyYAML
* Matplotlib (for plotting)
* SciPy

These are installed automatically with the package.

Optional acceleration (Numba)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

MNiShed's daily time loop has a `Numba <https://numba.pydata.org>`_
just-in-time (JIT) compiled implementation that runs roughly two orders of
magnitude faster than the pure-Python loop.  This is a large saving for
calibration, where the model is run thousands of times.  Numba is *not* a
core dependency; install it with the ``jit`` extra, which pulls in Numba
together with a compatible NumPy:

.. code-block:: bash

    pip install mnished[jit]

The extra caps NumPy at ``<2.3`` (the newest release the current Numba
supports); a bare ``pip install numba`` also works in an environment that
already has ``numpy < 2.3``.  Without a compatible NumPy, ``import numba``
fails; MNiShed then uses the pure-Python loop and emits a one-time warning so
the slowdown is not silent.

When Numba is present the JIT loop is used automatically for every supported
configuration (including PDM saturation-excess and ``et_water_stress``); when
Numba is absent, MNiShed falls back to the pure-Python loop, with identical
results.

From PyPI
~~~~~~~~~

The easiest way to install MNiShed is via pip:

.. code-block:: bash

    pip install mnished

From Source (Development)
~~~~~~~~~~~~~~~~~~~~~~~~~

For development or contributing to the project:

.. code-block:: bash

    git clone https://github.com/MNiMORPH/MNiShed.git
    cd MNiShed
    pip install -e '.[bmi,jit]'   # optional extras: bmi, jit, docs, lint

This installs the package in "editable" mode, so changes to the source code are
reflected immediately without reinstalling.  The optional extras are ``bmi``
(CSDMS BMI wrapper), ``jit`` (Numba acceleration), ``docs`` (build the
documentation), and ``lint`` (Ruff).

Installing Documentation Dependencies
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To build documentation locally:

.. code-block:: bash

    pip install -r docs/requirements.txt
    cd docs
    make html

The HTML documentation will be in ``docs/_build/html/``.

Verification
~~~~~~~~~~~~

To verify your installation:

.. code-block:: python

    >>> import mnished
    >>> print(mnished.__version__)
    
    >>> # Create and initialize a model
    >>> model = mnished.Buckets()
    >>> print("Installation successful!")

Troubleshooting
~~~~~~~~~~~~~~~

**ImportError when importing mnished**

  Ensure you've run ``pip install -e .`` from the repository root.

**YAML parsing errors**

  Update PyYAML: ``pip install --upgrade pyyaml``

**Plotting doesn't work**

  Install Matplotlib: ``pip install matplotlib``

Next Steps
~~~~~~~~~~

Head to the :doc:`quickstart` guide to run your first model!
