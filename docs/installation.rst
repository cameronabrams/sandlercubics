.. _installation:

Installation
============

Requirements
------------

sandlercubics requires:

* Python 3.7 or later
* NumPy
* SciPy
* sandlermisc (for some general utilities)
* sandlerprops (for compound property data, only for command-line usage)

Install from PyPI
-----------------

The easiest way to install sandlercubics is using pip:

.. code-block:: bash

   pip install sandlercubics

This will automatically install all required dependencies.

Install with Property Database
-------------------------------

To use the compound property database (recommended for most users):

.. code-block:: bash

   pip install sandlercubics sandlerprops

Install from Source
-------------------

To install the latest development version from GitHub:

.. code-block:: bash

   git clone https://github.com/cameronabrams/sandlercubics.git
   cd sandlercubics
   pip install -e .

This installs the package in "editable" mode, which is useful for development.

Verify Installation
-------------------

To verify that sandlercubics is installed correctly, run:

.. code-block:: bash

   sandlercubics --help

You should see the help message for the command-line interface.

Alternatively, from Python:

.. code-block:: python

   import sandlercubics
   print(sandlercubics.__version__)

Troubleshooting
---------------

Import Errors
~~~~~~~~~~~~~

If you encounter import errors, make sure that:

1. You have activated the correct Python environment
2. All dependencies are installed: ``pip install numpy scipy``
3. Your Python version is 3.7 or later: ``python --version``

NumPy/SciPy Installation Issues
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

On some systems, NumPy and SciPy may require additional system libraries. If pip installation fails:

**Linux (Debian/Ubuntu):**

.. code-block:: bash

   sudo apt-get install python3-numpy python3-scipy

**macOS (with Homebrew):**

.. code-block:: bash

   brew install openblas
   pip install numpy scipy

**Windows:**

Consider using Anaconda/Miniconda, which includes precompiled NumPy and SciPy:

.. code-block:: bash

   conda install numpy scipy
   pip install sandlercubics

Updating
--------

To update to the latest version:

.. code-block:: bash

   pip install --upgrade sandlercubics

Uninstalling
------------

To remove sandlercubics:

.. code-block:: bash

   pip uninstall sandlercubics
