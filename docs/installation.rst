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
* sandlerprops (for compound critical properties and heat-capacity data)

Install from PyPI
-----------------

The easiest way to install sandlercubics is using pip:

.. code-block:: bash

   pip install sandlercubics

This will automatically install all required dependencies.

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
