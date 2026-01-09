# sandlercubics Documentation

This directory contains the source files for the sandlercubics documentation, built with [Sphinx](https://www.sphinx-doc.org/).

## Building the Documentation

### Prerequisites

Install the documentation requirements:

```bash
pip install -r requirements.txt
```

Or install with the package:

```bash
cd ..
pip install -e ".[docs]"
```

### Building HTML Documentation

On Linux/macOS:

```bash
make html
```

On Windows:

```bash
make.bat html
```

The built documentation will be in `_build/html/`. Open `_build/html/index.html` in a web browser to view.

### Other Build Formats

```bash
make latexpdf  # Build PDF (requires LaTeX)
make epub      # Build EPUB
make man       # Build man pages
make help      # Show all available formats
```

### Cleaning Build Files

```bash
make clean
```

## Documentation Structure

```
docs/
├── conf.py              # Sphinx configuration
├── index.rst            # Main documentation page
├── installation.rst     # Installation guide
├── quickstart.rst       # Quick start guide
├── cli.rst              # CLI documentation
├── examples.rst         # Usage examples
├── theory.rst           # Theoretical background
├── contributing.rst     # Contributing guide
├── changelog.rst        # Version history
├── api/                 # API reference
│   ├── index.rst        # API overview
│   ├── eos.rst          # EOS module docs
│   ├── cli.rst          # CLI module docs
│   └── utils.rst        # Utils module docs
├── Makefile             # Build script (Unix)
├── make.bat             # Build script (Windows)
└── requirements.txt     # Documentation dependencies
```

## Writing Documentation

### reStructuredText Syntax

The documentation uses reStructuredText (reST) format. Key syntax:

**Headers:**
```rst
Page Title
==========

Section
-------

Subsection
~~~~~~~~~~
```

**Code blocks:**
```rst
.. code-block:: python

   from sandlercubics import PengRobinsonEOS
   eos = PengRobinsonEOS(Tc=190.4, Pc=4.6, omega=0.011)
```

**Links:**
```rst
:doc:`cli`           # Link to another doc
:ref:`api_eos`       # Link to a label
`text <url>`_        # External link
```

**Math:**
```rst
.. math::

   P = \frac{RT}{v-b} - \frac{a}{v^2}
```

**Autodoc (from docstrings):**
```rst
.. autoclass:: sandlercubics.eos.PengRobinsonEOS
   :members:
   :show-inheritance:
```

### Adding New Pages

1. Create a new `.rst` file in the appropriate location
2. Add it to a `.. toctree::` directive in the parent page
3. Build and verify

### Style Guidelines

- Use clear, concise language
- Include code examples for API documentation
- Cross-reference related sections
- Test all code examples
- Keep line length reasonable (~80 chars for prose)

## Updating Documentation

After modifying docstrings in the source code:

1. Rebuild the documentation: `make html`
2. Verify changes in browser
3. Commit both source changes and documentation updates

## Troubleshooting

**"sphinx-build command not found"**
- Install sphinx: `pip install sphinx`

**Import errors during build**
- Install the package: `pip install -e ..`
- Or install dependencies: `pip install -r requirements.txt`

**Autodoc not finding modules**
- Check `sys.path` configuration in `conf.py`
- Verify package is installed: `pip list | grep sandlercubics`

**Build warnings**
- Fix warnings before committing
- Some warnings are expected initially if source code isn't present

## Live Preview (Optional)

For automatic rebuilding while editing:

```bash
pip install sphinx-autobuild
sphinx-autobuild . _build/html
```

Then open http://127.0.0.1:8000 in your browser.

## Additional Resources

- [Sphinx Documentation](https://www.sphinx-doc.org/)
- [reStructuredText Primer](https://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html)
- [Napoleon Extension](https://www.sphinx-doc.org/en/master/usage/extensions/napoleon.html) (NumPy/Google docstrings)
- [Read the Docs Theme](https://sphinx-rtd-theme.readthedocs.io/)
