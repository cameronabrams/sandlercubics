# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys

# Add the package to the path for autodoc
sys.path.insert(0, os.path.abspath('../src'))

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'sandlercubics'
copyright = '2025-2026, Cameron F. Abrams'
author = 'Cameron F. Abrams'
# get version from the package

release = sandlercubics.__version__
version = release  # short X.Y.Z version

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.intersphinx',
    'sphinx.ext.mathjax',
    'sphinx.ext.napoleon',
    'sphinx.ext.todo',
    'sphinx.ext.viewcode',
    'sphinx_copybutton',  # Adds copy button to code blocks
    'sphinxcontrib.mermaid',
]

# Napoleon settings for NumPy/Google style docstrings
napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = True
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = True
napoleon_use_admonition_for_examples = False
napoleon_use_admonition_for_notes = False
napoleon_use_admonition_for_references = False
napoleon_use_ivar = False
napoleon_use_param = True
napoleon_use_rtype = True
napoleon_preprocess_types = False
napoleon_type_aliases = None
napoleon_attr_annotations = True

# Autosummary settings
autosummary_generate = True
autosummary_imported_members = False

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# The master toctree document.
master_doc = 'index'

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'furo' #sphinx_rtd_theme'  # or 'alabaster', 'furo', etc.
html_static_path = ['_static']

# Theme options
# In conf.py

html_theme = 'furo'

html_theme_options = {
    # ===== Color Scheme =====
    "light_css_variables": {
        "color-brand-primary": "#2962ff",  # Blue for links/accents
        "color-brand-content": "#2962ff",
        "font-stack": "system-ui, -apple-system, sans-serif",
        "font-stack--monospace": "Consolas, Monaco, 'Courier New', monospace",
    },
    "dark_css_variables": {
        "color-brand-primary": "#4fc3f7",  # Lighter blue for dark mode
        "color-brand-content": "#4fc3f7",
    },
    
    # ===== Sidebar =====
    "sidebar_hide_name": False,  # Show project name in sidebar
    
    # ===== Navigation =====
    "navigation_with_keys": True,  # Use arrow keys to navigate
    
    # ===== Footer =====
    "footer_icons": [
        {
            "name": "GitHub",
            "url": "https://github.com/cameronabrams/sandlercubics",
            "html": """
                <svg stroke="currentColor" fill="currentColor" stroke-width="0" viewBox="0 0 16 16">
                    <path fill-rule="evenodd" d="M8 0C3.58 0 0 3.58 0 8c0 3.54 2.29 6.53 5.47 7.59.4.07.55-.17.55-.38 0-.19-.01-.82-.01-1.49-2.01.37-2.53-.49-2.69-.94-.09-.23-.48-.94-.82-1.13-.28-.15-.68-.52-.01-.53.63-.01 1.08.58 1.23.82.72 1.21 1.87.87 2.33.66.07-.52.28-.87.51-1.07-1.78-.2-3.64-.89-3.64-3.95 0-.87.31-1.59.82-2.15-.08-.2-.36-1.02.08-2.12 0 0 .67-.21 2.2.82.64-.18 1.32-.27 2-.27.68 0 1.36.09 2 .27 1.53-1.04 2.2-.82 2.2-.82.44 1.1.16 1.92.08 2.12.51.56.82 1.27.82 2.15 0 3.07-1.87 3.75-3.65 3.95.29.25.54.73.54 1.48 0 1.07-.01 1.93-.01 2.2 0 .21.15.46.55.38A8.013 8.013 0 0 0 16 8c0-4.42-3.58-8-8-8z"></path>
                </svg>
            """,
            "class": "",
        },
        {
            "name": "LinkedIn",
            "url": "https://linkedin.com/in/cameron-abrams-b0143398",
            "html": """
                <svg role="img" width="24" height="24" viewBox="0 0 24 24" xmlns="http://www.w3.org/2000/svg" fill="currentColor">
                    <title>LinkedIn</title>
                    <path d="M20.447 20.452h-3.554v-5.403c0-1.288-.025-2.945-1.796-2.945-1.796 0-2.071 1.4-2.071 2.847v5.501h-3.554V9.001h3.414v1.561h.05c.475-.9 1.637-1.797 3.368-1.797 3.599 0 4.262 2.368 4.262 5.446v6.241zM5.337 7.433c-1.144 0-2.072-.93-2.072-2.072 0-1.142.928-2.07 2.072-2.07 1.142 0 2.07.928 2.07 2.07 0 1.144-.928 2.072-2.07 2.072zM6.814 20.452H3.859V9.001h2.955v11.451zM22.225 0H1.771C.792 0 0 .774 0 1.729v20.543C0 23.225.792 24 1.771 24h20.451C23.208 24 24 23.225 24 22.272V1.729C24 .774 23.208 0 22.225 0z"/>
                </svg>
            """,
            "class": "linkedin-icon",
        },
    ],
    
    # ===== Source Links =====
    "source_repository": "https://github.com/cameronabrams/sandlercubics",
    "source_branch": "main",
    "source_directory": "docs/",
    
    # ===== Top Announcement =====
    # "announcement": "<em>⚠️ This is educational software only - not for production use!</em>",
}

# Optional: Add logo
# html_logo = "_static/logo.png"  # If you have one

# Optional: Add favicon
# html_favicon = "_static/favicon.ico"

# Additional HTML options
html_title = f"{project} v{version}"
html_short_title = project
html_logo = None  # Add path to logo if you have one: '_static/logo.png'
html_favicon = None  # Add path to favicon if you have one

# If true, "Created using Sphinx" is shown in the HTML footer. Default is True.
html_show_sphinx = True

# If true, "(C) Copyright ..." is shown in the HTML footer. Default is True.
html_show_copyright = True

# -- Options for LaTeX output ------------------------------------------------

latex_elements = {
    'papersize': 'letterpaper',
    'pointsize': '10pt',
    'preamble': r'''
        \usepackage{amsmath}
        \usepackage{amssymb}
    ''',
}

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title, author, documentclass [howto, manual, or own class]).
latex_documents = [
    (master_doc, 'sandlercubics.tex', 'sandlercubics Documentation',
     author, 'manual'),
]

# -- Options for intersphinx extension ---------------------------------------
# https://www.sphinx-doc.org/en/master/usage/extensions/intersphinx.html#configuration

intersphinx_mapping = {
    'python': ('https://docs.python.org/3', None),
    'numpy': ('https://numpy.org/doc/stable/', None),
    'scipy': ('https://docs.scipy.org/doc/scipy/', None),
}

# -- Options for todo extension ----------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/extensions/todo.html#configuration

todo_include_todos = True
