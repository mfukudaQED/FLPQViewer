# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys
sys.path.insert(0, os.path.abspath(".."))

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'FLPQViewer'
copyright = '2025, Masahiro FUKUDA'
author = 'Masahiro FUKUDA'
release = '0.1.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.autodoc",  # docstring $B$+$i%I%-%e%a%s%H$r<+F0@8@.(B
    "sphinx.ext.napoleon",  # Google-style docstring $B$r%5%]!<%H(B
    "sphinx.ext.viewcode",  # $B%=!<%9%3!<%I$X$N%j%s%/$rDI2C(B
]

templates_path = ['_templates']
exclude_patterns = []



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'alabaster'
html_static_path = ['_static']
