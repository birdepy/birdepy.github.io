# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
# once BirDePy is on PyPI the below can be replaced with 'import birdepy'

import os
import sys
sys.path.insert(0, os.path.abspath('.'))


# -- Project information -----------------------------------------------------

project = 'BirDePy'
copyright = '2021, Sophie Hautphenne and Brendan Patch'
author = 'Sophie Hautphenne and Brendan Patch'

# The full version, including alpha/beta/rc tags
release = '0.0.10'


html_favicon = 'Birdepy_favicon_3.png'

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = ['sphinx_rtd_theme',
'sphinx.ext.autodoc',
'sphinx.ext.napoleon',
'sphinx.ext.autosectionlabel',
'sphinx.ext.githubpages',
'sphinx_copybutton'
]

# The master toctree document.
master_doc = 'index'

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'
html_logo = 'Birdepy_logo_1_resized.png'
html_theme_options = {
'style_nav_header_background': 'white',
'collapse_navigation': False,
'navigation_depth': 4,
'logo_only': True,
'titles_only': False,
}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

html_css_files = [
    'css/custom.css',
]

