# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'zonoopt'
copyright = '2025, Joshua Robbins'
author = 'Joshua Robbins'
release = '2.0.1'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon', # Recommended for NumPy/Google style docstrings
    'sphinx.ext.autosummary', # Optional: for generating summary tables
    'sphinx.ext.mathjax',
    'myst_parser'
]

myst_enable_extensions = [
    "dollarmath",
    "amsmath",
]

mathjax3_config = {
    'tex': {
        # Force the 'ams' package to load
        'packages': {'[+]': ['ams']},
        # Ensure inline math works
        'inlineMath': [['$', '$'], ['\\(', '\\)']]
    },
    'loader': {
        # Load the ams component
        'load': ['[tex]/ams']
    }
}

templates_path = ['_templates']
exclude_patterns = []

autoclass_content = "both"



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['../../../images']

def skip_special_members(app, what, name, obj, skip, options):
    if name.startswith('__') and name.endswith('__') and name not in ['__add__', '__sub__', '__mul__', '__truediv__', '__getitem__', '__setitem__']:
        return True
    return skip

def setup(app):
    app.connect('autodoc-skip-member', skip_special_members)