import importlib.metadata
import os
import sys
from datetime import datetime

sys.path.insert(0, os.path.abspath("../.."))

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "SpatialLeiden"
copyright = f"""
{datetime.now():%Y}, Niklas Müller-Bötticher, Shashwat Sahay, Naveed Ishaque, Roland Eils,
Berlin Institute of Health @ Charité"""
author = "Niklas Müller-Bötticher & Shashwat Sahay"
version = importlib.metadata.version("spatialleiden")
release = version

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration


extensions = [
    "sphinx_copybutton",
    "sphinx.ext.autodoc",
    "sphinx.ext.intersphinx",
    "sphinx.ext.napoleon",
    "sphinx.ext.autosummary",
    "sphinx.ext.mathjax",
    "myst_nb",
]


autodoc_typehints = "none"
autodoc_typehints_format = "short"

python_use_unqualified_type_names = True  # still experimental

autosummary_generate = True
autosummary_imported_members = True

nitpicky = True
nitpick_ignore = [("py:class", "optional")]

# MyST-NB config
nb_execution_timeout = 3 * 60

exclude_patterns: list[str] = []

intersphinx_mapping = dict(
    anndata=("https://anndata.readthedocs.io/en/stable/", None),
    matplotlib=("https://matplotlib.org/stable/", None),
    mudata=("https://mudata.readthedocs.io/en/stable/", None),
    numpy=("https://numpy.org/doc/stable/", None),
    python=("https://docs.python.org/3", None),
    scipy=("https://docs.scipy.org/doc/scipy/", None),
    squidpy=("https://squidpy.readthedocs.io/en/stable/", None),
)

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "sphinx_rtd_theme"
html_static_path: list[str] = []
