"""Sphinx configuration for TatooineMesher documentation."""

import os
import sys

import matplotlib

matplotlib.use("Agg")

sys.path.insert(0, os.path.abspath("../.."))

from tatooinemesher import VERSION  # noqa: E402

project = "TatooineMesher"
author = "Luc Duron"
copyright = "2018-2026, CNR"
release = VERSION
version = ".".join(VERSION.split(".")[:2])

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.mathjax",
    "sphinx.ext.intersphinx",
    "sphinx.ext.viewcode",
    "matplotlib.sphinxext.plot_directive",
    "myst_parser",
    "sphinx_copybutton",
    "sphinxcontrib.mermaid",
]

mermaid_version = "11.4.1"
mermaid_init_js = "mermaid.initialize({startOnLoad:true, theme:'default', flowchart:{useMaxWidth:true}});"

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "scipy": ("https://docs.scipy.org/doc/scipy/", None),
    "shapely": ("https://shapely.readthedocs.io/en/stable/", None),
    "matplotlib": ("https://matplotlib.org/stable/", None),
}

html_theme = "furo"
html_static_path = ["_static"]
html_title = f"TatooineMesher {version}"

templates_path = ["_templates"]
exclude_patterns = []

plot_include_source = True
plot_html_show_source_link = False
plot_html_show_formats = False
plot_formats = [("png", 100)]

autodoc_default_options = {
    "members": True,
    "show-inheritance": True,
}
autodoc_typehints = "description"
autodoc_mock_imports = ["pyteltools", "crue10", "osgeo", "shapefile", "triangle", "jinja2"]

napoleon_google_docstring = True
napoleon_numpy_docstring = True

mathjax3_config = {"tex": {"macros": {}}}

suppress_warnings = ["autodoc.import_object"]
