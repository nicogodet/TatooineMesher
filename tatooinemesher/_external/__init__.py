"""
Vendored third-party code.

This package contains code copied from external projects to keep TatooineMesher's
runtime dependency footprint small. Each subpackage carries the license and
attribution of its upstream source.

Currently vendored:
- pyteltools.geom and pyteltools.slf — copied from
  https://github.com/CNR-Engineering/PyTelTools (originally under GPL-3.0).
  Only the modules actually used by TatooineMesher are vendored, which lets
  us avoid pulling in PyQt5, descartes, pyproj and Rtree as transitive
  dependencies. See `pyteltools/__init__.py` for the upstream commit/version.
"""
