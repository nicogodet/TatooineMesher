TatooineMesher
==============

Mesh and interpolate rivers from 1D cross-section profiles and constraint lines.

This documentation focuses on the central algorithm, :func:`mesh_and_interpolate`,
which builds a 2D triangular mesh from a hydraulic axis, a sequence of
cross-sections, and an optional set of (2D or 3D) constraint lines, and
interpolates bathymetry and other variables onto the resulting nodes.

.. toctree::
   :maxdepth: 2
   :caption: User guide

   installation
   quickstart
   algorithm
   parameters

.. toctree::
   :maxdepth: 2
   :caption: Reference

   math
   api
   references

Indices and tables
------------------

* :ref:`genindex`
* :ref:`search`
