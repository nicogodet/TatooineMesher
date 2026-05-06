``mesh_crue10_run``
===================

.. currentmodule:: tatooinemesher.algorithms.mesh_crue10_run_alg

Builds a 2D mesh and interpolates the time-varying results of a Crue10
1D simulation onto it, for visualization in TELEMAC's Serafin format.
Same principle as :doc:`mesh_and_interpolate`, with a Crue10 model
(branches + cross-sections + casiers) replacing the standalone shapefile
inputs and with hydraulic results (steady states or transient frames)
written as time-dependent variables on the mesh.

Requires the ``crue10`` extra (see :doc:`../installation`).

.. note::

   Detailed pipeline, parameter list, and worked examples are not yet
   documented. Track progress in
   `GitHub issues <https://github.com/CNR-Engineering/TatooineMesher/issues>`_.

   In the meantime:

   * Source: :file:`tatooinemesher/algorithms/mesh_crue10_run_alg.py`
   * CLI help: ``python cli/mesh_crue10_run.py --help``
   * Floodplain (casier) meshing is supported; provide a raster (``--infile_dem``)
     for floodplain bottom elevation. Requires the ``gdal`` extra.
