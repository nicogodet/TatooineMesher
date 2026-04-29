``mesh_mascaret_run``
=====================

.. currentmodule:: tatooinemesher.algorithms.mesh_mascaret_run_alg

Builds a 2D mesh and interpolates the time-varying results of a MASCARET
1D simulation onto it, for visualization in TELEMAC's Serafin format.
Mirrors :doc:`mesh_crue10_run` for MASCARET geometry / Opthyca or Listing
result files.

Requires the ``crue10`` extra (its ``Mascaret`` parser is bundled with
``Crue10_tools``); see :doc:`../installation`.

.. note::

   Detailed pipeline, parameter list, and worked examples are not yet
   documented. Track progress in
   `GitHub issues <https://github.com/CNR-Engineering/TatooineMesher/issues>`_.

   In the meantime:

   * Source: :file:`tatooinemesher/algorithms/mesh_mascaret_run_alg.py`
   * CLI help: ``python cli/mesh_mascaret_run.py --help``
   * Variables are interpolated identically to
     :doc:`mesh_and_interpolate` — the lateral / longitudinal trade-offs
     in :doc:`../parameters` apply.
