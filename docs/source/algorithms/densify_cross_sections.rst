``densify_cross_sections``
==========================

.. currentmodule:: tatooinemesher.algorithms.densify_cross_sections_alg

Refines a 1D set of cross-sections without producing a 2D mesh, by
inserting intermediate sections between the user-provided ones. Used as a
pre-treatment for 1D simulations (MASCARET, Crue10) where the cross-section
spacing is too coarse for the desired resolution.

The algorithm reuses the same skeleton as :doc:`mesh_and_interpolate`
(steps 1–3) but stops before the triangulation: only the intermediate
sections, with their interpolated bathymetry, are exported.

.. note::

   Detailed pipeline, parameter list, and worked examples are not yet
   documented. Track progress in
   `GitHub issues <https://github.com/CNR-Engineering/TatooineMesher/issues>`_.

   In the meantime:

   * Source: :file:`tatooinemesher/algorithms/densify_cross_sections_alg.py`
   * CLI help: ``python cli/densify_cross_sections.py --help``
   * The interpolation modes ``LINEAR``, ``B-SPLINE``, ``AKIMA``, ``PCHIP``,
     ``CUBIC_SPLINE`` (see :ref:`values-interp`) all apply here too.
