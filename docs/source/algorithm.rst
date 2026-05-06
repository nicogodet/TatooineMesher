Algorithms overview
===================

TatooineMesher exposes four algorithms, each with its own CLI script under
:file:`cli/` and a corresponding Python entry point under
:mod:`tatooinemesher.algorithms`. They share the same underlying primitives
(cross-section reading, constraint-line interpolation, lateral / longitudinal
sampling, anisotropic bathymetry interpolation).

Algorithm catalog
-----------------

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Algorithm
     - Purpose
   * - :doc:`algorithms/mesh_and_interpolate`
     - Build a 2D triangular mesh between cross-sections and interpolate
       bathymetry (and any other per-section variables) onto its nodes.
       Handles 2D and 3D constraint lines.
   * - :doc:`algorithms/densify_cross_sections`
     - Refine 1D cross-sections for further 1D simulations (no mesh).
   * - :doc:`algorithms/mesh_crue10_run`
     - Mesh and interpolate the geometry **and** time-varying results of a
       Crue10 1D simulation for 2D visualization.
   * - :doc:`algorithms/mesh_mascaret_run`
     - Same as above, for MASCARET 1D simulations.

This documentation currently covers :doc:`algorithms/mesh_and_interpolate`
in depth. The other three are presented as stubs and will be filled in as
the surrounding work lands.

.. toctree::
   :hidden:
   :maxdepth: 2

   algorithms/mesh_and_interpolate
   algorithms/densify_cross_sections
   algorithms/mesh_crue10_run
   algorithms/mesh_mascaret_run

Common pipeline
---------------

All four algorithms share a four-step skeleton inherited from the original
TatooineMesher specification [TUC2019]_:

.. mermaid::

   flowchart LR
       A[Order cross-sections] --> B[Intersect with constraint lines]
       B --> C[Generate nodes per submesh]
       C --> D[Constrained Delaunay triangulation]

* **Step 1 — Order**: read cross-sections, project them on the hydraulic
  axis, drop those that do not intersect it, sort by curvilinear abscissa.
* **Step 2 — Intersect**: identify the limits where each constraint line
  crosses each cross-section. Each pair (line :math:`L_j`, line :math:`L_{j+1}`)
  defines a *bed*.
* **Step 3 — Generate nodes**: for every (zone, bed) pair, sample
  intermediate cross-sections at :math:`x_l \in (0, 1)` and lateral nodes at
  :math:`x_t \in [0, 1]`. The flow-oriented coordinates :math:`(x_l, x_t)`
  make the bathymetry interpolation anisotropic — see :doc:`math`.
* **Step 4 — Triangulate**: a constrained Delaunay triangulation is run
  (via :mod:`triangle` [Shewchuk1996]_) with the cross-sections and
  constraint lines forced as hard segments. The resulting "planar
  straight-line graph" satisfies the standard mesh-quality criteria
  reviewed in [Merkel2013]_; alternative meshers like Gmsh
  [Geuzaine2009]_ may be used downstream when adapting node density to
  hydraulic structures (bridge piers, groynes).

Where each algorithm differs is in the input data type (raw shapefile vs.
1D model), in the output produced (densified profiles vs. mesh + values),
and in whether time-dependent variables are advected onto the mesh.

For the math behind every interpolation step, follow the links in
:doc:`algorithms/mesh_and_interpolate` or jump straight to :doc:`math`.
