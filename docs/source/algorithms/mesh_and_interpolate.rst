``mesh_and_interpolate``
========================

.. currentmodule:: tatooinemesher.algorithms.mesh_and_interpolate_alg

The central algorithm of TatooineMesher. From a hydraulic axis, a sequence
of 1D cross-sections, and an optional set of (2D and/or 3D) constraint
lines, it builds a 2D triangular mesh aligned with the flow direction and
interpolates bathymetry — and any other per-section variable — onto every
node.

End-to-end pipeline
-------------------

.. mermaid::

   flowchart TD
       subgraph Inputs
           direction LR
           AX[axis.shp]
           CS[cross_sections.shp]
           CL2D[constraint_lines.shp<br/>2D]
           CL3D[constraint_3D_lines.shp<br/>POLYLINEZ]
       end

       subgraph Step1[Read & project]
           R1[get_hydraulic_axis]
           R2[CrossSectionSequence.from_file]
           R3[compute_dist_proj_axe]
           R4[check_intersections]
           R5[sort_by_dist]
           R1 --> R3 --> R4 --> R5
           R2 --> R3
       end

       subgraph Step2[Constraint lines]
           CL[ConstraintLine.get_lines_from_file<br/>or .get_lines_and_set_limits_from_sections]
           FL[find_and_add_limits]
           CL --> FL
       end

       subgraph Step3[Mesh interpolation]
           BI[MeshConstructor.build_interp]
           BIP[build_initial_profiles]
           SAM[coord_sampling_along_line]
           ZC[_compute_line_z_correction]
           ICL[interp_coord_linear<br/>per bed]
           AP[add_points]
           BI --> BIP --> SAM --> ZC --> ICL --> AP
       end

       subgraph Step4[Triangulate & values]
           BM[build_mesh<br/>via triangle]
           IV[interp_values_from_geom<br/>1D or 2D mode]
           BM --> IV
       end

       subgraph Outputs
           direction LR
           ON[outfile_nodes<br/>.shp .xyz]
           OM[outfile_mesh<br/>.slf .t3s .xml]
       end

       Inputs --> Step1 --> Step2 --> Step3 --> Step4 --> Outputs

The four steps below match the canonical mesh-generation skeleton from
[TUC2019]_, with steps 3 and 4 fused into the ``MeshConstructor`` machinery.

Step 1 — Read inputs and project
--------------------------------

* :func:`tatooinemesher.utils.get_hydraulic_axis` loads ``infile_axis``
  (single open polyline) and returns a :class:`shapely.geometry.LineString`
  augmented with a curvilinear abscissa :math:`X_t` (see
  :ref:`curvilinear-abscissa`).
* :meth:`tatooinemesher.section.CrossSectionSequence.from_file` parses
  ``infile_cross_sections`` (POLYLINEZ shapefile or ``.i3s``). Each
  cross-section is a :class:`~tatooinemesher.section.CrossSection` with its
  own structured points array; the first per-vertex attribute is treated as
  Z (bathymetry).
* :meth:`~tatooinemesher.section.CrossSectionSequence.compute_dist_proj_axe`
  projects each cross-section onto the axis and stores its
  ``dist_proj_axe`` (:math:`X_l`). Sections farther than ``dist_max`` from
  the axis raise.
* :meth:`~tatooinemesher.section.CrossSectionSequence.check_intersections`
  validates that successive sections do not cross each other (catches a
  common shapefile authoring mistake).
* :meth:`~tatooinemesher.section.CrossSectionSequence.sort_by_dist` orders
  the sequence by increasing :math:`X_l`. Sections that do not cross the
  axis are dropped here — see [TUC2019]_ §II.B step 1.

Step 2 — Build constraint lines
-------------------------------

Three modes, mutually compatible inside the same call:

1. **No file given** —
   :meth:`~tatooinemesher.constraint_line.ConstraintLine.get_lines_and_set_limits_from_sections`
   builds two lines from the section endpoints (left bank, right bank).
2. **2D shapefile** (``--infile_constraint_lines``) — POLYLINE / POLYLINEZ /
   POLYLINEM, parsed via
   :meth:`~tatooinemesher.constraint_line.ConstraintLine.get_lines_from_file`.
3. **3D shapefile** (``--infile_constraint_3D_lines``) — POLYLINEZ only,
   parsed with ``has_z=True``. The Z profile along each line is interpolated
   and applied as a longitudinal correction to the bathymetry between
   sections (see :ref:`z-correction`).

When both 2D and 3D files are given, the lists are concatenated and the
IDs are renumbered contiguously. The merged set must contain exactly two
lines when ``--nb_pts_lat`` is set, or when ``--interp_values`` is a 2D
mode (``BILINEAR``, ``BICUBIC``, ``BIVARIATE_SPLINE``).

:meth:`~tatooinemesher.section.CrossSectionSequence.find_and_add_limits`
then attaches each line as a limit on every section it crosses (snap
distance ``dist_max``).

Step 3 — Mesh interpolation
---------------------------

:class:`~tatooinemesher.mesh_constructor.MeshConstructor` walks the
sequence in
:meth:`~tatooinemesher.mesh_constructor.MeshConstructor.build_interp`:

* :meth:`~tatooinemesher.mesh_constructor.MeshConstructor.build_initial_profiles`
  emits the nodes of the first and last cross-sections without
  inter-section interpolation.
* For each *zone* (consecutive section pair at :math:`i, i+1`):

  * The longitudinal step ``long_step`` (or a constant per-zone count when
    ``constant_long_disc=True``) yields the intermediate :math:`x_l`
    values between the two sections.
  * For each *bed* (pair of constraint lines bounding a sub-channel),
    coordinates and Z correction are sampled along the two bordering lines
    via
    :meth:`~tatooinemesher.constraint_line.ConstraintLine.coord_sampling_along_line`.
  * When at least one of the two lines is 3D (``has_z=True``),
    :meth:`~tatooinemesher.mesh_constructor.MeshConstructor._compute_line_z_correction`
    yields the Z correction at the boundary nodes (see :ref:`z-correction`).
  * :meth:`~tatooinemesher.section.Bed.interp_coord_linear` interpolates
    transversal node coordinates (X, Y) inside the bed using
    :math:`x_l, x_t \in [0, 1]` (see :ref:`bed-linear-interp`).
  * :meth:`~tatooinemesher.mesh_constructor.MeshConstructor.add_points`
    appends the resulting nodes to the structured ``points`` array,
    carrying their ``z_correction`` along.

Step 4 — Triangulate and interpolate values
-------------------------------------------

* :meth:`~tatooinemesher.mesh_constructor.MeshConstructor.build_mesh`
  delegates to the :mod:`triangle` library on the nodes and constrained
  segments, producing a planar straight-line graph with cross-sections
  and constraint lines forced as hard segments [TUC2019]_ §II.B step 4.
* Bathymetry and other variables are sampled in
  :meth:`~tatooinemesher.mesh_constructor.MeshConstructor.interp_values_from_geom`,
  which calls either
  :meth:`~tatooinemesher.mesh_constructor.MeshConstructor.interp_1d_values_from_profiles`
  (1D modes) or
  :meth:`~tatooinemesher.mesh_constructor.MeshConstructor.interp_2d_values_from_profiles`
  (2D modes, ``BI``-prefix). The per-node ``z_correction`` is added to the
  first variable (Z / elevation) at the very end.
* Outputs:

  * ``outfile_nodes`` — nodes only (``.shp`` or ``.xyz``);
  * ``outfile_mesh`` — full mesh (``.slf`` Serafin via PyTelTools, ``.t3s``
    BlueKenue, ``.xml`` LandXML).

For the mathematical content of each interpolation step, jump to
:doc:`../math`. For the parameter cheat-sheet, see :doc:`../parameters`.

.. seealso::

   - :doc:`../parameters` — full reference and a parameter choice guide.
   - :doc:`../math` — every interpolation formula with LaTeX and plots.
   - :doc:`../api` — autodoc signature of :func:`mesh_and_interpolate`.
