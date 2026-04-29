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

The function executes as a linear chain of six stages:

.. mermaid::

   flowchart LR
       I["Inputs<br/>axis · sections<br/>· constraint lines"] --> S1["Step 1<br/>Read &amp; project"]
       S1 --> S2["Step 2<br/>Constraint lines<br/>+ limits"]
       S2 --> S3["Step 3<br/>build_interp"]
       S3 --> S4["Step 4<br/>Triangulate<br/>+ values"]
       S4 --> O["Outputs<br/>nodes · mesh"]

The four numbered steps match the canonical mesh-generation skeleton
from [TUC2019]_, with steps 3 and 4 fused into the ``MeshConstructor``
machinery. Step 3 itself is a triple loop over zones, beds, and
intermediate intermediate cross-sections — see
:ref:`mai-build-interp-loop` below.

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

.. _mai-build-interp-loop:

Step 3 — Mesh interpolation
---------------------------

:meth:`~tatooinemesher.mesh_constructor.MeshConstructor.build_interp`
is a triple loop:

.. mermaid::

   flowchart TD
       BIP[build_initial_profiles<br/>nodes of first &amp; last sections]
       BIP --> Z["For each zone (section i, i+1)"]
       Z --> B["For each bed (L_j, L_j+1)"]
       B --> SAM[coord_sampling_along_line<br/>on L_j and L_j+1]
       SAM --> ZC[_compute_line_z_correction<br/>if has_z]
       ZC --> XL["For each intermediate x_l"]
       XL --> ICL[Bed.interp_coord_linear<br/>lateral nodes]
       ICL --> AP[add_points<br/>append to MeshConstructor]
       AP --> XL
       XL -. next bed .-> B
       B -. next zone .-> Z

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
