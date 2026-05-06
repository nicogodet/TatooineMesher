Parameters reference
====================

All parameters of
:func:`~tatooinemesher.algorithms.mesh_and_interpolate_alg.mesh_and_interpolate`,
grouped by role.

Input files
-----------

.. list-table::
   :header-rows: 1
   :widths: 25 15 60

   * - Name
     - Type
     - Description
   * - ``infile_axis``
     - ``str``
     - Hydraulic axis file (``.shp``, ``.i2s``). Must contain a single open
       polyline used to project sections.
   * - ``infile_cross_sections``
     - ``str``
     - Cross-sections file (``.shp`` POLYLINEZ, ``.i3s``). The first per-vertex
       attribute is treated as Z (bathymetry).
   * - ``attr_cross_sections``
     - ``str | None``
     - Optional shapefile attribute to identify each cross-section. When
       absent, sections are numbered in file order.
   * - ``infile_constraint_lines``
     - ``str | None``
     - 2D constraint lines file (``.shp`` POLYLINE / POLYLINEZ / POLYLINEM,
       ``.i2s``). When ``None`` and no 3D file is provided, lines are
       inferred from the section endpoints.
   * - ``infile_constraint_3D_lines``
     - ``str | None``
     - 3D Z-aware constraint lines (``.shp`` POLYLINEZ only). The Z profile
       along the line is propagated as a longitudinal correction to the bed
       between sections. See :ref:`z-correction`.

Discretization
--------------

.. list-table::
   :header-rows: 1
   :widths: 25 15 60

   * - Name
     - Type
     - Description
   * - ``long_step``
     - ``float`` (m)
     - Longitudinal step between intermediate cross-sections. The number of
       intermediate sections in each zone is :math:`\lceil X_l / \text{long\_step} \rceil`.
   * - ``constant_long_disc``
     - ``bool``
     - When ``True``, every zone gets the same number of intermediate
       sections (independent of its length) — useful when downstream tools
       require a constant longitudinal density.
   * - ``lat_step``
     - ``float | None`` (m)
     - Lateral step inside each bed. Mutually exclusive with ``nb_pts_lat``.
   * - ``nb_pts_lat``
     - ``int | None``
     - Override of the lateral count: every bed gets exactly this many
       transversal points. Only compatible with **two** constraint lines.
   * - ``dist_max``
     - ``float`` (m)
     - Snap distance for ``compute_dist_proj_axe`` and ``find_and_add_limits``.
       Default ``0.01``.

Interpolation methods
---------------------

.. list-table::
   :header-rows: 1
   :widths: 25 20 55

   * - Name
     - Allowed values
     - Description
   * - ``interp_constraint_lines``
     - ``LINEAR``,
       ``CARDINAL``,
       ``FINITE_DIFF``
     - Coordinate (and Z, when ``has_z=True``) interpolation along
       constraint lines.
       ``LINEAR`` uses :meth:`shapely.geometry.LineString.interpolate`;
       ``CARDINAL`` and ``FINITE_DIFF`` rely on a Cubic Hermite Spline (see
       :ref:`chs`).
   * - ``interp_values``
     - ``LINEAR``,
       ``B-SPLINE``,
       ``AKIMA``,
       ``PCHIP``,
       ``CUBIC_SPLINE``,
       ``BILINEAR``,
       ``BICUBIC``,
       ``BIVARIATE_SPLINE``
     - Bathymetry and per-node variable interpolation between cross-sections.
       The first five are 1D longitudinal blends; the last three (``BI``-prefix)
       are bivariate :math:`(u, v)` interpolators. See :ref:`values-interp`.
   * - ``project_straight_line``
     - ``bool``
     - When ``True``, each cross-section is projected onto the straight line
       joining its endpoints before interpolation (smoother bathy on curved
       reaches).

3D constraint line correction
-----------------------------

These two parameters only matter when ``infile_constraint_3D_lines`` is set.
See :ref:`z-correction` for the full math.

.. list-table::
   :header-rows: 1
   :widths: 25 15 60

   * - Name
     - Type
     - Description
   * - ``z_line_strength``
     - ``float`` ∈ [0, 1]
     - Pull factor toward the 3D line at the mid-distance between two
       sections. ``0`` keeps cross-section Z exactly between sections (only
       the line *shape deviation* is applied);
       ``1`` makes the boundary node match the line Z at :math:`x_l = 0.5`.
       Section Z is preserved at :math:`x_l \in \{0, 1\}` regardless. See
       :ref:`the pull bell curve <pull-bell>`.
   * - ``z_line_gap_scale``
     - ``float`` ≥ 0 (m\ :sup:`-1`)
     - Lateral falloff scale. With ``0`` (default), the correction is spread
       linearly through the bed:
       :math:`(1-x_t) \Delta Z_{L_1} + x_t \Delta Z_{L_2}`.
       Larger values raise the exponents, localizing the correction to the
       bed boundary where the line Z diverges sharply from the bathymetry.
       See :ref:`lateral-spread`.

Outputs
-------

.. list-table::
   :header-rows: 1
   :widths: 25 15 60

   * - Name
     - Type
     - Description
   * - ``outfile_nodes``
     - ``str | None``
     - Optional path for the node set (``.shp``, ``.xyz``). Only the
       triangulation vertices are exported, no triangles.
   * - ``outfile_mesh``
     - ``str | None``
     - Optional path for the full mesh (``.slf`` Serafin, ``.t3s`` BlueKenue,
       ``.xml`` LandXML). When set, ``build_mesh()`` is called.
   * - ``lang``
     - ``"en"`` | ``"fr"``
     - Variable name language for the Serafin (``.slf``) header.
   * - ``verbose``
     - ``bool``
     - Sets the logger level to ``DEBUG`` when ``True``.

.. _parameter-guide:

Parameter choice guide
----------------------

The defaults are reasonable starting points but the right values depend
heavily on the river geometry and on the modelling goal. The advice below
is distilled from the TUC2019 benchmark on the Vaugris reach (10 km of the
Rhône, 27 cross-sections, ≈200 m wide, [TUC2019]_ §IV) and from common
practice on CNR projects.

Discretization
^^^^^^^^^^^^^^

**Lateral discretization** — prefer ``lat_step`` over ``nb_pts_lat`` as
soon as the river width varies along the reach. Forcing a constant
lateral count produces elements with poor minimum-angle statistics
(Mesh #4 in [TUC2019]_ Fig. 11) — the worst quality among all tested
configurations. Reserve ``nb_pts_lat`` for fully prismatic test cases or
when a downstream tool requires a constant transversal count.

**Element elongation** — pick ``long_step`` ≳ 1.5 × ``lat_step`` to
elongate elements along the flow direction. Mesh #1 vs #2 in [TUC2019]_
shows that going from a 1:1 ratio to 1:1.5 cuts the node count by ~25 %
with only a marginal hit on the minimum angle distribution. The
CFL-bound time step depends on the smallest edge, so elongation is a
free win when the flow is mostly longitudinal.

**Longitudinal discretization mode** — enable ``constant_long_disc=True``
when downstream visualization benefits from a structured-looking mesh
(e.g. side-by-side comparison with a 1D longitudinal profile). Otherwise
leave it off — the unstructured per-bed discretization adapts to local
widening / narrowing of the channel.

Constraint line interpolation (``interp_constraint_lines``)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* No constraint line file → prefer ``CARDINAL`` (Cubic Hermite Spline,
  see :ref:`chs`). When TatooineMesher synthesizes the two banks from
  the section endpoints, ``CARDINAL`` keeps the river width more
  conservative than ``LINEAR`` on curved reaches ([TUC2019]_ §II.C.1,
  Fig. 3).
* User-provided constraint lines with sparse vertices → ``CARDINAL``
  smooths the polyline between vertices without overshoot.
* User-provided constraint lines with dense vertices → ``LINEAR`` is
  cheap and visually indistinguishable.
* ``FINITE_DIFF`` behaves like ``CARDINAL`` with slightly tighter
  tangents at sharp bends; it is rarely needed in practice.

Bathymetry interpolation (``interp_values``)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The TUC2019 benchmark on Vaugris reports ([TUC2019]_ Fig. 13):

.. list-table::
   :header-rows: 1
   :widths: 25 15 15 15 15

   * - Mode
     - MSD [m]
     - MAD [m]
     - RMSD [m]
     - Note
   * - ``LINEAR`` (1D)
     - ≈ +0.025
     - 0.50
     - 0.69
     - Robust baseline
   * - ``AKIMA`` (1D)
     - ≈ −0.005
     - **0.47**
     - **0.65**
     - Best across all metrics
   * - ``B-SPLINE`` (1D)
     - ≈ −0.015
     - 0.49
     - 0.69
     - Slightly oscillates on sparse profiles
   * - ``PCHIP`` (1D)
     - ≈ +0.005
     - 0.48
     - 0.66
     - Good middle ground, monotone
   * - ``BILINEAR``
     - ≈ +0.015
     - 0.49
     - 0.69
     - Comparable to LINEAR
   * - ``BICUBIC``
     - ≈ −0.025
     - 0.49
     - 0.69
     - Comparable to LINEAR

Concrete advice:

* Densely sampled cross-sections (≥ 50 points per profile) → ``LINEAR``
  is enough. The lateral interpolator has marginal influence in that
  regime.
* Sparse cross-sections (≤ 25 points per profile, as in Vaugris with
  20) → ``AKIMA``. The Akima spline removes the staircase artefacts of
  ``LINEAR`` without introducing the wiggles of ``CUBIC_SPLINE``.
* Avoid ``CUBIC_SPLINE`` on bathymetry: classical cubic splines
  oscillate around discontinuities (banks, levees) — see [TUC2019]_
  Fig. 9 and the comparison plot in
  :ref:`the lateral interpolation comparison <lateral-interp-comparison>`.
* 2D modes (``BILINEAR``, ``BICUBIC``, ``BIVARIATE_SPLINE``) require
  exactly two constraint lines along the entire reach: their global
  :math:`(x_l, x_t)` parameterization is undefined when the bed topology
  changes. Use them only on prismatic reaches with continuous banks.
  ``BILINEAR`` is the safest 2D choice; ``BIVARIATE_SPLINE`` smoothes
  more aggressively and can mask narrow features.
* ``PCHIP`` is a good compromise when monotonicity matters (e.g.
  preserving a thalweg low point without overshoot).

Cross-section projection
^^^^^^^^^^^^^^^^^^^^^^^^

* Enable ``project_straight_line=True`` when input cross-sections are
  noticeably squiggly (sloppy digitizing or natural sinuosity at the
  scale of the section). Elements adjacent to the section become more
  regular ([TUC2019]_ Fig. 4) at the cost of shifting the original
  bathymetry points to the chord between the two banks.
* Disable it when cross-sections are clean lines and you want to
  preserve the original bathymetry XY positions exactly.

3D constraint line correction
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

These two parameters are extensions specific to TatooineMesher (post
TUC2019) and steer how the Z profile of a 3D constraint line bleeds
into the bed bathymetry between sections.

* Start with ``z_line_strength = 0``. With strength 0 only the line
  *shape* (its non-linear deviation from a chord) is applied;
  cross-section Z is preserved at section locations and the long-term
  trend is unchanged. This is usually enough when the 3D line tracks a
  thalweg whose absolute Z is already encoded in the cross-sections.
* Raise ``z_line_strength`` to 0.25 – 0.5 when the 3D line carries Z
  information not covered by the cross-sections (e.g. a bank crest
  measured by LiDAR but digitized only roughly in the cross-sections).
  Avoid values close to 1.0 unless you are confident in the absolute Z
  of the line — it will overrule the bathymetry at the mid-distance
  between sections.
* Keep ``z_line_gap_scale = 0`` as long as the line Z is consistent
  with the bed: the correction will spread linearly through the bed,
  like a standard bilinear blend.
* Increase ``z_line_gap_scale`` (start at 0.5–1.0 m\ :sup:`-1`) when
  the 3D line diverges sharply from the bed — for instance a narrow
  thalweg in a wide cross-section. Larger values keep the correction
  localized close to the line. See :ref:`lateral-spread` for a plot of
  the lateral weight :math:`(1-x_t)^p` for several gap magnitudes.

When in doubt, sweep ``z_line_strength`` ∈ {0, 0.25, 0.5} and
``z_line_gap_scale`` ∈ {0, 1} and visually inspect the result on a
representative reach before committing to a value.
