Mathematical reference
======================

This chapter spells out every interpolation formula used by
:func:`~tatooinemesher.algorithms.mesh_and_interpolate_alg.mesh_and_interpolate`,
with cross-references back to the implementation.

Notations
---------

.. list-table::
   :header-rows: 1
   :widths: 15 85

   * - Symbol
     - Meaning
   * - :math:`X_t`
     - Curvilinear abscissa along a polyline (m), measured from its first
       vertex.
   * - :math:`x_t \in [0, 1]`
     - Dimensionless lateral position inside a bed,
       :math:`x_t = (X_t - X_{t,0}) / (X_{t,\max} - X_{t,0})`.
   * - :math:`X_l`
     - Longitudinal curvilinear abscissa along the hydraulic axis (m).
   * - :math:`x_l \in [0, 1]`
     - Dimensionless longitudinal position inside a *zone* (between two
       consecutive cross-sections); ``0`` is the upstream section, ``1`` the
       downstream section.
   * - :math:`Z`
     - Bathymetry / elevation (m).
   * - *zone*
     - Region between two consecutive cross-sections.
   * - *bed*
     - Sub-channel of a cross-section, bounded by two constraint lines.

.. _curvilinear-abscissa:

Curvilinear abscissa along a polyline
-------------------------------------

Implemented in
:meth:`tatooinemesher.constraint_line.ConstraintLine.__init__` and
:class:`tatooinemesher.coord.Coord`.

For an open polyline of vertices :math:`(x_i, y_i)_{i=0..N}`,

.. math::

   X_{t,i} = \sum_{j=1}^{i} \sqrt{(x_j - x_{j-1})^2 + (y_j - y_{j-1})^2},
   \quad X_{t,0} = 0.

The dimensionless form normalizes against the bounds of the segment of
interest:

.. math::

   x_{t,i} = \frac{X_{t,i} - X_{t,0}}{X_{t,\max} - X_{t,0}} \in [0, 1].

Linear interpolation on constraint lines
----------------------------------------

Implemented in
:meth:`tatooinemesher.constraint_line.ConstraintLine.build_interp_linear`.

X and Y are interpolated by walking the polyline with
:meth:`shapely.geometry.LineString.interpolate`. When ``has_z=True``, Z is
interpolated linearly between the polyline's vertex Z values:

.. math::

   Z(X_t) = \mathrm{interp}\bigl(X_t,\ \{X_{t,i}\}_{i=0..N},\ \{Z_i\}_{i=0..N}\bigr)

i.e. piecewise-linear with breakpoints at the polyline vertices
(:func:`numpy.interp`).

.. _chs:

Cubic Hermite Spline
--------------------

Implemented in :class:`tatooinemesher.interp.cubic_hermite_spline.CubicHermiteSpline`.
Selected via ``interp_constraint_lines = "CARDINAL"`` or ``"FINITE_DIFF"``.

Hermite basis functions on the unit interval :math:`t \in [0, 1]`:

.. math::

   \begin{aligned}
   h_{00}(t) &= 2t^3 - 3t^2 + 1, \\
   h_{10}(t) &= t^3 - 2t^2 + t, \\
   h_{01}(t) &= -2t^3 + 3t^2, \\
   h_{11}(t) &= t^3 - t^2.
   \end{aligned}

For a query :math:`x \in [x_i, x_{i+1}]`, with
:math:`t = (x - x_i)/\Delta x_i` and :math:`\Delta x_i = x_{i+1} - x_i`:

.. math::

   p(x) = h_{00}(t)\, p_i
        + h_{10}(t)\, \Delta x_i\, m_i
        + h_{01}(t)\, p_{i+1}
        + h_{11}(t)\, \Delta x_i\, m_{i+1}.

The interior tangents :math:`m_i` are set by the chosen ``tan_method``:

**Cardinal spline** (``CARDINAL``, default, with ``c = 0`` Catmull-Rom):

.. math::

   m_i = (1 - c)\,\frac{p_{i+1} - p_{i-1}}{x_{i+1} - x_{i-1}}.

**Finite-difference** (``FINITE_DIFF``):

.. math::

   m_i = \tfrac{1}{2}\!\left(
     \frac{p_{i+1} - p_i}{x_{i+1} - x_i}
     + \frac{p_i - p_{i-1}}{x_i - x_{i-1}}
   \right).

End tangents (``end_tan = GRAD``, default) use the slope of the first /
last segment scaled by ``m`` (default ``1``).

.. plot::
   :caption: Tangent methods on a 1D test polyline (left) and shape
             preservation on a 2D river-bank-like polyline (right). On the
             2D case, ``CARDINAL`` produces a smoother, less distorted
             curve — the river width is more conservative than with
             ``LINEAR`` ([TUC2019]_ §II.C.1).
   :include-source: false

   import numpy as np
   import matplotlib.pyplot as plt

   def hermite_eval(x_query, X, Y, M):
       i = np.clip(np.searchsorted(X, x_query) - 1, 0, len(X) - 2)
       dx = X[i + 1] - X[i]
       t = (x_query - X[i]) / dx
       h00 = 2 * t**3 - 3 * t**2 + 1
       h10 = t**3 - 2 * t**2 + t
       h01 = -2 * t**3 + 3 * t**2
       h11 = t**3 - t**2
       return h00 * Y[i] + h10 * dx * M[i] + h01 * Y[i + 1] + h11 * dx * M[i + 1]

   def cardinal_tangents(X, Y, c=0.0):
       M = np.zeros_like(Y, dtype=float)
       M[1:-1] = (1 - c) * (Y[2:] - Y[:-2]) / (X[2:] - X[:-2])
       M[0] = (Y[1] - Y[0]) / (X[1] - X[0])
       M[-1] = (Y[-1] - Y[-2]) / (X[-1] - X[-2])
       return M

   def finite_diff_tangents(X, Y):
       M = np.zeros_like(Y, dtype=float)
       d_fwd = (Y[2:] - Y[1:-1]) / (X[2:] - X[1:-1])
       d_bwd = (Y[1:-1] - Y[:-2]) / (X[1:-1] - X[:-2])
       M[1:-1] = 0.5 * (d_fwd + d_bwd)
       M[0] = (Y[1] - Y[0]) / (X[1] - X[0])
       M[-1] = (Y[-1] - Y[-2]) / (X[-1] - X[-2])
       return M

   fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))

   ax = axes[0]
   xs = np.array([0.0, 1.0, 2.0, 3.5, 5.0])
   ys = np.array([0.0, 2.0, 1.0, 3.0, 2.5])
   xq = np.linspace(xs[0], xs[-1], 400)
   ax.plot(xs, ys, "ko", label="vertices")
   ax.plot(xq, hermite_eval(xq, xs, ys, cardinal_tangents(xs, ys)),
           label="CARDINAL (c=0, Catmull-Rom)")
   ax.plot(xq, hermite_eval(xq, xs, ys, finite_diff_tangents(xs, ys)),
           "--", label="FINITE_DIFF")
   ax.set_title("Tangent methods on a 1D polyline")
   ax.set_xlabel("x"); ax.set_ylabel("y")
   ax.legend(fontsize=9); ax.grid(alpha=0.3)

   ax = axes[1]
   theta = np.linspace(0.4, 2.6, 7)
   bx = 4 * np.cos(theta)
   by = 4 * np.sin(theta) - np.linspace(0, 0.6, 7)
   t = np.array([0.0] + list(np.cumsum(np.hypot(np.diff(bx), np.diff(by)))))
   tq = np.linspace(t[0], t[-1], 400)
   bx_lin = np.interp(tq, t, bx)
   by_lin = np.interp(tq, t, by)
   bx_card = hermite_eval(tq, t, bx, cardinal_tangents(t, bx))
   by_card = hermite_eval(tq, t, by, cardinal_tangents(t, by))
   ax.plot(bx_lin, by_lin, label="LINEAR")
   ax.plot(bx_card, by_card, label="CARDINAL")
   ax.plot(bx, by, "ko", ms=5, label="vertices")
   ax.set_title("Shape preservation on a curved bank")
   ax.set_xlabel("X"); ax.set_ylabel("Y")
   ax.set_aspect("equal")
   ax.legend(fontsize=9); ax.grid(alpha=0.3)

   fig.tight_layout()

.. _bed-linear-interp:

Linear bed interpolation between profiles
-----------------------------------------

Implemented in :meth:`tatooinemesher.section.Bed.interp_coord_linear`.

Inside one bed bounded by two constraint lines, at intermediate longitudinal
position :math:`x_l \in [0, 1]` and lateral position :math:`x_t \in [0, 1]`,
node coordinates are blended linearly between the two profiles:

.. math::

   x(x_l, x_t) = (1 - x_l)\, x_{\text{us}}(x_t) + x_l\, x_{\text{ds}}(x_t),

.. math::

   y(x_l, x_t) = (1 - x_l)\, y_{\text{us}}(x_t) + x_l\, y_{\text{ds}}(x_t).

The lateral functions :math:`x_{\text{us}}, y_{\text{us}}, x_{\text{ds}},
y_{\text{ds}}` are themselves piecewise linear in :math:`X_t` along their
respective bed.

.. plot::
   :caption: Schematic of a single submesh in flow-oriented coordinates
             :math:`(x_l, x_t)`. Constraint lines are the lateral edges
             :math:`x_t \in \{0, 1\}`; cross-sections are the longitudinal
             edges :math:`x_l \in \{0, 1\}`. Each interior node :math:`P` is
             interpolated from its 4 neighbours on the unit square.
   :include-source: false

   import numpy as np
   import matplotlib.pyplot as plt

   fig, ax = plt.subplots(figsize=(6, 5))

   nx, nt = 6, 5
   xl = np.linspace(0, 1, nx)
   xt = np.linspace(0, 1, nt)
   Xl, Xt = np.meshgrid(xl, xt)

   for i in range(nx):
       ax.plot([xl[i]] * nt, xt, color="0.7", lw=0.6)
   for j in range(nt):
       ax.plot(xl, [xt[j]] * nx, color="0.7", lw=0.6)
   ax.plot(Xl.ravel(), Xt.ravel(), "ko", ms=3)

   for x in [0, 1]:
       ax.plot([x, x], [0, 1], "g-", lw=2.5)
   for y in [0, 1]:
       ax.plot([0, 1], [y, y], "r-", lw=2.5)

   ax.plot(0.55, 0.65, "bs", ms=11)
   ax.annotate("P", (0.55, 0.65), xytext=(0.61, 0.71), fontsize=13, color="b")

   ax.plot(0.55, 0, "k^", ms=8)
   ax.annotate("A", (0.55, 0), xytext=(0.58, 0.04), fontsize=11)
   ax.plot(0.55, 1, "k^", ms=8)
   ax.annotate("B", (0.55, 1), xytext=(0.58, 0.93), fontsize=11)

   ax.set_xlabel(r"$x_l$ (longitudinal)")
   ax.set_ylabel(r"$x_t$ (lateral)")
   ax.set_xlim(-0.05, 1.05); ax.set_ylim(-0.05, 1.05)
   ax.set_aspect("equal")
   ax.set_title("Flow-oriented coordinate system (single bed)")

   from matplotlib.lines import Line2D
   handles = [
       Line2D([0], [0], color="g", lw=2.5, label="Constraint lines"),
       Line2D([0], [0], color="r", lw=2.5, label="Cross-sections"),
       Line2D([0], [0], marker="o", color="0.4", lw=0, label="Mesh nodes"),
       Line2D([0], [0], marker="s", color="b", lw=0, label="Query node P"),
       Line2D([0], [0], marker="^", color="k", lw=0, label="A, B (lateral interp)"),
   ]
   ax.legend(handles=handles, fontsize=8, loc="upper left")
   fig.tight_layout()

.. _values-interp:

Cross-section bathymetry interpolation
--------------------------------------

Implemented in
:meth:`tatooinemesher.mesh_constructor.MeshConstructor.interp_1d_values_from_profiles`
and
:meth:`tatooinemesher.mesh_constructor.MeshConstructor.interp_2d_values_from_profiles`.

1D longitudinal blend
^^^^^^^^^^^^^^^^^^^^^

For variable :math:`v` (Z, Manning, etc.) at a node of dimensionless
position :math:`(x_l, x_t)`:

.. math::

   v(x_l, x_t) = (1 - x_l)\, v_{\text{us}}\!\bigl(X_{t,\text{us}}(x_t)\bigr)
              + x_l\, v_{\text{ds}}\!\bigl(X_{t,\text{ds}}(x_t)\bigr),

where :math:`v_{\text{us}}, v_{\text{ds}}` are the cross-section profiles
of :math:`v` along their respective :math:`X_t`. The lateral interpolator
applied to :math:`v_{\text{us}}` and :math:`v_{\text{ds}}` depends on
``interp_values``:

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Mode
     - Lateral interpolator
   * - ``LINEAR``
     - :func:`numpy.interp`
   * - ``B-SPLINE``
     - :func:`scipy.interpolate.splrep` + :func:`scipy.interpolate.splev`
       (cubic B-spline)
   * - ``AKIMA``
     - :class:`scipy.interpolate.Akima1DInterpolator`
   * - ``PCHIP``
     - :class:`scipy.interpolate.PchipInterpolator`
   * - ``CUBIC_SPLINE``
     - :class:`scipy.interpolate.CubicSpline`

.. _lateral-interp-comparison:

.. plot::
   :caption: Comparison of the five 1D lateral interpolators on a synthetic
             cross-section with a low thalweg and a sharp inner-bank break.
             ``CUBIC_SPLINE`` overshoots near the discontinuity; ``AKIMA``
             and ``PCHIP`` track the data without wiggles.
   :include-source: false

   import numpy as np
   import matplotlib.pyplot as plt
   from scipy.interpolate import (
       Akima1DInterpolator,
       CubicSpline,
       PchipInterpolator,
       splev,
       splrep,
   )

   xt_data = np.array([0.0, 0.10, 0.20, 0.25, 0.35, 0.55, 0.75, 0.85, 1.0])
   z_data = np.array([305.0, 304.5, 303.5, 299.0, 299.2, 299.5, 301.0, 303.0, 305.0])

   xt = np.linspace(0, 1, 400)
   tck = splrep(xt_data, z_data, k=3)

   fig, ax = plt.subplots(figsize=(8, 4))
   ax.plot(xt_data, z_data, "kx", ms=10, mew=2, label="Cross-section data", zorder=10)
   ax.plot(xt, np.interp(xt, xt_data, z_data), label="LINEAR")
   ax.plot(xt, Akima1DInterpolator(xt_data, z_data)(xt), "--", label="AKIMA")
   ax.plot(xt, PchipInterpolator(xt_data, z_data)(xt), "-.", label="PCHIP")
   ax.plot(xt, CubicSpline(xt_data, z_data)(xt), ":", label="CUBIC_SPLINE")
   ax.plot(xt, splev(xt, tck), label="B-SPLINE")
   ax.set_xlabel(r"$x_t$ (lateral, dimensionless)")
   ax.set_ylabel(r"$Z$ (m)")
   ax.legend(loc="upper center", ncol=3, fontsize=9)
   ax.grid(alpha=0.3)
   ax.set_title("Lateral interpolation modes — sample cross-section")
   fig.tight_layout()

The decomposition into a lateral interpolator followed by a longitudinal
linear blend is widely used in the river-bed reconstruction literature
([Caviedes2014]_, [Schäppi2009]_); TatooineMesher exposes the choice of
the lateral kernel to the modeller while keeping the longitudinal blend
linear by default.

2D bivariate interpolation
^^^^^^^^^^^^^^^^^^^^^^^^^^

Selected by an ``interp_values`` value starting with ``BI``. The
cross-section samples of every zone are mapped to the local
:math:`(u, v) = (x_l, x_t)` square and a 2D interpolator is fit:

* ``BILINEAR`` / ``BICUBIC`` — :func:`scipy.interpolate.griddata` with
  ``method = "linear"`` or ``"cubic"``.
* ``BIVARIATE_SPLINE`` —
  :class:`scipy.interpolate.SmoothBivariateSpline` with ``kx = ky = 3``.

These modes require **exactly two** constraint lines (a single bed across
the river width).

.. _z-correction:

3D constraint line Z correction
-------------------------------

Implemented in
:meth:`tatooinemesher.mesh_constructor.MeshConstructor._compute_line_z_correction`
and applied in
:meth:`tatooinemesher.mesh_constructor.MeshConstructor.build_interp`.

Activated when at least one constraint line is loaded with ``has_z=True``.
The goal is to correct the longitudinal bathymetry between two cross-sections
so that it follows the Z profile of a 3D constraint line (typically the
talweg or a sharp bank), without altering Z at the section locations.

For a bed bounded by lines :math:`L_1, L_2`, at intermediate longitudinal
position :math:`x_l \in [0, 1]`, define along each line:

**Line Z evaluated at the intermediate** :math:`x_l`:
:math:`Z_{\text{line}}(x_l) = Z_{L_j}(X_{p,\text{us}} + x_l \, \Delta X_p)`
(via the line's own interpolator —
:meth:`~tatooinemesher.constraint_line.ConstraintLine.coord_sampling_along_line`).

**Linear cord** of the line Z between the two sections:

.. math::

   Z_{\text{line,lin}}(x_l) = (1 - x_l)\, Z_{L_j}(X_{p,\text{us}})
                            + x_l\, Z_{L_j}(X_{p,\text{ds}}).

**Cross-section Z** at the line's intersection with each section, linearly
blended:

.. math::

   Z_{\text{sec,lin}}^{(L_j)}(x_l) = (1 - x_l)\, Z_{\text{sec,us}}^{(L_j)}
                                   + x_l\, Z_{\text{sec,ds}}^{(L_j)}.

**Longitudinal gap** between the line and the section bathymetry:

.. math::

   \mathrm{gap}_j(x_l)
       = Z_{\text{line,lin}}(x_l) - Z_{\text{sec,lin}}^{(L_j)}(x_l).

**Shape deviation** (line departure from its own cord — vanishes at
:math:`x_l \in \{0, 1\}` so cross-section Z is preserved at the sections):

.. math::

   \mathrm{dev}_j(x_l)
       = Z_{\text{line}}(x_l) - Z_{\text{line,lin}}(x_l).

.. _pull-bell:

**Pull bell curve** (vanishes at :math:`x_l \in \{0, 1\}`,
:math:`= \text{strength}` at :math:`x_l = 0.5`):

.. math::

   \alpha(x_l) = \text{strength} \cdot 4\, x_l\, (1 - x_l).

**Boundary correction** at the bed edge (along line :math:`L_j`):

.. math::

   \Delta Z_j(x_l) = \mathrm{dev}_j(x_l) + \alpha(x_l)\, \mathrm{gap}_j(x_l).

.. plot::
   :caption: Pull factor :math:`\alpha(x_l) = \text{strength} \cdot 4\,
             x_l\,(1 - x_l)` for several ``z_line_strength`` values.
   :include-source: false

   import numpy as np
   import matplotlib.pyplot as plt

   xl = np.linspace(0, 1, 200)
   fig, ax = plt.subplots(figsize=(7, 3.5))
   for s in (0.0, 0.25, 0.5, 1.0):
       ax.plot(xl, s * 4 * xl * (1 - xl), label=f"strength = {s}")
   ax.set_xlabel(r"$x_l$"); ax.set_ylabel(r"$\alpha(x_l)$")
   ax.set_title("Pull bell curve")
   ax.legend(); ax.grid(alpha=0.3)
   fig.tight_layout()

.. _lateral-spread:

Lateral propagation through the bed
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The boundary corrections :math:`\Delta Z_1, \Delta Z_2` must now be spread
laterally through the bed at every :math:`x_t \in [0, 1]` between
:math:`L_1` and :math:`L_2`. The blend is a power blend whose exponents
grow with the local Z gap (controlled by ``z_line_gap_scale``):

.. math::

   p_{L_j} = 1 + \text{z\_line\_gap\_scale} \cdot |\mathrm{gap}_j(x_l)|,

.. math::

   \Delta Z(x_l, x_t) = (1 - x_t)^{p_{L_1}} \Delta Z_1(x_l)
                      + x_t^{p_{L_2}} \Delta Z_2(x_l).

With ``z_line_gap_scale = 0`` (default) we recover a plain linear blend
:math:`(1 - x_t) \Delta Z_1 + x_t \Delta Z_2`. With a positive value, the
exponent rises wherever the line Z diverges from the bed bathymetry,
keeping the correction localized near the offending line.

.. plot::
   :caption: Effect of ``z_line_gap_scale`` on the lateral weight
             :math:`(1 - x_t)^p` for various longitudinal gap magnitudes.
   :include-source: false

   import numpy as np
   import matplotlib.pyplot as plt

   xt = np.linspace(0, 1, 200)
   fig, axes = plt.subplots(1, 2, figsize=(9, 3.5), sharey=True)

   for ax, scale in zip(axes, (0.0, 1.0)):
       for gap in (0.0, 1.0, 5.0, 10.0):
           p = 1.0 + scale * abs(gap)
           ax.plot(xt, (1 - xt) ** p, label=f"|gap| = {gap}")
       ax.set_xlabel(r"$x_t$")
       ax.set_title(f"z_line_gap_scale = {scale}")
       ax.grid(alpha=0.3)
       ax.legend(fontsize=8)
   axes[0].set_ylabel(r"$(1 - x_t)^{p_{L_1}}$")
   fig.tight_layout()

The final correction is added to the first variable (Z) of every node in the
bed, in
:meth:`~tatooinemesher.mesh_constructor.MeshConstructor.interp_values_from_geom`:

.. math::

   Z_{\text{final}}(x_l, x_t) = Z_{\text{interp}}(x_l, x_t)
                              + \Delta Z(x_l, x_t),

where :math:`Z_{\text{interp}}` is the bathymetry interpolated by the chosen
:ref:`values-interp` mode.

End-to-end illustration
^^^^^^^^^^^^^^^^^^^^^^^

The plot below stitches the three terms together on a synthetic
thalweg case:

* :math:`L_1` is a 3D constraint line on the left bank with a 2 m sag at
  the mid-distance (:math:`Z_{L_1}` = 0 m at both cross-sections, but
  drops to −2 m at :math:`x_l = 0.5`);
* the cross-section bed at the :math:`L_1` intersection sits at
  :math:`Z_{\text{sec}}` = −1 m in both cross-sections (the cross-sections
  underestimate the thalweg depth);
* :math:`L_2` is a 2D line on the right bank (no Z information);
* the bed bathymetry between sections is a flat linear blend at −1 m.

At the mid-distance node :math:`x_l = 0.5`, the constants of the
correction become :math:`\mathrm{gap}_1 = +1\ \text{m}`,
:math:`\mathrm{dev}_1 = -2\ \text{m}`,
:math:`\Delta Z_1 = -2 + \alpha(0.5)`, and :math:`\Delta Z_2 = 0`.

.. plot::
   :caption: Lateral profile of the Z correction at :math:`x_l = 0.5`.
             **Top:** varying ``z_line_strength`` with ``z_line_gap_scale``
             = 0. At ``strength = 1`` the boundary node reaches the line
             Z exactly (Z_final − Z_interp = −1). **Bottom:** varying
             ``z_line_gap_scale`` with ``strength = 0.5``; larger values
             localize the correction near :math:`L_1`.
   :include-source: false

   import numpy as np
   import matplotlib.pyplot as plt

   xt = np.linspace(0, 1, 400)
   xl = 0.5

   z_line_us, z_line_ds = 0.0, 0.0
   z_line_at_xl = -2.0
   z_sec_us_L1, z_sec_ds_L1 = -1.0, -1.0

   z_line_lin = z_line_us * (1 - xl) + z_line_ds * xl
   z_sec_lin = z_sec_us_L1 * (1 - xl) + z_sec_ds_L1 * xl
   gap_L1 = z_line_lin - z_sec_lin
   dev_L1 = z_line_at_xl - z_line_lin
   gap_L2 = 0.0
   dev_L2 = 0.0

   def correction(strength, gap_scale):
       alpha = strength * 4 * xl * (1 - xl)
       dZ_L1 = dev_L1 + alpha * gap_L1
       dZ_L2 = dev_L2 + alpha * gap_L2
       p_L1 = 1.0 + gap_scale * abs(gap_L1)
       p_L2 = 1.0 + gap_scale * abs(gap_L2)
       return (1 - xt) ** p_L1 * dZ_L1 + xt ** p_L2 * dZ_L2

   fig, axes = plt.subplots(2, 1, figsize=(8, 6), sharex=True)

   ax = axes[0]
   for strength in (0.0, 0.25, 0.5, 1.0):
       ax.plot(xt, correction(strength, 0.0), label=f"strength = {strength}")
   ax.set_ylabel(r"$\Delta Z(x_l=0.5, x_t)$  [m]")
   ax.set_title("Effect of z_line_strength (gap_scale = 0)")
   ax.axhline(0, color="0.5", lw=0.5)
   ax.legend(fontsize=8, ncol=4)
   ax.grid(alpha=0.3)

   ax = axes[1]
   for gap_scale in (0.0, 0.5, 1.0, 2.0):
       ax.plot(xt, correction(0.5, gap_scale), label=f"gap_scale = {gap_scale}")
   ax.set_xlabel(r"$x_t$ (lateral)")
   ax.set_ylabel(r"$\Delta Z(x_l=0.5, x_t)$  [m]")
   ax.set_title(r"Effect of z_line_gap_scale (strength = 0.5)")
   ax.axhline(0, color="0.5", lw=0.5)
   ax.legend(fontsize=8, ncol=4)
   ax.grid(alpha=0.3)

   fig.tight_layout()
