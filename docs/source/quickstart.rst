Quickstart
==========

The ``mesh_and_interpolate`` algorithm is exposed both as a Python function
(:func:`tatooinemesher.algorithms.mesh_and_interpolate_alg.mesh_and_interpolate`)
and as a CLI script (``cli/mesh_and_interpolate.py``).

2D case (no Z on constraint lines)
----------------------------------

.. code-block:: bash

   python cli/mesh_and_interpolate.py \
       axis.shp \
       cross_sections.shp \
       --infile_constraint_lines banks.shp \
       --long_step 5.0 \
       --lat_step 1.0 \
       --interp_values LINEAR \
       --outfile_mesh mesh.slf

* ``axis.shp``: hydraulic axis (single open polyline)
* ``cross_sections.shp``: cross-sections (POLYLINEZ or POLYLINE)
* ``banks.shp``: two open polylines defining the river bed boundaries

3D case — bathymetry correction from Z-aware constraint lines
-------------------------------------------------------------

When the bed bottom (talweg) or a sharp lateral feature is poorly captured
between cross-sections, you can supply a **POLYLINEZ** shapefile whose Z values
guide the longitudinal bathymetry between sections:

.. code-block:: bash

   python cli/mesh_and_interpolate.py \
       axis.shp \
       cross_sections.shp \
       --infile_constraint_lines banks.shp \
       --infile_constraint_3D_lines talweg.shp \
       --z_line_strength 0.5 \
       --z_line_gap_scale 0.0 \
       --long_step 5.0 \
       --lat_step 1.0 \
       --outfile_mesh mesh.slf

* ``--z_line_strength`` (default 0): blends the cross-section Z (0) toward the
  3D line Z (1) at the mid-distance between two consecutive sections. Section
  Z is preserved at the section locations regardless. See
  :ref:`the pull bell curve <pull-bell>`.
* ``--z_line_gap_scale`` (default 0): localizes the lateral propagation of the
  correction when the line Z diverges sharply from the bed bathymetry. See
  :ref:`lateral-spread`.

Python API
----------

.. code-block:: python

   from tatooinemesher.algorithms import mesh_and_interpolate

   mesh_constr = mesh_and_interpolate(
       infile_axis="axis.shp",
       infile_cross_sections="cross_sections.shp",
       attr_cross_sections=None,
       long_step=5.0,
       infile_constraint_lines="banks.shp",
       infile_constraint_3D_lines="talweg.shp",
       z_line_strength=0.5,
       lat_step=1.0,
       interp_values="LINEAR",
   )
   mesh_constr.export_mesh("mesh.slf", lang="en")

For the full parameter list, see :doc:`parameters`.
