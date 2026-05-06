Installation
============

Requirements
------------

* Python ≥ 3.10 (tested up to 3.13)
* System libraries: ``libgdal-dev``, ``libspatialindex-dev`` (for the optional
  GDAL and shapely-rtree integrations)

Core install
------------

From a fresh clone:

.. code-block:: bash

   git clone https://github.com/CNR-Engineering/TatooineMesher.git
   cd TatooineMesher
   pip install -e .

This pulls only the core dependencies (numpy, scipy, shapely, matplotlib,
pyshp, jinja2, triangle).

Optional extras
---------------

Several integrations are gated behind extras to keep the core install light:

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Extra
     - Purpose
   * - ``crue10``
     - Pulls ``Crue10_tools`` (CNR-Engineering) for the ``mesh_crue10_run`` CLI
   * - ``pyteltools``
     - Pulls ``PyTelTools`` for ``.i2s`` / ``.i3s`` / ``.slf`` (Serafin) IO
   * - ``gdal``
     - Pulls the ``gdal`` python bindings for raster floor reading
   * - ``all``
     - Combines the three above
   * - ``docs``
     - Sphinx, Furo theme, ``myst-parser``, ``sphinx-copybutton``

Install one or several:

.. code-block:: bash

   pip install -e .[pyteltools]
   pip install -e .[all]
   pip install -e .[docs]

The ``mesh_and_interpolate`` algorithm only requires the core install plus
``pyteltools`` if you read shapefile constraint lines or write Serafin output.
