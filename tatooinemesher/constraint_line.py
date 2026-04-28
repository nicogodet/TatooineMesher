import numpy as np
import shapefile
from shapely.geometry import LineString

try:
    from pyteltools.geom import BlueKenue as bk
    from pyteltools.geom import Shapefile as shp
except ImportError:
    bk = None
    shp = None

from tatooinemesher.interp.cubic_hermite_spline import CubicHermiteSpline
from tatooinemesher.utils import TatooineException, float_vars


class ConstraintLine:
    """
    ConstraintLine: open polyline used to separate beds

    ### Attributes
    - id <integer>:  unique identifier (automatic numbering starting from 0)
    - nb_points <int>: number of points
    - coord <2D-array float>: coordinates (X and Y)
    - Z <1D-array float | None>: elevation along the line (None when has_z is False)
    - has_z <bool>: True when an independent Z profile is attached to the line
    - geom <shapely.geometry.LineString>: 2D geometry
    - interp: coordinates (and optionally Z) interpolator

    ### Methods
    - get_lines_from_file
    - get_lines_and_set_limits_from_sections
    - build_interp_linear
    - build_interp_chs
    - coord_sampling_along_line
    """

    def __init__(self, id, coord, interp_coord="LINEAR", has_z=False):
        """
        Create a ConstraintLine instance from X and Y coordinates (and optionally Z)
        @id <integer>: unique identifier (automatic numbering starting from 0)
        @coord: coordinates; 2D tuples (x, y) or 3D tuples (x, y, z) when has_z is True
        @interp_coord <str>: coordinates interpolation method (among: 'LINEAR', 'CARDINAL' and 'FINITE_DIFF')
        @has_z <bool>: if True, coord must provide a Z value that is interpolated along the line
        """
        self.id = id
        self.nb_points = len(coord)
        self.has_z = has_z
        coord_arr = np.asarray(coord, dtype=float)
        if has_z:
            if coord_arr.ndim != 2 or coord_arr.shape[1] < 3:
                raise TatooineException(
                    "ConstraintLine #%s: 3D coordinates (x, y, z) are required when has_z=True" % id
                )
            self.Z = coord_arr[:, 2].copy()
            self.coord = coord_arr[:, :2].copy()
        else:
            self.Z = None
            self.coord = coord_arr[:, :2].copy() if coord_arr.ndim == 2 else coord_arr
        Xt = np.sqrt(
            np.power(np.ediff1d(self.coord[:, 0], to_begin=0.0), 2)
            + np.power(np.ediff1d(self.coord[:, 1], to_begin=0.0), 2)
        )
        self.Xt = Xt.cumsum()
        self.geom = LineString(self.coord)

        if interp_coord == "LINEAR":
            self.interp = self.build_interp_linear()
        else:
            if interp_coord == "CARDINAL":
                tan_method = CubicHermiteSpline.CARDINAL
            elif interp_coord == "FINITE_DIFF":
                tan_method = CubicHermiteSpline.FINITE_DIFF
            else:
                raise NotImplementedError
            self.interp = self.build_interp_chs(tan_method)

    def __repr__(self):
        return f"ConstraintLine #{self.id} ({self.nb_points} points)"

    @staticmethod
    def get_lines_from_file(filename, interp_coord="LINEAR", has_z=False):
        """
        Returns a list of ConstraintLine from an input file
        @param has_z <bool>: if True, read Z along the line (POLYLINEZ shapefile required)
        TODO 1: Value is ignored in i2s file format
        """
        lines = []
        if filename is not None:
            if filename.endswith(".i2s"):
                if has_z:
                    raise TatooineException("i2s format does not carry Z; use a POLYLINEZ shapefile instead")
                if bk is None:
                    raise ImportError(
                        "PyTelTools is required for this feature. "
                        "Install it with: pip install TatooineMesher[pyteltools]"
                    )
                with bk.Read(filename) as in_i2s:
                    in_i2s.read_header()
                    for i, line in enumerate(in_i2s.get_open_polylines()):
                        lines.append(ConstraintLine(i, list(line.polyline().coords), interp_coord))

            elif filename.endswith(".shp"):
                if shp is None:
                    raise ImportError(
                        "PyTelTools is required for this feature. "
                        "Install it with: pip install TatooineMesher[pyteltools]"
                    )
                shp_type = shp.get_shape_type(filename)
                if has_z:
                    if shp_type != shapefile.POLYLINEZ:
                        raise TatooineException(
                            f"Z-aware constraint lines require POLYLINEZ; file {filename} is of a different type"
                        )
                    for i, line in enumerate(shp.get_open_polylines(filename)):
                        lines.append(ConstraintLine(i, list(line.polyline().coords), interp_coord, has_z=True))
                else:
                    if shp_type not in (shapefile.POLYLINE, shapefile.POLYLINEZ, shapefile.POLYLINEM):
                        raise TatooineException(f"The type of file {filename} is not POLYLINEZ[M]")
                    for i, line in enumerate(shp.get_open_polylines(filename)):
                        lines.append(ConstraintLine(i, list(line.polyline().coords), interp_coord))

            else:
                raise NotImplementedError("Only shp and i2s formats are supported for constraint lines")

        return lines

    @staticmethod
    def get_lines_and_set_limits_from_sections(section_seq, interp_coord="LINEAR"):
        """
        @brief: Returns a list of ConstraintLine from an sequence of cross-sections
        @param section_seq <CrossSectionSequence>: sequence of cross-sections
        @param interp_coord <str>: interpolation method
        """
        # Build 2 constraint lines from cross-section bounds
        first_coords = []
        last_coords = []
        for section in section_seq:
            first_coords.append(section.geom.coords[0][:2])
            last_coords.append(section.geom.coords[-1][:2])
        lines = [ConstraintLine(0, first_coords, interp_coord), ConstraintLine(1, last_coords, interp_coord)]

        # Set limits
        for line_id, line in enumerate(lines):
            for section, Xt_line in zip(section_seq, line.Xt):
                Xt_section = section.coord.array["Xt"][0] if line_id == 0 else section.coord.array["Xt"][-1]
                intersection = section.geom.interpolate(Xt_section)
                section.add_limit(line_id, Xt_section, Xt_line, intersection)

        return lines

    def build_interp_linear(self):
        """
        @brief: Build a double linear interpolator for X and Y coordinates (and Z when has_z)
        """
        fields = ["X", "Y", "Z"] if self.has_z else ["X", "Y"]

        def interp_xy_linear(Xt_new):
            Xt_new = np.asarray(Xt_new)
            result = np.empty(len(Xt_new), dtype=float_vars(fields))
            for i, dist in enumerate(Xt_new):
                point = self.geom.interpolate(dist)
                result["X"][i] = point.coords[0][0]
                result["Y"][i] = point.coords[0][1]
            if self.has_z:
                result["Z"] = np.interp(Xt_new, self.Xt, self.Z)
            return result

        return interp_xy_linear

    def build_interp_chs(self, tan_method):
        """
        @brief: Build a Cubic Hermite Spline interpolator for X and Y coordinates (and Z when has_z)
        """
        spline_x = CubicHermiteSpline()  # x = spline_x(Xt)
        spline_y = CubicHermiteSpline()  # y = spline_y(Xt)

        spline_x.Initialize(np.vstack((self.Xt, self.coord[:, 0])).T, tan_method=tan_method)
        spline_y.Initialize(np.vstack((self.Xt, self.coord[:, 1])).T, tan_method=tan_method)

        spline_z = None
        if self.has_z:
            spline_z = CubicHermiteSpline()
            spline_z.Initialize(np.vstack((self.Xt, self.Z)).T, tan_method=tan_method)

        fields = ["X", "Y", "Z"] if self.has_z else ["X", "Y"]

        def interp_xy_chs(Xt_new):
            Xt_new = np.asarray(Xt_new)
            result = np.empty(len(Xt_new), dtype=float_vars(fields))
            for i, dist in enumerate(Xt_new):
                result["X"][i] = spline_x.evaluate(dist)
                result["Y"][i] = spline_y.evaluate(dist)
                if spline_z is not None:
                    result["Z"][i] = spline_z.evaluate(dist)
            return result

        return interp_xy_chs

    def coord_sampling_along_line(self, Xp1, Xp2, Xp_adm_int):
        """
        @brief: Extract coordinates interpolated along line between Xp1 and Xp2
            at requested dimensionless curvilinear distances
        @param: Xp_adm_int <1D-array>: dimensionless curvilinear distances between Xp1 and Xp2 (0=Xp1 and 1=Xp2)
        @param: Xp1 <float>: starting curvilinear distance
        @param: Xp2 <float>: ending curvilinear distance
        """
        # Building list of curvilinear distance in meters
        Xp = (1 - Xp_adm_int) * Xp1 + Xp_adm_int * Xp2
        # Use coordinate interpolator
        return self.interp(Xp)
