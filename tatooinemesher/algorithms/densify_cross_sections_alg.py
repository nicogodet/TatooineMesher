from typing import Optional

from tatooinemesher.constraint_line import ConstraintLine
from tatooinemesher.mesh_constructor import MeshConstructor
from tatooinemesher.section import CrossSectionSequence
from tatooinemesher.utils import TatooineException, get_hydraulic_axis, logger, set_logger_level


def densify_cross_sections(
    infile_axis: str,
    infile_cross_sections: str,
    attr_cross_sections: Optional[str],
    long_step: float,
    outfile_sections: str,
    infile_constraint_lines: Optional[str] = None,
    interp_constraint_lines: str = "LINEAR",
    interp_values: str = "LINEAR",
    project_straight_line: bool = False,
    nb_pts_lat: Optional[int] = None,
    lat_step: Optional[float] = None,
    dist_max: float = 0.01,
    verbose: bool = False,
) -> None:
    """Interpolate initial and intermediate cross-sections.

    Multiple variables are supported if input cross-sections file is a shapefile
    with POINTZ type.

    Parameters
    ----------
    infile_axis:
        Hydraulic axis file (*.shp, *.i2s).
    infile_cross_sections:
        Cross-sections file (*.shp, *.i3s).
    attr_cross_sections:
        Attribute to identify cross-sections.
    long_step:
        Longitudinal space step (in m).
    outfile_sections:
        Output file containing interpolated cross-sections (*.i3s, *.georefC, *.shp).
    infile_constraint_lines:
        Constraint lines file (*.shp, *.i2s).
    interp_constraint_lines:
        Interpolation method for X and Y coordinates of constraint lines.
    interp_values:
        Interpolation method (crosswise for 1D or global 2D) for values.
    project_straight_line:
        Project cross-sections along a straight line (linking cross-section bounds).
    nb_pts_lat:
        Number of nodes crosswise.
    lat_step:
        Lateral space step (in m).
    dist_max:
        Maximum search distance to rescue intersections for limits (in m).
    verbose:
        Increase output verbosity.
    """
    set_logger_level(verbose)

    logger.info("~> Reading input files")
    axe = get_hydraulic_axis(infile_axis)
    section_seq = CrossSectionSequence.from_file(
        infile_cross_sections,
        "Cross-section",
        field_id=attr_cross_sections,
        project_straight_line=project_straight_line,
    )

    section_seq.compute_dist_proj_axe(axe, dist_max)
    section_seq.check_intersections()
    section_seq.sort_by_dist()

    if infile_constraint_lines is None:
        constraint_lines = ConstraintLine.get_lines_and_set_limits_from_sections(section_seq, interp_constraint_lines)
    else:
        constraint_lines = ConstraintLine.get_lines_from_file(infile_constraint_lines, interp_constraint_lines)
        if nb_pts_lat is not None and len(constraint_lines) != 2:
            raise TatooineException("Argument `--nb_pts_lat` is only compatible with 2 constraint lines!")
        if interp_values.startswith("BI") and len(constraint_lines) != 2:
            raise TatooineException("A 2D interpolation is only compatible with 2 constraint lines!")
        section_seq.find_and_add_limits(constraint_lines, dist_max)

    # section_seq.export_sections_shp('export_cross-sections.shp')

    mesh_constr = MeshConstructor(
        section_seq=section_seq, lat_step=lat_step, nb_pts_lat=nb_pts_lat, interp_values=interp_values
    )
    mesh_constr.build_interp(constraint_lines, long_step, True)
    mesh_constr.export_sections(outfile_sections)
