#!/usr/bin/env python3
"""
densify_cross_sections.py

Interpolate initial and intermediate cross-sections
Multiple variables are supported if input cross-sections file is a shapefile with POINTZ type.
"""
from time import perf_counter

from tatooinemesher.algorithms.densify_cross_sections_alg import densify_cross_sections
from tatooinemesher.utils import TatooineException, logger, set_logger_level
from tatooinemesher.utils.arg_command_line import MyArgParse

parser = MyArgParse(description=__doc__)
parser.add_common_args(project_straight_line=True)
# Inputs
parser.infile_args.add_argument("infile_axis", help="hydraulic axis file (*.shp, *.i2s)")
parser.infile_args.add_argument("infile_cross_sections", help="cross-sections file (*.shp, *.i3s)")
parser.infile_args.add_argument("--infile_constraint_lines", help="constraint lines file (*.shp, *.i2s)")
parser.infile_args.add_argument("--attr_cross_sections", help="attribute to identify cross-sections")
# Outputs
parser.outfile_args.add_argument(
    "outfile_sections", help="output file containing interpolated cross-sections (*.i3s, *.georefC, *.shp)"
)


if __name__ == "__main__":
    args = parser.parse_args()
    try:
        set_logger_level(args.verbose)
        t1 = perf_counter()
        densify_cross_sections(
            args.infile_axis,
            args.infile_cross_sections,
            args.attr_cross_sections,
            args.long_step,
            args.outfile_sections,
            infile_constraint_lines=args.infile_constraint_lines,
            interp_constraint_lines=args.interp_constraint_lines,
            interp_values=args.interp_values,
            project_straight_line=args.project_straight_line,
            nb_pts_lat=args.nb_pts_lat,
            lat_step=args.lat_step,
            dist_max=args.dist_max,
            verbose=args.verbose,
        )
        t2 = perf_counter()
        logger.info("=> Execution time: {}s".format(t2 - t1))
    except TatooineException as e:
        logger.critical(e.message)
