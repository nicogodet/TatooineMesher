#!/usr/bin/env python3
"""
mesh_and_interpolate.py

Channel mesher and/or interpolator from cross-sections
"""
# TODO: integrate "seuils?" (local correction of the bathymetry)
from time import perf_counter

from tatooinemesher.algorithms.mesh_and_interpolate_alg import mesh_and_interpolate
from tatooinemesher.utils import TatooineException, logger, set_logger_level
from tatooinemesher.utils.arg_command_line import MyArgParse

parser = MyArgParse(description=__doc__)
parser.add_common_args(project_straight_line=True, constant_long_disc=True)
# Inputs
parser.infile_args.add_argument("infile_axis", help="hydraulic axis file (*.shp, *.i2s)")
parser.infile_args.add_argument("infile_cross_sections", help="cross-sections file (*.shp, *.i3s)")
parser.infile_args.add_argument("--infile_constraint_lines", help="constraint lines file (*.shp, *.i2s)")
parser.infile_args.add_argument("--attr_cross_sections", help="attribute to identify cross-sections")
# TODO: add groynes
# parser_epis = parser.add_argument_group('Parameters to define lateral groynes (optional)')
# parser_epis.add_argument("--infile_epis", help="input file for groynes (*.shp, *.i3s)")
# parser_epis.add_argument("--attr_epis", help="attribute to identify groynes")
# parser_epis.add_argument("--dist_corr_epi", type=float,
#                          help="distance around groynes to modify nodes close to them "
#                               "(should be less than lateral and longitudinal space step)")
# Outputs
parser.add_out_mesh_file(is_optional=True)
parser.outfile_args.add_argument("--outfile_nodes", help="output points set file with mesh nodes (*.shp, *.xyz)")


if __name__ == "__main__":
    args = parser.parse_args()
    try:
        set_logger_level(args.verbose)
        t1 = perf_counter()
        mesh_and_interpolate(
            args.infile_axis,
            args.infile_cross_sections,
            args.attr_cross_sections,
            args.long_step,
            infile_constraint_lines=args.infile_constraint_lines,
            interp_constraint_lines=args.interp_constraint_lines,
            interp_values=args.interp_values,
            project_straight_line=args.project_straight_line,
            nb_pts_lat=args.nb_pts_lat,
            lat_step=args.lat_step,
            constant_long_disc=args.constant_long_disc,
            dist_max=args.dist_max,
            outfile_nodes=args.outfile_nodes,
            outfile_mesh=args.outfile_mesh,
            lang=args.lang,
            verbose=args.verbose,
        )
        t2 = perf_counter()
        logger.info("=> Execution time: {}s".format(t2 - t1))
    except TatooineException as e:
        logger.critical(e.message)
