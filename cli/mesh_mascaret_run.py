"""
mesh_mascaret_run.py

Mesh Mascaret geometry and results if specified
"""

import sys
from time import perf_counter

from crue10.utils import ExceptionCrue10

from tatooinemesher.algorithms.mesh_mascaret_run_alg import mesh_mascaret_run
from tatooinemesher.utils import logger, set_logger_level
from tatooinemesher.utils.arg_command_line import MyArgParse

parser = MyArgParse(description=__doc__)
parser.add_common_args(constant_long_disc=True)
# Inputs
parser.infile_args.title = "~> Input Mascaret files arguments"
parser.infile_args.add_argument("infile_geo", help="Mascaret geometry file (*.georef, *.georefC)")
parser.infile_args.add_argument("--infile_res", help="Mascaret results file (*.opt, *.rub)")
# Outputs
parser.add_out_mesh_file()


if __name__ == "__main__":
    args = parser.parse_args()
    try:
        set_logger_level(args.verbose)
        t1 = perf_counter()
        mesh_mascaret_run(
            args.infile_geo,
            args.long_step,
            args.outfile_mesh,
            infile_res=args.infile_res,
            lat_step=args.lat_step,
            nb_pts_lat=args.nb_pts_lat,
            interp_constraint_lines=args.interp_constraint_lines,
            interp_values=args.interp_values,
            constant_long_disc=args.constant_long_disc,
            dist_max=args.dist_max,
            lang=args.lang,
            verbose=args.verbose,
        )
        t2 = perf_counter()
        logger.info("=> Execution time: {}s".format(t2 - t1))
    except FileNotFoundError as e:
        logger.critical(e)
        sys.exit(1)
    except ExceptionCrue10 as e:
        logger.critical(e)
        sys.exit(2)
