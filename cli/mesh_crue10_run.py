"""
mesh_crue10_run.py

Génération d'un fichier résultat 2D (format Telemac) avec plusieurs variables et éventuellement plusieurs temps.
La taille des éléments du maillage est controlable.

Si un fichier rcal est spécifié alors on peut traiter :
- tous les calculs permanents
- un calcul transitoire en spécifiant son nom dans l'argument `--calc_unsteady`

Les variables écrites dans le fichier de sortie sont :
* FOND
* IS LIT ACTIVE (0 = lit inactif, 1 = lit actif)
* FROTTEMENT (moyenne sur la verticale si plusieurs valeurs)
* Les variables supplémentaires si un résultat est lu sont :
    * HAUTEUR D'EAU (la variable 'Z' aux sections/casiers est nécessaire.
        Attention, il ne faut pas avoir sorti la charge 'H' aux sections)
    * VITESSE SCALAIRE (seulement si la variable 'Vact' est présente)

Seulement les branches et les casiers actifs sont traités.
"""

import sys
from time import perf_counter

from crue10.emh.branche import Branche
from crue10.utils import ExceptionCrue10

from tatooinemesher.algorithms.mesh_crue10_run_alg import mesh_crue10_run
from tatooinemesher.utils import logger, set_logger_level
from tatooinemesher.utils.arg_command_line import MyArgParse

parser = MyArgParse(description=__doc__)
parser.add_common_args(project_straight_line=False, constant_long_disc=True)
# Inputs
parser.infile_args.title = "~> Crue10 input model and run (and the optional DEM)"
parser.infile_args.add_argument("infile_etu", help="Crue10 study file (*.etu.xml)")
parser.infile_args.add_argument("model_name", help="model name")
parser.infile_args.add_argument("--infile_rcal", help="Crue10 results file (*.rcal.xml)")
parser.infile_args.add_argument(
    "--calc_unsteady", help="name of the unsteady file (otherwise considers all steady calculations)"
)
parser.infile_args.add_argument(
    "--infile_dem",
    help='Raster file (geoTIFF format) containing bottom elevation for the "casiers" in the floodplain (*.tif)',
)
# Parameters to select branches
parser_branches = parser.add_argument_group("Parameters to filter branches")
parser_branches.add_argument(
    "--branch_types_filter",
    type=int,
    nargs="+",
    default=Branche.TYPES_IN_MINOR_BED,
    help="types of branches to consider",
)
parser_branches.add_argument(
    "--branch_patterns",
    nargs="+",
    default=None,
    help="list of patterns to filter branches which name does not contain any pattern",
)
# Mesh parameters
parser.mesher_args.add_argument("--floodplain_step", type=float, default=None, help="floodplain space step (in m)")
# Outputs
parser.add_out_mesh_file()


if __name__ == "__main__":
    args = parser.parse_args()
    try:
        set_logger_level(args.verbose)
        t1 = perf_counter()
        mesh_crue10_run(
            args.infile_etu,
            args.model_name,
            args.long_step,
            args.outfile_mesh,
            infile_rcal=args.infile_rcal,
            calc_unsteady=args.calc_unsteady,
            infile_dem=args.infile_dem,
            branch_types_filter=args.branch_types_filter,
            branch_patterns=args.branch_patterns,
            lat_step=args.lat_step,
            nb_pts_lat=args.nb_pts_lat,
            interp_constraint_lines=args.interp_constraint_lines,
            interp_values=args.interp_values,
            constant_long_disc=args.constant_long_disc,
            dist_max=args.dist_max,
            floodplain_step=args.floodplain_step,
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
