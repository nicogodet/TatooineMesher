import os.path
from typing import List, Optional

import numpy as np
import triangle
from crue10.emh.branche import Branche
from crue10.emh.section import SectionProfil
from crue10.etude import Etude
from crue10.run.resultats_calcul import ResultatsCalcul
from crue10.utils import ExceptionCrue10
from numpy.lib.recfunctions import unstructured_to_structured
from pyteltools.slf import Serafin

from tatooinemesher.constraint_line import ConstraintLine
from tatooinemesher.interp.raster import interp_raster
from tatooinemesher.mesh_constructor import MeshConstructor
from tatooinemesher.section import CrossSection, CrossSectionSequence
from tatooinemesher.utils import TatooineException, logger, resample_2d_line, set_logger_level

VARIABLES_FROM_GEOMETRY = ["B", "IS LIT ACTIVE", "W"]


def mesh_crue10_run(
    infile_etu: str,
    model_name: str,
    long_step: float,
    outfile_mesh: str,
    infile_rcal: Optional[str] = None,
    calc_unsteady: Optional[str] = None,
    infile_dem: Optional[str] = None,
    branch_types_filter: Optional[List[int]] = None,
    branch_patterns: Optional[List[str]] = None,
    lat_step: Optional[float] = None,
    nb_pts_lat: Optional[int] = None,
    interp_constraint_lines: str = "LINEAR",
    interp_values: str = "LINEAR",
    constant_long_disc: bool = False,
    dist_max: float = 0.01,
    floodplain_step: Optional[float] = None,
    lang: str = "fr",
    verbose: bool = False,
) -> None:
    """Mesh a Crue10 model and interpolate results if specified.

    Generates a 2D result file (Telemac format) with multiple variables and
    optionally multiple time steps. Mesh element size is controllable.

    Parameters
    ----------
    infile_etu:
        Crue10 study file (*.etu.xml).
    model_name:
        Model name.
    long_step:
        Longitudinal space step (in m).
    outfile_mesh:
        2D results file (*.slf, *.t3s, *.xml).
    infile_rcal:
        Crue10 results file (*.rcal.xml).
    calc_unsteady:
        Name of the unsteady calculation (otherwise considers all steady calculations).
    infile_dem:
        Raster file (geoTIFF) containing bottom elevation for casiers in the floodplain (*.tif).
    branch_types_filter:
        Types of branches to consider. Defaults to Branche.TYPES_IN_MINOR_BED.
    branch_patterns:
        List of patterns to filter branches which name does not contain any pattern.
    lat_step:
        Lateral space step (in m).
    nb_pts_lat:
        Number of nodes crosswise.
    interp_constraint_lines:
        Interpolation method for X and Y coordinates of constraint lines.
    interp_values:
        Interpolation method (crosswise for 1D or global 2D) for values.
    constant_long_disc:
        Method to compute number of intermediate cross-sections per zone
        (identical between 2 consecutive cross-sections) instead of per bed/submesh.
    dist_max:
        Maximum search distance to rescue intersections for limits (in m).
    floodplain_step:
        Floodplain space step (in m).
    lang:
        Language for standard variables in output file.
    verbose:
        Increase output verbosity.
    """
    set_logger_level(verbose)

    if branch_types_filter is None:
        branch_types_filter = Branche.TYPES_IN_MINOR_BED

    # Read the model and its submodels from xml/shp files
    etude = Etude(infile_etu)
    modele = etude.get_modele(model_name)
    modele.read_all()
    logger.info(modele)
    for sous_modele in modele.liste_sous_modeles:
        sous_modele.remove_sectioninterpolee()
        sous_modele.normalize_geometry()
        logger.info(sous_modele.summary())
        # sous_modele.write_shp_limites_lits_numerotes('limites_lits.shp')  # DEBUG
    logger.info(modele)

    global_mesh_constr = MeshConstructor()

    # Handle branches in minor bed
    for i, branche in enumerate(modele.get_liste_branches()):
        # Ignore branch if branch_patterns is set and do not match with current branch name
        if branch_patterns is not None:
            ignore = True
            for pattern in branch_patterns:
                if pattern in branche.id:
                    ignore = False
                    break
        else:
            ignore = False

        if branche.type not in branch_types_filter or not branche.is_active:
            ignore = True

        if not ignore:
            logger.info("===== TRAITEMENT DE LA BRANCHE %s =====" % branche.id)
            axe = branche.geom
            try:
                section_seq = CrossSectionSequence()
                for crue_section in branche.liste_sections_dans_branche:
                    if isinstance(crue_section, SectionProfil):
                        coords = list(crue_section.get_coord(add_z=True))
                        section = CrossSection(crue_section.id, [(coord[0], coord[1]) for coord in coords], "Section")

                        # Determine some variables (constant over the simulation) from the geometry
                        z = np.array([coord[2] for coord in coords])
                        is_bed_active = crue_section.get_is_bed_active_array()
                        mean_strickler = crue_section.get_friction_coeff_array()
                        section.coord.values = np.core.records.fromarrays(
                            np.column_stack((z, is_bed_active, mean_strickler)).T, names=VARIABLES_FROM_GEOMETRY
                        )

                        section_seq.add_section(section)

                section_seq.compute_dist_proj_axe(axe, dist_max)
                if len(section_seq) >= 2:
                    section_seq.check_intersections()
                    # section_seq.sort_by_dist() is useless because profiles are already sorted
                    constraint_lines = ConstraintLine.get_lines_and_set_limits_from_sections(
                        section_seq, interp_constraint_lines
                    )

                    mesh_constr = MeshConstructor(
                        section_seq=section_seq,
                        lat_step=lat_step,
                        nb_pts_lat=nb_pts_lat,
                        interp_values=interp_values,
                    )
                    mesh_constr.build_interp(constraint_lines, long_step, constant_long_disc)
                    mesh_constr.build_mesh(in_floworiented_crs=True)

                    global_mesh_constr.append_mesh_constr(mesh_constr)
                else:
                    logger.warning("Branche ignorée par manque de sections")
            except TatooineException as e:
                logger.error("/!\\ Branche ignorée à cause d'une erreur bloquante :")
                logger.error(e.message)
            logger.info("\n")

    # Handle casiers in floodplain
    nb_casiers = len(modele.get_liste_casiers())
    if infile_dem and nb_casiers > 0:
        logger.info("===== TRAITEMENT DES CASIERS =====")

        if not os.path.exists(infile_dem):
            raise TatooineException("File not found: %s" % infile_dem)
        from osgeo.gdal import Open

        raster = Open(infile_dem)
        dem_interp = interp_raster(raster)

        floodplain_step = floodplain_step if not None else long_step
        max_elem_area = floodplain_step * floodplain_step / 2.0
        simplify_dist = floodplain_step / 2.0

        for i, casier in enumerate(modele.get_liste_casiers()):
            if casier.is_active:
                if casier.geom is None:
                    raise TatooineException("Geometry of %s could not be found" % casier)
                line = casier.geom.simplify(simplify_dist)
                if not line.is_closed:
                    raise RuntimeError
                coords = resample_2d_line(line.coords, floodplain_step)[1:]  # Ignore last duplicated node

                hard_nodes_xy = np.array(coords, dtype=float)
                hard_nodes_idx = np.arange(0, len(hard_nodes_xy), dtype=int)
                hard_segments = np.column_stack((hard_nodes_idx, np.roll(hard_nodes_idx, 1)))

                tri = {
                    "vertices": np.array(np.column_stack((hard_nodes_xy[:, 0], hard_nodes_xy[:, 1]))),
                    "segments": hard_segments,
                }
                triangulation = triangle.triangulate(tri, opts="qpa%f" % max_elem_area)

                nodes_xy = np.array(triangulation["vertices"], dtype=float)
                bottom = dem_interp(nodes_xy)
                points = unstructured_to_structured(np.column_stack((nodes_xy, bottom)), names=["X", "Y", "Z"])

                global_mesh_constr.add_floodplain_mesh(triangulation, points)

    if len(global_mesh_constr.points) == 0:
        raise ExceptionCrue10(
            "Aucun point à traiter, adaptez l'option `--branch_patterns` et/ou `--branch_types_filter`"
        )

    logger.info(global_mesh_constr.summary())  # General information about the merged mesh

    if infile_rcal:
        # Read rcal result file
        resultats = ResultatsCalcul(infile_rcal)
        logger.info(resultats.summary())

        # Check result consistency
        missing_sections = modele.get_missing_active_sections(resultats.emh["Section"])
        if missing_sections:
            raise ExceptionCrue10(
                "Sections actives dans le scénario mais manquantes dans le Run :\n%s" % missing_sections
            )

        # Subset results to get requested variables at active sections
        varnames_1d = resultats.variables["Section"]
        logger.info("Variables 1D disponibles aux sections: %s" % varnames_1d)
        try:
            pos_z = varnames_1d.index("Z")
        except ValueError:
            raise TatooineException("La variable Z doit être présente dans les résultats aux sections")
        if global_mesh_constr.has_floodplain:
            try:
                pos_z_fp = resultats.variables["Casier"].index("Z")
            except ValueError:
                raise TatooineException("La variable Z doit être présente dans les résultats aux casiers")
        else:
            pos_z_fp = None

        pos_variables = [resultats.variables["Section"].index(var) for var in varnames_1d]
        pos_sections_list = [resultats.emh["Section"].index(profil.id) for profil in global_mesh_constr.section_seq]
        if global_mesh_constr.has_floodplain:
            pos_casiers_list = [
                resultats.emh["Casier"].index(casier.id) for casier in modele.get_liste_casiers() if casier.is_active
            ]
        else:
            pos_casiers_list = []

        if "H" in varnames_1d:
            raise TatooineException(
                "La variable H (charge) de Crue10 entre en conflit avec celle de TatooineMesher. "
                "Veuillez supprimer cette variable et relancer le traitement"
            )
        additional_variables_id = ["H"]
        if "Vact" in varnames_1d:
            additional_variables_id.append("M")

        values_geom = global_mesh_constr.interp_values_from_geom()
        z_bottom = values_geom[0, :]
        with Serafin.Write(outfile_mesh, lang, overwrite=True) as resout:
            title = "%s (written by TatooineMesher)" % os.path.basename(outfile_mesh)
            output_header = Serafin.SerafinHeader(title=title, lang=lang)
            output_header.from_triangulation(
                global_mesh_constr.triangle["vertices"], global_mesh_constr.triangle["triangles"] + 1
            )
            for var_name in VARIABLES_FROM_GEOMETRY:
                if var_name in ["B", "W"]:
                    output_header.add_variable_from_ID(var_name)
                else:
                    output_header.add_variable_str(var_name, var_name, "")
            for var_id in additional_variables_id:
                output_header.add_variable_from_ID(var_id)
            for var_name in varnames_1d:
                output_header.add_variable_str(var_name, var_name, "")
            resout.write_header(output_header)

            if calc_unsteady is None:
                for i, calc_name in enumerate(resultats.res_calc_pseudoperm.keys()):
                    logger.info("~> Calcul permanent %s" % calc_name)
                    # Read a single *steady* calculation
                    res_steady = resultats.get_data_pseudoperm(calc_name)
                    variables_at_profiles = res_steady["Section"][pos_sections_list, :][:, pos_variables]
                    if global_mesh_constr.has_floodplain:
                        z_at_casiers = res_steady["Casier"][pos_casiers_list, pos_z_fp]
                    else:
                        z_at_casiers = None

                    # Interpolate between sections and set in casiers
                    values_res = global_mesh_constr.interp_values_from_res(variables_at_profiles, z_at_casiers, pos_z)

                    # Compute water depth: H = Z - Zf and clip below 0m (avoid negative values)
                    depth = np.clip(values_res[pos_z, :] - z_bottom, a_min=0.0, a_max=None)

                    # Merge and write values
                    if "Vact" in varnames_1d:
                        # Compute velocity magnitude from Vact and apply mask "is active bed"
                        velocity = values_res[varnames_1d.index("Vact"), :] * values_geom[1, :]
                        values = np.vstack((values_geom, depth, velocity, values_res))
                    else:
                        values = np.vstack((values_geom, depth, values_res))

                    resout.write_entire_frame(output_header, 3600.0 * i, values)

            else:
                calc_unsteady_result = resultats.get_res_calc_trans(calc_unsteady)
                logger.info("Calcul transitoire %s" % calc_unsteady)
                res_unsteady = resultats.get_data_trans(calc_unsteady)

                for i, (time, _) in enumerate(calc_unsteady_result.frame_list):
                    logger.info("~> %fs" % time)
                    res_at_sections = res_unsteady["Section"][i, :, :]
                    variables_at_profiles = res_at_sections[pos_sections_list, :][:, pos_variables]
                    if global_mesh_constr.has_floodplain:
                        z_at_casiers = res_unsteady["Casier"][i, pos_casiers_list, pos_z_fp]
                    else:
                        z_at_casiers = None

                    # Interpolate between sections
                    values_res = global_mesh_constr.interp_values_from_res(variables_at_profiles, z_at_casiers, pos_z)

                    # Compute water depth: H = Z - Zf and clip below 0m (avoid negative values)
                    depth = np.clip(values_res[pos_z, :] - z_bottom, a_min=0.0, a_max=None)

                    # Merge and write values
                    if "Vact" in varnames_1d:
                        # Compute velocity magnitude from Vact and apply mask "is active bed"
                        velocity = values_res[varnames_1d.index("Vact"), :] * values_geom[1, :]
                        values = np.vstack((values_geom, depth, velocity, values_res))
                    else:
                        values = np.vstack((values_geom, depth, values_res))

                    resout.write_entire_frame(output_header, time, values)

    else:
        # Write a single frame with only variables from geometry
        global_mesh_constr.export_mesh(outfile_mesh, lang=lang)
