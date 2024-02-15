from tatooinemesher.constraint_line import ConstraintLine
from tatooinemesher.mesh_constructor import MeshConstructor
from tatooinemesher.section import CrossSection, CrossSectionSequence
from tatooinemesher.utils import get_hydraulic_axis, logger, set_logger_level, TatooineException


def mesh_and_interpolate(
    infile_axis,
    infile_cross_sections,
    attr_cross_sections,
    long_step,
    infile_constraint_lines=None,
    interp_constraint_lines="LINEAR",
    interp_values="LINEAR",
    project_straight_line=False,
    nb_pts_lat=None,
    lat_step=None,
    constant_long_disc=False,
    dist_max=0.01,
    outfile_nodes=None,
    outfile_mesh=None,
    lang="en",
    verbose=False,
):
    set_logger_level(verbose)

    logger.info("~> Reading input files")
    axe = get_hydraulic_axis(infile_axis)
    section_seq = CrossSectionSequence.from_file(
        infile_cross_sections,
        "Cross-section",
        field_id=attr_cross_sections,
        project_straight_line=project_straight_line,
    )

    # if args.infile_epis is not None and args.dist_corr_epi is not None:
    #     epis = CrossSectionSequence.from_file(args.infile_epis, "Groynes", field_id=args.attr_epis,
    #                                           project_straight_line=args.project_straight_line)
    # else:
    #     epis = None

    section_seq.compute_dist_proj_axe(axe, dist_max)
    section_seq.check_intersections()
    section_seq.sort_by_dist()
    # section_seq.export_sections_shp('profiles_projected.shp')  # DEBUG

    # TODO : Add option to create this line even if CL provided
    if infile_constraint_lines is None:
        constraint_lines = ConstraintLine.get_lines_and_set_limits_from_sections(section_seq, interp_constraint_lines)
    else:
        constraint_lines = ConstraintLine.get_lines_from_file(infile_constraint_lines, interp_constraint_lines)
        if nb_pts_lat is not None and len(constraint_lines) != 2:
            raise TatooineException("Argument `--nb_pts_lat` is only compatible with 2 constraint lines!")
        if interp_values.startswith("BI") and len(constraint_lines) != 2:
            raise TatooineException("A 2D interpolation is only compatible with 2 constraint lines!")
        section_seq.find_and_add_limits(constraint_lines, dist_max)

    mesh_constr = MeshConstructor(
        section_seq=section_seq, lat_step=lat_step, nb_pts_lat=nb_pts_lat, interp_values=interp_values
    )
    logger.info("~> Building interpolation")
    mesh_constr.build_interp(constraint_lines, long_step, constant_long_disc)
    # mesh_constr.export_segments('check_segments.shp')  # DEBUG

    # if epis is not None:
    #     mesh_constr.corr_bathy_on_epis(epis, args.dist_corr_epi)

    logger.info("~> Export outputs")
    if outfile_nodes is not None:
        mesh_constr.export_points(outfile_nodes)

    if outfile_mesh is not None:
        mesh_constr.build_mesh()
        mesh_constr.export_mesh(outfile_mesh, lang=lang)
