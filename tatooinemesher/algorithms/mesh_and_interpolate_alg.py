from tatooinemesher.constraint_line import ConstraintLine
from tatooinemesher.mesh_constructor import MeshConstructor
from tatooinemesher.section import CrossSectionSequence
from tatooinemesher.utils import TatooineException, get_hydraulic_axis, logger, set_logger_level


def mesh_and_interpolate(
    infile_axis,
    infile_cross_sections,
    attr_cross_sections,
    long_step,
    infile_constraint_lines=None,
    infile_constraint_3D_lines=None,
    z_line_strength=0.0,
    z_line_gap_scale=0.0,
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
    if infile_constraint_lines is None and infile_constraint_3D_lines is None:
        constraint_lines = ConstraintLine.get_lines_and_set_limits_from_sections(section_seq, interp_constraint_lines)
    else:
        constraint_lines = []
        if infile_constraint_lines is not None:
            constraint_lines.extend(
                ConstraintLine.get_lines_from_file(infile_constraint_lines, interp_constraint_lines)
            )
        if infile_constraint_3D_lines is not None:
            constraint_lines_3d = ConstraintLine.get_lines_from_file(
                infile_constraint_3D_lines, interp_constraint_lines, has_z=True
            )
            constraint_lines.extend(constraint_lines_3d)
            logger.info(
                "~> %d Z-aware constraint line(s) loaded for longitudinal Z correction" % len(constraint_lines_3d)
            )
        # Renumber ids contiguously across the merged list
        for i, line in enumerate(constraint_lines):
            line.id = i
        if nb_pts_lat is not None and len(constraint_lines) != 2:
            raise TatooineException("Argument `--nb_pts_lat` is only compatible with 2 constraint lines!")
        if interp_values.startswith("BI"):
            if len(constraint_lines) != 2:
                raise TatooineException("A 2D interpolation is only compatible with 2 constraint lines!")
            if any(line.has_z for line in constraint_lines):
                logger.warning("2D interpolation mode (%s): Z profiles on constraint lines are ignored" % interp_values)
        section_seq.find_and_add_limits(constraint_lines, dist_max)

    mesh_constr = MeshConstructor(
        section_seq=section_seq,
        lat_step=lat_step,
        nb_pts_lat=nb_pts_lat,
        interp_values=interp_values,
        z_line_strength=z_line_strength,
        z_line_gap_scale=z_line_gap_scale,
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
