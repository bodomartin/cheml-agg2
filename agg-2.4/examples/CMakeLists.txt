
#+bm we need this because the demos look for the includes via "...", not < ... >
include_directories( ${antigrain_SOURCE_DIR}/include )

ADD_EXECUTABLE( aa_demo ${WIN32GUI}
    aa_demo.cpp
    )
target_link_libraries(aa_demo agg_platform agg_controls agg)
  

ADD_EXECUTABLE( aa_test ${WIN32GUI}
    aa_test.cpp
)
target_link_libraries(aa_test agg_platform agg_controls agg)

ADD_EXECUTABLE( alpha_gradient ${WIN32GUI}
    alpha_gradient.cpp
)
target_link_libraries(alpha_gradient agg_platform agg_controls agg)

ADD_EXECUTABLE( alpha_mask ${WIN32GUI}
    alpha_mask.cpp
    parse_lion.cpp
)
target_link_libraries(alpha_mask agg_platform agg_controls agg)

ADD_EXECUTABLE( alpha_mask2 ${WIN32GUI}
    alpha_mask.cpp
    parse_lion.cpp
)
target_link_libraries(alpha_mask2 agg_platform agg_controls agg)

ADD_EXECUTABLE( alpha_mask3 ${WIN32GUI}
    alpha_mask.cpp
    make_arrows.cpp
    make_gb_poly.cpp
    parse_lion.cpp
)
target_link_libraries(alpha_mask3 agg_platform agg_controls agg)

ADD_EXECUTABLE( bezier_div ${WIN32GUI}
    bezier_div.cpp
    interactive_polygon.cpp
)
target_link_libraries(bezier_div agg_platform agg_controls agg)

ADD_EXECUTABLE( blend_color ${WIN32GUI}
    blend_color.cpp
)
target_link_libraries(blend_color agg_platform agg_controls agg)

ADD_EXECUTABLE( blur ${WIN32GUI}
    blur.cpp
)
target_link_libraries(blur agg_platform agg_controls agg)

ADD_EXECUTABLE( bspline ${WIN32GUI}
    bspline.cpp
    interactive_polygon.cpp
)
target_link_libraries(bspline agg_platform agg_controls agg)

ADD_EXECUTABLE( circles ${WIN32GUI}
    circles.cpp
)
target_link_libraries(circles agg_platform agg_controls agg)

ADD_EXECUTABLE( component_rendering ${WIN32GUI}
    component_rendering.cpp
)
target_link_libraries(component_rendering agg_platform agg_controls agg)

ADD_EXECUTABLE( compositing ${WIN32GUI}
    compositing.cpp
)
target_link_libraries(compositing agg_platform agg_controls agg)

ADD_EXECUTABLE( compositing2 ${WIN32GUI}
    compositing2.cpp
)
target_link_libraries(compositing2 agg_platform agg_controls agg)

ADD_EXECUTABLE( conv_contour ${WIN32GUI}
    conv_contour.cpp
)
target_link_libraries(conv_contour agg_platform agg_controls agg)

ADD_EXECUTABLE( conv_dash_marker ${WIN32GUI}
    conv_dash_marker.cpp
)
target_link_libraries(conv_dash_marker agg_platform agg_controls agg)

ADD_EXECUTABLE( conv_stroke ${WIN32GUI}
    conv_stroke.cpp
)
target_link_libraries(conv_stroke agg_platform agg_controls agg)

ADD_EXECUTABLE( distortions ${WIN32GUI}
    distortions.cpp
)
target_link_libraries(distortions agg_platform agg_controls agg)

ADD_EXECUTABLE( flash_rasterizer ${WIN32GUI}
    flash_rasterizer.cpp
)
target_link_libraries(flash_rasterizer agg_platform agg_controls agg)

ADD_EXECUTABLE( flash_rasterizer2 ${WIN32GUI}
    flash_rasterizer2.cpp
)
target_link_libraries(flash_rasterizer2 agg_platform agg_controls agg)

IF ( agg_USE_FREETYPE )

  ADD_EXECUTABLE( freetype_test ${WIN32GUI}
    freetype_test.cpp
    make_arrows.cpp
    make_gb_poly.cpp
    )
  target_link_libraries(treetype_test agg_platform agg_controls agg freetype)
      
ENDIF()

ADD_EXECUTABLE( gamma_correction ${WIN32GUI}
    gamma_correction.cpp
)
target_link_libraries(gamma_correction agg_platform agg_controls agg)

ADD_EXECUTABLE( gamma_ctrl ${WIN32GUI}
    gamma_ctrl.cpp
)
target_link_libraries(gamma_ctrl agg_platform agg_controls agg)

ADD_EXECUTABLE( gamma_tuner ${WIN32GUI}
    gamma_tuner.cpp
)
target_link_libraries(gamma_tuner agg_platform agg_controls agg)

ADD_EXECUTABLE( gouraud ${WIN32GUI}
    gouraud.cpp
)
target_link_libraries(gouraud agg_platform agg_controls agg)

ADD_EXECUTABLE( gouraud_mesh ${WIN32GUI}
    gouraud_mesh.cpp
)
target_link_libraries(gouraud_mesh agg_platform agg_controls agg)

IF ( agg_USE_GPC )
    ADD_EXECUTABLE( gpc_test ${WIN32GUI}
        gpc_test.cpp
        make_arrows.cpp
        make_gb_poly.cpp
        )
      
      target_link_libraries(gourand_mesh agg_platform agg_controls gpc)
      
ENDIF()

ADD_EXECUTABLE( gradients ${WIN32GUI}
    gradients.cpp
)
target_link_libraries(gradients agg_platform agg_controls agg)

ADD_EXECUTABLE( gradient_focal ${WIN32GUI}
    gradient_focal.cpp
)
target_link_libraries(gradient_focal agg_platform agg_controls agg)

ADD_EXECUTABLE( gradients_contour ${WIN32GUI}
    gradients_contour.cpp
    make_arrows.cpp
    make_gb_poly.cpp
    parse_lion.cpp
)
target_link_libraries(gradients_contour agg_platform agg_controls agg)

ADD_EXECUTABLE( graph_test ${WIN32GUI}
    graph_test.cpp
)
target_link_libraries(graph_test agg_platform agg_controls agg)

ADD_EXECUTABLE( idea ${WIN32GUI}
    idea.cpp
)
target_link_libraries(idea agg_platform agg_controls agg)

ADD_EXECUTABLE( image1 ${WIN32GUI}
    image1.cpp
)
target_link_libraries(image1 agg_platform agg_controls agg)

ADD_EXECUTABLE( image_alpha ${WIN32GUI}
    image_alpha.cpp
)
target_link_libraries(image_alpha agg_platform agg_controls agg)

ADD_EXECUTABLE( image_filters ${WIN32GUI}
    image_filters.cpp
)
target_link_libraries(image_filters agg_platform agg_controls agg)

ADD_EXECUTABLE( image_filters2 ${WIN32GUI}
    image_filters2.cpp
)
target_link_libraries(image_filters2 agg_platform agg_controls agg)

ADD_EXECUTABLE( image_fltr_graph ${WIN32GUI}
    image_fltr_graph.cpp
)
target_link_libraries(image_fltr_graph agg_platform agg_controls agg)

ADD_EXECUTABLE( image_perspective ${WIN32GUI}
    image_perspective.cpp
    interactive_polygon.cpp
)
target_link_libraries(image_perspective agg_platform agg_controls agg)

ADD_EXECUTABLE( image_resample ${WIN32GUI}
    image_resample.cpp
    interactive_polygon.cpp
)
target_link_libraries(image_resample agg_platform agg_controls agg)

ADD_EXECUTABLE( image_transforms ${WIN32GUI}
    image_transforms.cpp
)
target_link_libraries(image_transforms agg_platform agg_controls agg)

ADD_EXECUTABLE( line_patterns ${WIN32GUI}
    line_patterns.cpp
)
target_link_libraries(line_patterns agg_platform agg_controls agg)

ADD_EXECUTABLE( line_patterns_clip ${WIN32GUI}
    line_patterns_clip.cpp
)
target_link_libraries(line_patterns_clip agg_platform agg_controls agg)

ADD_EXECUTABLE( lion ${WIN32GUI}
    lion.cpp
    parse_lion.cpp
)
target_link_libraries(lion agg_platform agg_controls agg)

ADD_EXECUTABLE( lion_lens ${WIN32GUI}
    lion_lens.cpp
    parse_lion.cpp
)
target_link_libraries(lion_lens agg_platform agg_controls agg)

ADD_EXECUTABLE( lion_outline ${WIN32GUI}
    lion_outline.cpp
    parse_lion.cpp
)
target_link_libraries(lion_outline agg_platform agg_controls agg)

ADD_EXECUTABLE( mol_view ${WIN32GUI}
    mol_view.cpp
)
target_link_libraries(mol_view agg_platform agg_controls agg)

# bm 2016-05-22 expanded example (3d rotation, picking ...)
# to be used for cheml project 
ADD_EXECUTABLE( mol_view2 ${WIN32GUI}
    mol_view2/mol_view.cpp
    mol_view2/sv2
    mol_view2/sv_math
    mol_view2/molecule
    
    )
target_link_libraries(mol_view2 agg_platform agg_controls agg)
add_custom_target( mol_view2_sdfiles
  COMMAND ${CMAKE_COMMAND} -E copy ${antigrain_SOURCE_DIR}/examples/mol_view2/ketofluorene2_v2.sdf ${PROJECT_BINARY_DIR}/examples
  COMMAND ${CMAKE_COMMAND} -E copy ${antigrain_SOURCE_DIR}/examples/mol_view2/HEM_ideal.sdf ${PROJECT_BINARY_DIR}/examples
  )
add_dependencies(mol_view2 mol_view2_sdfiles)

ADD_EXECUTABLE( multi_clip ${WIN32GUI}
    multi_clip.cpp
    parse_lion.cpp
)
target_link_libraries(multi_clip agg_platform agg_controls agg)

ADD_EXECUTABLE( pattern_fill ${WIN32GUI}
    pattern_fill.cpp
)
target_link_libraries(pattern_fill agg_platform agg_controls agg)

ADD_EXECUTABLE( pattern_perspective ${WIN32GUI}
    pattern_perspective.cpp
    interactive_polygon.cpp
)
target_link_libraries(pattern_perspective agg_platform agg_controls agg)

ADD_EXECUTABLE( pattern_resample ${WIN32GUI}
    pattern_resample.cpp
    interactive_polygon.cpp
)
target_link_libraries(pattern_resample agg_platform agg_controls agg)

ADD_EXECUTABLE( perspective ${WIN32GUI}
    perspective.cpp
    interactive_polygon.cpp
    parse_lion.cpp
)
target_link_libraries(perspective agg_platform agg_controls agg)

ADD_EXECUTABLE( polymorphic_renderer ${WIN32GUI}
    polymorphic_renderer.cpp
)
target_link_libraries(polymorphic_renderer agg_platform agg_controls agg)

ADD_EXECUTABLE( rasterizers ${WIN32GUI}
    rasterizers.cpp
)
target_link_libraries(rasterizers agg_platform agg_controls agg)

ADD_EXECUTABLE( rasterizers2 ${WIN32GUI}
    rasterizers2.cpp
)
target_link_libraries(rasterizers2 agg_platform agg_controls agg)

ADD_EXECUTABLE( rasterizer_compound ${WIN32GUI}
    rasterizer_compound.cpp
)
target_link_libraries(rasterizer_compound agg_platform agg_controls agg)

ADD_EXECUTABLE( raster_text ${WIN32GUI}
    raster_text.cpp
)
target_link_libraries(raster_text agg_platform agg_controls agg)

ADD_EXECUTABLE( rounded_rect ${WIN32GUI}
    rounded_rect.cpp
)
target_link_libraries(rounded_rect agg_platform agg_controls agg)

ADD_EXECUTABLE( scanline_boolean ${WIN32GUI}
    scanline_boolean.cpp
    interactive_polygon.cpp
)
target_link_libraries(scanline_boolean agg_platform agg_controls agg)

ADD_EXECUTABLE( scanline_boolean2 ${WIN32GUI}
    scanline_boolean2.cpp
    make_arrows.cpp
    make_gb_poly.cpp
)
target_link_libraries(scanline_boolean2 agg_platform agg_controls agg)

ADD_EXECUTABLE( simple_blur ${WIN32GUI}
    simple_blur.cpp
    parse_lion.cpp
)
target_link_libraries(simple_blur agg_platform agg_controls agg)

IF(WIN32)

  ADD_EXECUTABLE( trans_curve1 ${WIN32GUI}
    trans_curve1.cpp
    interactive_polygon.cpp
    )
  target_link_libraries(trans_curve1 agg_platform agg_controls agg)

    ADD_EXECUTABLE( trans_curve2 ${WIN32GUI}
        trans_curve2.cpp
        interactive_polygon.cpp
    )
  target_link_libraries(trans_curve2 agg_platform agg_controls agg)

    ADD_EXECUTABLE( truetype_test ${WIN32GUI}
        truetype_test.cpp
    )
  target_link_libraries(truetype_test agg_platform agg_controls agg)

ENDIF(WIN32)

ADD_EXECUTABLE( trans_polar ${WIN32GUI}
    trans_polar.cpp
)
target_link_libraries(trans_polar agg_platform agg_controls agg)

IF ( agg_USE_EXPAT )

  ADD_EXECUTABLE( svg_test ${WIN32GUI}
        ./svg_viewer/svg_test.cpp
        ./svg_viewer/agg_svg_exception.h
        ./svg_viewer/agg_svg_parser.cpp
        ./svg_viewer/agg_svg_parser.h
        ./svg_viewer/agg_svg_path_renderer.cpp
        ./svg_viewer/agg_svg_path_renderer.h
        ./svg_viewer/agg_svg_path_tokenizer.cpp
        ./svg_viewer/agg_svg_path_tokenizer.h
        )

      target_link_libraries(svg_test agg_platform agg_controls agg)

ENDIF()

IF(WIN32)
    ADD_EXECUTABLE( pure_api ${WIN32GUI}
        ./win32_api/pure_api/pure_api.h
        ./win32_api/pure_api/pure_api.cpp
        ./win32_api/pure_api/resource.h
        ./win32_api/pure_api/StdAfx.h
        ./win32_api/pure_api/StdAfx.cpp
        ./win32_api/pure_api/pure_api.rc
        parse_lion.cpp
        
        )

      target_link_libraries(pure_api agg_platform agg_controls agg)
     
ENDIF()

IF( agg_USE_AGG2D )
    ADD_EXECUTABLE( agg2_demo ${WIN32GUI}
        agg2d_demo.cpp
    )

  target_link_libraries(agg2d_demo agg_platform agg_controls agg)

ENDIF( agg_USE_AGG2D )
    

