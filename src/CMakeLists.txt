
# build target library



target_include_directories(NuXSSplMkr PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)
target_link_libraries(NuXSSplMkr 
    photospline
    ${LHAPDF_LIBRARIES} 
    ${APFEL_LIBRARIES} 
    ${GSL_LIBRARIES} 
    ${BOOST_LIBRARIES}
)

install(TARGETS NuXSSplMkr
    EXPORT NuXSSplMkrTargets
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    RUNTIME DESTINATION bin
    INCLUDES DESTINATION include
)

install(DIRECTORY ${PROJECT_SOURCE_DIR}/include/ 
    DESTINATION include
)