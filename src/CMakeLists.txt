add_executable(PunchTable)
target_include_directories(PunchTable PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

target_sources(PunchTable PRIVATE
    CubicSpline.h
    CubicSpline.cpp
    MassLookup.h
    MassLookup.cpp
    PunchTable.h
    PunchTable.cpp
    GenerateTable.h
    GenerateTable.cpp
    main.cpp
)

target_link_libraries(PunchTable catima)

set_target_properties(PunchTable PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PUNCHTABLE_BINARY_DIR})