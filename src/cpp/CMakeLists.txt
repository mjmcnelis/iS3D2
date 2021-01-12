set (SOURCES
    AnisoVariables.cpp
    Arsenal.cpp
    BinSampledParticle.cpp
    DeltafData.cpp
    EmissionFunction.cpp
    GaussThermal.cpp
    iS3D.cpp
    LocalRestFrame.cpp
    Main.cpp
    Momentum.cpp
    MomentumSpectra.cpp
    ParameterReader.cpp
    ParticleSampler.cpp
    Polarization.cpp
    readindata.cpp
    SpacetimeDistribution.cpp
    Table.cpp
    )

add_library(iS3D_lib SHARED ${SOURCES})
add_executable(iS3D ${SOURCES})


target_link_libraries (iS3D_lib ${GSL_LIBRARIES})
target_link_libraries (iS3D ${GSL_LIBRARIES})

install(TARGETS iS3D DESTINATION ${CMAKE_HOME_DIRECTORY})
