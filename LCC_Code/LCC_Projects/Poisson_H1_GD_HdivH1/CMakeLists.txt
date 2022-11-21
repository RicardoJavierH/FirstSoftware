file(GLOB header *.h)
file(GLOB sources *.cpp)

set(main_SRCS ${sources} ${headers} ${LCCCommon_SRC} ${LCCMaterial_SRC})

add_executable(Simulator_H1_GD_HdivH1 ${main_SRCS})

target_link_libraries(Simulator_H1_GD_HdivH1 ${PTHREAD_LIB})
target_link_libraries(Simulator_H1_GD_HdivH1 ${PZ_LIB})

