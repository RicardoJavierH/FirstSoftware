# common source files used to analysis classes 
#

file(GLOB all_sources *.cpp *.h)

set(LCCAnalysis_SRC ${all_sources} ${LCCMaterial_SRC} ${LCC_Estimator_SRC} PARENT_SCOPE)
   
include_directories(.)