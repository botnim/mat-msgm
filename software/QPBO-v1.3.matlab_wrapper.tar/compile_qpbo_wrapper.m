function compile_qpbo_wrapper()
%
% Compiling wrapper for QPBO
%



mex -O -largeArrayDims CXXFLAGS="\$CXXFLAGS -Wno-write-strings"...
    QPBO.cpp QPBO_extra.cpp QPBO_maxflow.cpp QPBO_postprocessing.cpp ...
    QPBO_wrapper_mex.cpp -output QPBO_wrapper_mex
