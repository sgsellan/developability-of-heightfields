% Build mex functions required to run our ADMM optimization. This assumes you have
% already built gptoolbox so that libigl has been cloned into gptoolbox/mex/external/


mex('argmin_X_full_mex.cpp','-I../gptoolbox/mex/external/libigl/include/','-I../gptoolbox/mex/external/libigl/external/eigen/','-V');
mex('argmin_X_full_L2_mex.cpp','-I../gptoolbox/mex/external/libigl/include/','-I../gptoolbox/mex/external/libigl/external/eigen/','-V');
