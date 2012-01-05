#pushd ~/scratch0/CMS/CMSSW_1_6_11/src
pushd /soft/lip-sw/cmssw/users/anayak/CMS/NewFWK/CMSSW_4_2_8/src/
setcms
eval `scramv1 runtime -csh`
popd
setenv MY_PATH ${PWD}/MiniTree/Selection/src/

if( $?LD_LIBRARY_PATH ) then
  setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${MY_PATH}
else
  setenv LD_LIBRARY_PATH ${MY_PATH}
endif

