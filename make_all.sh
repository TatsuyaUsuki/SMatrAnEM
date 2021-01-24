#!/bin/sh
cd ./NonUniform  && make && cd ../
cd ./Curvature  && make && cd ../
cd ./CoordTrans  && make && cd ../
cd ./OuterRegion  && make && cd ../
cd ./Structure  && make && cd ../
cd ./Medium3D  && make && cd ../
cd ./Medium3D2a  && make && cd ../
cd ./Medium3DSi  && make && cd ../
cd ./MixedData  && make && cd ../
cd ./CorrFactor  && make && cd ../
cd ./StrSlice  && make && cd ../
cd ./ModeCalc  && make && cd ../
cd ./ModeCalc2  && make && cd ../
cd ./ModeSlice  && make && cd ../
cd ./PMLset  && make && cd ../
cd ./Smat  && make && cd ../
cd ./FDTDcalc  && make && cd ../

