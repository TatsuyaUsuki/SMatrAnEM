#!/bin/sh
cd ./NonUniform  && sh ./comp_clang.sh && sh ./comp.sh && cd ../
cd ./Curvature  && sh ./comp_clang.sh && sh ./comp.sh && cd ../
cd ./CoordTrans  && sh ./comp_clang.sh && sh ./comp.sh && cd ../
cd ./OuterRegion  && sh ./comp_clang.sh && sh ./comp.sh && cd ../
cd ./Structure  && sh ./comp_clang.sh && sh ./comp.sh && cd ../
cd ./Medium3D  && sh ./comp_clang.sh && sh ./comp.sh && cd ../
cd ./Medium3D2a  && sh ./comp_clang.sh && sh ./comp.sh && cd ../
cd ./MixedData  && sh ./comp_clang.sh && sh ./comp.sh && cd ../
cd ./CorrFactor  && sh ./comp_clang.sh && sh ./comp.sh && cd ../
cd ./StrSlice  && sh ./comp_clang.sh && sh ./comp.sh && cd ../
cd ./ModeCalc  && sh ./comp_clang.sh && sh ./comp.sh && cd ../
cd ./ModeCalc2  && sh ./comp_clang.sh && sh ./comp.sh && cd ../
cd ./PMLset  && sh ./comp_clang.sh && sh ./comp.sh && cd ../
cd ./Smat  && sh ./comp_clang.sh && sh ./comp.sh && cd ../
cd ./FDTDcalc  && sh ./comp_clang.sh && sh ./comp.sh && cd ../

