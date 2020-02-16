#!/bin/sh
#------- 1 -------
./xi2u ./input.txt && \
#------- 2 -------
./mkCurv ./input.txt && \
#------- 3 -------
./u2r2x ./input.txt 2> stderr_urx.dat && \
#------- 4 -------
./outer ./input.txt && \
#------- 5 -------
./StrMap ./input.txt 2> stderr_Str.dat && \
#------- 6 -------
./Med3D ./input.txt 2> stderr_Med3D.dat && \
#------- 7 -------
./BindData ./input.txt 2> stderr_BindData.dat && \
#------- 8 -------
./CorrF ./input.txt 2> stderr_CorrF.dat && \
#------- 9 -------
./YeeSlice ./input.txt 2> stderr_YeeSlice.dat && \
#------ 10 -------
./Mcalc ./input.txt 2> stderr_Mcalc.dat && \
#------ 11 -------
./PMLgen ./input.txt 2> stderr_PMLgen.dat && \
#------ 12 -------
./Scalc ./input.txt 2> stderr_Scalc.dat && \
#------ 13 -------
./FDTD ./input.txt 2> stderr_FDTD.dat #&& \

