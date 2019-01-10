#!/bin/bash

for i in *.csv; do
 Name=`basename "$i" .csv`
 echo $Name
 ./sect <$i | grep -v "NaN" | sed -e "s/Normalized\ heat\ flow/Normalized_heat_flow/"> $Name.dat
# ./start-R.sh Deconvolution-Kinetic-Sequential.R $Name
done

#./start-R.sh Deconvolution-Kinetic-Sequential.R 2-NO2C6H4N2+_TfO-_75_Nitrogen_12-20-16  0.0033424 22.0 -414.0 0.0026 0.00026 0.001128 0.0001128
./start-R.sh Deconvolution-Kinetic-Sequential.R 2-NO2C6H4N2+_TfO-_75_Nitrogen_12-20-16  0.0033424 2.209046e+01 -4.580509e+02 9.838210e-01 7.673343e-02 7.263796e-04 1.580044e-03

#./start-R.sh Deconvolution-Kinetic-Sequential.R 2-NO2C6H4N2+_TfO-_80_Nitrogen_8-22-16   0.0033424 22.0 -386.0 0.0055 0.00055 0.001243 0.0001243
./start-R.sh Deconvolution-Kinetic-Sequential.R 2-NO2C6H4N2+_TfO-_80_Nitrogen_8-22-16   0.0033424 2.910810e+01 -4.269288e+02 5.796413e-01 3.670586e-02 1.155185e-03 3.207568e-03

#./start-R.sh Deconvolution-Kinetic-Sequential.R 2-NO2C6H4N2+_TfO-_85_Nitrogen_7-1-16    0.0033424 22.0 -396.0 0.0076 0.00076 0.001461 0.0001461
./start-R.sh Deconvolution-Kinetic-Sequential.R 2-NO2C6H4N2+_TfO-_85_Nitrogen_7-1-16    0.0033424 2.486793e+01 -4.133467e+02 9.934544e-01 5.962921e-02 7.021213e-04 2.755259e-03


#./start-R.sh Deconvolution-Kinetic-Sequential.R 4-CH3OC6H4N2+_TfO-_75_Nitrogen_12-25-15 0.0035185 29.1 -183.1 0.0210 0.00210 0.000109 0.0000109
./start-R.sh Deconvolution-Kinetic-Sequential.R 4-CH3OC6H4N2+_TfO-_75_Nitrogen_12-25-15 0.0035185 -9.870400e+01 -1.946068e+02 2.527550e+01 1.107108e-01 6.573825e-04 4.039706e-04

#./start-R.sh Deconvolution-Kinetic-Sequential.R 4-CH3OC6H4N2+_TfO-_80_Nitrogen_2-15-16  0.0035185 29.1 -183.2 0.0440 0.00440 0.000113 0.0000113
./start-R.sh Deconvolution-Kinetic-Sequential.R 4-CH3OC6H4N2+_TfO-_80_Nitrogen_2-15-16  0.0035185 2.964153e+01 -2.133650e+02 6.223540e-02 4.383461e-01 2.797659e+01 1.115711e-04

#./start-R.sh Deconvolution-Kinetic-Sequential.R 4-CH3OC6H4N2+_TfO-_85_Nitrogen_1-12-16  0.0035185 29.1 -106.0 0.1280 0.01280 0.000238 0.0000238
./start-R.sh Deconvolution-Kinetic-Sequential.R 4-CH3OC6H4N2+_TfO-_85_Nitrogen_1-12-16  0.0035185 1.312668e+01 -1.194273e+02 3.767194e-02 7.390249e-01 2.905274e+01 2.392651e-04
