#!/bin/bash

for i in *.csv; do
 Name=`basename "$i" .csv`
 echo $Name
# ./sect <$i | grep -v "NaN" | sed -e "s/Normalized\ heat\ flow/Normalized_heat_flow/"> $Name.dat
# ./start-R.sh Deconvolution-Kinetic-Concurent.R $Name
done

#./start-R.sh Deconvolution-Kinetic-Concurent.R 2-NO2C6H4N2+_TfO-_75_Nitrogen_12-20-16  0.0033424 22.0 -414.0 0.0026 0.00026 0.001128 0.0001128
./start-R.sh Deconvolution-Kinetic-Concurent.R 2-NO2C6H4N2+_TfO-_75_Nitrogen_12-20-16  0.0033424 4.987779e+01 -4.927423e+02 5.127006e-01 4.791654e-02 1.353114e-03 3.115688e-03

#./start-R.sh Deconvolution-Kinetic-Concurent.R 2-NO2C6H4N2+_TfO-_80_Nitrogen_8-22-16   0.0033424 22.0 -386.0 0.0055 0.00055 0.001243 0.0001243
##./start-R.sh Deconvolution-Kinetic-Concurent.R 2-NO2C6H4N2+_TfO-_80_Nitrogen_8-22-16   0.0033424 5.617100e+01 -5.071608e+02 5.138825e-01 3.294921e-03 1.169264e-03 4.203698e-02
./start-R.sh Deconvolution-Kinetic-Concurent.R 2-NO2C6H4N2+_TfO-_80_Nitrogen_8-22-16   0.0033424 35.17100 -450.1608 0.62 0.047 1.169264e-03 3.00e-03

#./start-R.sh Deconvolution-Kinetic-Concurent.R 2-NO2C6H4N2+_TfO-_85_Nitrogen_7-1-16    0.0033424 22.0 -396.0 0.0076 0.00076 0.001461 0.0001461
./start-R.sh Deconvolution-Kinetic-Concurent.R 2-NO2C6H4N2+_TfO-_85_Nitrogen_7-1-16    0.0033424 2.812087e+01 -4.259087e+02 9.342798e-01 6.351878e-02 7.149883e-04 2.537486e-03


#./start-R.sh Deconvolution-Kinetic-Concurent.R 4-CH3OC6H4N2+_TfO-_75_Nitrogen_12-25-15 0.0035185 29.1 -183.1 0.0210 0.00210 0.000109 0.0000109
./start-R.sh Deconvolution-Kinetic-Concurent.R 4-CH3OC6H4N2+_TfO-_75_Nitrogen_12-25-15 0.0035185 7.260436e+01 -4.374580e+02 6.255051e-01 5.758671e-02 1.933897e-06 2.064136e-03

#./start-R.sh Deconvolution-Kinetic-Concurent.R 4-CH3OC6H4N2+_TfO-_80_Nitrogen_2-15-16  0.0035185 29.1 -183.2 0.0440 0.00440 0.000113 0.0000113
./start-R.sh Deconvolution-Kinetic-Concurent.R 4-CH3OC6H4N2+_TfO-_80_Nitrogen_2-15-16  0.0035185 7.260436e+01 -4.374580e+02 3.227525e-01 3.150887e-01 1.005415e-03 1.016961e-03

#./start-R.sh Deconvolution-Kinetic-Concurent.R 4-CH3OC6H4N2+_TfO-_85_Nitrogen_1-12-16  0.0035185 29.1 -106.0 0.1280 0.01280 0.000238 0.0000238
./start-R.sh Deconvolution-Kinetic-Concurent.R 4-CH3OC6H4N2+_TfO-_85_Nitrogen_1-12-16  0.0035185 2.730354e+01 -1.380430e+02 1.655319e-01 4.406386e-01 3.662907e-02 5.621479e-04
