set terminal postscript eps mono
set key inside right top vertical Right noreverse enhanced autotitles box linetype -1 linewidth 0.200
#set title "Deconvolution for 85 C" 
set ylabel "Heat flow (P), mW/g" font "Helvetica-Bold,26"
set xlabel "Time, h" font "Helvetica-Bold,26"
set bars small
#set xrange [0:0.1]
#set yrange [-100:0]
#set size 0.5,0.5
#set terminal postscript enhanced "Courier" 20

set terminal svg size 1200,900 font "Helvetica,26"
set key autotitle columnhead
set datafile separator ","
set termoption dash

set linestyle 1 lt 1 lw 2 lc 0
set linestyle 2 lt 2 lw 2 lc 1
set linestyle 3 lt 3 lw 2 lc 2
set linestyle 4 lt 4 lw 2 lc 3
set linestyle 5 lt 5 lw 2 lc 4
set linestyle 6 lt 6 lw 2 lc 5
set linestyle 7 lt 7 lw 2 lc 6
set linestyle 8 lt 8 lw 2 lc 7

set output "test-deconvolution-85.svg"
plot "test-deconvolution-85.csv" using (($1)/3600):(1000*($8)) with lines linestyle 4 ti "Hf1", \
 "test-deconvolution-85.csv" using (($1)/3600):(1000*($9)) with lines linestyle 3 ti "Hf2", \
 "test-deconvolution-85.csv" using (($1)/3600):(1000*($10)) with lines linestyle 1 ti "Hf", \
 "o-NO2-Ph-N2-OTf-85.dat" using (($1)/3600):(1000*($5)) with lines linestyle 2 ti "Hf-exp"

set output "test-deconvolution-85-black.svg"
set style line 1 lt 1 lw 3 pt 5 linecolor rgb "black"
set style line 2 lt 2 lw 3 pt 6 linecolor rgb "black"
set style line 3 lt 3 lw 3 pt 7 linecolor rgb "gray"
set style line 4 lt 14 lw 3 pt 8 linecolor rgb "gray"
set label "P_1" at 50,0.0
set label "P_2" at 2,2.2
set label "P_{Sum}" at 50,1.7
set label "P_{Exp}" at 260,0.75
set nokey
#set key samplen -1 spacing 1.3 font ",26" box lw 0
plot "test-deconvolution-85.csv" using (($1)/3600):(1000*($8)) with lines linestyle 4 ti "P_1", \
 "test-deconvolution-85.csv" using (($1)/3600):(1000*($9)) with lines linestyle 3 ti "P_2", \
 "test-deconvolution-85.csv" using (($1)/3600):(1000*($10)) with lines linestyle 1 ti "P_{Sum}", \
 "o-NO2-Ph-N2-OTf-85.dat" using (($1)/3600):(1000*($5)) with lines linestyle 2 ti "P_{Exp}"



quit
