# Gnuplot script file for plotting data in file "force.dat"
      # This file is called   plot-type.p
      set   autoscale                        # scale axes automatically
      unset log                              # remove any log-scaling
      unset label                            # remove any previous labels
      set xtic auto                          # set xtics automatically
      set ytic auto                          # set ytics automatically
      set title "TEST WALL WITH MIRRORS, CASE n = 1.5"
      set xlabel "Time (hr)"
      set ylabel "Power (W)"
      set key bottom
      set style line 1 lt 1 lw 2
      set style line 2 lt 2 lw 2
      set style line 3 lt 3 lw 2
      set style line 4 lt 4 lw 2 
      set style line 5 lt 5 lw 2
      set style line 6 lt 6 lw 2
      set key spacing 1.5
      #set key 0.01,100
      #set label "Yield Point" at 0.003,260
      #set arrow from -3.55,0 to -3.55,50 nohead  
      #set arrow from -3.04,0 to -3.04,50 nohead  
      set xr [10:18]
      #set yr [3e+08:4e+08]
      #set object 1 rect from -3.55,0 to -3.04,50 
      plot    "difference.txt" using 1:2 title 'Difference w-wo Mirrors WITH FIXED CODE' with lines linestyle 1, \
"expected.txt" using 1:2 title 'Approx. expected mirrors contribution' with lines linestyle 2, \
"difference-beforefixing.txt" using 1:2 title 'Difference w-wo Mirrors BEFORE FIXING CODE' with lines linestyle 3 
