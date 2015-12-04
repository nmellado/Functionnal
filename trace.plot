set terminal png nocrop enhanced size 1000,1000 font "arial,8"
set output 'curves.png'

# Functions (1/0 means not defined)
a = 0.9
f(x) = abs(x)<2*pi ? a*sin(x)           : 1/0
g(x) = abs(x)<2*pi ? a*sin(x+pi/2)      : 1/0
h(x) = abs(x)<2*pi ? a*sin(x+pi)        : 1/0
k(x) = abs(x)<2*pi ? a*sin(x+3.0/2*pi)  : 1/0

set xrange [0:1]
set yrange [0:1]

### Start multiplot (2x2 layout)
set multiplot layout 2,2 rowsfirst
# --- GRAPH a
plot "./control0.txt" with lines ls 1 t "Control Polygon 0", \
     "./curves.txt" u 1:2 with lines ls 2  t "Bezier curve 0"
# --- GRAPH b
plot "./control1.txt" with lines ls 1 t "Control Polygon 1 ", \
     "./curves.txt" u 3:4 with lines ls 2  t "Bezier curve 1 "

set xrange [-1:1]
set yrange [-1:1]
# --- GRAPH c
plot "./hcontrol0.txt" with lines ls 1 t "Control Hodograph 0 ", \
     "./hodographs.txt" u 1:2 with lines ls 2  t "Hodograph 0 "
 # --- GRAPH d
 plot "./hcontrol1.txt" with lines ls 1 t "Control Hodograph 1 ", \
      "./hodographs.txt" u 3:4 with lines ls 2  t "Hodograph 1 "
unset multiplot
### End multiplot
