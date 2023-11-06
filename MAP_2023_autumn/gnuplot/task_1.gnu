count = 990 
step = 10
set term gif animate delay 10 size 2000, 1000
set yrange [0:1.5]

set output "data/task1.gif"

do for [n = 0 : count : step] {
    filename = sprintf("data/task1/out_%01d.dat", n)
    plot filename with lines title ""
}