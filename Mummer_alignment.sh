nucmer --prefix prefix1 -t 160 -g 1000 -l 40 -c 90 ref1.fa ref2.fa
delta-filter -1 -q -r prefix1.delta > prefix1.filter
show-coords -rclT prefix1.filter > prefix1.coords
