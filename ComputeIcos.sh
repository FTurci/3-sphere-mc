for (( i = 10000; i < 1000000; i+=20000 )); do
    python sphere_voronoi.py $i
    ./icos 120 bonds.txt 
done

# grep 'b)' logs | sed 's/[^0-9]//g' 