cat /dev/null > psolver_benchmarking_result.csv
TIMEFORMAT='%S'
for size in 4 8 16 32 64 128 256 512 1024
do
    for iter in 1 2 3 4 5 6 7 8
    do
        cp ../inputs/benchmarking.inp ./benchmarking_$size.inp
        sed -i "s/n_x = 256/n_x = $size/g" ./benchmarking_$size.inp
        sed -i "s/n_y = 256/n_y = $size/g" ./benchmarking_$size.inp
        echo $size >> psolver_benchmarking_result.csv
        time ( ./../build/GrossPitaevskiiFDM_run ./benchmarking_$size.inp ) 2>&1 1>/dev/null 2>> psolver_benchmarking_result.csv
    done
done
cat psolver_benchmarking_result.csv > sed 'N;s/\n/, /' > sed 'N;s/\n/, /' > sed 'N;s/\n/, /' > tmp.txt
cat tmp.txt > psolver_benchmarking_result.csv
rm tmp.txt
