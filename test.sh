export SIZE=1000
./serial -n $SIZE -o serial.txt
./openmp -n $SIZE -o openmp.txt
./pthreads -n $SIZE -o pthreads.txt
diff serial.txt openmp.txt | wc
diff serial.txt pthreads.txt | wc
