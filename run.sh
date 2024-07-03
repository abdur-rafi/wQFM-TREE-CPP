# rm -f a.out
g++ main.cpp -O3 -march=native
# ./a.out input.txt input.txt output.txt
./a.out ../run/astral2/estimated-gene-trees/model.50.2000000.0.000001/03/gt-cleaned ../run/astral2/estimated-consensus-trees/model.50.2000000.0.000001/03/cons-paup.tre output.txt 2> log.log