g++ -fopenmp main.cpp -O3 # -std=c++11

./a.out structure/test.xyz structure/ROI_test.txt -o output/test -s 0.5
./a.out structure/biomembrane_complete.pdb structure/ROI.txt -s 0.1 -o output/biomembrane_complete -vs 2
./a.out structure/biomembrane_missingCHL.pdb structure/ROI.txt -s 0.1 -o output/biomembrane_missingCHL -vs 2