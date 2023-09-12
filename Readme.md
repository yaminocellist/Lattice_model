# Based on Teif's Lattice model papers;

The partition function (normalization coefficient) is easy to calculate through:

Z = &Sigma;<sub>i</sub>(e<sup>-&beta;H<sub>i</sub></sup>)

where the summation is over all possible combinations of DNA and proteins;
Then the probability of i-th combination is given by:

p<sub>i</sub> = 1/Z * e<sup>-&beta;H<sub>i</sub></sup>

# Input format:
(Protein 1's) mg, Vg, binding constant, concentration (Protein 2's) mg, Vg, binding constant, concentration, ... (cycling), target n(th segment of DNA), target g(th protein);

# Run in Macos, use:
g++ -I/opt/homebrew/Cellar/boost/1.78.0_1/include -I/opt/homebrew/Cellar/eigen/3.4.0_1/include *.cpp -std=c++14 && ./a.out 4 7 12.8 10 12 13 11.3 10 2 1 59.7 10 3 0

# Run in WSL, first change the including phases of BOOST in "functions.h", then run the below magic:

g++ *.cpp -I /usr/include/eigen3 -std=c++14 -lstdc++fs && ./a.out 4 7 12.8 10 12 13 11.3 10 2 1 59.7 10 3 0

# To give the system permission of editing, use:
sudo chown -R yaminocellist:yaminocellist Lattice_model
