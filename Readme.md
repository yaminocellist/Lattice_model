# Run in Macos, use:
g++ -I/opt/homebrew/Cellar/boost/1.78.0_1/include -I/opt/homebrew/Cellar/eigen/3.4.0_1/include *.cpp -std=c++14 && ./a.out 4 7 12.8 1 12 13 11.3 1 2 1 59.7 1

# Run in WSL, first change the including phases of BOOST in "functions.h",
# then run the below magic:

g++ *.cpp -I /usr/include/eigen3 -std=c++14 -lstdc++fs && ./a.out 4 7 12.8 1 12 13 11.3 1 2 1 59.7 1
