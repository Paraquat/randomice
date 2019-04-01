# Library paths
LPATH = -L/Users/zamaan/Software/backward-cpp/
IPATH = -I/usr/local/Cellar/eigen/3.3.5/include/eigen3/ -I/usr/local/Cellar/gsl/2.4/include/ -I/usr/local/Cellar/boost/1.67.0_1/include/

# Debugging
DEBUG =  -g -Wall

# Optimisations
OPT = -O3

# Compiler
CC = g++ $(OPT) $(DEBUG) $(IPATH)
