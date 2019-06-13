# Library paths
LPATH = -L/Users/zamaan/Software/backward-cpp/
IPATH = -I/usr/local/Cellar/eigen/3.3.7/include/eigen3/ -I/usr/local/Cellar/gsl/2.4/include/ -I/usr/local/Cellar/boost/1.70.0/include/

# Debugging
#DEBUG = -g -Wall

# Optimisations
OPT = -O3

# Compiler
CC = g++ $(OPT) $(DEBUG) $(IPATH)
