# Library paths
LPATH = 
IPATH = -I/usr/local/Cellar/eigen/3.3.4/include/eigen3/ -I/usr/local/Cellar/gsl/2.4/include/ -I/usr/local/Cellar/boost/1.64.0_1/include/

# Debugging
DEBUG =  -g -Wall

# Optimisations
OPT = 

# Compiler
CC = g++ $(OPT) $(DEBUG) $(LPATH) $(IPATH)