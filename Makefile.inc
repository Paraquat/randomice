# Library paths
LPATH = 
IPATH = -I/usr/local/Cellar/eigen/3.3.4/include/eigen3/

# Debugging
DEBUG =  -g -Wall

# Optimisations
OPT = 

# Compiler
CC = g++ $(OPT) $(DEBUG) $(LPATH) $(IPATH)