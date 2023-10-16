# Overlap_Tool
Simple Program to display overlap integrals between STO shells on atom pairs.

Handles S-D shells with atomic numbers from H-Kr

--------------------------------------

All modern C++ compilers should successfully compile the program. Because the program uses header files, no Makefile is provided here.

Simply run (in the case of GNU g++)

g++ -c -Iinclude/ overlap_tool.cpp

g++ -o overlap_tool overlap_tool.o

--------------------------------------

This makes an executable called overlap_tool

To run, simply use:

./overlap_tool


