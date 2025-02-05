# UnexpectedDtourAhead
A critical evaluation of D statistics and the theory of Neanderthal introgression

This project critically examines the assumptions behind D-statistics, a widely used method to infer archaic introgression from Neanderthals into modern humans. The analysis reveals that key assumptions—constant mutation rates and negligible recurrent mutations—are invalid. 

C++ and Python codes to reproduce all the figures in the manuscript are available in the enclosed folders.
- The C++ script Source1\Source.cpp generates Tables 2 and 3 and data for Figures 7, 8, 9, 10 and S2. 
- The C++ script Source2\Source2n.cpp generates data for Figures 4, 5, 6, 11 and S3
- The Python script MutationCounts\PlotWindowPercents.py generates Figure 2. 
- The Python script MutationCounts\PlotWindowPercents_all_chrs.py generates Figure S1. 

The folders Source1 and Source2 include the source files in C++, as well as compiled files for Windows and Linux.
Each folder also includes the outcome of the code for chromosome 22 (chr22.zip).
Note, that Source2 requires the HGDP file. Due to its size, we provided a file with the first 20K rows.

