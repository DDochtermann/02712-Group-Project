# 02712-Group-Project

baseA.c and baseB.c are code to reproduce the results in Figure 2A and 2B of the paper 
**Evolutionary Rescue Through Partly HeritablePhenotypic Variability**. To run the code, use the following commands and replace the `variance`
and `memory_prob` with corresponding numbers:

```
gcc -o baseA.o baseA.c -lgsl -lgslcblas -lm 
./baseA.o variance memory_prob
```


