/* 
 * File:   main.cpp
 * Author: zhakov
 * Явная схема для трёхмерного уравнения теплопроводности
 *
 * Created on 26 Сентябрь 2012 г., 21:05
 */

#include <math.h>
#include <stdio.h>

#include "td.h"
/**
 * Entry point
 */
int main(int argc, char** argv) {
    
    td* solution = new td(100, 1, 0.0001, 0.003, 0.00001);
    
    solution->initMatrix();
    solution->enableDebugMode(100);
    solution->solve();
    
    printf("IterationsCount = %d\n", solution->getIterationsCount());
    
    
    solution->printMatrixToFileGnuPlotFormat("outfile.txt", 10);

    return 0;
}