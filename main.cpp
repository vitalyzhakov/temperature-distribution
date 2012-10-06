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
    
    td* solution = new td(10, 1, 0.0001, 0.003, 0.001);
    
    solution->initMatrix();
    solution->enableDebugMode(100);
    solution->solve();
    
    printf("IterationsCount = %d", solution->getIterationsCount());
    
    
    solution->printMatrixToFile("outfile.txt");

    return 0;
}