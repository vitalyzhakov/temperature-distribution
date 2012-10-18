/* 
 * File:   main.cpp
 * Author: zhakov
 * Явная схема для трёхмерного уравнения теплопроводности
 *
 * Created on 26 Сентябрь 2012 г., 21:05
 */

#include <math.h>
#include <stdio.h>
#include <omp.h>


#include "td.h"
/**
 * Entry point
 */
int main(int argc, char** argv) {
    
    //Замер времени
    double startTime = omp_get_wtime();
    
    //Класс-решение
    td* solution = new td(100, 1, 0.00001, 0.003, 0.000001);
    
    solution->initMatrix();
    solution->enableDebugMode(100);
    solution->solve();
    
    printf("IterationsCount = %d\n", solution->getIterationsCount());
    
    solution->printMatrixToFile("outfile-full.txt");
    
    
    double finishTime = omp_get_wtime();
    
    printf("Calculation time: %f seconds\n", finishTime - startTime);

    return 0;
}