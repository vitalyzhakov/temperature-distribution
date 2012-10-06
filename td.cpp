/* 
 * File:   td.cpp
 * Author: zhakov
 * 
 * Created on 26 Сентябрь 2012 г., 23:08
 */

#include "td.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>

//Конструктор

td::td(int N = 100, float dimensity = 1, float Eps = 0.001, float RhoSqr = 0.003, float tau = 0.00001) {


    this->N = N;
    this->dimencity = dimensity;
    this->Eps = Eps;
    this->RhoSqr = RhoSqr;
    this->tau = tau;

    this->h = this->dimencity / this->N;

    debugMode = false;
    maxIterations = 100000;
}

td::td(const td& orig) {
}

td::~td() {
}

/**
 * Включить отладочный режим (вывод максимальной разницы)
 * @param step - шаг, через который выводить отладочную информацию
 * @return 
 */
int td::enableDebugMode(unsigned step = 100) {
    debugMode = true;
    debugStep = step;
    return 0;
}

/**
 * Отменить режим отладки
 * @return 
 */
int td::disableDebugMode() {
    debugMode = false;
    return 0;
}

/**
 * Получить количество проделанных итераций
 * @return 
 */
int td::getIterationsCount() {
    return iterationsCount;
}

/**
 * Инициализируем нулями матрицы
 * Нужно перед проведением рассчётов
 * @return 
 */
int td::initMatrix() {

    //Инициализация массива для адекватной работы с памятью
    solution = new float**[N];
    tempSolution = new float**[N];
    for (int i = 0; i < N; ++i) {
        solution[i] = new float*[N];
        tempSolution[i] = new float*[N];

        for (int j = 0; j < N; ++j) {
            solution[i][j] = new float[N];
            tempSolution[i][j] = new float[N];
        }
    }

    //Заполнение нулями
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            for (int k = 0; k < N; k++) {
                solution[i][j][k] = 0;
                tempSolution[i][j][k] = 0;
            }
        }
    }

    //Считаем длину шага по каждой компненте
    xStep = dimencity / N;
    yStep = dimencity / N;
    zStep = dimencity / N;
    

    return 0;
}

/**
 * Подсчёт значений на границах куба
 * @param matrix
 * @return 
 */
int td::boundaryCalculate(float*** matrix) {
    //x = const
    for (int j = 0; j < N; j++) {
        for (int k = 0; k < N; k++) {
            float t = iterationsCount * tau;

            matrix[0][j][k] = t0yz(t, j, k, matrix);
            matrix[N - 1][j][k] = tNyz(t, j, k, matrix);
        }
    }

    //y = const
    for (int i = 0; i < N; i++) {
        for (int k = 0; k < N; k++) {
            float t = iterationsCount * tau;

            matrix[i][0][k] = tx0z(t, i, k, matrix);
            matrix[i][N - 1][k] = txNz(t, i, k, matrix);
        }
    }


    //z = const
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            float t = iterationsCount * tau;

            matrix[i][j][0] = txy0(t, i, j, matrix);
            matrix[i][j][N - 1] = txyN(t, i, j, matrix);
        }
    }
    
    
    return 0;
}

int td::solve() {

    float maxDifference = 1;
    iterationsCount = 0;
    
    //Считаем значения на границе
    boundaryCalculate(solution);
    
    while ((maxDifference > Eps) && (iterationsCount < maxIterations)) {
        iterationsCount++;
        maxDifference = 0;

        //Считаем значения на границе
        boundaryCalculate(tempSolution);
        
        //Считаем во внутренних точках
        for (int i = 1; i < N - 1; i++) {
            for (int j = 1; j < N - 1; j++) {
                for (int k = 1; k < N - 1; k++) {
                    //Разница
                    float difference = RhoSqr * tau / h / h * (solution[i + 1][j][k] + solution[i - 1][j][k] + solution[i][j + 1][k] + solution[i][j - 1][k] + solution[i][j][k + 1] + solution[i + 1][j][k - 1] - 6 * solution[i][j][k]);

                    //Новое значение
                    tempSolution[i][j][k] += difference;
                    //Условие по eps
                    if (difference > maxDifference) {
                        maxDifference = difference;
                    }
                }
            }
        }
        //Если включена отладка, выводим найденную разницу через заданное количество шагов
        if (debugMode) {
            if (iterationsCount % debugStep == 1) {
                printf("maxDiff %.7f on iteration %d\n", maxDifference, iterationsCount);
            }
        }

        //Обмениваем указатели массивов
        float ***tempArr;
        tempArr = solution;
        solution = tempSolution;
        tempSolution = tempArr;
    }

    return 0;
}

/**
 * Краевое условие при x = 0
 * @param t
 * @param y
 * @param z
 * @return 
 */
float td::t0yz(float t, int y, int z, float*** matrix) {
    return matrix[1][y][z];
}

/**
 * Краевое условие при x = N
 * @param t
 * @param y
 * @param z
 * @param matrix
 * @return 
 */
float td::tNyz(float t, int y, int z, float*** matrix) {
    return matrix[N - 1][y][z];
}

/**
 * Краевое условие при z = 0
 * @param t
 * @param x
 * @param y
 * @param matrix
 * @return 
 */
float td::txy0(float t, int x, int y, float*** matrix) {
    return matrix[x][y][1];
}

/**
 * Краевое условие при z = N
 * @param t
 * @param x
 * @param y
 * @param matrix
 * @return 
 */
float td::txyN(float t, int x, int y, float*** matrix) {
    return matrix[x][y][N - 1];
}

/**
 * Краевое условие при y = 0
 * @param t
 * @param x
 * @param z
 * @param matrix
 * @return 
 */
float td::tx0z(float t, int x, int z, float*** matrix) {
    float floatX = x * xStep;
    float floatZ = z * zStep;
    return -sin(4 * M_PI * floatX) * sin(2 * M_PI * floatZ);
}

/**
 * Краевое условие при y = N
 * @param t
 * @param x
 * @param z
 * @param matrix
 * @return 
 */
float td::txNz(float t, int x, int z, float*** matrix) {
    return 0;
}

/**
 * Распечатать матрицу
 * 
 * @param matrix
 * @return 
 */
int td::printMatrix(float ***matrix) {
    printf("x y z f\n");
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            for (int k = 0; k < N; k++) {
                printf("%d %d %d %.3f\n", i, j, k, matrix[i][j][k]);
            }
        }
    }

    return 0;
}

/**
 * Записать получившуюся матрицу в файл
 * @param path
 * @return 
 */
int td::printMatrixToFile(char* path) {
    //Печатаем в файл
    FILE *outfile = fopen(path, "w");

    fprintf(outfile, "x     y     z      f\n");
    for (int i = 0; i < N; i++) {
        float x = i * xStep;
        for (int j = 0; j < N; j++) {
            float y = j * yStep;
            for (int k = 0; k < N; k++) {
                float z = k * zStep;
                fprintf(outfile, "%.5f %.5f %.5f %.5f\n", x, y, z, solution[i][j][k]);
            }
        }
    }

    fclose(outfile);

    return 0;
}