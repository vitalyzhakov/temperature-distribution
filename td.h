/* 
 * File:   td.h
 * Author: zhakov
 *
 * Created on 26 Сентябрь 2012 г., 23:08
 */

#ifndef TD_H
#define	TD_H

class td {
public:
    td(int N, float dimensity, float Eps, float RhoSqr, float tau);
    td(const td& orig);
    virtual ~td();
    int initMatrix();
    int solve();
    int printMatrix(float ***matrix);
    int printMatrixToFile(char *path);
    int enableDebugMode(unsigned step);
    int disableDebugMode();
    int getIterationsCount();
    int printMatrixToFileGnuPlotFormat(char *path, int groupCount);
    
private:
    float ***solution;
    float ***tempSolution;
    float xStep; //Шаг по x
    float yStep; //Шаг по y
    float zStep; //Шаг по z
    int N; //Количество шагов
    float Eps; //Точность вычислений
    float RhoSqr; //ро квадрат
    float h; //Шаг по координате
    float tau; //Шаг по времени
    int iterationsCount; //Счётчик количества итераций
    float dimencity; //Длина ребра куба
    bool debugMode; //Включение отладки
    unsigned debugStep; //Шаг вывода для отладки
    unsigned maxIterations; //Максимальное количество итераций
    
    int boundaryCalculate(float ***matrix); //Подсчёт значений на границах
    
    float t0yz(float t, int y, int z, float*** matrix);
    float tNyz(float t, int y, int z, float*** matrix);
    float tx0z(float t, int x, int z, float*** matrix);
    float txNz(float t, int x, int z, float*** matrix);
    float txy0(float t, int x, int y, float*** matrix);
    float txyN(float t, int x, int y, float*** matrix);
};

#endif	/* TD_H */