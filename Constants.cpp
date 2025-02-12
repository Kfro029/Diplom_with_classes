#include "Constants.h"
#include<cmath>

double L = 0.04; // размер пространства
double eps = 1 * L * 1.0e-16; // сколько частица не долетает до правой стенки

int num_ceil = 100; // число ячеек
double dx = L / (num_ceil * 1.0 - 1.0); // шаг по сетке

double T = 3.5e-9; // всего времени
double dt = 3.65e-12; // шаг по времени

//Желаемая концентрация 1e17
double n_2 = 4e11; // сколько частиц моделирует 1 частица

double eps_0 = 8.85e-12;


//константа для перевода посчитанных rho в правую часть уравнения Пуассона
double denom = dx * dx * eps_0;

//число волн в распределении электронов
int k = 8;


// int NX = L / dx; // в какой ячейке частица

double E_0 = 1.0 * pow(10, -12);
double B_0 = 1.0;

int N = num_ceil * 50; // число частиц


double m_ion = 1.67 * pow(10, -27);
double m_el = 9.31 * pow(10, -31);
double q = 1.6 * pow(10, -19);