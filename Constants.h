#pragma once
#include<cmath>

extern double L; // размер пространства
extern double eps; // сколько частица не долетает до правой стенки

extern int num_ceil; // число ячеек
extern double dx; // шаг по сетке

extern double T; // всего времени
extern double dt; // шаг по времени

//Желаемая концентрация 1e17
extern double n_2; // сколько частиц моделирует 1 частица

extern double eps_0;


//константа для перевода посчитанных rho в правую часть уравнения Пуассона
extern double denom;

//число волн в распределении электронов
extern int k;


// int NX = L / dx; // в какой ячейке частица

extern double E_0;
extern double B_0;

extern int N; // число частиц

extern double v_t_ions; //тепловая скорость ионов
extern double v_t_el; //тепловая скорость электронов



extern double m_ion;
extern double m_el;
extern double q;
