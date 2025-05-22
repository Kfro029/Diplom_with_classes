#include "Constants.h"

double L = 0.04; // размер пространства
double eps = 0 * L * 1.0e-16; // сколько частица не долетает до правой стенки

int num_ceil = 100; // число ячеек
double dx = L / (num_ceil * 1.0 - 1.0); // шаг по сетке

double T = 3.5e-6; // всего времени
double dt = 3.5e-12; // шаг по времени
int shot = 1000;

//Желаемая концентрация 1e17
double n_2 = 4e11; // сколько частиц моделирует 1 частица

double n_ntr = 4e5; // сколько нейтронов моделирует 1 частица

double eps_0 = 8.85e-12;


//константа для перевода посчитанных rho в правую часть уравнения Пуассона
double denom = dx * dx * eps_0;

//число волн в распределении электронов
int k = 8;


// int NX = L / dx; // в какой ячейке частица

double E_0 = 1.0 * pow(10, -12);
//double B_0 = 1.0;

int N = num_ceil * 100; // число частиц


int N_ntr = num_ceil * 1;
//double v_t_el = L / (10 * dt);
//double v_t_ions = 0.0236 * v_t_el;


double m_ion = 1.67 * pow(10, -27);
double m_el = 9.31 * pow(10, -31);
double q = 1.6 * pow(10, -19);


double kB = 1.380649e-23;   // Постоянная Больцмана (Дж/К)
double eV_to_K = 11604.525; // 1 эВ = 11604.525 K


double T_eV = 10; // эВ

double v_t_el = std::sqrt(3 * kB * T_eV * eV_to_K / m_el);
double v_t_ions = std::sqrt(3 * kB * T_eV * eV_to_K / m_ion);


double BPHI0 = 0.00215;
double BPHIWDTH = 5.e-7;

double emission = 0.95 * L;
int ceil_emission = emission / dx;
double loc_emission = fmod(emission, dx);