#include "Constants.h"
#include<cmath>

double L = 0.04; // ������ ������������
double eps = 1 * L * 1.0e-16; // ������� ������� �� �������� �� ������ ������

int num_ceil = 100; // ����� �����
double dx = L / (num_ceil * 1.0 - 1.0); // ��� �� �����

double T = 3.5e-9; // ����� �������
double dt = 3.65e-12; // ��� �� �������

//�������� ������������ 1e17
double n_2 = 4e11; // ������� ������ ���������� 1 �������

double eps_0 = 8.85e-12;


//��������� ��� �������� ����������� rho � ������ ����� ��������� ��������
double denom = dx * dx * eps_0;

//����� ���� � ������������� ����������
int k = 8;


// int NX = L / dx; // � ����� ������ �������

double E_0 = 1.0 * pow(10, -12);
double B_0 = 1.0;

int N = num_ceil * 50; // ����� ������


double m_ion = 1.67 * pow(10, -27);
double m_el = 9.31 * pow(10, -31);
double q = 1.6 * pow(10, -19);