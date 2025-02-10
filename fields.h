#pragma once
#include <iostream>
#include <cmath>
#include <vector>
// #include "Constants.h"

class Particles;


class Field {
private:
    std::vector<double> E;
    std::vector<double> p; // ��� �������� ��������
    std::vector<double> q; // ��� �������� ��������

    std::vector<double> fi; // ��������� �������������� ����

    std::vector<double> rho; // ��������� ������
    std::vector<double> rho_el; // ��������� ������ ����������
    std::vector<double> rho_ions; // ��������� ������ ��������

public:
    Field(); // �����������

    double field_by_x(double x); //������� �������� ���� �� ���������� x
    void solve_field(Particles& el, Particles& ions); // ������� ��������� �������� -> ������� fi
    void calc_E();

    void fill_null_field(); //��������� ��� ������ ����� Move/SETV


};
