#pragma once
#include <iostream>
#include <cmath>
#include <vector>
#include "Constants.h"

class Field;


class Particles {
private:
public:
    double m;
    double q;
    std::vector<double> x;
    std::vector<double> v_x;
    std::vector<double> v_y;

    double v_x_new; //������� ��� ������� ���� (���. 18 �-�)

    // �������, ����� � ������� ��� ��������
    double t; 
    double s;
    double c;

    std::vector<double> rho;
    

    Particles(double m, double q, std::vector<double> want_x, std::vector<double> want_v_x, std::vector <double> want_v_y); // �����������


    void move(Field& field); // ��������
    void SETV(Field& field); //��������� �������� �� -t/2
    void CIC(); //������ ��������� ������
    std::vector<double> give_rho(); //������ ��������� ��� �������� � �����

    void fill_null_part(); //����� Move ��������� ��� ������� ��� �������� rho

    void swap_and_delete(int number);

    void ionization_first();
    // std::vector<double> give_x(); // ������� ����������, ��������, ��� ������� �����


};
