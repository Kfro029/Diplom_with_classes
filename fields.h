#pragma once
#include <iostream>
#include <cmath>
#include <vector>
// #include "Constants.h"

class Particles;


class Field {
private:
    std::vector<double> E;
    std::vector<double> p; // для прогонки Пуассона
    std::vector<double> q; // для прогонки Пуассона

    std::vector<double> fi; // потенциал электрического поля

    std::vector<double> rho; // плотность заряда
    std::vector<double> rho_el; // плотность заряда электронов
    std::vector<double> rho_ions; // плотность заряда протонов

public:
    Field(); // Конструктор

    double field_by_x(double x); //вернуть значение поля по координате x
    void solve_field(Particles& el, Particles& ions); // Решение уравнения Пуассона -> находим fi
    void calc_E();

    void fill_null_field(); //заполнить все нулями после Move/SETV


};
