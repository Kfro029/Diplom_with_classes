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

    double v_x_new; //костыль для решения СЛАУ (стр. 18 Б-Л)

    // тангенс, синус и косинус для вращения
    double t; 
    double s;
    double c;

    std::vector<double> rho;
    

    Particles(double m, double q, std::vector<double> want_x, std::vector<double> want_v_x, std::vector <double> want_v_y); // Конструктор


    void move(Field& field); // Движение
    void SETV(Field& field); //начальное даижение на -t/2
    void CIC(); //расчет плотности заряда
    std::vector<double> give_rho(); //выдаем плотность для Пуассона в полях

    void fill_null_part(); //после Move запускать это действо для массивов rho

    void swap_and_delete(int number);

    void ionization_first();
    // std::vector<double> give_x(); // верунть координаты, например, для расчета полей


};
