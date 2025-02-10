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

    std::vector<double> rho;
    

    Particles(double m, double q, std::vector<double> want_x, std::vector<double> want_v_x, std::vector <double> want_v_y); //  онструктор


    void move(Field& field); // ƒвижение
    void SETV(Field& field); //начальное даижение на -t/2
    void CIC(); //расчет плотности зар€да
    std::vector<double> give_rho(); //выдаем плотность дл€ ѕуассона в пол€х

    void fill_null_part(); //после Move запускать это действо дл€ массивов rho
    // std::vector<double> give_x(); // верунть координаты, например, дл€ расчета полей


};
