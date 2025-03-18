#pragma once
#include <iostream>
#include <cmath>
#include <vector>
#include "Constants.h"
#include <random>
#include "fields.h"
#include <fstream>

class Field;


class Particles {
private:
public:
    double m;
    double q;
    std::vector<double> x;
    std::vector<double> v_x;
    std::vector<double> v_y;

    double v_x_new; //костыль дл€ решени€ —Ћј” (стр. 18 Ѕ-Ћ)

    // тангенс, синус и косинус дл€ вращени€
    double temp_t; 
    double temp_s;
    double temp_c;

    // случайные числа дл€ скоростей
    double R_s;
    double R_theta;

    double v_t; //теплова€ скорость

    int out_an; //количество вылетевших частиц с анода
    int out_kat; //количество вылетевших частиц с катода

    std::vector<double> rho;

    int fly;
    std::vector<int> divergention;

    //генератор случ чисел
    std::mt19937_64 gen;
    std::uniform_real_distribution<double> dist;
    
    //вектор ионизации
    std::vector<double> ionaze;

    //дл€ интерпол€ции
    int idx;
    double alpha;
    double x_grid;

    //вектор веро€тностей
    std::vector<double> collision;
    int x_ceil;
    double x_loc;



    //методы
    // 
    //загрузка файла
    void loadFromFile(std::string filename);



    Particles(double m, double q, double v_t, int seed); //  онструктор

    void fill(Field& fieldE); // заполнение 

    void move(Field& fieldE); // ƒвижение
    void SETV(Field& fieldE); //начальное даижение на -t/2
    void CIC(); //расчет плотности зар€да
    //std::vector<double> give_rho(); //выдаем плотность дл€ ѕуассона в пол€х

    void fill_null_part(); //после Move запускать это действо дл€ массивов rho

    void swap_and_delete(int number);

    void ionization_first();
    void ionization();
    void emission_ionization(int out);
    // std::vector<double> give_x(); // верунть координаты, например, дл€ расчета полей

    double get_random();

    double if_collis(double x_part);

};
