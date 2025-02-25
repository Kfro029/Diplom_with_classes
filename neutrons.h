#pragma once
// #include "fieldB.h"
#include <random>
#include <cmath>
#include <fstream>


class Neutrons {
private:
public:
    double m;
    std::vector<double> x;
    std::vector<double> v_x;
    std::vector<double> v_y;

    std::vector<double> n;

    // тангенс, синус и косинус для вращения
    double temp_t;
    double temp_s;
    double temp_c;

    // случайные числа для скоростей
    double R_s;
    double R_theta;

    double v_t; //тепловая скорость

    int out; //количество вылетевших частиц

    std::vector<double> rho;

    //генератор случ чисел
    std::mt19937_64 gen;
    std::uniform_real_distribution<double> dist;

    //вектор ионизации
    std::vector<double> ionaze;

    //для интерполяции
    int idx;
    double alpha;
    double x_grid;

    //вектор вероятностей
    int x_ceil;
    double x_loc;



    //методы
    // 
    //загрузка файла
    //void loadFromFile(std::string filename);



    Neutrons(double m, double v_t, int seed); // Конструктор

    void fill(); // заполнение 

    void move(); // Движение
    void concentration(); //расчет плотности заряда
    //std::vector<double> give_rho(); //выдаем плотность для Пуассона в полях

    void swap_and_delete(int number);

    void ionization_first();
    void ionization();
    // std::vector<double> give_x(); // верунть координаты, например, для расчета полей

    double get_random();

    //double if_collis(double x_part);

};