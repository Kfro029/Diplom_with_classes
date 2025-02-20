#pragma once
#include <iostream>
#include <cmath>
#include <vector>
// #include "Constants.h"



class FieldB {
private:
public:
    std::vector<double> B;

    FieldB(); // Конструктор

    void loadFromFile(std::string& filename);

    double b_by_x(double x);

};
