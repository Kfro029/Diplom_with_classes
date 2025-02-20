#include "fieldB.h"
#include <fstream>
#include "cmath"
#include "Constants.h"


// –еализаци€ конструктора
FieldB::FieldB()
	: B(num_ceil, 1.0)
{}


// –еализаци€ методов
void FieldB::loadFromFile(std::string& filename) {
	std::ifstream file(filename);
	if (!file) {
		std::cerr << "ќшибка: не удалось открыть файл " << filename << std::endl;
		return;
	}

	std::vector<double> file_data;
	double value;
	while (file >> value) {
		file_data.push_back(value);
	}

	int file_size = file_data.size();
	int B_size = num_ceil;

	if (file_size == B_size) {
		B = file_data; // ≈сли размеры совпадают, просто копируем
	}
	else {
		std::cout << "»нтерпол€ци€ значений B..." << std::endl;

		double dx_file = 1.0 / (file_size - 1);  // Ўаг сетки в файле (нормализованный)
		double dx_grid = dx;     // Ўаг сетки границ €чеек

		for (int i = 0; i < B_size; i++) {
			double x_grid = i * dx_grid;  //  оордината в новой сетке

			// ќпредел€ем индексы в file_data дл€ интерпол€ции
			int idx = static_cast<int>(x_grid / dx_file);
			double alpha = (x_grid - idx * dx_file) / dx_file;  // ¬ес линейной интерпол€ции

			if (idx >= file_size - 1) {
				B[i] = file_data.back();  //  райний правый случай (x = L)
			}
			else {
				B[i] = (1 - alpha) * file_data[idx] + alpha * file_data[idx + 1];  // Ћинейна€ интерпол€ци€
			}
		}
	}
}

double FieldB::b_by_x(double x) {
	int ceil_x = x / dx; //номер €чейки дл€ частицы
	double x_loc = fmod(x, dx); //координата в €чейке

	double B_in_x = (B[ceil_x] * (dx - x_loc) + B[ceil_x + 1] * x_loc) / dx;
	//std::cout << E_in_x << '\t' << x << "\n";

	//return 1.0e9;
	return B_in_x;
}