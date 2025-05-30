#include "fields.h"

// –еализаци€ конструктора
Field::Field()
	: E(num_ceil, 0.0),
	p(num_ceil - 1, 0.0), q(num_ceil - 1, 0.0),
	fi(num_ceil, 0.0), rho(num_ceil, 0.0), rho_el(num_ceil, 0.0), rho_ions(num_ceil, 0.0),
	B(num_ceil, 0.0)
{}

// –еализаци€ методов
double Field::field_by_x(double x) {
	int ceil_x = x / dx; //номер €чейки дл€ частицы
	double x_loc = fmod(x, dx); //координата в €чейке
	
	double E_in_x = (E[ceil_x] * (dx - x_loc) + E[ceil_x + 1] * x_loc) / dx;
	//std::cout << E_in_x << '\t' << x << "\n";
	
	return E_in_x;
}


void Field::calc_E() {
	E[0] = (fi[0] - fi[1]) / dx;
	for (int i = 1; i + 1 < fi.size(); i++) {
		E[i] = (fi[i - 1] - fi[i + 1]) / (2.0 * dx);
		//std::cout << E[i] << "\ti" << "\n";
	}
	E[E.size() - 1] = (fi[fi.size() - 2] - fi[fi.size() - 1]) / dx;

};

void Field::fill_null_field() {
	std::fill(rho.begin(), rho.end(), 0.0);
	// std::fill(fi.begin(), fi.end(), 0.0);
}

void Field::solve_field(Particles& el, Particles& ions) {

	//rho_el = el.rho;
	//rho_ions = ions.rho;

	for (std::size_t i = 0; i < rho.size(); i++) {
		rho[i] = el.rho[i] + ions.rho[i];
	}

	//√”

	fi[0] = 200.0;
	fi[fi.size() - 1] = 0.0;

	//пр€мой ход
	p[2] = 1.0 / 2.0;
	q[2] = (-rho[1] + fi[0]) / (2.0);


	for (std::size_t i = 2; i + 2 < fi.size(); i++) {
		p[i + 1] = -1.0 / (p[i] - 2.0);
		q[i + 1] = (rho[i] - q[i]) / (p[i] - 2.0);
	}


	fi[fi.size() - 2] = (rho[fi.size() - 2] - fi[fi.size() - 1] - q[fi.size() - 2]) / (p[fi.size() - 2] - 2.0);

	//обратный ход
	for (int i = fi.size() - 3; i > 0; i--) {
		fi[i] = p[i + 1] * fi[i + 1] + q[i + 1];
	}
	
	
	fi_em = (fi[ceil_emission] * (dx - loc_emission) + fi[ceil_emission + 1] * loc_emission) / dx;

	for (std::size_t i = 0; i < fi.size(); i++) {
		fi[i] -= fi_em / emission * (i * dx);
	}
	
	calc_E();
}

void Field::loadFromFile(std::string filename) {
	//магнитное поле задаетс€ в √с, а мы потом переводим в “л
	std::ifstream file(filename);
	if (!file) {
		std::cerr << "Don't find field B in file: " << filename << std::endl;
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
		B = file_data;
	}
	else {
		std::cout << "start process of interpolation B..." << std::endl;
		// B.resize(B_size);  

		double dx_file = L / (file_size - 1);  // »справлено
		double dx_grid = dx;

		for (int i = 0; i < B_size; i++) {
			double x_grid = i * dx_grid;

			int idx = static_cast<int>(x_grid / dx_file);
			double alpha = (x_grid - idx * dx_file) / dx_file;

			if (idx >= file_size - 1) {
				B[i] = file_data.back();
			}
			else {
				B[i] = (1 - alpha) * file_data[idx] + alpha * file_data[idx + 1];
			}
		}
	}
	
	for (std::size_t i = 0; i < B.size(); i++) {
		B[i] = B[i] * 1. / 10000.;
	}
}

double Field::b_by_x(double x) {
	int ceil_x = x / dx; //????? ?????? ??? ???????
	double x_loc = fmod(x, dx); //?????????? ? ??????

	double B_in_x = (B[ceil_x] * (dx - x_loc) + B[ceil_x + 1] * x_loc) / dx;
	//std::cout << E_in_x << '\t' << x << "\n";

	//return 1.0e9;
	return B_in_x;
}


