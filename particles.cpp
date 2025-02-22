#include "particles.h" 
#include "fields.h"
// #include "fieldB.h"
#include <random>
#include <cmath>
#include <fstream>

// Реализация конструктора
Particles::Particles(double m, double q, double v_t, int seed)
	: m(m), q(q), rho(num_ceil, 0.0), temp_t(0.0), temp_s(0.0), temp_c(0.0), v_t(v_t), gen(seed), dist(0.0, 1.0)
{
	x.resize(N);
	v_x.resize(N);
	v_y.resize(N);
	//std::cout << "q in constructor = " << q << std::endl;
}

// Реализация методов

//ее надо будет запускать после Move
void Particles::fill_null_part() {
	std::fill(rho.begin(), rho.end(), 0.0);
};


void Particles::swap_and_delete(int number) {

	if (x.empty() || v_x.empty() || v_y.empty()) return;

	if (number < x.size() - 1) {
		x[number] = x.back();
		v_x[number] = v_x.back();
		v_y[number] = v_y.back();
	}

	x.resize(x.size() - 1);
	v_x.resize(v_x.size() - 1);
	v_y.resize(v_y.size() - 1);
};



void Particles::move(Field& fieldE) {
	
	for (int i = x.size() - 1; i >= 0; i--) {

		//расчет скорости на dt/2
		v_x[i] += fieldE.field_by_x(x[i]) * q / m * (dt / 2.0);
		

		//расчет поворота частицы
		double temp_t = q * 1.0 / m * dt / 2;
		double temp_s = 2 * temp_t / (1 + temp_t * temp_t);
		double temp_c = (1 - temp_t * temp_t) / (1 + temp_t * temp_t);

		double v_x_new = v_x[i] * temp_c + v_y[i] * temp_s;
		v_y[i] = -v_x[i] * temp_s + v_y[i] * temp_c;
		v_x[i] = v_x_new;

		//еще подвинулись на dt/2
		v_x[i] += fieldE.field_by_x(x[i]) * q / m * (dt / 2.0);
		x[i] += v_x[i] * dt;
		

		
		/*
		// ГУ вылет и все
		if (this->x[i] >= L || this->x[i] < 0) {
			swap_and_delete(i);
			//i--;
		}
		*/
	
		
		// граничные условия: жесткая стенка
		if (x[i] >= L - eps) {
			x[i] = (L - eps) - (x[i] - (L - eps));
			v_x[i] *= (-1);
		}
		if (x[i] - eps <= 0) {
			x[i] *= (-1);
			v_x[i] *= (-1);
		}
		
		

	}

	// fill_null_part();
	// field.fill_null_field();
}

void Particles::CIC() {
	for (std::size_t i = 0; i < x.size(); i++) {


		int ceil_left = (x[i])  / dx;
		int ceil_right = ceil_left + 1;		
		//int ceil_right = std::min(static_cast<size_t>(ceil_left + 1), x.size()-1);
		
		double x_loc = fmod(x[i], dx);
		//std::cout << "x[i]\t" << x[i] << "\ti\t" << i;
		//std::cout << "\t left \t" << ceil_left << "\tright \t " << ceil_right << "\n";
		
		rho[ceil_left] += q * (dx - x_loc);
		rho[ceil_right] += q * x_loc;
		
	}
	//std::cout << rho[1];

	for (std::size_t i = 0; i < rho.size(); i++) {
		rho[i] *= (-n_2) / denom;
		//std::cout << rho[i] << "\t" << i << "\n";
	}
	
}

void Particles::SETV(Field& fieldE) {
	for (int i = x.size() - 1; i >= 0; i--) {

		//расчет скорости на dt/2
		v_x[i] += fieldE.field_by_x(x[i]) * (-q/ 2.) / m * (dt / 2.0);


		//расчет поворота частицы
		double temp_t = q * 1.0 / m * dt / 2;
		double temp_s = 2 * temp_t / (1 + temp_t * temp_t);
		double temp_c = (1 - temp_t * temp_t) / (1 + temp_t * temp_t);

		double v_x_new = v_x[i] * temp_c + v_y[i] * temp_s;
		v_y[i] = -v_x[i] * temp_s + v_y[i] * temp_c;
		v_x[i] = v_x_new;

		//еще подвинулись на dt/2
		v_x[i] += fieldE.field_by_x(x[i]) * (-q/2.0) / m * (dt / 2.0);
		x[i] += v_x[i] * dt;



		/*
		// ГУ вылет и все
		if (x[i] >= L || x[i] < 0) {
			swap_and_delete(i);
			//i--;
		}
		*/


		// граничные условия: жесткая стенка
		if (x[i] >= L - eps) {
			x[i] = (L - eps) - (x[i] - (L - eps));
			v_x[i] *= (-1);
		}
		if (x[i] - eps <= 0) {
			x[i] *= (-1);
			v_x[i] *= (-1);
		}



	}
}

double Particles::get_random() {
	return dist(gen);
}



void Particles::ionization_first() {
	x.push_back(L * get_random());

	R_s = get_random();
	R_theta = get_random();


	v_x.push_back(v_t * std::sqrt(-2 * std::log(R_s)) * std::cos(2 * 3.1416 * R_theta));
	v_y.push_back(v_t * std::sqrt(-2 * std::log(R_s)) * std::sin(2 * 3.1416 * R_theta));
};


void Particles::fill() {
	for (std::size_t i = 0; i < N; i++) {
		ionization_first();
	}
}

void Particles::ionization() {
	for (std::size_t i = 0; i < ionaze.size(); i++) {

		// целая часть ионизации в ячейке
		for (int j = 0; j < int(ionaze[i]); j++) {
			x.push_back((i + get_random()) * dx);

			R_s = get_random();
			R_theta = get_random();

			v_x.push_back(v_t * std::sqrt(-2 * std::log(R_s)) * std::cos(2 * 3.1416 * R_theta));
			v_y.push_back(v_t * std::sqrt(-2 * std::log(R_s)) * std::sin(2 * 3.1416 * R_theta));
		}
		//дробная
		if (get_random() < fmod(ionaze[i], 1.0)) {
			x.push_back((i + get_random()) * dx);

			R_s = get_random();
			R_theta = get_random();

			v_x.push_back(v_t * std::sqrt(-2 * std::log(R_s)) * std::cos(2 * 3.1416 * R_theta));
			v_y.push_back(v_t * std::sqrt(-2 * std::log(R_s)) * std::sin(2 * 3.1416 * R_theta));
		}

	}

}


void Particles::loadFromFile(std::string filename) {
	std::ifstream file(filename);
	if (!file) {
		std::cerr << "Don't find ionization func in file: " << filename << std::endl;
		return;
	}

	std::vector<double> file_data;
	double value;
	while (file >> value) {
		file_data.push_back(value);
	}

	int file_size = file_data.size();
	int B_size = num_ceil;

	if (file_size == num_ceil) {
		ionaze = file_data; // 
	}
	else {
		std::cout << "start process of interpolation ionaze function..." << std::endl;

		double dx_file = 1.0 / (file_size - 1);  // ??? ????? ? ????? (???????????????)
		double dx_grid = dx;     // ??? ????? ?????? ?????

		for (int i = 0; i < B_size; i++) {
			double x_grid = i * dx_grid;  // ?????????? ? ????? ?????

			// ?????????? ??????? ? file_data ??? ????????????
			int idx = static_cast<int>(x_grid / dx_file);
			double alpha = (x_grid - idx * dx_file) / dx_file;  // ??? ???????? ????????????

			if (idx >= file_size - 1) {
				ionaze[i] = file_data.back();  // ??????? ?????? ?????? (x = L)
			}
			else {
				ionaze[i] = (1 - alpha) * file_data[idx] + alpha * file_data[idx + 1];  // ???????? ????????????
			}
		}
	}
}