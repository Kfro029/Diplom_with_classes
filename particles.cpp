#include "particles.h" 
#include "fields.h"
#include "fieldB.h"
#include <random>
#include <cmath>

// Реализация конструктора
Particles::Particles(double m, double q, double v_t, int seed)
	: m(m), q(q), rho(num_ceil, 0.0), t(0.0), s(0.0), c(0.0), v_t(v_t), seed(seed)
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
		std::swap(x[number], x.back());
		std::swap(v_x[number], v_x.back());
		std::swap(v_y[number], v_y.back());
	}

	x.resize(x.size() - 1);
	v_x.resize(v_x.size() - 1);
	v_y.resize(v_y.size() - 1);
};



void Particles::move(Field& fieldE, FieldB& B) {
	
	for (int i = x.size() - 1; i >= 0; i--) {

		//расчет скорости на dt/2
		v_x[i] += fieldE.field_by_x(x[i]) * q / m * (dt / 2.0);
		

		//расчет поворота частицы
		double temp_t = this->q * B.b_by_x(x[i]) / this->m * dt / 2;
		double temp_s = 2 * temp_t / (1 + temp_t * temp_t);
		double temp_c = (1 - temp_t * temp_t) / (1 + temp_t * temp_t);

		double v_x_new = this->v_x[i] * temp_c + this->v_y[i] * temp_s;
		this->v_y[i] = -this->v_x[i] * temp_s + this->v_y[i] * temp_c;
		this->v_x[i] = v_x_new;

		//еще подвинулись на dt/2
		this->v_x[i] += fieldE.field_by_x(x[i]) * this->q / this->m * (dt / 2.0);
		this->x[i] += this->v_x[i] * dt;
		

		
		/*
		// ГУ вылет и все
		if (this->x[i] >= L || this->x[i] < 0) {
			swap_and_delete(i);
			//i--;
		}
		*/
	
		
		// граничные условия: жесткая стенка
		if (this->x[i] >= L - eps) {
			this->x[i] = (L - eps) - (this->x[i] - (L - eps));
			this->v_x[i] *= (-1);
		}
		if (this->x[i] - eps <= 0) {
			this->x[i] *= (-1);
			this->v_x[i] *= (-1);
		}
		
		

	}

	// fill_null_part();
	// field.fill_null_field();
}

void Particles::CIC() {
	for (std::size_t i = 0; i < this->x.size(); i++) {


		int ceil_left = (this->x[i])  / dx;
		int ceil_right = ceil_left + 1;		
		//int ceil_right = std::min(static_cast<size_t>(ceil_left + 1), x.size()-1);
		
		double x_loc = fmod(this->x[i], dx);
		//std::cout << "x[i]\t" << x[i] << "\ti\t" << i;
		//std::cout << "\t left \t" << ceil_left << "\tright \t " << ceil_right << "\n";
		
		this->rho[ceil_left] += this->q * (dx - x_loc);
		this->rho[ceil_right] += this->q * x_loc;
		
	}
	//std::cout << rho[1];

	for (std::size_t i = 0; i < this->rho.size(); i++) {
		this->rho[i] *= (-n_2) / denom;
		//std::cout << rho[i] << "\t" << i << "\n";
	}
	
}

void Particles::SETV(Field& fieldE, FieldB& B) {
	for (int i = x.size() - 1; i >= 0; i--) {

		//расчет скорости на dt/2
		this->v_x[i] += fieldE.field_by_x(x[i]) * (-this->q/ 2.) / this->m * (dt / 2.0);


		//расчет поворота частицы
		double temp_t = this->q * B.b_by_x(x[i]) / this->m * dt / 2;
		double temp_s = 2 * temp_t / (1 + temp_t * temp_t);
		double temp_c = (1 - temp_t * temp_t) / (1 + temp_t * temp_t);

		double v_x_new = this->v_x[i] * temp_c + this->v_y[i] * temp_s;
		this->v_y[i] = -this->v_x[i] * temp_s + this->v_y[i] * temp_c;
		this->v_x[i] = v_x_new;

		//еще подвинулись на dt/2
		this->v_x[i] += fieldE.field_by_x(x[i]) * (-this->q/2.0) / this->m * (dt / 2.0);
		this->x[i] += this->v_x[i] * dt;



		/*
		// ГУ вылет и все
		if (x[i] >= L || x[i] < 0) {
			swap_and_delete(i);
			//i--;
		}
		*/


		// граничные условия: жесткая стенка
		if (this->x[i] >= L - eps) {
			this->x[i] = (L - eps) - (this->x[i] - (L - eps));
			this->v_x[i] *= (-1);
		}
		if (this->x[i] - eps <= 0) {
			this->x[i] *= (-1);
			this->v_x[i] *= (-1);
		}



	}
}

double Particles::get_random() {
	static std::mt19937_64 gen(seed); // Фиксированный seed
	static std::uniform_real_distribution<double> dist(0.0, 1.0);

	return dist(gen);
}



void Particles::ionization_first() {
	x.push_back(L * get_random());

	R_s = get_random();
	R_theta = get_random();


	v_x.push_back(v_t * std::sqrt(-2 * std::log(R_s)) * std::cos(2 * 3.1416 * R_theta));
	v_y.push_back(v_t * std::sqrt(-2 * std::log(R_s)) * std::sin(2 * 3.1416 * R_theta));
};

std::vector<double> Particles::give_rho() {
	return rho;
}

void Particles::fill() {
	for (std::size_t i = 0; i < N; i++) {
		x[i] = L * get_random();

	}
	for (std::size_t i = 0; i < N; i++) {
		R_s = get_random();
		R_theta = get_random();

		v_x[i] = v_t * std::sqrt(-2 * std::log(R_s)) * std::cos(2 * 3.1416 * R_theta);
		v_y[i] = v_t * std::sqrt(-2 * std::log(R_s)) * std::sin(2 * 3.1416 * R_theta);
	}

} 