#include "neutrons.h"

// –еализаци€ конструктора
Neutrons::Neutrons(double m, double v_t, int seed)
	: m(m), 
	n(num_ceil, 0.0),
	temp_t(0.0), temp_s(0.0), temp_c(0.0), v_t(v_t),
	gen(seed), dist(0.0, 1.0),
	out(0),
	ionaze(num_ceil, 0.0)
	
{
	x.resize(N_ntr);
	v_x.resize(N_ntr);
	v_y.resize(N_ntr);
	ionaze[0] = 1.;
	//std::cout << "q in constructor = " << q << std::endl;
}



void Neutrons::concentration() {
	std::fill(n.begin(), n.end(), 0.0);

	for (std::size_t i = 0; i < x.size(); i++) {
		x_ceil = x[i] / dx;
		x_loc = fmod(x[i], dx);

		n[x_ceil] += n_ntr * (dx - x_loc) / dx;
		n[x_ceil + 1] += n_ntr * x_loc / dx;
	}

};

void Neutrons::swap_and_delete(int number) {
	out++;

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


void Neutrons::move() {

	for (int i = x.size() - 1; i >= 0; i--) {
		/*
		//расчет скорости на dt/2
		//v_x[i] += ;


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

		*/
		x[i] += v_x[i] * dt;

		// √” вылет и все
		if (this->x[i] >= L || this->x[i] < 0) {
			swap_and_delete(i);
			continue;
			//i--;
		}

		/*
		// граничные услови€: жестка€ стенка
		if (x[i] >= L - eps) {
			x[i] = (L - eps) - (x[i] - (L - eps));
			v_x[i] *= (-1);
		}
		if (x[i] - eps <= 0) {
			x[i] *= (-1);
			v_x[i] *= (-1);
		}
		*/

		/*
		if (get_random() < if_collis(x[i])) {
			R_s = get_random();
			R_theta = get_random();

			v_x[i] = (v_t * std::sqrt(-2 * std::log(R_s)) * std::cos(2 * 3.1416 * R_theta));
			v_y[i] = (v_t * std::sqrt(-2 * std::log(R_s)) * std::sin(2 * 3.1416 * R_theta));
		}
		*/
	}
}


double Neutrons::get_random() {
	return dist(gen);
}


//по€вление одной частицы в случайном месте с ћаксвелловским распределением
void Neutrons::ionization_first() {
	x.push_back(L * get_random());

	R_s = get_random();
	R_theta = get_random();

	v_x.push_back(v_t * std::sqrt(-2 * std::log(R_s)) * std::cos(2 * 3.1416 * R_theta));
	v_y.push_back(v_t * std::sqrt(-2 * std::log(R_s)) * std::sin(2 * 3.1416 * R_theta));
};


// начальное заполнение частицами
void Neutrons::fill() {
	for (std::size_t i = 0; i < N_ntr; i++) {
		ionization_first();
	}
}


void Neutrons::ionization() {
	for (std::size_t i = 0; i + 1 < num_ceil; i++) {

		// цела€ часть ионизации в €чейке
		for (int j = 0; j < int(ionaze[i]); j++) {
			x.push_back(((i * 1.0) + get_random()) * dx);

			R_s = get_random();
			R_theta = get_random();

			v_x.push_back(v_t * std::sqrt(-2 * std::log(R_s)) * std::cos(2 * 3.1416 * R_theta));
			v_y.push_back(v_t * std::sqrt(-2 * std::log(R_s)) * std::sin(2 * 3.1416 * R_theta));
		}
		//дробна€
		if (get_random() < fmod(ionaze[i], 1.0)) {
			x.push_back((i + get_random()) * dx);

			R_s = get_random();
			R_theta = get_random();

			v_x.push_back(v_t * std::sqrt(-2 * std::log(R_s)) * std::cos(2 * 3.1416 * R_theta));
			v_y.push_back(v_t * std::sqrt(-2 * std::log(R_s)) * std::sin(2 * 3.1416 * R_theta));
		}

	}

}