#include "particles.h" 

// Реализация конструктора
Particles::Particles(double m, double q, double v_t, int seed)
	: m(m), q(q),
	rho(num_ceil, 0.0), 
	temp_t(0.0), temp_s(0.0), temp_c(0.0), v_t(v_t),
	gen(seed), dist(0.0, 1.0),
	out_kat(0),
	out_an(0),
	ionaze(num_ceil -1, 0.0),
	collision(num_ceil, 0.00),
	fly(0),
	divergention(num_ceil, 0),
	coll(0),
	theor_coll(0.),
	born(0)
{
	x.reserve(N);
	v_x.reserve(N);
	v_y.reserve(N);
	//std::cout << "q in constructor = " << q << std::endl;
}

// Реализация методов

//вероятность столкновения
double Particles::if_collis(double x_part) {
	x_ceil = x_part / dx;
	x_loc = fmod(x_part, dx);
	return (collision[x_ceil] * (dx - x_loc) + collision[x_ceil + 1] * x_loc) / dx;
}


//ее надо будет запускать после Move
void Particles::fill_null_part() {
	std::fill(rho.begin(), rho.end(), 0.0);
};

//одно из граничных условий
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


void Particles::zeros() {
	std::fill(divergention.begin(), divergention.end(), 0);
	coll = 0;
	theor_coll = 0.;
	born = 0;
}

// движение
void Particles::move(Field& fieldE) {
	out_kat = 0;
	out_an = 0;
	
	for (int i = x.size() - 1; i >= 0; i--) {
		x_ceil = x[i] / dx;
		//расчет скорости на dt/2
		v_x[i] += fieldE.field_Ex_by_x(x[i]) * q / m * (dt / 2.0);
		v_y[i] += fieldE.field_Ey_by_x(x[i]) * q / m * (dt / 2.0);


		//расчет поворота частицы
		double temp_t = q * fieldE.b_by_x(x[i]) / m * dt / 2.;
		double temp_s = 2. * temp_t / (1. + temp_t * temp_t);
		double temp_c = (1. - temp_t * temp_t) / (1. + temp_t * temp_t);

		double v_x_new = v_x[i] * temp_c + v_y[i] * temp_s;
		v_y[i] = -v_x[i] * temp_s + v_y[i] * temp_c;
		v_x[i] = v_x_new;

		//еще подвинулись на dt/2
		v_x[i] += fieldE.field_Ex_by_x(x[i]) * q / m * (dt / 2.0);
		v_y[i] += fieldE.field_Ey_by_x(x[i]) * q / m * (dt / 2.0);
		
		x[i] += v_x[i] * dt;
		

		
		
		// ГУ вылет и все
		//подсчет вылета для разных частиц
		if ((x[i] / dx) > x_ceil + 1) {
			//divergention[x_ceil]--;
			divergention[x_ceil + 1]++;
		}
		if ((x[i] / dx) < x_ceil) {
			//divergention[x_ceil]--;
			divergention[x_ceil]--;
		}

		if (x[i] < 0) {
			swap_and_delete(i);
			out_an++;
			continue;
			//i--;
		}

		if (x[i] >= L) {
			swap_and_delete(i);
			out_kat++;
			continue;
			//i--;
		}
		/*
		// граничные условия: жесткая стенка
		if (x[i] >= L - eps) {
			x[i] = (L - eps) - (x[i] - (L - eps));
			v_x[i] *= (-1);
		}
		if (x[i] - eps <= 0) {
			x[i] *= (-1);
			v_x[i] *= (-1);
		}
		*/
		

		//СТОЛКНОВЕНИЯ
		/*
		theor_coll += if_collis(x[i]);
		if (get_random() < if_collis(x[i])) {
			R_s = get_random();
			R_theta = get_random();

			v_x[i] = (v_t * std::sqrt(-2 * std::log(R_s)) * std::cos(2 * 3.1416 * R_theta));
			v_y[i] = (v_t * std::sqrt(-2 * std::log(R_s)) * std::sin(2 * 3.1416 * R_theta));

			coll++;
		}
		*/
		
		
	}
	

	// fill_null_part();
	// field.fill_null_field();
}

// подсчет плотности заряда
void Particles::CIC() {
	for (std::size_t i = 0; i < x.size(); i++) {


		int ceil_left = (x[i])  / dx;
		int ceil_right = ceil_left + 1;
		
		double x_loc = fmod(x[i], dx);
		rho[ceil_left] += (dx - x_loc) / dx;
		rho[ceil_right] += x_loc / dx;
		
	}

	for (std::size_t i = 0; i < rho.size(); i++) {
		rho[i] *= q * (-n_2) * denom / dx;
	}
}

// начальная подвижка
void Particles::SETV(Field& fieldE) {
	for (int i = x.size() - 1; i >= 0; i--) {
		x_ceil = x[i] / dx;

		//расчет скорости на dt/2
		v_x[i] += fieldE.field_Ex_by_x(x[i]) * (-q/ 2.) / m * (dt / 2.0);
		v_y[i] += fieldE.field_Ey_by_x(x[i]) * (-q / 2.) / m * (dt / 2.0);


		//расчет поворота частицы
		double temp_t = (-q/2.) * fieldE.b_by_x(x[i]) / m * dt / 2;
		double temp_s = 2 * temp_t / (1 + temp_t * temp_t);
		double temp_c = (1 - temp_t * temp_t) / (1 + temp_t * temp_t);

		double v_x_new = v_x[i] * temp_c + v_y[i] * temp_s;
		v_y[i] = -v_x[i] * temp_s + v_y[i] * temp_c;
		v_x[i] = v_x_new;

		//еще подвинулись на dt/2
		v_x[i] += fieldE.field_Ex_by_x(x[i]) * (-q / 2.) / m * (dt / 2.0);
		v_y[i] += fieldE.field_Ey_by_x(x[i]) * (-q / 2.) / m * (dt / 2.0);

		//Столкновения
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

//рандомное число
double Particles::get_random() {
	return dist(gen);
}


//появление одной частицы в случайном месте с Максвелловским распределением
void Particles::ionization_first() {
	x.push_back(L * get_random());

	R_s = get_random();
	R_theta = get_random();


	v_x.push_back(v_t * std::sqrt(-2 * std::log(R_s)) * std::cos(2 * 3.1416 * R_theta));
	v_y.push_back(v_t * std::sqrt(-2 * std::log(R_s)) * std::sin(2 * 3.1416 * R_theta));
};


//появление нескольких частиц в эмиссионной плоскости
void Particles::emission_ionization(int out) {
	for (int i = 0; i < out; i++) {
		x.push_back(emission);

		R_s = get_random();
		R_theta = get_random();


		v_x.push_back(v_t * std::sqrt(-2 * std::log(R_s)) * std::cos(2 * 3.1416 * R_theta));
		v_y.push_back(v_t * std::sqrt(-2 * std::log(R_s)) * std::sin(2 * 3.1416 * R_theta));
	}
};

// начальное заполнение частицами и массива столкновения
void Particles::fill(Field& fieldE) {
	for (std::size_t i = 0; i < N; i++) {
		ionization_first();
	}

	//заполнение массива столкновений
	for (std::size_t i = 0; i < num_ceil; i++) {
		if (fieldE.B[i] > BPHI0) {
			collision[i] = abs(q) * fieldE.B[i] / m * (1. - 0.997 * exp(-((fieldE.B[i] - BPHI0) * (fieldE.B[i] - BPHI0) / BPHIWDTH / 7.) * ((fieldE.B[i] - BPHI0) * (fieldE.B[i] - BPHI0) / BPHIWDTH / 7.)));
		}
		else {
			collision[i] = abs(q) * fieldE.B[i] / m * (1. - 0.997 * exp(-((fieldE.B[i] - BPHI0) * (fieldE.B[i] - BPHI0) / BPHIWDTH) * ((fieldE.B[i] - BPHI0) * (fieldE.B[i] - BPHI0) / BPHIWDTH)));
		}
	}

	for (std::size_t i = 0; i < num_ceil; i++) {
		collision[i] = 1 - exp(-collision[i] * dt);
	}

}

// ионизация
void Particles::ionization() {

	
	for (std::size_t i = 0; i + 1 < num_ceil; i++) {

		// целая часть ионизации в ячейке
		for (int j = 0; j < int(ionaze[i]); j++) {
			x.push_back(((i * 1.0) + get_random()) * dx);

			R_s = get_random();
			R_theta = get_random();

			v_x.push_back(v_t * std::sqrt(-2 * std::log(R_s)) * std::cos(2 * 3.1416 * R_theta));
			v_y.push_back(v_t * std::sqrt(-2 * std::log(R_s)) * std::sin(2 * 3.1416 * R_theta));

			std::cout << "born ";
			born++;
		}
		//дробная
		if (get_random() < fmod(ionaze[i], 1.0)) {
			x.push_back((i  + get_random()) * dx);

			R_s = get_random();
			R_theta = get_random();

			v_x.push_back(v_t * std::sqrt(-2 * std::log(R_s)) * std::cos(2 * 3.1416 * R_theta));
			v_y.push_back(v_t * std::sqrt(-2 * std::log(R_s)) * std::sin(2 * 3.1416 * R_theta));
			
			std::cout << "born ";
			born++;
		}


	}

}

// подгрузка функции ионизации из файла
void Particles::loadFromFile_ionization(const std::string filename) {
	std::ifstream file(filename);

	if (!file.is_open()) {
		std::cerr << "Не удалось открыть файл: " << filename << std::endl;
		return;
	}

	for (std::size_t i = 0; i < ionaze.size() && file >> ionaze[i]; ++i);

	for (std::size_t i = 0; i < ionaze.size(); ++i) {
		ionaze[i] = ionaze[i] / n_2 / std::abs(q) * dt * dx;
	}
}
