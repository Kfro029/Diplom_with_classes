#include "particles.h" 
#include "fields.h"

// Реализация конструктора
Particles::Particles(double m, double q, std::vector <double> want_x, std::vector <double> want_v_x, std::vector <double> want_v_y)
	: m(m), q(q), x(want_x), v_x(want_v_x), v_y(want_v_y), rho(num_ceil, 0.0)
{
	//std::cout << "q in constructor = " << q << std::endl;
}

// Реализация методов

//ее надо будет запускать после Move
void Particles::fill_null_part() {
	std::fill(rho.begin(), rho.end(), 0.0);
};



void Particles::move(Field& field) {
	
	for (std::size_t i = 0; i < x.size(); i++) {

		//расчет скорости на dt/2
		v_x[i] += field.field_by_x(x[i]) * q / m * (dt / 2.0);
		

		//расчет поворота частицы
		t = q * B_0 / m * dt / 2;
		s = 2 * t / (1 + t * t);
		c = (1 - t * t) / (1 + t * t);

		v_x_new = v_x[i] * c + v_y[i] * s;
		v_y[i] = -v_x[i] * s + v_y[i] * c;
		v_x[i] = v_x_new;

		//еще подвинулись на dt/2
		v_x[i] += field.field_by_x(x[i]) * q / m * (dt / 2.0);
		x[i] += v_x[i] * dt;
	

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

void Particles::SETV(Field& field) {
	for (std::size_t i = 0; i < x.size(); i++) {
		

		v_x[i] = field.field_by_x(x[i]) * (-q / 2.0) / m * dt;
		x[i] += v_x[i] * dt;
		// граничные условия: жесткая стенка
		if (x[i] >= L - eps) {
			x[i] = (L - eps) - (x[i] - (L - eps));
			v_x[i] *= (-1);
		}
		if (x[i] <= eps) {
			x[i] *= (-1);
			v_x[i] *= (-1);
		}

	}
	// fill_null_part();
	// field.fill_null_field();
}

std::vector<double> Particles::give_rho() {
	return rho;
}