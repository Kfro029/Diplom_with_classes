// test_cmake.cpp: определяет точку входа для приложения.
//

#include "particles.h"
#include "fields.h"
#include <random>
#include <fstream>
#include <random>
#include <cmath>
#include "neutrons.h"

#include <algorithm>


double el_x = 0;
void ionization(Particles& ions, Particles& el) {


	for (std::size_t i = 0; i < el.ionaze.size(); i++) {

		// целая часть ионизации в ячейке
		for (int j = 0; j < int(el.ionaze[i]); j++) {
			el_x = el.get_random();
			el.x.push_back(((i * 1.0) + el_x) * dx);
			ions.x.push_back(((i * 1.0) + el_x) * dx);


			el.R_s = el.get_random();
			el.R_theta = el.get_random();

			ions.R_s = ions.get_random();
			ions.R_theta = ions.get_random();


			el.v_x.push_back(el.v_t * std::sqrt(-2 * std::log(el.R_s)) * std::cos(2 * 3.1416 * el.R_theta));
			el.v_y.push_back(el.v_t * std::sqrt(-2 * std::log(el.R_s)) * std::sin(2 * 3.1416 * el.R_theta));


			ions.v_x.push_back(ions.v_t * std::sqrt(-2 * std::log(ions.R_s)) * std::cos(2 * 3.1416 * ions.R_theta));
			ions.v_y.push_back(ions.v_t * std::sqrt(-2 * std::log(ions.R_s)) * std::sin(2 * 3.1416 * ions.R_theta));
			
			
			//std::cout << "born ";
			//el.born++;
			//ions.born++;
		}
		//дробная
		if (el.get_random() < fmod(el.ionaze[i], 1.0)) {
			
			el_x = el.get_random();
			el.x.push_back((i + el_x) * dx);
			ions.x.push_back((i + el_x) * dx);


			el.R_s = el.get_random();
			el.R_theta = el.get_random();

			ions.R_s = ions.get_random();
			ions.R_theta = ions.get_random();


			el.v_x.push_back(el.v_t * std::sqrt(-2 * std::log(el.R_s)) * std::cos(2 * 3.1416 * el.R_theta));
			el.v_y.push_back(el.v_t * std::sqrt(-2 * std::log(el.R_s)) * std::sin(2 * 3.1416 * el.R_theta));


			ions.v_x.push_back(ions.v_t * std::sqrt(-2 * std::log(ions.R_s)) * std::cos(2 * 3.1416 * ions.R_theta));
			ions.v_y.push_back(ions.v_t * std::sqrt(-2 * std::log(ions.R_s)) * std::sin(2 * 3.1416 * ions.R_theta));
			
			
			//std::cout << "born ";
			//el.born++;
			//ions.born++;
		}


	}

}



int main() {

	Field fields;
	fields.loadFromFile_B("B.txt");
	fields.loadFromFile_Ey("Ey.txt");


	//Создаются электроны и ионы
	Particles electrons(m_el, -q, v_t_el, 1);
	Particles ions(m_ion, q, v_t_ions, 2);

	//Neutrons neutrons(m_ion, std::sqrt(3 * kB * 400 / m_ion), 3);
	
	//заполняются координаты и скорости соответственно заданным функциям
	ions.fill(fields);
	electrons.fill(fields);

	//neutrons.fill();

	//ions.loadFromFile("ionization_ion_values.txt");
	electrons.loadFromFile_ionization("ionization.txt");



	//std::cout << (ions.ionaze.size() == num_ceil) << std::endl;

	ions.fill_null_part();
	electrons.fill_null_part();

	electrons.CIC();
	ions.CIC();
	
	


	fields.solve_field(electrons, ions);
		
	electrons.SETV(fields);
	ions.SETV(fields);



	std::ofstream ne1;
	ne1.open("n_ntr.txt");

	std::ofstream fi1;
	fi1.open("fi.txt");
	
	std::ofstream rho_ions1;
	rho_ions1.open("rho_ions.txt");

	std::ofstream rho_el1;
	rho_el1.open("rho_el.txt");

	std::ofstream out_el1;
	out_el1.open("out_el.txt");

	std::ofstream out_ions1;
	out_ions1.open("out_ions.txt");

	std::ofstream coll_el1;
	coll_el1.open("coll_el.txt");

	std::ofstream born_el1;
	born_el1.open("born_el.txt");

	std::ofstream born_ions1;
	born_ions1.open("born_ions.txt");


	



	//std::vector<double> rho_el = electrons.rho;
	//std::vector<double> rho_ions = ions.rho;
	for (std::size_t p = 0; p < electrons.ionaze.size(); p++) {

		
		
		//rho_ions1 << ions.rho[p] / denom << " ";
		
		
		//std::cout << electrons.ionaze[p] << " ";
		
		
		//ne1 << neutrons.n[p] << " ";

	}

	
	
	for (std::size_t p = 0; p < num_ceil; p++) {

		rho_ions1 << ions.rho[p] / denom << " ";
		rho_el1 << electrons.rho[p] / denom << " ";
		//ne1 << neutrons.n[p] << " ";

	}
	rho_ions1 << std::endl;
	rho_el1 << std::endl;
	ne1 << std::endl;
	

	electrons.fill_null_part();
	ions.fill_null_part();

	electrons.CIC();
	ions.CIC();
	


	fields.solve_field(electrons, ions);

	std::vector<double> fi = fields.fi;

	
	//std::cout << "coll_want = "  << coll_want << std::endl;

	for (int i = 1; i <= (T / dt); i++) {
		electrons.move(fields);
		ions.move(fields);
		
		electrons.emission_ionization(electrons.out_an - ions.out_kat);

		//neutrons.move();

		//electrons.ionization_first();
		//std::cout << rho_el.x[0] << "\t" << i << "\n";

		electrons.fill_null_part();
		ions.fill_null_part();

		electrons.CIC();
		ions.CIC();

		//neutrons.concentration();


		//ПОТОМ ДОБАВИТЬ
		ionization(ions, electrons);

		//neutrons.ionization();



		fields.solve_field(electrons, ions);

		//fi = fields.fi;

		if (i % shot == 0) {

			
			/*
			std::ofstream coll_el1("coll_el.txt", std::ios::app);
			std::ofstream out_el1("out_el.txt", std::ios::app);
			std::ofstream out_ions1("out_ions.txt", std::ios::app);
			std::ofstream rho_ions1("rho_ions.txt", std::ios::app);
			std::ofstream rho_el1("rho_el.txt", std::ios::app);
			std::ofstream fi1("fi.txt", std::ios::app);
			std::ofstream ne1("n_ntr.txt", std::ios::app);
			std::ofstream born_el1("born_el.txt", std::ios::app);
			std::ofstream born_ions1("born_ions.txt", std::ios::app);
			*/
			std::cout << "writing i = " << i << std::endl;
			
			coll_el1 << electrons.coll << std::endl;
			
			for (std::size_t p = 0; p < num_ceil; p++) {
				out_el1 << electrons.divergention[p] * n_2 / (dt * shot) * q << " "; // запись тока в СИ (А/м^2)
				// out_el1 << electrons.divergention[p] << " "; Запись в штуках
				out_ions1 << ions.divergention[p] * n_2 / (dt * shot) * q << " ";


				//rho_ions1 << ions.rho[p] / denom / (q * (-n_2) * denom / dx) << " "; //в частицах
				//rho_el1 << electrons.rho[p] / denom / (q * (-n_2) * denom / dx) << " ";

				rho_ions1 << ions.rho[p] / denom << " "; // в Кл/м3
				rho_el1 << electrons.rho[p] / denom << " "; // в Кл/м3
				fi1 << fields.fi[p] << " ";
				//ne1 << neutrons.n[p] << " ";

			}
			born_el1 << electrons.born << std::endl;
			born_ions1 << ions.born << std::endl;

			out_el1 << std::endl;
			out_ions1 << std::endl;
			rho_ions1 << std::endl;
			rho_el1 << std::endl;
			fi1 << std::endl;
			ne1 << std::endl;
			
			ions.zeros();
			electrons.zeros();
			// break;
		}
		
	}
}
