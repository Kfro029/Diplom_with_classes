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


int main() {

	Field fields;
	fields.loadFromFile("B_values.txt");


	//Создаются электроны и ионы
	Particles electrons(m_el, -q, v_t_el, 1);
	Particles ions(m_ion, q, v_t_ions, 2);

	//Neutrons neutrons(m_ion, std::sqrt(3 * kB * 400 / m_ion), 3);
	
	//заполняются координаты и скорости соответственно заданным функциям
	ions.fill(fields);
	electrons.fill(fields);

	//neutrons.fill();

	ions.loadFromFile("ionization_ion_values.txt");
	electrons.loadFromFile("ionization_el_values.txt");

	//std::cout << (ions.ionaze.size() == num_ceil) << std::endl;

	ions.fill_null_part();
	electrons.fill_null_part();

	electrons.CIC();
	ions.CIC();
	
	


	fields.solve_field(electrons, ions);
		
	electrons.SETV(fields);
	ions.SETV(fields);

	electrons.emission_ionization(electrons.out_an - ions.out_kat);


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


	



	//std::vector<double> rho_el = electrons.rho;
	//std::vector<double> rho_ions = ions.rho;


	
	
	for (std::size_t p = 0; p < num_ceil; p++) {

		rho_ions1 << ions.rho[p] * denom << " ";
		rho_el1 << electrons.rho[p] * denom << " ";
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


	double coll_want = 0.;
	for (std::size_t i = 0; i < num_ceil; i++) {
		coll_want += electrons.collision[i];
	}
	double max_collis = 0;
	

	/*
	for (std::size_t i = 0; i < num_ceil; i++) {
		std::cout << electrons.collision[i] << " ";
	}
	*/
	std::cout << "max collis = " << max_collis << std::endl;
	 
	//std::cout << "coll_want = "  << coll_want << std::endl;

	for (int i = 1; i <= (T / dt); i++) {
		electrons.move(fields);
		ions.move(fields);
		
		
		//electrons.emission_ionization(electrons.out_an - ions.out_kat);

		//neutrons.move();

		//electrons.ionization_first();
		//std::cout << rho_el.x[0] << "\t" << i << "\n";

		electrons.fill_null_part();
		ions.fill_null_part();

		electrons.CIC();
		ions.CIC();

		//neutrons.concentration();


		//ПОТОМ ДОБАВИТЬ
		electrons.ionization();
		ions.ionization();
		
		//neutrons.ionization();



		fields.solve_field(electrons, ions);

		//fi = fields.fi;

		if (i % shot == 0) {

			std::ofstream coll_el1("coll_el.txt", std::ios::app);
			std::ofstream out_el1("out_el.txt", std::ios::app);
			std::ofstream out_ions1("out_ions.txt", std::ios::app);
			std::ofstream rho_ions1("rho_ions.txt", std::ios::app);
			std::ofstream rho_el1("rho_el.txt", std::ios::app);
			std::ofstream fi1("fi.txt", std::ios::app);
			std::ofstream ne1("n_ntr.txt", std::ios::app);

			std::cout << "writing i = " << i << std::endl;
			
			coll_el1 << electrons.coll << std::endl;
			
			for (std::size_t p = 0; p < num_ceil; p++) {
				out_el1 << electrons.divergention[p] * n_2 * dt * shot * q << " ";
				out_ions1 << ions.divergention[p] * n_2 * dt * shot * q << " ";
				rho_ions1 << ions.rho[p] * denom / dx << " ";
				rho_el1 << electrons.rho[p] * denom / dx << " ";
				fi1 << fields.fi[p] << " ";
				//ne1 << neutrons.n[p] << " ";

			}
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
