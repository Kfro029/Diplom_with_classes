// test_cmake.cpp: определяет точку входа для приложения.
//

#include "particles.h"
#include "fields.h"
#include <random>
#include <fstream>
#include <random>
#include <cmath>
#include "fieldB.h"







int main() {

	
	
	/*
	for (std::size_t i = 0; i < X_ions.size(); i++) {
		X_ions[i] = (i * 1.0) / X_ions.size() * L;
		X_el[i] = (i * 1.0) / X_el.size() * L;
	}
	for (std::size_t i = 0; i < X_el.size(); i++) {
		X_el[i] += 0.0005 * L * sin(X_el[i] / L * 2 * 3.1415926536 * k);
	}
	*/

	//Создаются электроны и ионы
	Particles electrons(m_el, -q, v_t_el, 1);
	Particles ions(m_ion, q, v_t_ions, 2);

	
	//B.loadFromFile("B_values.txt");


	//заполняются координаты и скорости соответственно заданным функциям
	ions.fill();
	electrons.fill();

	ions.fill_null_part();
	electrons.fill_null_part();

	electrons.CIC();
	ions.CIC();
	
	Field fields;
	fields.loadFromFile("B_values.txt");
	fields.solve_field(electrons, ions);
		
	electrons.SETV(fields);
	ions.SETV(fields);



	std::ofstream fi1;
	fi1.open("fi.txt");
	
	std::ofstream rho_ions1;
	rho_ions1.open("rho_ions.txt");

	std::ofstream rho_el1;
	rho_el1.open("rho_el.txt");
	



	std::vector<double> rho_el = electrons.rho;
	std::vector<double> rho_ions = ions.rho;


	
	
	for (std::size_t p = 0; p < rho_el.size(); p++) {

		rho_ions1 << rho_ions[p] << " ";
		rho_el1 << rho_el[p] << " ";

	}
	rho_ions1 << std::endl;
	rho_el1 << std::endl;
	

	electrons.fill_null_part();
	ions.fill_null_part();

	electrons.CIC();
	ions.CIC();
	


	fields.solve_field(electrons, ions);

	std::vector<double> fi = fields.fi;

	for (int i = 1; i <= (T / dt); i++) {
		electrons.move(fields);
		ions.move(fields);

		//electrons.ionization_first();
		//std::cout << rho_el.x[0] << "\t" << i << "\n";

		electrons.fill_null_part();
		ions.fill_null_part();

		electrons.CIC();
		ions.CIC();

		rho_el = electrons.rho;
		rho_ions = ions.rho;

		fields.solve_field(electrons, ions);

		fi = fields.fi;

		
		for (std::size_t p = 0; p < rho_el.size(); p++) {

			rho_ions1 << rho_ions[p] << " ";
			rho_el1 << rho_el[p] << " ";
			fi1 << fi[p] << " ";

		}
		rho_ions1 << std::endl;
		rho_el1 << std::endl;
		fi1 << std::endl;
		// break;
		
		
	}
}
