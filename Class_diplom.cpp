// test_cmake.cpp: определяет точку входа для приложения.
//

#include "particles.h"
#include "fields.h"
#include <random>
#include <fstream>



int main() {

	//Надо будет сначала создать массивы, потом передавать в констукторы класса
	std::vector<double> X_ions(N, 0.0);
	std::vector<double> X_el(N, 0.0);


	std::vector<double> V_ions(N, 0.0);
	std::vector<double> V_el(N, L / (1000 * dt));



	for (std::size_t i = 0; i < X_ions.size(); i++) {
		X_ions[i] = (i * 1.0) / X_ions.size() * L;
		X_el[i] = (i * 1.0) / X_el.size() * L;
	}
	for (std::size_t i = 0; i < X_el.size(); i++) {
		X_el[i] += 0.0005 * L * sin(X_el[i] / L * 2 * 3.1415926536 * k);
	}


	Particles electrons(m_el, -q, X_el, V_el, V_el);
	Particles ions(m_ion, q, X_ions, V_ions, V_ions);

	ions.fill_null_part();
	electrons.fill_null_part();

	electrons.CIC();
	ions.CIC();
	
	Field E;
	E.solve_field(electrons, ions);
		
	electrons.SETV(E);
	ions.SETV(E);




	
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
	


	E.solve_field(electrons, ions);


	for (int i = 1; i <= (T / dt); i++) {
		electrons.move(E);
		ions.move(E);

		electrons.ionization_first();
		//std::cout << rho_el.x[0] << "\t" << i << "\n";

		electrons.fill_null_part();
		ions.fill_null_part();

		electrons.CIC();
		ions.CIC();

		rho_el = electrons.rho;
		rho_ions = ions.rho;

		E.solve_field(electrons, ions);

		
		for (std::size_t p = 0; p < rho_el.size(); p++) {

			rho_ions1 << rho_ions[p] << " ";
			rho_el1 << rho_el[p] << " ";

		}
		rho_ions1 << std::endl;
		rho_el1 << std::endl;
		// break;
		
		
	}
}
