#include <iostream>
#include "AAfrag.h"
#include "kamae.h"
#include "kelner-aharonian.h"

int main() {
	try {
		Kamae kamae(Particle::neutrinos);
		KelnerAharonian kelner(Particle::neutrinos);
		AAfrag aafrag(Particle::neutrinos);

		for (double E_gamma = 1; E_gamma < 3e6; E_gamma *= 1.01) {
			std::cout << std::scientific;
			std::cout << E_gamma << " ";
			std::cout << kamae.get(1e1, E_gamma) << " ";
			std::cout << kamae.get(1e2, E_gamma) << " ";
			std::cout << kamae.get(1e3, E_gamma) << " ";
			std::cout << kamae.get(1e5, E_gamma) << " ";
			std::cout << kelner.get(1e1, E_gamma) << " ";
			std::cout << kelner.get(1e2, E_gamma) << " ";
			std::cout << kelner.get(1e3, E_gamma) << " ";
			std::cout << kelner.get(1e5, E_gamma) << " ";
			std::cout << aafrag.get(1e1, E_gamma) << " ";
			std::cout << aafrag.get(1e2, E_gamma) << " ";
			std::cout << aafrag.get(1e3, E_gamma) << " ";
			std::cout << aafrag.get(1e5, E_gamma) << " ";
			std::cout << "\n";
		}
	}
	catch (const std::exception& e) {
		std::cout << "\n" << e.what() << "\n";
	}
	return 0;
}
