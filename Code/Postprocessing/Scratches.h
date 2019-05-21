/*
 * Scratches.h
 *
 *  Created on: Apr 6, 2016
 *      Author: iwtm84
 */

#ifndef CODE_POSTPROCESSING_SCRATCHES_H_
#define CODE_POSTPROCESSING_SCRATCHES_H_

#include "Extractor.h"
#include <functional>

using namespace std;
using namespace dealii;

template<int dim, int nslip>
struct gpData {
	Vector<double> gamma; // e.g. to calculate plastic strain
	Vector<double> dGamma_actualStep; // e.g. to calculate pi
	Vector<double> current_gamma_abs; // for hardening
	vector<Tensor<1, dim> > gradGamma; // for xi

	SymmetricTensor<2, dim> strain; // to calculate stress and Schmid stress

	double temperature;
	double dt;

	PerCrystalData<dim, nslip> *perCrystalData;

	gpData(){
		gamma = Vector<double>(nslip);
		dGamma_actualStep = Vector<double>(nslip);
		current_gamma_abs = Vector<double>(nslip);
		gradGamma = vector<Tensor<1, dim>>(nslip);
		strain = SymmetricTensor<2, dim>();
		temperature = 0;
		dt = 0;
		perCrystalData = nullptr;
	};
};

template<int dim, int nslip>
struct ScalarScratch {
	string name; // Name of the postprocessed data

	// vector with one Vector per Crystal with solution to build the patches
	vector<Vector<double> > solution;

	// functionpointer to a function that return the solution for the patches
	std_cxx11::function<double(gpData<dim, nslip>&)> post_function;

	// old construtor with functionpointer
	ScalarScratch(string n, const unsigned int ndof_solution,
			const unsigned int nDomains,
			double (*function)(gpData<dim, nslip>&)) :
			name(n), solution(nDomains, Vector<double>(ndof_solution)), post_function(
					function) {
	}

	// new constructor with function from c++11
	ScalarScratch(string n, const unsigned int ndof_solution,
			const unsigned int nDomains,
			std_cxx11::function<double(gpData<dim, nslip>&)> function) :
			name(n), solution(nDomains, Vector<double>(ndof_solution)), post_function(
					function) {
	}

};

//template<int dim, int nslip>
//struct ScalarScratch {
//	string name; // Name of the postprocessed data
//
//	// vector with one Vector per Crystal with solution to build the patches
//	vector<Vector<double> > solution;
//
//	// functionpointer to a function that return the solution for the patches
//	std_cxx11::function<double(Material<dim, nslip>&)> post_function;
//
//	// old construtor with functionpointer
//	ScalarScratch(string n, const unsigned int ndof_solution,
//			const unsigned int nDomains,
//			double (*function)(Material<dim, nslip>&)) :
//			name(n), solution(nDomains, Vector<double>(ndof_solution)), post_function(
//					function) {
//	}
//
//	// new constructor with function from c++11
//	ScalarScratch(string n, const unsigned int ndof_solution,
//			const unsigned int nDomains,
//			std_cxx11::function<double(Material<dim, nslip>&)> function) :
//			name(n), solution(nDomains, Vector<double>(ndof_solution)), post_function(
//					function) {
//	}
//};

#endif /* CODE_POSTPROCESSING_SCRATCHES_H_ */
