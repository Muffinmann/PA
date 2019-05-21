/*
 * Materialclass_cristalplasticity_norton_creep.h
 *
 *  Created on: Jan 19, 2016
 *      Author: iwtm84
 */

#ifndef MATERIALCLASS_CRISTALPLASTICITY_NORTON_CREEP_H_
#define MATERIALCLASS_CRISTALPLASTICITY_NORTON_CREEP_H_

#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/lac/vector.h>
#include <deal.II/base/tensor.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <deal.II/lac/householder.h>
#include <deal.II/base/parameter_handler.h>

namespace Material {
using namespace dealii;
using namespace std;

template<int dim>
SymmetricTensor<4, dim> stress_strain_cubic_symmetric(double C11, double C12,
		double C44) {
	SymmetricTensor<4, dim> tmp;

	for (int i = 0; i < dim; i++) {
		for (int j = 0; j < dim; j++) {
			for (int k = 0; k < dim; k++) {
				for (int l = 0; l < dim; l++) {
					for (int a = 0; a < dim; a++) {
						if (a == i && a == j && a == k && a == l) {
							tmp[i][j][k][l] = C11;
						}
						for (int b = a + 1; b < dim; b++) {
							if (a == i && a == j && b == k && b == l && a < b) {
								tmp[i][j][k][l] = C12;
							}
							if (b == i && b == j && a == k && a == l && a < b) {
								tmp[i][j][k][l] = C12;
							}
							if (a == i && b == j && a == k && b == l && a < b) {
								tmp[i][j][k][l] = C44;
							}
						}
					}
				}
			}
		}
	}
	return tmp;
}

template<int dim>
Tensor<2, dim> rotmat(double angle, int normal) {
	Assert(dim==3, ExcMessage("rotmat only for dim=3 implemented"));
	Assert(normal==1 || normal==2,
			ExcMessage("rotmat only for normal 1 and 2 implemented"));
	angle = angle / 180 * 3.1415926535;
	Tensor<2, dim> rotmat;
	if (normal == 0) {
		rotmat[TableIndices<2>(0, 0)] = 1;
		rotmat[TableIndices<2>(1, 1)] = cos(angle);
		rotmat[TableIndices<2>(2, 2)] = cos(angle);
		rotmat[TableIndices<2>(1, 2)] = sin(angle);
		rotmat[TableIndices<2>(2, 1)] = -sin(angle);
	} else if (normal == 1) {
		rotmat[TableIndices<2>(0, 0)] = cos(angle);
		rotmat[TableIndices<2>(1, 1)] = 1;
		rotmat[TableIndices<2>(2, 2)] = cos(angle);
		rotmat[TableIndices<2>(0, 2)] = -sin(angle);
		rotmat[TableIndices<2>(2, 0)] = sin(angle);
	} else if (normal == 2) {
		rotmat[TableIndices<2>(0, 0)] = cos(angle);
		rotmat[TableIndices<2>(1, 1)] = cos(angle);
		rotmat[TableIndices<2>(2, 2)] = 1;
		rotmat[TableIndices<2>(0, 1)] = sin(angle);
		rotmat[TableIndices<2>(1, 0)] = -sin(angle);
	}
	return rotmat;
}

template<int dim>
std::vector<Tensor<1, dim>> load_s() {
	//load s
	std::vector<Tensor<1, dim>> s(24);
	{
		ifstream in("Input/slipdirections_fcc");
		if (in) {
			int i = 0, j = 0;
			string line;
			while (getline(in, line)) {
				stringstream sep(line);
				string field;
				j = 0;
				Assert(i<s.size(),
						ExcMessage("Too many slipdirections in slipdirections_fcc."));
				while (getline(sep, field, ',')) {
					Assert(j<dim,
							ExcMessage("Too many vectorcomponents in slipdirections_fcc."));
					s[i][j++] = stod(field);
				}
				i++;
			}
		} else {
			Assert(false, ExcMessage("File slipdirections_fcc not found."));
		}
	}
	return s;
}

template<int dim>
std::vector<Tensor<1, dim>> load_m() {
	//load m
	std::vector<Tensor<1, dim>> m(24);
	{
		ifstream in("Input/slipplanenormals_fcc");
		if (in) {
			int i = 0, j = 0;
			string line;
			while (getline(in, line)) {
				stringstream sep(line);
				string field;
				j = 0;
				Assert(i<m.size(),
						ExcMessage("Too many normals in slipplanenormals_fcc."));
				while (getline(sep, field, ',')) {
					Assert(j<dim,
							ExcMessage("Too many normals in slipplanenormals_fcc."));
					m[i][j++] = stod(field);
				}
				i++;
			}
		} else {
			Assert(false, ExcMessage("File slipplanenormals_fcc not found."));
		}
	}
	return m;
}

template<int dim>
class Material {
public:
	Material();
	~Material();

	SymmetricTensor<2, dim> get_old_stress() {
		return stress;
	}
	SymmetricTensor<2, dim> get_plastic_strain() {
		return old_plastic_strain;
	}
	double get_equivalent_plastic_strain() {
		return old_plastic_strain.norm();
	}
	int get_n_active_systems() {
		return n_active_systems_last_step;
	}
	int get_n_changed_active_systems() {
		return n_changed_active_systems_last_step;
	}
	int get_niterations() {
		return niterations;
	}

	SymmetricTensor<4, dim> calculateEPTangent(SymmetricTensor<2, dim> strain,
			double dt);
	void update_pointhistorydata();
	void init(ParameterHandler &parameter, double Euler1, double Euler2,
			double Euler3);
	Vector<double> get_direction(int dir);

private:
	void load_pointhistorydata();
	void init_P();
	void init_C_e();

	static const SymmetricTensor<2, dim> I;
	static const std::vector<Tensor<1, dim>> m;
	static const std::vector<Tensor<1, dim>> s;
	static const SymmetricTensor<4, dim> IxI;
	static const SymmetricTensor<4, dim> deviatoric_identity;

	double p; //exponent of evolutionequation
//	double kappa;
//	double mu;
	double eta;
	double C11;
	double C12;
	double C44;
	double h0;
	double tau0;
	double tauS;
	double hardeningtype;
	SymmetricTensor<4, dim> stress_strain_cubic;

	double eulerangle1;
	double eulerangle2;
	double eulerangle3;
	vector<SymmetricTensor<2, dim>> P;

	//Historyvariables
	SymmetricTensor<2, dim> stress;
	SymmetricTensor<2, dim> old_plastic_strain;
	SymmetricTensor<2, dim> updated_plastic_strain;
	Vector<double> old_g;
	Vector<double> updated_g;
	Vector<double> dlambda_old;
	int n_active_systems_last_step;
	int n_changed_active_systems_last_step;
	double old_A; //Average of hardening
	double updated_A;
	int niterations;		//to store niterations
	int maxNRIterations;	//read form parameterfile
	double maxErrorNR;		//read form parameterfile
};
//------------------------------------------------
//Intitialize Classvariables

template<int dim>
const SymmetricTensor<2, dim> Material<dim>::I = unit_symmetric_tensor<dim>();

template<int dim>
const SymmetricTensor<4, dim> Material<dim>::IxI = outer_product(I, I);

template<int dim>
const SymmetricTensor<4, dim> Material<dim>::deviatoric_identity =
		deviator_tensor<dim>();

template<int dim>
const std::vector<Tensor<1, dim>> Material<dim>::s = load_s<dim>();

template<int dim>
const std::vector<Tensor<1, dim>> Material<dim>::m = load_m<dim>();

//------------------------------------------------
//Constructors, Destrucctors
template<int dim> Material<dim>::Material() {
	//Materialparameter todo define with call of constructor
	//lambda:1125, mu:562.5, E:1500, nu:1/3, kappa:1500
//	kappa = 1500 * 200;
//	mu = 562.5 * 200;

	//http://solidmechanics.org/text/Chapter3_2/Chapter3_2.htm
	//in MPa (N/mm^2), Daten von Nickel

	//Pointhistory
	stress = 0;
	old_plastic_strain = 0;
	updated_plastic_strain = 0;
	old_g.reinit(24);
	updated_g.reinit(24);
	old_A = 0;
	updated_A = 0;
	dlambda_old.reinit(24);

	niterations = 0;
	n_active_systems_last_step = 0;
	n_changed_active_systems_last_step = 0;

}

template<int dim> Material<dim>::~Material() {
}
//------------------------------------------------

template<int dim>
SymmetricTensor<4, dim> Material<dim>::calculateEPTangent(
		SymmetricTensor<2, dim> strain, double dt) {
	niterations = 0; //to store the maximum number of inner NR Iterations
	SymmetricTensor<4, dim> C_el;
	C_el = stress_strain_cubic;
	int ndof = P.size();

	Vector<double> delta_lambda;
	delta_lambda.reinit(24);
	delta_lambda = dlambda_old;

	Vector<double> tau;

	Vector<double> residuum;

	Householder<double> solver;
	int iterationcounter = 0;
	double first_residuum = 0;
	do {

		load_pointhistorydata(); //copy old historydata into updatevariables
		for (int i = 0; i < ndof; i++) {
			updated_A += delta_lambda(i);
			updated_plastic_strain += delta_lambda(i) * P[i];
//			updated_g(i) += h0 * delta_lambda(i);
		}
		// sinh and cosh terms for hardening law
		double cosh_term = (exp(h0 * updated_A / (tauS - tau0))
				+ exp(-h0 * updated_A / (tauS - tau0))) / 2;
		double sinh_term = (exp(h0 * updated_A / (tauS - tau0))
				- exp(-h0 * updated_A / (tauS - tau0))) / 2;
		for (int i = 0; i < ndof; i++) {
			for (int j = 0; j < ndof; j++) {
				updated_g(i) += h0 / pow(cosh_term, 2)
						* (hardeningtype
								+ (1 - hardeningtype) * (i == j ? 1 : 0))
										* delta_lambda(j);
			}
		}
		stress = C_el * (strain - updated_plastic_strain);
		tau.reinit(ndof);
		for (int i = 0; i < tau.size(); i++) {
			tau(i) = stress * P[i];
			tau(i) = max(0.0, tau(i));
		}
		residuum.reinit(ndof);
		for (int i = 0; i < residuum.size(); i++) {
			// Attention -res is calculated
			residuum(i) = -delta_lambda(i)
					+ dt / eta * pow((tau(i) / updated_g(i)), p);
		}

		if (iterationcounter == 0) {
			first_residuum = residuum.l2_norm();
			if (first_residuum < 1e-5) {
				if (first_residuum < maxErrorNR) {
					return C_el;
				}
				first_residuum = 1e-5;
			}
		}
		if (residuum.l2_norm() < maxErrorNR * first_residuum) {
			niterations = iterationcounter;
//			if (iterationcounter == 0) {
//				return C_el;
//			}
			break;
		}
		if (iterationcounter++ >= maxNRIterations) {
			cout << "Iterationscounter inner NewtonRaphson >= "
					<< maxNRIterations << endl;
			cout << "First Residuum: " << first_residuum << " actual Residuum: "
					<< residuum.l2_norm() << endl;
			cout << "Strain: " << strain << endl;
			cout << "plastic Strain: " << old_plastic_strain << endl;
			cout << "Abs(Strain - plastic Strain): "
					<< (strain - old_plastic_strain).norm() << endl;
			cout << "Stress: " << stress << endl;
			throw bad_exception();
		}
		FullMatrix<double> J(ndof, ndof);
		for (int i = 0; i < ndof; i++) {
			double latent_hardening = 0;
			for (int gamma = 0; gamma < ndof; gamma++) {
				latent_hardening += delta_lambda(gamma)
						* (hardeningtype
								+ (1 - hardeningtype) * (i == gamma ? 1 : 0))
						* 2 * h0 * h0 / (tauS - tau0)
						* (-sinh_term / pow(cosh_term, 3));
			}
			double h_A = h0 / pow(cosh_term, 2);
			J(i, i) = 1;
			if (tau(i) > 0) {
				for (int j = 0; j < ndof; j++) {
					double h_ij = h_A
							* (hardeningtype
									+ (1 - hardeningtype) * (i == j ? 1 : 0));

					J(i, j) += -p * dt / eta * pow(tau(i) / updated_g(i), p - 1)
							/ pow(updated_g(i), 2)
							* (-P[i] * C_el * P[j] * updated_g(i)
									- tau(i) * (h_ij + latent_hardening));
//					J(i, j) +=
//							p * dt / eta * pow(tau(i) / updated_g(i), p - 1)
//									* (P[i] * C_el * P[j] * updated_g(i)
//											+ tau(i) * h_ij)
//									/ pow(updated_g(i), 2);
				}
			}
		}

		//solve J x ddlambda = -res
		solver.initialize(J);
		Vector<double> ddlambda;
		ddlambda.reinit(ndof);
		solver.least_squares(ddlambda, residuum);

		delta_lambda += ddlambda;
		if (ddlambda.l2_norm() < maxErrorNR*maxErrorNR){
			deallog << "ddlambda very small" << endl;
		}
//		cout << "Residuum of iteration " << iterationcounter << ": "
//				<< residuum.l2_norm() << endl;

	} while (true);

	n_active_systems_last_step = 0;
	for (int i = 0; i < delta_lambda.size(); i++) {
		if (abs(delta_lambda(i)) > 1e-10) {
			n_active_systems_last_step++;
		}
	}
	n_changed_active_systems_last_step = 0;
	for (int i = 0; i < delta_lambda.size(); i++) {
		if ((abs(delta_lambda(i)) > 1e-10) != (abs(dlambda_old(i)) > 1e-10)) {
			n_changed_active_systems_last_step++;
		}
	}
	dlambda_old = delta_lambda;

	// assemble rhs partial Ra / partial eps_n+1
	// -d_R/d_eps = d_R/d_dlambda x d_dlambda/d_eps
	// Stored in voigt like format
	// Attention -d_R/d_eps is calculated
	vector<Vector<double>> dr_deps(6, Vector<double>(ndof));

	for (int i = 0; i < ndof; i++) {
		SymmetricTensor<2, dim> dRi_deps;
		dRi_deps = dt * p / eta / updated_g(i)
				* pow(tau(i) / updated_g(i), p - 1) * C_el * P[i];
		dr_deps[0](i) = dRi_deps(TableIndices<2>(0, 0));
		dr_deps[1](i) = dRi_deps(TableIndices<2>(1, 1));
		dr_deps[2](i) = dRi_deps(TableIndices<2>(2, 2));
		dr_deps[3](i) = dRi_deps(TableIndices<2>(0, 1));
		dr_deps[4](i) = dRi_deps(TableIndices<2>(0, 2));
		dr_deps[5](i) = dRi_deps(TableIndices<2>(1, 2));
	}	// END ndof-loop

	// calculate d_dlambda / d_eps solving
	// -d_R/d_eps = d_R/d_dlambda x d_dlambda/d_eps
	Vector<double> ddlambda_deps11;
	Vector<double> ddlambda_deps22;
	Vector<double> ddlambda_deps33;
	Vector<double> ddlambda_deps12;
	Vector<double> ddlambda_deps13;
	Vector<double> ddlambda_deps23;
	ddlambda_deps11.reinit(ndof);
	ddlambda_deps22.reinit(ndof);
	ddlambda_deps33.reinit(ndof);
	ddlambda_deps12.reinit(ndof);
	ddlambda_deps13.reinit(ndof);
	ddlambda_deps23.reinit(ndof);
	solver.least_squares(ddlambda_deps11, dr_deps[0]);
	solver.least_squares(ddlambda_deps22, dr_deps[1]);
	solver.least_squares(ddlambda_deps33, dr_deps[2]);
	solver.least_squares(ddlambda_deps12, dr_deps[3]);
	solver.least_squares(ddlambda_deps13, dr_deps[4]);
	solver.least_squares(ddlambda_deps23, dr_deps[5]);
	SymmetricTensor<4, dim> C_pl;
	for (int i = 0; i < ndof; i++) {
		SymmetricTensor<2, dim> ddlambda_deps;
		ddlambda_deps(TableIndices<2>(0, 0)) = ddlambda_deps11(i);
		ddlambda_deps(TableIndices<2>(1, 1)) = ddlambda_deps22(i);
		ddlambda_deps(TableIndices<2>(2, 2)) = ddlambda_deps33(i);
		ddlambda_deps(TableIndices<2>(0, 1)) = ddlambda_deps12(i);
		ddlambda_deps(TableIndices<2>(0, 2)) = ddlambda_deps13(i);
		ddlambda_deps(TableIndices<2>(1, 2)) = ddlambda_deps23(i);
		C_pl += outer_product(C_el * P[i], ddlambda_deps);
	}

	return C_el - C_pl;
}

template<int dim> void Material<dim>::update_pointhistorydata() {
	old_g = updated_g;
	old_plastic_strain = updated_plastic_strain;
	old_A = updated_A;
}

template<int dim> void Material<dim>::load_pointhistorydata() {
	updated_g = old_g;
	updated_plastic_strain = old_plastic_strain;
	updated_A = old_A;
}

template<int dim> void Material<dim>::init(ParameterHandler &parameter,
		double Euler1, double Euler2, double Euler3) {
	double E = parameter.get_double("EModul");
	double nu = parameter.get_double("nu");
	double mu = parameter.get_double("mu");
	C11 = E * (1 - nu) / (1 - nu - 2 * nu * nu);
	C12 = E * nu / (1 - nu - 2 * nu * nu);
	C44 = mu;
	h0 = parameter.get_double("h0");
	tau0 = parameter.get_double("tau0");
	old_g = tau0;
	tauS = parameter.get_double("tauS");
	eta = parameter.get_double("eta");
	hardeningtype = parameter.get_double("hardeningtype");

	p = parameter.get_double("p");

	maxNRIterations = parameter.get_integer("maxInnerIterations");
	maxErrorNR = parameter.get_double("innerTolNewton");

	eulerangle1 = Euler1;
	eulerangle2 = Euler2;
	eulerangle3 = Euler3;
	init_P();
	init_C_e();
}

template<int dim>
void Material<dim>::init_P() {
	Tensor<2, dim> rotationmatrix;
	rotationmatrix = rotmat<dim>(eulerangle3, 2) * rotmat<dim>(eulerangle2, 1)
			* rotmat<dim>(eulerangle1, 2);

	P.resize(24);
	for (int i = 0; i < P.size(); i++) {
		P[i] = 0;
		Tensor<2, dim> tmp;
		outer_product(tmp, rotationmatrix * s[i] / s[i].norm(),
				rotationmatrix * m[i] / m[i].norm());
//		cout << "Norm vor Sym:  " << tmp.norm() << endl;
		P[i] = 0.5 * (tmp + transpose(tmp));
//		cout << "Norm nach Sym: " << P[i].norm() << endl;
//		cout << i<<": "<<P[i] << endl;
	}
}

template<int dim>
Vector<double> Material<dim>::get_direction(int dir) {
	Assert(dir>=0 && dir < 3, ExcMessage("Direction wrong"));
	Tensor<2, dim> rotationmatrix;
	rotationmatrix = rotmat<dim>(eulerangle3, 2) * rotmat<dim>(eulerangle2, 1)
			* rotmat<dim>(eulerangle1, 2);
	Vector<double> tmp;
	tmp.reinit(dim);
	tmp(0) = rotationmatrix[TableIndices<2>(0, dir)];
	tmp(1) = rotationmatrix[TableIndices<2>(1, dir)];
	tmp(2) = rotationmatrix[TableIndices<2>(2, dir)];
	return tmp;
}

template<int dim> void Material<dim>::init_C_e() {
//	SymmetricTensor<4, dim> tmp = stress_strain_cubic_symmetric<dim>(
//			kappa + 4.0 / 3 * mu, kappa - 2.0 / 3 * mu, mu);
	SymmetricTensor<4, dim> tmp = stress_strain_cubic_symmetric<dim>(C11, C12,
			C44);
//	cout << eulerangle1 << "/" << eulerangle2 << "/" << eulerangle3 << endl;
	Tensor<2, dim> R;
	R = rotmat<dim>(eulerangle3, 2) * rotmat<dim>(eulerangle2, 1)
			* rotmat<dim>(eulerangle1, 2);
//	cout << "R:  "<< R << endl;
	for (int m = 0; m < dim; m++) {
		for (int n = m; n < dim; n++) {
			for (int o = 0; o < dim; o++) {
				for (int p = o; p < dim; p++) {
					for (int i = 0; i < dim; i++) {
						for (int j = 0; j < dim; j++) {
							for (int k = 0; k < dim; k++) {
								for (int l = 0; l < dim; l++) {
									stress_strain_cubic[m][n][o][p] += R[m][i]
											* R[n][j] * R[o][k] * R[p][l]
											* tmp[i][j][k][l];
//									if(abs(stress_strain_cubic[m][n][o][p]) < 1e-20){
//										stress_strain_cubic[m][n][o][p]=0;
//									}
								}
							}
						}
					}
				}
			}
		}
	}
//	cout << "Norm vor Drehung:  " << tmp.norm() << endl;
//	cout << "Norm nach Drehung: " << stress_strain_cubic.norm() << endl;
//	cout << stress_strain_cubic - tmp << endl;
}

} //End Namespace

#endif /* MATERIALCLASS_CRISTALPLASTICITY_NORTON_CREEP_H_ */
