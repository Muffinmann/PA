/*
 * localCP.h
 *
 *  Created on: Apr 12, 2018
 *      Author: iwtm84
 */

#ifndef CODE_MATERIAL_LOCALCP_H_
#define CODE_MATERIAL_LOCALCP_H_

#include "PerCrystalData.h"
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/lac/vector.h>

namespace localCP {

template<int dim>
struct localCPdata{
	SymmetricTensor<2,dim> stress;
	SymmetricTensor<4,dim> tangent;
};

/*
 * returns the latent hardening factor
 */
template<int dim, int nslip>
inline double hardeningTerm(int alpha, int beta,
		PerCrystalData<dim, nslip> *perCrystalData) {
	if (alpha == beta) {
		return perCrystalData->get_hardeningtype();
	} else {
		return 1;
	}
}

template<int dim, int nslip>
localCPdata<dim> getLocalCPdata(const SymmetricTensor<2, dim> strain,
		double dt, double /*temperature*/,
		PerCrystalData<dim, nslip> *perCrystalData){
	Assert(perCrystalData->isLocalFormulation(),
				ExcMessage("localCP namespace used but is non local formulation."));
	localCPdata<dim> retval;

		double maxErrorNR = 1e-8;
		int maxNRIterations = 50;
		int iterationcounter = 0;
		double first_residuum = 0;

		double K = perCrystalData->get_h0();
		double tau0 = perCrystalData->get_pi0();
		double eta = perCrystalData->get_gamma0dot();
		double p = 1 / perCrystalData->get_p();

		SymmetricTensor<4, dim> C_el = perCrystalData->get_Cel();

		vector<SymmetricTensor<2, dim>> P = perCrystalData->get_P();

		int ndof = P.size();

		Vector<double> delta_lambda;
		delta_lambda.reinit(ndof);
		delta_lambda = 0;

		Vector<double> tau;
		tau.reinit(ndof);

		Vector<double> residuum;
		residuum.reinit(ndof);

		Householder<double> solver;

		Vector<double> g;
		g.reinit(ndof);
		SymmetricTensor<2, dim> plastic_strain;
		SymmetricTensor<2, dim> stress;

		do {
			/*
			 * reset A and plastic_strain because only one loadstep is allowed
			 * values are always 0 in the beginning
			 */
			plastic_strain = 0;
			for (unsigned int i = 0; i < delta_lambda.size(); i++) {
				plastic_strain += delta_lambda(i) * P[i];
			}

			g = tau0;
			for (unsigned int i = 0; i < g.size(); i++) {
				g(i) += K * abs(delta_lambda(i));
			}
			stress = C_el * (strain - plastic_strain);

			for (unsigned int i = 0; i < tau.size(); i++) {
				tau(i) = stress * P[i];
			}

			residuum = 0;
			for (unsigned int i = 0; i < residuum.size(); i++) {
				// Attention -res is calculated
				residuum(i) = -delta_lambda(i)
						+ dt / eta * pow((abs(tau(i)) / g(i)), p) * sgn(tau(i));
			}

			if (iterationcounter == 0) {
				first_residuum = residuum.l2_norm();
				if (first_residuum < maxErrorNR) {
					retval.stress = stress;
					retval.tangent = C_el;
					return retval;
				}
				if (first_residuum < 1) {
					first_residuum = 1;
				}
			}

			if (residuum.l2_norm() < maxErrorNR * first_residuum) {
				break;
			}

			if (iterationcounter++ >= maxNRIterations) {
				deallog << "Locally not converging res: " << residuum.l2_norm()
						<< " dt: " << dt << endl << "strain: " << strain << endl;
				throw bad_exception();
			}

			FullMatrix<double> J(ndof, ndof);
			for (int i = 0; i < ndof; i++) {

				J(i, i) = 1;
				for (int j = 0; j < ndof; j++) {
					double DtauDlambda = -sgn(tau(i)) * (P[i] * C_el * P[j]);
					double DgDlambda = (i == j ? 1 : 0) * sgn(delta_lambda(j)) * K;

					J(i, j) += -p * dt / eta * pow(abs(tau(i)) / g(i), p - 1)
							* sgn(tau(i))
							* (DtauDlambda * g(i) - abs(tau(i)) * DgDlambda)
							/ pow(g(i), 2);
				}
			}

			//solve J x ddlambda = -res
			solver.initialize(J);
			Vector<double> ddlambda;
			ddlambda.reinit(ndof);
			solver.least_squares(ddlambda, residuum);

			delta_lambda += ddlambda;
			if (ddlambda.l2_norm() < maxErrorNR * maxErrorNR) {
				deallog << "ddlambda very small" << endl;
			}
			//		cout << "Residuum of iteration " << iterationcounter << ": "
			//				<< residuum.l2_norm() << endl;

		} while (true);

	//	deallog << "Local problem solved" << endl;

		// assemble rhs partial Ra / partial eps_n+1
		// -d_R/d_eps = d_R/d_dlambda x d_dlambda/d_eps
		// Stored in voigt like format
		// Attention -d_R/d_eps is calculated
		vector<Vector<double>> dr_deps(6, Vector<double>(ndof));

		for (int i = 0; i < ndof; i++) {
			SymmetricTensor<2, dim> dRi_deps;
			dRi_deps = dt * p / eta / g(i) * pow(abs(tau(i)) / g(i), p - 1)
					* sgn(tau(i)) * C_el * P[i];
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

		retval.stress = stress;
		retval.tangent = C_el - C_pl;
		return retval;
}

template<int dim, int nslip>
SymmetricTensor<2, dim> getStress(const SymmetricTensor<2, dim> strain,
		double dt, double temperature,
		PerCrystalData<dim, nslip> *perCrystalData) {

	localCPdata<dim> retval = getLocalCPdata(strain, dt, temperature, perCrystalData);

	return retval.stress;

}

template<int dim, int nslip> //checked 27.6.16
SymmetricTensor<4, dim> getDsigDeps(const SymmetricTensor<2, dim> strain,
		double dt, double temperature,
		PerCrystalData<dim, nslip> *perCrystalData) {



	localCPdata<dim> retval = getLocalCPdata(strain, dt, temperature, perCrystalData);

	return retval.tangent;

}

} // end namespace localCP

#endif /* CODE_MATERIAL_LOCALCP_H_ */
