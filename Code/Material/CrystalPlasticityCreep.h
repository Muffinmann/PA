/*
 * CrystalPlasticityCreep.h
 *
 *  Created on: Apr 1, 2016
 *      Author: iwtm84
 */

#ifndef CRYSTALPLASTICITYCREEP_H_
#define CRYSTALPLASTICITYCREEP_H_

using namespace dealii;
using namespace std;

#include "PerCrystalData.h"
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
#include <boost/lexical_cast.hpp>
#include "../Various/Tensor8.h"
#include "../Various/AuxilaryFunctions.h"
#include "localCP.h"

namespace GradientCP {

//----------------------------------
//---------Supportfunctions---------
//----------------------------------

/*
 * returns the latent hardening factor
 */
template<int dim, int nslip>
inline double hardeningTerm(int alpha, int beta,
		PerCrystalData<dim, nslip> *perCrystalData) {
	Assert(!perCrystalData->isLocalFormulation(),
			ExcMessage("gradientCP namespace used but is local formulation."));
	//1.58416e-06
	if (alpha == beta) {
		return perCrystalData->get_hardeningtype();
	} else {
		return 1;
	}

	//1.99888e-06
//	return perCrystalData->get_hardeningtype()
//						+ (1 - perCrystalData->get_hardeningtype())
//								* (alpha == beta ? 1 : 0);

//2.30975e-07
//	return perCrystalData->get_hardeningtype() + alpha == beta ?
//			(1 - perCrystalData->get_hardeningtype()) : 0;
}

/*
 * returns the plStrain due to the slips in gamma
 */
template<int dim, int nslip>
SymmetricTensor<2, dim> get_plStrain(const Vector<double> &gamma,
		PerCrystalData<dim, nslip> *perCrystalData) {
	Assert(!perCrystalData->isLocalFormulation(),
			ExcMessage("gradientCP namespace used but is local formulation."));
	SymmetricTensor<2, dim> plStrain;
	for (unsigned int i = 0; i < gamma.size(); ++i) {
		plStrain += gamma(i) * perCrystalData->get_P(i);
	}
	return plStrain;
}

/*
 * returns the thermal strain to a reference temperature of 23 degree celsius
 */
template<int dim, int nslip>
SymmetricTensor<2, dim> get_thermalStrain(double temperature,
		PerCrystalData<dim, nslip> *perCrystalData) {
	Assert(!perCrystalData->isLocalFormulation(),
			ExcMessage("gradientCP namespace used but is local formulation."));
	double referenceTemperature = 23;
	SymmetricTensor<2, dim> thermalStrain = perCrystalData->get_alpha()
			* (temperature - referenceTemperature)
			* unit_symmetric_tensor<dim>();
	return thermalStrain;
}

/*
 * return the cummulated hardening value
 */
double getA(const Vector<double> &current_gamma_abs) {
	double retVal = 0;
	for (unsigned int i = 0; i < current_gamma_abs.size(); i++) {
		Assert(current_gamma_abs(i)>=0, ExcMessage("Negative hardening value"));
		retVal += current_gamma_abs(i);
	}
	return retVal;
}

/*
 * return the temperature dependent slip resistance
 */
template<int dim, int nslip>
double getPiTemp(const double temperature,
		PerCrystalData<dim, nslip> *perCrystalData) {
	Assert(!perCrystalData->isLocalFormulation(),
			ExcMessage("gradientCP namespace used but is local formulation."));
	double t_0_stiffness = 3500;
	double roomTemperature = 23;
	return perCrystalData->get_pi0() * (t_0_stiffness - temperature)
			/ (t_0_stiffness - roomTemperature);
}

/*
 * returns the hardening for system alpha
 */
template<int dim, int nslip>
double getGAlpha(const int alpha, double temperature,
		const Vector<double> &dGamma, const Vector<double> &current_gamma_abs,
		const Vector<double> &prior_hardening,
		PerCrystalData<dim, nslip> *perCrystalData) {
	Assert(!perCrystalData->isLocalFormulation(),
			ExcMessage("gradientCP namespace used but is local formulation."));
	double retVal = prior_hardening(alpha); // hardening of last converged step

	switch (perCrystalData->get_regularisation()) {
	case tanh_lin: { // calculate linear hardening
		for (unsigned int beta = 0; beta < dGamma.size(); beta++) {
			retVal += perCrystalData->get_h0()
					* hardeningTerm(alpha, beta, perCrystalData)
					* abs(dGamma(beta));
		}
		break;
	}
	case powerlaw_sech: { //calculate sech hardening
		double updated_A = getA(current_gamma_abs);
		double piTemp = getPiTemp(temperature, perCrystalData);

		double h0A_taus0 = perCrystalData->get_h0() * updated_A
				/ (perCrystalData->get_piS() - piTemp);
		double cosh_term = (exp(h0A_taus0) + exp(-h0A_taus0)) / 2;
		double sechPow2 = 1.0 / pow(cosh_term, 2);

		for (unsigned int beta = 0; beta < dGamma.size(); beta++) {
			retVal += perCrystalData->get_h0() * sechPow2
					* hardeningTerm(alpha, beta, perCrystalData)
					* abs(dGamma(beta));
		}
		break;
	}//
	case bardella_lin: { 
		for (unsigned int beta = 0; beta < dGamma.size(); beta++) {
			retVal += perCrystalData->get_h0()
					* hardeningTerm(alpha,beta,perCrystalData)
					* abs(dGamma(beta)) //TODO
					* perCrystalData->get_Nh()
					* pow(dGamma(beta),perCrystalData->get_Nh()-1); //Nh = 1 here
		}
		break;
	}//
	case unknown: {
		Assert(false, ExcMessage("Unknown regularisation"))
		break;
	}
	};
//	if (retVal != retVal) {
//		deallog << endl << "Nan in getDAlpha" << endl << "alpha: " << alpha
//				<< " temperature: " << temperature << " dGamma: " << dGamma
//				<< " current_gamma_abs: " << current_gamma_abs
//				<< " prior_hardening: " << prior_hardening << endl;
//	}
	return retVal;
}

//----------------------------------
//---------------RHS----------------
//----------------------------------
/*
 * returns the macroscopic stress due to elastic strains
 */
template<int dim, int nslip>
SymmetricTensor<2, dim> get_stress(const SymmetricTensor<2, dim> strain,
		const Vector<double> &gamma, double temperature,
		PerCrystalData<dim, nslip> *perCrystalData) {
	Assert(!perCrystalData->isLocalFormulation(),
			ExcMessage("gradientCP namespace used but is local formulation."));
	SymmetricTensor<2, dim> plStrain = get_plStrain(gamma, perCrystalData);
	SymmetricTensor<2, dim> thermalStrain = get_thermalStrain(temperature,
			perCrystalData);

	SymmetricTensor<2, dim> stress = perCrystalData->get_Cel()
			* (strain - plStrain - thermalStrain);

	return stress;
}

/*
 * returns the Schmidstress on system alpha
 * based on the stress
 */
template<int dim, int nslip>
double get_TauAlpha(int alpha, SymmetricTensor<2, dim> stress,
		PerCrystalData<dim, nslip> *perCrystalData) {
	Assert(!perCrystalData->isLocalFormulation(),
			ExcMessage("gradientCP namespace used but is local formulation."));
	return stress * perCrystalData->get_P(alpha);
}

/*
 * returns the Schmidstress on system alpha
 * based on strains and so on
 */
template<int dim, int nslip>
double get_TauAlpha(int alpha, SymmetricTensor<2, dim> strain,
		const Vector<double> &gamma, double temperature,
		PerCrystalData<dim, nslip> *perCrystalData) {
	Assert(!perCrystalData->isLocalFormulation(),
			ExcMessage("gradientCP namespace used but is local formulation."));
	SymmetricTensor<2, dim> stress = get_stress(strain, gamma, temperature);
	return stress * perCrystalData->get_P(alpha);
}

//----------------------------------
//-------------TANGENT--------------
//----------------------------------
template<int dim, int nslip> //checked 27.6.16
SymmetricTensor<4, dim> getDsigDeps(
		PerCrystalData<dim, nslip> *perCrystalData) {
	Assert(!perCrystalData->isLocalFormulation(),
			ExcMessage("gradientCP namespace used but is local formulation."));
	return perCrystalData->get_Cel();
}

template<int dim, int nslip> //checked 27.6.16
SymmetricTensor<2, dim> getDsigDgammaAlpha(int alpha,
		PerCrystalData<dim, nslip> *perCrystalData) {
	return -1.0 * perCrystalData->get_Cel_P(alpha);
}

template<int dim, int nslip> //checked 27.6.16
SymmetricTensor<2, dim> getDxiAlphaDgradgammaBeta(int alpha, int beta,
		double temperature, PerCrystalData<dim, nslip> *perCrystalData) {
	Assert(!perCrystalData->isLocalFormulation(),
			ExcMessage("gradientCP namespace used but is local formulation."));
	SymmetricTensor<2, dim> retval;
	if (alpha != beta) {
		return retval;
	}

	double piTemp = getPiTemp(temperature, perCrystalData);
	return piTemp * perCrystalData->get_DxiAlphaDgradgammaAlpha(alpha);
}

template<int dim, int nslip> //checked 27.6.16
SymmetricTensor<2, dim> getDtauAlphaDeps(int alpha,
		PerCrystalData<dim, nslip> *perCrystalData) {
	Assert(!perCrystalData->isLocalFormulation(),
			ExcMessage("gradientCP namespace used but is local formulation."));
	return perCrystalData->get_Cel_P(alpha);
}

template<int dim, int nslip> //checked 27.6.16
double getDtauAlphaDgammaBeta(int alpha, int beta,
		PerCrystalData<dim, nslip> *perCrystalData) {
	Assert(!perCrystalData->isLocalFormulation(),
			ExcMessage("gradientCP namespace used but is local formulation."));
	return perCrystalData->get_DtauAlphaDgammaBeta(alpha, beta);
}

template<int dim, int nslip> //checked 27.6.16
double getDpiAlphaDgammaBeta_Reddy_sech(int alpha, int beta, double dt,
		double temperature, Vector<double> dGamma,
		const Vector<double> &current_gamma_abs,
		const Vector<double> &current_hardening,
		PerCrystalData<dim, nslip> *perCrystalData) {
	Assert(!perCrystalData->isLocalFormulation(),
			ExcMessage("gradientCP namespace used but is local formulation."));
	double retval = 0;
	for (unsigned int i = 0; i < dGamma.size(); i++) {
		if (dGamma[i] == 0) {
			dGamma[i] = 1e-10;
		}
	}

	double piTemp = getPiTemp(temperature, perCrystalData);
	double updated_A = getA(current_gamma_abs);
	double updated_g_alpha = current_hardening(alpha);
//			getGAlpha(alpha, temperature, dGamma,
//			prior_gamma_abs, perCrystalData);
	double hardening_term = hardeningTerm(alpha, beta, perCrystalData);

	if (alpha == beta) {		//first summand
		retval = (piTemp + updated_g_alpha) * perCrystalData->get_p()
				/ perCrystalData->get_gamma0dot() / dt
				* pow(abs(dGamma[alpha]) / perCrystalData->get_gamma0dot() / dt,
						perCrystalData->get_p() - 1);
	}

	// sinh and cosh terms for hardening law
	double h0A_taus0 = perCrystalData->get_h0() * updated_A
			/ (perCrystalData->get_piS() - piTemp);
	double cosh_term = (exp(h0A_taus0) + exp(-h0A_taus0)) / 2;
	double sinh_term = (exp(h0A_taus0) - exp(-h0A_taus0)) / 2;

	double dgAlpha_dgammaBeta = perCrystalData->get_h0() / pow(cosh_term, 2)
			* hardening_term * (dGamma[beta] < 0 ? -1 : 1);
	for (int Theta = 0; Theta < nslip; Theta++) {
		double dhAlphaA_dgammaTheta = hardening_term * 2
				* perCrystalData->get_h0() * (-sinh_term / pow(cosh_term, 3))
				* perCrystalData->get_h0() * updated_A
				/ (perCrystalData->get_piS() - piTemp)
				* (dGamma[beta] < 0 ? -1 : 1);
		dgAlpha_dgammaBeta += dhAlphaA_dgammaTheta * abs(dGamma[Theta]);
	}
	//second summand linear hardening, only selfhardening
	retval += dgAlpha_dgammaBeta
			* pow(abs(dGamma[alpha]) / perCrystalData->get_gamma0dot() / dt,
					perCrystalData->get_p()) * (dGamma[alpha] < 0 ? -1 : 1);

	return retval;
}

template<int dim, int nslip> //checked 27.6.16
double getDpiAlphaDgammaBeta_Miehe_lin(int alpha, int beta, double dt,
		double temperature, const Vector<double> &dGamma,
		const Vector<double> &current_hardening,
		PerCrystalData<dim, nslip> *perCrystalData) {
	Assert(!perCrystalData->isLocalFormulation(),
			ExcMessage("gradientCP namespace used but is local formulation."));
	perCrystalData->set_dt(dt);
	double retval = 0;

	double piTemp = getPiTemp(temperature, perCrystalData);
	double updated_g_alpha = current_hardening(alpha);

	double tanh_term = tanh(
			dGamma[alpha] / dt / perCrystalData->get_gamma0dot());
	double hardening_term = hardeningTerm(alpha, beta, perCrystalData);

	if (alpha == beta) {		//first summand
		retval = ((piTemp + updated_g_alpha) / perCrystalData->get_gamma0dot()
				* (1 - pow(tanh_term, 2)) + perCrystalData->get_eta()) / dt;
	}

	// dg -> dh-at_dgamma-b 
	retval += tanh_term * perCrystalData->get_h0() * hardening_term
			* (dGamma[beta] < 0 ? -1 : 1);

	return retval;
}
//TODO
template<int dim, int nslip> 
double getDpiAlphaDgammaBeta_Bardella_lin(int alpha, int beta, double dt,
		double temperature, const Vector<double> &dGamma,
		const Vector<double> &current_hardening,
		PerCrystalData<dim, nslip> *perCrystalData) {
	Assert(!perCrystalData->isLocalFormulation(),
			ExcMessage("gradientCP namespace used but is local formulation."));
	perCrystalData->set_dt(dt);
	double retval = 0;
	
	double piTemp = getPiTemp(temperature, perCrystalData);
	double updated_g_alpha = current_hardening(alpha);
	double hardening_term = hardeningTerm(alpha, beta, perCrystalData);

	double AbsGammaAlphaDot = abs(dGamma[alpha]) / dt;
	//double Gamma0Dot = perCrystalData->get_gamma0dot;
	double p = perCrystalData->get_p();
	if(AbsGammaAlphaDot <= perCrystalData->get_gamma0dot()){
		if(alpha == beta){
			retval = (piTemp + updated_g_alpha)/ 2 / perCrystalData->get_gamma0dot()/ dt ;
		}
		retval += AbsGammaAlphaDot/2/perCrystalData->get_gamma0dot() 
		* perCrystalData->get_h0() * hardening_term;
	}else{
		if(alpha == beta){
			retval = (piTemp + updated_g_alpha)
			* pow(AbsGammaAlphaDot/perCrystalData->get_gamma0dot(),p-2)/2/perCrystalData->get_gamma0dot()/dt
			;
		}
		retval += (p - 2 + pow(AbsGammaAlphaDot / perCrystalData->get_gamma0dot(), p - 1))/ 2 / (p - 1)
		* perCrystalData->get_h0() * hardening_term;
	}
	return retval;
}	
	
template<int dim, int nslip> //checked 27.6.16
double getDpiAlphaDgammaBeta(int alpha, int beta, double dt, double temperature,
		const Vector<double> &dGamma, const Vector<double> &current_gamma_abs,
		const Vector<double> &current_hardening,
		PerCrystalData<dim, nslip> *perCrystalData) {
	Assert(!perCrystalData->isLocalFormulation(),
			ExcMessage("gradientCP namespace used but is local formulation."));
	switch (perCrystalData->get_regularisation()) {
	case tanh_lin:
		return getDpiAlphaDgammaBeta_Miehe_lin(alpha, beta, dt, temperature,
				dGamma, current_hardening, perCrystalData);
		break;
	case powerlaw_sech:
		return getDpiAlphaDgammaBeta_Reddy_sech(alpha, beta, dt, temperature,
				dGamma, current_gamma_abs, current_hardening, perCrystalData);
		break;
	case bardella_lin:
		return getDpiAlphaDgammaBeta_Bardella_lin(alpha, beta, dt, temperature,
				dGamma,current_hardening, perCrystalData);
		break;
	default:
		cerr << "Formulation unknown" << endl;
		return -1;
	};
}

template<int dim, int nslip> //checked 28.6.16
double get_PiAlpha_Reddy(int alpha, double dt, double temperature,
		const Vector<double> &dGamma, const Vector<double> &current_hardening,
		PerCrystalData<dim, nslip> *perCrystalData) {
	Assert(!perCrystalData->isLocalFormulation(),
			ExcMessage("gradientCP namespace used but is local formulation."));
	double dgammaAlpha = dGamma[alpha];
	if (dgammaAlpha == 0) {
		dgammaAlpha = 1e-10;
	}

	double piTemp = getPiTemp(temperature, perCrystalData);
	double updated_g_alpha = current_hardening(alpha);

	double retval = (piTemp + updated_g_alpha)
			* pow(abs(dgammaAlpha) / perCrystalData->get_gamma0dot() / dt,
					perCrystalData->get_p()) * (dgammaAlpha > 0 ? 1 : -1);

	return retval;
}

template<int dim, int nslip> //checked 28.6.16
double get_PiAlpha_Miehe(int alpha, double dt, double temperature,
		const Vector<double> &dGamma, const Vector<double> &current_hardening,
		PerCrystalData<dim, nslip> *perCrystalData) {
	Assert(!perCrystalData->isLocalFormulation(),
			ExcMessage("gradientCP namespace used but is local formulation."));

	double piTemp = getPiTemp(temperature, perCrystalData);
	double updated_g_alpha = current_hardening(alpha);

	double tanh_term = tanh(
			dGamma[alpha] / dt / perCrystalData->get_gamma0dot());
	double retval = (piTemp + updated_g_alpha) * tanh_term
			+ perCrystalData->get_eta() * dGamma[alpha] / dt;
 
	return retval;
}
//TODO
template<int dim, int nslip> 
double get_PiAlpha_Bardella_lin(int alpha, double dt, double temperature,
		const Vector<double> &dGamma, const Vector<double> &current_hardening,
		PerCrystalData<dim, nslip> *perCrystalData) {
	Assert(!perCrystalData->isLocalFormulation(),
			ExcMessage("gradientCP namespace used but is local formulation."));

	double piTemp = getPiTemp(temperature, perCrystalData);
	double updated_g_alpha = current_hardening(alpha);
	double AbsGammaAlphaDot = abs(dGamma[alpha]/dt);
	double p = perCrystalData->get_p();
	//deallog<<"p="<<p<<endl;
	//deallog<<"p0="<<perCrystalData->get_gamma0dot()<<endl;
	if (AbsGammaAlphaDot <= perCrystalData->get_gamma0dot()){
		double retval = (piTemp + updated_g_alpha)
		*AbsGammaAlphaDot/perCrystalData->get_gamma0dot()/2* (dGamma[alpha] > 0 ? 1 : -1);
		return retval;
	}
	else{
		double retval = (piTemp + updated_g_alpha)
		*(p-2+pow((AbsGammaAlphaDot/perCrystalData->get_gamma0dot()),p-1))/2/(p-1)* (dGamma[alpha] > 0 ? 1 : -1);
		return retval;
	}

}
//
template<int dim, int nslip> //checked 28.6.16
double get_PiAlpha(int alpha, double dt, double temperature,
		const Vector<double> &dGamma, const Vector<double> &current_hardening,
		PerCrystalData<dim, nslip> *perCrystalData) {
	Assert(!perCrystalData->isLocalFormulation(),
			ExcMessage("gradientCP namespace used but is local formulation."));
	switch (perCrystalData->get_regularisation()) {
	case tanh_lin:
		return get_PiAlpha_Miehe(alpha, dt, temperature, dGamma,
				current_hardening, perCrystalData);
		break;
	case powerlaw_sech:
		return get_PiAlpha_Reddy(alpha, dt, temperature, dGamma,
				current_hardening, perCrystalData);
		break;
	case bardella_lin:
		return get_PiAlpha_Bardella_lin(alpha, dt, temperature, dGamma,
				current_hardening, perCrystalData);
		break;
	case unknown:
		Assert(false, ExcMessage("Unknown regularisation"))
		break;
	};
	return INFINITY;
}

template<int dim, int nslip> //
double get_gammaDot_Miehe(const double dt, const double temperature,
		const double pi, const double current_hardening,
		PerCrystalData<dim, nslip> *perCrystalData) {
	double piTemp = getPiTemp(temperature, perCrystalData);
	double updated_g_alpha = current_hardening;

	double retval = dt * perCrystalData->get_gamma0dot()
			* atanh(pi / (piTemp + updated_g_alpha));
	return retval;
}

template<int dim, int nslip>//TODO
double get_gammaDot(const double dt, const double temperature,
		const double pi, const double current_hardening,
		PerCrystalData<dim, nslip> *perCrystalData) {
	Assert(!perCrystalData->isLocalFormulation(),
			ExcMessage("gradientCP namespace used but is local formulation."));
	switch (perCrystalData->get_regularisation()) {
	case tanh_lin:
		return get_gammaDot_Miehe(dt, temperature, pi,
				current_hardening, perCrystalData);
		break;
	case powerlaw_sech:
		Assert(false, ExcMessage("Not implemented yet"))
		break;
	case bardella_lin:
		Assert(false, ExcMessage("Not implemented yet"))
		break;
	case unknown:
		Assert(false, ExcMessage("Unknown regularisation"))
		break;
	};
	return INFINITY;
}

template<int dim, int nslip> //checked 28.6.16
Tensor<1, dim> get_XiAlpha(vector<Tensor<1, dim> > &gradGamma, int alpha,
		double temperature, PerCrystalData<dim, nslip> *perCrystalData) {
	Assert(!perCrystalData->isLocalFormulation(),
			ExcMessage("gradientCP namespace used but is local formulation."));
	double piTemp = getPiTemp(temperature, perCrystalData);

	return piTemp * perCrystalData->get_le() * perCrystalData->get_le()
			* (perCrystalData->get_c1() * perCrystalData->get_sxs(alpha)
					* gradGamma[alpha]
					+ perCrystalData->get_c2() * perCrystalData->get_lxl(alpha)
							* gradGamma[alpha]);
}

}

//end namespace

#endif /* CRYSTALPLASTICITYCREEP_H_ */
