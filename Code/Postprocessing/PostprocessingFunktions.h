/*
 * PostprocessingFunktions.h
 *
 *  Created on: Apr 6, 2016
 *      Author: iwtm84
 */

#ifndef CODE_POSTPROCESSING_POSTPROCESSINGFUNKTIONS_H_
#define CODE_POSTPROCESSING_POSTPROCESSINGFUNKTIONS_H_

using namespace std;
using namespace dealii;
#include <algorithm>
#include "../Material/CrystalPlasticityCreep.h"
#include "../Material/PerCrystalData.h"
#include "Scratches.h"

template<int dim, int nslip>
double calculateTiltangle(gpData<dim, nslip>& gpd) {
	Vector<double> tmp1;
	tmp1 = gpd.perCrystalData->get_direction(0);
	double tilt1 = atan(sqrt(tmp1(0) * tmp1(0) + tmp1(1) * tmp1(1)) / tmp1(2))
			/ 3.14159265 * 180;
	tilt1 = abs(tilt1);

	Vector<double> tmp2;
	tmp2 = gpd.perCrystalData->get_direction(1);
	double tilt2 = atan(sqrt(tmp2(0) * tmp2(0) + tmp2(1) * tmp2(1)) / tmp2(2))
			/ 3.14159265 * 180;
	tilt2 = abs(tilt2);

	Vector<double> tmp3;
	tmp3 = gpd.perCrystalData->get_direction(2);
	double tilt3 = atan(sqrt(tmp3(0) * tmp3(0) + tmp3(1) * tmp3(1)) / tmp3(2))
			/ 3.14159265 * 180;
	tilt3 = abs(tilt3);
//	deallog << "tilts: " << tilt1 << " " << tilt2 << " " << tilt3 << endl;

	double tiltangle;
	if (tilt1 < tilt2) {
		if (tilt1 < tilt3) {
			// tilt1 is the smallest
			tiltangle = tilt1;
//			deallog << "tilt1 is the smallest: " << tilt1 << endl;
		} else {
			// tilt3 is the smallest
			tiltangle = tilt3;
//			deallog << "tilt3 is the smallest: " << tilt3 << endl;
		}
	} else {
		if (tilt2 < tilt3) {
			// tilt2 is the smallest
			tiltangle = tilt2;
//			deallog << "tilt2 is the smallest: " << tilt2 << endl;
		} else {
			// tilt3 is the smallest
			tiltangle = tilt3;
//			deallog << "tilt3 is the smallest: " << tilt3 << endl;
		}
	}

//	if (tiltangle > 20) {
//		deallog << "tiltangle: " << tiltangle << endl;
//	}

	return tiltangle;
}

template<int dim, int nslip>
double calculateRotationangle(gpData<dim, nslip>& gpd) {
	// first look for the direction which is pointing in positive or negative z direction
	Vector<double> tmp1;
	tmp1 = gpd.perCrystalData->get_direction(0);
	double tilt1 = atan(sqrt(tmp1(0) * tmp1(0) + tmp1(1) * tmp1(1)) / tmp1(2))
			/ 3.14159265 * 180;
	tilt1 = abs(tilt1);

	Vector<double> tmp2;
	tmp2 = gpd.perCrystalData->get_direction(1);
	double tilt2 = atan(sqrt(tmp2(0) * tmp2(0) + tmp2(1) * tmp2(1)) / tmp2(2))
			/ 3.14159265 * 180;
	tilt2 = abs(tilt2);

	Vector<double> tmp3;
	tmp3 = gpd.perCrystalData->get_direction(2);
	double tilt3 = atan(sqrt(tmp3(0) * tmp3(0) + tmp3(1) * tmp3(1)) / tmp3(2))
			/ 3.14159265 * 180;
	tilt3 = abs(tilt3);

	int zdir;
	if (tilt1 < tilt2) {
		if (tilt1 < tilt3) {
			// tilt1 is the smallest
			zdir = 0;
			//			deallog << "tilt1 is the smallest: " << tilt1 << endl;
		} else {
			// tilt3 is the smallest
			zdir = 2;
			//			deallog << "tilt3 is the smallest: " << tilt3 << endl;
		}
	} else {
		if (tilt2 < tilt3) {
			// tilt2 is the smallest
			zdir = 1;
			//			deallog << "tilt2 is the smallest: " << tilt2 << endl;
		} else {
			// tilt3 is the smallest
			zdir = 2;
			//			deallog << "tilt3 is the smallest: " << tilt3 << endl;
		}
	}

	// calculate the respective roationangle
	Vector<double> tmp;
	tmp = gpd.perCrystalData->get_direction((zdir + 1) % dim);
	double angle = atan(tmp(1) / tmp(0)) / 3.14159265 * 180;
	while (abs(angle) > 45) {
		angle -= 45 * sgn(angle);
	}
	return angle;

////	tmp = history.get_direction(0);
//	tmp = gpd.perCrystalData->get_direction(0);
//	double angle = atan(tmp(1) / tmp(0)) / 3.14159265 * 180;
//	while (abs(angle) > 45) {
//		angle -= 45 * sgn(angle);
//	}
//	return angle;
}

template<int dim, int nslip>
double calculateSig11(gpData<dim, nslip>& gpd) {
	SymmetricTensor<2, dim> stress;

	if (gpd.perCrystalData->isLocalFormulation()) {
		stress = localCP::getStress(gpd.strain, gpd.dt, gpd.temperature,
				gpd.perCrystalData);
	} else {
		stress = GradientCP::get_stress(gpd.strain, gpd.gamma, gpd.temperature,
				gpd.perCrystalData);
	}
	return stress[0][0];
}

template<int dim, int nslip>
double calculateSig22(gpData<dim, nslip>& gpd) {
	SymmetricTensor<2, dim> stress;

	if (gpd.perCrystalData->isLocalFormulation()) {
		stress = localCP::getStress(gpd.strain, gpd.dt, gpd.temperature,
				gpd.perCrystalData);
	} else {
		stress = GradientCP::get_stress(gpd.strain, gpd.gamma, gpd.temperature,
				gpd.perCrystalData);
	}
	return stress[1][1];
}

template<int dim, int nslip>
double calculateSig33(gpData<dim, nslip>& gpd) {
	SymmetricTensor<2, dim> stress;

	if (gpd.perCrystalData->isLocalFormulation()) {
		stress = localCP::getStress(gpd.strain, gpd.dt, gpd.temperature,
				gpd.perCrystalData);
	} else {
		stress = GradientCP::get_stress(gpd.strain, gpd.gamma, gpd.temperature,
				gpd.perCrystalData);
	}
	return stress[2][2];
}

template<int dim, int nslip>
double calculateSig12(gpData<dim, nslip>& gpd) {
	SymmetricTensor<2, dim> stress;

	if (gpd.perCrystalData->isLocalFormulation()) {
		stress = localCP::getStress(gpd.strain, gpd.dt, gpd.temperature,
				gpd.perCrystalData);
	} else {
		stress = GradientCP::get_stress(gpd.strain, gpd.gamma, gpd.temperature,
				gpd.perCrystalData);
	}
	return stress[0][1];
}

template<int dim, int nslip>
double calculateSig13(gpData<dim, nslip>& gpd) {
	SymmetricTensor<2, dim> stress;

	if (gpd.perCrystalData->isLocalFormulation()) {
		stress = localCP::getStress(gpd.strain, gpd.dt, gpd.temperature,
				gpd.perCrystalData);
	} else {
		stress = GradientCP::get_stress(gpd.strain, gpd.gamma, gpd.temperature,
				gpd.perCrystalData);
	}
	return stress[0][2];
}

template<int dim, int nslip>
double calculateSig23(gpData<dim, nslip>& gpd) {
	SymmetricTensor<2, dim> stress;

	if (gpd.perCrystalData->isLocalFormulation()) {
		stress = localCP::getStress(gpd.strain, gpd.dt, gpd.temperature,
				gpd.perCrystalData);
	} else {
		stress = GradientCP::get_stress(gpd.strain, gpd.gamma, gpd.temperature,
				gpd.perCrystalData);
	}
	return stress[1][2];
}

template<int dim, int nslip>
double calculateEpsP12(gpData<dim, nslip>& /*gpd*/) {
	return 0;
//	return history.get_plastic_strain()[0][1];
}

template<int dim, int nslip>
double calculateEquivalentPlasticStrain(gpData<dim, nslip>& gpd) {
	SymmetricTensor<2, dim> plasticStrain;

	if (gpd.perCrystalData->isLocalFormulation()) {
		return -1;
	} else {
		plasticStrain = GradientCP::get_plStrain(gpd.gamma, gpd.perCrystalData);

	}

	return sqrt(
			0.5
					* (pow(plasticStrain[0][0] - plasticStrain[1][1], 2)
							+ pow(plasticStrain[1][1] - plasticStrain[2][2], 2)
							+ pow(plasticStrain[2][2] - plasticStrain[0][0], 2)
							+ 6
									* (pow(plasticStrain[0][1], 2)
											+ pow(plasticStrain[1][2], 2)
											+ pow(plasticStrain[0][2], 2))));
}

template<int dim, int nslip>
double calculateEquivalentPlasticStrainRate(gpData<dim, nslip>& gpd) {

	SymmetricTensor<2, dim> deltaPlasticStrain;
	if (gpd.perCrystalData->isLocalFormulation()) {
		return -1;
	} else {
		deltaPlasticStrain = GradientCP::get_plStrain(gpd.dGamma_actualStep,
				gpd.perCrystalData);
	}

	return sqrt(
			0.5
					* (pow(deltaPlasticStrain[0][0] - deltaPlasticStrain[1][1],
							2)
							+ pow(
									deltaPlasticStrain[1][1]
											- deltaPlasticStrain[2][2], 2)
							+ pow(
									deltaPlasticStrain[2][2]
											- deltaPlasticStrain[0][0], 2)
							+ 6
									* (pow(deltaPlasticStrain[0][1], 2)
											+ pow(deltaPlasticStrain[1][2], 2)
											+ pow(deltaPlasticStrain[0][2], 2))));
}

template<int dim, int nslip>
double calculateEquivalentThermalStrain(gpData<dim, nslip>& gpd) {
	SymmetricTensor<2, dim> thermalStrain;
	if (gpd.perCrystalData->isLocalFormulation()) {
		return -1;
	} else {
		thermalStrain = GradientCP::get_thermalStrain(gpd.temperature,
				gpd.perCrystalData);
	}
	return sqrt(
			0.5
					* (pow(thermalStrain[0][0] - thermalStrain[1][1], 2)
							+ pow(thermalStrain[1][1] - thermalStrain[2][2], 2)
							+ pow(thermalStrain[2][2] - thermalStrain[0][0], 2)
							+ 6
									* (pow(thermalStrain[0][1], 2)
											+ pow(thermalStrain[1][2], 2)
											+ pow(thermalStrain[0][2], 2))));
}

template<int dim, int nslip>
double calculateTemperatur(gpData<dim, nslip>& gpd) {
	return gpd.temperature;
}

template<int dim, int nslip>
double calculateA(gpData<dim, nslip>& gpd) {
	if (gpd.perCrystalData->isLocalFormulation()) {
		return -1;
	} else {
		return GradientCP::getA(gpd.current_gamma_abs);
	}
}

template<int dim, int nslip>
double calculateEdgeDislocationDensity(gpData<dim, nslip>& gpd) {
	double edgeDislocationDensity = 0;
	for (unsigned int i = 0; i < gpd.gradGamma.size(); ++i) {
		edgeDislocationDensity += pow(
				gpd.perCrystalData->get_s(i) * gpd.gradGamma[i], 2);
	}
	return pow(edgeDislocationDensity, 0.5);
}

template<int dim, int nslip>
double calculateScrewDislocationDensity(gpData<dim, nslip>& gpd) {
	double screwDislocationDensity = 0;
	for (unsigned int i = 0; i < gpd.gradGamma.size(); ++i) {
		screwDislocationDensity += pow(
				gpd.perCrystalData->get_l(i) * gpd.gradGamma[i], 2);
	}
	return pow(screwDislocationDensity, 0.5);
}

template<int dim, int nslip>
double calculateDislocationDensity(gpData<dim, nslip>& gpd) {
	double dislocationDensity = 0;
	for (unsigned int i = 0; i < gpd.gradGamma.size(); ++i) {
		dislocationDensity += pow(
				gpd.perCrystalData->get_m(i) * gpd.gradGamma[i], 2);
		dislocationDensity += pow(
				gpd.perCrystalData->get_s(i) * gpd.gradGamma[i], 2);
	}
	return pow(dislocationDensity, 0.5);
}

template<int dim, int nslip>
double calculateDomainID(gpData<dim, nslip>& gpd) {
	return gpd.perCrystalData->getDomainID();
}

template<int dim, int nslip>
double calculateMises(gpData<dim, nslip>& gpd) {
	SymmetricTensor<2, dim> stress;

	if (gpd.perCrystalData->isLocalFormulation()) {
		stress = localCP::getStress(gpd.strain, gpd.dt, gpd.temperature,
				gpd.perCrystalData);
	} else {
		stress = GradientCP::get_stress(gpd.strain, gpd.gamma, gpd.temperature,
				gpd.perCrystalData);
	}

	return sqrt(
			0.5
					* (pow(stress[0][0] - stress[1][1], 2)
							+ pow(stress[1][1] - stress[2][2], 2)
							+ pow(stress[2][2] - stress[0][0], 2)
							+ 6
									* (pow(stress[0][1], 2)
											+ pow(stress[1][2], 2)
											+ pow(stress[0][2], 2))));
}

#endif /* CODE_POSTPROCESSING_POSTPROCESSINGFUNKTIONS_H_ */
