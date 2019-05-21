/*
 * Materialclass_elastoplastic.h
 *
 *  Created on: Jan 19, 2016
 *      Author: iwtm84
 */

#ifndef MATERIALCLASS_ELASTOPLASTIC_H_
#define MATERIALCLASS_ELASTOPLASTIC_H_

#include <deal.II/base/symmetric_tensor.h>

namespace Material {
using namespace dealii;

template<int dim>
class Material {
public:
	Material();
	~Material();
	static SymmetricTensor<4, dim> getElasticStressStrainTensor() {
		return elastic_stress_strain_tensor;
	}
	SymmetricTensor<2, dim> old_stress() {
		return stress;
	}

	void output(void);
	SymmetricTensor<4, dim> calculateEPTangentElastoPlastischMises(
			SymmetricTensor<2, dim> strain);
	void update_pointhistorydata();

private:
	static const SymmetricTensor<2, dim> I;
	static const SymmetricTensor<4, dim> IxI;
	static const SymmetricTensor<4, dim> elastic_stress_strain_tensor;
	static const SymmetricTensor<4, dim> deviatoric_identity;

	SymmetricTensor<2, dim> stress;
	SymmetricTensor<2, dim> old_plastic_strain;
	double old_g;
	SymmetricTensor<2, dim> updated_plastic_strain;
	double updated_g;
};

using namespace dealii;

template<int dim>
SymmetricTensor<4, dim> get_stress_strain_tensor(const double lambda,
		const double mu) {
	SymmetricTensor<4, dim> tmp;
	for (unsigned int i = 0; i < dim; ++i)
		for (unsigned int j = 0; j < dim; ++j)
			for (unsigned int k = 0; k < dim; ++k)
				for (unsigned int l = 0; l < dim; ++l)
					tmp[i][j][k][l] = (((i == k) && (j == l) ? mu : 0.0)
							+ ((i == l) && (j == k) ? mu : 0.0)
							+ ((i == j) && (k == l) ? lambda : 0.0));
	return tmp;
}

template<int dim>
const SymmetricTensor<2, dim> Material<dim>::I = unit_symmetric_tensor<dim>();

template<int dim>
const SymmetricTensor<4, dim> Material<dim>::IxI = outer_product(I, I);

template<int dim>
const SymmetricTensor<4, dim> Material<dim>::elastic_stress_strain_tensor =
		get_stress_strain_tensor<dim>(
		/*lambda = */1125,
		/*mu     = */562.5); //E:1500, nu:1/3, kappa:1500

template<int dim>
const SymmetricTensor<4, dim> Material<dim>::deviatoric_identity =
		deviator_tensor<dim>();

template<int dim> Material<dim>::Material() {
	stress = 0;
	old_plastic_strain = 0;
	old_g = 0;
	updated_g = 0;
	updated_plastic_strain = 0;
//	std::cout << "Material erstellt" << std::endl;
}

template<int dim> Material<dim>::~Material() {
//	std::cout << "Material zerstört" << std::endl;
}
template<int dim>
SymmetricTensor<4, dim> Material<dim>::calculateEPTangentElastoPlastischMises(
		SymmetricTensor<2, dim> strain) {

	// Parameter des Werkstoffgesetzes
	double kappa = 1500.0;
	double mu = 562.5;

	double K = 0.0;
	double sig_y = 10.0;

	SymmetricTensor<4, dim> Cep_vol = kappa * IxI; // volumetric part of EPTangente
	// compute trial values
	SymmetricTensor<2, dim> sigma_vol_np1 = kappa * (strain * I) * I;
	SymmetricTensor<2, dim> eps_p_n = old_plastic_strain;
	SymmetricTensor<2, dim> sigma_dev_trial = 2 * mu
			* (strain * deviatoric_identity - eps_p_n);
	double alpha_n = old_g;
	double R_trial = -K * alpha_n;
	double yield_fun = sigma_dev_trial.norm()
			- std::sqrt(2.0 / 3.0) * (sig_y - R_trial);
	if (yield_fun <= 0) {
		stress = sigma_vol_np1 + sigma_dev_trial;
		return Cep_vol + 2.0 * mu * deviatoric_identity;
	}

	double delta_lambda = yield_fun / (2.0 * mu + 2.0 / 3.0 * K);
	SymmetricTensor<2, dim> normal = sigma_dev_trial / sigma_dev_trial.norm();
	SymmetricTensor<4, dim> Cep_dev = 2.0 * mu * deviatoric_identity
			- 4 * mu * mu / (2 * mu + 2.0 / 3.0 * K)
					* outer_product(normal, normal)
			- 4 * mu * mu * delta_lambda / sigma_dev_trial.norm()
					* (deviatoric_identity - outer_product(normal, normal));

//	std::cout << "davor a: " << point_history.old_alpha << "\nsig:      "
//			<< point_history.old_stress << "\nstrain_p: "
//			<< point_history.old_plastic_strain << std::endl;

	stress = sigma_vol_np1 + sigma_dev_trial
			- delta_lambda * 2 * mu * normal;
	updated_g = old_g + delta_lambda * std::sqrt(2.0 / 3.0);
	updated_plastic_strain = old_plastic_strain + delta_lambda * normal;

//	std::cout << "danach a: " << point_history.old_alpha << "\nsig:      "
//			<< point_history.old_stress << "\nstrain_p: "
//			<< point_history.old_plastic_strain << std::endl;
//	std::cout << "plastisch" << std::endl;
	return Cep_vol + Cep_dev;
}

template<int dim> void Material<dim>::update_pointhistorydata() {
	old_g = updated_g;
	old_plastic_strain = updated_plastic_strain;
}

} //End Namespace
#endif /* MATERIALCLASS_ELASTOPLASTIC_H_ */
