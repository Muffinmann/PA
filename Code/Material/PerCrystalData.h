#ifndef PERCRYSTALDATA_H_
#define PERCRYSTALDATA_H_

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
#include "CrystalPlasticityCreep.h"
#include <boost/lexical_cast.hpp>
#include "../Various/Tensor8.h"
#include "../Various/AuxilaryFunctions.h"

template<int dim, int nslip>
class Material;
//TODO
enum regularisations {
	unknown = -1, powerlaw_sech = 0, tanh_lin = 1, bardella_lin = 2,
};
//
DeclException1(ExcIndexOutOfBounce, double,
		<<"Index out of bounce, Index was: "<<arg1);

using namespace dealii;
using namespace std;

template<int dim, int nslip>
class PerCrystalData {
public:
	PerCrystalData();
	~PerCrystalData();
	void init(ParameterHandler &parameter, double Euler1, double Euler2,
			double Euler3, int ID, int layer);

	// GETTER
	static double get_C_AB_ab(const PerCrystalData<dim, nslip> *crystalA,
			const PerCrystalData<dim, nslip> *crystalB,
			const Tensor<1, dim> normal_A, const Tensor<1, dim> normal_B,
			const int a, const int b);
	static double get_C_AB_ab(Material<dim, nslip> *crystalA,
			Material<dim, nslip> *crystalB, const Tensor<1, dim> normal_A,
			const Tensor<1, dim> normal_B, const int a, const int b);

	int get_NSlipSystems() const {
		return P.size();
	}
	SymmetricTensor<4, dim> get_Cel() const {
		return stress_strain_cubic;
	}
	vector<SymmetricTensor<2, dim>> get_P() const {
		vector<SymmetricTensor<2, dim>> retval = P;
//		for (unsigned int i=0;i<retval.size();++i){
//			retval[i] *= 0.01;
//		}
		return retval;
	}
	SymmetricTensor<2, dim> get_Cel_P(int a) const {
		return Cel_P[a];
	}
	SymmetricTensor<2, dim> get_DxiAlphaDgradgammaAlpha(int a) const {
		return DxiAlphaDgradgammaAlpha[a];
	}
	double get_DtauAlphaDgammaBeta(int a, int b) const {
		return DtauAlphaDgammaBeta[a][b];
	}
	SymmetricTensor<2, dim> get_sxs(int a) const {
		Assert(a < sxs.size() && a >= 0, ExcIndexOutOfBounce(a));
		return sxs[a];
	}
	SymmetricTensor<2, dim> get_lxl(int a) const {
		Assert(a < lxl.size() && a >= 0, ExcIndexOutOfBounce(a));
		return lxl[a];
	}
	Tensor<1, dim> get_s(int a) const {
		Assert(a < s_rotated.size() && a >= 0, ExcIndexOutOfBounce(a));
		return s_rotated[a];
	}
	Tensor<1, dim> get_m(int a) const {
		Assert(a < m_rotated.size() && a >= 0, ExcIndexOutOfBounce(a));
		return m_rotated[a];
	}
	Tensor<1, dim> get_l(int a) const {
		Assert(a < l_rotated.size() && a >= 0, ExcIndexOutOfBounce(a));
		return l_rotated[a];
	}
	SymmetricTensor<2, dim> get_P(int a) const {
		Assert(a < P.size() && a >= 0, ExcIndexOutOfBounce(a));
		return get_P()[a];
	}
	Vector<double> get_direction(int dir);
	int getDomainID() const {
		return domainID;
	}
	double get_alpha() const {
		return alpha;
	}
	double get_p() const {//TODO
		if (regularisation == tanh_lin) {
			return p;
		} else if (regularisation == powerlaw_sech) {
			return p / regularisationFactor;
		} else if(regularisation == bardella_lin){
			return p / regularisationFactor;
		}
		deallog << "--!--unknown regularization type--!--";
		throw std::bad_exception();
	}
	double get_Nh() const {
		return Nh;
	}
	double get_h0() const {
		return h0;
	}
	double get_pi0() const {
		return pi0;
	}
	double get_piS() const {
		return piS;
	}
	double get_hardeningtype() const {
		return hardeningtype;
	}
	double get_gamma0dot() const {//TODO
		if (regularisation == tanh_lin) {
			return gamma0dot / regularisationFactor;
		} else if (regularisation == powerlaw_sech) {
			return gamma0dot;
		} else if (regularisation == bardella_lin){
			return gamma0dot/regularisationFactor;
		}
		deallog << "--!--unknown regularization type--!--";
		throw std::bad_exception();
	}
	double get_regularisationFactor() {
		return regularisationFactor;
	}
	void set_regularisationFactor(double val) {
		regularisationFactor = val;
	}
	double get_eta() const {
		return eta;
	}
	double get_le() const {
		return le;
	}
	double get_c1() const {
		return c1;
	}
	double get_c2() const {
		return c2;
	}
	bool is_initiated() {
		return isInitiated;
	}
	double get_dt() {
		return dt;
	}
	void set_dt(double value) {
		dt = value;
	}
	regularisations get_regularisation() {
		return regularisation;
	}
	bool isLocalFormulation() {
		return localFormulation;
	}
	vector<double> getEulerAngles() {
		return {eulerangle1, eulerangle2, eulerangle3};
	}

private:
	void init_C_e();
	void init_P();
	void init_sxs();
	void init_lxl();
	void init_Cel_P();
	void init_DxiAlphaDgradgammaAlpha();
	void init_DtauAlphaDgammaBeta();

	Tensor<2, dim> getRotationMatrix();

	static const std::vector<Tensor<1, dim> > m;
	static const std::vector<Tensor<1, dim> > s;
	// elastic materialdata
	double C11;
	double C12;
	double C44;
	SymmetricTensor<4, dim> stress_strain_cubic;

	// plastic materialdata
	double p; //exponent of evolutionequation
	double h0;
	double pi0;
	double piS;
	double hardeningtype;
	double gamma0dot;
	double regularisationFactor;
	double eta;
	double le;
	double c1;
	double c2;
	regularisations regularisation;
	double Nh;

	// Thermal parameter
	double alpha; // themal expansion

	// orientation
	double eulerangle1;
	double eulerangle2;
	double eulerangle3;
	Tensor<2, dim> correctionB4;
	vector<SymmetricTensor<2, dim> > P;
	vector<SymmetricTensor<2, dim> > sxs;
	vector<SymmetricTensor<2, dim> > lxl;
	vector<Tensor<1, dim> > s_rotated;
	vector<Tensor<1, dim> > m_rotated;
	vector<Tensor<1, dim> > l_rotated;
	vector<SymmetricTensor<2, dim> > Cel_P;
	vector<SymmetricTensor<2, dim> > DxiAlphaDgradgammaAlpha;
	vector<vector<double> > DtauAlphaDgammaBeta;

	// Geometry
	int layer;
	double csAngle;
	int domainID;

	double dt;

	bool isInitiated;
	bool localFormulation;
};

//-------------------------------------------
//PerCellData
//-------------------------------------------
//Intitialize Classvariables

template<int dim, int nslip>
const std::vector<Tensor<1, dim>> PerCrystalData<dim, nslip>::s = load_s<dim,
		nslip>();

template<int dim, int nslip>
const std::vector<Tensor<1, dim>> PerCrystalData<dim, nslip>::m = load_m<dim,
		nslip>();

//-------------------------------------------

template<int dim, int nslip> PerCrystalData<dim, nslip>::PerCrystalData() {
// elastic materialdata
	C11 = 0;
	C12 = 0;
	C44 = 0;
	stress_strain_cubic = 0;

	c1 = 0;
	c2 = 0;
	le = 0;
	pi0 = 0;
	piS = 0;
	gamma0dot = 0;
	regularisationFactor = 0;
	h0 = 0;
	eta = 0;
	p = 0;
	hardeningtype = 0;
	Nh = 0;

	alpha = 0;

// orientation
	eulerangle1 = 0;
	eulerangle2 = 0;
	eulerangle3 = 0;
	domainID = 0;
	csAngle = 0;
	layer = 0;

	dt = 0;
	regularisation = unknown;

	isInitiated = false;
	localFormulation = false;
}

template<int dim, int nslip> PerCrystalData<dim, nslip>::~PerCrystalData() {
}

template<int dim, int nslip> void PerCrystalData<dim, nslip>::init(
		ParameterHandler &parameter, double Euler1, double Euler2,
		double Euler3, int ID, int layer) {
	double E = parameter.get_double("EModul");
	double nu = parameter.get_double("nu");
	double mu = parameter.get_double("mu");
// Quelle Anandarajah Computational Methods in Elasticity and Plasicity S. 87
	// C8.2.2-81
	C12 = E * nu / (1 + nu) / (1 - 2 * nu);
	C11 = E * (1 - nu) / (1 - nu - 2 * nu * nu);
	C44 = mu;

	c1 = parameter.get_double("c1");
	c2 = parameter.get_double("c2");
	le = parameter.get_double("le");
	p = parameter.get_double("p");
	h0 = parameter.get_double("h0");
	pi0 = parameter.get_double("pi0");
	piS = parameter.get_double("piS");
	gamma0dot = parameter.get_double("gamma0dot");
	regularisationFactor = 1;
	eta = parameter.get_double("eta");
	hardeningtype = parameter.get_double("hardeningtype");
	Nh = parameter.get_double("Nh");

	eulerangle1 = Euler1;
	eulerangle2 = Euler2;
	eulerangle3 = Euler3;

	if (parameter.get("orientations") == "orientations_B4") {
		correctionB4 = rotmat<dim>(parameter.get_double("euler1")-parameter.get_double("euler3"), 2)
				* rotmat<dim>(parameter.get_double("euler2"), 0)
				* rotmat<dim>(parameter.get_double("euler3"), 2);
	} else {
		correctionB4 = rotmat<dim>(0, 0);
	}

	domainID = ID;
	this->layer = layer;
	if (parameter.get_bool("crosssnake")) {
		csAngle = parameter.get_double("csAngle");
	} else {
		csAngle = 0;
	}
	init_P();
	init_C_e();
	init_sxs();
	init_lxl();
	init_Cel_P();
	init_DxiAlphaDgradgammaAlpha();
	init_DtauAlphaDgammaBeta();

	alpha = parameter.get_double("alpha");

	if (parameter.get("regularisation").compare("tanh_lin") == 0) {
		regularisation = tanh_lin;
	} else if (parameter.get("regularisation").compare("powerlaw_sech") == 0) {
		regularisation = powerlaw_sech;
	} else if (parameter.get("regularisation").compare("bardella_lin" ) == 0) { //TODO
		regularisation = bardella_lin;
	} else {
		Assert(false, ExcMessage("Unknown Regularisation"));
	}

	isInitiated = true;
	localFormulation = parameter.get_bool("localCP");
}

template<int dim, int nslip>
void PerCrystalData<dim, nslip>::init_P() {
	Tensor<2, dim> rotationmatrix = getRotationMatrix();

	// nslip == 0 for local formulation, set to 12
	int ns = nslip;
	if (nslip == 0) {
		ns = 12;
	}

	P.resize(ns);
	for (unsigned int i = 0; i < P.size(); i++) {
		Tensor<2, dim> tmp;
		tmp = outer_product(rotationmatrix * s[i] / s[i].norm(),
				rotationmatrix * m[i] / m[i].norm());
		P[i] = 0.5 * (tmp + transpose(tmp));

	}
}

template<int dim, int nslip>
void PerCrystalData<dim, nslip>::init_sxs() {
	Tensor<2, dim> rotationmatrix = getRotationMatrix();

	// nslip == 0 for local formulation, set to 12
	int ns = nslip;
	if (nslip == 0) {
		ns = 12;
	}
	sxs.resize(ns);
	s_rotated.resize(ns);
	m_rotated.resize(ns);
	for (unsigned int i = 0; i < sxs.size(); i++) {
		Tensor<2, dim> tmp;
		s_rotated[i] = rotationmatrix * s[i] / s[i].norm();
		m_rotated[i] = rotationmatrix * m[i] / m[i].norm();
		sxs[i] = outer_product(s_rotated[i], s_rotated[i]);
	}
}

template<int dim, int nslip>
void PerCrystalData<dim, nslip>::init_lxl() {
	Tensor<2, dim> rotationmatrix = getRotationMatrix();

	lxl.resize(nslip);
	l_rotated.resize(nslip);
	for (unsigned int i = 0; i < lxl.size(); i++) {
		l_rotated[i] = cross_product_3d(rotationmatrix * m[i] / m[i].norm(),
				rotationmatrix * s[i] / s[i].norm());
		lxl[i] = outer_product(l_rotated[i], l_rotated[i]);
	}
}

template<int dim, int nslip>
void PerCrystalData<dim, nslip>::init_Cel_P() {

	// nslip == 0 for local formulation, set to 12
	int ns = nslip;
	if (nslip == 0) {
		ns = 12;
	}
	Cel_P.resize(ns);
	for (unsigned int i = 0; i < Cel_P.size(); i++) {
		double_contract(Cel_P[i], get_Cel(), get_P(i));
	}
}

template<int dim, int nslip>
void PerCrystalData<dim, nslip>::init_DxiAlphaDgradgammaAlpha() {
	DxiAlphaDgradgammaAlpha.resize(nslip);
	for (unsigned int i = 0; i < DxiAlphaDgradgammaAlpha.size(); i++) {
		DxiAlphaDgradgammaAlpha[i] = le * le
				* (c1 * get_sxs(i) + c2 * get_lxl(i));
	}
}

template<int dim, int nslip>
void PerCrystalData<dim, nslip>::init_DtauAlphaDgammaBeta() {
	DtauAlphaDgammaBeta.resize(nslip);
	for (unsigned int i = 0; i < DxiAlphaDgradgammaAlpha.size(); i++) {
		DtauAlphaDgammaBeta[i].resize(nslip);
		for (unsigned int j = 0; j < DtauAlphaDgammaBeta[i].size(); j++) {
			SymmetricTensor<2, dim> tmp;
			double_contract(tmp, get_P(i), get_Cel());
			DtauAlphaDgammaBeta[i][j] = -tmp * get_P(j);
		}
	}
}

template<int dim, int nslip> void PerCrystalData<dim, nslip>::init_C_e() {
//	SymmetricTensor<4, dim> tmp = stress_strain_cubic_symmetric<dim>(
//			kappa + 4.0 / 3 * mu, kappa - 2.0 / 3 * mu, mu);
	SymmetricTensor<4, dim> tmp = stress_strain_cubic_symmetric<dim>(C11, C12,
			C44);
//	cout << eulerangle1 << "/" << eulerangle2 << "/" << eulerangle3 << endl;
	Tensor<2, dim> R = getRotationMatrix();
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

template<int dim, int nslip>
Vector<double> PerCrystalData<dim, nslip>::get_direction(int dir) {
	Assert(dir>=0 && dir < 3, ExcMessage("Direction wrong"));
	Tensor<2, dim> rotationmatrix = getRotationMatrix();
	Vector<double> tmp;
	tmp.reinit(dim);
	tmp(0) = rotationmatrix[TableIndices<2>(0, dir)];
	tmp(1) = rotationmatrix[TableIndices<2>(1, dir)];
	tmp(2) = rotationmatrix[TableIndices<2>(2, dir)];
	return tmp;
}

//template<int dim, int nslip>
//double PerCrystalData<dim, nslip>::get_C_AB_ab(Material<dim, nslip> *crystalA,
//		Material<dim, nslip> *crystalB, const Tensor<1, dim> normal_A,
//		const Tensor<1, dim> normal_B, const int a, const int b) {
//	return get_C_AB_ab(crystalA->get_per_cell_data(),
//			crystalB->get_per_cell_data(), normal_A, normal_B, a, b);
//}

template<int dim, int nslip>
double PerCrystalData<dim, nslip>::get_C_AB_ab(
		const PerCrystalData<dim, nslip> *crystalA,
		const PerCrystalData<dim, nslip> *crystalB,
		const Tensor<1, dim> normal_A, const Tensor<1, dim> normal_B,
		const int a, const int b) {
	Assert(
			(normal_A+normal_B).norm() < 1e-6 || (normal_A-normal_B).norm() < 1e-6,
			StandardExceptions::ExcMessage("Normal_A and Normal_B not parallel."));
	const Tensor<1, dim> s_A_a = crystalA->get_s(a);
	const Tensor<1, dim> m_A_a = crystalA->get_m(a);
	const Tensor<1, dim> s_B_b = crystalB->get_s(b);
	const Tensor<1, dim> m_B_b = crystalB->get_m(b);
	return contract<0, 0>(s_A_a, s_B_b)
			* contract<0, 0>(cross_product_3d(m_A_a, normal_A),
					cross_product_3d(m_B_b, normal_B));
}

/*
 * rotationmatrix to rotate a body not a coordinate system
 */
template<int dim, int nslip>
Tensor<2, dim> PerCrystalData<dim, nslip>::getRotationMatrix() {
	Tensor<2, dim> rotationmatrix;
	rotationmatrix = rotmat<dim>(eulerangle1, 2) * rotmat<dim>(eulerangle2, 0)
			* rotmat<dim>(eulerangle3, 2);

	Tensor<2, dim> rotationCS;
	/*
	 * Eulerangles in degree
	 */
	if (layer == 0) {
		rotationCS = rotmat<dim>(-csAngle, 1);
	} else if (layer == 1) {
		rotationCS = rotmat<dim>(csAngle, 0);
	} else if (layer == 2) {
		rotationCS = rotmat<dim>(csAngle, 1);
	} else if (layer == 3) {
		rotationCS = rotmat<dim>(-csAngle, 0);
	} else {
		deallog << "getRotationMatrix: Layer not in specified range of [1; 4]"
				<< endl;
		throw bad_exception();
	}

	return correctionB4 * rotationCS * rotationmatrix;
}

#endif /* PERCRYSTALDATA_H_ */
