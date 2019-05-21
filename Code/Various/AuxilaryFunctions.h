/*
 * AuxilaryFunctions.h
 *
 *  Created on: Mar 17, 2017
 *      Author: iwtm84
 */

#ifndef CODE_VARIOUS_AUXILARYFUNCTIONS_H_
#define CODE_VARIOUS_AUXILARYFUNCTIONS_H_

#include <algorithm>
#include "../Material/CrystalPlasticityCreep.h"
#include <deal.II/base/symmetric_tensor.h>
#include <iostream>
#include <fstream>
#include <string>

using namespace dealii;
using namespace std;

template<int dim>
Tensor<2, dim> rotmat(double angle, int normal);

namespace Meshgeneration {

vector<double> calculateEulerAngles(const int domainID,
		const ParameterHandler &parameter) {

	//Map for Eulerangles
	std::vector<std::vector<double>> domains(4096, std::vector<double>(3));
	if (parameter.get("orientations") == "orientations_B4") {
#include"../../Input/orientations_B4"

		int offset_domain = parameter.get_integer("offset_domain") + domainID;
		if (parameter.get_bool("randomly_oriented")) {
			double tmp = rand() / double(RAND_MAX);
			offset_domain += domains.size() * tmp;
		}
		offset_domain %= domains.size();
//		deallog << "offset_domain: " << offset_domain << endl;
		double euler1 = domains[offset_domain][0];
		double euler2 = domains[offset_domain][1];
		double euler3 = domains[offset_domain][2];

		vector<double> retval = { euler1, euler2, euler3 };
		return retval;
	} else if (parameter.get("orientations") == "orientations") {
#include"../../Input/orientations"
	} else {
		deallog << "No orientationfile loaded" << endl;
		exit(1);
	}

	// Rot = RotZ(Euler2) * RotY(Euler1) * RotZ(Euler0)
	double factorEulerangle0 = 1;
	double factorEulerangle1 = 1;
	double factorEulerangle2 = 1;

	// Isotrop means random Orientation
	if (!parameter.get_bool("isotrop")) {
		// Anisotrop is the setting for columnar grains with building direction in Z-Axis
		// Tilted {0;tiltangle} around the Y-Axis, the rotated randomly arround new Z
			factorEulerangle0 = 1;
			factorEulerangle1 = parameter.get_double("tiltangle") / 360;
			factorEulerangle2 = 1;
	}

	if (parameter.get("loadtype").compare("testcase_andrew_microhard") == 0
			|| parameter.get("loadtype").compare("testcase_andrew_microfree")
					== 0
			|| parameter.get("loadtype").compare("testcase_grain_boundary") == 0
			|| parameter.get("loadtype").compare("testcase_gottschalk") == 0) {
		factorEulerangle0 = 1;
		factorEulerangle1 = 0;
		factorEulerangle2 = 0;
	}

	double euler1, euler2, euler3;
	if (parameter.get_bool("randomly_oriented")) {
			double alpha = 360 * (rand() / double(RAND_MAX));
			double beta = acos(2 * (rand() / double(RAND_MAX)) - 1)
					/ 3.14159265358979323 * 180;
			double gamma = 360 * (rand() / double(RAND_MAX));
			euler1 = alpha * factorEulerangle0;
			euler2 = beta * factorEulerangle1;
			euler3 = (gamma - 180) * factorEulerangle2;
			//		deallog << "Gamma: " << gamma << " Euler1: " << euler1 << " Euler2: "
			//				<< euler2 << " Euler3: " << euler3 << endl;
	} else { // not random oriented
		if (parameter.get_double("n_crystals") > 1) {

			int offset_domain = parameter.get_integer("offset_domain")
					+ domainID;
			offset_domain %= domains.size();
			euler1 = domains[offset_domain][0] * factorEulerangle0;
			euler2 = domains[offset_domain][1] * factorEulerangle1;
			euler3 = (domains[offset_domain][2] - 180) * factorEulerangle2;

		} else { // if 1 crystal, use orientation given in parameterfile

			euler1 = parameter.get_double("euler1");
			euler2 = parameter.get_double("euler2");
			euler3 = parameter.get_double("euler3");
		}
	}

	// Spezial cases for the testcases andrew, boundary and gottschalk
	if (parameter.get("loadtype").compare("testcase_andrew_microhard") == 0
			|| parameter.get("loadtype").compare("testcase_andrew_microfree")
					== 0
			|| parameter.get("loadtype").compare("testcase_grain_boundary") == 0
			|| parameter.get("loadtype").compare("testcase_gottschalk") == 0) {
		if (parameter.get_double("n_crystals") == 2 && domainID == 0) {
			euler1 = -parameter.get_double("euler2");
			euler2 = 0;
			euler3 = 0;
		}
		if (parameter.get_double("n_crystals") == 2 && domainID == 1) {
			euler1 = -parameter.get_double("euler1");
			euler2 = 0;
			euler3 = 0;
		}
	}

	vector<double> retval = { euler1, euler2, euler3 };
	return retval;
}

/*
 * returns the roationmatrix to rotate to unitcells in the direction of grain growths
 */
template<int dim>
Tensor<2, dim> getCSRotation(ParameterHandler &parameter, int layer) {
	Tensor<2, dim> retval;

	if (parameter.get_bool("crosssnake")) {
		/*
		 * Eulerangles in degree
		 */
		double eulerangle = parameter.get_double("csAngle");
		if (layer == 1) {
			retval = rotmat<dim>(eulerangle, 0);
		} else if (layer == 2) {
			retval = rotmat<dim>(eulerangle, 1);
		} else if (layer == 3) {
			retval = rotmat<dim>(-eulerangle, 0);
		} else if (layer == 4) {
			retval = rotmat<dim>(-eulerangle, 1);
		} else {
			deallog << "getCSRotation: Layer not in specified range of [1;4]"
					<< endl;
			throw bad_exception();
		}
	}
	return retval;
}
} // End Namespace Mesh Generation

/*
 * read string that contains a symmetric tensor in one line
 * used for reading of the written stress_strain files
 * and strains for prescribed load
 */
template<int dim>
SymmetricTensor<2, dim> readTensorLine(string line) {

	int index = line.find_first_of("0123456789-");

	SymmetricTensor<2, dim> retval;
	for (int i = 0; i < dim; ++i) {
		for (int j = 0; j < dim; ++j) {
			int length = line.find_first_of(" ", index) - index;
			string subString = line.substr(index, length);
			//			deallog << "line: '" << line << "' Index: " << index << "/"
			//					<< lastIndex << " subString: '" << subString << "'" << endl;
			double val = boost::lexical_cast<double>(subString);
			index = index + length + 1;

			retval[i][j] = val;
		}
	}
	return retval;
}

/*
 * read string that contains a double in one line
 * used for reading of the written dissipation
 */
double readDoubleLine(string line) {

	int index = line.find_first_of("0123456789-");

	double retval;
	int length = line.find_first_of(" ", index) - index;
	string subString = line.substr(index, length);
	//			deallog << "line: '" << line << "' Index: " << index << "/"
	//					<< lastIndex << " subString: '" << subString << "'" << endl;
	double val = boost::lexical_cast<double>(subString);
	index = index + length + 1;

	retval = val;

	return retval;
}

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

	//	Tensor<4, dim> tmp2;
	//	for (int i = 0; i < dim; i++) {
	//		for (int j = 0; j < dim; j++) {
	//			for (int k = 0; k < dim; k++) {
	//				for (int l = 0; l < dim; l++) {
	//					if (i == j && k == l) {
	//						tmp2[i][j][k][l] += C12;
	//					}
	//
	//					if (i == k && j == l) {
	//						tmp2[i][j][k][l] += C44;
	//					}
	//					if (i == l && j == k) {
	//						tmp2[i][j][k][l] += C44;
	//					}
	//
	//					for (int a = 0; a < dim; a++) {
	//						if (a == i && a == j && a == k && a == l) {
	//							tmp2[i][j][k][l] += C11-C12-2*C44;
	//						}
	//					}
	//				}
	//			}
	//		}
	//	}
	//	print(cout, tmp2);
	//
	//	for (int i = 0; i < dim; i++) {
	//		for (int j = 0; j < dim; j++) {
	//			for (int k = 0; k < dim; k++) {
	//				for (int l = 0; l < dim; l++) {
	//					if (tmp[i][j][k][l] != tmp2[i][j][k][l]) {
	//						cout << i << "/" << j << "/" << k << "/" << l << "/"
	//								<< endl;
	//					}
	//				}
	//			}
	//		}
	//	}
	//
	//	print(cout, tmp);
	//	exit(0);
	return tmp;
}

/*
 * returns a rotationmatrix for a rotation of angle in degree around normal
 * 0: x-Axis
 * 1: y-Axis
 * 2: z-Axis
 */
template<int dim>
Tensor<2, dim> rotmat(double angle, int normal) {
	Assert(dim == 3, ExcMessage("rotmat only for dim = 3 implemented"));
	Assert(normal == 0 || normal == 1 || normal == 2,
			ExcMessage("rotmat only for normal 1, 2 and 3 implemented"));
	angle = angle / 180 * 3.1415926535;
	double cosval = cos(angle);
	double sinval = sin(angle);

	//	deallog << "Cos("<<angle/3.1415926535*180<<"): " << cosval << "   Sin("<<angle/3.1415926535*180<<"): " << sinval << endl;

	// checked 17.08.2018
	Tensor<2, dim> rotmat;
	if (normal == 0) {
		rotmat[TableIndices<2>(0, 0)] = 1;
		rotmat[TableIndices<2>(1, 1)] = cosval;
		rotmat[TableIndices<2>(2, 2)] = cosval;
		rotmat[TableIndices<2>(1, 2)] = -sinval;
		rotmat[TableIndices<2>(2, 1)] = sinval;
	} else if (normal == 1) {
		rotmat[TableIndices<2>(0, 0)] = cosval;
		rotmat[TableIndices<2>(1, 1)] = 1;
		rotmat[TableIndices<2>(2, 2)] = cosval;
		rotmat[TableIndices<2>(0, 2)] = sinval;
		rotmat[TableIndices<2>(2, 0)] = -sinval;
	} else if (normal == 2) {
		rotmat[TableIndices<2>(0, 0)] = cosval;
		rotmat[TableIndices<2>(1, 1)] = cosval;
		rotmat[TableIndices<2>(2, 2)] = 1;
		rotmat[TableIndices<2>(0, 1)] = -sinval;
		rotmat[TableIndices<2>(1, 0)] = sinval;
	}

//	Tensor<1, dim> x;
//	x[0] = 1;
//	Tensor<1, dim> y;
//	y[1] = 1;
//	Tensor<1, dim> z;
//	z[2] = 1;
//	deallog << "Angle: " << angle << " Rotmat (normal:"<<normal<< "): " << rotmat << endl <<
//			" Rotation of [1 0 0]: " << rotmat * x << endl <<
//			" Rotation of [0 1 0]: " << rotmat * y << endl <<
//			" Rotation of [0 0 1]: " << rotmat * z << endl;

	return rotmat;
}

template<int dim, int nslip>
std::vector<Tensor<1, dim>> load_s() {
	// nslip == 0 for local formulation, set to 12
	int ns = nslip;
	if (nslip == 0) {
		ns = 12;
	}
	//load s
	string filename;
	std::vector<Tensor<1, dim>> s(ns);
	switch (ns) {
	case 0:
		return s;
	case 1:
		filename = "Input/Slipsystems/slipdirections_one_slip";
		break;
	case 2:
		filename = "Input/Slipsystems/slipdirections_testcase_andrew";
		break;
	case 12:
		filename = "Input/Slipsystems/slipdirections_fcc";
		break;
	default:
		Assert(false, ExcMessage("Unknown number of Slipsystems"))
		;
	}
	ifstream in(filename);
	if (in) {
		int i = 0, j = 0;
		string line;
		while (getline(in, line)) {
			stringstream sep(line);
			string field;
			j = 0;
			Assert(i < s.size(),
					ExcMessage( "Too many slipdirections in slipdirections_fcc."));
			while (getline(sep, field, ',')) {
				Assert(j < dim,
						ExcMessage( "Too many vectorcomponents in slipdirections_fcc."));
				s[i][j++] = boost::lexical_cast<double>(field);
			}
			i++;
		}
	} else {
		Assert(false, ExcMessage("File slipdirections_fcc not found."));
	}
	return s;
}

template<int dim, int nslip>
std::vector<Tensor<1, dim>> load_m() {
	// nslip == 0 for local formulation, set to 12
	int ns = nslip;
	if (nslip == 0) {
		ns = 12;
	}
	//load m
	string filename;
	std::vector<Tensor<1, dim>> m(ns);
	switch (ns) {
	case 0:
		return m;
	case 1:
		filename = "Input/Slipsystems/slipplanenormals_one_slip";
		break;
	case 2:
		filename = "Input/Slipsystems/slipplanenormals_testcase_andrew";
		break;
	case 12:
		filename = "Input/Slipsystems/slipplanenormals_fcc";
		break;
	default:
		Assert(false, ExcMessage("Unknown number of Slipsystems"))
		;
	}
	ifstream in(filename);
	if (in) {
		int i = 0, j = 0;
		string line;
		while (getline(in, line)) {
			stringstream sep(line);
			string field;
			j = 0;
			Assert(i < m.size(),
					ExcMessage("Too many normals in slipplanenormals_fcc."));
			while (getline(sep, field, ',')) {
				Assert(j < dim,
						ExcMessage( "Too many normals in slipplanenormals_fcc."));
				m[i][j++] = boost::lexical_cast<double>(field);
			}
			i++;
		}
	} else {
		Assert(false, ExcMessage("File slipplanenormals_fcc not found."));
	}
	return m;
}

/*
 * Function to read in textfile with doublevalues, one value per line, one loadstep per line
 * is needed for prescribed loadpath
 */
std::vector<double> loadDoubleArchive(string filename) {
	std::vector<double> retval;
	ifstream in(filename);
	if (in) {
		string line;
		while (getline(in, line)) {
			stringstream sep(line);
			string field;
			while (getline(sep, field, ' ')) {
				retval.push_back(boost::lexical_cast<double>(field));
			}
		}
	}
	return retval;
}

/*
 * Function to read in textfile with doublevalues, one value per line, one loadstep per line
 */
template<int dim>
std::vector<SymmetricTensor<2, dim>> loadSymmetricTensorArchive(
		string filename) {
	std::vector<SymmetricTensor<2, dim>> retval;
	ifstream in(filename);
	if (in) {
		string line;
		while (getline(in, line)) {
			SymmetricTensor<2, dim> tmp = readTensorLine<dim>(line);
			retval.push_back(tmp);
		}
	}
	return retval;
}

template<typename T> int sgn(T x) {
	if (x > 0)
		return 1;
	if (x < 0)
		return -1;
	return 0;
}

template<int dim>
bool writeYieldSurfaceFile(SymmetricTensor<2, dim> macroscopicStrain,
		SymmetricTensor<2, dim> macroscopicStress) {
	std::string filename;

	filename = "Results/YieldSurface.txt";
	std::ofstream File(filename, std::ios::out | std::ios::app);
	if (File.is_open()) {
		File << macroscopicStrain << " " << macroscopicStress << endl;
		return true;
	}
	return false;
}

/*
 * calculates yieldingpoint based on last loadpath in stress_strain_file
 * interpolates linearly between given points
 * returns a vector with yieldStress, yieldStrain, yieldPlasticStrain, yieldPlasticStrainRate
 */
template<int dim>
vector<SymmetricTensor<2, dim>> getYieldingPoint(ParameterHandler &para,
		double &yieldDissipation) {
	vector<double> dissipation;

	// filename
	string dir = para.get("straindirection");
	std::stringstream sstm;
	sstm << "Results/stress_strain_" << dir << ".txt";
	string filename = sstm.str();

	vector<vector<SymmetricTensor<2, dim> > > loadpath;
	ifstream myReadFile;
	myReadFile.open(filename);
	string output;
	if (myReadFile.is_open()) {
		while (!myReadFile.eof()) {
			getline(myReadFile, output);
			if (output.find("Begin step 0") != string::npos) {
				loadpath.clear();
				dissipation.clear();
			}

			if (output.find("Begin step") != string::npos) {
				vector<SymmetricTensor<2, dim>> tmp;
				getline(myReadFile, output); // stress
				tmp.push_back(readTensorLine<dim>(output));

				getline(myReadFile, output); // macroStrain from boundary object
				getline(myReadFile, output); // strain mean in rve
				tmp.push_back(readTensorLine<dim>(output));

				getline(myReadFile, output); // plastic_strain_mean
				tmp.push_back(readTensorLine<dim>(output));

				getline(myReadFile, output); // plastic_strain_rate_mean
				tmp.push_back(readTensorLine<dim>(output));

				getline(myReadFile, output); // dissipation value, delta in actual step
				double dissipationValue = readDoubleLine(output);
				if (dissipation.size() > 0) {
					dissipationValue += dissipation[dissipation.size() - 1];
				}
				dissipation.push_back(dissipationValue);

				loadpath.push_back(tmp);
			}
		}
	}
	myReadFile.close();

	// search yieldpoint
	double lastYieldCriterion = 0;
	//	deallog << "Yieldcriterion:";
	for (int i = 1; i < loadpath.size(); ++i) {
		SymmetricTensor<2, dim> macroStress = loadpath[i][0];
		SymmetricTensor<2, dim> macroPlasticStrain = loadpath[i][2];
		double yieldCriterion;
		if (yieldDissipation < 0) {
			yieldCriterion = macroStress * macroPlasticStrain
					/ macroStress.norm();
		} else {
			yieldCriterion = dissipation[i - 1];
		}

		//		SymmetricTensor<2, dim> macroStrainRate = loadpath[i][1]
		//				- loadpath[i - 1][1];
		//		SymmetricTensor<2, dim> macroPlasticStrainRate = loadpath[i][2]
		//				- loadpath[i - 1][2];
		//		double yieldCriterion = (macroStrainRate * macroPlasticStrainRate)
		//				/ (macroStrainRate * macroStrainRate);
		//		deallog << " " << yieldCriterion;

		double criticalYC;
		if (yieldDissipation < 0) {
			criticalYC = para.get_double("criticalYC");
		} else {
			criticalYC = yieldDissipation;
		}

		if (yieldCriterion >= criticalYC) {
			double fraction = (criticalYC - lastYieldCriterion)
					/ (yieldCriterion - lastYieldCriterion);
			vector<SymmetricTensor<2, dim>> tmp;
			for (int j = 0; j < loadpath[i].size(); ++j) {
				tmp.push_back(
						loadpath[i - 1][j]
								+ fraction
										* (loadpath[i][j] - loadpath[i - 1][j]));
			}
			//			deallog << "Macro stress: " << tmp[0] << endl;
			//			deallog << "Macro plastic strain " << tmp[2] << endl;
			if (yieldDissipation < 0) {
				yieldDissipation = dissipation[i - 2]
						+ fraction * (dissipation[i - 1] - dissipation[i - 2]);
			}
			return tmp;
		}

		lastYieldCriterion = yieldCriterion;
	}
	//	deallog << endl;
	vector<SymmetricTensor<2, dim>> tmp;
	//	deallog << "Last Yieldcriterion: " << lastYieldCriterion << endl;
	return tmp;
}

/*
 * update solution_history, is needed for mechanical prediction,
 * dynamic relaxation and other things
 */
void update_solutionhistory(vector<Vector<double>> &solution_history,
		Vector<double> &solution, vector<double> &times_history, double time,
		ParameterHandler &parameter) {
	// insert last solutiondelta
	solution_history.insert(solution_history.begin(), solution);
	times_history.insert(times_history.begin(), time);

	// necessary number of old solutiondeltas
	// for mechanical prediction
	unsigned int necessary_length = parameter.get_integer("MP_order") + 1;

	//	deallog << "Length of solution history before removing: "
	//			<< solution_history.size() << endl;
	// remove no longer needed solutions from history
	while (solution_history.size() > necessary_length) {
		solution_history.pop_back();
	}
	while (times_history.size() > necessary_length + 1) {
		times_history.pop_back();
	}
	//	deallog << "Length of solution history after  removing: "
	//			<< solution_history.size() << endl;

}

/*
 * 	Polynomial extrapolation according to
 An adaptive polynomial based forward prediction algorithm for
 multi-actuator real-time dynamic substructuring
 M.I Wallace, D.J Wagg, S.A Neild

 returns a vector with coefficients for polynomial extrapolation with maximal polynomial degree
 for solutions given at the points times_history
 -> polydegree = times_history.size() - 2
 */
Vector<double> get_MP_coefficients(vector<double> &times_history) {
	if (times_history.size() == 1) {
		return Vector<double>(0);
	}
	/*
	 * first entry is the time of the next timestep
	 * has to be erased from history
	 */
	double t_star = times_history[0];
	int polyDegree = times_history.size() - 2;

	/* Assemble X Matrix from eq 3.2, 3.3
	 *
	 * n = times_history.size() - 1
	 t_1^0	t_1^1 ... t_1^polyDegree
	 t_2^0	t_2^1 ... t_2^polyDegree
	 ...
	 t_n^0	t_n^1 ... t_n^polyDegree
	 */
	FullMatrix<double> X(times_history.size() - 1, polyDegree + 1);
	for (unsigned int i = 0; i < X.m(); ++i) {
		for (unsigned int j = 0; j < X.n(); ++j) {
			X(i, j) = pow(times_history[i + 1], j);
		}
	}

	/*
	 * Assemble X_p eq. 3.11, like additional line in eq. 3.3
	 * [t_star^0, t_star^1, ..., t_star^polyDegree
	 */
	Vector<double> X_p(polyDegree + 1);
	for (unsigned int i = 0; i < X_p.size(); ++i) {
		X_p[i] = pow(t_star, i);
	}

	// The calculation being performed is (X'*X)-1 *X'.
	FullMatrix<double> X_txX_inv = X;
	if (X.m() == 1) {
		X_txX_inv(0, 0) = 1 / X(0, 0);
	} else {
		X_txX_inv.invert(X);
	}

	Vector<double> a(X_p.size());
	X_txX_inv.Tvmult(a, X_p, false);
	return a;

}

template<int dim>
void writeOutSolutionVector(Vector<double> &solution,
		hp::DoFHandler<dim> dof_handler) {
	// write stress strain to file
	std::stringstream sstm;
	sstm << "Results/SolutionVectors.txt";
	string filename = sstm.str();

	std::ofstream File(filename, std::ios::out | std::ios::app);
	//	if (File.is_open()) {
	//		solution.print(File);
	//	}

}

template<int dim>
void print(std::ostream &o, const Tensor<2, dim, double> &t) {
	int w = 14;
	for (unsigned int i = 0; i < dim; ++i) {
		for (unsigned int j = 0; j < dim; ++j) {
			o << setw(w) << t[i][j];
		}
		o << endl;
	}
}

template<int dim>
void print(std::ostream &o, const SymmetricTensor<2, dim, double> &t) {
	int w = 14;
	for (unsigned int i = 0; i < dim; ++i) {
		for (unsigned int j = 0; j < dim; ++j) {
			o << setw(w) << t[i][j];
		}
		o << endl;
	}
}

/*
 * Attention, is just for a Tensor of rank 1 with
 * [sig_xx, sig_yy, sig_zz] and sig_xy = sig_xz = sig_yz = 0
 */
template<int dim>
Tensor<1, dim> getDeviatoricPart(Tensor<1, dim> &inp) {
	Tensor<1, dim> retVal = inp;
	for (int i = 0; i < dim; ++i) {
		for (int j = 0; j < dim; ++j) {
			retVal[j] -= inp[i] / dim;
		}
	}
	return retVal;
}

/*
 * create macroscopic tangent from parameterfile entry
 */
template<int dim>
SymmetricTensor<4, dim> getMacroTangentFromParameterfile(
		ParameterHandler &parameter) {
	string inp = parameter.get("macroTangent");
	SymmetricTensor<4, dim> retval;

	const std::vector<std::string> splitted_entries =
			Utilities::split_string_list(inp, ',');

	if (splitted_entries.size() != 36) {
		return retval;
	}

	Tensor<2, dim> normalEntries;
	Tensor<2, dim> couplingEntries;
	Tensor<2, dim> shearEntries;

	for (unsigned int i = 0; i < splitted_entries.size(); ++i) {
		int row = i / 6;
		int column = i % 6;

		if (row < dim && column < dim) {
			normalEntries[row % dim][column % dim] = std::stod(
					splitted_entries[i]);
		} else if (row < dim && column >= dim) {
			couplingEntries[row % dim][column % dim] = std::stod(
					splitted_entries[i]);
		} else if (row >= dim && column >= dim) {
			shearEntries[row % dim][column % dim] = std::stod(
					splitted_entries[i]);
		}
	}

	for (unsigned int i = 0; i < dim; ++i) {
		for (unsigned int j = 0; j < dim; ++j) {
			retval[i][i][j][j] = normalEntries[i][j];
		}
	}

	for (unsigned int i = 0; i < dim; ++i) {
		for (unsigned int j = 0; j < dim; ++j) {
			retval[i][i][(j + 1) % dim][(j + 2) % dim] = couplingEntries[i][j];
			retval[i][i][(j + 2) % dim][(j + 1) % dim] = couplingEntries[i][j];

			retval[(j + 2) % dim][(j + 1) % dim][i][i] = couplingEntries[i][j];
			retval[(j + 1) % dim][(j + 2) % dim][i][i] = couplingEntries[i][j];
		}
	}

	for (unsigned int i = 0; i < dim; ++i) {
		for (unsigned int j = 0; j < dim; ++j) {
			int k = (i + 1) % dim;
			int l = (i + 2) % dim;
			int m = (j + 1) % dim;
			int n = (j + 2) % dim;

			retval[k][l][m][n] = shearEntries[i][j];
			retval[l][k][m][n] = shearEntries[i][j];
			retval[k][l][n][m] = shearEntries[i][j];
			retval[l][k][n][m] = shearEntries[i][j];
		}
	}

	return retval;
}

template<int dim>
SymmetricTensor<2, dim> getUnitLoad(string load) {
	SymmetricTensor<2, dim> retval;

	if (load == "11") {
		retval[0][0] = 1;
	} else if (load == "22") {
		retval[1][1] = 1;
	} else if (load == "33") {
		retval[2][2] = 1;
	} else if (load == "12") {
		retval[0][1] = 0.5;
	} else if (load == "13") {
		retval[0][2] = 0.5;
	} else if (load == "23") {
		retval[1][2] = 0.5;
	}

	return retval;

}

#endif /* CODE_VARIOUS_AUXILARYFUNCTIONS_H_ */
