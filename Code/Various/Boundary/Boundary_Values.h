/*
 * Boundary_Values.h
 *
 *  Created on: May 31, 2017
 *      Author: iwtm84
 */

#ifndef CODE_VARIOUS_BOUNDARY_VALUES_H_
#define CODE_VARIOUS_BOUNDARY_VALUES_H_

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

using namespace dealii;

template<int dim>
class BoundaryValues: public Function<dim> {
public:
	BoundaryValues();
	virtual ~BoundaryValues() {

	}

	virtual void vector_value_list(const std::vector<Point<dim> > &points,
			std::vector<Vector<double> > &value_list) const = 0;

	// values at point p
	virtual void vector_value(const Point<dim> &p,
			Vector<double> &values) const = 0;

	// init class and set timestep to the first step
	virtual void first_timestep() = 0;

	// indicate new timestep, load will be applied in next step
	virtual void next_timestep() = 0;

	// indicate that first newtonstep is done, in the next step the load will not be applied
	virtual void deactivate_inhomogenious_bc() = 0;

	// get actual full macroscopic strain in actual loadstep
	virtual SymmetricTensor<2, dim> get_full_macroscopic_strain() const = 0;

	// get delta x between last and actual loadstep
	virtual double get_dx() const = 0;

	// get delta t between last and actual loadstep
	virtual double get_dt() const = 0;

	// get theta in actual loadstep
	virtual double get_theta(const Point<dim> &p) const = 0;

	// get initial theta, boundary object will be initialized here
	virtual double get_initial_theta(const Point<dim> &p) const = 0;

	// returns false if last loadstep is reached, otherwise true
	virtual bool hasNextIncrement() const = 0;

	// returns false if last loadstep is reached, otherwise true
	virtual double get_time() const = 0;

	// returns true if the loaddirection of the acutal loadstep is equal to the one before
	virtual bool isMonotonicLoad() const = 0;

//	// sets the macroscopic strain
//	virtual void setNormalizedMacroscopicStrain(
//			SymmetricTensor<2, dim> macroscopic_strain) const = 0;

};

template<int dim>
BoundaryValues<dim>::BoundaryValues() :
		Function<dim>(dim) {
}

/*
 * ----------------------------------------------------
 * -------------- Boundary Values Linear --------------
 * ----------------------------------------------------
 */

template<int dim>
class BoundaryValuesLinear: public BoundaryValues<dim> {
public:
	BoundaryValuesLinear(ParameterHandler *para);
	~BoundaryValuesLinear() {
	}
	;

	void vector_value_list(const std::vector<Point<dim> > &points,
			std::vector<Vector<double> > &value_list) const;

	// values at point p
	void vector_value(const Point<dim> &p, Vector<double> &values) const;

	// init class and set timestep to the first step
	void first_timestep();

	// indicate new timestep, load will be applied in next step
	void next_timestep();

	// indicate that first newtonstep is done, in the next step the load will not be applied
	void deactivate_inhomogenious_bc();

	// get actual full macroscopic strain in actual loadstep
	SymmetricTensor<2, dim> get_full_macroscopic_strain() const;

	// get delta x between last and actual loadstep
	double get_dx() const;

	// get delta t between last and actual loadstep
	double get_dt() const;

	// get theta in actual loadstep
	double get_theta(const Point<dim> &p) const;

	// get initial theta, boundary object will be initialized here
	double get_initial_theta(const Point<dim> &p) const;

	// returns false if last loadstep is reached, otherwise true
	bool hasNextIncrement() const;

	//returns time
	double get_time() const;

	// returns true if the loaddirection of the acutal loadstep is equal to the one before
	bool isMonotonicLoad() const {
		return true;
	}

	// sets the normalized macroscopic strain, should be in the linear strain class
	void setNormalizedMacroscopicStrain(
			SymmetricTensor<2, dim> macroscopic_strain) {
		macro_strain = macroscopic_strain / macroscopic_strain.norm();
	}

private:
	// get actual full macroscopic strain in actual loadstep
	SymmetricTensor<2, dim> get_macroscopic_strain() const;
	double get_displacement(double t) const;
	double get_duration() const;
	double get_theta(const double t, const double z) const;

	ParameterHandler *parameter;
	SymmetricTensor<2, dim> macro_strain;
	vector<double> times_loadsteps;
	bool load_active;

	/*
	 * Range is [0, n_loadsteps]
	 */
	unsigned int loadstep;
};

template<int dim>
BoundaryValuesLinear<dim>::BoundaryValuesLinear(ParameterHandler *para) :
		BoundaryValues<dim>(), load_active(false), loadstep(0) {
	parameter = para;
	macro_strain = SymmetricTensor<2, dim>();
	times_loadsteps = vector<double>(0);
	load_active = false;
	loadstep = -1;
}

/*
 * uses vector_value to calculate values at all points in vector points
 */
template<int dim>
void BoundaryValuesLinear<dim>::vector_value_list(
		const std::vector<Point<dim> > &points,
		std::vector<Vector<double> > &value_list) const {
	const unsigned int n_points = points.size();
	Assert(value_list.size() == n_points,
			ExcDimensionMismatch(value_list.size(), n_points));
	for (unsigned int p = 0; p < n_points; ++p) {
		vector_value(points[p], value_list[p]);
	}
}

/*
 * calculates value at point p
 */
template<int dim>
void BoundaryValuesLinear<dim>::vector_value(const Point<dim> &p,
		Vector<double> &values) const {
	Assert(values.size() == dim, ExcDimensionMismatch(values.size(), dim));
	SymmetricTensor<2, dim> macroscopic_strain = get_macroscopic_strain();
	Tensor<1, dim> tmp = macroscopic_strain * p;
	for (unsigned int i = 0; i < dim; i++) {
		values[i] = tmp[i];
	}
}

template<int dim>
void BoundaryValuesLinear<dim>::first_timestep() {
	// Initialize times_loadsteps for whole loadpath
	times_loadsteps.clear();

	for (unsigned int i = 0; i <= parameter->get_double("n_loadsteps"); ++i) {
		double tmp = get_duration() / parameter->get_double("n_loadsteps") * i;
		times_loadsteps.push_back(tmp);
	}
	loadstep = 0;
	load_active = false;

	// Initialize macroscopic strain
	if (parameter->get("straindirection").compare("detailed") == 0) {
		if (macro_strain.norm() < 1e-8) {
			deallog << "macro_strain not set for straindirection detailed"
					<< endl;
		}
	} else {
		macro_strain = 0;
		if (parameter->get("straindirection").compare("11") == 0) {
			macro_strain[0][0] = 1;
		} else if (parameter->get("straindirection").compare("22") == 0) {
			macro_strain[1][1] = 1;
		} else if (parameter->get("straindirection").compare("33") == 0) {
			macro_strain[2][2] = 1;
		} else if (parameter->get("straindirection").compare("12") == 0) {
			macro_strain[0][1] = 0.5;
		} else if (parameter->get("straindirection").compare("13") == 0) {
			macro_strain[0][2] = 0.5;
		} else if (parameter->get("straindirection").compare("23") == 0) {
			macro_strain[1][2] = 0.5;
		}
	}
}

template<int dim>
void BoundaryValuesLinear<dim>::next_timestep() {
	load_active = true;
	loadstep++;
}

template<int dim>
void BoundaryValuesLinear<dim>::deactivate_inhomogenious_bc() {
	load_active = false;
}

/*
 * get full macroscopic strain at actual loadstep
 */
template<int dim>
SymmetricTensor<2, dim> BoundaryValuesLinear<dim>::get_full_macroscopic_strain() const {
	SymmetricTensor<2, dim> tmp = macro_strain
			* get_displacement(times_loadsteps[loadstep]);

	return tmp;
}

/*
 * get dx of actual timestep, if load inactive dx = 0
 */
template<int dim>
double BoundaryValuesLinear<dim>::get_dx() const {
	Assert(loadstep > 0 && loadstep < times_loadsteps.size(),
			ExcMessage("Loadstep out of Range"));

	if (load_active) {
		double t_1 = times_loadsteps[loadstep - 1];
		double t_2 = times_loadsteps[loadstep];

		return get_displacement(t_2) - get_displacement(t_1);
	}
	return 0;
}

/*
 * get delta t of actual timestep
 */
template<int dim>
double BoundaryValuesLinear<dim>::get_dt() const {
	if (!(loadstep > 0 && loadstep < times_loadsteps.size())) {
		return 1;
	}
	double t_n_1 = times_loadsteps[loadstep - 1];
	double t_n = times_loadsteps[loadstep];
	double dt = t_n - t_n_1;
	return dt;
}

/*
 * returns temperature at height z and time t
 */
template<int dim>
double BoundaryValuesLinear<dim>::get_theta(const double t,
		const double z) const {
	double theta_top = parameter->get_double("t_start")
			+ t / get_duration()
					* (parameter->get_double("t_end")
							- parameter->get_double("t_start"));
	double theta_bottom = parameter->get_double("t_end");
	return theta_bottom + z * (theta_top - parameter->get_double("t_end"));
}

/*
 * returns temperature at point p in actual timestep
 */
template<int dim>
double BoundaryValuesLinear<dim>::get_theta(const Point<dim>& p) const {
	return get_theta(times_loadsteps[loadstep], p[2]);
}

/*
 * returns temperature at point p in initial timestep
 */
template<int dim>
double BoundaryValuesLinear<dim>::get_initial_theta(const Point<dim> &p) const {
	return get_theta(0, p[2]);
}

/*
 * returns true if
 */
template<int dim>
bool BoundaryValuesLinear<dim>::hasNextIncrement() const {
	return loadstep < parameter->get_double("n_loadsteps");
}

/*
 * return time of actual timestep
 */
template<int dim>
double BoundaryValuesLinear<dim>::get_time() const {
	if (times_loadsteps.size() == 0) {
		return 0;
	}
	return times_loadsteps[loadstep];
}

/*
 * get actual delta in macroscopic strain, if load inactive delta = 0
 */
template<int dim>
SymmetricTensor<2, dim> BoundaryValuesLinear<dim>::get_macroscopic_strain() const {
	SymmetricTensor<2, dim> tmp;
	if (load_active) {
		tmp = macro_strain * get_dx();
	}
	return tmp;
}

/*
 * returns the displacement at time t
 */
template<int dim>
double BoundaryValuesLinear<dim>::get_displacement(double t) const {
	if (parameter->get_double("displacement") == 0) {
		return 0;
	}
	return parameter->get_double("displacement") * t / get_duration();
}

/*
 * returns the over all duration with respect to displacement and strainrate
 */
template<int dim>
double BoundaryValuesLinear<dim>::get_duration() const {
	double duration = parameter->get_double("displacement")
			/ parameter->get_double("strainrate");
	if (duration < 1e-16) {
		return 1;
	}
	return duration;
}

/*
 * ----------------------------------------------------
 * ------- Boundary Values Prescribed form File -------
 * ----------------------------------------------------
 */

template<int dim>
class BoundaryValuesPrescribed: public BoundaryValues<dim> {
public:
	BoundaryValuesPrescribed(ParameterHandler *para);
	~BoundaryValuesPrescribed() {
	}
	;

	void vector_value_list(const std::vector<Point<dim> > &points,
			std::vector<Vector<double> > &value_list) const;

	// values at point p
	void vector_value(const Point<dim> &p, Vector<double> &values) const;

	// init class and set timestep to the first step
	void first_timestep();

	// indicate new timestep, load will be applied in next step
	void next_timestep();

	// indicate that first newtonstep is done, in the next step the load will not be applied
	void deactivate_inhomogenious_bc();

	// get actual full macroscopic strain in actual loadstep
	SymmetricTensor<2, dim> get_full_macroscopic_strain() const;

	// get delta x between last and actual loadstep
	double get_dx() const;

	// get delta t between last and actual loadstep
	double get_dt() const;

	// get theta in actual loadstep
	double get_theta(const Point<dim> &p) const;

	// get initial theta, boundary object will be initialized here
	double get_initial_theta(const Point<dim> &p) const;

	// returns false if last loadstep is reached, otherwise true
	bool hasNextIncrement() const;

	//returns time
	double get_time() const;

	// returns true if the loaddirection of the acutal loadstep is equal to the one before
	bool isMonotonicLoad() const;

	// sets the normalized macroscopic strain, should be in the linear strain class
	void setNormalizedMacroscopicStrain(
			SymmetricTensor<2, dim> macroscopic_strain) {
		deallog
				<< "The function setNormalizedMacroscopicStrain should never been called for BoundaryValuesPrescribed"
				<< endl;
		exit(0);
	}

private:
	// get actual full macroscopic strain in actual loadstep
	SymmetricTensor<2, dim> get_macroscopic_strain() const;
	double get_theta(const double t, const double z) const;

	ParameterHandler *parameter;
	vector<double> times_loadsteps;
	vector<double> temperatures_loadsteps;
	vector<SymmetricTensor<2, dim>> strains_loadsteps;
	bool load_active;

	/*
	 * Range is [0, times_loadsteps.size()-1]
	 */
	unsigned int loadstep;
};

template<int dim>
BoundaryValuesPrescribed<dim>::BoundaryValuesPrescribed(ParameterHandler *para) :
		BoundaryValues<dim>(), load_active(false), loadstep(0) {
	parameter = para;
	times_loadsteps.clear();
	temperatures_loadsteps.clear();
	strains_loadsteps.clear();

	times_loadsteps = loadDoubleArchive("Input/Loadpath/times");
	temperatures_loadsteps = loadDoubleArchive("Input/Loadpath/temperatures");
	strains_loadsteps = loadSymmetricTensorArchive<dim>(
			"Input/Loadpath/strains");
}

/*
 * uses vector_value to calculate values at all points in vector points
 */
template<int dim>
void BoundaryValuesPrescribed<dim>::vector_value_list(
		const std::vector<Point<dim> > &points,
		std::vector<Vector<double> > &value_list) const {
	const unsigned int n_points = points.size();
	Assert(value_list.size() == n_points,
			ExcDimensionMismatch(value_list.size(), n_points));
	for (unsigned int p = 0; p < n_points; ++p) {
		vector_value(points[p], value_list[p]);
	}
}

/*
 * calculates value at point p
 */
template<int dim>
void BoundaryValuesPrescribed<dim>::vector_value(const Point<dim> &p,
		Vector<double> &values) const {
	Assert(values.size() == dim, ExcDimensionMismatch(values.size(), dim));
	SymmetricTensor<2, dim> macroscopic_strain = get_macroscopic_strain();
	Tensor<1, dim> tmp = macroscopic_strain * p;
	for (unsigned int i = 0; i < dim; i++) {
		values[i] = tmp[i];
	}
}

template<int dim>
void BoundaryValuesPrescribed<dim>::first_timestep() {
	// Initialize times_loadsteps for whole loadpath
	times_loadsteps.clear();
	temperatures_loadsteps.clear();
	strains_loadsteps.clear();
	loadstep = 0;
	load_active = false;

	times_loadsteps = loadDoubleArchive("Input/Loadpath/times");
	temperatures_loadsteps = loadDoubleArchive("Input/Loadpath/temperatures");
	strains_loadsteps = loadSymmetricTensorArchive<dim>(
			"Input/Loadpath/strains");

	if (times_loadsteps.size() != temperatures_loadsteps.size()
			|| times_loadsteps.size() != strains_loadsteps.size()) {
		deallog
				<< "Times, Temperatures, Strains with different numbers of timesteps."
				<< endl;
		exit(0);
	}

//	for (unsigned int i = 0; i < times_loadsteps.size(); ++i) {
//		deallog << "Time: " << setw(4) << times_loadsteps[i] << " Temperature: "
//				<< setw(4) << temperatures_loadsteps[i] << " Strain: "
//				<< strains_loadsteps[i] << endl;
//	}

}

template<int dim>
void BoundaryValuesPrescribed<dim>::next_timestep() {
	load_active = true;
	loadstep++;
}

template<int dim>
void BoundaryValuesPrescribed<dim>::deactivate_inhomogenious_bc() {
	load_active = false;
}

/*
 * get full macroscopic strain at actual loadstep
 */
template<int dim>
SymmetricTensor<2, dim> BoundaryValuesPrescribed<dim>::get_full_macroscopic_strain() const {
	if (loadstep < 0 || loadstep >= strains_loadsteps.size()) {
		Assert(loadstep > 0 && loadstep < times_loadsteps.size(),
				ExcMessage("Loadstep out of Range"));
	}
	SymmetricTensor<2, dim> tmp = strains_loadsteps[loadstep];
	return tmp;
}

/*
 * get dx of actual timestep, if load inactive dx = 0
 */
template<int dim>
double BoundaryValuesPrescribed<dim>::get_dx() const {
	if (parameter->get("loadtype").compare("testcase_andrew_microhard") == 0
			|| parameter->get("loadtype").compare("testcase_grain_boundary")
					== 0) {
		if (load_active) {
			return (strains_loadsteps[loadstep][0][1]
					- strains_loadsteps[loadstep - 1][0][1]) * 2;
		}
		return 0;
	}
	deallog << "The function get_dx should only been called for the testcases"
			<< endl;
	exit(0);
}

/*
 * get delta t of actual timestep
 */
template<int dim>
double BoundaryValuesPrescribed<dim>::get_dt() const {
	if (loadstep < 1){
		return 1;
	}

	double t_n_1 = times_loadsteps[loadstep - 1];
	double t_n = times_loadsteps[loadstep];
	double dt = t_n - t_n_1;
	if (dt < 1e-16) {
		deallog << "dt = " << dt << " returning 1e-16" << endl;
		return 1e-16;
	}
	return t_n - t_n_1;
}

/*
 * returns temperature at height z and time t
 */
template<int dim>
double BoundaryValuesPrescribed<dim>::get_theta(const double t,
		const double /*z*/) const {
	Assert(loadstep > 0 && loadstep < temperatures_loadsteps.size(),
			ExcMessage("Loadstep out of Range"));
	return temperatures_loadsteps[loadstep];
}

/*
 * returns temperature at point p in actual timestep
 */
template<int dim>
double BoundaryValuesPrescribed<dim>::get_theta(const Point<dim>& /*p*/) const {
	Assert(loadstep > 0 && loadstep < temperatures_loadsteps.size(),
			ExcMessage("Loadstep out of Range"));
	return temperatures_loadsteps[loadstep];
}

/*
 * returns temperature at point p in initial timestep
 */
template<int dim>
double BoundaryValuesPrescribed<dim>::get_initial_theta(
		const Point<dim> &/*p*/) const {
	return temperatures_loadsteps[0];
}

/*
 * returns true if
 */
template<int dim>
bool BoundaryValuesPrescribed<dim>::hasNextIncrement() const {
	return loadstep < times_loadsteps.size() - 1;
}

/*
 * return time of actual timestep
 */
template<int dim>
double BoundaryValuesPrescribed<dim>::get_time() const {
	if (!(loadstep > 0 && loadstep < times_loadsteps.size())) {
		return 1;
	}
	return times_loadsteps[loadstep];
}

/*
 * get actual delta in macroscopic strain, if load inactive delta = 0
 */
template<int dim>
SymmetricTensor<2, dim> BoundaryValuesPrescribed<dim>::get_macroscopic_strain() const {
	SymmetricTensor<2, dim> tmp;
	if (load_active) {
		tmp = strains_loadsteps[loadstep];
	}
	return tmp;
}

template<int dim>
bool BoundaryValuesPrescribed<dim>::isMonotonicLoad() const {
	if (loadstep <= 1) {
		return true;
	}
	SymmetricTensor<2, dim> lastStep = strains_loadsteps[loadstep - 2];
	lastStep -= strains_loadsteps[loadstep - 1];
	SymmetricTensor<2, dim> nextStep = strains_loadsteps[loadstep - 1];
	nextStep -= strains_loadsteps[loadstep];

	double tmp1 = lastStep * nextStep;
	double tmp2 = lastStep * lastStep;
	double tmp3 = nextStep * nextStep;

	double result = tmp1 / pow(tmp2 * tmp3, 0.5);
//	deallog << "Result monotonic load: " << result << endl;
	return abs(result - 1.0) < 1e-8;

}

#endif /* CODE_VARIOUS_BOUNDARY_VALUES_H_ */
