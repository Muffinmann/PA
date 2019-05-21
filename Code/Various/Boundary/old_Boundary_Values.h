/*
 * Boundary_Values.h
 *
 *  Created on: Jun 30, 2016
 *      Author: iwtm84
 */

#ifndef CODE_VARIOUS_OLD_BOUNDARY_VALUES_H_
#define CODE_VARIOUS_OLD_BOUNDARY_VALUES_H_

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

using namespace dealii;

template<int dim>
class BoundaryValues: public Function<dim> {
public:
	BoundaryValues(ParameterHandler *para);
//	BoundaryValues(const BoundaryValues<dim>& source);
	BoundaryValues<dim> &operator=(const BoundaryValues<dim>& source);
	virtual ~BoundaryValues() {

	}
	void vector_value(const Point<dim> &p, Vector<double> &values) const;
	void vector_value(const Point<dim> &p, Vector<double> &values,
			bool apply_load);
	virtual
	void vector_value_list(const std::vector<Point<dim> > &points,
			std::vector<Vector<double> > &value_list) const;
	void init();

	SymmetricTensor<2, dim> get_full_macroscopic_strain() const;
	void next_timestep();
	void first_timestep();
	void deactivate_inhomogenious_bc();
	double get_dt() const;
	double get_dx() const;
	bool apply_load() const;
	double get_theta(const Point<dim> &p) const;
	double get_initial_theta(const Point<dim> &p) const;
	double apply_bc_via_newtonguess(double factor);
	bool hasNextIncrement();
	int NumberOfLoadsteps();
	double get_time();

	void setMacroscopicStrain(SymmetricTensor<2, dim> macroscopic_strain) {
		macro_strain = macroscopic_strain / macroscopic_strain.norm();
	}
private:
	double get_displacement(double t) const;
	double get_theta(double t, double z) const;
	SymmetricTensor<2, dim> get_macroscopic_strain() const;
	ParameterHandler *para;
	bool load;
	vector<double> time_loadsteps;
	double propotion_applied_via_newtonguess;

	vector<double> times_macro;
	vector<double> temperatures_macro;
	vector<SymmetricTensor<2, dim> > strains_macro;
	int timestep;

	enum modes {
		unknown = -1, prescribed_from_file = 0, linear_displacement = 1
	};
	modes mode;

	SymmetricTensor<2, dim> macro_strain;

	double duration() const {
		return para->get_double("displacement") / para->get_double("strainrate");
	}
};

template<int dim>
BoundaryValues<dim>::BoundaryValues(ParameterHandler *para) :
		Function<dim>(dim), para(para), load(true), propotion_applied_via_newtonguess(
				0), timestep(0), mode(unknown), macro_strain() {
}

//template<int dim>
//BoundaryValues<dim>::BoundaryValues(BoundaryValues<dim>& source) {
//	 			Function<dim>(dim);
//	 			para(source.para);
//	 			load(source.load);
//	 			propotion_applied_via_newtonguess(source.propotion_applied_via_newtonguess);
//	 			timestep(source.timestep);
//	 			mode(source.mode);
//	 			macro_strain(source.macro_strain);
//}

template<int dim>
BoundaryValues<dim> &BoundaryValues<dim>::operator=(
		const BoundaryValues<dim>& source) {
	if (this == &source) {
		return *this;
	}
	Function<dim>::operator =(source);
	para = source.para;
	load = source.load;
	propotion_applied_via_newtonguess = source.propotion_applied_via_newtonguess;
	timestep = source.timestep;
	mode = source.mode;
	macro_strain = source.macro_strain;

	return *this;
}

template<int dim>
void BoundaryValues<dim>::init() {
	time_loadsteps.clear();
	load = true;
	propotion_applied_via_newtonguess = 0;
	timestep = 0;
	if (para->get("strainpath").compare("prescribed_from_file") == 0) {
		mode = prescribed_from_file;

		std::ifstream ifstream_boost_archive_temp(
				"Input/data_macrocalculation/macro_temperatures.txt");
		std::ifstream ifstream_boost_archive_time(
				"Input/data_macrocalculation/macro_times.txt");
		std::ifstream ifstream_boost_archive_strain(
				"Input/data_macrocalculation/macro_strains.txt");

		if (!ifstream_boost_archive_temp.is_open()) {
			deallog << "Trouble with temp archive" << endl;
		}
		if (!ifstream_boost_archive_time.is_open()) {
			deallog << "Trouble with time archive" << endl;
		}
		if (!ifstream_boost_archive_strain.is_open()) {
			deallog << "Trouble with strain archive" << endl;
		}

		boost::archive::text_iarchive ia_temp(ifstream_boost_archive_temp);
		boost::archive::text_iarchive ia_time(ifstream_boost_archive_time);
		boost::archive::text_iarchive ia_strain(ifstream_boost_archive_strain);

		vector<FullMatrix<double> > strains_macro_FM;
		ia_temp >> temperatures_macro;
		ia_time >> times_macro;
		ia_strain >> strains_macro_FM;

		for (int i = 0; i < strains_macro_FM.size(); i++) {
			Tensor<2, dim> strain_T;
			strains_macro_FM[i].copy_to(strain_T);
			strains_macro.push_back(strain_T);
		}

		ifstream_boost_archive_temp.close();
		ifstream_boost_archive_time.close();
		ifstream_boost_archive_strain.close();

//		for (int i=0;i<temperatures_macro.size();i++){
//			cout << "Time: " << times_macro[i] << " Temperature: " << temperatures_macro[i] << " Strain: " << strains_macro[i] << endl;
//		}
	} else if (para->get("strainpath").compare("linear_displacement") == 0) {
		mode = linear_displacement;

	} else {
		mode = unknown;
		cerr << "Unknown loadtype (Constructor)" << endl;
//		exit(0);
	}
}

template<int dim>
void BoundaryValues<dim>::vector_value(const Point<dim> &p,
		Vector<double> &values) const {
	Assert(values.size() == dim, ExcDimensionMismatch(values.size(), dim));
	SymmetricTensor<2, dim> macroscopic_strain = get_macroscopic_strain();
	Tensor<1, dim> tmp = macroscopic_strain * p;
	values[0] = tmp[0];
	values[1] = tmp[1];
	values[2] = tmp[2];
}

template<int dim>
void BoundaryValues<dim>::vector_value(const Point<dim> &p,
		Vector<double> &values, bool apply_load) {
	bool tmp_load = load;
	if (apply_load) {
		load = true;
	}
	Assert(values.size() == dim, ExcDimensionMismatch(values.size(), dim));
	SymmetricTensor<2, dim> macroscopic_strain = get_macroscopic_strain();
	Tensor<1, dim> tmp = macroscopic_strain * p;
	values[0] = tmp[0];
	values[1] = tmp[1];
	values[2] = tmp[2];
	if (apply_load) {
		load = tmp_load;
	}
}

template<int dim>
void BoundaryValues<dim>::vector_value_list(
		const std::vector<Point<dim> > &points,
		std::vector<Vector<double> > &value_list) const {
	const unsigned int n_points = points.size();
	Assert(value_list.size() == n_points,
			ExcDimensionMismatch(value_list.size(), n_points));
	for (unsigned int p = 0; p < n_points; ++p) {
		vector_value(points[p], value_list[p]);
	}
}
template<int dim>
SymmetricTensor<2, dim> BoundaryValues<dim>::get_macroscopic_strain() const {
	SymmetricTensor<2, dim> macroscopic_strain;

	if (mode == linear_displacement) {
		double displacementvalue = get_dx();
		if (para->get("straindirection").compare("11") == 0) {
			macroscopic_strain[0][0] = displacementvalue;
		} else if (para->get("straindirection").compare("22") == 0) {
			macroscopic_strain[1][1] = displacementvalue;
		} else if (para->get("straindirection").compare("33") == 0) {
			macroscopic_strain[2][2] = displacementvalue;
		} else if (para->get("straindirection").compare("12") == 0) {
			macroscopic_strain[0][1] = displacementvalue / 2;
		} else if (para->get("straindirection").compare("13") == 0) {
			macroscopic_strain[0][2] = displacementvalue / 2;
		} else if (para->get("straindirection").compare("23") == 0) {
			macroscopic_strain[1][2] = displacementvalue / 2;
		} else if (para->get("straindirection").compare("detailed") == 0) {
			macroscopic_strain = macro_strain * displacementvalue;
		}
	} else if (mode == prescribed_from_file) {
		// only if load, else 0
		// get actual timestep
		// return delta eps (actual step - last step)
		if (load) {
			macroscopic_strain = strains_macro[timestep];
			macroscopic_strain -= strains_macro[timestep - 1];
		}
	} else {
		cerr << "Unknown loadtype (get macroscopic strain)" << endl;
		exit(0);
	}
//	deallog << "MacroStrain: " << macroscopic_strain << endl;
	return macroscopic_strain;
}
template<int dim>
SymmetricTensor<2, dim> BoundaryValues<dim>::get_full_macroscopic_strain() const {
	SymmetricTensor<2, dim> macroscopic_strain;

	if (mode == linear_displacement) {
		double displacementvalue = para->get_double("displacement")
				* time_loadsteps[time_loadsteps.size() - 1] / duration();
		if (para->get("straindirection").compare("11") == 0) {
			macroscopic_strain[0][0] = displacementvalue;
		} else if (para->get("straindirection").compare("22") == 0) {
			macroscopic_strain[1][1] = displacementvalue;
		} else if (para->get("straindirection").compare("33") == 0) {
			macroscopic_strain[2][2] = displacementvalue;
		} else if (para->get("straindirection").compare("12") == 0) {
			macroscopic_strain[0][1] = displacementvalue / 2;
		} else if (para->get("straindirection").compare("13") == 0) {
			macroscopic_strain[0][2] = displacementvalue / 2;
		} else if (para->get("straindirection").compare("23") == 0) {
			macroscopic_strain[1][2] = displacementvalue / 2;
		} else if (para->get("straindirection").compare("detailed") == 0) {
			macroscopic_strain = macro_strain * displacementvalue;
		}
	} else if (mode == prescribed_from_file) {
		// return eps(aktuell) - eps(0)
		macroscopic_strain = strains_macro[timestep];
		macroscopic_strain -= strains_macro[0];
	} else {
		cerr << "Unknown loadtype (get full macroscopic strain)" << endl;
	}

	return macroscopic_strain;
}
template<int dim>
void BoundaryValues<dim>::next_timestep() {
//	cout << "Next timestep" << endl;
	if (mode == linear_displacement) {
		propotion_applied_via_newtonguess = 0;
		double dt;
		if (load) {
			// dt = para.get_double("duration") / para.get_double("n_loadsteps");
			dt = duration() / para->get_double("n_loadsteps");
			deallog << "delta_t: " << dt << endl;
		} else {
			dt = 0;
		}
		time_loadsteps.push_back(
				time_loadsteps[time_loadsteps.size() - 1] + dt);

		if ((time_loadsteps[time_loadsteps.size() - 1] - duration()) < 1e-12) {
		} else {
			time_loadsteps[time_loadsteps.size() - 1] = duration();
		}
	} else if (mode == prescribed_from_file) {
		// counter ++
		timestep++;
		load = true;
		SymmetricTensor<2, dim> strain = get_macroscopic_strain();
		double dt = get_dt();

		deallog << "Equivalent strainrate: " << strain.norm() / dt << endl;

//		cout << "Macro strain: " << get_macroscopic_strain() << endl;
//		cout << "Full macro strain: " << get_full_macroscopic_strain() << endl;
	} else {
		cerr << "Unknown loadtype (next timestep)" << endl;
	}

}
template<int dim>
void BoundaryValues<dim>::first_timestep() {
//	cout << "First timestep" << endl;
	if (mode == linear_displacement) {
		load = true;
		time_loadsteps.push_back(0);
	} else if (mode == prescribed_from_file) {
		// counter = 0
		timestep = 0;
		load = true;
	} else {
		cerr << "Unknown loadtype (first timestep)" << endl;
	}
}

template<int dim>
void BoundaryValues<dim>::deactivate_inhomogenious_bc() {
//	cout << "next newtonstep" << endl;
	load = false;
}

template<int dim>
bool BoundaryValues<dim>::apply_load() const {
//	cout << "apply load" << endl;
	return load;
}

template<int dim>
double BoundaryValues<dim>::get_dt() const {
//	cout << "get_dt" << endl;
	if (mode == linear_displacement) {
		return time_loadsteps[time_loadsteps.size() - 1]
				- time_loadsteps[time_loadsteps.size() - 2];
	} else if (mode == prescribed_from_file) {
		// delta t
		return times_macro[timestep] - times_macro[timestep - 1];
	} else {
		cerr << "Unknown loadtype (get dt in boundary)" << endl;
		exit(0);
	}
}

template<int dim>
double BoundaryValues<dim>::apply_bc_via_newtonguess(double factor) {
//	cout << "apply via newtonguess" << endl;
	if (mode == linear_displacement) {
		double t_actual_timestep = time_loadsteps[time_loadsteps.size() - 1];
		double t_previous_timestep = time_loadsteps[time_loadsteps.size() - 2];
		double delta_t_actual_timestep = time_loadsteps[time_loadsteps.size()
				- 1] - time_loadsteps[time_loadsteps.size() - 2];
		double delta_t_previous_timestep = time_loadsteps[time_loadsteps.size()
				- 2] - time_loadsteps[time_loadsteps.size() - 3];

		// alte Zeit plus alte Zeitschrittweite
		double t_with_newtonguess = t_actual_timestep
				+ delta_t_previous_timestep * factor;
		if (t_actual_timestep - duration() < 1e-12) {
			if (t_with_newtonguess - duration() < 1e-12) {
				time_loadsteps[time_loadsteps.size() - 1] = t_with_newtonguess;
				deallog << "Added newtonguess, time: " << t_actual_timestep
						<< " -> " << t_with_newtonguess << endl;
				propotion_applied_via_newtonguess = delta_t_previous_timestep
						/ (delta_t_previous_timestep + delta_t_actual_timestep);
				return factor;
			} else {
				time_loadsteps[time_loadsteps.size() - 1] = duration();
				deallog << "Factor is reduced" << endl;
				propotion_applied_via_newtonguess = (duration()
						- t_actual_timestep)
						/ (duration() - t_previous_timestep);
				return (duration() - t_actual_timestep)
						/ delta_t_previous_timestep;
			}
		} else {
			// should never be reached
			deallog << "allready finished, this point should never be reached" << endl;
			return 0;
		}
	} else if (mode == prescribed_from_file) {
		return 0;
	} else {
		cerr << "Unknown loadtype (apply bc via newtonguess)" << endl;
		exit(0);
	}
}

template<int dim>
double BoundaryValues<dim>::get_time() {
	return time_loadsteps[time_loadsteps.size() - 1];
}

template<int dim>
double BoundaryValues<dim>::get_dx() const {
//	cout << "get_dx" << endl;
	if (mode == linear_displacement) {

		if (load) {
			double t_1 = time_loadsteps[time_loadsteps.size() - 2];
			double t_2 = time_loadsteps[time_loadsteps.size() - 1];
			t_1 = t_1 + (t_2 - t_1) * propotion_applied_via_newtonguess;

			return get_displacement(t_2) - get_displacement(t_1);
		} else {
			return 0;
		}
	} else if (mode == prescribed_from_file) {
		cerr << "get_dx in boundary should not be called" << endl;
		exit(0);
	} else {
		cerr << "Unknown loadtype (get dx in boundary values)" << endl;
		exit(0);
	}
}

template<int dim>
double BoundaryValues<dim>::get_displacement(double t) const {
//	cout << "get displacement" << endl;
	if (mode == linear_displacement) {
// linear loading
		return para->get_double("displacement") * t / duration();
	} else if (mode == prescribed_from_file) {
		cerr << "should not be called (get displacement)" << endl;
		exit(0);
	} else {
		cerr << "Unknown loadtype (get displacement)" << endl;
		exit(0);
	}
}

template<int dim>
double BoundaryValues<dim>::get_theta(double t, double z) const {
//	cout << "get theta tz" << endl;
	if (mode == linear_displacement) {
		double theta_top =
				para->get_double("t_start")
						+ t / duration()
								* (para->get_double("t_end")
										- para->get_double("t_start"));
		double theta_bottom = para->get_double("t_end");
		return theta_bottom + z * (theta_top - para->get_double("t_end"));
	} else if (mode == prescribed_from_file) {
		// return actual temperature
		if (t == 0) {
			return temperatures_macro[0];
		} else {
			return temperatures_macro[timestep];
		}
	} else {
		cerr << "Unknown loadtype (get theta(t,z))" << endl;
		return -1;
	}
}

template<int dim>
double BoundaryValues<dim>::get_theta(const Point<dim> &p) const {
//	cout << "get theta p" << endl;
	if (mode == linear_displacement) {
		return get_theta(time_loadsteps[time_loadsteps.size() - 1], p[2]);
	} else if (mode == prescribed_from_file) {
		return get_theta(times_macro[timestep], p[2]);
	} else {
		cerr << "Unknown loadtype (get theta(p))" << endl;
		return -1;
	}
}

template<int dim>
double BoundaryValues<dim>::get_initial_theta(const Point<dim> &p) const {
//	cout << "get initial theta" << endl;
	return get_theta(0, p[2]);
}

template<int dim>
bool BoundaryValues<dim>::hasNextIncrement() {
//	cout << "has next increment" << endl;
	if (mode == linear_displacement) {
		return ((time_loadsteps[time_loadsteps.size() - 1] - duration())
				< -1e-12);
	} else if (mode == prescribed_from_file) {
		// has next step?
		return timestep + 1 < temperatures_macro.size();
	} else {
		cerr << "Unknown loadtype (get theta (p))" << endl;
		return false;
	}
}

template<int dim>
int BoundaryValues<dim>::NumberOfLoadsteps() {
//	cout << "get theta p" << endl;
	if (mode == linear_displacement) {
		return para->get_double("n_loadsteps");
	} else if (mode == prescribed_from_file) {
		return times_macro.size() - 1;
	} else {
		cerr << "Unknown loadtype (get theta(p))" << endl;
		return -1;
	}
}

#endif /* CODE_VARIOUS_OLD_BOUNDARY_VALUES_H_ */
