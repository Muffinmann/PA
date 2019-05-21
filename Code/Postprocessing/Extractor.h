/*
 * Extractor.h
 *
 *  Created on: Apr 5, 2016
 *      Author: iwtm84
 */

#ifndef CODE_POSTPROCESSING_EXTRACTOR_H_
#define CODE_POSTPROCESSING_EXTRACTOR_H_
#include <deal.II/lac/vector.h>
#include <deal.II/numerics/data_postprocessor.h>

using namespace std;
using namespace dealii;

template<int dim, int nslip>
class Extract_scalar_value: public DataPostprocessorScalar<dim> {
public:
	Extract_scalar_value();
	Extract_scalar_value(std::string name);
	Extract_scalar_value(std::string name, int component);
	virtual
	~Extract_scalar_value();
	void
	compute_derived_quantities_vector(const std::vector<Vector<double> > &uh,
			const std::vector<std::vector<Tensor<1, dim> > > &duh,
			const std::vector<std::vector<Tensor<2, dim> > > &dduh,
			const std::vector<Point<dim> > &normals,
			const std::vector<Point<dim> > &evaluation_points,
			std::vector<Vector<double> > &computed_quantities) const;
	enum postprocessingmodes {
		unknown = -1, scalar_data = 0, slips = 1
	};
private:
	int component;
	postprocessingmodes postprocessingmode;
};
template<int dim, int nslip>
Extract_scalar_value<dim, nslip>::Extract_scalar_value() :
		DataPostprocessorScalar<dim>("Standard", update_values), component(0), postprocessingmode(
				postprocessingmodes::scalar_data) {
}
template<int dim, int nslip>
Extract_scalar_value<dim, nslip>::Extract_scalar_value(std::string name) :
		DataPostprocessorScalar<dim>(name, update_values), component(0), postprocessingmode(
				postprocessingmodes::scalar_data) {
}
template<int dim, int nslip>
Extract_scalar_value<dim, nslip>::Extract_scalar_value(std::string name,
		int component) :
		DataPostprocessorScalar<dim>(name, update_values), component(component), postprocessingmode(
				postprocessingmodes::slips) {
}

template<int dim, int nslip>
Extract_scalar_value<dim, nslip>::~Extract_scalar_value() {

}
template<int dim, int nslip>
void Extract_scalar_value<dim, nslip>::compute_derived_quantities_vector(
		const std::vector<Vector<double> > &uh,
		const std::vector<std::vector<Tensor<1, dim> > > & /*duh*/,
		const std::vector<std::vector<Tensor<2, dim> > > & /*dduh*/,
		const std::vector<Point<dim> > & /*normals*/,
		const std::vector<Point<dim> > & /*evaluation_points*/,
		std::vector<Vector<double> > &computed_quantities) const {
	Assert(computed_quantities.size() == uh.size(),
			ExcDimensionMismatch(computed_quantities.size(), uh.size()));
	for (unsigned int i = 0; i < computed_quantities.size(); i++) {
		Assert(computed_quantities[i].size() == 1,
				ExcDimensionMismatch(computed_quantities[i].size(), 1));
		if (postprocessingmode == postprocessingmodes::scalar_data) {
			computed_quantities[i](0) = uh[i](component);
		} else if (postprocessingmode == postprocessingmodes::slips) {
			computed_quantities[i](0) = 0;
			for (unsigned int j = component; j < uh[i].size(); j += nslip) {
				computed_quantities[i](0) += uh[i](j);
			}
		}
	}
}

template<int dim, int nslip>
class Extract_displacement_value: public DataPostprocessorVector<dim> {
public:
	Extract_displacement_value();
	virtual
	~Extract_displacement_value();
	void
	compute_derived_quantities_vector(const std::vector<Vector<double> > &uh,
			const std::vector<std::vector<Tensor<1, dim> > > &duh,
			const std::vector<std::vector<Tensor<2, dim> > > &dduh,
			const std::vector<Point<dim> > &normals,
			const std::vector<Point<dim> > &evaluation_points,
			std::vector<Vector<double> > &computed_quantities) const;
};
template<int dim, int nslip>
Extract_displacement_value<dim, nslip>::Extract_displacement_value() :
		DataPostprocessorVector<dim>("Displacement", update_values) {
}

template<int dim, int nslip>
Extract_displacement_value<dim, nslip>::~Extract_displacement_value() {

}
template<int dim, int nslip>
void Extract_displacement_value<dim, nslip>::compute_derived_quantities_vector(
		const std::vector<Vector<double> > &uh,
		const std::vector<std::vector<Tensor<1, dim> > > & /*duh*/,
		const std::vector<std::vector<Tensor<2, dim> > > & /*dduh*/,
		const std::vector<Point<dim> > & /*normals*/,
		const std::vector<Point<dim> > & /*evaluation_points*/,
		std::vector<Vector<double> > &computed_quantities) const {
	Assert(computed_quantities.size() == uh.size(),
			ExcDimensionMismatch(computed_quantities.size(), uh.size()));
	for (unsigned int i = 0; i < computed_quantities.size(); i++) {
		Assert(computed_quantities[i].size() == 3,
				ExcDimensionMismatch(computed_quantities[i].size(), 1));
//		Assert(uh[i].size() == dim+nslip, ExcDimensionMismatch(uh[i].size(), 2));
		for (int j = 0; j < dim; j++) {
			computed_quantities[i](j) = uh[i](j);
		}
	}
}

#endif /* CODE_POSTPROCESSING_EXTRACTOR_H_ */
