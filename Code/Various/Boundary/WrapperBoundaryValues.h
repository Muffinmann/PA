/*
 * WrapperBoundaryValues.h
 *
 *  Created on: Mar 14, 2018
 *      Author: iwtm84
 */

#ifndef CODE_VARIOUS_BOUNDARY_WRAPPERBOUNDARYVALUES_H_
#define CODE_VARIOUS_BOUNDARY_WRAPPERBOUNDARYVALUES_H_

using namespace dealii;

template<int dim>
class WrapperBoundaryValues: public Function<dim> {
public:
	WrapperBoundaryValues();
	WrapperBoundaryValues(BoundaryValues<dim> *boundary,
			unsigned int n_components) :
			Function<dim>(n_components), boundary(boundary), solution_dim(
					n_components) {
	}
	;
	virtual ~WrapperBoundaryValues() {

	}

	virtual void vector_value_list(const std::vector<Point<dim> > &points,
			std::vector<Vector<double> > &value_list) const;

	// values at point p
	virtual void vector_value(const Point<dim> &p,
			Vector<double> &values) const;

	BoundaryValues<dim> *boundary;
	unsigned int solution_dim;
};

/*
 * uses vector_value to calculate values at all points in vector points
 */
template<int dim>
void WrapperBoundaryValues<dim>::vector_value_list(
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
void WrapperBoundaryValues<dim>::vector_value(const Point<dim> &p,
		Vector<double> &values) const {
	Assert(values.size() == solution_dim,
			ExcDimensionMismatch(values.size(), solution_dim));

	Vector<double> values_dim(dim);
	boundary->vector_value(p, values_dim);

	values = 0;
	for (unsigned int i = 0; i < dim; ++i) {
		values(i) = values_dim(i);
	}
}

#endif /* CODE_VARIOUS_BOUNDARY_WRAPPERBOUNDARYVALUES_H_ */
