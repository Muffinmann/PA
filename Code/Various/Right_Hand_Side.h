/*
 * Right_Hand_Side.h
 *
 *  Created on: Jun 30, 2016
 *      Author: iwtm84
 */

#ifndef CODE_VARIOUS_RIGHT_HAND_SIDE_H_
#define CODE_VARIOUS_RIGHT_HAND_SIDE_H_

using namespace dealii;



template<int dim>
class RightHandSide: public Function<dim> {
public:
	RightHandSide();
	virtual ~RightHandSide() {
	}
	;
	virtual void vector_value(const Point<dim> &p,
			Vector<double> &values) const;

	virtual void vector_value_list(const std::vector<Point<dim> > &points,
			std::vector<Vector<double> > &value_list) const;
};
template<int dim>
RightHandSide<dim>::RightHandSide() :
		Function<dim>(dim) {
}
template<int dim>

inline
void RightHandSide<dim>::vector_value(
		const Point<dim> &p __attribute__((unused)),
		Vector<double> &values) const {
	Assert(values.size() == dim, ExcDimensionMismatch (values.size(), dim));
	Assert(dim >= 2, ExcNotImplemented());

	for (unsigned int i = 0; i < dim; i++) {
		values(i) = 0;
	}
}

template<int dim>
void RightHandSide<dim>::vector_value_list(
		const std::vector<Point<dim> > &points,
		std::vector<Vector<double> > &value_list) const {
	Assert(value_list.size() == points.size(),
			ExcDimensionMismatch (value_list.size(), points.size()));

	const unsigned int n_points = points.size();

	for (unsigned int p = 0; p < n_points; ++p)
		RightHandSide<dim>::vector_value(points[p], value_list[p]);
}




#endif /* CODE_VARIOUS_RIGHT_HAND_SIDE_H_ */
