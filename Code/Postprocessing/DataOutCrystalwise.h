/*
 * DataOutCrystalwise.h
 *
 *  Created on: Apr 7, 2016
 *      Author: iwtm84
 */

#ifndef CODE_POSTPROCESSING_DATAOUTCRYSTALWISE_H_
#define CODE_POSTPROCESSING_DATAOUTCRYSTALWISE_H_

using namespace dealii;

template<int dim, int nslip, typename DoFHandlerType=DoFHandler<dim> >
class DataOutCrystalwise: public DataOut<dim, DoFHandlerType> {
public:
	DataOutCrystalwise();
	virtual ~DataOutCrystalwise();
	typename DataOut<dim>::cell_iterator next_cell(
			const typename DataOut<dim>::cell_iterator & cell);
	typename DataOut<dim>::cell_iterator first_cell();
	void  setID(int i){
		ID=i;
	}
private:
	int ID;

};

template<int dim, int nslip, typename DoFHandlerType>
DataOutCrystalwise<dim, nslip, DoFHandlerType>::DataOutCrystalwise(){
	ID=0;
}

template<int dim, int nslip, typename DoFHandlerType>
DataOutCrystalwise<dim, nslip, DoFHandlerType>::~DataOutCrystalwise() {

}

template<int dim, int nslip, typename DoFHandlerType>
typename DataOut<dim>::cell_iterator DataOutCrystalwise<dim, nslip, DoFHandlerType>::first_cell() {
	typename Triangulation<dim>::active_cell_iterator active_cell =
			this->triangulation->begin_active(), last_cell =
			this->triangulation->last_active();

	do {
		if(active_cell->material_id() == ID){
			break;
		}
		if (active_cell == last_cell) {
			active_cell++;
			break;
		}
		active_cell++;
	} while (true);
	return active_cell;
}

template<int dim, int nslip, typename DoFHandlerType>
typename DataOut<dim>::cell_iterator DataOutCrystalwise<dim, nslip, DoFHandlerType>::next_cell(
		const typename DataOut<dim>::cell_iterator & cell) {
	// convert the iterator to an active_iterator and advance this to the next
	// active cell
	typename Triangulation<dim>::active_cell_iterator active_cell = cell,
			last_cell = this->triangulation->last_active();

	do {
		if (active_cell == last_cell) {
			active_cell++;
			break;
		}
		++active_cell;
	} while (active_cell->material_id() != ID);

	return active_cell;
}

#endif /* CODE_POSTPROCESSING_DATAOUTCRYSTALWISE_H_ */
