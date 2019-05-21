/*
 * Postprocessing.h
 *
 *  Created on: Apr 5, 2016
 *      Author: iwtm84
 */

#ifndef CODE_POSTPROCESSING_POSTPROCESSING_H_
#define CODE_POSTPROCESSING_POSTPROCESSING_H_

#include <string.h>
#include "Extractor.h"
#include "Scratches.h"
#include "DataOutCrystalwise.h"
using namespace std;
using namespace dealii;

template<int dim>
struct ScratchDataPost {
	hp::FEValues<dim> hp_fe_values;

	Vector<double> &current_solution;
	Vector<double> &prior_solution; // solution of the previous loadstep
	Vector<double> &current_solution_abs; // absolute solutions
	double dt;

	ScratchDataPost(const hp::FECollection<dim> &fe_collection,
			const hp::QCollection<dim> &quadrature_collection,
			const UpdateFlags update_flags, Vector<double> &current_solution,
			Vector<double> &prior_solution, Vector<double> &current_solution_abs, double dt, int /*check*/) :

			hp_fe_values(fe_collection, quadrature_collection, update_flags), current_solution(
					current_solution), prior_solution(prior_solution), current_solution_abs(
							current_solution_abs), dt(dt) {
	}
	ScratchDataPost(const ScratchDataPost &scratch) :
			hp_fe_values(scratch.hp_fe_values.get_fe_collection(),
					scratch.hp_fe_values.get_quadrature_collection(),
					scratch.hp_fe_values.get_update_flags()), current_solution(
					scratch.current_solution), prior_solution(
					scratch.prior_solution), current_solution_abs(
					scratch.current_solution_abs), dt(scratch.dt) {
	}
};

template<int dim>
struct PerTaskDataPost {
	std::vector<types::global_dof_index> local_dof_indices;
	int crystalID;
	vector<Vector<double>> postSolutions;

	PerTaskDataPost(const hp::FECollection<dim> &fe_collection, int ntasks) :
			local_dof_indices(fe_collection.max_dofs_per_cell()), crystalID(-1), postSolutions(
					ntasks, Vector<double>(fe_collection.max_dofs_per_cell())) {
	}
};

template<int dim, int nslip>
class Postprocessing {
public:
	Postprocessing(const hp::DoFHandler<dim> &DH,
			const hp::QCollection<dim> &QF, const hp::FECollection<dim> &fe,
			std::vector<PerCrystalData<dim, nslip> > &crystal_data,
			FEValuesExtractors::Vector &fe_values_extractor_displacement,
			vector<vector<FEValuesExtractors::Scalar>> &fe_values_extractor_slips);
	~Postprocessing();
	void writeResults(unsigned int loadstep, double dt, Vector<double> &current_solution,
			Vector<double> &prior_solution, Vector<double> &current_solution_abs);
	void init();
	void save_double_slip_diagrams(
			vector<vector<FEValuesExtractors::Scalar>> extractor_slips,
			const Vector<double> &solution);

private:
	void addNodewiseByDomainData(
			DataOutCrystalwise<dim, nslip, hp::DoFHandler<dim> > &data_out,
			std::map<string, Extract_scalar_value<dim, nslip> > &extractors,
			Extract_displacement_value<dim, nslip> &extractor_displacement,
			vector<Extract_scalar_value<dim, nslip>> &extractor_slips,
			/*vector<ScalarScratch<dim, nslip>> &tasks,*/
			double dt,
			Vector<double> &current_solution, Vector<double> &prior_solution,
			Vector<double> &current_solution_abs);

	void calculateMassMatrix();

	void postprocessing_on_one_cell(
			const typename hp::DoFHandler<dim>::active_cell_iterator &cell,
			ScratchDataPost<dim> &scratch, PerTaskDataPost<dim> &data);
	void copy_local_to_global(const PerTaskDataPost<dim> &data);

	const hp::DoFHandler<dim> &hp_dof_handler;
	const hp::QCollection<dim> &quadrature_collection;
	const hp::FECollection<dim> &fe_collection;
	std::vector<PerCrystalData<dim, nslip> > &crystal_data;
	FEValuesExtractors::Vector &fe_values_extractor_displacement;
	vector<vector<FEValuesExtractors::Scalar>> &fe_values_extractor_slips;

	vector<int> domainIDs;

	vector<Vector<double>> mass_matrix_postprocessing;
	vector<ScalarScratch<dim, nslip>> *tasks;
};

template<int dim, int nslip> Postprocessing<dim, nslip>::Postprocessing(
		const hp::DoFHandler<dim> &DH, const hp::QCollection<dim> &QF,
		const hp::FECollection<dim> &fe,
		std::vector<PerCrystalData<dim, nslip> > &crystal_data,
		FEValuesExtractors::Vector &fe_values_extractor_displacement,
		vector<vector<FEValuesExtractors::Scalar>> &fe_values_extractor_slips) :
		hp_dof_handler(DH), quadrature_collection(QF), fe_collection(fe), crystal_data(
				crystal_data), fe_values_extractor_displacement(
				fe_values_extractor_displacement), fe_values_extractor_slips(
				fe_values_extractor_slips), tasks(nullptr) {
}

template<int dim, int nslip> Postprocessing<dim, nslip>::~Postprocessing() {
}

template<int dim, int nslip> void Postprocessing<dim, nslip>::init() {
	calculateMassMatrix();
}

#include "Postprocessing.cc"

#endif /* CODE_POSTPROCESSING_POSTPROCESSING_H_ */
