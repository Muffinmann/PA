/*
 * writeResults.cc
 *
 *  Created on: Apr 6, 2016
 *      Author: iwtm84
 */

#include "Postprocessing.h"
#include "Scratches.h"
#include "PostprocessingFunktions.h"
#include "Extractor.h"
#include "DataOutCrystalwise.h"
#include <string.h>

using namespace std;
using namespace dealii;

template<int dim, int nslip>
void Postprocessing<dim, nslip>::writeResults(unsigned int loadstep, double dt,
		Vector<double> &current_solution, Vector<double> &prior_solution,
		Vector<double> &current_solution_abs) {

	std::map<string, Extract_scalar_value<dim, nslip>> extractors_bydomain;
	Extract_displacement_value<dim, nslip> extractor_displacement;
	vector<Extract_scalar_value<dim, nslip>> extractor_slips;
	for (int i = 0; i < nslip; i++) {
		std::stringstream sstm;
		if (i < 9) {
			sstm << "Slipsystem_0" << i + 1;
		} else {
			sstm << "Slipsystem_" << i + 1;
		}
		extractor_slips.push_back(
				Extract_scalar_value<dim, nslip>(sstm.str(), i + dim));
	}
	vector<ScalarScratch<dim, nslip>> tasks_nodewise_bydomain;
	tasks = &tasks_nodewise_bydomain;

	DataOutCrystalwise<dim, nslip, hp::DoFHandler<dim> > data_out_bydomain;
	data_out_bydomain.attach_dof_handler(hp_dof_handler);

	tasks_nodewise_bydomain.insert(tasks_nodewise_bydomain.begin(),
			ScalarScratch<dim, nslip>("Sig11", hp_dof_handler.n_dofs(),
					domainIDs.size(), calculateSig11));

	tasks_nodewise_bydomain.insert(tasks_nodewise_bydomain.begin(),
			ScalarScratch<dim, nslip>("Sig22", hp_dof_handler.n_dofs(),
					domainIDs.size(), calculateSig22));

	tasks_nodewise_bydomain.insert(tasks_nodewise_bydomain.begin(),
			ScalarScratch<dim, nslip>("Sig33", hp_dof_handler.n_dofs(),
					domainIDs.size(), calculateSig33));

	tasks_nodewise_bydomain.insert(tasks_nodewise_bydomain.begin(),
			ScalarScratch<dim, nslip>("Sig12", hp_dof_handler.n_dofs(),
					domainIDs.size(), calculateSig12));

	tasks_nodewise_bydomain.insert(tasks_nodewise_bydomain.begin(),
			ScalarScratch<dim, nslip>("Sig13", hp_dof_handler.n_dofs(),
					domainIDs.size(), calculateSig13));

	tasks_nodewise_bydomain.insert(tasks_nodewise_bydomain.begin(),
			ScalarScratch<dim, nslip>("Sig23", hp_dof_handler.n_dofs(),
					domainIDs.size(), calculateSig23));

	tasks_nodewise_bydomain.insert(tasks_nodewise_bydomain.begin(),
			ScalarScratch<dim, nslip>("VonMisesStress", hp_dof_handler.n_dofs(),
					domainIDs.size(), calculateMises));

	tasks_nodewise_bydomain.insert(tasks_nodewise_bydomain.begin(),
			ScalarScratch<dim, nslip>("EquivalentPlasticStrain",
					hp_dof_handler.n_dofs(), domainIDs.size(),
					calculateEquivalentPlasticStrain));

	tasks_nodewise_bydomain.insert(tasks_nodewise_bydomain.begin(),
			ScalarScratch<dim, nslip>("EquivalentPlasticStrainRate",
					hp_dof_handler.n_dofs(), domainIDs.size(),
					calculateEquivalentPlasticStrainRate));

	tasks_nodewise_bydomain.insert(tasks_nodewise_bydomain.begin(),
			ScalarScratch<dim, nslip>("EquivalentThermalStrain",
					hp_dof_handler.n_dofs(), domainIDs.size(),
					calculateEquivalentThermalStrain));

	tasks_nodewise_bydomain.insert(tasks_nodewise_bydomain.begin(),
			ScalarScratch<dim, nslip>("Tiltangles", hp_dof_handler.n_dofs(),
					domainIDs.size(), calculateTiltangle));

	tasks_nodewise_bydomain.insert(tasks_nodewise_bydomain.begin(),
			ScalarScratch<dim, nslip>("Rotationangle", hp_dof_handler.n_dofs(),
					domainIDs.size(), calculateRotationangle));

	tasks_nodewise_bydomain.insert(tasks_nodewise_bydomain.begin(),
			ScalarScratch<dim, nslip>("Temperatur", hp_dof_handler.n_dofs(),
					domainIDs.size(), calculateTemperatur));

	tasks_nodewise_bydomain.insert(tasks_nodewise_bydomain.begin(),
			ScalarScratch<dim, nslip>("A", hp_dof_handler.n_dofs(),
					domainIDs.size(), calculateA));

	tasks_nodewise_bydomain.insert(tasks_nodewise_bydomain.begin(),
			ScalarScratch<dim, nslip>("ScrewDislocationDensity",
					hp_dof_handler.n_dofs(), domainIDs.size(),
					calculateScrewDislocationDensity));

	tasks_nodewise_bydomain.insert(tasks_nodewise_bydomain.begin(),
			ScalarScratch<dim, nslip>("EdgeDislocationDensity",
					hp_dof_handler.n_dofs(), domainIDs.size(),
					calculateEdgeDislocationDensity));

	tasks_nodewise_bydomain.insert(tasks_nodewise_bydomain.begin(),
			ScalarScratch<dim, nslip>("DislocationDensity",
					hp_dof_handler.n_dofs(), domainIDs.size(),
					calculateDislocationDensity));

	tasks_nodewise_bydomain.insert(tasks_nodewise_bydomain.begin(),
			ScalarScratch<dim, nslip>("DomainID", hp_dof_handler.n_dofs(),
					domainIDs.size(), calculateDomainID));

//	for (int i = 0; i < nslip; i++) {
//		tasks_nodewise_bydomain.insert(tasks_nodewise_bydomain.begin(),
//				ScalarScratch<dim, nslip>("Xi_" + to_string(i),
//						hp_dof_handler.n_dofs(), domainIDs.size(),
//						[i](Material<dim, nslip> &history) {
//							return history.get_Xi(i).norm();
//						}));
//	}
//
//	for (int i = 0; i < dim; i++) {
//		tasks_nodewise_bydomain.insert(tasks_nodewise_bydomain.begin(),
//				ScalarScratch<dim, nslip>("100_" + to_string(i),
//						hp_dof_handler.n_dofs(), domainIDs.size(),
//						[i](Material<dim, nslip> &history) {
//							return history.get_direction(2)[i];
//						}));
//	}

	addNodewiseByDomainData(data_out_bydomain, extractors_bydomain,
			extractor_displacement, extractor_slips, /*tasks_nodewise_bydomain,*/
			dt, current_solution, prior_solution, current_solution_abs);
	// Offstream for vtu output
	//-------------------------------
	string filename;
	std::stringstream sstm;
	{
		sstm << "Results/solution-bydomain" << loadstep << ".vtu";
		filename = sstm.str();
	}
	std::ofstream output_bydomain(filename.c_str());
	data_out_bydomain.write_vtu(output_bydomain);
	data_out_bydomain.clear();
	//-------------------------------------------
	tasks = NULL;
}

template<int dim, int nslip>
void Postprocessing<dim, nslip>::postprocessing_on_one_cell(
		const typename hp::DoFHandler<dim>::active_cell_iterator &cell,
		ScratchDataPost<dim> &scratch, PerTaskDataPost<dim> &data) {

//	Material<dim, nslip> *local_quadrature_points_data =
//			reinterpret_cast<Material<dim, nslip>*>(cell->user_pointer());
	gpData<dim, nslip> gpData;
	gpData.perCrystalData = &crystal_data[cell->active_fe_index()];

	data.crystalID = cell->active_fe_index();
	scratch.hp_fe_values.reinit(cell);
	const FEValues<dim> &fe_values =
			scratch.hp_fe_values.get_present_fe_values();
	unsigned int n_q_points = fe_values.n_quadrature_points;

	const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;
	data.local_dof_indices.resize(dofs_per_cell);
	cell->get_dof_indices(data.local_dof_indices);

	for (unsigned int t = 0; t < data.postSolutions.size(); t++) {
		data.postSolutions[t] = 0;
	}

	// calculate quantities
	vector<SymmetricTensor<2, dim> > strains(n_q_points);
	fe_values[fe_values_extractor_displacement].get_function_symmetric_gradients(
			scratch.current_solution, strains);

	vector<vector<double> > slips(nslip, vector<double>(n_q_points));
	vector<vector<double> > prior_slips(nslip, vector<double>(n_q_points));
	vector<vector<double> > current_gamma_abs_cell(nslip, vector<double>(n_q_points));
	vector<vector<Tensor<1, dim> > > slip_grads(nslip,
			vector<Tensor<1, dim> >(n_q_points));
	for (unsigned int i = 0;
			i < fe_values_extractor_slips[cell->active_fe_index()].size();
			i++) {
		fe_values[fe_values_extractor_slips[cell->active_fe_index()][i]].get_function_values(
				scratch.current_solution, slips[i]);
		fe_values[fe_values_extractor_slips[cell->active_fe_index()][i]].get_function_values(
				scratch.prior_solution, prior_slips[i]);
		fe_values[fe_values_extractor_slips[cell->active_fe_index()][i]].get_function_values(
				scratch.current_solution_abs, current_gamma_abs_cell[i]);
		fe_values[fe_values_extractor_slips[cell->active_fe_index()][i]].get_function_gradients(
				scratch.current_solution, slip_grads[i]);
	}

	Vector<double> dGamma(slips.size());
	Vector<double> gamma(slips.size());
	Vector<double> current_gamma_abs_gp(slips.size());
	vector<Tensor<1, dim> > gradGamma(slip_grads.size());

	for (unsigned int q_point = 0; q_point < n_q_points; ++q_point) {
		double JxW = fe_values.JxW(q_point);

		// calculate values for gpData
		for (unsigned int i = 0; i < dGamma.size(); i++) {
			gamma[i] = slips[i][q_point];
			current_gamma_abs_gp[i] = current_gamma_abs_cell[i][q_point];
			dGamma[i] = slips[i][q_point] - prior_slips[i][q_point];
			gradGamma[i] = slip_grads[i][q_point];
		}

		//set values in gpData
		gpData.dGamma_actualStep = dGamma;
		gpData.dt = scratch.dt;
		gpData.gamma = gamma;
		gpData.gradGamma = gradGamma;
		gpData.current_gamma_abs = current_gamma_abs_gp;
		gpData.strain = strains[q_point];
		gpData.temperature = 23;

		for (unsigned int i = 0; i < dofs_per_cell; ++i) {
			const unsigned int component_i =
					fe_collection[cell->active_fe_index()].system_to_component_index(
							i).first;
			if (component_i >= 1) { // only have to calculate the postprocessing for the first dof TODO try that
				continue;
			}
			double N_i = fe_values.shape_value(i, q_point);
			for (unsigned int t = 0; t < tasks->size(); t++) {
				data.postSolutions[t](i) += N_i
						* tasks->at(t).post_function(gpData) * JxW;
			}
		}
	}
}
template<int dim, int nslip>
void Postprocessing<dim, nslip>::copy_local_to_global(
		const PerTaskDataPost<dim> &data) {

	for (unsigned int t = 0; t < tasks->size(); t++) {
		for (unsigned int i = 0; i < data.local_dof_indices.size(); ++i) {
			tasks->at(t).solution[data.crystalID](data.local_dof_indices[i]) +=
					data.postSolutions[t](i);
		}
	}
}

template<int dim, int nslip>
void Postprocessing<dim, nslip>::addNodewiseByDomainData(
		DataOutCrystalwise<dim, nslip, hp::DoFHandler<dim> > &data_out,
		std::map<string, Extract_scalar_value<dim, nslip> > &extractors,
		Extract_displacement_value<dim, nslip> &extractor_displacement,
		vector<Extract_scalar_value<dim, nslip>> &extractor_slips,
		/*vector<ScalarScratch<dim, nslip>> &tasks,*/double dt,
		Vector<double> &current_solution, Vector<double> &prior_solution,
		Vector<double> &current_solution_abs) {

	typename hp::DoFHandler<dim>::active_cell_iterator cell =
			hp_dof_handler.begin_active(), endc = hp_dof_handler.end();

	PerTaskDataPost<dim> per_task_data(fe_collection, tasks->size());
	ScratchDataPost<dim> scratch(fe_collection, quadrature_collection,
			update_values | update_gradients | update_JxW_values,
			current_solution, prior_solution, current_solution_abs, dt, 0);
	WorkStream::run(cell, endc, *this,
			&Postprocessing<dim, nslip>::postprocessing_on_one_cell,
			&Postprocessing<dim, nslip>::copy_local_to_global, scratch,
			per_task_data);

	DataOutCrystalwise<dim, nslip, hp::DoFHandler<dim> > data_out_tmp;
	data_out_tmp.attach_dof_handler(hp_dof_handler);
	const int subdivisions = 0;

	for (vector<int>::const_iterator ID = domainIDs.begin();
			ID != domainIDs.end(); ++ID) {

		// Add solutionvector
		if (ID == domainIDs.begin()) {
			data_out.add_data_vector(current_solution, extractor_displacement);
			for (unsigned int i = 0; i < extractor_slips.size(); i++) {
				data_out.add_data_vector(current_solution, extractor_slips[i]);
			}
		} else {
			data_out_tmp.add_data_vector(current_solution,
					extractor_displacement);
			for (unsigned int i = 0; i < extractor_slips.size(); i++) {
				data_out_tmp.add_data_vector(current_solution,
						extractor_slips[i]);
			}
		}

		// add results from taskarray
		for (unsigned int t = 0; t < tasks->size(); t++) {
			tasks->at(t).solution[*ID].scale(mass_matrix_postprocessing[*ID]);
			if (ID == domainIDs.begin()) {
				Assert(extractors.find(tasks->at(t).name) == extractors.end(),
						ExcInternalError());
				extractors.insert(
						pair<string, Extract_scalar_value<dim, nslip>>(
								tasks->at(t).name,
								Extract_scalar_value<dim, nslip>(
										tasks->at(t).name)));
				data_out.add_data_vector(tasks->at(t).solution[*ID],
						extractors.at(tasks->at(t).name));
			} else {
				data_out_tmp.add_data_vector(tasks->at(t).solution[*ID],
						extractors.at(tasks->at(t).name));
			}
		}

		//build and merge patches
		if (ID == domainIDs.begin()) {
			data_out.setID(*ID);
			data_out.build_patches(subdivisions);
		} else {
			data_out_tmp.setID(*ID);
			data_out_tmp.build_patches(subdivisions);
			data_out.merge_patches(data_out_tmp);
			data_out_tmp.clear_data_vectors();
		}
	}
//	debugFileOut.close();
}

template<int dim, int nslip> void Postprocessing<dim, nslip>::calculateMassMatrix() {
	hp::FEValues<dim> hp_fe_values(fe_collection, quadrature_collection,
			update_values | update_JxW_values);

	// search all DomainIDs
	typename hp::DoFHandler<dim>::active_cell_iterator cell =
			hp_dof_handler.begin_active(), endc = hp_dof_handler.end();
	for (; cell != endc; ++cell) {
		int domainID = cell->active_fe_index();
		if (find(domainIDs.begin(), domainIDs.end(), domainID)
				== domainIDs.end()) {
			domainIDs.push_back(domainID);
		}
	}

	mass_matrix_postprocessing.resize(domainIDs.size(),
			Vector<double>(hp_dof_handler.n_dofs()));
	for (vector<int>::const_iterator ID = domainIDs.begin();
			ID != domainIDs.end(); ++ID) {

		typename hp::DoFHandler<dim>::active_cell_iterator cell =
				hp_dof_handler.begin_active(), endc = hp_dof_handler.end();
		for (; cell != endc; ++cell) {
			if ((int) cell->active_fe_index() != *ID) {
				continue;
			}
			hp_fe_values.reinit(cell);
			const FEValues<dim> &fe_values =
					hp_fe_values.get_present_fe_values();
			const unsigned int dofs_per_cell = fe_values.dofs_per_cell;
			std::vector<types::global_dof_index> local_dof_indices(
					dofs_per_cell);
			cell->get_dof_indices(local_dof_indices);
			for (unsigned int q_point = 0;
					q_point < fe_values.n_quadrature_points; ++q_point) {
				double JxW = fe_values.JxW(q_point);
				for (unsigned int i = 0; i < dofs_per_cell; ++i) {
					double N_i = fe_values.shape_value(i, q_point);
					for (unsigned int j = 0; j < dofs_per_cell; j++) {
						if (fe_collection[cell->active_fe_index()].system_to_component_index(
								i).first
								== fe_collection[cell->active_fe_index()].system_to_component_index(
										j).first) {
							double N_j = fe_values.shape_value(j, q_point);
							mass_matrix_postprocessing[*ID](
									local_dof_indices[i]) += (N_i * N_j * JxW);
						}
					}
				}
			}
		}
		for (unsigned int i = 0; i < mass_matrix_postprocessing[*ID].size();
				i++) {
			mass_matrix_postprocessing[*ID](i) = 1
					/ mass_matrix_postprocessing[*ID](i);
		}

	}
}

template<int dim, int nslip> void Postprocessing<dim, nslip>::save_double_slip_diagrams(
		vector<vector<FEValuesExtractors::Scalar>> extractor_slips,
		const Vector<double> &solution) {

	if (nslip == 0){
		return;
	}

	Quadrature<dim> q(
			fe_collection[0].base_element(dim).get_unit_support_points());
	hp::QCollection<dim> q_collection(q);
	hp::FEValues<dim> hp_fe_values(fe_collection, q_collection,
			update_quadrature_points | update_values);

	// post_data[i][j] is a map<y-coord, slip-value system i> for crystal j
	vector<vector<map<double, double>>> post_data(nslip, vector<map<double, double>>(fe_collection.size()));
	vector<map<double, double>> equiv_pl_st(fe_collection.size(),
			map<double, double>());

	typename hp::DoFHandler<dim>::active_cell_iterator cell =
			hp_dof_handler.begin_active(), endc = hp_dof_handler.end();
	for (; cell != endc; ++cell) {
		hp_fe_values.reinit(cell);
		const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values();

//		Material<dim, nslip> *local_quadrature_points_data =
//				reinterpret_cast<Material<dim, nslip>*>(cell->user_pointer());

		const vector<Point<dim> > node_coords =
				fe_values.get_quadrature_points();

		vector<vector<double> > node_slips(nslip,
				vector<double>(fe_values.n_quadrature_points));

		for (unsigned int i = 0; i < extractor_slips[cell->active_fe_index()].size();
				i++) {
			fe_values[extractor_slips[cell->active_fe_index()][i]].get_function_values(
					solution, node_slips[i]);
		}

		for (unsigned int q_point = 0; q_point < node_coords.size(); q_point++) {
			if (node_coords[q_point][0] * node_coords[q_point][0]
					+ node_coords[q_point][2] * node_coords[q_point][2]
					< 1e-5) {
				SymmetricTensor<2, dim> equiv_pl_strain;
				for (int i = 0; i < nslip; i++) {
					post_data[i][cell->active_fe_index()].insert(
							pair<double, double>(node_coords[q_point][1],
									node_slips[i][q_point]));
//					equiv_pl_strain += node_slips[i][q_point]
					//							* local_quadrature_points_data->get_P(i);
					equiv_pl_strain += node_slips[i][q_point]
							* SymmetricTensor<2, dim>();

				}
				equiv_pl_st[cell->active_fe_index()].insert(
						pair<double, double>(node_coords[q_point][1],
								equiv_pl_strain.norm()));
			}
		}

	}

	const char* args[] = { "blue", "red", "green" };
	std::vector<std::string> colors(args, args + 3);

	{ // write out the values of the slipsystems
		std::stringstream sstm;
		sstm << "Results/double_slip_diagram.txt";
		string filename = sstm.str();

		std::ofstream File(filename, std::ios::out/* | std::ios::app*/);
		if (!File.is_open()) {
			return;
		}

		for (unsigned int i = 0; i < post_data.size(); i++) { // Slipsystem i
			for (unsigned int j = 0; j < post_data[i].size(); j++) { // Crystal j
				File << "\\addplot [color="
						<< colors[(i < colors.size() ? i : colors.size() - 1)]
						<< ", mark = none, ";
				if (j == 0) {
					File << "forget plot";
				}
				File << "]\n  table[row sep=newline]{%" << endl;
				map<double, double>::iterator it;
				for (it = post_data[i][j].begin(); it != post_data[i][j].end();
						it++) {
					File << setw(15) << it->second  // string (key)
							<< " " << setw(15) << it->first   // string's value
							<< endl;
				}
				File << "};\n" << endl;
				if (j == post_data[i].size() - 1) {
					File << "\\addlegendentry{Slipsystem " << i + 1
							<< ", $\\lambda_{gb}$ = xx};\n\n" << endl;
				}
			}
		}
	}   // END: write out the values of the slipsystems

	{   // write out the values of the equivalent plastic strain
		std::stringstream sstm;
		sstm << "Results/double_slip_equivalent_plastic_strain_diagram.txt";
		string filename = sstm.str();

		std::ofstream File(filename, std::ios::out/* | std::ios::app*/);
		if (!File.is_open()) {
			return;
		}

		for (unsigned int j = 0; j < equiv_pl_st.size(); j++) { // Crystal j
			File << "\\addplot [color=blue, mark = none, ";
			if (j != equiv_pl_st.size() - 1) {
				File << "forget plot";
			}
			File << "]\n  table[row sep=newline]{%" << endl;
			map<double, double>::iterator it;
			for (it = equiv_pl_st[j].begin(); it != equiv_pl_st[j].end();
					it++) {
				File << setw(15) << it->second  // string (key)
						<< " " << setw(15) << it->first   // string's value
						<< endl;
			}
			File << "};\n" << endl;
		}
		File
				<< "\\addlegendentry{$|$\\mbox{\\boldmath$\\varepsilon$}$^p|$, $\\lambda_{gb}$ = xx};\n\n"
				<< endl;
	}   // END: write out the values of the equivalent plastic strain
}
