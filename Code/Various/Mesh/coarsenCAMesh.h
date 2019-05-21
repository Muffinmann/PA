/*
 * coarsenCAMesh.h
 *
 *  Created on: Mar 29, 2018
 *      Author: iwtm84
 */

#ifndef CODE_VARIOUS_MESH_COARSENCAMESH_H_
#define CODE_VARIOUS_MESH_COARSENCAMESH_H_

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/grid_out.h>

#include <iostream>
#include <fstream>
#include <cmath>
#include <boost/lexical_cast.hpp>

using namespace dealii;
using namespace std;

namespace readCA {

template<int dim>
bool checkRefinementConstraintOnGrainBoundaries(
		Triangulation<dim> &triangulation) {
	bool refinedSomething = false;
	typename Triangulation<dim>::active_cell_iterator cell =
			triangulation.begin_active(), endc = triangulation.end();

	for (; cell != endc; ++cell) {
		int cell_domainID = cell->material_id();
		for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; ++f) {
			if (cell->at_boundary(f)) {
				continue; // If the face is at the boundary it cant be a grainboundary
			}
			int neighbor_domainID = cell->neighbor(f)->material_id();
			if (cell_domainID == neighbor_domainID) {
				continue; // If both have the same MaterialID the face is inside a grain
			}

			if (!cell->neighbor(f)->active()) {
				if (!cell->refine_flag_set()) {
//					cell->set_refine_flag();
					cell->set_refine_flag(RefinementCase<dim>::cut_xy);
					refinedSomething = true;
				}
			}
		}
	}
	return refinedSomething;
}

template<int dim>
bool setRefinementAndCoarseningFlags(Triangulation<dim> &triangulation) {
	int numberOfChangedFlags = 0;
	typename Triangulation<dim>::active_cell_iterator cell =
			triangulation.begin_active(), endc = triangulation.end();

	for (; cell != endc; ++cell) {
		// if cell is refined, loop over all childs of its parent
		if (cell->level() > 0) {
			bool cell_can_be_coarsened = true;
			// if all childs have the same materialID -> set all coarsening flags for childs
			unsigned int material_id_of_current_cell = cell->material_id();
			typename Triangulation<dim>::cell_iterator parent = cell->parent();

			/*Check for parents children if the all share the same material id as this cell*/
			for (unsigned int i = 0; i < parent->n_children(); ++i) {
				typename Triangulation<dim>::cell_iterator child =
						parent->child(i);
				if (child->active()) {
					if (material_id_of_current_cell != child->material_id()) {
						/*if they don't share the same material id the id is cleared*/
						cell_can_be_coarsened = false;
						break;
						/*at this point it is clear that "cell" can not be coarsened, i.e.
						 jump to next cell*/
					}
				}
			}

			// else delete all coarsening flags for childs
			if (!cell_can_be_coarsened) {
				for (unsigned int i = 0; i < parent->n_children(); ++i) {
					typename Triangulation<dim>::cell_iterator ptr_to_child =
							parent->child(i);
					if (ptr_to_child->active())
						ptr_to_child->clear_coarsen_flag();

				}
			} else {
				/*If this point is reached all children of the parent share the same material id. Therefore
				 * it is set here to that value*/
				parent->set_material_id(material_id_of_current_cell);
				cell->set_coarsen_flag();
				++numberOfChangedFlags;
			}
		}
	}
	return numberOfChangedFlags != 0;
}

template<int dim>
int countRefinementAndCoarseningFlags(Triangulation<dim> &triangulation) {
	int counter = 0;

	typename Triangulation<dim>::active_cell_iterator cell =
			triangulation.begin_active(), endc = triangulation.end();

	for (; cell != endc; ++cell) {
		if (cell->coarsen_flag_set() || cell->refine_flag_set()) {
			++counter;
		}
	}
	return counter;
}

template<int dim>
void changeRefinementFlagsToAniso(Triangulation<dim> &triangulation) {
	typename Triangulation<dim>::active_cell_iterator cell =
			triangulation.begin_active(), endc = triangulation.end();

	for (; cell != endc; ++cell) {
		if (cell->refine_flag_set()
				== RefinementCase<dim>::isotropic_refinement) {
			cell->set_refine_flag(RefinementCase<dim>::cut_xy);
		}
	}
}

template<int dim>
void outputMesh(Triangulation<dim> &triangulation, int refinementLoops) {
	std::string filename = "grid-";
	filename += ('0' + refinementLoops);
	filename += ".vtk";

	std::ofstream out(filename);
	GridOut grid_out;
	grid_out.write_vtk(triangulation, out);
}

template<int dim>
void coarsen_mesh(Triangulation<dim> &triangulation) {
	// coarsen as many cells as possible
	int maxRefinementLevel = triangulation.n_levels();
	int refinementLoops = 1;

//	outputMesh(triangulation, refinementLoops);

	// merge elements with the same material ID
	while (true) {
		bool a = setRefinementAndCoarseningFlags(triangulation);
		bool c = triangulation.prepare_coarsening_and_refinement();
		changeRefinementFlagsToAniso(triangulation);
		int numberOfFlags = countRefinementAndCoarseningFlags(triangulation);

		deallog << "Flags are set due to same MaterialID: " << a << endl;
		deallog << "Flags are changed by prepare_coarsening_and_refinement: "
				<< c << endl;
		deallog << "Number of refinement/coarsening flags: " << numberOfFlags
				<< endl;
		deallog << "-------------------------------------------" << endl;

		if (numberOfFlags == 0) {
			deallog << endl << "Coarsening due to material ID done" << endl;
			deallog
					<< "----------------------------------------------------------------------------------"
					<< endl;
			break;
		}

		triangulation.execute_coarsening_and_refinement();
		++refinementLoops;

//		outputMesh(triangulation, refinementLoops);

		if (refinementLoops > maxRefinementLevel) {
			break;
		}
	}

	// refine to hold the constraint of same refinement level on grain boundaries
	maxRefinementLevel *= 2;
	while (true) {
		bool b = checkRefinementConstraintOnGrainBoundaries(triangulation);
		bool c = triangulation.prepare_coarsening_and_refinement();
		changeRefinementFlagsToAniso(triangulation);
		int numberOfFlags = countRefinementAndCoarseningFlags(triangulation);

		deallog << "Flags are changed due to grain boundary constraint: " << b
				<< endl;
		deallog << "Flags are changed by prepare_coarsening_and_refinement: "
				<< c << endl;
		deallog << "Number of refinement/coarsening flags: " << numberOfFlags
				<< endl;
		deallog << "-------------------------------------------" << endl;

		if (numberOfFlags == 0) {
			deallog << endl << "Refinement due to grain boundaries done"
					<< endl;
			deallog
					<< "----------------------------------------------------------------------------------"
					<< endl;
			break;
		}

		triangulation.execute_coarsening_and_refinement();
		++refinementLoops;

//		outputMesh(triangulation, refinementLoops);

		if (refinementLoops > maxRefinementLevel) {
			deallog << endl << "Too many iterations" << endl;
			break;
		}
	}
}

vector<vector<int>> read_orientation_file(string filename) {
	string pathname = "Input/Meshs/caData/" + filename + ".txt";
	ifstream in(pathname);
	if (in) {
		int n_x = -1; // cells in x direction
		int n_y = -1; // cells in y direction
		int n_c = -1; // number of crystals
		string line;

		// read n_x, n_y, n_c
		while (getline(in, line)) {
			if (line.compare(0, 1, "%") == 0) {
				continue;
			}
			if (line.find("n_x") != string::npos) {
				n_x = boost::lexical_cast<int>(
						line.substr(4, line.length() - 4));
			}
			if (line.find("n_y") != string::npos) {
				n_y = boost::lexical_cast<int>(
						line.substr(4, line.length() - 4));
			}
			if (line.find("n_c") != string::npos) {
				n_c = boost::lexical_cast<int>(
						line.substr(4, line.length() - 4));
			}
			if (line.find("IDs:") != string::npos) {
				break;
			}
		}
		if (n_x < 0 || n_y < 0 || n_c < 0) {
			deallog << "n_x: " << n_x << endl;
			deallog << "n_y: " << n_y << endl;
			deallog << "n_c: " << n_c << endl;
			Assert(false, ExcMessage("ID_file header incorrect"));
		}

		// Read crystal IDs
		vector<vector<int>> ids(n_x, vector<int>(n_y));
		int i = 0;
		while (getline(in, line)) {
			if (line.compare(0, 1, "%") == 0) {
				continue;
			}
			int y = i % n_x;
			int x = i / n_x;
			ids[x][y] = boost::lexical_cast<int>(line);
			++i;
		}
		return ids;

	} else {
		deallog << "Filename: " << pathname << endl;
		Assert(false, ExcMessage("ID file not found."));
	}
	return vector<vector<int>>();
}

template<int dim>
void refine_global_xy(Triangulation<dim> &triangulation, int nRefinements) {
	for (int i = 0; i < nRefinements; ++i) {
		typename Triangulation<dim>::active_cell_iterator cell =
				triangulation.begin_active(), endc = triangulation.end();
		for (; cell != endc; ++cell) {
			cell->set_refine_flag(RefinementCase<dim>::cut_xy);
		}
		triangulation.execute_coarsening_and_refinement();
	}
}

//template<int dim>
//void Triangulation<dim>::refine_global_xy(int nRefinements){
//	for(int i=0;i<nRefinements;++i){
//	typename Triangulation<dim>::active_cell_iterator cell =
//				triangulation.begin_active(), endc = triangulation.end();
//		for (; cell != endc; ++cell) {
//			cell->set_refine_flag(RefinementCase<dim>::cut_xy);
//		}
//	}
//}

template<int dim>
void make_grid(Triangulation<dim> & triangulation, string filename,
		vector<double> &geometric_range, ParameterHandler &parameter) {
	deallog.push("MeshGen");

	vector<vector<int>> ids = read_orientation_file(filename);
	int refinementLevelCA = log2(ids.size());
	Assert(refinementLevelCA == log2(ids.size()),
			ExcMessage("Read id file has no pow of 2 as size"));

	// If linear CG, 5 elements in thinkness direction
	int rep = 5;
	if (parameter.get("loadtype").find("linear") != string::npos
			&& parameter.get("loadtype").find("CG") == string::npos) {
		rep = 5;
	} else {
		rep = 1;
	}

	geometric_range[0] = 1;
	geometric_range[1] = 1;
	int refinementLevel = refinementLevelCA + parameter.get_integer("n_refinements");
	geometric_range[2] = pow(0.5, refinementLevel) * rep;
	Point<dim> p1(0, 0, 0);
	Point<dim> p2(geometric_range[0], geometric_range[1], geometric_range[2]);
	const vector<unsigned int> repetitions = { (unsigned int) pow(2.0,
			refinementLevel), (unsigned int) pow(2.0, refinementLevel),
			(unsigned int) (rep) };
	GridGenerator::subdivided_hyper_rectangle(triangulation, repetitions, p1,
			p2, /* colorize =*/
			true);

//	triangulation.begin(0)->recursively_set_material_id(254);

//	outputMesh(triangulation, 0);

	int n_x = pow(2.0, refinementLevelCA);
	int n_y = n_x;
	AssertIndexRange(n_x, ids.size() + 1);
	AssertIndexRange(n_y, ids[0].size() + 1);

	// set material IDs
	typename Triangulation<dim>::active_cell_iterator cell =
			triangulation.begin_active(), endc = triangulation.end();
	for (; cell != endc; ++cell) {
		Point<dim> barycenter = cell->barycenter();
		int x = barycenter[0] * n_x;
		int y = barycenter[1] * n_y;
		//		deallog<< "x/y: " << x << "/" << y << endl;
		cell->set_material_id(ids[x][y]);
	}

//	std::ofstream out("ca_mesh.msh");
//	GridOut grid_out;
//	grid_out.write_msh(triangulation, out);

//	coarsen_mesh(triangulation);

	deallog << "Number of active cells: " << triangulation.n_active_cells()
			<< endl;
	deallog.pop();
}

} // end namespace readCA

#endif /* CODE_VARIOUS_MESH_COARSENCAMESH_H_ */
