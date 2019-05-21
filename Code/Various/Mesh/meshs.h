/*
 * meshs.h
 *
 *  Created on: Sep 29, 2017
 *      Author: iwtm84
 */

#ifndef CODE_VARIOUS_MESH_MESHS_H_
#define CODE_VARIOUS_MESH_MESHS_H_

#include "coarsenCAMesh.h"

/*
 * return layer of point in range [0; 4]
 * t = z/zMax
 * 0 for points with t<0.5 || t>3.5
 * 1 for points with t>0.5 && t<1.5
 * linear in between
 */
double getLayer(double z, double zMax) {
	double layer = (z / zMax) * 4 + 0.5;
	if (layer > 4) {
		layer -= 4;
	}
//	deallog << "z: " << z / zMax << " Layer: " << layer << endl;
	return layer;
}
template<int dim>
double getLayer(Point<dim> vertex, vector<double> geometric_range) {
	return getLayer(vertex[2], geometric_range[2]);
}

template<int dim>
void generateCrossSnake(ParameterHandler &parameter,
		hp::DoFHandler<dim> &dof_handler, vector<double> geometric_range) {

	if (!parameter.get_bool("crosssnake")) {
		return;
	}
	deallog << "Generate Crosssnake" << endl;

	double dx = geometric_range[2] / 4
			* tan(parameter.get_double("csAngle") / 180 * 3.14159265359);

	std::set<unsigned int> vertex_indices;
	typename Triangulation<dim>::active_cell_iterator cell =
			dof_handler.get_triangulation().begin_active(), endc =
			dof_handler.get_triangulation().end();

	for (; cell != endc; ++cell) {
		for (unsigned int i = 0; i < GeometryInfo<dim>::vertices_per_cell;
				++i) {

			unsigned int vertex_index = cell->vertex_index(i);
			if (vertex_indices.count(vertex_index) == 0) {
				vertex_indices.insert(vertex_index);
				Point<dim> &v = cell->vertex(i);
				double layer = getLayer(v, geometric_range);

				if (layer < 1 - 1e-8) {
					v(0) += dx * layer;
				} else {
					v(0) += dx;
					if (layer < 2 - 1e-8) {
						v(1) += dx * (layer - 1);
					} else {
						v(1) += dx;
						if (layer < 3 - 1e-8) {
							v(0) -= dx * (layer - 2);
						} else {
							v(0) -= dx;
							v(1) -= dx * (layer - 3);
						}
					}
				}
				v(0) -= dx / 2;
			}
		}
	}
}

namespace testcases {

template<int dim>
void setMaterialIDsBD(ParameterHandler &parameter,
		hp::DoFHandler<dim> &dof_handler) {

	typename hp::DoFHandler<dim>::active_cell_iterator cell =
			dof_handler.begin_active(), endc = dof_handler.end();

	for (; cell != endc; ++cell) {
		Point<dim> center = cell->barycenter();
		int numberOfCrystals = parameter.get_integer("n_crystals");
		types::material_id domainID = center(1) * numberOfCrystals;
		cell->set_material_id(domainID);
	}
}

template<int dim>
void setMaterialIDs(ParameterHandler &parameter,
		hp::DoFHandler<dim> &dof_handler, double zMax) {

	typename hp::DoFHandler<dim>::active_cell_iterator cell =
			dof_handler.begin_active(), endc = dof_handler.end();
	for (; cell != endc; ++cell) {
		Point<dim> center = cell->barycenter();
		int numberOfCrystals = parameter.get_integer("n_crystals");
		// Calculate domainID based on coordinates
		// Use modulo operator and a shift of half a crystal to avoid
		// grain boundaries on the geometric boundary of the cube
		int x = center(0) * numberOfCrystals + 0.5;
		int y = center(1) * numberOfCrystals + 0.5;
		int z = center(2) * numberOfCrystals + 0.5;
		types::material_id domainID;

		if (parameter.get_bool("isotrop")) {
			domainID = pow(numberOfCrystals, 2) * (x % numberOfCrystals)
					+ pow(numberOfCrystals, 1) * (y % numberOfCrystals)
					+ pow(numberOfCrystals, 0) * (z % numberOfCrystals);
		} else {
			if (parameter.get_bool("crosssnake")) {
				/*
				 * like for the else case but 4 consecutive IDs for the 4 parts of a crystal
				 */
				int layer = getLayer(center(2), zMax);
				domainID = (pow(numberOfCrystals, 1) * (x % numberOfCrystals)
						+ pow(numberOfCrystals, 0) * (y % numberOfCrystals)) * 4
						+ layer;
//				deallog << "Layer: " << layer << " x: " << x << " y: " << y
//						<< " DomainID: " << (int) domainID << endl;
			} else {
				/*
				 * x=0, y=0 -> ID 0
				 * x=0, y=1 -> ID 1
				 * x=1, y=0 -> ID numberOfCrystals
				 */
				domainID = pow(numberOfCrystals, 1) * (x % numberOfCrystals)
						+ pow(numberOfCrystals, 0) * (y % numberOfCrystals);
			}
		}
		cell->set_material_id(domainID);
	}
}

namespace doubleSlip {

template<int dim>
void createMesh(ParameterHandler &parameter, Triangulation<dim> &triangulation,
		hp::DoFHandler<dim> &dof_handler) {
	unsigned int cells_in_other_directions = 1;
	double edgeLenght = cells_in_other_directions
			* pow(0.5, parameter.get_integer("n_refinements"));
	Point<dim> p1(0, 0, 0);
	Point<dim> p2(edgeLenght, 1, edgeLenght);
	const vector<unsigned int> repetitions = { cells_in_other_directions,
			(unsigned int) pow(2.0, parameter.get_integer("n_refinements")),
			cells_in_other_directions };
	GridGenerator::subdivided_hyper_rectangle(triangulation, repetitions, p1,
			p2, /* colorize =*/
			true);
	// Set Material IDs
	setMaterialIDsBD(parameter, dof_handler);
}

} // end namespace doubleSlip

namespace gottschalk {

template<int dim>
void createMesh(ParameterHandler &parameter, Triangulation<dim> &triangulation,
		hp::DoFHandler<dim> &dof_handler) {
	Point<dim> p1(0, 0, 0);
	Point<dim> p2(1, 1, 1);
	const vector<unsigned int> repetitions = { 10, 40, 8 };
	GridGenerator::subdivided_hyper_rectangle(triangulation, repetitions, p1,
			p2, /* colorize =*/
			true);

	// Set Material IDs
	setMaterialIDsBD(parameter, dof_handler);
}

} // end namespace gottschalk

namespace GAMM16 {

template<int dim>
void createMesh(ParameterHandler &parameter, Triangulation<dim> &triangulation,
		hp::DoFHandler<dim> &dof_handler, vector<double> &geometric_range) {
	// generate a cube like the one for gamm16 but just one layer

	if (parameter.get_bool("crosssnake")) {
		int elementsPerLayer = parameter.get_integer("elementsCSLayer");
		double totalHeight = 4 * parameter.get_double("csHeight");
		geometric_range[2] = totalHeight;
		Point<dim> p1(0, 0, 0);
		Point<dim> p2(1, 1, geometric_range[2]);
		const vector<unsigned int> repetitions = { (unsigned int) pow(2.0,
				parameter.get_integer("n_refinements")), (unsigned int) pow(2.0,
				parameter.get_integer("n_refinements")), (unsigned int) (4
				* elementsPerLayer) };
		GridGenerator::subdivided_hyper_rectangle(triangulation, repetitions,
				p1, p2, /* colorize =*/
				true);
	} else {
		geometric_range[2] = pow(0.5, parameter.get_integer("n_refinements"));
		Point<dim> p1(0, 0, 0);
		Point<dim> p2(1, 1, geometric_range[2]);
		const vector<unsigned int> repetitions = { (unsigned int) pow(2.0,
				parameter.get_integer("n_refinements")), (unsigned int) pow(2.0,
				parameter.get_integer("n_refinements")), (unsigned int) (1) };
		GridGenerator::subdivided_hyper_rectangle(triangulation, repetitions,
				p1, p2, /* colorize =*/
				true);
	}

	setMaterialIDs(parameter, dof_handler, geometric_range[2]);
	generateCrossSnake(parameter, dof_handler, geometric_range);
}

} // end namespace GAMM16

} // end namespace testcases

namespace RVE {

template<int dim>
void createMesh(ParameterHandler &parameter, Triangulation<dim> &triangulation,
		hp::DoFHandler<dim> &dof_handler, vector<double> &geometric_range) {

	Triangulation<dim - 1> tmp_triangulation;

	if (parameter.get("mesh").find("ca") != string::npos) { // generate mesh from CA data
		readCA::make_grid(triangulation, parameter.get("mesh"), geometric_range,
				parameter);
	} else { // read 2 dimensional .msh file
		GridIn<dim - 1> gridin;
		gridin.attach_triangulation(tmp_triangulation);

		if (parameter.get("mesh").find(".inp") != string::npos) {
			// read abaqus file
			string filename;
			{
				std::stringstream sstm;
				sstm << "Input/Meshs/" << parameter.get("mesh");
				filename = sstm.str();
			}
			std::ifstream f(filename);
			gridin.read_abaqus(f);
		} else {
			// Read 2 dimensional mesh
			string filename;
			{
				std::stringstream sstm;
				sstm << "Input/Meshs/" << parameter.get("mesh") << ".msh";
				filename = sstm.str();
			}
			std::ifstream f(filename);
			gridin.read_msh(f);
		}

		tmp_triangulation.refine_global(parameter.get_double("n_refinements"));

		//extrude in third dimension
		if (parameter.get_bool("crosssnake")) {
			int elementsPerLayer = parameter.get_integer("elementsCSLayer");
			double totalHeight = 4 * parameter.get_double("csHeight");
			geometric_range[2] = totalHeight;
			GridGenerator::extrude_triangulation(tmp_triangulation,
					4 * elementsPerLayer + 1, geometric_range[2],
					triangulation);
		} else {
			geometric_range[2] = 1 / sqrt(tmp_triangulation.n_active_cells());
			GridGenerator::extrude_triangulation(tmp_triangulation, 2,
					geometric_range[2], triangulation);
		}
	}

	// find min max of geometry
	{
		double xmin = numeric_limits<double>::max();
		double ymin = numeric_limits<double>::max();
		double zmin = numeric_limits<double>::max();
		double xmax = -numeric_limits<double>::max();
		double ymax = -numeric_limits<double>::max();
		double zmax = -numeric_limits<double>::max();
		typename Triangulation<dim>::active_cell_iterator cell =
				triangulation.begin_active(), endc = triangulation.end();
		for (; cell != endc; ++cell) {
			for (unsigned int i = 0; i < GeometryInfo<dim>::faces_per_cell; i++) {
				if (cell->face(i)->at_boundary()) {
					Point<dim> center = cell->face(i)->center();
					if (center(0) < xmin) {
						xmin = center(0);
					}
					if (center(0) > xmax) {
						xmax = center(0);
					}
					if (center(1) < ymin) {
						ymin = center(1);
					}
					if (center(1) > ymax) {
						ymax = center(1);
					}
					if (center(2) < zmin) {
						zmin = center(2);
					}
					if (center(2) > zmax) {
						zmax = center(2);
					}
				}
			}
		}
//		deallog << "x [" << xmin << " - " << xmax << "], ";
//		deallog << "y [" << ymin << " - " << ymax << "], ";
//		deallog << "z [" << zmin << " - " << zmax << "]" << endl;
//		deallog << geometric_range[0] << " " << geometric_range[1] << " " << geometric_range[2] << endl;

		GridTools::shift(Tensor<1, dim>({-xmin, -ymin, -zmin}), triangulation);

		if (abs(xmax - xmin - 1) > 1e-6 || abs(ymax - ymin - 1) > 1e-6) {
			deallog << "Wrong size of the geometry" << endl;
			deallog << "x [" << xmin << " - " << xmax << "], ";
			deallog << "y [" << ymin << " - " << ymax << "], ";
			deallog << "z [" << zmin << " - " << zmax << "]" << endl;
		}
	}

	// set BoundaryIDs
	typename Triangulation<dim>::active_cell_iterator cell =
			triangulation.begin_active(), endc = triangulation.end();
	for (; cell != endc; ++cell) {
		for (unsigned int i = 0; i < GeometryInfo<dim>::faces_per_cell; i++) {
			if (cell->face(i)->at_boundary()) {
				Point<dim> center = cell->face(i)->center();
				if (abs(center(0) - 0) < 1e-5) {
//					deallog << "id0 " << center(0) << endl;
					cell->face(i)->set_boundary_id(0);

				} else if (abs(center(0) - geometric_range[0]) < 1e-5) {
//					deallog << "id1 " << center(0) << endl;
					cell->face(i)->set_boundary_id(1);

				} else if (abs(center(1) - 0) < 1e-5) {
//					deallog << "id2 " << center(1) << endl;
					cell->face(i)->set_boundary_id(2);

				} else if (abs(center(1) - geometric_range[1]) < 1e-5) {
//					deallog << "id3 " << center(1) << endl;
					cell->face(i)->set_boundary_id(3);

				} else if (abs(center(2) - 0) < 1e-5) {
//					deallog << "id4 " << center(2) << endl;
					cell->face(i)->set_boundary_id(4);

				} else if (abs(center(2) - geometric_range[2]) < 1e-5) {
//					deallog << "id5 " << center(2) << endl;
					cell->face(i)->set_boundary_id(5);

				}
//				else {
//					if (cell->face(i)->boundary_id()
//							!= numbers::internal_face_boundary_id) {
//						deallog << "internal " << center(0) << " " << center(1) << " " << center(2) << endl;
//						cell->face(i)->set_boundary_id(
//								numbers::internal_face_boundary_id);
//					}
//				}
			}
		}
		if (parameter.get_bool("crosssnake")) {
			int layer = getLayer(cell->center(), geometric_range);
			int newMatID = cell->material_id() * 4 + layer;
			cell->set_material_id(newMatID);
		}
	}

	// generate Crosssnake by moving vertices
	generateCrossSnake(parameter, dof_handler, geometric_range);
	deallog << "z: " << geometric_range[2] << endl;
}

} // end namespace RVE

#endif /* CODE_VARIOUS_MESH_MESHS_H_ */
