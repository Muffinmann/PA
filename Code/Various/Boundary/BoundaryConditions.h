/*
 * BoundaryConditions.h
 *
 *  Created on: Mar 14, 2017
 *      Author: iwtm84
 */

#ifndef CODE_VARIOUS_BOUNDARYCONDITIONS_H_
#define CODE_VARIOUS_BOUNDARYCONDITIONS_H_

#include "Boundary_Values.h"
#include "WrapperBoundaryValues.h"
using namespace std;
using namespace dealii;

template<int dim>
inline bool isReference(int dof, vector<vector<int> > &dofsRef) {
	Assert(dofsRef.size() == dim + 1, ExcInternalError());
	for (unsigned int i = 0; i < dofsRef.size(); i++) {
		Assert(dofsRef[i].size() == dim, ExcInternalError());
		for (unsigned int j = 0; j < dofsRef[i].size(); j++) {
			if (dofsRef[i][j] == dof) {
				return true;
			}
		}
	}
	return false;
}

inline void printConstraintsOnDof(ConstraintMatrix &constraints,
		unsigned int dof) {
	if (constraints.is_constrained(dof)) {
		const vector<pair<unsigned int, double> > *entries =
				constraints.get_constraint_entries(dof);
//		deallog << " dof: " << dof << ":";
		for (unsigned int k = 0; k < entries->size(); ++k) {
			deallog << " " << entries->at(k).first << "/"
					<< entries->at(k).second;
		}
		if (constraints.is_inhomogeneously_constrained(dof)) {
			deallog << " + " << constraints.get_inhomogeneity(dof);
		} else {
			deallog << " + 0";
		}
	}
}

void printConstraintsOnReferenceNodes(vector<vector<int> > &dof_reference_node,
		ConstraintMatrix &constraints) {

	// Print out constraints on reference nodes
	deallog << endl;
	for (unsigned int i = 0; i < dof_reference_node.size(); ++i) {
		for (unsigned int j = 0; j < dof_reference_node[i].size(); ++j) {
			deallog << "i/j: " << (int) i - 1 << "/" << j << " dof: "
					<< dof_reference_node[i][j] << ":";
			printConstraintsOnDof(constraints, dof_reference_node[i][j]);
			deallog << endl;
		}
	}
}

namespace BC {

/*
 * Sets diriclet BC for the node origin in dof_reference_node with the value from boundary
 */
template<int dim>
void constrainOrigin(BoundaryValues<dim> *boundary,
		std::map<int, std::pair<Point<dim, double>, int> > &dofdata,
		vector<vector<int> > &dof_reference_node,
		ConstraintMatrix &constraints) {

	// calculate value for the diriclet bc
	Vector<double> values(dim);
	boundary->vector_value(dofdata[dof_reference_node[0][0]].first, values);

	// add diriclet bc to constraint matrix
	for (unsigned int l = 0; l < dof_reference_node[0].size(); l++) {
		constraints.add_line(dof_reference_node[0][l]);
		constraints.set_inhomogeneity(dof_reference_node[0][l], values[l]);
	}
}

/*
 * Sets diriclet BC for reference Node X to the values in values in the directions specified in dirToConstrain
 */
void constrainRefX(vector<vector<int> > &dof_reference_node,
		ConstraintMatrix &constraints, Vector<double> values,
		vector<bool> &dirToConstrain) {

//	deallog << "Constraint RefX" << endl;

// add diriclet bc to constraint matrix
	for (unsigned int l = 0; l < dof_reference_node[1].size(); l++) {
//		deallog << "Schleife: " << l << " bool is: " << dirToConstrain[l] << endl;
		if (dirToConstrain[l]) {
			constraints.add_line(dof_reference_node[1][l]);
			constraints.set_inhomogeneity(dof_reference_node[1][l], values[l]);
//			deallog << "Constrain X in: " << l << " with value: " << values[l] << endl;
		}
	}
}

/*
 * Sets diriclet BC for reference Node X to the values in values in the directions specified in dirToConstrain
 */
template<int dim>
void constrainRefX(BoundaryValues<dim> *boundary,
		std::map<int, std::pair<Point<dim, double>, int> > &dofdata,
		vector<vector<int> > &dof_reference_node,
		ConstraintMatrix &constraints) {

	// calculate value for the diriclet bc
	Vector<double> values(dim);
	boundary->vector_value(dofdata[dof_reference_node[1][0]].first, values);

	vector<bool> dirToConstrain(dim, true);

	constrainRefX(dof_reference_node, constraints, values, dirToConstrain);
}

/*
 * Sets diriclet BC for reference Node X to the values in values in x direction
 */
template<int dim>
void constrainRefX_x(BoundaryValues<dim> *boundary,
		std::map<int, std::pair<Point<dim, double>, int> > &dofdata,
		vector<vector<int> > &dof_reference_node,
		ConstraintMatrix &constraints) {

	// calculate value for the diriclet bc
	Vector<double> values(dim);
	boundary->vector_value(dofdata[dof_reference_node[1][0]].first, values);

	vector<bool> dirToConstrain(false, dim);
	dirToConstrain[0] = true;
//	deallog << "Bool is " << dirToConstrain[0] << "/" <<dirToConstrain[1] << "/" <<dirToConstrain[2] << endl;

	constrainRefX(dof_reference_node, constraints, values, dirToConstrain);
}

/*
 * Sets diriclet BC for reference Node X to the values in values in y direction
 */
template<int dim>
void constrainRefX_y(BoundaryValues<dim> *boundary,
		std::map<int, std::pair<Point<dim, double>, int> > &dofdata,
		vector<vector<int> > &dof_reference_node,
		ConstraintMatrix &constraints) {

	// calculate value for the diriclet bc
	Vector<double> values(dim);
	boundary->vector_value(dofdata[dof_reference_node[1][0]].first, values);

	vector<bool> dirToConstrain(false, dim);
	dirToConstrain[1] = true;

	constrainRefX(dof_reference_node, constraints, values, dirToConstrain);
}

/*
 * Sets diriclet BC for reference Node X to the values in values in z direction
 */
template<int dim>
void constrainRefX_z(BoundaryValues<dim> *boundary,
		std::map<int, std::pair<Point<dim, double>, int> > &dofdata,
		vector<vector<int> > &dof_reference_node,
		ConstraintMatrix &constraints) {

	// calculate value for the diriclet bc
	Vector<double> values(dim);
	boundary->vector_value(dofdata[dof_reference_node[1][0]].first, values);

	vector<bool> dirToConstrain(false, dim);
	dirToConstrain[2] = true;

	constrainRefX(dof_reference_node, constraints, values, dirToConstrain);
}

///*
// * Sets diriclet BC for reference Node X to the values in values in y and z direction
// */
//template<int dim>
//void constrainRefX_yz(BoundaryValues<dim> *boundary,
//		std::map<int, std::pair<Point<dim, double>, int> > &dofdata,
//		vector<vector<int> > &dof_reference_node,
//		ConstraintMatrix &constraints) {
//
//	// calculate value for the diriclet bc
//	Vector<double> values(dim);
//	boundary->vector_value(dofdata[dof_reference_node[1][0]].first, values);
//
//	vector<bool> dirToConstrain(false, dim);
//	dirToConstrain[1] = true;
//	dirToConstrain[2] = true;
//
//	constrainRefX(dof_reference_node, constraints, values, dirToConstrain);
//}

/*
 * Sets diriclet BC for reference Node Y to the values in values in the directions specified in dirToConstrain
 */
void constrainRefY(vector<vector<int> > &dof_reference_node,
		ConstraintMatrix &constraints, Vector<double> values,
		vector<bool> &dirToConstrain) {

	// add diriclet bc to constraint matrix
	for (unsigned int l = 0; l < dof_reference_node[2].size(); l++) {
		if (dirToConstrain[l]) {
			constraints.add_line(dof_reference_node[2][l]);
			constraints.set_inhomogeneity(dof_reference_node[2][l], values[l]);
		}
	}
}

/*
 * Sets diriclet BC for reference Node Y to the values got from boundary
 */
template<int dim>
void constrainRefY(BoundaryValues<dim> *boundary,
		std::map<int, std::pair<Point<dim, double>, int> > &dofdata,
		vector<vector<int> > &dof_reference_node,
		ConstraintMatrix &constraints) {

	// calculate value for the diriclet bc
	Vector<double> values(dim);
	boundary->vector_value(dofdata[dof_reference_node[2][0]].first, values);

	vector<bool> dirToConstrain(dim, true);

	constrainRefY(dof_reference_node, constraints, values, dirToConstrain);
}

/*
 * Sets diriclet BC for reference Node Y to the values in values in y direction
 */
template<int dim>
void constrainRefY_y(BoundaryValues<dim> *boundary,
		std::map<int, std::pair<Point<dim, double>, int> > &dofdata,
		vector<vector<int> > &dof_reference_node,
		ConstraintMatrix &constraints) {

	// calculate value for the diriclet bc
	Vector<double> values(dim);
	boundary->vector_value(dofdata[dof_reference_node[2][0]].first, values);

	vector<bool> dirToConstrain(false, dim);
	dirToConstrain[1] = true;

	constrainRefY(dof_reference_node, constraints, values, dirToConstrain);
}

/*
 * Sets diriclet BC for reference Node Y to the values in values in z direction
 */
template<int dim>
void constrainRefY_z(BoundaryValues<dim> *boundary,
		std::map<int, std::pair<Point<dim, double>, int> > &dofdata,
		vector<vector<int> > &dof_reference_node,
		ConstraintMatrix &constraints) {

	// calculate value for the diriclet bc
	Vector<double> values(dim);
	boundary->vector_value(dofdata[dof_reference_node[2][0]].first, values);

	vector<bool> dirToConstrain(false, dim);
	dirToConstrain[2] = true;

	constrainRefY(dof_reference_node, constraints, values, dirToConstrain);
}

/*
 * Sets diriclet BC for reference Node Z to the values in values in the directions specified in dirToConstrain
 */
void constrainRefZ(vector<vector<int> > &dof_reference_node,
		ConstraintMatrix &constraints, Vector<double> values,
		vector<bool> &dirToConstrain) {

//	deallog << "Constraint RefZ" << endl;

// add diriclet bc to constraint matrix
	for (unsigned int l = 0; l < dof_reference_node[3].size(); l++) {
		if (dirToConstrain[l]) {
			constraints.add_line(dof_reference_node[3][l]);
			constraints.set_inhomogeneity(dof_reference_node[3][l], values[l]);
//			deallog << "Constrain Z in: " << l << " with value: " << values[l] << endl;
		}
	}
}

/*
 * Sets diriclet BC for reference Node Z to the values in values in the directions specified in dirToConstrain
 */
template<int dim>
void constrainRefZ(BoundaryValues<dim> *boundary,
		std::map<int, std::pair<Point<dim, double>, int> > &dofdata,
		vector<vector<int> > &dof_reference_node,
		ConstraintMatrix &constraints) {

	// calculate value for the diriclet bc
	Vector<double> values(dim);
	boundary->vector_value(dofdata[dof_reference_node[3][0]].first, values);

	vector<bool> dirToConstrain(dim, true);

	constrainRefZ(dof_reference_node, constraints, values, dirToConstrain);
}


/*
 * Sets diriclet BC for reference Node Z to the values in values in z direction
 */
template<int dim>
void constrainRefZ_z(BoundaryValues<dim> *boundary,
		std::map<int, std::pair<Point<dim, double>, int> > &dofdata,
		vector<vector<int> > &dof_reference_node,
		ConstraintMatrix &constraints) {

	// calculate value for the diriclet bc
	Vector<double> values(dim);
	boundary->vector_value(dofdata[dof_reference_node[3][0]].first, values);

	vector<bool> dirToConstrain(false, dim);
	dirToConstrain[2] = true;

	constrainRefZ(dof_reference_node, constraints, values, dirToConstrain);
}

template<int dim>
void microhardTopBottom(hp::DoFHandler<dim> &dof_handler,
		ConstraintMatrix &periodic_constraints, int nslip, int maxCrystals) {
	std::vector<bool> components_DBC_slips(dim + nslip * maxCrystals, true);
	for (int i = 0; i < dim; i++) {
		components_DBC_slips[i] = false;
	}
	ComponentMask mask_DBC_slips(components_DBC_slips);
	for (int i = 2; i < 4; i++) { // Boundary ID 2&3 are top and bottom
		DoFTools::make_zero_boundary_constraints(dof_handler, i,
				periodic_constraints, mask_DBC_slips);
	}
}

template<int dim>
void constrainMacroscopicallySymmetric(
		std::map<int, std::pair<Point<dim, double>, int> > &dofdata,
		vector<vector<int> > &dof_reference_node,
		ConstraintMatrix &constraints) {

	vector<double> geometricRange(dim);
	for (unsigned int i = 0; i < dim; ++i) {
		geometricRange[i] = (dofdata[dof_reference_node[i + 1][0]].first
				- dofdata[dof_reference_node[0][0]].first)[i];
//		deallog << "Size[" << i << "]: " << geometricRange[i] << endl;
	}

	// add diriclet bc to constraint matrix
	// u_x(ref_y) = u_y(ref_x) * L_y/L_x
	constraints.add_line(dof_reference_node[2][0]);
	constraints.add_entry(dof_reference_node[2][0], dof_reference_node[1][1],
			geometricRange[1] / geometricRange[0]);
	// u_x(ref_z) = u_z(ref_x) * L_z/L_x
	constraints.add_line(dof_reference_node[3][0]);
	constraints.add_entry(dof_reference_node[3][0], dof_reference_node[1][2],
			geometricRange[2] / geometricRange[0]);
	// u_y(ref_z) = u_z(ref_y) * L_z/L_y
	constraints.add_line(dof_reference_node[3][1]);
	constraints.add_entry(dof_reference_node[3][1], dof_reference_node[2][2],
			geometricRange[2] / geometricRange[1]);
}

namespace Periodic {

/*
 * Apply constraints on one face
 */
template<int dim>
void constrainFace(std::map<int, std::pair<Point<dim, double>, int> > &dofdata,
		vector<vector<int> > &dof_reference_node,
		ConstraintMatrix &periodic_constraints,
		hp::DoFHandler<dim> &dof_handler, int constraint_facedir,
		vector<bool> &constraint_dir) {

	ConstraintMatrix tmp_constraints;
	tmp_constraints.clear();

	//make periodicity constraints in one direction for displacement and slip-dofs in tmp_constraints
	DoFTools::make_periodicity_constraints<hp::DoFHandler<dim> >(dof_handler,
			2 * constraint_facedir + 1, 2 * constraint_facedir,
			constraint_facedir, tmp_constraints);

	// loop over all dof
	for (unsigned int i = 0; i < dof_handler.n_dofs(); i++) {
		// if dof is not constrained -> next dof
		if (!tmp_constraints.is_constrained(i)) {
			continue;
		}
		if (isReference<dim>(i, dof_reference_node)) {
			continue; // Referenzknoten werden nicht constrained
		}
		// Add displacement of reference node as inhomogeneity to the constraint
		// u = u_periodic_partner + u_ref
		if (dofdata.count(i) != 0) { // constrain dof wit dof_ref
			if (!constraint_dir[dofdata[i].second]) {
				continue;
			}
			tmp_constraints.add_entry(i,
					dof_reference_node[constraint_facedir + 1][dofdata[i].second],
					1);
		}

		// Merge the constraints from tmp_constraints in periodic constraints
		// dealii function does not work because dofs on edges are constraint
		// from different other edges
		const std::vector<std::pair<unsigned int, double> > *periodic_vector =
				tmp_constraints.get_constraint_entries(i);
		if (!periodic_constraints.is_constrained(i)) {
			for (unsigned int j = 0; j < periodic_vector->size(); j++) {
				periodic_constraints.add_line(i);
				periodic_constraints.add_entry(i, periodic_vector->at(j).first,
						periodic_vector->at(j).second);
			}
		}
	}
}

/*
 * Apply constraints on all relevant faces
 */
template<int dim>
void constrainFaces(std::map<int, std::pair<Point<dim, double>, int> > &dofdata,
		vector<vector<int> > &dof_reference_node, ConstraintMatrix &constraints,
		hp::DoFHandler<dim> &dof_handler, vector<bool> &faces,
		vector<bool> &directions) {

	ConstraintMatrix periodic_constraints;
	periodic_constraints.clear();
	for (int constraint_dir = 0; constraint_dir < dim; constraint_dir++) {
		if (faces[constraint_dir]) {
			constrainFace(dofdata, dof_reference_node, periodic_constraints,
					dof_handler, constraint_dir, directions);
		}
	}
	periodic_constraints.close();
	constraints.merge(periodic_constraints, ConstraintMatrix::left_object_wins);
}

/*
 *  Apply constraints for periodic BC
 */
template<int dim>
void constrain(BoundaryValues<dim> *boundary,
		std::map<int, std::pair<Point<dim, double>, int> > &dofdata,
		vector<vector<int> > &dof_reference_node, ConstraintMatrix &constraints,
		hp::DoFHandler<dim> &dof_handler) {

	BC::constrainOrigin(boundary, dofdata, dof_reference_node, constraints);
	BC::constrainRefX(boundary, dofdata, dof_reference_node, constraints);
	BC::constrainRefY(boundary, dofdata, dof_reference_node, constraints);
	BC::constrainRefZ(boundary, dofdata, dof_reference_node, constraints);

	// Constrain faces in all directions
	vector<bool> faces(dim, true);
	vector<bool> directions(dim, true);
	Periodic::constrainFaces(dofdata, dof_reference_node, constraints,
			dof_handler, faces, directions);

}

/*
 *  Apply constraints for periodic uniaxial tension BC
 */
template<int dim>
void constrainUniaxialTension(BoundaryValues<dim> *boundary,
		std::map<int, std::pair<Point<dim, double>, int> > &dofdata,
		vector<vector<int> > &dof_reference_node, ConstraintMatrix &constraints,
		hp::DoFHandler<dim> &dof_handler, string load) {

	BC::constrainOrigin(boundary, dofdata, dof_reference_node, constraints);

	if (load.compare("11") == 0) {
		BC::constrainRefX_x(boundary, dofdata, dof_reference_node, constraints);
	} else if (load.compare("22") == 0) {
		BC::constrainRefY_y(boundary, dofdata, dof_reference_node, constraints);
	} else if (load.compare("33") == 0) {
		BC::constrainRefZ_z(boundary, dofdata, dof_reference_node, constraints);
	} else if (load.compare("12") == 0) {
		BC::constrainRefX_y(boundary, dofdata, dof_reference_node, constraints);
	} else if (load.compare("13") == 0) {
		BC::constrainRefX_z(boundary, dofdata, dof_reference_node, constraints);
	} else if (load.compare("23") == 0) {
		BC::constrainRefY_z(boundary, dofdata, dof_reference_node, constraints);
	} else {
		AssertThrow(false, ExcMessage("Wrong straindir."));
	}

	BC::constrainMacroscopicallySymmetric(dofdata, dof_reference_node,
			constraints);

	// Constrain faces in all directions
	vector<bool> faces(dim, true);
	vector<bool> directions(dim, true);
	Periodic::constrainFaces(dofdata, dof_reference_node, constraints,
			dof_handler, faces, directions);

}
/*
 *  Apply constraints for periodic plane stress BC
 */
template<int dim>
void constrainPlaneStress(BoundaryValues<dim> *boundary,
		std::map<int, std::pair<Point<dim, double>, int> > &dofdata,
		vector<vector<int> > &dof_reference_node, ConstraintMatrix &constraints,
		hp::DoFHandler<dim> &dof_handler, vector<string> &loads) {

	BC::constrainOrigin(boundary, dofdata, dof_reference_node, constraints);

//	BC::constrainRefX_y(boundary, dofdata, dof_reference_node, constraints);
//	BC::constrainRefX_z(boundary, dofdata, dof_reference_node, constraints);
//	BC::constrainRefZ_y(boundary, dofdata, dof_reference_node, constraints);

//	BC::constrainRefX_x(boundary, dofdata, dof_reference_node, constraints);
//	BC::constrainRefZ_z(boundary, dofdata, dof_reference_node, constraints);

	for (unsigned int i = 0; i < loads.size(); ++i) {
		if (loads[i].compare("11") == 0) {
			BC::constrainRefX_x(boundary, dofdata, dof_reference_node,
					constraints);
		} else if (loads[i].compare("22") == 0) {
			BC::constrainRefY_y(boundary, dofdata, dof_reference_node,
					constraints);
		} else if (loads[i].compare("33") == 0) {
			BC::constrainRefZ_z(boundary, dofdata, dof_reference_node,
					constraints);
		} else if (loads[i].compare("12") == 0) {
			BC::constrainRefX_y(boundary, dofdata, dof_reference_node,
					constraints);
		} else if (loads[i].compare("13") == 0) {
			BC::constrainRefX_z(boundary, dofdata, dof_reference_node,
					constraints);
		} else if (loads[i].compare("23") == 0) {
			BC::constrainRefY_z(boundary, dofdata, dof_reference_node,
					constraints);
		} else {
			AssertThrow(false, ExcMessage("Wrong straindir."));
		}
	}

	BC::constrainMacroscopicallySymmetric(dofdata, dof_reference_node,
			constraints);

// Constrain faces in all directions
	vector<bool> faces(dim, true);
	vector<bool> directions(dim, true);
	Periodic::constrainFaces(dofdata, dof_reference_node, constraints,
			dof_handler, faces, directions);

}
///*
// *  Apply constraints for periodic plane stress BC in XY plane
// */
//template<int dim>
//void constrainPlaneStressXY(BoundaryValues<dim> *boundary,
//		std::map<int, std::pair<Point<dim, double>, int> > &dofdata,
//		vector<vector<int> > &dof_reference_node, ConstraintMatrix &constraints,
//		hp::DoFHandler<dim> &dof_handler) {
//
//	BC::constrainOrigin(boundary, dofdata, dof_reference_node, constraints);
//
////	BC::constrainRefX_y(boundary, dofdata, dof_reference_node, constraints);
////	BC::constrainRefX_z(boundary, dofdata, dof_reference_node, constraints);
////	BC::constrainRefZ_y(boundary, dofdata, dof_reference_node, constraints);
//
//	BC::constrainRefX_x(boundary, dofdata, dof_reference_node, constraints);
//	BC::constrainRefY_y(boundary, dofdata, dof_reference_node, constraints);
//
//	BC::constrainMacroscopicallySymmetric(dofdata, dof_reference_node,
//			constraints);
//
//// Constrain faces in all directions
//	vector<bool> faces(dim, true);
//	vector<bool> directions(dim, true);
//	Periodic::constrainFaces(dofdata, dof_reference_node, constraints,
//			dof_handler, faces, directions);
//
//}

/*
 *  Apply constraints for periodic stressfree BC
 */
template<int dim>
void constrainStressfree(BoundaryValues<dim> *boundary,
		std::map<int, std::pair<Point<dim, double>, int> > &dofdata,
		vector<vector<int> > &dof_reference_node, ConstraintMatrix &constraints,
		hp::DoFHandler<dim> &dof_handler) {

	BC::constrainOrigin(boundary, dofdata, dof_reference_node, constraints);
	BC::constrainMacroscopicallySymmetric(dofdata, dof_reference_node,
			constraints);
//	BC::constrainRefX_yz(boundary, dofdata, dof_reference_node, constraints);
//	BC::constrainRefY_z(boundary, dofdata, dof_reference_node, constraints);

// Constrain faces in all directions
	vector<bool> faces(dim, true);
	vector<bool> directions(dim, true);
	Periodic::constrainFaces(dofdata, dof_reference_node, constraints,
			dof_handler, faces, directions);

}

/*
 *  Apply constraints for periodic stressfree BC
 */
template<int dim>
void constrainStressfreeCooling(BoundaryValues<dim> *boundary,
		std::map<int, std::pair<Point<dim, double>, int> > &dofdata,
		vector<vector<int> > &dof_reference_node, ConstraintMatrix &constraints,
		hp::DoFHandler<dim> &dof_handler) {

	BC::constrainOrigin(boundary, dofdata, dof_reference_node, constraints);
	BC::constrainMacroscopicallySymmetric(dofdata, dof_reference_node,
			constraints);
//	BC::constrainRefX_yz(boundary, dofdata, dof_reference_node, constraints);
//	BC::constrainRefY_z(boundary, dofdata, dof_reference_node, constraints);

// Constrain faces in all directions
	vector<bool> faces(dim, true);
	faces[2] = false;
	vector<bool> directions(dim, true);
	Periodic::constrainFaces(dofdata, dof_reference_node, constraints,
			dof_handler, faces, directions);

}

template<int dim>
void displacementShearPlaneStrain(hp::DoFHandler<dim> &dof_handler,
		ConstraintMatrix &periodic_constraints, BoundaryValues<dim> *boundary,
		int nslip, int maxCrystals) {

// encastre x-y-z on bottom
	std::vector<bool> components_DBC_XYZ(dim + nslip * maxCrystals, false);
	for (int i = 0; i < dim; i++) {
		components_DBC_XYZ[i] = true;
	}
	ComponentMask mask_DBC_XYZ(components_DBC_XYZ);
	VectorTools::interpolate_boundary_values(dof_handler, 2,
			ZeroFunction<dim>(dim + nslip * maxCrystals), periodic_constraints,
			mask_DBC_XYZ);

// encastre y-z on top
	std::vector<bool> components_DBC_YZ(dim + nslip * maxCrystals, false);
	components_DBC_YZ[1] = true;
	components_DBC_YZ[2] = true;
	ComponentMask mask_DBC_YZ(components_DBC_YZ);
	VectorTools::interpolate_boundary_values(dof_handler, 3,
			ZeroFunction<dim>(dim + nslip * maxCrystals), periodic_constraints,
			mask_DBC_YZ);

// shear x on top
	std::vector<bool> components_DBC_X(dim + nslip * maxCrystals, false);
	components_DBC_X[0] = true;
	ComponentMask mask_DBC_X(components_DBC_X);
	VectorTools::interpolate_boundary_values(dof_handler, 3,
			ConstantFunction<dim>(boundary->get_dx(),
					dim + nslip * maxCrystals), periodic_constraints,
			mask_DBC_X);

//plane strain
	std::vector<bool> components_DBC_plane_strain(dim + nslip * maxCrystals,
			false);
	components_DBC_plane_strain[2] = true;
	ComponentMask mask_DBC_slips(components_DBC_plane_strain);
	DoFTools::make_zero_boundary_constraints(dof_handler, 4,
			periodic_constraints, mask_DBC_slips);
	DoFTools::make_zero_boundary_constraints(dof_handler, 5,
			periodic_constraints, mask_DBC_slips);

//---------------- Periodic Constraints -----------------
// x direction, +x, -x surfaces
	DoFTools::make_periodicity_constraints<hp::DoFHandler<dim> >(dof_handler,
			2 * 0 + 1, 2 * 0, 0, periodic_constraints);
// x direction, +y, -y surfaces
//	DoFTools::make_periodicity_constraints<hp::DoFHandler<dim> >(dof_handler,
//			2 * 1 + 1, 2 * 1, 0, periodic_constraints);
}

template<int dim>
void displacementShearlikeGottschalk(hp::DoFHandler<dim> &dof_handler,
		ConstraintMatrix &periodic_constraints, BoundaryValues<dim> *boundary,
		int nslip, int maxCrystals) {

// encastre x-y-z on bottom
	std::vector<bool> components_DBC_XYZ(dim + nslip * maxCrystals, false);
	for (int i = 0; i < dim; i++) {
		components_DBC_XYZ[i] = true;
	}
	ComponentMask mask_DBC_XYZ(components_DBC_XYZ);
	VectorTools::interpolate_boundary_values(dof_handler, 2,
			ZeroFunction<dim>(dim + nslip * maxCrystals), periodic_constraints,
			mask_DBC_XYZ);

// shear x on top
	std::vector<bool> components_DBC_X(dim + nslip * maxCrystals, false);
	components_DBC_X[0] = true;
	ComponentMask mask_DBC_X(components_DBC_X);
	VectorTools::interpolate_boundary_values(dof_handler, 3,
			ConstantFunction<dim>(boundary->get_dx(),
					dim + nslip * maxCrystals), periodic_constraints,
			mask_DBC_X);

}

} // Namespace Periodic

// BC linear in plane and periodic in building direction
// CG: columnar grains
namespace Linear {

template<int dim>
void constrainFaces(std::map<int, std::pair<Point<dim, double>, int> > &dofdata,
		vector<vector<int> > &dof_reference_node, ConstraintMatrix &constraints,
		hp::DoFHandler<dim> &dof_handler, vector<bool> &faces,
		vector<bool> &directions) {

	// calculate size of RVE
	Point<dim> origin = dofdata[dof_reference_node[0][0]].first;
	vector<Point<dim>> referenceNodes;
	for (unsigned int i = 0; i < dim; ++i) {
		referenceNodes.push_back(
				Point<dim>(dofdata[dof_reference_node[i + 1][0]].first));
	}

	vector<double> sizeRVE(dim);
	for (unsigned int i = 0; i < dim; ++i) {
		sizeRVE[i] = referenceNodes[i][i] - origin[i];
	}
//	deallog << "Size RVE: " << sizeRVE[0] << " " << sizeRVE[1] << " " << sizeRVE[2] << endl;

	// make linear constraints based on reference nodes
	ConstraintMatrix tmp_constraints;
	tmp_constraints.clear();

	// component mask for displacement dofs
	unsigned int n_components = dof_handler.get_fe().n_components();
	std::vector<bool> u_dofs(n_components, false);
	u_dofs[0] = true;
	u_dofs[1] = true;
	u_dofs[2] = true;
	ComponentMask u_mask(u_dofs);

	// fill tmp_constraints with zero boundary conditions
	// and add later on the other constraints
	for (unsigned int i = 0; i < 2 * dim; ++i) {
		if (faces[i / 2]) {
			DoFTools::make_zero_boundary_constraints(dof_handler, i,
					tmp_constraints, u_mask);
		}
	}

	// loop over all constrained dofs and add the other constraints
	for (unsigned int i = 0; i < dof_handler.n_dofs(); i++) {
		// if dof is not constrained -> next dof
		if (!tmp_constraints.is_constrained(i)) {
			continue;
		}
		// reference nodes get there constraints in the main function if necessary
		if (isReference<dim>(i, dof_reference_node)) {
			continue;
		}

		// u_x = ( u_x(ref_y) - u_x(origin) ) / L_y * x_y
		if (dofdata.count(i) != 0) { // otherwise its no displacement dof
			Point<dim> position = dofdata[i].first;
			unsigned int component = dofdata[i].second; // =x
//			unsigned int dof_origin = dof_reference_node[0][component];

			if (!directions[component]) {
				continue;
			}

			if (!constraints.is_constrained(i)) {
				constraints.add_line(i);
			}
			for (unsigned y = 0; y < dim; ++y) {
				unsigned int dof_referenceNode =
						dof_reference_node[y + 1][component];
				double l_y = sizeRVE[y];
				double factor = position[y] / l_y;

				constraints.add_entry(i, dof_referenceNode, factor);
//				tmp_constraints.add_entry(i, dof_origin, -factor); // this is needed if origin is not fixed
			}
		} else {
			Assert(false, ExcInternalError());
		}
	}
}
template<int dim>
void constrain(BoundaryValues<dim> *boundary, ConstraintMatrix &constraints,
		hp::DoFHandler<dim> &dof_handler) {

	unsigned int n_components = dof_handler.get_fe().n_components();
	std::vector<bool> u_dofs(n_components, false);
	u_dofs[0] = true;
	u_dofs[1] = true;
	u_dofs[2] = true;

	WrapperBoundaryValues<dim> wrapperBoundaryValues(boundary, n_components);
	ComponentMask u_mask(u_dofs);

	ConstraintMatrix tmp_constraints;
	tmp_constraints.clear();

	for (unsigned int i = 0; i < 2 * dim; ++i) {
		VectorTools::interpolate_boundary_values(dof_handler, i,
				wrapperBoundaryValues, tmp_constraints, u_mask);
	}

	tmp_constraints.close();
	constraints.merge(tmp_constraints, ConstraintMatrix::left_object_wins);

}

/*
 *  Apply constraints for linear uniaxial tension BC
 */
template<int dim>
void constrainUniaxialTension(BoundaryValues<dim> *boundary,
		std::map<int, std::pair<Point<dim, double>, int> > &dofdata,
		vector<vector<int> > &dof_reference_node, ConstraintMatrix &constraints,
		hp::DoFHandler<dim> &dof_handler, string load) {

	BC::constrainOrigin(boundary, dofdata, dof_reference_node, constraints);

	if (load.compare("11") == 0) {
		BC::constrainRefX_x(boundary, dofdata, dof_reference_node, constraints);
	} else if (load.compare("22") == 0) {
		BC::constrainRefY_y(boundary, dofdata, dof_reference_node, constraints);
	} else if (load.compare("33") == 0) {
		BC::constrainRefZ_z(boundary, dofdata, dof_reference_node, constraints);
	} else if (load.compare("12") == 0) {
		BC::constrainRefX_y(boundary, dofdata, dof_reference_node, constraints);
	} else if (load.compare("13") == 0) {
		BC::constrainRefX_z(boundary, dofdata, dof_reference_node, constraints);
	} else if (load.compare("23") == 0) {
		BC::constrainRefY_z(boundary, dofdata, dof_reference_node, constraints);
	} else {
		AssertThrow(false, ExcMessage("Wrong straindir."));
	}

	BC::constrainMacroscopicallySymmetric(dofdata, dof_reference_node,
			constraints);

	// Constrain faces in all directions
	vector<bool> directions_linear(dim, true);
	vector<bool> faces_linear(dim, true);
	Linear::constrainFaces(dofdata, dof_reference_node, constraints,
			dof_handler, faces_linear, directions_linear);

}

/*
 *  Apply constraints for linear plane stress
 */
template<int dim>
void constrainPlaneStress(BoundaryValues<dim> *boundary,
		std::map<int, std::pair<Point<dim, double>, int> > &dofdata,
		vector<vector<int> > &dof_reference_node, ConstraintMatrix &constraints,
		hp::DoFHandler<dim> &dof_handler, vector<string> &loads) {

	BC::constrainOrigin(boundary, dofdata, dof_reference_node, constraints);

//	BC::constrainRefX_x(boundary, dofdata, dof_reference_node, constraints);
//	BC::constrainRefZ_z(boundary, dofdata, dof_reference_node, constraints);

	for (unsigned int i = 0; i < loads.size(); ++i) {
		if (loads[i].compare("11") == 0) {
			BC::constrainRefX_x(boundary, dofdata, dof_reference_node,
					constraints);
		} else if (loads[i].compare("22") == 0) {
			BC::constrainRefY_y(boundary, dofdata, dof_reference_node,
					constraints);
		} else if (loads[i].compare("33") == 0) {
			BC::constrainRefZ_z(boundary, dofdata, dof_reference_node,
					constraints);
		} else if (loads[i].compare("12") == 0) {
			BC::constrainRefX_y(boundary, dofdata, dof_reference_node,
					constraints);
		} else if (loads[i].compare("13") == 0) {
			BC::constrainRefX_z(boundary, dofdata, dof_reference_node,
					constraints);
		} else if (loads[i].compare("23") == 0) {
			BC::constrainRefY_z(boundary, dofdata, dof_reference_node,
					constraints);
		} else {
			AssertThrow(false, ExcMessage("Wrong straindir."));
		}
	}

	BC::constrainMacroscopicallySymmetric(dofdata, dof_reference_node,
			constraints);

// Constrain faces in all directions
	vector<bool> directions_linear(dim, true);
	vector<bool> faces_linear(dim, true);
	Linear::constrainFaces(dofdata, dof_reference_node, constraints,
			dof_handler, faces_linear, directions_linear);

}

///*
// *  Apply constraints for linear plane stress in XY BC
// */
//template<int dim>
//void constrainPlaneStressXY(BoundaryValues<dim> *boundary,
//		std::map<int, std::pair<Point<dim, double>, int> > &dofdata,
//		vector<vector<int> > &dof_reference_node, ConstraintMatrix &constraints,
//		hp::DoFHandler<dim> &dof_handler) {
//
//	BC::constrainOrigin(boundary, dofdata, dof_reference_node, constraints);
//
//	BC::constrainRefX_x(boundary, dofdata, dof_reference_node, constraints);
//	BC::constrainRefY_y(boundary, dofdata, dof_reference_node, constraints);
//
//	BC::constrainMacroscopicallySymmetric(dofdata, dof_reference_node,
//			constraints);
//
//// Constrain faces in all directions
//	vector<bool> directions_linear(dim, true);
//	vector<bool> faces_linear(dim, true);
//	Linear::constrainFaces(dofdata, dof_reference_node, constraints,
//			dof_handler, faces_linear, directions_linear);
//
//}

namespace CG {
/*
 *  Apply linear constraints in isotropic plane and periodic in BD
 */
template<int dim>
void constrain(BoundaryValues<dim> *boundary,
		std::map<int, std::pair<Point<dim, double>, int> > &dofdata,
		vector<vector<int> > &dof_reference_node, ConstraintMatrix &constraints,
		hp::DoFHandler<dim> &dof_handler) {

	BC::constrainOrigin(boundary, dofdata, dof_reference_node, constraints);
	BC::constrainRefX(boundary, dofdata, dof_reference_node, constraints);
	BC::constrainRefY(boundary, dofdata, dof_reference_node, constraints);
	BC::constrainRefZ(boundary, dofdata, dof_reference_node, constraints);

	// Constrain faces in all directions
	vector<bool> directions_linear(dim, true);
	vector<bool> faces_linear(dim, true);
	faces_linear[dim - 1] = false;
	Linear::constrainFaces(dofdata, dof_reference_node, constraints,
			dof_handler, faces_linear, directions_linear);

	vector<bool> directions_periodic(dim, true);
//	directions_periodic[dim] = true;
	vector<bool> faces(dim, false);
	faces[dim - 1] = true;
	Periodic::constrainFaces(dofdata, dof_reference_node, constraints,
			dof_handler, faces, directions_periodic);

}

/*
 *  Apply constraints for linear uniaxial tension BC with periodicity in BD
 */
template<int dim>
void constrainUniaxialTension(BoundaryValues<dim> *boundary,
		std::map<int, std::pair<Point<dim, double>, int> > &dofdata,
		vector<vector<int> > &dof_reference_node, ConstraintMatrix &constraints,
		hp::DoFHandler<dim> &dof_handler, string load) {

	BC::constrainOrigin(boundary, dofdata, dof_reference_node, constraints);

	if (load.compare("11") == 0) {
		BC::constrainRefX_x(boundary, dofdata, dof_reference_node, constraints);
	} else if (load.compare("22") == 0) {
		BC::constrainRefY_y(boundary, dofdata, dof_reference_node, constraints);
	} else if (load.compare("33") == 0) {
		BC::constrainRefZ_z(boundary, dofdata, dof_reference_node, constraints);
	} else if (load.compare("12") == 0) {
		BC::constrainRefX_y(boundary, dofdata, dof_reference_node, constraints);
	} else if (load.compare("13") == 0) {
		BC::constrainRefX_z(boundary, dofdata, dof_reference_node, constraints);
	} else if (load.compare("23") == 0) {
		BC::constrainRefY_z(boundary, dofdata, dof_reference_node, constraints);
	} else {
		AssertThrow(false, ExcMessage("Wrong straindir."));
	}

	BC::constrainMacroscopicallySymmetric(dofdata, dof_reference_node,
			constraints);

	// Constrain faces in all directions
	vector<bool> directions_linear(dim, true);
	vector<bool> faces_linear(dim, true);
	faces_linear[2] = false;
	Linear::constrainFaces(dofdata, dof_reference_node, constraints,
			dof_handler, faces_linear, directions_linear);

	vector<bool> directions_periodic(dim, true);
	vector<bool> faces_periodic(dim, false);
	faces_periodic[2] = true;

	Periodic::constrainFaces(dofdata, dof_reference_node, constraints,
			dof_handler, faces_periodic, directions_periodic);

}

/*
 *  Apply constraints for linear plane stress BC with periodicity in BD
 */
template<int dim>
void constrainPlaneStress(BoundaryValues<dim> *boundary,
		std::map<int, std::pair<Point<dim, double>, int> > &dofdata,
		vector<vector<int> > &dof_reference_node, ConstraintMatrix &constraints,
		hp::DoFHandler<dim> &dof_handler, vector<string> &loads) {

	BC::constrainOrigin(boundary, dofdata, dof_reference_node, constraints);

//	BC::constrainRefX_x(boundary, dofdata, dof_reference_node, constraints);
//	BC::constrainRefZ_z(boundary, dofdata, dof_reference_node, constraints);

	for (unsigned int i = 0; i < loads.size(); ++i) {
		if (loads[i].compare("11") == 0) {
			BC::constrainRefX_x(boundary, dofdata, dof_reference_node,
					constraints);
		} else if (loads[i].compare("22") == 0) {
			BC::constrainRefY_y(boundary, dofdata, dof_reference_node,
					constraints);
		} else if (loads[i].compare("33") == 0) {
			BC::constrainRefZ_z(boundary, dofdata, dof_reference_node,
					constraints);
		} else if (loads[i].compare("12") == 0) {
			BC::constrainRefX_y(boundary, dofdata, dof_reference_node,
					constraints);
		} else if (loads[i].compare("13") == 0) {
			BC::constrainRefX_z(boundary, dofdata, dof_reference_node,
					constraints);
		} else if (loads[i].compare("23") == 0) {
			BC::constrainRefY_z(boundary, dofdata, dof_reference_node,
					constraints);
		} else {
			AssertThrow(false, ExcMessage("Wrong straindir."));
		}
	}

	BC::constrainMacroscopicallySymmetric(dofdata, dof_reference_node,
			constraints);

	// Constrain faces in all directions
	vector<bool> directions_linear(dim, true);
	vector<bool> faces_linear(dim, true);
	faces_linear[2] = false;
	Linear::constrainFaces(dofdata, dof_reference_node, constraints,
			dof_handler, faces_linear, directions_linear);

	vector<bool> directions_periodic(dim, true);
	vector<bool> faces_periodic(dim, false);
	faces_periodic[2] = true;

	Periodic::constrainFaces(dofdata, dof_reference_node, constraints,
			dof_handler, faces_periodic, directions_periodic);

}

///*
// *  Apply constraints for linear plane stress in XY BC with periodicity in BD
// */
//template<int dim>
//void constrainPlaneStressXY(BoundaryValues<dim> *boundary,
//		std::map<int, std::pair<Point<dim, double>, int> > &dofdata,
//		vector<vector<int> > &dof_reference_node, ConstraintMatrix &constraints,
//		hp::DoFHandler<dim> &dof_handler) {
//
//	BC::constrainOrigin(boundary, dofdata, dof_reference_node, constraints);
//
//	BC::constrainRefX_x(boundary, dofdata, dof_reference_node, constraints);
//	BC::constrainRefY_y(boundary, dofdata, dof_reference_node, constraints);
//
//	BC::constrainMacroscopicallySymmetric(dofdata, dof_reference_node,
//			constraints);
//
//	// Constrain faces in all directions
//	vector<bool> directions_linear(dim, true);
//	vector<bool> faces_linear(dim, true);
//	faces_linear[2] = false;
//	Linear::constrainFaces(dofdata, dof_reference_node, constraints,
//			dof_handler, faces_linear, directions_linear);
//
//	vector<bool> directions_periodic(dim, true);
//	vector<bool> faces_periodic(dim, false);
//	faces_periodic[2] = true;
//
//	Periodic::constrainFaces(dofdata, dof_reference_node, constraints,
//			dof_handler, faces_periodic, directions_periodic);
//
//}

}// end Namespace CG

} // end Namespace Linear

namespace TestcaseShear {

template<int dim>
void microHardGottschalk(hp::DoFHandler<dim> &dof_handler,
		ConstraintMatrix &constraints, BoundaryValues<dim> *boundary, int nslip,
		int maxCrystals) {

	ConstraintMatrix periodic_constraints;
	periodic_constraints.clear();

	BC::microhardTopBottom(dof_handler, periodic_constraints, nslip,
			maxCrystals);
	BC::Periodic::displacementShearlikeGottschalk(dof_handler,
			periodic_constraints, boundary, nslip, maxCrystals);

	periodic_constraints.close();

	constraints.merge(periodic_constraints, ConstraintMatrix::left_object_wins);

}

template<int dim>
void microHard(hp::DoFHandler<dim> &dof_handler, ConstraintMatrix &constraints,
		BoundaryValues<dim> *boundary, int nslip, int maxCrystals) {

	ConstraintMatrix periodic_constraints;
	periodic_constraints.clear();

	BC::microhardTopBottom(dof_handler, periodic_constraints, nslip,
			maxCrystals);
	BC::Periodic::displacementShearPlaneStrain(dof_handler,
			periodic_constraints, boundary, nslip, maxCrystals);

	periodic_constraints.close();

	constraints.merge(periodic_constraints, ConstraintMatrix::left_object_wins);

}

template<int dim>
void microFree(hp::DoFHandler<dim> &dof_handler, ConstraintMatrix &constraints,
		BoundaryValues<dim> *boundary, int nslip, int maxCrystals) {

	ConstraintMatrix periodic_constraints;
	periodic_constraints.clear();

	BC::Periodic::displacementShearPlaneStrain(dof_handler,
			periodic_constraints, boundary, nslip, maxCrystals);

	periodic_constraints.close();

	constraints.merge(periodic_constraints, ConstraintMatrix::left_object_wins);

}

} // Namespace TestcaseShear

} // Namespace BC

#endif /* CODE_VARIOUS_BOUNDARYCONDITIONS_H_ */
