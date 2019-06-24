/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2000 - 2015 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the deal.II distribution.
 *
 * ---------------------------------------------------------------------

 *
 * Author: Wolfgang Bangerth, University of Heidelberg, 2000
 */

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/work_stream.h>
#include <deal.II/base/parameter_handler.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/constraint_matrix.h>

#include <deal.II/lac/sparse_direct.h>
#include <deal.II/dofs/dof_renumbering.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_nothing.h>

#include <deal.II/hp/dof_handler.h>
#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/fe_values.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>

#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include "Code/Material/CrystalPlasticityCreep.h"
#include "Code/Material/localCP.h"
#include "Code/Various/parameter.h"
#include "Code/Postprocessing/Postprocessing.h"
#include "Code/Various/Tensor8.h"
#include "Code/Various/Right_Hand_Side.h"
#include "Code/Various/Boundary/BoundaryConditions.h"
#include "Code/Various/AuxilaryFunctions.h"
#include "Code/Various/Boundary/Boundary_Values.h"
#include "Code/Various/Mesh/meshs.h"

namespace SFB814_C5 {
using namespace dealii;

inline void setConstrainedDofsToZero(Vector<double> &newton_update,
                                     ConstraintMatrix &constraints) {
    for (unsigned int i = 0; i < newton_update.size(); i++) {
        if (constraints.is_inhomogeneously_constrained(i)) {
            newton_update(i) = 0;
            deallog << i << " ";
        }
    }
    deallog << endl;
}

template<int dim>
inline bool isReference(int dof, vector<vector<int> > &dofsRef) {
    Assert(dofsRef.size() == dim+1, ExcInternalError());
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

template<int dim>
struct PerTaskData {
    FullMatrix<double> cell_matrix;
    Vector<double> cell_rhs;
    std::vector<unsigned int> local_dof_indices;

    vector<FullMatrix<double>> cell_coupling_matrices;
    vector<vector<unsigned int>> local_dof_indices_neighbours;
    vector<bool> is_grain_boundary;

    //	PerTaskData(const FESystem<dim> &fe) :
    //			cell_matrix(fe.dofs_per_cell, fe.dofs_per_cell), cell_rhs(
    //					fe.dofs_per_cell), local_dof_indices(fe.dofs_per_cell) {
    //	}
    PerTaskData(const hp::FECollection<dim> fe_collection) :
        cell_matrix(fe_collection.max_dofs_per_cell(),
                    fe_collection.max_dofs_per_cell()), cell_rhs(
                        fe_collection.max_dofs_per_cell()), local_dof_indices(
                            fe_collection.max_dofs_per_cell()), cell_coupling_matrices(
                                GeometryInfo<dim>::faces_per_cell,
                                FullMatrix<double>(fe_collection.max_dofs_per_cell(),
                                        fe_collection.max_dofs_per_cell())), local_dof_indices_neighbours(
                                            GeometryInfo<dim>::faces_per_cell,
                                            vector<unsigned int>(fe_collection.max_dofs_per_cell())), is_grain_boundary(
                                                GeometryInfo<dim>::faces_per_cell) {
    }
};

template<int dim>
struct ScratchData {
    hp::FEValues<dim> hp_fe_values;
    hp::FEFaceValues<dim> hp_fe_face_values;
    hp::FEFaceValues<dim> hp_fe_face_values_neighbor;
    ScratchData(const hp::FECollection<dim> &fe_collection,
                const hp::QCollection<dim> &quadrature_collection_cell,
                const UpdateFlags update_flags_cell,
                const hp::QCollection<dim - 1> &quadrature_collection_interface,
                const UpdateFlags update_flags_interface) :
        hp_fe_values(fe_collection, quadrature_collection_cell,
                     update_flags_cell), hp_fe_face_values(fe_collection,
                             quadrature_collection_interface, update_flags_interface), hp_fe_face_values_neighbor(
                                 fe_collection, quadrature_collection_interface,
                                 update_flags_interface) {
    }
    ScratchData(const ScratchData &scratch) :
        hp_fe_values(scratch.hp_fe_values.get_fe_collection(),
                     scratch.hp_fe_values.get_quadrature_collection(),
                     scratch.hp_fe_values.get_update_flags()), hp_fe_face_values(
                         scratch.hp_fe_face_values.get_fe_collection(),
                         scratch.hp_fe_face_values.get_quadrature_collection(),
                         scratch.hp_fe_face_values.get_update_flags()), hp_fe_face_values_neighbor(
                             scratch.hp_fe_face_values_neighbor.get_fe_collection(),
                             scratch.hp_fe_face_values_neighbor.get_quadrature_collection(),
                             scratch.hp_fe_face_values_neighbor.get_update_flags()) {
    }
};

template<int dim, int nslip>
class ElastoplasticProblem {
    friend class BoundaryValues<dim> ;
public:
    ElastoplasticProblem(/*int feq_u, int feq_gamma*/int number_gp);
    ~ElastoplasticProblem();
    void run();

private:
    // stuff for setup system
    void createMesh();
    void setup_system();
    void setup_sparsity_pattern();
    //	void update_quadrature_point_history();
    //	void quadrature_point_history_copy_newtonstep_to_history();
    void setup_quadrature_point_history();
    void identify_grainboundarys();

    // stuff for assembling the linear system
    void assemble_system();
    void apply_constraints();
    void sparseindices_on_one_cell(
        const typename hp::DoFHandler<dim>::active_cell_iterator &cell,
        PerTaskData<dim> &data);
    void assemble_on_one_cell(
        const typename hp::DoFHandler<dim>::active_cell_iterator &cell,
        ScratchData<dim> &scratch, PerTaskData<dim> &data);
    void copy_local_to_global(const PerTaskData<dim> &data);
    void update_hardening();

    // other calculation modes
    void calculateYieldSurface();
    void calculateYieldPoint();

    // ----------------------------------------
    // stuff for the solver
    // ----------------------------------------
    void solve();
    Vector<double> mechanical_prediction_get_delta(
        vector<Vector<double>> &solution_history,
        vector<double> &times_history);
    vector<double> calculateResiduals(Vector<double> &res);
    void slotConditionNumber(double d);

    // ----------------------------------------
    // Some postprocessing
    // ----------------------------------------
    void postprocessing(unsigned int loadstep);
    void calculateMacroscopicTangent();
    void saveSystemMatrix(int loadstep/*Vector<double> &dSolution*/);
    SymmetricTensor<2, dim> getMacroscopicStress();
    SymmetricTensor<2, dim> getMacroscopicStrain();
    SymmetricTensor<2, dim> getMacroscopicPlasticStrain();
    SymmetricTensor<2, dim> getMacroscopicPlasticStrainRate();
    double calculate_dissipation();

    int newtonCorrection(Vector<double> &newtonStep);
    //	double maxValueEvolutionEquation();
    //	int count_changed_active_systems();

    // ----------------------------------------
    // Some output for the console
    // ----------------------------------------
    void print_conv_footer();
    void print_conv_header();
    void print_detailed_crystal_informations();

    static const Tensor<8, dim> P_transversal_isotropic;
    int maxCrystals;

    TimerOutput timer;
    ParameterHandler parameter;

    Triangulation<dim> triangulation;

    hp::DoFHandler<dim> dof_handler;
    hp::FECollection<dim> fe_collection;
    hp::QCollection<dim> quadrature_collection;
    hp::QCollection<dim - 1> face_quadrature_collection;

    //	DoFHandler<dim> dof_handler;
    //	FESystem<dim> fe;

    SparsityPattern sparsity_pattern;
    SparseMatrix<double> system_matrix;

    Vector<double> current_solution;
    Vector<double> prior_solution; // solution of the previous loadstep
    Vector<double> current_solution_abs; // absolute solutions for hardening
    Vector<double> prior_solution_abs; // absolute solutions for hardening
    Vector<double> current_hardening; // current hardening in slipsystems
    Vector<double> prior_hardening; // prior hardening in slipsystems
    Vector<double> system_rhs;

    //	QGauss<dim> quadrature_formula;

    //	std::vector<Material<dim, nslip> > quadrature_point_history;
    std::vector<PerCrystalData<dim, nslip> > crystal_data; // TODO use a Map here, Problem is that this
    //implementation only works if MaterialIDs are ongoing

    //std::vector<int> domainIDs;
    vector<bool> dof_is_grainboundary;

    ConstraintMatrix constraints;

    FEValuesExtractors::Vector extractor_displacement;
    vector<vector<FEValuesExtractors::Scalar>> extractor_slips;
    BoundaryValues<dim>* boundary;
    vector<double> geometric_range; //vector<double(3), assume that 0/0/0 is one edge, than the x/y and z limits are stored

    // dof_reference_node[i][j] = dof of reference node i in direction j
    // i=0 -> origin, i=1 -> refX, ...
    vector<vector<int> > dof_reference_node; // list of dof of the reference nodes for displacement bc

    vector<int> compToCPSystem;

    // map for coordinates, dofs, and comp(dof)
    //coords[key]= <dof, <Point(x,y,z), comp(dof)> >
    //only displacement-dofs are stored
    std::map<int, std::pair<Point<dim, double>, int> > dofdata;   // relate key----value
    double conditionNumber = 0;

    bool reCreateSparsityPattern;

    int n_iter_remove_me;
    Postprocessing<dim, nslip> post;

};

template<int dim, int nslip>
ElastoplasticProblem<dim, nslip>::ElastoplasticProblem(
    /*int feq_u, int feq_gamma,*/int number_gp) :
    maxCrystals(-1), // number of crystals will be stored here
    timer(deallog.get_file_stream(), TimerOutput::summary,
          TimerOutput::wall_times), // set options for timer
    dof_handler(triangulation), //
    quadrature_collection(QGauss<dim>(number_gp)), //QGauss<dim>(n) = Generate a formula with n quadrature points (in each space direction), exact for polynomials of degree 2n-1. 
    face_quadrature_collection(QGauss<dim - 1>(number_gp)), //
    sparsity_pattern(), //
    extractor_displacement(0), //
    reCreateSparsityPattern(true), //
    n_iter_remove_me(-1), //
    post(dof_handler, quadrature_collection, fe_collection, crystal_data,
         extractor_displacement, extractor_slips) {

    // set standard values for geomatric range (size of rve)
    geometric_range.resize(dim);
    for (int i = 0; i < dim; i++) {
        geometric_range[i] = 1;
    }

    // read in parameterfile
    read_Parameterfile("Input/parameter", parameter);

    // decide on which loadtype is needed and set appropriate Boundary object
    if (parameter.get("strainpath").compare("prescribed_from_file") == 0) {
        boundary = new BoundaryValuesPrescribed<dim>(&parameter);
    } else if (parameter.get("strainpath").compare("linear_displacement")
               == 0) {
        boundary = new BoundaryValuesLinear<dim>(&parameter);
    }

    // set limit for shared memory parallelization
    unsigned int n_cores = parameter.get_integer("nCoresAssembly");
    MultithreadInfo::set_thread_limit(n_cores);
}

template<int dim, int nslip>
const Tensor<8, dim> ElastoplasticProblem<dim, nslip>::P_transversal_isotropic =
    calculate_transversal_isotropic_projector<dim>();

template<int dim, int nslip>
ElastoplasticProblem<dim, nslip>::~ElastoplasticProblem() {
    dof_handler.clear();
}

template<int dim, int nslip>
void ElastoplasticProblem<dim, nslip>::setup_system() {
    timer.enter_section("Setup System");

    /*
     * initialize extractor_slip
     * first dimension (nc) stands for the fe_element with
     * dim x fe_u, nc x nslip x fe_nothing, nslip x fe_g, (cmax - 1 - nc) x nslip x fe_nothing
     */
    for (int c = 0; c < maxCrystals; c++) {
        vector<FEValuesExtractors::Scalar> extractor_slips_cn; 
        for (int i = 0; i < nslip; i++) {
            extractor_slips_cn.push_back(
                FEValuesExtractors::Scalar(dim + c * nslip + i)); //It's an index corresponding to dim + c*nslip +i, can be used as
        }                                                          //fe_values[extractor_slips_cn[n]] 
        extractor_slips.push_back(extractor_slips_cn);
    }

    /*
     * initialize compToCPSystem, a vector with the information which
     * local dof [0, nslip-1] a component of the FESystem belongs to
     */
    compToCPSystem.resize(dim + nslip * maxCrystals);
    int comp_vec = 0;
    // first [0; dim-1]
    for (; comp_vec < dim; comp_vec++) {
        compToCPSystem[comp_vec] = comp_vec;
    }
    // then the slipsystems [dim; dim + nslip -1]
    for (int i = 0; i < maxCrystals; i++) {
        for (int k = 0; k < nslip; k++) {
            compToCPSystem[comp_vec] = k + dim;    //compToCPSystem = [0,1,2,[3,4]*2_Crystals]= [0,1,2,3,4,3,4]
            comp_vec++;                                              
        }
    }

    // determine the degree of ansatzfunctions for displacments and slips
    int feq_u = 0;
    if (parameter.get("el_type_u").compare("Q1") == 0) {
        feq_u = 1;
    } else if (parameter.get("el_type_u").compare("Q2") == 0) {
        feq_u = 2;
    } else if (parameter.get("el_type_u").compare("Q3") == 0) {
        feq_u = 3;
    }
    int feq_gamma = 0;
    if (parameter.get("el_type_gamma").compare("Q1") == 0) {
        feq_gamma = 1;
    } else if (parameter.get("el_type_gamma").compare("Q2") == 0) {
        feq_gamma = 2;
    } else if (parameter.get("el_type_gamma").compare("Q3") == 0) {
        feq_gamma = 3;
    }

    /*
     * generate the fe_collection
     * for each crystal one entry with
     * [u1, u2, u3, fe_nothing..., slip1, slip2, ... slip_nslip, fe_nothing...]
     */
    for (int i = 0; i < maxCrystals; i++) {
        vector<const FiniteElement<dim>*> fe_list;
        for (int k = 0; k < dim; k++) {
            fe_list.push_back(new FE_Q<dim>(feq_u));       
        }
        for (int j = 0; j < maxCrystals; j++) {
            for (int k = 0; k < nslip; k++) {
                if (i == j) {
                    fe_list.push_back(new FE_Q<dim>(feq_gamma));   
                } else {
                    fe_list.push_back(new FE_Nothing<dim>());        
                }
            }
        }         //fe_list = [F1ux,F1uy,F1uz,F1ge1,F1ge2,0,0]
        fe_collection.push_back(
            FESystem<dim>(fe_list,
                          vector<unsigned int>(fe_list.size(), 1))); 
        for (unsigned int j = 0; j < fe_list.size(); j++) {
            delete (fe_list[j]);
        }
    }
	//fe_collection = [[F1ux,F1uy,F1uz,F1ge1,F1ge2,0,0],[F2ux,F2uy,F2uz,0,0,F2g1,F2g2]]
    dof_handler.distribute_dofs(fe_collection);          //renumbering all vertices and quadrature points

    //	DoFRenumbering::Cuthill_McKee (dof_handler); //condition number -30percent
    DoFRenumbering::boost::king_ordering(dof_handler); //condition number -40percent
    //	DoFRenumbering::boost::minimum_degree (dof_handler); //condition number -10percent

    deallog << "Ndof: " << dof_handler.n_dofs() << std::endl;
    deallog << "Nel: " << triangulation.n_active_cells() << std::endl;

    current_solution.reinit(dof_handler.n_dofs());
    prior_solution.reinit(dof_handler.n_dofs());
    current_solution_abs.reinit(dof_handler.n_dofs());
    prior_solution_abs.reinit(dof_handler.n_dofs());
    current_hardening.reinit(dof_handler.n_dofs());
    prior_hardening.reinit(dof_handler.n_dofs());
    system_rhs.reinit(dof_handler.n_dofs());

    dof_is_grainboundary.resize(dof_handler.n_dofs());
    dof_reference_node.resize(dim + 1, vector<int>(dim));      //dof_reference_node = [[v0,v1,v2],[],[],[]]

    // identify Reference Nodes and fill list of dofData
    // Coordinates of reference points
    Point<dim> origin(0, 0, 0); // origin
    Point<dim> ref_x(geometric_range[0], 0, 0); // reference x
    Point<dim> ref_y(0, geometric_range[1], 0); // reference y
    Point<dim> ref_z(0, 0, geometric_range[2]); // reference z     

    Quadrature<dim> q(fe_collection[0].get_unit_support_points()); //In Unit cell,support point p_j = where shape function v_i(p_j) = delta_ij
    //Construct a dummy quadrature formula from a list of points, with weights set to infinity
	hp::QCollection<dim> q_collection(q);
    hp::FEValues<dim> hp_fe_values(fe_collection, q_collection,    // collection of fe values at vertex and quadrature points
                                   update_quadrature_points);

    typename hp::DoFHandler<dim>::active_cell_iterator cell =
        dof_handler.begin_active(), endc = dof_handler.end();
    for (; cell != endc; ++cell) {
        const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;
        vector<types::global_dof_index> local_dof_indices(dofs_per_cell);    // global_dof_index = an unsigned integer type.
        // local_dof_indices = [x1, x2, ... x_cellDof]
        //---------------------------
        hp_fe_values.reinit(cell);            //reinitialize values and gradients
        const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values();
        const vector<Point<dim> > node_coords =
            fe_values.get_quadrature_points();
        cell->get_dof_indices(local_dof_indices);    //write the global indices of the DoF located on this cell in the standard ordering defined by the finite element to vector local_dof_indices
        for (unsigned int i = 0; i < dofs_per_cell; i++) {
            int dof = local_dof_indices[i];
            const unsigned int component_dof =
                compToCPSystem[fe_collection[cell->active_fe_index()].system_to_component_index(   //active_fe_index return the index inside the hp::FECollection of the FiniteElement used for this cell
                                   i).first]; // i.e. fe_collection[0].s_t_c_i(0).first = 0 ---> compToCPSystem[0]=0
            if (component_dof < dim) {
                // map for coordinates, dofs, and comp(dof)
                //coords[key]= <dof, <Point(x,y,z), comp(dof)> >
                //only displacement-dofs are stored
                dofdata.insert(             //<dof, <Point(x,y,z), comp(dof)> >
                    std::pair<int, std::pair<Point<dim, double>, int> >(dof, //pair <int, char> PAIR1, PAIR1.first = 100, PAIR1.second=G
                            std::pair<Point<dim, double>, int>(         //
                                node_coords[i], component_dof)));
                if ((origin - node_coords[i]).norm() < 1e-5) { //check distance
                    dof_reference_node[0][component_dof] = dof;
//										deallog << "Found origin" << endl;
                }
                if ((ref_x - node_coords[i]).norm() < 1e-5) { //check distance
                    dof_reference_node[1][component_dof] = dof;
//										deallog << "Found reference node x " << endl;
                }
                if ((ref_y - node_coords[i]).norm() < 1e-5) { //check distance
                    dof_reference_node[2][component_dof] = dof;
//										deallog << "Found reference node y " << endl;
                }
                if ((ref_z - node_coords[i]).norm() < 1e-5) { //check distance
                    dof_reference_node[3][component_dof] = dof;
//										deallog << "Found reference node z " << endl;
                }
                //				deallog << "dof: " << dof << " coords: " << node_coords[i] << endl;
            }
        }
    }
    timer.exit_section("Setup System");
}

template<int dim, int nslip>
void ElastoplasticProblem<dim, nslip>::setup_sparsity_pattern() {
    timer.enter_section("Make sparsity pattern");
    if (!reCreateSparsityPattern) {
        //		system_matrix.reinit(sparsity_pattern);
        system_matrix = 0;
        timer.exit_section("Make sparsity pattern");
        return;
    }
    reCreateSparsityPattern = false;
    DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
    typename hp::DoFHandler<dim>::active_cell_iterator cell =
        dof_handler.begin_active(), end_cell = dof_handler.end();

    for (; cell != end_cell; ++cell) {
        PerTaskData<dim> data(fe_collection);

        sparseindices_on_one_cell(cell, data);

        constraints.add_entries_local_to_global(data.local_dof_indices, dsp, /*keep constrained dofs*/
                                                false);

        if (parameter.get("GB_behaviour").compare("detailed") == 0) {
            for (unsigned int i = 0; i < data.is_grain_boundary.size(); i++) {
                if (data.is_grain_boundary[i]) {
                    constraints.add_entries_local_to_global(
                        data.local_dof_indices,
                        data.local_dof_indices_neighbours[i], dsp, /*keep constrained dofs*/
                        false);
                }
            }
        }
    }
    dsp.compress();
    sparsity_pattern.copy_from(dsp);
    sparsity_pattern.compress();
    system_matrix.reinit(sparsity_pattern);
    timer.exit_section("Make sparsity pattern");
}

template<int dim, int nslip>
void ElastoplasticProblem<dim, nslip>::assemble_on_one_cell(
    const typename hp::DoFHandler<dim>::active_cell_iterator &cell,
    ScratchData<dim> &scratch, PerTaskData<dim> &data) {

    const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    scratch.hp_fe_values.reinit(cell);
    const FEValues<dim> &fe_values =
        scratch.hp_fe_values.get_present_fe_values();
    const unsigned int n_q_points = fe_values.n_quadrature_points;
    if (data.cell_rhs.size() == dofs_per_cell) {
        data.cell_matrix = 0;                // initalize all components to 0
        data.cell_rhs = 0;
    } else {
        data.cell_matrix.reinit(dofs_per_cell, dofs_per_cell);
        data.cell_rhs.reinit(dofs_per_cell);
    }

    vector<SymmetricTensor<2, dim> > strains(n_q_points);        //
    fe_values[extractor_displacement].get_function_symmetric_gradients(   // return values within index range [extractor_displacement] = [a:b] of current_solution to strains
        current_solution, strains);

    vector<vector<double> > slips(nslip, vector<double>(n_q_points));   // slips = [[vector with size n_q_points]*nslip]
    vector<vector<double> > prior_slips(nslip, vector<double>(n_q_points));
    vector<vector<double> > prior_slips_abs(nslip, vector<double>(n_q_points));
    vector<vector<double> > current_hardening_cell(nslip,
            vector<double>(n_q_points));
    vector<vector<double> > current_gamma_abs_cell(nslip,
            vector<double>(n_q_points));
    vector<vector<Tensor<1, dim> > > slip_grads(nslip,              //slip_grads = [[Rank_1_Tensor*n_q_points]*nslip]
            vector<Tensor<1, dim> >(n_q_points));
    for (unsigned int i = 0;
            i < extractor_slips[cell->active_fe_index()].size(); i++) {
        fe_values[extractor_slips[cell->active_fe_index()][i]].get_function_values(    // [in]: current_solution, [out]: slips[i]
            current_solution, slips[i]);
        fe_values[extractor_slips[cell->active_fe_index()][i]].get_function_values(
            prior_solution, prior_slips[i]);
        fe_values[extractor_slips[cell->active_fe_index()][i]].get_function_values(
            prior_solution_abs, prior_slips_abs[i]);

        fe_values[extractor_slips[cell->active_fe_index()][i]].get_function_values(
            current_hardening, current_hardening_cell[i]);
        fe_values[extractor_slips[cell->active_fe_index()][i]].get_function_values(
            current_solution_abs, current_gamma_abs_cell[i]);

        fe_values[extractor_slips[cell->active_fe_index()][i]].get_function_gradients(
            current_solution, slip_grads[i]);
    }
    data.local_dof_indices.resize(dofs_per_cell);
    cell->get_dof_indices(data.local_dof_indices);

    Vector<double> dGamma(slips.size());
    Vector<double> gamma(slips.size());
    Vector<double> prior_gamma(slips.size());
    Vector<double> prior_gamma_abs(slips.size());
    vector<Tensor<1, dim> > gradGamma(slip_grads.size());
    Vector<double> current_hardening_gp(slips.size());
    Vector<double> current_gamma_abs_gp(slips.size());
    SymmetricTensor<4, dim> dsig_deps;

    // Elementsteifigkeitsmatrix
    for (unsigned int q_point = 0; q_point < n_q_points; ++q_point) {
        PerCrystalData<dim, nslip> *perCrystalData =
            &crystal_data[cell->active_fe_index()];
        for (unsigned int i = 0; i < dGamma.size(); i++) {
            gamma[i] = slips[i][q_point];
            prior_gamma[i] = prior_slips[i][q_point];
            prior_gamma_abs[i] = prior_slips_abs[i][q_point];

            current_hardening_gp[i] = current_hardening_cell[i][q_point];
            current_gamma_abs_gp[i] = current_gamma_abs_cell[i][q_point];

            dGamma[i] = slips[i][q_point] - prior_slips[i][q_point];
            gradGamma[i] = slip_grads[i][q_point];
        }

        // constant things on the gp
        SymmetricTensor<2, dim> stress;
        if (perCrystalData->isLocalFormulation()) {
            stress = localCP::getStress(strains[q_point], boundary->get_dt(),
                                        23 /*temperature*/, perCrystalData);
            dsig_deps = localCP::getDsigDeps(strains[q_point],
                                             boundary->get_dt(), 23 /*temperature*/, perCrystalData);
        } else {
            stress = GradientCP::get_stress(strains[q_point], gamma, /*Temperature*/
                                            23, perCrystalData);
            dsig_deps = GradientCP::getDsigDeps(perCrystalData);
        }
        double JxW = fe_values.JxW(q_point);

        for (unsigned int i = 0; i < dofs_per_cell; ++i) {
            const unsigned int component_i =
                compToCPSystem[fe_collection[cell->active_fe_index()].system_to_component_index(   
                                   i).first];                    //index of the only nonzero vector component of shape function i using the fe.system_to_component_index(i)
            if (component_i < dim) { //u - x dof

                SymmetricTensor<2, dim> grad_phi_i_u =
                    fe_values[extractor_displacement].symmetric_gradient(i,
                            q_point);
                for (unsigned int j = i; j < dofs_per_cell; ++j) {
                    const unsigned int component_j =
                        compToCPSystem[fe_collection[cell->active_fe_index()].system_to_component_index(
                                           j).first];
                    if (component_j < dim) { //uu dof

                        SymmetricTensor<2, dim> grad_phi_j_u =
                            fe_values[extractor_displacement].symmetric_gradient(
                                j, q_point);
                        data.cell_matrix(i, j) += (grad_phi_i_u * dsig_deps
                                                   * grad_phi_j_u * JxW);
                    } else { // u -Gamma dof
                        int slipsystem_j = component_j - dim;

                        const Tensor<0, dim> phi_j =
                            fe_values[extractor_slips[cell->active_fe_index()][slipsystem_j]].value(
                                j, q_point);

                        SymmetricTensor<2, dim> dsig_dgammaJ =
                            GradientCP::getDsigDgammaAlpha(slipsystem_j,
                                                           perCrystalData);

                        data.cell_matrix(i, j) += (grad_phi_i_u * dsig_dgammaJ
                                                   * phi_j * JxW);
                    }
                }
                // u dof rhs

                data.cell_rhs(i) += (/*rhs_values[q_point](component_i)
				 * scratch.fe_values.shape_value(i, q_point)*/
                                        -grad_phi_i_u * stress) * JxW;

                //				if (data.cell_rhs(i) != data.cell_rhs(i)) {
                //					deallog << "Nan in u-rhs(" << i << ")" << endl;
                //					throw std::bad_exception();
                //				}
            } else { // gamma - x dof
                int slipsystem_i = component_i - dim;

                Tensor<1, dim> grad_phi_i_g =
                    fe_values[extractor_slips[cell->active_fe_index()][slipsystem_i]].gradient(
                        i, q_point);
                const double phi_i_g =
                    fe_values[extractor_slips[cell->active_fe_index()][slipsystem_i]].value(
                        i, q_point);
                for (unsigned int j = i; j < dofs_per_cell; ++j) {
                    const unsigned int component_j =
                        compToCPSystem[fe_collection[cell->active_fe_index()].system_to_component_index(
                                           j).first];

                    if (component_j < dim) { //Gamma - u dof
                        SymmetricTensor<2, dim> grad_phi_j_u =
                            fe_values[extractor_displacement].symmetric_gradient(
                                j, q_point);

                        SymmetricTensor<2, dim> dtauI_deps =
                            GradientCP::getDtauAlphaDeps(slipsystem_i,
                                                         perCrystalData);

                        data.cell_matrix(i, j) += -(phi_i_g * dtauI_deps
                                                    * grad_phi_j_u) * JxW;
                    } else { // Gamma -Gamma dof, variable part
                        int slipsystem_j = component_j - dim;

                        Tensor<1, dim> grad_phi_j_g =
                            fe_values[extractor_slips[cell->active_fe_index()][slipsystem_j]].gradient(
                                j, q_point);
                        const double phi_j =
                            fe_values[extractor_slips[cell->active_fe_index()][slipsystem_j]].value(
                                j, q_point);

                        Tensor<2, dim> dxiI_dgradGammaJ =
                            GradientCP::getDxiAlphaDgradgammaBeta(
                                slipsystem_i, slipsystem_j, /*Temperature*/
                                23, perCrystalData);
                        double dtauI_dGammaJ =
                            GradientCP::getDtauAlphaDgammaBeta(slipsystem_i,
                                                               slipsystem_j, perCrystalData);
                        double dpiI_dGammaJ = GradientCP::getDpiAlphaDgammaBeta(
                                                  slipsystem_i, slipsystem_j, boundary->get_dt(), /*temperatur*/
                                                  23, dGamma, current_gamma_abs_gp,
                                                  current_hardening_gp, perCrystalData);
                        data.cell_matrix(i, j) += (contract3(grad_phi_i_g,
                                                             dxiI_dgradGammaJ, grad_phi_j_g)
                                                   - (phi_i_g * (dtauI_dGammaJ - dpiI_dGammaJ)
                                                      * phi_j)) * JxW;
                    }
                }
                // Gamma rhs
                Tensor<1, dim> xiI = GradientCP::get_XiAlpha(gradGamma,
                                     slipsystem_i, /*Temperature*/23, perCrystalData);
                double tauI = GradientCP::get_TauAlpha(slipsystem_i, stress,
                                                       perCrystalData);
                double piI = GradientCP::get_PiAlpha(slipsystem_i,
                                                     boundary->get_dt(), /*temperature*/23, dGamma,
                                                     current_hardening_gp, perCrystalData);
                //				deallog << "In assemble system gp(" << q_point << ") pi^"
                //						<< slipsystem_i << "(dt: " << boundary->get_dt()
                //						<< ", dGamma: " << dGamma << ") = " << piI << endl;
                //				deallog << (grad_phi_i * xiI) << " - " << phi_i * tauI << " + " << phi_i * piI  << endl;

                data.cell_rhs(i) += -((grad_phi_i_g * xiI)
                                      - phi_i_g * (tauI - piI)) * JxW;
                //				if (data.cell_rhs(i) != data.cell_rhs(i)) {
                //					deallog << "Nan in g-rhs(" << i << ")" << endl;
                //					throw std::bad_exception();
                //				}
            }
        }
    }
    // Just the lower triangle and diagonal of the elementmatrices are assembled
    // copy the entries in die upper triangle
    for (unsigned int i = 0; i < dofs_per_cell; ++i) {
        for (unsigned int j = 0; j < i; ++j) {
            data.cell_matrix[i][j] = data.cell_matrix[j][i];
        }
    }

    // Interfaceterm
    if (parameter.get("GB_behaviour").compare("detailed") == 0) {
        // just if face is on interface
        // just once per face, identify if neighbour cell allready visited. One has to identify the faces before assembling them
        PerCrystalData<dim, nslip> *cell_crystaldata =
            &crystal_data[cell->active_fe_index()];
        int domainID = cell->active_fe_index();

        for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; ++f) {
            data.is_grain_boundary[f] = false;
            if (cell->at_boundary(f)) {
                continue;// If the cell is at the boundary it cant be a grainboundary
            }
            PerCrystalData<dim, nslip> *neighbor_crystaldata =
                &crystal_data[cell->neighbor(f)->material_id()];
            int neighbor_domainID = cell->neighbor(f)->material_id();
            if (domainID == neighbor_domainID) {
                continue;// If both have the same MaterialID the face is inside a grain
            }
            // f is a Grainboundary Face, reinit FEFaceValues
            scratch.hp_fe_face_values.reinit(cell, f);
            const FEFaceValues<dim> &fe_face_values =
                scratch.hp_fe_face_values.get_present_fe_values();
            scratch.hp_fe_face_values_neighbor.reinit(cell->neighbor(f),
                    cell->neighbor_of_neighbor(f));
            const FEFaceValues<dim> &fe_face_values_neighbor =
                scratch.hp_fe_face_values_neighbor.get_present_fe_values();
            Assert(
                fe_face_values.n_quadrature_points == fe_face_values_neighbor.n_quadrature_points,
                ExcInternalError());

            if (data.cell_coupling_matrices[f].m() == dofs_per_cell) {
                data.cell_coupling_matrices[f] = 0;
            } else {
                data.cell_coupling_matrices[f].reinit(dofs_per_cell,
                                                      dofs_per_cell);
            }

            //Values of slips on both sides of the interface
            vector<vector<double> > slips_cellfaces(nslip,
                                                    vector<double>(fe_face_values.n_quadrature_points));
            vector<vector<double> > prior_slips_cellfaces(nslip,
                    vector<double>(fe_face_values.n_quadrature_points));
            vector<vector<double> > slips_neighborfaces(nslip,
                    vector<double>(
                        fe_face_values_neighbor.n_quadrature_points));

            for (unsigned int i = 0;
                    i < extractor_slips[cell->active_fe_index()].size(); i++) {
                fe_face_values[extractor_slips[cell->active_fe_index()][i]].get_function_values(
                    current_solution, slips_cellfaces[i]);

                fe_face_values[extractor_slips[cell->active_fe_index()][i]].get_function_values(
                    prior_solution, prior_slips_cellfaces[i]);

                fe_face_values_neighbor[extractor_slips[cell->neighbor(f)->active_fe_index()][i]].get_function_values(
                    current_solution, slips_neighborfaces[i]);
            }

            // Find correlating Gauspoints
            const vector<Point<dim> > gp_coords_cell =
                fe_face_values.get_quadrature_points();
            const vector<Point<dim> > gp_coords_neighbor =
                fe_face_values_neighbor.get_quadrature_points();
            vector<unsigned int> gp_map_cell_to_neighbor(gp_coords_cell.size());
            for (unsigned int i = 0; i < gp_coords_cell.size(); i++) {
                for (unsigned int j = 0; j < gp_coords_cell.size(); j++) {
                    Tensor<1, dim> distance_tmp = gp_coords_cell[i]
                                                  - gp_coords_neighbor[j];
                    if (distance_tmp.norm() < 1e-8) {
                        gp_map_cell_to_neighbor[i] = j;
                        break;
                    }
                }
                if (i != gp_map_cell_to_neighbor[i]) {
                    deallog << "Map gp cell i: " << i << " to neighbor: "
                            << gp_map_cell_to_neighbor[i] << endl;
                }
            }

            Tensor<1, dim> normal_cell;
            Tensor<1, dim> normal_neighbor;
            //Loop over all Gauspoints
            for (unsigned int q_point = 0;
                    q_point < fe_face_values.n_quadrature_points; ++q_point) {
                double JxW_cell = fe_face_values.JxW(q_point);
                //			double JxW_neighbor = fe_face_values_neighbor.JxW(q_point);
                //Loop over Dofs of cell
                for (unsigned int i = 0; i < dofs_per_cell; ++i) {
                    const unsigned int component_i =
                        compToCPSystem[fe_collection[cell->active_fe_index()].system_to_component_index(
                                           i).first];
                    if (component_i >= dim) {
                        int slipsystem_i = component_i - dim;
                        const double phi_i =
                            fe_face_values[extractor_slips[cell->active_fe_index()][slipsystem_i]].value(
                                i, q_point);
                        normal_cell = fe_face_values.normal_vector(q_point);
                        normal_neighbor = fe_face_values_neighbor.normal_vector(
                                              q_point);
                        for (unsigned int j = 0; j < dofs_per_cell; ++j) {
                            const unsigned int component_j =
                                compToCPSystem[fe_collection[cell->active_fe_index()].system_to_component_index(
                                                   j).first];
                            const unsigned int component_j_neighbor =
                                compToCPSystem[fe_collection[cell->neighbor(
                                                                 f)->active_fe_index()].system_to_component_index(
                                                   j).first];
                            Assert(component_j==component_j_neighbor,
                                   StandardExceptions::ExcMessage("component_j of cell and neighbor unequal"));
                            if (component_j >= dim) { //C_AA and C_AB
                                int slipsystem_j = component_j - dim;
                                int slipsystem_j_neighbor = component_j_neighbor
                                                            - dim;
                                const double phi_j_cell =
                                    fe_face_values[extractor_slips[cell->active_fe_index()][slipsystem_j]].value(
                                        j, q_point);
                                const double phi_j_neighbor =
                                    fe_face_values_neighbor[extractor_slips[cell->neighbor(
                                                                f)->active_fe_index()][slipsystem_j_neighbor]].value(
                                        j, q_point);
                                const double C_CellCell = PerCrystalData<dim,
                                             nslip>::get_C_AB_ab(cell_crystaldata,
                                                                 cell_crystaldata, normal_cell,
                                                                 normal_cell, slipsystem_i,
                                                                 slipsystem_j); //C_AA
                                const double C_CellNeighbor = PerCrystalData<
                                                              dim, nslip>::get_C_AB_ab(
                                                                  cell_crystaldata, neighbor_crystaldata,
                                                                  normal_cell, normal_neighbor,
                                                                  slipsystem_i, slipsystem_j_neighbor); //C_AB

                                const double d_xi_iA_x_n_dgamma_jA =
                                    parameter.get_double("lambda_gb")
                                    * C_CellCell;
                                const double d_xi_iA_x_n_dgamma_jB =
                                    parameter.get_double("lambda_gb")
                                    * C_CellNeighbor;

                                data.cell_matrix(i, j) += -phi_i
                                                          * (-d_xi_iA_x_n_dgamma_jA) * phi_j_cell
                                                          * JxW_cell;
                                data.cell_coupling_matrices[f](i, j) += -phi_i
                                                                        * (-d_xi_iA_x_n_dgamma_jB)
                                                                        * phi_j_neighbor * JxW_cell;

                                // Tangent of dissipative part
                                if (slipsystem_i == slipsystem_j) {
                                    const double dgamma_cellface =
                                        prior_slips_cellfaces[slipsystem_i][q_point]
                                        - slips_cellfaces[slipsystem_i][q_point];
                                    const double d_xi_iA_x_n_dgamma_jA_diss =
                                        parameter.get_double("S_gb")
                                        * parameter.get_double(
                                            "p_gb")
                                        * pow(abs(dgamma_cellface),
                                              parameter.get_double(
                                                  "p_gb")
                                              - 1);
                                    data.cell_matrix(i, j) += -phi_i
                                                              * (-d_xi_iA_x_n_dgamma_jA_diss)
                                                              * phi_j_cell * JxW_cell;
                                }
                            }
                        }

                        // Interface RHS
                        double xi_x_n_cell = 0;
                        for (int beta = 0; beta < nslip; beta++) {
                            const double C_CellCell =
                                PerCrystalData<dim, nslip>::get_C_AB_ab(
                                    cell_crystaldata, cell_crystaldata,
                                    normal_cell, normal_cell,
                                    slipsystem_i, beta); //C_AA
                            const double C_CellNeighbor = PerCrystalData<dim,
                                         nslip>::get_C_AB_ab(cell_crystaldata,
                                                             neighbor_crystaldata, normal_cell,
                                                             normal_neighbor, slipsystem_i, beta); //C_AB

                            xi_x_n_cell +=
                                -parameter.get_double("lambda_gb")
                                * (C_CellCell
                                   * slips_cellfaces[beta][q_point]
                                   + C_CellNeighbor
                                   * slips_neighborfaces[beta][q_point]);

                            //						deallog << "C_" << slipsystem_i << slipsystem_i << ": "
                            //								<< C_CellCell << "             C_"
                            //								<< slipsystem_i << beta << ": "
                            //								<< C_CellNeighbor << endl;
                        }
                        const double dgamma_cellface =
                            prior_slips_cellfaces[slipsystem_i][q_point]
                            - slips_cellfaces[slipsystem_i][q_point];
                        const double rhs_gb_diss = -parameter.get_double("S_gb")
                                                   * pow(abs(dgamma_cellface),
                                                         parameter.get_double("p_gb"))
                                                   * sgn(dgamma_cellface);
                        xi_x_n_cell += -rhs_gb_diss;
                        // -rhs
                        data.cell_rhs(i) += -(-phi_i * xi_x_n_cell * JxW_cell);
                        //					deallog << "Xi contracted with n: " << xi_x_n_cell << endl;
                    }
                }
            }
            data.local_dof_indices_neighbours[f].resize(dofs_per_cell);
            cell->neighbor(f)->get_dof_indices(
                data.local_dof_indices_neighbours[f]);
            data.is_grain_boundary[f] = true;
        }
    }
}

template<int dim, int nslip>
void ElastoplasticProblem<dim, nslip>::sparseindices_on_one_cell(
    const typename hp::DoFHandler<dim>::active_cell_iterator &cell,
    PerTaskData<dim> &data) {

    const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    data.local_dof_indices.resize(dofs_per_cell);
    cell->get_dof_indices(data.local_dof_indices);

    // Interfaceindices

    // just if face is on interface
    int domainID = cell->material_id();

    for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; ++f) {
        data.is_grain_boundary[f] = false;
        if (cell->at_boundary(f)) {
            continue; // If the cell is at the boundary it cant be a grainboundary
        }
        int neighbor_domainID = cell->neighbor(f)->material_id();
        if (domainID == neighbor_domainID) {
            continue; // If both have the same MaterialID the face is inside a grain
        }

        data.local_dof_indices_neighbours[f].resize(dofs_per_cell);
        cell->neighbor(f)->get_dof_indices(
            data.local_dof_indices_neighbours[f]);
        data.is_grain_boundary[f] = true;
    }

}

template<int dim, int nslip>
void ElastoplasticProblem<dim, nslip>::copy_local_to_global(
    const PerTaskData<dim> &data) {
    //	const unsigned int dofs_per_cell = fe.dofs_per_cell;
    //	for (unsigned int i = 0; i < dofs_per_cell; ++i) {
    //		for (unsigned int j = 0; j < dofs_per_cell; ++j)
    //			system_matrix.add(data.local_dof_indices[i],
    //					data.local_dof_indices[j], data.cell_matrix(i, j));
    //
    //		system_rhs(data.local_dof_indices[i]) += data.cell_rhs(i);
    //	}
    constraints.distribute_local_to_global(data.cell_matrix, data.cell_rhs,
                                           data.local_dof_indices, system_matrix, system_rhs);
    if (parameter.get("GB_behaviour").compare("detailed") == 0) {
        for (unsigned int i = 0; i < data.is_grain_boundary.size(); i++) {
            if (data.is_grain_boundary[i]) {
                constraints.distribute_local_to_global(
                    data.cell_coupling_matrices[i], data.local_dof_indices,
                    data.local_dof_indices_neighbours[i], system_matrix);
            }
        }
    }
}

template<int dim, int nslip>
void ElastoplasticProblem<dim, nslip>::apply_constraints() {
    // make hanging node constraints
    DoFTools::make_hanging_node_constraints(dof_handler, constraints);

    // Micro hard grain boundaries
    if (parameter.get("GB_behaviour").compare("microhard") == 0) {
        //Zero microhard-grainboundarys
        for (unsigned int i = 0; i < dof_is_grainboundary.size(); i++) {
            if (dof_is_grainboundary[i]) {
                constraints.add_line(i);
                constraints.set_inhomogeneity(i, 0);
            }
        }
    }

    // If Tangent is calculated no plasticity is considered
    if (parameter.get("mode").compare("calculateTangent") == 0) {
        //		deallog << "Set slip dofs to 0, just elasticity" << endl;
        typename hp::DoFHandler<dim>::active_cell_iterator cell =
            dof_handler.begin_active(), endc = dof_handler.end();
        for (; cell != endc; ++cell) {
            const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;
            vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

            cell->get_dof_indices(local_dof_indices);
            for (unsigned int i = 0; i < dofs_per_cell; i++) {
                const unsigned int component_dof =
                    compToCPSystem[fe_collection[cell->active_fe_index()].system_to_component_index(
                                       i).first];
                if (component_dof >= dim) {
                    int dof = local_dof_indices[i];
                    constraints.add_line(dof);
                    constraints.set_inhomogeneity(dof, 0);
                }
            }
        }
    }

    if (parameter.get("loadtype").compare("periodic") == 0) {
        BC::Periodic::constrain(boundary, dofdata, dof_reference_node,
                                constraints, dof_handler);
    }
    if (parameter.get("loadtype").compare("periodic_uniaxial_tension") == 0) {
        BC::Periodic::constrainUniaxialTension(boundary, dofdata,
                                               dof_reference_node, constraints, dof_handler,
                                               parameter.get("straindirection"));
    }
    if (parameter.get("loadtype").compare("periodic_plane_stress") == 0) {
        vector<string> loadDirs = { parameter.get("planeStressDir1"),
                                    parameter.get("planeStressDir2")
                                  };
        BC::Periodic::constrainPlaneStress(boundary, dofdata,
                                           dof_reference_node, constraints, dof_handler, loadDirs);
    }
    if (parameter.get("loadtype").compare("periodic_stressfree") == 0
            || parameter.get("loadtype").compare("periodic_stressfree_cooling")
            == 0) {
        BC::Periodic::constrainStressfree(boundary, dofdata, dof_reference_node,
                                          constraints, dof_handler);
    }
    /*
     * linear
     */
    if (parameter.get("loadtype").compare("linear") == 0) {
        BC::Linear::constrain(boundary, constraints, dof_handler);
    }

    if (parameter.get("loadtype").compare("linear_uniaxial_tension") == 0) {
        BC::Linear::constrainUniaxialTension(boundary, dofdata,
                                             dof_reference_node, constraints, dof_handler,
                                             parameter.get("straindirection"));
    }

    if (parameter.get("loadtype").compare("linear_plane_stress") == 0) {
        vector<string> loadDirs = { parameter.get("planeStressDir1"),
                                    parameter.get("planeStressDir2")
                                  };
        BC::Linear::constrainPlaneStress(boundary, dofdata, dof_reference_node,
                                         constraints, dof_handler, loadDirs);
    }

    /*
     * linear in isotropic plane and periodic in buildingdirection
     */

    if (parameter.get("loadtype").compare("linear_CG") == 0) {
        BC::Linear::CG::constrain(boundary, dofdata, dof_reference_node,
                                  constraints, dof_handler);
    }

    if (parameter.get("loadtype").compare("linear_CG_uniaxial_tension") == 0) {
        BC::Linear::CG::constrainUniaxialTension(boundary, dofdata,
                dof_reference_node, constraints, dof_handler,
                parameter.get("straindirection"));
    }

    if (parameter.get("loadtype").compare("linear_CG_plane_stress") == 0) {
        vector<string> loadDirs = { parameter.get("planeStressDir1"),
                                    parameter.get("planeStressDir2")
                                  };
        BC::Linear::CG::constrainPlaneStress(boundary, dofdata,
                                             dof_reference_node, constraints, dof_handler, loadDirs);
    }

    if (parameter.get("loadtype").compare("testcase_andrew_microhard") == 0
            || parameter.get("loadtype").compare("testcase_grain_boundary")
            == 0) {
        BC::TestcaseShear::microHard(dof_handler, constraints, boundary, nslip,
                                     maxCrystals);
    }
    if (parameter.get("loadtype").compare("testcase_gottschalk") == 0) {
        BC::TestcaseShear::microHardGottschalk(dof_handler, constraints,
                                               boundary, nslip, maxCrystals);
    }
    if (parameter.get("loadtype").compare("testcase_andrew_microfree") == 0) {
        BC::TestcaseShear::microFree(dof_handler, constraints, boundary, nslip,
                                     maxCrystals);
    }

}

template<int dim, int nslip>
void ElastoplasticProblem<dim, nslip>::assemble_system() {
    timer.enter_section("Assemble");

    constraints.clear();
    apply_constraints();
    constraints.close();

    setup_sparsity_pattern();

    // Assemble Systemmatrix constant
    typename hp::DoFHandler<dim>::active_cell_iterator cell =
        dof_handler.begin_active(), end_cell = dof_handler.end();

    //------------------------------------------------------
    ScratchData<dim> scratch(fe_collection, quadrature_collection,
                             update_values | update_gradients | update_JxW_values
                             | update_quadrature_points, face_quadrature_collection,
                             update_values | update_JxW_values | update_normal_vectors
                             | update_quadrature_points);
    //parallel with workstream
    if (parameter.get_integer("nCoresAssembly") > 1) {
        PerTaskData<dim> per_task_data(fe_collection);
        //		ScratchData<dim> scratch(fe_collection, quadrature_collection,
        //				update_values | update_gradients | update_JxW_values
        //						| update_quadrature_points, face_quadrature_collection,
        //				update_values | update_JxW_values | update_normal_vectors
        //						| update_quadrature_points);
        WorkStream::run(cell, end_cell, *this,
                        &ElastoplasticProblem<dim, nslip>::assemble_on_one_cell,
                        &ElastoplasticProblem<dim, nslip>::copy_local_to_global,
                        scratch, per_task_data);
        //------------------------------------------------------
    } else {
        //------------------------------------------------------
        // without workstream
        //		ScratchData<dim> scratch(fe_collection, quadrature_collection,
        //				update_values | update_gradients | update_JxW_values
        //						| update_quadrature_points, face_quadrature_collection,
        //				update_values | update_JxW_values | update_normal_vectors
        //						| update_quadrature_points);
        for (; cell != end_cell; ++cell) {
            PerTaskData<dim> per_task_data(fe_collection);
            assemble_on_one_cell(cell, scratch, per_task_data);
            copy_local_to_global(per_task_data);
        }
        //------------------------------------------------------
    }
    timer.exit_section("Assemble");
}

template<int dim, int nslip>
void ElastoplasticProblem<dim, nslip>::solve() {
    timer.enter_section("Solve without Assemble");
    SolverControl solver_control(
        std::max(
            parameter.get_double("maxCGIterationsRel")
            * dof_handler.n_dofs(), 1000.0),
        parameter.get_double("tolCGsolver"));
    SolverCG<> cg(solver_control);
    cg.connect_condition_number_slot(
        std_cxx11::bind(
            &ElastoplasticProblem<dim, nslip>::slotConditionNumber,
            this, std_cxx11::_1));
    PreconditionSSOR<> preconditioner;

    double norm_res = -1;
    double last_residuum = -1;
    double regIncreaseFactor = parameter.get_double("regularIncreaseF");

    int loadstep = 0;
    boundary->first_timestep();
    timer.exit_section("Solve without Assemble");
    postprocessing(loadstep);
    timer.enter_section("Solve without Assemble");

    // vector for old solutions, used for mechanical predictor
    vector<Vector<double>> solution_history; //TODO alloc
    vector<double> times_history;
    update_solutionhistory(solution_history, current_solution, times_history,
                           boundary->get_time(), parameter);
    Vector<double> newton_update(dof_handler.n_dofs()); // Verschiebungsvektor des Newtonschritts
    Vector<double> d_newtonupdate(current_solution.size()); //just for displaying the newtonstepsize
    while (boundary->hasNextIncrement()) { //loadstep loop
        int total_iterations = 0;
        boundary->next_timestep();
        loadstep++;

        //Newton Raphson Verfahren
        unsigned int newton_iteration = 1;
        double first_residuum = 0;

        // save last converged solution
        prior_hardening = current_hardening;
        prior_solution_abs = current_solution_abs;
        prior_solution = current_solution;

        // insert delta solution in history, remove unnecessary entries
        update_solutionhistory(solution_history, current_solution,
                               times_history, boundary->get_time(), parameter);
        current_solution += mechanical_prediction_get_delta(solution_history,
                            times_history);
        update_hardening();

        deallog << "Loadstep " << loadstep << " @ " << boundary->get_time()
                << endl;
        bool isConverged;
        double regularizationFactor = 1;
        Vector<double> initial_current_solution(current_solution);
        Vector<double> initial_prior_solution(prior_solution);
        Vector<double> initial_current_solution_abs(current_solution_abs);
        Vector<double> initial_prior_solution_abs(prior_solution_abs);
        Vector<double> initial_current_hardening(current_hardening);
        Vector<double> initial_prior_hardening(prior_solution_abs);
        int regLowerLoops = 1;
        int changed_slips = 0;
        do {
            n_iter_remove_me = newton_iteration;
            isConverged = true; // will be set to false if not converging
            try {
                //				Vector<double> d_newtonupdate(current_solution.size()); //just for displaying the newtonstepsize
                d_newtonupdate = 0;
                print_conv_header();
                //				double relaxation_factor = 1;
                //				vector<Vector<double>> newton_history;
                do { //Newton Raphson Loop
                    deallog << "Iteration " << setw(2) << regLowerLoops << "/"
                            << setw(2) << newton_iteration << " " << flush;

                    // set relaxationFactors
                    for (unsigned int i = 0; i < crystal_data.size(); i++) {
                        crystal_data[i].set_regularisationFactor(
                            regularizationFactor);
                    }

                    //compute rhs, system_matrix
                    //					Vector<double> newton_update(dof_handler.n_dofs()); // Verschiebungsvektor des Newtonschritts
                    newton_update = 0;
                    //					system_rhs.reinit(dof_handler.n_dofs());
                    system_rhs = 0;
                    timer.exit_section("Solve without Assemble");
                    deallog << "Assemble " << flush;
                    assemble_system();
                    timer.enter_section("Solve without Assemble");

                    //					saveSystemMatrix(newton_iteration);

                    vector<double> res = calculateResiduals(system_rhs);
                    vector<double> res_dSolution = calculateResiduals(
                                                       d_newtonupdate);
                    last_residuum = norm_res;

                    /*
                     * The rhs of constrained DOFs are 0 because I set the update_vector
                     * at the first Newton step to 0
                     */
                    norm_res = system_rhs.l2_norm();
                    if (newton_iteration == 1 && regLowerLoops == 1) { // erster Newtonschritt
                        first_residuum = norm_res;
                        if (first_residuum < 1) {
                            first_residuum = 1;
                        }
                    }

                    if (regularizationFactor > 1 - 1e-6
                            && norm_res
                            < parameter.get_double("globalTolNewton")
                            * first_residuum) {
                        static const unsigned int c_width = 12;
                        deallog << setw(14) << "Converged " << "| "
                                << setw(c_width) << norm_res << setw(c_width)
                                << res[0] << setw(c_width) << res[1]
                                << setw(c_width) << res_dSolution[0]
                                << setw(c_width) << res_dSolution[1]

                                //								<< setw(c_width) << maxValueEvolutionEquation()
                                << setw(c_width) << solver_control.last_step()
                                << setw(c_width) << conditionNumber

                                //								<< setw(c_width) << relaxation_factor
                                << setw(c_width) << regularizationFactor
                                << setw(c_width)
                                << log10(norm_res) / log10(last_residuum)
                                << setw(c_width) << regIncreaseFactor
                                << setw(c_width) << changed_slips << endl;
                        double d_dissipation = calculate_dissipation();
                        deallog << "Dissipation actual step: " << d_dissipation
                                << " Total iteration: " << total_iterations
                                << endl;
                        print_conv_footer();
                        break;
                    }
                    if (norm_res != norm_res) {
                        deallog << "--!--Res == NAN--!--";
                        throw std::bad_exception();
                    }

                    //			int depth = deallog.depth_console(0);

                    // iterative solver (CG)
                    total_iterations++;
                    deallog << "Solve " << flush;
                    timer.enter_section("Only Solve");
                    preconditioner.initialize(system_matrix, 1);
                    cg.solve(system_matrix, newton_update, system_rhs,
                             preconditioner);
                    timer.exit_section("Only Solve");
                    constraints.distribute(newton_update);
                    //			deallog.depth_console(depth);
                    deallog.pop();
                    deallog << "Solved " << flush;

                    //					update_newtonhistory(newton_history, newton_update);
                    //					relaxation_factor = calculate_relaxation_factor(
                    //							newton_history, relaxation_factor);
                    static const unsigned int c_width = 12;
                    deallog << setw(1) << " " << "| " << setw(c_width)
                            << norm_res << setw(c_width) << res[0]
                            << setw(c_width) << res[1] << setw(c_width)
                            << res_dSolution[0] << setw(c_width)
                            << res_dSolution[1] << setw(c_width)

                            //							<< maxValueEvolutionEquation() << setw(c_width)
                            << solver_control.last_step() << setw(c_width)
                            << conditionNumber << setw(c_width)

                            //							<< relaxation_factor << setw(c_width)
                            << regularizationFactor << setw(c_width)
                            << log10(norm_res) / log10(last_residuum)
                            << setw(c_width) << regIncreaseFactor
                            << setw(c_width) << changed_slips << endl;
                    deallog.push("DEAL");
                    //					newton_update *= relaxation_factor;

                    changed_slips = newtonCorrection(newton_update);
                    // update solution
                    current_solution += newton_update;
                    // update global current_hardening vector and abs_solution
                    update_hardening();

                    boundary->deactivate_inhomogenious_bc();
                    newton_iteration++;
                    d_newtonupdate = newton_update;

                    if (newton_iteration
                            > parameter.get_integer("maxGlobalIterations")) {
                        deallog << "Iterationcounter NewtonRaphson > "
                                << parameter.get_integer("maxGlobalIterations");
                        throw std::bad_exception();
                    }

                    // increase regularisation factor
                    if (norm_res
                            < parameter.get_double("NewtonTolIncReF")
                            * first_residuum
                            && regularizationFactor < 1) {
                        double convergF_last_step = log10(norm_res)
                                                    / log10(last_residuum);
                        if (newton_iteration > 2) {
                            double targetConvergence = parameter.get_double(
                                                           "targetConvergence");
                            regIncreaseFactor *= pow(
                                                     convergF_last_step / targetConvergence,
                                                     0.5);
                            regIncreaseFactor = max(regIncreaseFactor,
                                                    pow(
                                                        parameter.get_double(
                                                            "regularIncreaseF"), 0.3));
                        }
                        regularizationFactor = min(
                                                   regularizationFactor * regIncreaseFactor, 1.0);
                        newton_iteration = 1;
                    }
                } while (true);
            } catch (...) {
                if (parameter.get_double("regularisationF") == 1) {
                    deallog
                            << "Not converging and regularization softening disabled."
                            << endl;
                    throw std::bad_exception();
                }
                if (regLowerLoops > 6) {
                    deallog
                            << "Not converging too many regularization lowering loops."
                            << endl;
                    throw std::bad_exception();
                }
                deallog << " - Lower Regularisation Factor" << endl;
                regLowerLoops++;
                regularizationFactor = pow(
                                           parameter.get_double("regularisationF"),
                                           regLowerLoops - 1);
                current_solution = initial_current_solution;
                prior_solution = initial_prior_solution;
                current_solution_abs = initial_current_solution_abs;
                prior_solution_abs = initial_prior_solution_abs;
                current_hardening = initial_current_hardening;
                prior_hardening = initial_prior_hardening;
                isConverged = false;
                newton_iteration = 1;
            } // End NR
        } while (!isConverged); // runs until coverged

        //		quadrature_point_history_copy_newtonstep_to_history(); // renew internal variables
        timer.exit_section("Solve without Assemble");
        postprocessing(loadstep);
        timer.enter_section("Solve without Assemble");
    }
    timer.exit_section("Solve without Assemble");
}

/*
 * mechanical prediction according to paper Erbts, Duester
 * Accelerated staggered coupling schemes for problems of thermoelasticity at finite strains
 */
template<int dim, int nslip>
Vector<double> ElastoplasticProblem<dim, nslip>::mechanical_prediction_get_delta(
    vector<Vector<double>> &solution_history,
    vector<double> &times_history) {

    Vector<double> a = get_MP_coefficients(times_history);
    deallog << setprecision(3);
    deallog << " x(t=" << times_history[0] << ") = ";
    for (unsigned int i = 1; i < times_history.size(); ++i) {
        deallog << a[i - 1] << " x(t=" << times_history[i] << "); ";
    }
    deallog << setprecision(6);

    if (!boundary->isMonotonicLoad()) {
        a.reinit(0);
        deallog << "no MP because of changing load ";
    }

    if (a.size() == 0) {
        a.reinit(1);
        a[0] = 1;
    }
    a[0] -= 1;

    Vector<double> delta_solution;
    delta_solution.reinit(current_solution.size());
    Vector<double> tmp;
    tmp.reinit(current_solution.size());
    for (unsigned int i = 0; i < a.size(); ++i) {
        tmp = solution_history[i];
        tmp *= a[i];
        delta_solution += tmp;
    }

    constraints.clear();
    apply_constraints();
    constraints.close();
    constraints.distribute(delta_solution);
    boundary->deactivate_inhomogenious_bc();
    //	current_solution += delta_solution;
    return delta_solution;
}

template<int dim, int nslip>
void ElastoplasticProblem<dim, nslip>::saveSystemMatrix(
    int loadstep/*Vector<double> &dSolution*/) {
    if (dof_handler.n_dofs() > 100000) {
        return;
    }

    {
        string filename;
        std::stringstream sstm;
        sstm << "Results/SystemMatrix-" << loadstep << ".txt";
        filename = sstm.str();
        std::ofstream File(filename.c_str(), std::ios::out);
        system_matrix.print(File);
    }

    Vector<double> dofType(dof_handler.n_dofs());

    for (unsigned int i = 0; i < dof_handler.n_dofs(); i++) {
        if (dofdata.count(i) != 0) { // u_dof
            dofType[i] = 0;
        } else { // gammadof
            dofType[i] = 1;
        }
    }

    {
        string filename;
        std::stringstream sstm;
        sstm << "Results/dofType-" << loadstep << ".txt";
        filename = sstm.str();
        std::ofstream File(filename.c_str(), std::ios::out);
        dofType.print(File);
    }

    {
        string filename;
        std::stringstream sstm;
        sstm << "Results/Solution-" << loadstep << ".txt";
        filename = sstm.str();
        std::ofstream File(filename.c_str(), std::ios::out);
        current_solution.print(File);
    }

    {
        string filename;
        std::stringstream sstm;
        sstm << "Results/rhs-" << loadstep << ".txt";
        filename = sstm.str();
        std::ofstream File(filename.c_str(), std::ios::out);
        system_rhs.print(File);
    }
}

//template<int dim, int nslip>
//int ElastoplasticProblem<dim, nslip>::count_changed_active_systems() {
//	typename hp::DoFHandler<dim>::active_cell_iterator cell =
//			dof_handler.begin_active(), endc = dof_handler.end();
////	const unsigned int n_q_points = quadrature_formula.size();
////	int n_active_systems = 0;
////	for (; cell != endc; ++cell) {
////		Material<dim> *local_quadrature_points_data = reinterpret_cast<Material<
////				dim>*>(cell->user_pointer());
////		for (unsigned int q_point = 0; q_point < n_q_points; ++q_point) {
////			n_active_systems +=
////					local_quadrature_points_data[q_point].get_n_changed_active_systems();
////		}
////	}
//	return -1;
//}

template<int dim, int nslip>
void ElastoplasticProblem<dim, nslip>::createMesh() {
    if (parameter.get("mesh").compare("simple_mesh") == 0) {

        // One stripe of elements
        if (parameter.get("loadtype").compare("testcase_andrew_microhard") == 0
                || parameter.get("loadtype").compare(
                    "testcase_andrew_microfree") == 0
                || parameter.get("loadtype").compare("testcase_grain_boundary")
                == 0) {
            testcases::doubleSlip::createMesh(parameter, triangulation,
                                              dof_handler);

        } else if (parameter.get("loadtype").compare("testcase_gottschalk")
                   == 0) {
            testcases::gottschalk::createMesh(parameter, triangulation,
                                              dof_handler);

        } else {
            // a cube like the one for GAMM16
            testcases::GAMM16::createMesh(parameter, triangulation, dof_handler,
                                          geometric_range);
        }
    } else { // no simple mesh (hypercube), read in a .msh file
        RVE::createMesh(parameter, triangulation, dof_handler, geometric_range);
    }

    // count crystals and set active_fe_indices
    set<int> active_fe_indices;
    typename hp::DoFHandler<dim>::active_cell_iterator cell =
        dof_handler.begin_active(), endc = dof_handler.end();
    for (; cell != endc; ++cell) {
        cell->set_active_fe_index(cell->material_id());

        if (active_fe_indices.count(cell->material_id()) == 0) {
            active_fe_indices.insert(cell->material_id());
        }
    }
    maxCrystals = active_fe_indices.size();
    deallog << "Number of Crystals: " << maxCrystals << endl;

}

template<int dim, int nslip>
void ElastoplasticProblem<dim, nslip>::print_detailed_crystal_informations() {
    if (!parameter.get_bool("print_crystal_data")) {
        return;
    }

    for (unsigned int i = 0; i < crystal_data.size(); i++) {
        for (int j = 0; j < nslip; j++) {
            vector<double> eulerAngles = crystal_data[i].getEulerAngles();
            deallog << "Eulerangles: " << eulerAngles[0] << "/"
                    << eulerAngles[1] << "/" << eulerAngles[2] << " " << "ID: "
                    << crystal_data[i].getDomainID() << "   s[" << i << "] = "
                    << crystal_data[i].get_s(j) << "   m[" << i << "] = "
                    << crystal_data[i].get_m(j) << "   P[" << i << "] = "
                    << crystal_data[i].get_P(j) << endl;
        }
    }

    typename hp::DoFHandler<dim>::active_cell_iterator cell =
        dof_handler.begin_active(), end_cell = dof_handler.end();
    for (; cell != end_cell; cell++) {

        // Interfaceindices

        // just if face is on interface
        // just once per face, identify if neighbour cell already visited. One has to identify the faces before assembling them
        PerCrystalData<dim, nslip> *cell_crystaldata =
            &crystal_data[cell->active_fe_index()];
        int domainID = cell->active_fe_index();

        for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; ++f) {
            if (cell->at_boundary(f)) {
                continue;// If the cell is at the boundary it cant be a grainboundary
            }
            PerCrystalData<dim, nslip> *neighbor_crystaldata =
                &crystal_data[cell->neighbor(f)->active_fe_index()];
            int neighbor_domainID = cell->neighbor(f)->active_fe_index();
            if (domainID == neighbor_domainID) {
                continue;// If both have the same MaterialID the face is inside a grain
            }

            // Is a grainboundary, write data to console
            deallog << "C_" << domainID << "-" << neighbor_domainID << endl;
            deallog
                    << "---------------------------------------------------------"
                    << endl;
            hp::FEFaceValues<dim> hp_fe_face_values(fe_collection,
                                                    face_quadrature_collection, update_normal_vectors);
            hp_fe_face_values.reinit(cell, f);
            const FEFaceValues<dim> &fe_face_values =
                hp_fe_face_values.get_present_fe_values();
            const Tensor<1, dim> normal_cell = fe_face_values.normal_vector(0);
            for (int alpha = 0; alpha < nslip; alpha++) {
                for (int beta = 0; beta < nslip; beta++) {
                    deallog << setw(12)
                            << PerCrystalData<dim, nslip>::get_C_AB_ab(
                                cell_crystaldata, cell_crystaldata,
                                normal_cell, normal_cell, alpha, beta);
                }
                for (int beta = 0; beta < nslip; beta++) {
                    deallog << setw(12)
                            << PerCrystalData<dim, nslip>::get_C_AB_ab(
                                cell_crystaldata, neighbor_crystaldata,
                                normal_cell, -normal_cell, alpha, beta);
                }
                deallog << endl;
            }
            for (int alpha = 0; alpha < nslip; alpha++) {
                for (int beta = 0; beta < nslip; beta++) {
                    deallog << setw(12)
                            << -PerCrystalData<dim, nslip>::get_C_AB_ab(
                                neighbor_crystaldata, cell_crystaldata,
                                -normal_cell, -normal_cell, alpha, beta);
                }
                for (int beta = 0; beta < nslip; beta++) {
                    deallog << setw(12)
                            << -PerCrystalData<dim, nslip>::get_C_AB_ab(
                                neighbor_crystaldata, neighbor_crystaldata,
                                -normal_cell, normal_cell, alpha, beta);
                }
                deallog << endl;
            }
            deallog
                    << "---------------------------------------------------------"
                    << endl;
        }
    }

}

/*! \brief Identify the dofs which are on grain boundarys
 *
 *  The grainboundary dofs are marked with an true in the vector
 *  dof_is_grainboundary and can be used to apply micro-hard conditions
 */

template<int dim, int nslip>
void ElastoplasticProblem<dim, nslip>::identify_grainboundarys() {

    typename hp::DoFHandler<dim>::active_cell_iterator cell =
        dof_handler.begin_active(), endc = dof_handler.end();
    for (; cell != endc; ++cell) {

        for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; f++) {

            // If the face is at the boundary it cant be a grainboundary
            if (cell->at_boundary(f)) {
                continue;
            }

            //			int neighbor_domainID = cell->neighbor(f)->material_id();
            //			if (cell_domainID == neighbor_domainID) {
            //				continue; // If both have the same MaterialID the face is inside a grain
            //			}

            // If both have the same MaterialID the face is inside a grain
            if (cell->material_id() == cell->neighbor(f)->material_id()) {
                continue;
            }

            // f is a Grainboundary Face, all slip_dofs are relevant
            //		//Loop over vertices
            for (unsigned int j = 0; j < GeometryInfo<dim>::vertices_per_face;
                    j++) {
                for (int k = dim; k < dim + nslip; k++) { // only slipdofs
                    int dof = cell->face(f)->vertex_dof_index(j, k,
                              cell->active_fe_index());

                    dof_is_grainboundary[dof] = true;

                } //End of dof loop
            } //End of vertex loop
        } //End of face loop
    } // End of cell loop
}

template<int dim, int nslip>
void ElastoplasticProblem<dim, nslip>::calculateMacroscopicTangent() {
    if (parameter.get_double("displacement") > 0) {
        deallog << "Displacement > 0; Not implemented yet, for the calculation "
                << "of the macro tangent, slipdofs are constrained to 0"
                << endl;
        AssertThrow(false, ExcMessage("Not implemented yet."));

        deallog << "Run calculation to inital loading" << endl;
        solve();
    }
    Vector<double> initial_solution(current_solution);
    //	vector<Material<dim, nslip> > initial_quadrature_point_history =
    //			quadrature_point_history;
    vector<PerCrystalData<dim, nslip>> initial_crystal_data = crystal_data;
    SymmetricTensor<2, dim> initialmacroscopicStress = getMacroscopicStress();
    SymmetricTensor<2, dim> initialmacroscopicStrain = getMacroscopicStrain();
    Tensor<4, dim> macroscopicTangent;

    //	double pertubation = 1e-6;
    double pertubation = parameter.get_double("pertubation");
    parameter.set("displacement", pertubation);
    if (!(parameter.get("loadtype").compare("periodic") == 0
            || parameter.get("loadtype").compare("linear") == 0
            || parameter.get("loadtype").compare("linear_CG") == 0)) {
        parameter.set("loadtype", "periodic");
    }
    parameter.set("n_loadsteps", "1");

    vector<string> loadcases = { "11", "22", "33", "23", "13", "12" };
    for (unsigned int i = 0; i < loadcases.size(); i++) {
        parameter.set("straindirection", loadcases[i]);
        deallog << "Loadcase: " << parameter.get("straindirection")
                << " / displacement: " << parameter.get_double("displacement")
                << " / with" << parameter.get("loadtype") << "bc" << endl;
        current_solution = initial_solution;
        //		quadrature_point_history = initial_quadrature_point_history;
        crystal_data = initial_crystal_data;
        solve();

        //		SymmetricTensor<2, dim> macroscopic_strain;
        //		if (parameter.get("straindirection").compare("11") == 0) {
        //			macroscopic_strain[0][0] = 1;
        //		} else if (parameter.get("straindirection").compare("22") == 0) {
        //			macroscopic_strain[1][1] = 1;
        //		} else if (parameter.get("straindirection").compare("33") == 0) {
        //			macroscopic_strain[2][2] = 1;
        //		} else if (parameter.get("straindirection").compare("12") == 0) {
        //			macroscopic_strain[0][1] = 1;
        //		} else if (parameter.get("straindirection").compare("13") == 0) {
        //			macroscopic_strain[0][2] = 1;
        //		} else if (parameter.get("straindirection").compare("23") == 0) {
        //			macroscopic_strain[1][2] = 1;
        //		}
        SymmetricTensor<2, dim> macroscopic_strain = getMacroscopicStrain()
                - initialmacroscopicStrain;

        SymmetricTensor<2, dim> incrementalMacroStress = getMacroscopicStress()
                - initialmacroscopicStress;
        //		deallog << "Stress" << parameter.get("straindirection") << ": "
        //				<< incrementalMacroStress << endl;
        //		deallog << "Strain" << parameter.get("straindirection") << ": "
        //				<< macroscopic_strain << endl;
        for (int j = 0; j < dim; j++) {
            for (int k = 0; k < dim; k++) {
                for (int l = 0; l < dim; l++) {
                    for (int m = 0; m < dim; m++) {
                        if (abs(macroscopic_strain[l][m])
                                > pertubation * 1e-4) {
                            if (l == m) {
                                macroscopicTangent[j][k][l][m] +=
                                    incrementalMacroStress[j][k]
                                    / macroscopic_strain[l][m];
                            } else {
                                macroscopicTangent[j][k][l][m] += 0.5
                                                                  * incrementalMacroStress[j][k]
                                                                  / macroscopic_strain[l][m];
                            }

                        }

                    }
                }
            }
        }
    }

    //	macroscopicTangent[1][1][0][0] = macroscopicTangent[0][0][0][0]-2*macroscopicTangent[0][1][0][1];
    //	macroscopicTangent[0][0][1][1] = macroscopicTangent[0][0][0][0]-2*macroscopicTangent[0][1][0][1];
    //	macroscopicTangent[0][1][0][1] = (macroscopicTangent[0][0][0][0]
    //			- macroscopicTangent[0][1][0][1])/2;

    Tensor<4, dim> projectedMacroscopicTangent = contract4(
                P_transversal_isotropic, macroscopicTangent);

    double C11 = macroscopicTangent[0][0][0][0];
    double C33 = macroscopicTangent[2][2][2][2];
    double C44 = macroscopicTangent[1][2][1][2];
    double C66 = macroscopicTangent[0][1][0][1];
    double C12 = macroscopicTangent[1][1][0][0];
    double C13 = macroscopicTangent[2][2][0][0];

    double E_Z = C33 - 2 * C13 * C13 / (C11 + C12);
    double E_X = (C11 - C12) * (C11 * C33 + C12 * C33 - 2 * C13 * C13)
                 / (C11 * C33 - C13 * C13);
    double G_xy = C44;
    double G_xz = C66;
    double nu_xz = C13 / (C11 + C12);
    double nu_yz = (C13 * C13 - C12 * C33) / (C13 * C13 - C11 * C33);

    deallog << "E_perpendicular = " << E_X << ", E_parallel = " << E_Z
            << ", G_perpendicular = " << G_xy << ", G_parallel = " << G_xz
            << ", nu_perpendicular = " << nu_yz << ", nu_parallel = " << nu_xz
            << endl;

    print(deallog, macroscopicTangent);

    C11 = projectedMacroscopicTangent[0][0][0][0];
    C33 = projectedMacroscopicTangent[2][2][2][2];
    C44 = projectedMacroscopicTangent[1][2][1][2];
    C66 = projectedMacroscopicTangent[0][1][0][1];
    C12 = projectedMacroscopicTangent[1][1][0][0];
    C13 = projectedMacroscopicTangent[2][2][0][0];

    E_Z = C33 - 2 * C13 * C13 / (C11 + C12);
    E_X = (C11 - C12) * (C11 * C33 + C12 * C33 - 2 * C13 * C13)
          / (C11 * C33 - C13 * C13);
    G_xy = C44;
    G_xz = C66;
    nu_xz = C13 / (C11 + C12);
    nu_yz = (C13 * C13 - C12 * C33) / (C13 * C13 - C11 * C33);

    deallog << "E_perpendicular = " << E_X << ", E_parallel = " << E_Z
            << ", G_perpendicular = " << G_xy << ", G_parallel = " << G_xz
            << ", nu_perpendicular = " << nu_yz << ", nu_parallel = " << nu_xz
            << endl;

    print(deallog, projectedMacroscopicTangent);
    {
        std::string filename;
        filename = "Results/Tangent/Tangent_periodic.txt";
        std::ofstream File(filename, std::ios::out | std::ios::app);
        if (File.is_open()) {
            //		File << E_Z << "," << E_X << "," << G_xy << "," << nu_xz << "," << nu_yz
            //				<< ";" << endl;
            print(File, macroscopicTangent);
        }
    }
    {
        std::string filename;
        filename = "Results/Tangent/Tangent_projected_periodic.txt";
        std::ofstream File(filename, std::ios::out | std::ios::app);
        if (File.is_open()) {
            //		File << E_Z << "," << E_X << "," << G_xy << "," << nu_xz << "," << nu_yz
            //				<< ";" << endl;
            print(File, projectedMacroscopicTangent);
        }
    }
}

template<int dim, int nslip>
void ElastoplasticProblem<dim, nslip>::calculateYieldSurface() {

    string loadtype = parameter.get("loadtype");
    string straindirection = parameter.get("straindirection");

    /*
     * read in macroStiffnessTensor and invert it
     * calculate equivalent strain rate for loadings as
     * average of 11 and 33 uniaxial loadcases
     */
    SymmetricTensor<4, dim> macroComplianceMatrix = invert(
                getMacroTangentFromParameterfile<dim>(parameter));
    SymmetricTensor<2, dim> sig_perpendicular, sig_parallel;
    sig_perpendicular[0][0] = sig_parallel[2][2] = 1;

    SymmetricTensor<2, dim> eps_perpendicular = macroComplianceMatrix
            * sig_perpendicular;
    SymmetricTensor<2, dim> eps_parallel = macroComplianceMatrix * sig_parallel;
    eps_perpendicular /= eps_perpendicular[0][0];
    eps_parallel /= eps_parallel[2][2];
    double eps_equivalent_soll = (deviator(eps_perpendicular).norm()
                                  + deviator(eps_parallel).norm()) / 2;

    /*
     * Loop over points in quadrant and over quadrants
     * each quadrant contains its starting point but not the end point
     * phi = [0; 90[ + j * 90
     */
    int n_point = 0; // counter for points on ys
    vector<int> n_points_quadrants(2);
    n_points_quadrants[0] = (parameter.get_integer("nPointsYS") - 1) / 2.0;
    n_points_quadrants[1] = (parameter.get_integer("nPointsYS") - 1) / 2.0;

    // loop over quadrants
    for (unsigned int j = 0; j < n_points_quadrants.size(); ++j) {
        // angle in bogenmass [0; 90[ with nPointsYS
        double shift_in_quadrant = j * 3.141592653589793238462643383279 / 2.0;
        int n_points_in_quadrant = n_points_quadrants[j];

        // loop over all points in a quadrant, in addition a uniaxial loading in the beginning
        for (double i = -1.0 / n_points_in_quadrant; i < 1 - 1e-8;
                i += 1.0 / n_points_in_quadrant) {
            if (n_point < parameter.get_integer("startPointYS")
                    || n_point > parameter.get_integer("endPointYS")) {
                ++n_point;
                continue;
            }

            ++n_point;
            double angle = shift_in_quadrant
                           + i * 90 / 180.0 * 3.141592653589793238462643383279;

            /*
             * calculate macroscopic strain for BC
             */
            SymmetricTensor<2, dim> sigN = getUnitLoad<dim>(
                                               parameter.get("planeStressDir1")) * cos(angle)
                                           + getUnitLoad<dim>(parameter.get("planeStressDir2"))
                                           * sin(angle);
            // norm of deviatoric part of the resulting strain, is used to calculate the corrected sigma
            double epsEquivN = deviator(macroComplianceMatrix * sigN).norm();
            SymmetricTensor<2, dim> sigCorrectedN = sigN
                                                    / (epsEquivN / eps_equivalent_soll);
            SymmetricTensor<2, dim> strain = macroComplianceMatrix
                                             * sigCorrectedN;

            /*
             * check if first (uniaxial loading) or other
             * run (detailed with calculated stress)
             * set necessary BC
             */
            if (i < 0 - 1e-8) {
                // first point in quadrant is uniaxial loading
                if (loadtype == "periodic_plane_stress") {
                    parameter.set("loadtype", "periodic_uniaxial_tension");
                } else if (loadtype == "linear_plane_stress") {
                    parameter.set("loadtype", "linear_uniaxial_tension");
                } else if (loadtype == "linear_CG_plane_stress") {
                    parameter.set("loadtype", "linear_CG_uniaxial_tension");
                }

                if (j % 2 == 0) { // first quadrant
                    // uniaxial tension in first direction
                    parameter.set("straindirection", "11");
                    deallog << "Angle: uniaxial 11" << endl;
                } else if (j % 2 == 1) { // second quadrant
                    // no uniaxial loading
                    --n_point;
                    continue;
                }

            } else {
                // detailed strain in direction strain
                parameter.set("straindirection", "detailed");
                parameter.set("loadtype", loadtype);

                BoundaryValuesLinear<dim> *boundLinPtr = (BoundaryValuesLinear<
                        dim>*) boundary;
                boundLinPtr->setNormalizedMacroscopicStrain(strain);

                deallog << "Angle: "
                        << angle * 180 / 3.141592653589793238462643383279
                        << endl;
            }
            calculateYieldPoint();
        }
    }
}

template<int dim, int nslip>
void ElastoplasticProblem<dim, nslip>::calculateYieldPoint() {
    // reinit solutionstate
    current_solution = 0;
    prior_solution = 0;
    current_solution_abs = 0; // absolute solutions for hardening
    prior_solution_abs = 0; // absolute solutions for hardening
    current_hardening = 0; // current hardening in slipsystems
    prior_hardening = 0; // prior hardening in slipsystems

    setup_quadrature_point_history();
    reCreateSparsityPattern = true;

    // solve with parameters in parameterhandler
    solve();
}

template<int dim, int nslip>
void ElastoplasticProblem<dim, nslip>::run() {

    createMesh();

    setup_system();
    setup_quadrature_point_history();
    identify_grainboundarys();
    print_detailed_crystal_informations();
    timer.enter_section("Post init");
    post.init();
    timer.exit_section("Post init");

    if (parameter.get("mode").compare("solve") == 0) {
        deallog << "Solve" << endl;
        solve();
    } else if (parameter.get("mode").compare("calculateTangent") == 0) {
        deallog << "Calculate tangent" << endl;
        calculateMacroscopicTangent();
    } else if (parameter.get("mode").compare("calculateYieldSurface") == 0) {
        deallog << "Calculate yield surface" << endl;
        calculateYieldSurface();
    }
    parameter.print_parameters(deallog.get_file_stream(),
                               ParameterHandler::OutputStyle::Text);

}

template<int dim, int nslip>
int ElastoplasticProblem<dim, nslip>::newtonCorrection(
    Vector<double> &newtonStep) {

    if (parameter.get("newtonCorrection") == "none") {
        return 0;
    }

    int n_components = dof_handler.get_fe().n_components();
    std::vector<bool> mask(n_components, true);
    for (int i = 0; i < dim; ++i) {
        mask[i] = false;
    }
    ComponentMask compMask = ComponentMask(mask);

    std::vector<bool> slip_dofs(dof_handler.n_dofs());
    DoFTools::extract_dofs(dof_handler, compMask, slip_dofs);

    PerCrystalData<dim, nslip> *perCrystalData = &crystal_data[0];
    double dt = boundary->get_dt();
    unsigned int n = 0;

    for (unsigned int i = 0; i < dof_handler.n_dofs(); ++i) {
        if (!slip_dofs[i]) {
            continue;
        }

        double threshold = perCrystalData->get_gamma0dot() * dt
                           *1.088659492482653;                           
        int old_section;
        int new_section;

        // check in which section the old solution was
        if (current_solution[i] - prior_solution[i] < -threshold) {
            old_section = -1;
        } else if (current_solution[i] - prior_solution[i] > threshold) {
            old_section = 1;
        } else {
            old_section = 0;
        }

        // check in which section the new solution is
        if (current_solution[i] + newtonStep[i] - prior_solution[i]
                < -threshold) {
            new_section = -1;
        } else if (current_solution[i] + newtonStep[i] - prior_solution[i]
                   > threshold) {
            new_section = 1;
        } else {
            new_section = 0;
        }

        if (parameter.get("newtonCorrection") == "naive") {
            if (abs(old_section - new_section) == 2) {
                newtonStep[i] = -(current_solution[i] - prior_solution[i]);
                ++n;
            }
        }

        /*
         * TODO
         * calculate pi_1 at (current_solution - prior_solution)
         * calculate tangent at (current_solution - prior_solution)
         * calculate pi_2 at (current_solution + newtonStep - prior_solution) using the tangent at (current_solution - prior_solution)
         * calculate calculate dGamma at pi_2
         * calculate decide which dGamma to use and change the newtonStep
         */

        if (parameter.get("newtonCorrection") == "improved"
                && abs(old_section - new_section) == 2) {
            Vector<double> localDGamma(1);
            localDGamma[0] = current_solution[i] - prior_solution[i];
            Vector<double> localHardening(1);
            localHardening[0] = current_hardening[i];
            Vector<double> localGammaAbs(1);
            localGammaAbs[0] = prior_solution_abs[i] + abs(localDGamma[0]);

//			calculate pi_1 at (current_solution - prior_solution)
            double pi_1 = GradientCP::get_PiAlpha(0, dt, 23, localDGamma,
                                                  localHardening, perCrystalData);

//			calculate tangent at (current_solution - prior_solution)
            double dPidGamma = GradientCP::getDpiAlphaDgammaBeta(0, 0, dt, 23,
                               localDGamma, localGammaAbs, localHardening, perCrystalData);

            Vector<double> localDGamma_2(1);
            localDGamma_2[0] = current_solution[i] + newtonStep[i]
                               - prior_solution[i];
//			calculate pi_2 at (current_solution + newtonStep - prior_solution)
//			using the tangent at (current_solution - prior_solution)
            double pi_2 = pi_1 + dPidGamma * newtonStep[i];

//			calculate dGamma for pi_2
            double dGamma_2 = GradientCP::get_gammaDot(dt, 23, pi_2,
                              localHardening[0], perCrystalData);

            double delta_2 = dGamma_2 - current_solution[i] + prior_solution[i];

//			decide which dGamma to use and change the newtonStep
            if (abs(delta_2) < newtonStep[i]) {
                if (abs(delta_2) > 0) {
                    newtonStep[i] = delta_2;
                    ++n;
                }
            }
        }
    }
    constraints.distribute(newtonStep);
    return n;
}

template<int dim, int nslip>
void ElastoplasticProblem<dim, nslip>::update_hardening() {

// update current_solution_abs
    for (unsigned int i = 0; i < current_solution_abs.size(); ++i) {
        current_solution_abs(i) = prior_solution_abs(i)
                                  + abs(current_solution(i) - prior_solution(i));
    }

    std::vector<unsigned int> local_dof_indices;
    Vector<double> prior_gamma_cell;
    Vector<double> current_gamma_cell;
    Vector<double> dGamma_cell;
    Vector<double> current_gamma_abs_cell;
    Vector<double> current_hardening_cell;
    Vector<double> prior_hardening_cell;
    typename hp::DoFHandler<dim>::active_cell_iterator cell =
        dof_handler.begin_active(), endc = dof_handler.end();

    for (; cell != endc; ++cell) {
        const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;
        PerCrystalData<dim, nslip> *perCrystalData =
            &crystal_data[cell->active_fe_index()];

        local_dof_indices.resize(dofs_per_cell);
        cell->get_dof_indices(local_dof_indices);

        prior_gamma_cell.reinit(dofs_per_cell);
        current_gamma_cell.reinit(dofs_per_cell);
        dGamma_cell.reinit(dofs_per_cell);
        current_gamma_abs_cell.reinit(dofs_per_cell);
        current_hardening_cell.reinit(dofs_per_cell);
        prior_hardening_cell.reinit(dofs_per_cell);

        cell->get_dof_values(prior_solution, prior_gamma_cell);
        cell->get_dof_values(current_solution, current_gamma_cell);
        cell->get_dof_values(current_solution_abs, current_gamma_abs_cell);
        cell->get_dof_values(prior_hardening, prior_hardening_cell);

        for (unsigned int i = 0; i < dGamma_cell.size(); ++i) {
            dGamma_cell[i] = current_gamma_cell[i] - prior_gamma_cell[i];
        }

        // loop over all dofs
        for (unsigned int i = 0; i < dofs_per_cell;) {
            const unsigned int component =
                compToCPSystem[fe_collection[cell->active_fe_index()].system_to_component_index(
                                   i).first];
            // find first slip-dof of the next node
            if (component >= dim) {

                // generate node_value_vectors by copying from the cell-ones
                Vector<double> dGamma_node(dGamma_cell.begin() + i,
                                           dGamma_cell.begin() + i + nslip);
                Vector<double> current_gamma_abs_node(
                    current_gamma_abs_cell.begin() + i,
                    current_gamma_abs_cell.begin() + i + nslip);
                Vector<double> prior_hardening_node(
                    prior_hardening_cell.begin() + i,
                    prior_hardening_cell.begin() + i + nslip);

                // calculate hardening for all slipdofs on that node
                for (unsigned int j = 0; j < nslip; ++j) {
                    const unsigned int slipsystem =
                        compToCPSystem[fe_collection[cell->active_fe_index()].system_to_component_index(
                                           i + j).first] - dim;
                    current_hardening_cell[i + j] = GradientCP::getGAlpha(
                                                        slipsystem, 23, dGamma_node, current_gamma_abs_node,
                                                        prior_hardening_node, perCrystalData);
                }
                i += nslip;
            } else {
                ++i;
            }
        }
        //		if (current_hardening_cell != current_hardening_cell) {
        //			deallog << endl << current_hardening_cell << endl;
        //		}
        cell->set_dof_values(current_hardening_cell, current_hardening);
    }
}

//template<int dim, int nslip>
//void ElastoplasticProblem<dim, nslip>::quadrature_point_history_copy_newtonstep_to_history() {
//	hp::FEValues<dim> hp_fe_values(fe_collection, quadrature_collection,
//			update_default);
//	typename hp::DoFHandler<dim>::active_cell_iterator cell =
//			dof_handler.begin_active(), endc = dof_handler.end();
//	for (; cell != endc; ++cell) {
//		Material<dim, nslip> *local_quadrature_points_history =
//				reinterpret_cast<Material<dim, nslip> *>(cell->user_pointer());
//		Assert(
//				local_quadrature_points_history >= &quadrature_point_history.front(),
//				ExcInternalError());
//		Assert(
//				local_quadrature_points_history <= &quadrature_point_history.back(),
//				ExcInternalError());
//
//		hp_fe_values.reinit(cell);
//		const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values();
//		for (unsigned int q = 0; q < fe_values.n_quadrature_points; ++q) {
//			local_quadrature_points_history[q].update_pointhistorydata();
//		}
//	}
//}

// setup just der per Crystal data vector
template<int dim, int nslip>
void ElastoplasticProblem<dim, nslip>::setup_quadrature_point_history() {
    timer.enter_section("Setup QuadraturPointHistory");

// clear user pointer, gauspoint history and cell data
    triangulation.clear_user_data();
    {
        std::vector<PerCrystalData<dim, nslip> > tmp;
        tmp.swap(crystal_data);	// clear Per Crystal data
    }
    crystal_data.resize(256);
    typename hp::DoFHandler<dim>::active_cell_iterator cell =
        dof_handler.begin_active(), endc = dof_handler.end();

    hp::FEValues<dim> hp_fe_values(fe_collection, quadrature_collection,
                                   update_quadrature_points);
    /* initialize random seed: */
    srand(time(NULL));
    std::map<int, std::vector<double>> angles;
    for (; cell != endc; ++cell) {
        int domainID = cell->material_id();

        if (angles.count(domainID) == 0) {

            /*
             * x=0, y=0 -> ID 0-3
             * x=0, y=1 -> ID 4-7
             * -> domainID % 4 gives the layer
             */
            int csFaktor = 1;
            if (parameter.get_bool("crosssnake")) {
                csFaktor = 4;
            }
            vector<double> eulerangles = Meshgeneration::calculateEulerAngles(
                                             domainID / csFaktor, parameter);
            for (int i = 0; i < csFaktor; ++i) {
                // domainID - domainID%csFaktor gives the ID of the crystals bottom
                angles[domainID - domainID%csFaktor + i] = eulerangles;
                crystal_data[domainID - domainID%csFaktor + i].init(parameter, eulerangles[0],
                        eulerangles[1], eulerangles[2], domainID - domainID%csFaktor + i,
                        (domainID - domainID%csFaktor + i) % csFaktor);
            }
        }

//		//initialize GP data
//		hp_fe_values.reinit(cell);
//		const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values();
//		for (unsigned int i = 0; i < fe_values.n_quadrature_points; i++) {
//		}
    }

// Erase empty PerCrystalData
    for (typename vector<PerCrystalData<dim, nslip> >::iterator it =
                crystal_data.begin(); it != crystal_data.end();) {
        if (!it->is_initiated()) {
            crystal_data.erase(it);
        } else {
            it++;
        }
    }

//	Assert(history_index == quadrature_point_history.size(), ExcInternalError());
    timer.exit_section("Setup QuadraturPointHistory");
}

// Setup quadrature point history and per crystal data vector
//template<int dim, int nslip>
//void ElastoplasticProblem<dim, nslip>::setup_quadrature_point_history() {
//	timer.enter_section("Setup QuadraturPointHistory");
//
//	// clear user pointer, gauspoint history and cell data
//	triangulation.clear_user_data();
//	{
//		std::vector<Material<dim, nslip> > tmp;
//		tmp.swap(quadrature_point_history);	// clear quadrature point history
//		std::vector<PerCrystalData<dim, nslip> > tmp2;
//		tmp2.swap(crystal_data);	// clear Per Crystal data
//	}
//	quadrature_point_history.resize(
//			triangulation.n_active_cells()
//					* quadrature_collection.max_n_quadrature_points());
//	crystal_data.resize(256);
//	unsigned int history_index = 0;
//	typename hp::DoFHandler<dim>::active_cell_iterator cell =
//			dof_handler.begin_active(), endc = dof_handler.end();
//
//	hp::FEValues<dim> hp_fe_values(fe_collection, quadrature_collection,
//			update_quadrature_points);
//	/* initialize random seed: */
//	srand(time(NULL));
//	std::map<int, std::vector<double>> angles;
//	for (; cell != endc; ++cell) {
//		cell->set_user_pointer(&quadrature_point_history[history_index]);
//		int domainID = cell->material_id();
//
//		if (angles.count(domainID) == 0) {
//
//			/*
//			 * x=0, y=0 -> ID 0-3
//			 * x=0, y=1 -> ID 4-7
//			 * -> domainID % 4 gives the layer
//			 */
//			int csFaktor = 1;
//			if (parameter.get_bool("crosssnake")) {
//				csFaktor = 4;
//			}
//			vector<double> eulerangles = Meshgeneration::calculateEulerAngles(
//					domainID / csFaktor, parameter);
//			angles[domainID] = eulerangles;
//			crystal_data[domainID].init(parameter, eulerangles[0],
//					eulerangles[1], eulerangles[2], domainID,
//					domainID % csFaktor);
//		}
//
//		//initialize GP data
//		hp_fe_values.reinit(cell);
//		const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values();
//		for (unsigned int i = 0; i < fe_values.n_quadrature_points; i++) {
//			quadrature_point_history[history_index + i].init(parameter,
//					&crystal_data[domainID],
//					boundary->get_initial_theta(fe_values.quadrature_point(i)));
//		}
//		//end set material ID
//		history_index += fe_values.n_quadrature_points;
//	}
//
//	// Erase empty PerCrystalData
//	for (typename vector<PerCrystalData<dim, nslip> >::iterator it =
//			crystal_data.begin(); it != crystal_data.end();) {
//		if (!it->is_initiated()) {
//			crystal_data.erase(it);
//		} else {
//			it++;
//		}
//	}
//
//	Assert(history_index == quadrature_point_history.size(), ExcInternalError());
//	timer.exit_section("Setup QuadraturPointHistory");
//}

template<int dim, int nslip>
void ElastoplasticProblem<dim, nslip>::postprocessing(unsigned int loadstep) {
    timer.enter_section("Postprocessing");

    if (parameter.get("postprocessing").find("full,") != string::npos) {
        deallog << "Full Postprocessing" << endl;
        post.writeResults(loadstep, boundary->get_dt(), current_solution,
                          prior_solution, current_solution_abs);
        post.save_double_slip_diagrams(extractor_slips, current_solution);
    }
    {
        // write stress strain to file
        std::stringstream sstm;
        sstm << "Results/stress_strain_" /*<< parameter.get("straindirection")*/
			 << parameter.get("regularisation")
			 <<"_"<<parameter.get("gamma0dot")
			 <<"_"<<parameter.get("p")
			 <<"_"<<parameter.get("strainrate")
             << ".txt";
        string filename = sstm.str();

        std::ofstream File(filename, std::ios::out | std::ios::app);
        if (File.is_open()) {
            File << "Begin step " << loadstep << endl;
            File << "Stress " << getMacroscopicStress() << endl;
            File << "Strain " << boundary->get_full_macroscopic_strain()
                 << endl;
            File << "Strain mean " << getMacroscopicStrain() << endl;
            File << "Plastic strain_mean " << getMacroscopicPlasticStrain()
                 << endl;
            File << "Plastic strain_rate_mean "
                 << getMacroscopicPlasticStrainRate() << endl;
            File << "Dissipation " << calculate_dissipation() << endl;
            File << "End step " << loadstep << endl;
        }
    }
    if (parameter.get_bool("writeOutVectors")) {
        Vector<double> dofType(dof_handler.n_dofs());

        for (unsigned int i = 0; i < dof_handler.n_dofs(); i++) {
            if (dofdata.count(i) != 0) { // u_dof
                dofType[i] = 0;
            } else { // gammadof
                dofType[i] = 1;
            }
        }

        {
            string filename;
            std::stringstream sstm;
            sstm << "Results/dofType.txt";
            filename = sstm.str();

            std::ofstream File(filename.c_str(), std::ios::out | std::ios::app);
            dofType.print(File);
        }

        {
            string filename;
            std::stringstream sstm;
            sstm << "Results/Solution.txt";
            filename = sstm.str();
            std::ofstream File(filename.c_str(), std::ios::out | std::ios::app);
            current_solution.print(File);
        }
    }
    timer.exit_section("Postprocessing");
}

template<int dim, int nslip>
SymmetricTensor<2, dim> ElastoplasticProblem<dim, nslip>::getMacroscopicStress() {
//	hp::FEValues<dim> hp_fe_values(fe_collection, quadrature_collection,
//			update_JxW_values);
//	SymmetricTensor<2, dim> macroscopicStress;
//	double volume = 0;
//	typename Triangulation<dim>::active_cell_iterator cell =
//			triangulation.begin_active(), endc = triangulation.end();
//	for (; cell != endc; ++cell) {
//		Material<dim, nslip> *local_quadrature_points_data =
//				reinterpret_cast<Material<dim, nslip>*>(cell->user_pointer());
//		hp_fe_values.reinit(cell);
//		const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values();
//		for (unsigned int q = 0; q < fe_values.n_quadrature_points; ++q) {
//			macroscopicStress += local_quadrature_points_data[q].get_stress()
//					* fe_values.JxW(q);
//			volume += fe_values.JxW(q);
//		}
//	}
//	return macroscopicStress / volume;

    hp::FEValues<dim> hp_fe_values(fe_collection, quadrature_collection,
                                   update_values | update_gradients | update_JxW_values);
    SymmetricTensor<2, dim> macroscopicStress;
    double volume = 0;

    typename hp::DoFHandler<dim>::active_cell_iterator cell =
        dof_handler.begin_active(), endc = dof_handler.end();
    for (; cell != endc; ++cell) {

        hp_fe_values.reinit(cell);
        const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values();

        // get strains and slips from solutionvector
        const unsigned int n_q_points = fe_values.n_quadrature_points;
        vector<SymmetricTensor<2, dim> > strains(n_q_points);
        fe_values[extractor_displacement].get_function_symmetric_gradients(
            current_solution, strains);
        vector<vector<double> > slips(nslip, vector<double>(n_q_points));
        for (unsigned int i = 0;
                i < extractor_slips[cell->active_fe_index()].size(); i++) {
            fe_values[extractor_slips[cell->active_fe_index()][i]].get_function_values(
                current_solution, slips[i]);
        }
        Vector<double> gamma(slips.size());

        // Loop over all qp and integrate volume and stress
        for (unsigned int q = 0; q < n_q_points; ++q) {
            //			PerCrystalData<dim, nslip> *perCrystalData =
            //					local_quadrature_points_data[q].get_per_cell_data();
            PerCrystalData<dim, nslip> *perCrystalData =
                &crystal_data[cell->active_fe_index()];
            for (unsigned int i = 0; i < gamma.size(); i++) {
                gamma[i] = slips[i][q];
            }
            SymmetricTensor<2, dim> stress;

            if (perCrystalData->isLocalFormulation()) {

            } else {
                macroscopicStress += GradientCP::get_stress(strains[q], gamma, /*Temperature*/
                                     23, perCrystalData) * fe_values.JxW(q);
            }
            volume += fe_values.JxW(q);
        }
    }
//	deallog << "stress based volume average: " << endl;
//	print(deallog.get_console(), Tensor<2, dim>(macroscopicStress) / volume);

// return averaged stress
    return macroscopicStress / volume;
}

template<int dim, int nslip>
SymmetricTensor<2, dim> ElastoplasticProblem<dim, nslip>::getMacroscopicStrain() {
//	hp::FEValues<dim> hp_fe_values(fe_collection, quadrature_collection,
//			update_values | update_gradients | update_JxW_values
//					| update_quadrature_points);
//	SymmetricTensor<2, dim> macroscopicStrain;
//	double volume = 0;
//	typename hp::DoFHandler<dim>::active_cell_iterator cell =
//			dof_handler.begin_active(), end_cell = dof_handler.end();
//	for (; cell != end_cell; ++cell) {
//		hp_fe_values.reinit(cell);
//		const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values();
//
//		vector<SymmetricTensor<2, dim> > displacement_grads(
//				fe_values.n_quadrature_points);
//		fe_values[extractor_displacement].get_function_symmetric_gradients(
//				current_solution, displacement_grads);
//
//		for (unsigned int q = 0; q < fe_values.n_quadrature_points; ++q) {
//			macroscopicStrain += displacement_grads[q] * fe_values.JxW(q);
//			volume += fe_values.JxW(q);
//		}
//	}
//	return macroscopicStrain / volume;

    hp::FEValues<dim> hp_fe_values(fe_collection, quadrature_collection,
                                   update_gradients | update_JxW_values);
    SymmetricTensor<2, dim> macroscopicStrain;
    double volume = 0;

    typename hp::DoFHandler<dim>::active_cell_iterator cell =
        dof_handler.begin_active(), endc = dof_handler.end();
    for (; cell != endc; ++cell) {

        hp_fe_values.reinit(cell);
        const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values();

        // get strains from solutionvector
        const unsigned int n_q_points = fe_values.n_quadrature_points;
        vector<SymmetricTensor<2, dim> > strains(n_q_points);
        fe_values[extractor_displacement].get_function_symmetric_gradients(
            current_solution, strains);

        // Loop over all qp and integrate volume and stress
        for (unsigned int q = 0; q < n_q_points; ++q) {
            macroscopicStrain += strains[q] * fe_values.JxW(q);
            volume += fe_values.JxW(q);
        }
    }

// calculate macrostrain with reference nodes
    Point<dim> u_origin;
    Point<dim> u_ref_x;
    Point<dim> u_ref_y;
    Point<dim> u_ref_z;

    for (unsigned int i = 0; i < dim; ++i) {
        u_origin[i] = current_solution[dof_reference_node[0][i]];
        u_ref_x[i] = current_solution[dof_reference_node[1][i]];
        u_ref_y[i] = current_solution[dof_reference_node[2][i]];
        u_ref_z[i] = current_solution[dof_reference_node[3][i]];
    }

    Tensor<2, dim> macro_strain_uRef;
    for (unsigned int i = 0; i < dim; ++i) {
        macro_strain_uRef[0][i] = (u_ref_x[i] - u_origin[i])
                                  / geometric_range[0];
        macro_strain_uRef[1][i] = (u_ref_y[i] - u_origin[i])
                                  / geometric_range[1];
        macro_strain_uRef[2][i] = (u_ref_z[i] - u_origin[i])
                                  / geometric_range[2];
    }

    if ((macro_strain_uRef - Tensor<2, dim>(macroscopicStrain) / volume).norm()
            > 1e-6) {
        deallog << "strain based on uRef: " << endl;
        print(deallog.get_console(), macro_strain_uRef);
        deallog << "strain based on volume average: " << endl;
        print(deallog.get_console(),
              Tensor<2, dim>(macroscopicStrain) / volume);
    }

// return averaged strain
    return macroscopicStrain / volume;
}

template<int dim, int nslip>
SymmetricTensor<2, dim> ElastoplasticProblem<dim, nslip>::getMacroscopicPlasticStrain() {
//	hp::FEValues<dim> hp_fe_values(fe_collection, quadrature_collection,
//			update_JxW_values);
//	SymmetricTensor<2, dim> macroscopicPlasticStrain;
//	double volume = 0;
//	typename Triangulation<dim>::active_cell_iterator cell =
//			triangulation.begin_active(), endc = triangulation.end();
//	for (; cell != endc; ++cell) {
//		Material<dim, nslip> *local_quadrature_points_data =
//				reinterpret_cast<Material<dim, nslip>*>(cell->user_pointer());
//		hp_fe_values.reinit(cell);
//		const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values();
//		for (unsigned int q = 0; q < fe_values.n_quadrature_points; ++q) {
//			macroscopicPlasticStrain +=
//					local_quadrature_points_data[q].get_plastic_strain()
//							* fe_values.JxW(q);
//			volume += fe_values.JxW(q);
//		}
//	}
//	return macroscopicPlasticStrain / volume;
    hp::FEValues<dim> hp_fe_values(fe_collection, quadrature_collection,
                                   update_values | update_JxW_values);
    SymmetricTensor<2, dim> macroscopicPlasticStrain;
    double volume = 0;

    typename hp::DoFHandler<dim>::active_cell_iterator cell =
        dof_handler.begin_active(), endc = dof_handler.end();
    for (; cell != endc; ++cell) {

        hp_fe_values.reinit(cell);
        const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values();

        // get slips from solutionvector
        const unsigned int n_q_points = fe_values.n_quadrature_points;
        vector<vector<double> > slips(nslip, vector<double>(n_q_points));
        for (unsigned int i = 0;
                i < extractor_slips[cell->active_fe_index()].size(); i++) {
            fe_values[extractor_slips[cell->active_fe_index()][i]].get_function_values(
                current_solution, slips[i]);
        }
        Vector<double> gamma(slips.size());

        // Loop over all qp and integrate volume and plastic strains
        for (unsigned int q = 0; q < n_q_points; ++q) {
            PerCrystalData<dim, nslip> *perCrystalData =
                &crystal_data[cell->active_fe_index()];
            for (unsigned int i = 0; i < gamma.size(); i++) {
                gamma[i] = slips[i][q];
            }

            if (perCrystalData->isLocalFormulation()) {

            } else {
                macroscopicPlasticStrain += GradientCP::get_plStrain(gamma,
                                            perCrystalData) * fe_values.JxW(q);
            }
            volume += fe_values.JxW(q);
        }
    }

// return averaged plastic strain
    return macroscopicPlasticStrain / volume;
}

template<int dim, int nslip>
SymmetricTensor<2, dim> ElastoplasticProblem<dim, nslip>::getMacroscopicPlasticStrainRate() {

//	hp::FEValues<dim> hp_fe_values(fe_collection, quadrature_collection,
//			update_JxW_values);
//	SymmetricTensor<2, dim> macroscopicPlasticStrainRate;
//	double volume = 0;
//	typename Triangulation<dim>::active_cell_iterator cell =
//			triangulation.begin_active(), endc = triangulation.end();
//	for (; cell != endc; ++cell) {
//		Material<dim, nslip> *local_quadrature_points_data =
//				reinterpret_cast<Material<dim, nslip>*>(cell->user_pointer());
//		hp_fe_values.reinit(cell);
//		const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values();
//		for (unsigned int q = 0; q < fe_values.n_quadrature_points; ++q) {
//			macroscopicPlasticStrainRate +=
//					local_quadrature_points_data[q].get_plastic_strain_rate()
//							* fe_values.JxW(q);
//			volume += fe_values.JxW(q);
//		}
//	}
//	return macroscopicPlasticStrainRate / volume;

//	hp::FEValues<dim> hp_fe_values(fe_collection, quadrature_collection,
//			update_values | update_JxW_values);
//	SymmetricTensor<2, dim> macroscopicPlasticStrainRate;
//	double volume = 0;
//
//	typename hp::DoFHandler<dim>::active_cell_iterator cell =
//			dof_handler.begin_active(), endc = dof_handler.end();
//	for (; cell != endc; ++cell) {
//
//		Material<dim, nslip> *local_quadrature_points_data =
//				reinterpret_cast<Material<dim, nslip>*>(cell->user_pointer());
//		hp_fe_values.reinit(cell);
//		const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values();
//
//		// get slips from solutionvector
//		const unsigned int n_q_points = fe_values.n_quadrature_points;
//		vector<vector<double> > slips(nslip, vector<double>(n_q_points));
//		vector<vector<double> > prior_slips(nslip, vector<double>(n_q_points));
//		for (unsigned int i = 0;
//				i < extractor_slips[cell->active_fe_index()].size(); i++) {
//			fe_values[extractor_slips[cell->active_fe_index()][i]].get_function_values(
//					current_solution, slips[i]);
//			fe_values[extractor_slips[cell->active_fe_index()][i]].get_function_values(
//					prior_solution, prior_slips[i]);
//		}
//		Vector<double> gamma(slips.size());
//		Vector<double> prior_gamma(slips.size());
//
//		// Loop over all qp and integrate volume and plastic strain rates
//		for (unsigned int q = 0; q < n_q_points; ++q) {
//			PerCrystalData<dim, nslip> *perCrystalData =
//					local_quadrature_points_data[q].get_per_cell_data();
//			for (unsigned int i = 0; i < gamma.size(); i++) {
//				gamma[i] = slips[i][q];
//				prior_gamma[i] = prior_slips[i][q];
//			}
//			SymmetricTensor<2, dim> priorPlasticStrain =
//					local_quadrature_points_data[q].get_plStrain(prior_gamma,
//							perCrystalData) * fe_values.JxW(q);
//			;
//			SymmetricTensor<2, dim> plasticStrain =
//					local_quadrature_points_data[q].get_plStrain(gamma,
//							perCrystalData) * fe_values.JxW(q);
//			;
//			macroscopicPlasticStrainRate += (plasticStrain - priorPlasticStrain)
//					/ boundary->get_dt();
//			volume += fe_values.JxW(q);
//		}
//	}
//
//	// return averaged plastic strain rate
//	return macroscopicPlasticStrainRate / volume;
    return SymmetricTensor<2, dim>();
}

// dissipationdensity
template<int dim, int nslip>
double ElastoplasticProblem<dim, nslip>::calculate_dissipation() {
    double d_dissipation = 0;
    if (boundary->get_time() == 0) {
        return d_dissipation;
    }
    double volume = 0;
    hp::FEValues<dim> hp_fe_values(fe_collection, quadrature_collection,
                                   update_values | update_JxW_values);
//	double volume = 0;
    typename hp::DoFHandler<dim>::active_cell_iterator cell =
        dof_handler.begin_active(), endc = dof_handler.end();
    for (; cell != endc; ++cell) {

        PerCrystalData<dim, nslip> *perCrystalData =
            &crystal_data[cell->active_fe_index()];

        hp_fe_values.reinit(cell);
        const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values();

        vector<vector<double> > slips(nslip,
                                      vector<double>(fe_values.n_quadrature_points));
        vector<vector<double> > prior_slips(nslip,
                                            vector<double>(fe_values.n_quadrature_points));
        vector<vector<double> > hardening_cell(nslip,
                                               vector<double>(fe_values.n_quadrature_points));
        for (unsigned int i = 0;
                i < extractor_slips[cell->active_fe_index()].size(); i++) {
            fe_values[extractor_slips[cell->active_fe_index()][i]].get_function_values(
                current_solution, slips[i]);
            fe_values[extractor_slips[cell->active_fe_index()][i]].get_function_values(
                prior_solution, prior_slips[i]);
            fe_values[extractor_slips[cell->active_fe_index()][i]].get_function_values(
                current_hardening, hardening_cell[i]);
        }

        for (unsigned int q = 0; q < fe_values.n_quadrature_points; ++q) {
            Vector<double> dGamma(slips.size());
            Vector<double> hardening_gp(slips.size());
            for (unsigned int i = 0; i < dGamma.size(); i++) {
                dGamma[i] = slips[i][q] - prior_slips[i][q];
                hardening_gp[i] = hardening_cell[i][q];
            }
            for (unsigned int alpha = 0; alpha < nslip; ++alpha) {
                double tmp = GradientCP::get_PiAlpha(alpha, boundary->get_dt(), /*temperature*/
                                                     23, dGamma, hardening_gp, perCrystalData) * dGamma[alpha]
                             * fe_values.JxW(q);
                if (tmp < 0) {
                    deallog << "negative dissipation" << endl;
                    exit(0);
                }
                d_dissipation += tmp;
                volume += fe_values.JxW(q);
            }
        }
    }
    return d_dissipation / volume;
}

//template<int dim, int nslip>
//double ElastoplasticProblem<dim, nslip>::maxValueEvolutionEquation() {
//	double retval = 0;
//	for (unsigned int i = 0; i < crystal_data.size(); i++) {
//		if (crystal_data[i].get_possible_convergence_problems() > retval) {
//			retval = crystal_data[i].get_possible_convergence_problems();
//		}
//	}
//	return retval;
//}

template<int dim, int nslip>
vector<double> ElastoplasticProblem<dim, nslip>::calculateResiduals(
    Vector<double> &res) {
    vector<double> retval(2);

    for (unsigned int i = 0; i < dof_handler.n_dofs(); i++) {
        if (dofdata.count(i) != 0) { // u_dof
            retval[0] += pow(res[i], 2);
        } else { // gammadof
            retval[1] += pow(res[i], 2);
        }
    }
    retval[0] = sqrt(retval[0]);
    retval[1] = sqrt(retval[1]);
    return retval;
}

template<int dim, int nslip>
void ElastoplasticProblem<dim, nslip>::print_conv_header() {
    static const unsigned int l_width = 170;
    static const unsigned int c_width = 12;
    for (unsigned int i = 0; i < l_width; ++i)
        deallog << "_";
    deallog << std::endl;
    deallog << setw(39) << "SOLVER STEP " << "| " << setw(c_width) << "RES_NORM"
            << setw(c_width) << "RES_U" << setw(c_width) << "RES_G"
            << setw(c_width) << "D_U" << setw(c_width) << "D_G" << setw(c_width)
            << "LIN_IT" << setw(c_width) << "COND_NR" << setw(c_width)
            /*<< "Relax_F"*/<< setw(c_width) << "Regular_F" << setw(c_width)
            << "ConvergF" << setw(c_width) << "regIncF" << setw(c_width)
            << "NewtonCor" << std::endl;
    for (unsigned int i = 0; i < l_width; ++i)
        deallog << "_";
    deallog << std::endl;
}

template<int dim, int nslip>
void ElastoplasticProblem<dim, nslip>::print_conv_footer() {
    static const unsigned int l_width = 170;
    deallog << std::endl;
    for (unsigned int i = 0; i < l_width; ++i)
        deallog << "_";
    deallog << std::endl;
}

template<int dim, int nslip>
void ElastoplasticProblem<dim, nslip>::slotConditionNumber(double d) {
    conditionNumber = d;
}

}		//End Namespace
int main(int argc, char *argv[]) {

    std::cerr << "Test std::cerr" << endl;
    std::cout << "Test std::out" << endl;
    if (system("mkdir -p Results") + system("mkdir -p Results/Tangent") != 0) {
        cout << "Error when creating directorys Result" << endl;
    }

    std::string path;
    if (argc > 1) {
        if (std::strcmp("path", argv[1])) {
            std::cout << "Path: " << argv[2] << std::endl;
            path = argv[2];
        }
    } else {
        //		std::cout << argv[0] << std::endl;
    }

    using namespace dealii;
	ParameterHandler parameter;
	read_Parameterfile("Input/parameter", parameter);
    deallog.depth_console(1);
    deallog.depth_file(1);
	
	std::string test(parameter.get("regularisation")),test2(parameter.get("gamma0dot")),test3(parameter.get("p")),
	test4("_"), test5(parameter.get("strainrate"));
	std::string details = test4 + test + test4 + test2 + test4 + test3 + test4 + test5;
	std::cout<<"!!!!!!!!!!!!!!"<<details<<endl;
    std::ofstream logfile("Results/deallog"+details);/*+test+test4+test2+test4+test3);*/
    deallog.attach(logfile);

    //ParameterHandler parameter;
    //read_Parameterfile("Input/parameter", parameter);

//	int feq_u = 0;
//	if (parameter.get("el_type_u").compare("Q1") == 0) {
//		feq_u = 1;
//	} else if (parameter.get("el_type_u").compare("Q2") == 0) {
//		feq_u = 2;
//	} else if (parameter.get("el_type_u").compare("Q3") == 0) {
//		feq_u = 3;
//	}
//	int feq_gamma = 0;
//	if (parameter.get("el_type_gamma").compare("Q1") == 0) {
//		feq_gamma = 1;
//	} else if (parameter.get("el_type_gamma").compare("Q2") == 0) {
//		feq_gamma = 2;
//	} else if (parameter.get("el_type_gamma").compare("Q3") == 0) {
//		feq_gamma = 3;
//	}
    try {
        using namespace dealii;
        if (parameter.get_integer("n_slip") == 0) {
        	SFB814_C5::ElastoplasticProblem<3, 0> elastic_problem_3d(
        			parameter.get_integer("number_gp"));
        	elastic_problem_3d.run();
        } else if (parameter.get_integer("n_slip") == 1) {
        	SFB814_C5::ElastoplasticProblem<3, 1> elastic_problem_3d(
        			parameter.get_integer("number_gp"));
        	elastic_problem_3d.run();
        }else
        if (parameter.get_integer("n_slip") == 2) {
            SFB814_C5::ElastoplasticProblem<3, 2> elastic_problem_3d(
                parameter.get_integer("number_gp"));
            elastic_problem_3d.run();
        } else if (parameter.get_integer("n_slip") == 12) {
			SFB814_C5::ElastoplasticProblem<3, 12> elastic_problem_3d(
					parameter.get_integer("number_gp"));
			elastic_problem_3d.run();
		}  
    } catch (std::exception &exc) {
        std::cerr << std::endl << std::endl
                  << "----------------------------------------------------"
                  << std::endl;
        std::cerr << "Exception on processing: " << std::endl << exc.what()
                  << std::endl << "Aborting!" << std::endl
                  << "----------------------------------------------------"
                  << std::endl;

        return 1;
    } catch (...) {
        std::cerr << std::endl << std::endl
                  << "----------------------------------------------------"
                  << std::endl;
        std::cerr << "Unknown exception!" << std::endl << "Aborting!"
                  << std::endl
                  << "----------------------------------------------------"
                  << std::endl;
        return 2;
    }

    return 0;
}





