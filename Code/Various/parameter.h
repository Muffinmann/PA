/*
 * parameter.h
 *
 *  Created on: Apr 5, 2016
 *      Author: iwtm84
 */

#ifndef CODE_VARIOUS_PARAMETER_H_
#define CODE_VARIOUS_PARAMETER_H_


#include <deal.II/base/parameter_handler.h>

void read_Parameterfile(std::string filename, ParameterHandler &parameter) {
	parameter.declare_entry("EModul", "136300", Patterns::Double(),
			"EModul in MPa");
	parameter.declare_entry("nu", "0.37", Patterns::Double(), "nu");
	parameter.declare_entry("mu", "127400", Patterns::Double(), "mu in MPa");
	parameter.declare_entry("h0", "1", Patterns::Double(),
			"Hardening in Slipsystem");
	parameter.declare_entry("pi0", "60", Patterns::Double(),
			"Initial Slipresistance");
	parameter.declare_entry("piS", "120", Patterns::Double(),
			"Saturated Slipresistance");
	parameter.declare_entry("le", "0", Patterns::Double(),
			"energetic length for gradientplasticity");
	parameter.declare_entry("c1", "120", Patterns::Double(),
			"constant for skrewdislocations");
	parameter.declare_entry("c2", "120", Patterns::Double(),
			"constant for edgedislocations");
	parameter.declare_entry("hardeningtype", "1", Patterns::Double(),
			"Hardeningtype");
	parameter.declare_entry("eta", "1000", Patterns::Double(),
			"Eta of Nortonslaw");
	parameter.declare_entry("gamma0dot", "0.0001", Patterns::Double(),
			"Reference sliprate");
	parameter.declare_entry("p", "10", Patterns::Double(),
			"Exponent of Nortonslaw");
	parameter.declare_entry("Nh","1",Patterns::Double(),
			"Exponent of gammaAlpha in Bardella");
	parameter.declare_entry("regularisation", "tanh_lin",
			Patterns::Selection("tanh_lin|powerlaw_sech|bardella_lin"), "Regularsation");//TODO
	parameter.declare_entry("alpha", "16.8e-6", Patterns::Double(),
			"Thermal expansion coefficient");

	parameter.declare_entry("GB_behaviour", "detailed",
			Patterns::Selection("detailed|microhard|microfree"),
			"Grain boundary behaviour");
	parameter.declare_entry("lambda_gb", "1", Patterns::Double(),
			"Coefficient for Grain Boundary Law");
	parameter.declare_entry("p_gb", "1", Patterns::Double(),
			"Rate sensitivity exponent for Grain Boundary Law");
	parameter.declare_entry("S_gb", "1", Patterns::Double(),
			"Hardening parameter for dissipative part of the Grain Boundary Law");

	parameter.declare_entry("n_slip", "12", Patterns::Selection("0|1|2|12"),
			"Number of Slipsystems");
	parameter.declare_entry("n_loadsteps", "20", Patterns::Integer(),
			"Number of Loadsteps");
	parameter.declare_entry("displacement", "0.01", Patterns::Double(),
			"Displacement in Tensildirection");
	parameter.declare_entry("strainrate", "0.02", Patterns::Double(),
			"Rate of norm of the strain");
	parameter.declare_entry("n_refinements", "2", Patterns::Integer(),
			"Number of Meshrefinements");
	parameter.declare_entry("isotrop", "false", Patterns::Bool(),
			"Isotrop or long grains");

	parameter.declare_entry("crosssnake", "false", Patterns::Bool(),
			"zigzag RVE or not");
	parameter.declare_entry("csAngle", "5", Patterns::Double(),
			"Tiltangle of grains in direction of meltpool movement");
	parameter.declare_entry("csHeight", "0.2", Patterns::Double(),
			"height of one layer");
	parameter.declare_entry("elementsCSLayer", "2", Patterns::Integer(),
			"number of elements per csLayer");

	parameter.declare_entry("n_crystals", "4", Patterns::Integer(),
			"Number of crystals per direction");
	parameter.declare_entry("loadtype", "periodic",
			Patterns::Selection(
					"periodic|periodic_uniaxial_tension|periodic_plane_stress|linear|linear_uniaxial_tension|linear_plane_stress|linear_CG|linear_CG_uniaxial_tension|linear_CG_plane_stress|testcase_andrew_microhard|testcase_andrew_microfree|testcase_grain_boundary|testcase_gottschalk"),
			"Loadtype (periodic|periodic_uniaxial_tension|periodic_plane_stress|linear|linear_uniaxial_tension|linear_plane_stress|linear_CG|linear_CG_uniaxial_tension|linear_CG_plane_stress|testcase_andrew_microhard|testcase_andrew_microfree|testcase_grain_boundary|testcase_gottschalk)");
	parameter.declare_entry("straindirection", "11",
			Patterns::Selection("11|22|33|12|13|23|detailed"), "Loaddirection");
	parameter.declare_entry("planeStressDir1", "11",
			Patterns::Selection("11|22|33|12|13|23"),
			"Direction one for biaxial loading");
	parameter.declare_entry("planeStressDir2", "22",
			Patterns::Selection("11|22|33|12|13|23"),
			"Direction two for biaxial loading");
	parameter.declare_entry("strainpath", "linear_displacement",
			Patterns::Selection("prescribed_from_file|linear_displacement"),
			"Type of prescribed strainpath");
	parameter.declare_entry("tiltangle", "15", Patterns::Double(),
			"Tilt against Buildingdirection");
	parameter.declare_entry("mesh", "simple_mesh", Patterns::Anything(),
			"simple_mesh|out7|out7periodic|out12|out12periodic|out38|out38periodic");
	parameter.declare_entry("orientations", "orientations", Patterns::Selection("orientations_B4|orientations"),
			"Filenames: orientations_B4|orientations");
	parameter.declare_entry("t_start", "1000", Patterns::Double(),
			"Temperature start");
	parameter.declare_entry("t_end", "23", Patterns::Double(),
			"Temperature end");
	parameter.declare_entry("el_type_u", "Q2", Patterns::Selection("Q1|Q2|Q3"),
			"Elementtype (Q1|Q2|Q3)");
	parameter.declare_entry("el_type_gamma", "Q1",
			Patterns::Selection("Q1|Q2|Q3"), "Elementtype (Q1|Q2|Q3)");
	parameter.declare_entry("number_gp", "3", Patterns::Integer(),
			"Number of GP per direction");

	parameter.declare_entry("localCP", "false", Patterns::Bool(),
			"local on nonlocal crystal plasticity");
	parameter.declare_entry("newtonCorrection", "naive", Patterns::Anything(),
			"none|naive|improved");
	parameter.declare_entry("maxCGIterationsRel", "1", Patterns::Double(),
			"Max Iterations for CG Solver is ndof*maxCGIterationsRel");
	parameter.declare_entry("maxGlobalIterations", "10", Patterns::Integer(),
			"Maximum global Newtoniterations");
	parameter.declare_entry("tolCGsolver", "1e-10", Patterns::Double(),
			"Tolerance for CG Solver");
	parameter.declare_entry("globalTolNewton", "1e-8", Patterns::Double(),
			"Terminationcriterion for global NewtonRaphson");

	parameter.declare_entry("regularisationF", "1", Patterns::Double(),
			"RegularisationFaktor for tanh");
	parameter.declare_entry("regularIncreaseF", "1.1", Patterns::Double(),
			"Faktor to increase RegularisationFaktor for tanh");
	parameter.declare_entry("NewtonTolIncReF", "1e-3", Patterns::Double(),
			"Tolerance to increase RegularisationFaktor for tanh");
	parameter.declare_entry("targetConvergence", "2", Patterns::Double(),
			"Target convergence for last newtonstep in each reg increase loop, 2 for quadratic");

//	parameter.declare_entry("DampingFactorNewton", "1", Patterns::Double(),
//			"Damping factor for NewtonRaphson");
//	parameter.declare_entry("Relaxation_Type", "static", Patterns::Anything(),
//			"static|dynamic");
//	parameter.declare_entry("lb_DR", "0.9", Patterns::Double(),
//			"Lower Bound for Dynamic Relaxation");
//	parameter.declare_entry("ub_DR", "1.0", Patterns::Double(),
//			"Upper Bound for Dynamic Relaxation");
	parameter.declare_entry("MP_order", "-1", Patterns::Integer(),
			"Order of mechanical prediction, -1 for guess=0, 0 for guess = oldNewtonUpdate, ...");
	parameter.declare_entry("mode", "solve",
			Patterns::Selection(
					"solve|calculateTangent|calculateYieldSurface"),
			"Calculationmode");
	parameter.declare_entry("pertubation", "1e-6", Patterns::Double(),
			"Pertubation for the calculation of macroTangent");
	parameter.declare_entry("nCoresAssembly", "8", Patterns::Integer(),
			"Number of cores for assembly with Workstream");
	parameter.declare_entry("writeOutVectors", "false", Patterns::Bool(),
			"write solution vectors for analysis in file");

	parameter.declare_entry("macroTangent", "", Patterns::List(Patterns::Double(), 0, 36), "Macroscopic Tangent rowwise");
	parameter.declare_entry("nPointsYS", "1", Patterns::Integer(),
			"number of points on yield surface including the uniaxial directions");
	parameter.declare_entry("startPointYS", "0", Patterns::Integer(),
			"First point to calculate on Yieldsurface");
	parameter.declare_entry("endPointYS", "999", Patterns::Integer(),
			"Last point to calculate on Yieldsurface");

	parameter.declare_entry("postprocessing", "full,", Patterns::Anything(),
			"1:full, 2:EModul, 3:stress_strain_curve");
	parameter.declare_entry("randomly_oriented", "false", Patterns::Bool(),
			"Orientation of Crystalls with prescribed list or randomly");
	parameter.declare_entry("print_crystal_data", "false", Patterns::Bool(),
			"Print interaction moduli and slipsystem informations");

	parameter.declare_entry("euler1", "0", Patterns::Double(), "Eulerangle 1");
	parameter.declare_entry("euler2", "0", Patterns::Double(), "Eulerangle 2");
	parameter.declare_entry("euler3", "0", Patterns::Double(), "Eulerangle 3");
	parameter.declare_entry("offset_domain", "0", Patterns::Integer(),
			"offset for angles");

	parameter.parse_input(filename);

	//-------------------------------------------------------------------
	// set degree of Ansatzfunktions depending on Elementtype
	string feq_u = "";
	if (parameter.get("el_type_u").compare("Q1") == 0) {
		feq_u = "1";
	} else if (parameter.get("el_type_u").compare("Q2") == 0) {
		feq_u = "2";
	} else if (parameter.get("el_type_u").compare("Q3") == 0) {
		feq_u = "3";
	}
	string feq_gamma = "";
	if (parameter.get("el_type_gamma").compare("Q1") == 0) {
		feq_gamma = "1";
	} else if (parameter.get("el_type_gamma").compare("Q2") == 0) {
		feq_gamma = "2";
	} else if (parameter.get("el_type_gamma").compare("Q3") == 0) {
		feq_gamma = "3";
	}
	parameter.declare_entry("feq_u", feq_u, Patterns::Double(),
			"Polynomial Degree of Displacement Ansatzfunction");
	parameter.declare_entry("feq_gamma", feq_gamma, Patterns::Double(),
			"Polynomial Degree of Slip Ansatzfunction");
	if (parameter.get_double("feq_u") - parameter.get_double("feq_gamma")
			!= 1) {
		deallog
				<< "Parameterhandler missmatch of degree of Ansatzfunktions u and gamma"
				<< endl;
	}
	//-------------------------------------------------------------------

	//-------------------------------------------------------------------
	// check loaddirection
	// has to be 12 for testcase_andrew and 11|22|33 for uniaxial tension
	if (parameter.get("loadtype").compare("testcase_andrew_microhard") == 0
			|| parameter.get("loadtype").compare("testcase_andrew_microfree")
					== 0) {
		if (parameter.get("straindirection").compare("12") != 0) {
			parameter.set("straindirection", "12");
			deallog
					<< "Parameterhandler changed straindirection to 12 for testcase_andrew"
					<< endl;
		}

		if (parameter.get("mesh").compare("simple_mesh") != 0) {
			parameter.set("mesh", "simple_mesh");
			deallog
					<< "Parameterhandler changed mesh to simple_mesh for testcase_andrew"
					<< endl;
		}
		if (parameter.get_double("n_crystals") == 1) { // One crystal in 45 degree
			parameter.set("euler1", "0");
			parameter.set("euler2", "0");
			parameter.set("euler3", "0");
		} else {
			parameter.set("offset_domain", "4000");
		}
	}

	// Testcase Grain boundary singe slip
	if (parameter.get("loadtype").compare("testcase_grain_boundary") == 0) {
		if (parameter.get("straindirection").compare("12") != 0) {
			parameter.set("straindirection", "12");
			deallog
					<< "Parameterhandler changed straindirection to 12 for testcase_grain_boundary"
					<< endl;
		}

		if (parameter.get("mesh").compare("simple_mesh") != 0) {
			parameter.set("mesh", "simple_mesh");
			deallog
					<< "Parameterhandler changed mesh to simple_mesh for testcase_grain_boundary"
					<< endl;
		}
		if (parameter.get_double("n_crystals") == 1) { // One crystal in 45 degree
			parameter.set("euler1", "0");
			parameter.set("euler2", "0");
			parameter.set("euler3", "0");
		} else {
			parameter.set("offset_domain", "4000");
		}
	}

	if (parameter.get("loadtype").compare("periodic_uniaxial_tension") == 0) {
		if (parameter.get("straindirection").compare("11") != 0
				&& parameter.get("straindirection").compare("22") != 0
				&& parameter.get("straindirection").compare("33") != 0) {
			deallog
					<< "Parameterhandler straindirection not 11, 22 or 33 for periodic uniaxial tension, think if this is what you want."
					<< endl;
		}
	}

	if (!parameter.get_bool("crosssnake")) {
		parameter.set("csAngle", "0");
	}else{
//		if (parameter.get_bool("randomly_oriented")){
//			deallog
//					<< "Crosssnake and random orientation not implemented yet."
//					<< endl;
//			throw std::bad_exception();
//		}
	}

	if (parameter.get("mode") == "calculateYieldSurface") {
		if (parameter.get_integer("nPointsYS") % 2 == 1) {
			deallog
					<< "Yieldsurface: (nPointsYS - 1) not completely divisible by 2"
					<< endl;
			throw std::bad_exception();
		}
	}

}

#endif /* CODE_VARIOUS_PARAMETER_H_ */
