// $Id: FESystemStaticMemberDefinitions.C,v 1.25.4.9 2008-08-21 00:37:07 manav Exp $


#include "Loads/load.h"
#include "Loads/LoadCase.h"
#include "Loads/LoadSet.h"

#include "Discipline/AnalysisDisciplineBase.h"
#include "Discipline/StructuralAnalysis.h"
#include "Discipline/ThermalAnalysis.h"
#include "Discipline/PistonTheory.h"
#include "Discipline/ThermalDisciplineInfo.h"
#include "Discipline/StructuralDisciplineInfo.h"
#include "Discipline/PistonTheoryInfo.h"
#include "Radiation/RadiationCavityAnalysis.h"


#include "AnalysisDriver/AnalysisDriver.h"
#include "AnalysisDriver/LinearAnalysisDriver.h"
#include "AnalysisDriver/NonLinearAnalysisDriver.h"
#include "AnalysisDriver/EigenProblemAnalysisDriver.h"
#include "AnalysisDriver/TransientAnalysisDriver.h"
#include "AnalysisDriver/RogerApproximationAeroelasticity.h"

#include "Solutions/SolutionBase.h"
#include "Solutions/ThermalSolution.h"
#include "Solutions/TransientThermalSolution.h"
#include "Solutions/TransientStructuralSolution.h"
#include "Solutions/LinearStressSolution.h"
#include "Solutions/StructuralVibrationEigenSolution.h"
#include "Solutions/LinearizedBucklingEigenSolution.h"
#include "Solutions/AeroelasticitySolution.h"


#include "DesignData/DesignParameter.h"
#include "DesignData/ShapeParameter.h"
#include "DesignData/PropertyParameter.h"

#include "FESystem/FESystemElem.h"

#include "Interpolation/FEInterpolationElem.h"
#include "Interpolation/InterpolationBar2.h"
#include "Interpolation/InterpolationQuad4.h"
#include "Interpolation/InterpolationTri3.h"
#include "Interpolation/InterpolationHex8.h"

#include "Numerics/FunctionBase.h"
#include "Numerics/MultilinearFunction.h"

#include "OutputProcessor/OutputProcessor.h"
#include "OutputProcessor/GMVOutputProcessor.h"
#include "OutputProcessor/GmshOutputProcessor.h"
#include "OutputProcessor/TecPlotOutputProcessor.h"

#include "Properties/ElemDataCard.h"
#include "Properties/MaterialPropertyNameEnums.h"
#include "Properties/PropertyCardParameter.h"
#include "Properties/Isotropic1DElemDataCard.h"
#include "Properties/Isotropic2DElemDataCard.h"
#include "Properties/Isotropic3DElemDataCard.h"
#include "Properties/IsotropicMaterialPropertyCard.h"
#include "Properties/ForcedConvection1DElemDataCard.h"


#include "Solvers/FESystemSolverBase.h"
#include "Solvers/SolverInfo.h"
#include "Solvers/LinearSolverInfo.h"
#include "Solvers/LinearSolver.h"
#include "Solvers/PetscLinearSolver.h"
#include "Solvers/NonlinearSolver.h"
#include "Solvers/NonlinearSolverInfo.h"
#include "Solvers/NewtonNonlinearSolver.h"
#include "Solvers/PetscNonlinearSolver.h"
#include "Solvers/EigenSolver.h"
#include "Solvers/SlepcEigenSolver.h"
#include "Solvers/ArpackEigenSolver.h"
#include "Solvers/EigenSolverInfo.h"
#include "Solvers/TransientSolver.h"
#include "Solvers/TransientSolverInfo.h"
#include "Solvers/NewmarkTransientSolver.h"
#include "Solvers/NewmarkTransientSolverInfo.h"
#include "Solvers/EulerTransientSolver.h"
#include "Solvers/EulerTransientSolverInfo.h"

#include "StructuralElems/structural_elem.h"
#include "StructuralElems/plate_MITC4.h"
#include "StructuralElems/membrane_tri3.h"
#include "StructuralElems/membrane_quad4.h"
#include "StructuralElems/bar.h"
#include "StructuralElems/beam2.h"
#include "StructuralElems/brick_hex8.h"
#include "StructuralElems/linear_spring.h"
#include "StructuralElems/PlateDKT.h"
#include "StructuralElems/Tri3VonKarman.h"


#include "ThermalElems/thermal_elem.h"
#include "ThermalElems/conduction_1d.h"
#include "ThermalElems/conduction_hex8.h"
#include "ThermalElems/conduction_prism6.h"
#include "ThermalElems/conduction_quad4.h"
#include "ThermalElems/conduction_tri3.h"
#include "ThermalElems/ForcedConvection1D.h"

#include "FESystem/FESystemElemTypeEnumHandler.h"

#include "FESystem/FESystemController.h"



// ********************* Packages being used  ******************************
NameEnumerationHandler FESystem::PackageNameEnum::enum_handler;

FESystem::PackageNameEnum::EnumRegistration<FESystem::PETSC_PACKAGE>
FESystem::PETSC_PACKAGE::registration; 


// ********************* Loads  ******************************
NameEnumerationHandler LoadSetKindEnum::enum_handler;
NameEnumerationHandler LoadCaseKindEnum::enum_handler;
NameEnumerationHandler LoadClassEnum::enum_handler;
NameEnumerationHandler LoadNameEnum::enum_handler;


LoadNameEnum::EnumRegistration<NODAL_FORCE> NODAL_FORCE::registration;
LoadNameEnum::EnumRegistration<NODAL_TEMPERATURE> NODAL_TEMPERATURE::registration;
LoadNameEnum::EnumRegistration<SURFACE_PRESSURE> SURFACE_PRESSURE::registration;
LoadNameEnum::EnumRegistration<DISPLACEMENT_BOUNDARY_CONDITION> 
DISPLACEMENT_BOUNDARY_CONDITION::registration;
LoadNameEnum::EnumRegistration<SELF_WEIGHT> SELF_WEIGHT::registration;


LoadNameEnum::EnumRegistration<TEMPERATURE_BOUNDARY_CONDITION> 
TEMPERATURE_BOUNDARY_CONDITION::registration;
LoadNameEnum::EnumRegistration<NODAL_HEAT_LOAD> NODAL_HEAT_LOAD::registration;
LoadNameEnum::EnumRegistration<VOLUME_HEAT_LOAD> VOLUME_HEAT_LOAD::registration;
LoadNameEnum::EnumRegistration<SURFACE_HEAT_LOAD> SURFACE_HEAT_LOAD::registration;
LoadNameEnum::EnumRegistration<SURFACE_RADIATION_HEAT_LOAD> SURFACE_RADIATION_HEAT_LOAD::registration;
LoadNameEnum::EnumRegistration<SURFACE_CONVECTION_HEAT_LOAD> SURFACE_CONVECTION_HEAT_LOAD::registration;

LoadNameEnum::EnumRegistration<PISTON_THEORY_SURFACE> PISTON_THEORY_SURFACE::registration;


LoadCaseKindEnum::EnumRegistration<STATIC_LOAD_CASE> STATIC_LOAD_CASE::registration;
LoadCaseKindEnum::EnumRegistration<TRANSIENT_LOAD_CASE> TRANSIENT_LOAD_CASE::registration;

LoadSetKindEnum::EnumRegistration<VOLUME_LOAD_SET> VOLUME_LOAD_SET::registration;
LoadSetKindEnum::EnumRegistration<SURFACE_LOAD_SET> SURFACE_LOAD_SET::registration;
LoadSetKindEnum::EnumRegistration<NODAL_LOAD_SET> NODAL_LOAD_SET::registration;
LoadSetKindEnum::EnumRegistration<BOUNDARY_CONDITION_LOAD_SET> BOUNDARY_CONDITION_LOAD_SET::registration;


LoadSetKindEnum::EnumRegistration<AEROELASTIC_LOAD_SET> AEROELASTIC_LOAD_SET::registration;


LoadClassEnum::EnumRegistration<NODAL_POINT_LOAD> NODAL_POINT_LOAD::registration;
LoadClassEnum::EnumRegistration<DIRICHLET_BOUNDARY_CONDITION> DIRICHLET_BOUNDARY_CONDITION::registration;
LoadClassEnum::EnumRegistration<SCALAR_VOLUME_LOAD> SCALAR_VOLUME_LOAD::registration;
LoadClassEnum::EnumRegistration<VECTOR_VOLUME_LOAD> VECTOR_VOLUME_LOAD::registration;
LoadClassEnum::EnumRegistration<SURFACE_RADIATION_LOAD> SURFACE_RADIATION_LOAD::registration;
LoadClassEnum::EnumRegistration<SCALAR_SURFACE_LOAD> SCALAR_SURFACE_LOAD::registration;
LoadClassEnum::EnumRegistration<SURFACE_CONVECTION_LOAD> SURFACE_CONVECTION_LOAD::registration;



// **********************  Discipline ***************************
NameEnumerationHandler Discipline::AnalysisDisciplineEnum::enum_handler;
NameEnumerationHandler Discipline::DisciplineAnalysisTypeEnum::enum_handler;
NameEnumerationHandler Discipline::DisciplineInfoTypeEnum::enum_handler;


Discipline::DisciplineAnalysisTypeEnum::EnumRegistration<Discipline::LINEAR_ANALYSIS> 
Discipline::LINEAR_ANALYSIS::registration;
Discipline::DisciplineAnalysisTypeEnum::EnumRegistration<Discipline::NONLINEAR_ANALYSIS> 
Discipline::NONLINEAR_ANALYSIS::registration;

Discipline::AnalysisDisciplineEnum::EnumRegistration<Discipline::STRUCTURAL_DISCIPLINE> 
Discipline::STRUCTURAL_DISCIPLINE::registration;

Discipline::AnalysisDisciplineEnum::EnumRegistration<Discipline::RADIATION_DISCIPLINE> 
Discipline::RADIATION_DISCIPLINE::registration;



Discipline::AnalysisDisciplineEnum::EnumRegistration<Discipline::THERMAL_DISCIPLINE> 
Discipline::THERMAL_DISCIPLINE::registration;

Discipline::AnalysisDisciplineEnum::EnumRegistration<Discipline::PISTON_THEORY> 
Discipline::PISTON_THEORY::registration;


Discipline::DisciplineInfoTypeEnum::EnumRegistration<Discipline::THERMAL_DISCIPLINE_INFO> 
Discipline::THERMAL_DISCIPLINE_INFO::registration;
Discipline::DisciplineInfoTypeEnum::EnumRegistration<Discipline::STRUCTURAL_DISCIPLINE_INFO> 
Discipline::STRUCTURAL_DISCIPLINE_INFO::registration;
Discipline::DisciplineInfoTypeEnum::EnumRegistration<Discipline::PISTON_THEORY_INFO> 
Discipline::PISTON_THEORY_INFO::registration;

// ********************* Driver ***********************************
NameEnumerationHandler Driver::AnalysisDriverTypeEnum::enum_handler;
NameEnumerationHandler Driver::AnalysisTypeEnum::enum_handler;
NameEnumerationHandler Driver::AnalysisDriverQtyEnum::enum_handler;

Driver::AnalysisDriverTypeEnum::EnumRegistration<Driver::LINEAR_ANALYSIS_DRIVER> 
Driver::LINEAR_ANALYSIS_DRIVER::registration;
Driver::AnalysisDriverTypeEnum::EnumRegistration<Driver::NONLINEAR_ANALYSIS_DRIVER> 
Driver::NONLINEAR_ANALYSIS_DRIVER::registration;
Driver::AnalysisDriverTypeEnum::EnumRegistration<Driver::EIGENPROBLEM_ANALYSIS_DRIVER> 
Driver::EIGENPROBLEM_ANALYSIS_DRIVER::registration;
Driver::AnalysisDriverTypeEnum::EnumRegistration<Driver::LINEAR_TRANSIENT_ANALYSIS_DRIVER> 
Driver::LINEAR_TRANSIENT_ANALYSIS_DRIVER::registration;
Driver::AnalysisDriverTypeEnum::EnumRegistration<Driver::NONLINEAR_TRANSIENT_ANALYSIS_DRIVER> 
Driver::NONLINEAR_TRANSIENT_ANALYSIS_DRIVER::registration;
Driver::AnalysisDriverTypeEnum::EnumRegistration<Driver::ROGER_APPROXIMATION_AEROELASTICITY_DRIVER> 
Driver::ROGER_APPROXIMATION_AEROELASTICITY_DRIVER::registration;

Driver::AnalysisTypeEnum::EnumRegistration<Driver::ANALYSIS> Driver::ANALYSIS::registration;
Driver::AnalysisTypeEnum::EnumRegistration<Driver::SENSITIVITY_ANALYSIS> 
Driver::SENSITIVITY_ANALYSIS::registration;

Driver::AnalysisDriverQtyEnum::EnumRegistration<Driver::SYSTEM_MATRIX> 
Driver::SYSTEM_MATRIX::registration;
Driver::AnalysisDriverQtyEnum::EnumRegistration<Driver::FORCE_VECTOR> 
Driver::FORCE_VECTOR::registration;
//Driver::AnalysisDriverQtyEnum::EnumRegistration<Driver::RESIDUAL_VECTOR> 
//Driver::RESIDUAL_VECTOR::registration;
Driver::AnalysisDriverQtyEnum::EnumRegistration<Driver::JACOBIAN_MATRIX> 
Driver::JACOBIAN_MATRIX::registration;
Driver::AnalysisDriverQtyEnum::EnumRegistration<Driver::SYSTEM_MATRIX_SENSITIVITY> 
Driver::SYSTEM_MATRIX_SENSITIVITY::registration;
Driver::AnalysisDriverQtyEnum::EnumRegistration<Driver::FORCE_VECTOR_SENSITIVITY> 
Driver::FORCE_VECTOR_SENSITIVITY::registration;
Driver::AnalysisDriverQtyEnum::EnumRegistration<Driver::EIGENPROBLEM_A_MATRIX> 
Driver::EIGENPROBLEM_A_MATRIX::registration;
Driver::AnalysisDriverQtyEnum::EnumRegistration<Driver::EIGENPROBLEM_A_MATRIX_SENSITIVITY> 
Driver::EIGENPROBLEM_A_MATRIX_SENSITIVITY::registration;
Driver::AnalysisDriverQtyEnum::EnumRegistration<Driver::EIGENPROBLEM_B_MATRIX> 
Driver::EIGENPROBLEM_B_MATRIX::registration;
Driver::AnalysisDriverQtyEnum::EnumRegistration<Driver::EIGENPROBLEM_B_MATRIX_SENSITIVITY> 
Driver::EIGENPROBLEM_B_MATRIX_SENSITIVITY::registration;
Driver::AnalysisDriverQtyEnum::EnumRegistration<Driver::MODEL_MASS> 
Driver::MODEL_MASS::registration;
Driver::AnalysisDriverQtyEnum::EnumRegistration<Driver::MODEL_MASS_SENSITIVITY> 
Driver::MODEL_MASS_SENSITIVITY::registration;



// ********************* Design Data *****************************
NameEnumerationHandler DesignData::DesignParameterTypeEnum::enum_handler;
NameEnumerationHandler DesignData::SensitivityMethodEnum::enum_handler;

DesignData::SensitivityMethodEnum::EnumRegistration<DesignData::EULER_FD_SENSITIVITY> 
DesignData::EULER_FD_SENSITIVITY::registration;
DesignData::SensitivityMethodEnum::EnumRegistration<DesignData::ANALYTIC_SENSITIVITY> 
DesignData::ANALYTIC_SENSITIVITY::registration;

DesignData::DesignParameterTypeEnum::EnumRegistration<DesignData::SHAPE_PARAMETER> 
DesignData::SHAPE_PARAMETER::registration;
DesignData::DesignParameterTypeEnum::EnumRegistration<DesignData::PROPERTY_PARAMETER> 
DesignData::PROPERTY_PARAMETER::registration;



// ********************* Functions *****************************
NameEnumerationHandler FunctionTypeEnum::enum_handler;

FunctionTypeEnum::EnumRegistration<MULTILINEAR_FUNCTION> MULTILINEAR_FUNCTION::registration;

// ********************** Output *****************************
NameEnumerationHandler OutputFormatEnum::enum_handler;

OutputFormatEnum::EnumRegistration<GMV_OUTPUT_PROCESSOR> GMV_OUTPUT_PROCESSOR::registration;
OutputFormatEnum::EnumRegistration<TECPLOT_OUTPUT_PROCESSOR> 
TECPLOT_OUTPUT_PROCESSOR::registration;
OutputFormatEnum::EnumRegistration<GMSH_OUTPUT_PROCESSOR> GMSH_OUTPUT_PROCESSOR::registration;



// ********************** Properties *************************
NameEnumerationHandler ElemDataCardEnum::enum_handler;
NameEnumerationHandler Isotropic_1D_ElemDataEnum::enum_handler;
NameEnumerationHandler Isotropic_2D_ElemDataEnum::enum_handler;
NameEnumerationHandler Isotropic_3D_ElemDataEnum::enum_handler;

NameEnumerationHandler ElemDataCardFactorsEnum::enum_handler;

MaterialPropertyNameEnumerationHandler::PropertyKindMapType 
MaterialPropertyNameEnumerationHandler::property_kind_map;

NameEnumerationHandler MaterialPropertyCardEnum::enum_handler;
NameEnumerationHandler Property::PropertyCardParameterType::enum_handler;

NameEnumerationHandler Property::LocalParameterType::enum_handler;



// ******************  Solutions ************************* 
NameEnumerationHandler Solution::SolutionEnum::enum_handler;

Solution::SolutionEnum::EnumRegistration<Solution::THERMAL_SOLUTION> 
Solution::THERMAL_SOLUTION::registration;
Solution::SolutionEnum::EnumRegistration<Solution::TRANSIENT_THERMAL_SOLUTION> 
Solution::TRANSIENT_THERMAL_SOLUTION::registration;
Solution::SolutionEnum::EnumRegistration<Solution::TRANSIENT_STRUCTURAL_SOLUTION> 
Solution::TRANSIENT_STRUCTURAL_SOLUTION::registration;
Solution::SolutionEnum::EnumRegistration<Solution::LINEAR_STRESS_SOLUTION> 
Solution::LINEAR_STRESS_SOLUTION::registration;
Solution::SolutionEnum::EnumRegistration<Solution::STRUCTURAL_VIBRATION_EIGEN_SOLUTION> 
Solution::STRUCTURAL_VIBRATION_EIGEN_SOLUTION::registration;
Solution::SolutionEnum::EnumRegistration<Solution::LINEARIZED_BUCKLING_EIGEN_SOLUTION> 
Solution::LINEARIZED_BUCKLING_EIGEN_SOLUTION::registration;
Solution::SolutionEnum::EnumRegistration<Solution::AEROELASTICITY_SOLUTION> 
Solution::AEROELASTICITY_SOLUTION::registration;


// ******************  SOLVER ************************* 
NameEnumerationHandler Solver::SolverClassEnum::enum_handler;

NameEnumerationHandler Solver::SolverInfoEnum::enum_handler;


// ****************** Linear Solver ************************
NameEnumerationHandler Solver::LinearSolverKindEnum::enum_handler;
NameEnumerationHandler Solver::PCTypeEnum::enum_handler;
NameEnumerationHandler Solver::KSPTypeEnum::enum_handler;

Solver::SolverInfoEnum::EnumRegistration<Solver::LINEAR_SOLVER_INFO> 
Solver::LINEAR_SOLVER_INFO::registration;

Solver::SolverClassEnum::EnumRegistration<Solver::LINEAR_SOLVER> 
Solver::LINEAR_SOLVER::registration;

Solver::LinearSolverKindEnum::EnumRegistration<Solver::PETSC_LINEAR_SOLVER> 
Solver::PETSC_LINEAR_SOLVER::registration;


Solver::PCTypeEnum::EnumRegistration<Solver::PC_LU> 
Solver::PC_LU::registration;
Solver::PCTypeEnum::EnumRegistration<Solver::PC_JACOBI> 
Solver::PC_JACOBI::registration;
Solver::PCTypeEnum::EnumRegistration<Solver::PC_BJACOBI> 
Solver::PC_BJACOBI::registration;
Solver::PCTypeEnum::EnumRegistration<Solver::PC_SOR> 
Solver::PC_SOR::registration;
Solver::PCTypeEnum::EnumRegistration<Solver::PC_EISENSTAT> 
Solver::PC_EISENSTAT::registration;
Solver::PCTypeEnum::EnumRegistration<Solver::PC_ICC> 
Solver::PC_ICC::registration;
Solver::PCTypeEnum::EnumRegistration<Solver::PC_ILU> 
Solver::PC_ILU::registration;
Solver::PCTypeEnum::EnumRegistration<Solver::PC_ASM> 
Solver::PC_ASM::registration;
Solver::PCTypeEnum::EnumRegistration<Solver::PC_KSP> 
Solver::PC_KSP::registration;
Solver::PCTypeEnum::EnumRegistration<Solver::PC_CHOLESKY> 
Solver::PC_CHOLESKY::registration;
Solver::PCTypeEnum::EnumRegistration<Solver::PC_NONE> 
Solver::PC_NONE::registration;

Solver::KSPTypeEnum::EnumRegistration<Solver::KSP_RICHARDSON> 
Solver::KSP_RICHARDSON::registration;
Solver::KSPTypeEnum::EnumRegistration<Solver::KSP_CHEBYCHEV> 
Solver::KSP_CHEBYCHEV::registration;
Solver::KSPTypeEnum::EnumRegistration<Solver::KSP_CG> 
Solver::KSP_CG::registration;
Solver::KSPTypeEnum::EnumRegistration<Solver::KSP_GMRES> 
Solver::KSP_GMRES::registration;
Solver::KSPTypeEnum::EnumRegistration<Solver::KSP_TCQMR> 
Solver::KSP_TCQMR::registration;
Solver::KSPTypeEnum::EnumRegistration<Solver::KSP_BCG> 
Solver::KSP_BCG::registration;
Solver::KSPTypeEnum::EnumRegistration<Solver::KSP_CGS> 
Solver::KSP_CGS::registration;
Solver::KSPTypeEnum::EnumRegistration<Solver::KSP_TFQMR> 
Solver::KSP_TFQMR::registration;
Solver::KSPTypeEnum::EnumRegistration<Solver::KSP_CR> 
Solver::KSP_CR::registration;
Solver::KSPTypeEnum::EnumRegistration<Solver::KSP_LSQR> 
Solver::KSP_LSQR::registration;
Solver::KSPTypeEnum::EnumRegistration<Solver::KSP_BICG> 
Solver::KSP_BICG::registration;
Solver::KSPTypeEnum::EnumRegistration<Solver::KSP_PRE_ONLY> 
Solver::KSP_PRE_ONLY::registration;



// ****************** Nonlinear Solver ************************
NameEnumerationHandler Solver::NonlinearSolverKindEnum::enum_handler;

Solver::SolverInfoEnum::EnumRegistration<Solver::NONLINEAR_SOLVER_INFO> 
Solver::NONLINEAR_SOLVER_INFO::registration;

Solver::SolverClassEnum::EnumRegistration<Solver::NONLINEAR_SOLVER> 
Solver::NONLINEAR_SOLVER::registration;
Solver::NonlinearSolverKindEnum::EnumRegistration<Solver::FESYSTEM_NEWTON_NONLINEAR_SOLVER> 
Solver::FESYSTEM_NEWTON_NONLINEAR_SOLVER::registration;
Solver::NonlinearSolverKindEnum::EnumRegistration<Solver::NEWTON_NONLINEAR_SOLVER_CUBIC_LINE_SEARCH> 
Solver::NEWTON_NONLINEAR_SOLVER_CUBIC_LINE_SEARCH::registration;
Solver::NonlinearSolverKindEnum::EnumRegistration<Solver::NEWTON_NONLINEAR_SOLVER_QUADRATIC_LINE_SEARCH> 
Solver::NEWTON_NONLINEAR_SOLVER_QUADRATIC_LINE_SEARCH::registration;
Solver::NonlinearSolverKindEnum::EnumRegistration<Solver::NEWTON_NONLINEAR_SOLVER_NO_LINE_SEARCH> 
Solver::NEWTON_NONLINEAR_SOLVER_NO_LINE_SEARCH::registration;
Solver::NonlinearSolverKindEnum::EnumRegistration<Solver::NEWTON_NONLINEAR_SOLVER_NO_NORMS_LINE_SEARCH> 
Solver::NEWTON_NONLINEAR_SOLVER_NO_NORMS_LINE_SEARCH::registration;


// ****************** Eigen Solver ************************
NameEnumerationHandler Solver::EigenProblemKindEnum::enum_handler;
NameEnumerationHandler Solver::EigenSolutionSpectrumKindEnum::enum_handler;
NameEnumerationHandler Solver::EigenSolutionShiftKindEnum::enum_handler;



Solver::SolverClassEnum::EnumRegistration<Solver::EIGEN_SOLVER> 
Solver::EIGEN_SOLVER::registration;
Solver::EigenProblemKindEnum::EnumRegistration<Solver::HERMITIAN_EIGENPROBLEM> 
Solver::HERMITIAN_EIGENPROBLEM::registration;
Solver::EigenProblemKindEnum::EnumRegistration<Solver::GENERALIZED_HERMITIAN_EIGENPROBLEM> 
Solver::GENERALIZED_HERMITIAN_EIGENPROBLEM::registration;
Solver::EigenProblemKindEnum::EnumRegistration<Solver::NON_HERMITIAN_EIGENPROBLEM> 
Solver::NON_HERMITIAN_EIGENPROBLEM::registration;
Solver::EigenProblemKindEnum::EnumRegistration<Solver::GENERALIZED_NON_HERMITIAN_EIGENPROBLEM> 
Solver::GENERALIZED_NON_HERMITIAN_EIGENPROBLEM::registration;

Solver::EigenSolutionSpectrumKindEnum::EnumRegistration<Solver::LARGEST_MAGNITUDE> 
Solver::LARGEST_MAGNITUDE::registration;
Solver::EigenSolutionSpectrumKindEnum::EnumRegistration<Solver::SMALLEST_MAGNITUDE> 
Solver::SMALLEST_MAGNITUDE::registration;
Solver::EigenSolutionSpectrumKindEnum::EnumRegistration<Solver::LARGEST_REAL> 
Solver::LARGEST_REAL::registration;
Solver::EigenSolutionSpectrumKindEnum::EnumRegistration<Solver::SMALLEST_REAL> 
Solver::SMALLEST_REAL::registration;
Solver::EigenSolutionSpectrumKindEnum::EnumRegistration<Solver::LARGEST_IMAGINARY>
Solver::LARGEST_IMAGINARY::registration;
Solver::EigenSolutionSpectrumKindEnum::EnumRegistration<Solver::SMALLEST_IMAGINARY> 
Solver::SMALLEST_IMAGINARY::registration;


Solver::EigenSolutionShiftKindEnum::EnumRegistration<Solver::NO_SHIFT> 
Solver::NO_SHIFT::registration;
Solver::EigenSolutionShiftKindEnum::EnumRegistration<Solver::ORIGIN_SHIFT> 
Solver::ORIGIN_SHIFT::registration;
Solver::EigenSolutionShiftKindEnum::EnumRegistration<Solver::SPECTRUM_FOLD> 
Solver::SPECTRUM_FOLD::registration;
Solver::EigenSolutionShiftKindEnum::EnumRegistration<Solver::SHIFT_AND_INVERT> 
Solver::SHIFT_AND_INVERT::registration;
Solver::EigenSolutionShiftKindEnum::EnumRegistration<Solver::CAYLEY_SHIFT> 
Solver::CAYLEY_SHIFT::registration;


Solver::SolverInfoEnum::EnumRegistration<Solver::EIGEN_SOLVER_INFO> 
Solver::EIGEN_SOLVER_INFO::registration;


// ************************ Slepc Solver enum static declerations *************************

NameEnumerationHandler Solver::EigenSolverKindEnum::enum_handler;

Solver::EigenSolverKindEnum::EnumRegistration<Solver::EPS_LAPACK_EIGEN_SOLVER> 
Solver::EPS_LAPACK_EIGEN_SOLVER::registration;
Solver::EigenSolverKindEnum::EnumRegistration<Solver::EPS_POWER_EIGEN_SOLVER> 
Solver::EPS_POWER_EIGEN_SOLVER::registration;
Solver::EigenSolverKindEnum::EnumRegistration<Solver::EPS_SUBSPACE_EIGEN_SOLVER> 
Solver::EPS_SUBSPACE_EIGEN_SOLVER::registration;
Solver::EigenSolverKindEnum::EnumRegistration<Solver::EPS_ARNOLDI_EIGEN_SOLVER> 
Solver::EPS_ARNOLDI_EIGEN_SOLVER::registration;
Solver::EigenSolverKindEnum::EnumRegistration<Solver::EPS_LANCZOS_EIGEN_SOLVER> 
Solver::EPS_LANCZOS_EIGEN_SOLVER::registration;
Solver::EigenSolverKindEnum::EnumRegistration<Solver::EPS_KRYLOV_SCHUR_EIGEN_SOLVER> 
Solver::EPS_KRYLOV_SCHUR_EIGEN_SOLVER::registration;

Solver::EigenSolverKindEnum::EnumRegistration<Solver::ARPACK_EIGEN_SOLVER> 
Solver::ARPACK_EIGEN_SOLVER::registration;


// ****************** Transient Solver ************************
NameEnumerationHandler Solver::TransientSolverKindEnum::enum_handler;

Solver::SolverClassEnum::EnumRegistration<Solver::LINEAR_TRANSIENT_SOLVER> 
Solver::LINEAR_TRANSIENT_SOLVER::registration;
Solver::SolverClassEnum::EnumRegistration<Solver::NONLINEAR_TRANSIENT_SOLVER> 
Solver::NONLINEAR_TRANSIENT_SOLVER::registration;

Solver::SolverInfoEnum::EnumRegistration<Solver::LINEAR_TRANSIENT_SOLVER_INFO> 
Solver::LINEAR_TRANSIENT_SOLVER_INFO::registration;
Solver::SolverInfoEnum::EnumRegistration<Solver::NONLINEAR_TRANSIENT_SOLVER_INFO> 
Solver::NONLINEAR_TRANSIENT_SOLVER_INFO::registration;
Solver::SolverInfoEnum::EnumRegistration<Solver::NEWMARK_LINEAR_TRANSIENT_SOLVER_INFO> 
Solver::NEWMARK_LINEAR_TRANSIENT_SOLVER_INFO::registration;
Solver::SolverInfoEnum::EnumRegistration<Solver::NEWMARK_NONLINEAR_TRANSIENT_SOLVER_INFO> 
Solver::NEWMARK_NONLINEAR_TRANSIENT_SOLVER_INFO::registration;
Solver::SolverInfoEnum::EnumRegistration<Solver::EULER_LINEAR_TRANSIENT_SOLVER_INFO> 
Solver::EULER_LINEAR_TRANSIENT_SOLVER_INFO::registration;
Solver::SolverInfoEnum::EnumRegistration<Solver::EULER_NONLINEAR_TRANSIENT_SOLVER_INFO> 
Solver::EULER_NONLINEAR_TRANSIENT_SOLVER_INFO::registration;


Solver::TransientSolverKindEnum::EnumRegistration<Solver::NEWMARK_LINEAR_TRANSIENT_SOLVER> 
Solver::NEWMARK_LINEAR_TRANSIENT_SOLVER::registration;
Solver::TransientSolverKindEnum::EnumRegistration<Solver::NEWMARK_NONLINEAR_TRANSIENT_SOLVER> 
Solver::NEWMARK_NONLINEAR_TRANSIENT_SOLVER::registration;
Solver::TransientSolverKindEnum::EnumRegistration<Solver::EULER_LINEAR_TRANSIENT_SOLVER> 
Solver::EULER_LINEAR_TRANSIENT_SOLVER::registration;
Solver::TransientSolverKindEnum::EnumRegistration<Solver::EULER_NONLINEAR_TRANSIENT_SOLVER> 
Solver::EULER_NONLINEAR_TRANSIENT_SOLVER::registration;


// ************************ FESystem Element Base *************************
FESystemElem::FESystemElemNameEnumerationHandler FESystemElem::FESystemElemTypeEnum::enum_handler;
FESystemElem::FESystemElemNameEnumerationHandler::ElemKindMapType
FESystemElem::FESystemElemNameEnumerationHandler::elem_kind_map;
NameEnumerationHandler FESystemElem::FESystemElemQtyEnum::enum_handler;
NameEnumerationHandler FESystemElem::IntegrationDomainEnum::enum_handler;
NameEnumerationHandler FESystemElem::DesignPointElemEnum::enum_handler;

FESystemElem::FESystemElemQtyEnum::EnumRegistration<FESystemElem::TRANSFORM_MATRIX> 
FESystemElem::TRANSFORM_MATRIX::registration;

FESystemElem::FESystemElemQtyEnum::EnumRegistration<FESystemElem::N_N_FACTOR> 
FESystemElem::N_N_FACTOR::registration;

FESystemElem::FESystemElemQtyEnum::EnumRegistration<FESystemElem::N_X_N_X_FACTOR> 
FESystemElem::N_X_N_X_FACTOR::registration;

FESystemElem::FESystemElemQtyEnum::EnumRegistration<FESystemElem::N_Y_N_Y_FACTOR> 
FESystemElem::N_Y_N_Y_FACTOR::registration;

FESystemElem::FESystemElemQtyEnum::EnumRegistration<FESystemElem::N_Z_N_Z_FACTOR> 
FESystemElem::N_Z_N_Z_FACTOR::registration;

FESystemElem::FESystemElemQtyEnum::EnumRegistration<FESystemElem::N_X_N_Y_FACTOR> 
FESystemElem::N_X_N_Y_FACTOR::registration;

FESystemElem::FESystemElemQtyEnum::EnumRegistration<FESystemElem::N_Y_N_Z_FACTOR> 
FESystemElem::N_Y_N_Z_FACTOR::registration;

FESystemElem::FESystemElemQtyEnum::EnumRegistration<FESystemElem::N_Z_N_X_FACTOR> 
FESystemElem::N_Z_N_X_FACTOR::registration;

FESystemElem::FESystemElemQtyEnum::EnumRegistration<FESystemElem::N_X_N_FACTOR> 
FESystemElem::N_X_N_FACTOR::registration;

FESystemElem::FESystemElemQtyEnum::EnumRegistration<FESystemElem::N_Y_N_FACTOR> 
FESystemElem::N_Y_N_FACTOR::registration;

FESystemElem::FESystemElemQtyEnum::EnumRegistration<FESystemElem::N_Z_N_FACTOR> 
FESystemElem::N_Z_N_FACTOR::registration;

FESystemElem::FESystemElemQtyEnum::EnumRegistration<FESystemElem::N_FACTOR> 
FESystemElem::N_FACTOR::registration;

FESystemElem::FESystemElemQtyEnum::EnumRegistration<FESystemElem::SURFACE_NORMAL> 
FESystemElem::SURFACE_NORMAL::registration;


FESystemElem::IntegrationDomainEnum::EnumRegistration<FESystemElem::ELEM_VOLUME> 
FESystemElem::ELEM_VOLUME::registration;

FESystemElem::IntegrationDomainEnum::EnumRegistration<FESystemElem::SIDE_ZERO> 
FESystemElem::SIDE_ZERO::registration;

FESystemElem::IntegrationDomainEnum::EnumRegistration<FESystemElem::SIDE_ONE> 
FESystemElem::SIDE_ONE::registration;

FESystemElem::IntegrationDomainEnum::EnumRegistration<FESystemElem::SIDE_TWO> 
FESystemElem::SIDE_TWO::registration;

FESystemElem::IntegrationDomainEnum::EnumRegistration<FESystemElem::SIDE_THREE> 
FESystemElem::SIDE_THREE::registration;

FESystemElem::IntegrationDomainEnum::EnumRegistration<FESystemElem::SIDE_FOUR> 
FESystemElem::SIDE_FOUR::registration;

FESystemElem::IntegrationDomainEnum::EnumRegistration<FESystemElem::SIDE_FIVE> 
FESystemElem::SIDE_FIVE::registration;

FESystemElem::DesignPointElemEnum::EnumRegistration<FESystemElem::BASE_ELEM> 
FESystemElem::BASE_ELEM::registration;

FESystemElem::DesignPointElemEnum::EnumRegistration<FESystemElem::BASE_PLUS_DELTA_ELEM> 
FESystemElem::BASE_PLUS_DELTA_ELEM::registration;

// ******************  Structural Elements ************************* 
NameEnumerationHandler FESystemElem::PlateMITC4QtyEnum::enum_handler;
NameEnumerationHandler FESystemElem::StructuralElemQtyEnum::enum_handler;
NameEnumerationHandler FESystemElem::PlateDKBatozQtyEnum::enum_handler;
NameEnumerationHandler FESystemElem::PlateVonKarmanQtyEnum::enum_handler;

FESystemElem::StructuralElemQtyEnum::EnumRegistration<FESystemElem::STRUCTURAL_M_MATRIX> 
FESystemElem::STRUCTURAL_M_MATRIX::registration;
FESystemElem::StructuralElemQtyEnum::EnumRegistration<FESystemElem::STRUCTURAL_C_MATRIX> 
FESystemElem::STRUCTURAL_C_MATRIX::registration;
FESystemElem::StructuralElemQtyEnum::EnumRegistration<FESystemElem::STRUCTURAL_K_MATRIX> 
FESystemElem::STRUCTURAL_K_MATRIX::registration;
FESystemElem::StructuralElemQtyEnum::EnumRegistration<FESystemElem::STRUCTURAL_K_G_MATRIX> 
FESystemElem::STRUCTURAL_K_G_MATRIX::registration;
FESystemElem::StructuralElemQtyEnum::EnumRegistration<FESystemElem::STRUCTURAL_F_T_VECTOR> 
FESystemElem::STRUCTURAL_F_T_VECTOR::registration;
FESystemElem::StructuralElemQtyEnum::EnumRegistration<FESystemElem::STRUCTURAL_F_PRESSURE_VECTOR> 
FESystemElem::STRUCTURAL_F_PRESSURE_VECTOR::registration;
FESystemElem::StructuralElemQtyEnum::EnumRegistration<FESystemElem::STRUCTURAL_STRAIN_OPERATOR> 
FESystemElem::STRUCTURAL_STRAIN_OPERATOR::registration;


FESystemElem::FESystemElemTypeEnum::EnumRegistration<FESystemElem::STRUCTURAL_BAR_EDGE2> 
FESystemElem::STRUCTURAL_BAR_EDGE2::registration;
FESystemElem::FESystemElemTypeEnum::EnumRegistration<FESystemElem::STRUCTURAL_BEAM_EDGE2> 
FESystemElem::STRUCTURAL_BEAM_EDGE2::registration;
FESystemElem::FESystemElemTypeEnum::EnumRegistration<FESystemElem::STRUCTURAL_LINEAR_SPRING_EDGE2> 
FESystemElem::STRUCTURAL_LINEAR_SPRING_EDGE2::registration;
FESystemElem::FESystemElemTypeEnum::EnumRegistration<FESystemElem::STRUCTURAL_MEMBRANE_QUAD4> 
FESystemElem::STRUCTURAL_MEMBRANE_QUAD4::registration;
FESystemElem::FESystemElemTypeEnum::EnumRegistration<FESystemElem::STRUCTURAL_MEMBRANE_TRI3> 
FESystemElem::STRUCTURAL_MEMBRANE_TRI3::registration;
FESystemElem::FESystemElemTypeEnum::EnumRegistration<FESystemElem::STRUCTURAL_PLATE_MITC4_QUAD4> 
FESystemElem::STRUCTURAL_PLATE_MITC4_QUAD4::registration;
FESystemElem::FESystemElemTypeEnum::EnumRegistration<FESystemElem::STRUCTURAL_PLATE_DKT> 
FESystemElem::STRUCTURAL_PLATE_DKT::registration;
FESystemElem::FESystemElemTypeEnum::EnumRegistration<FESystemElem::STRUCTURAL_TRI3_VON_KARMAN> 
FESystemElem::STRUCTURAL_TRI3_VON_KARMAN::registration;


FESystemElem::PlateDKBatozQtyEnum::EnumRegistration<FESystemElem::PLATE_DK_BATOZ_HX_HX_FACTOR> 
FESystemElem::PLATE_DK_BATOZ_HX_HX_FACTOR::registration;
FESystemElem::PlateDKBatozQtyEnum::EnumRegistration<FESystemElem::PLATE_DK_BATOZ_HY_HY_FACTOR> 
FESystemElem::PLATE_DK_BATOZ_HY_HY_FACTOR::registration;
FESystemElem::PlateDKBatozQtyEnum::EnumRegistration<FESystemElem::PLATE_DK_BATOZ_HXX_HXX_FACTOR> 
FESystemElem::PLATE_DK_BATOZ_HXX_HXX_FACTOR::registration;
FESystemElem::PlateDKBatozQtyEnum::EnumRegistration<FESystemElem::PLATE_DK_BATOZ_HYY_HYY_FACTOR> 
FESystemElem::PLATE_DK_BATOZ_HYY_HYY_FACTOR::registration;
FESystemElem::PlateDKBatozQtyEnum::EnumRegistration<FESystemElem::PLATE_DK_BATOZ_HXX_HYY_FACTOR> 
FESystemElem::PLATE_DK_BATOZ_HXX_HYY_FACTOR::registration;
FESystemElem::PlateDKBatozQtyEnum::EnumRegistration<FESystemElem::PLATE_DK_BATOZ_HXY_HXY_FACTOR> 
FESystemElem::PLATE_DK_BATOZ_HXY_HXY_FACTOR::registration;
FESystemElem::PlateDKBatozQtyEnum::EnumRegistration<FESystemElem::PLATE_DK_BATOZ_HYX_HXY_FACTOR> 
FESystemElem::PLATE_DK_BATOZ_HYX_HXY_FACTOR::registration;
FESystemElem::PlateDKBatozQtyEnum::EnumRegistration<FESystemElem::PLATE_DK_BATOZ_HYX_HYX_FACTOR> 
FESystemElem::PLATE_DK_BATOZ_HYX_HYX_FACTOR::registration;

FESystemElem::PlateVonKarmanQtyEnum::EnumRegistration<FESystemElem::PLATE_VON_KARMAN_NONLINEAR_STIFFNESS_COMPONENT> 
FESystemElem::PLATE_VON_KARMAN_NONLINEAR_STIFFNESS_COMPONENT::registration;



// ******************  Thermal Elements ************************* 
NameEnumerationHandler FESystemElem::ThermalElemQtyEnum::enum_handler;

FESystemElem::FESystemElemTypeEnum::EnumRegistration<FESystemElem::THERMAL_CONDUCTION_EDGE2> 
FESystemElem::THERMAL_CONDUCTION_EDGE2::registration;
FESystemElem::FESystemElemTypeEnum::EnumRegistration<FESystemElem::THERMAL_CONDUCTION_HEX8> 
FESystemElem::THERMAL_CONDUCTION_HEX8::registration;
FESystemElem::FESystemElemTypeEnum::EnumRegistration<FESystemElem::THERMAL_CONDUCTION_PRISM6> 
FESystemElem::THERMAL_CONDUCTION_PRISM6::registration;
FESystemElem::FESystemElemTypeEnum::EnumRegistration<FESystemElem::THERMAL_CONDUCTION_QUAD4> 
FESystemElem::THERMAL_CONDUCTION_QUAD4::registration;
FESystemElem::FESystemElemTypeEnum::EnumRegistration<FESystemElem::THERMAL_CONDUCTION_TRI3> 
FESystemElem::THERMAL_CONDUCTION_TRI3::registration;
FESystemElem::FESystemElemTypeEnum::EnumRegistration<FESystemElem::THERMAL_FORCED_CONVECTION_EDGE2> 
FESystemElem::THERMAL_FORCED_CONVECTION_EDGE2::registration;

FESystemElem::ThermalElemQtyEnum::EnumRegistration<FESystemElem::THERMAL_C_MATRIX> 
FESystemElem::THERMAL_C_MATRIX::registration;
FESystemElem::ThermalElemQtyEnum::EnumRegistration<FESystemElem::THERMAL_K_C_MATRIX> 
FESystemElem::THERMAL_K_C_MATRIX::registration;
FESystemElem::ThermalElemQtyEnum::EnumRegistration<FESystemElem::THERMAL_K_H_MATRIX> 
FESystemElem::THERMAL_K_H_MATRIX::registration;
FESystemElem::ThermalElemQtyEnum::EnumRegistration<FESystemElem::THERMAL_F_H_VECTOR> 
FESystemElem::THERMAL_F_H_VECTOR::registration;
FESystemElem::ThermalElemQtyEnum::EnumRegistration<FESystemElem::THERMAL_F_VOL_VECTOR> 
FESystemElem::THERMAL_F_VOL_VECTOR::registration;
FESystemElem::ThermalElemQtyEnum::EnumRegistration<FESystemElem::THERMAL_F_Q_SURF_VECTOR> 
FESystemElem::THERMAL_F_Q_SURF_VECTOR::registration;
FESystemElem::ThermalElemQtyEnum::EnumRegistration<FESystemElem::THERMAL_F_EMITTED_RAD_VECTOR> 
FESystemElem::THERMAL_F_EMITTED_RAD_VECTOR::registration;
FESystemElem::ThermalElemQtyEnum::EnumRegistration<FESystemElem::THERMAL_EMITTED_RAD_JAC_MATRIX> 
FESystemElem::THERMAL_EMITTED_RAD_JAC_MATRIX::registration;
FESystemElem::ThermalElemQtyEnum::EnumRegistration<FESystemElem::THERMAL_F_SIGMA_FACTOR> 
FESystemElem::THERMAL_F_SIGMA_FACTOR::registration;
FESystemElem::ThermalElemQtyEnum::EnumRegistration<FESystemElem::THERMAL_K_C_JAC_MATRIX> 
FESystemElem::THERMAL_K_C_JAC_MATRIX::registration;
FESystemElem::ThermalElemQtyEnum::EnumRegistration<FESystemElem::THERMAL_C_JAC_MATRIX> 
FESystemElem::THERMAL_C_JAC_MATRIX::registration;



// ********************* Interpolation Element ******************* 
NameEnumerationHandler FEInterpolationElemQtyEnum::enum_handler;

FESystemElem::FESystemElemTypeEnum::EnumRegistration<FESystemElem::FE_INTERPOLATION_EDGE2> 
FESystemElem::FE_INTERPOLATION_EDGE2::registration;
FESystemElem::FESystemElemTypeEnum::EnumRegistration<FESystemElem::FE_INTERPOLATION_HEX8> 
FESystemElem::FE_INTERPOLATION_HEX8::registration;
FESystemElem::FESystemElemTypeEnum::EnumRegistration<FESystemElem::FE_INTERPOLATION_QUAD4> 
FESystemElem::FE_INTERPOLATION_QUAD4::registration;
FESystemElem::FESystemElemTypeEnum::EnumRegistration<FESystemElem::FE_INTERPOLATION_TRI3> 
FESystemElem::FE_INTERPOLATION_TRI3::registration;



// ****************** Property ************************************
MaterialPropertyNameEnumerationHandler MaterialPropertyNameBase::enum_handler;

ElemDataCardFactorsEnum::EnumRegistration<MASS_FACTOR> MASS_FACTOR::registration;
ElemDataCardFactorsEnum::EnumRegistration<MEMBRANE_MASS_FACTOR> MEMBRANE_MASS_FACTOR::registration;
ElemDataCardFactorsEnum::EnumRegistration<PLATE_MASS_FACTOR> PLATE_MASS_FACTOR::registration;
ElemDataCardFactorsEnum::EnumRegistration<THERMAL_EXPANSION_FACTOR> 
THERMAL_EXPANSION_FACTOR::registration;
ElemDataCardFactorsEnum::EnumRegistration<THERMAL_CAPACITANCE_FACTOR> 
THERMAL_CAPACITANCE_FACTOR::registration;
ElemDataCardFactorsEnum::EnumRegistration<THERMAL_CONDUCTANCE_FACTOR> 
THERMAL_CONDUCTANCE_FACTOR::registration;
ElemDataCardFactorsEnum::EnumRegistration<THERMAL_EMITTED_LOAD_FACTOR> 
THERMAL_EMITTED_LOAD_FACTOR::registration;
ElemDataCardFactorsEnum::EnumRegistration<STIFFNESS_A_MATRIX_FACTOR> 
STIFFNESS_A_MATRIX_FACTOR::registration;
ElemDataCardFactorsEnum::EnumRegistration<STIFFNESS_B_MATRIX_FACTOR> 
STIFFNESS_B_MATRIX_FACTOR::registration;
ElemDataCardFactorsEnum::EnumRegistration<STIFFNESS_D_MATRIX_FACTOR> 
STIFFNESS_D_MATRIX_FACTOR::registration;
ElemDataCardFactorsEnum::EnumRegistration<RADIATION_EPSILON_FACTOR_1> 
RADIATION_EPSILON_FACTOR_1::registration;
ElemDataCardFactorsEnum::EnumRegistration<RADIATION_EPSILON_FACTOR_2> 
RADIATION_EPSILON_FACTOR_2::registration;
ElemDataCardFactorsEnum::EnumRegistration<STRESS_STRAIN_FACTOR> STRESS_STRAIN_FACTOR::registration;

MaterialPropertyNameBase::EnumRegistration<DENSITY> DENSITY::registration;
MaterialPropertyNameBase::EnumRegistration<SPECIFIC_HEAT> SPECIFIC_HEAT::registration;
MaterialPropertyNameBase::EnumRegistration<ALPHA_EXPANSION> ALPHA_EXPANSION::registration;
MaterialPropertyNameBase::EnumRegistration<THERMAL_CONDUCTIVITY> THERMAL_CONDUCTIVITY::registration;
MaterialPropertyNameBase::EnumRegistration<YOUNGS_MODULUS> YOUNGS_MODULUS::registration;
MaterialPropertyNameBase::EnumRegistration<POISSONS_RATIO> POISSONS_RATIO::registration;
MaterialPropertyNameBase::EnumRegistration<SPRING_CONSTANT> SPRING_CONSTANT::registration;
MaterialPropertyNameBase::EnumRegistration<EMISSIVITY> EMISSIVITY::registration;
//MaterialPropertyNameBase::EnumRegistration<TEMP_REF> TEMP_REF::registration;

MaterialPropertyCardEnum::EnumRegistration<ISOTROPIC_MATERIAL_PROPERTY_CARD> 
ISOTROPIC_MATERIAL_PROPERTY_CARD::registration;

ElemDataCardEnum::EnumRegistration<ISOTROPIC_1D_ELEM_DATA_CARD> 
ISOTROPIC_1D_ELEM_DATA_CARD::registration;
ElemDataCardEnum::EnumRegistration<ISOTROPIC_2D_ELEM_DATA_CARD> 
ISOTROPIC_2D_ELEM_DATA_CARD::registration;
ElemDataCardEnum::EnumRegistration<ISOTROPIC_3D_ELEM_DATA_CARD> 
ISOTROPIC_3D_ELEM_DATA_CARD::registration;


ElemDataCardFactorsEnum::EnumRegistration<EA_FACTOR> EA_FACTOR::registration;
ElemDataCardFactorsEnum::EnumRegistration<GA_FACTOR> GA_FACTOR::registration;
ElemDataCardFactorsEnum::EnumRegistration<EIZZ_FACTOR> EIZZ_FACTOR::registration;
ElemDataCardFactorsEnum::EnumRegistration<EIYY_FACTOR> EIYY_FACTOR::registration;
ElemDataCardFactorsEnum::EnumRegistration<SPRING_STIFFNESS_FACTOR> 
SPRING_STIFFNESS_FACTOR::registration;
ElemDataCardFactorsEnum::EnumRegistration<SOLID_STIFFNESS_MATRIX_FACTOR> 
SOLID_STIFFNESS_MATRIX_FACTOR::registration;

Isotropic_1D_ElemDataEnum::EnumRegistration<AREA> AREA::registration;
Isotropic_1D_ElemDataEnum::EnumRegistration<THICKNESS_1D_ELEM> THICKNESS_1D_ELEM::registration;
Isotropic_1D_ElemDataEnum::EnumRegistration<IZZ> IZZ::registration;
Isotropic_1D_ElemDataEnum::EnumRegistration<IYY> IYY::registration;
Isotropic_1D_ElemDataEnum::EnumRegistration<ROTATION_AXIS_X> ROTATION_AXIS_X::registration;
Isotropic_1D_ElemDataEnum::EnumRegistration<ROTATION_AXIS_Y> ROTATION_AXIS_Y::registration;
Isotropic_1D_ElemDataEnum::EnumRegistration<ROTATION_AXIS_Z> ROTATION_AXIS_Z::registration;

Isotropic_2D_ElemDataEnum::EnumRegistration<THICKNESS_2D_ELEM> 
THICKNESS_2D_ELEM::registration;


Property::PropertyCardParameterType::EnumRegistration<Property::LOCAL_PROPERTY_PARAMETER> 
Property::LOCAL_PROPERTY_PARAMETER::registration;
Property::PropertyCardParameterType::EnumRegistration<Property::GLOBAL_PROPERTY_PARAMETER> 
Property::GLOBAL_PROPERTY_PARAMETER::registration;

Property::LocalParameterType::EnumRegistration<Property::TEMPERATURE> 
Property::TEMPERATURE::registration;


// ****************** Property ************************************

ElemDataCardEnum::EnumRegistration<FORCED_CONVECTION_1D_ELEM_DATA_CARD> 
FORCED_CONVECTION_1D_ELEM_DATA_CARD::registration;

NameEnumerationHandler ForcedConvection_1D_ElemDataEnum::enum_handler;
NameEnumerationHandler ForcedConvection_1D_ElemDataCardFactorsEnum::enum_handler;

ForcedConvection_1D_ElemDataEnum::EnumRegistration<WALL_AREA>  WALL_AREA::registration; 
ForcedConvection_1D_ElemDataEnum::EnumRegistration<FLUID_AREA>  FLUID_AREA::registration; 
ForcedConvection_1D_ElemDataEnum::EnumRegistration<CONTACT_PERIMETER>  
CONTACT_PERIMETER::registration; 
ForcedConvection_1D_ElemDataEnum::EnumRegistration<MASS_FLOW_RATE>  MASS_FLOW_RATE::registration; 
ForcedConvection_1D_ElemDataEnum::EnumRegistration<CONVECTION_COEFFICIENT>  
CONVECTION_COEFFICIENT::registration; 


ForcedConvection_1D_ElemDataCardFactorsEnum::EnumRegistration<FLUID_CAPACITANCE_FACTOR>  
FLUID_CAPACITANCE_FACTOR::registration; 
ForcedConvection_1D_ElemDataCardFactorsEnum::EnumRegistration<WALL_CAPACITANCE_FACTOR>  
WALL_CAPACITANCE_FACTOR::registration; 
ForcedConvection_1D_ElemDataCardFactorsEnum::EnumRegistration<FLUID_CONVECTIVE_FACTOR>  
FLUID_CONVECTIVE_FACTOR::registration; 
ForcedConvection_1D_ElemDataCardFactorsEnum::EnumRegistration<CONVECTIVE_EXCANGE_FACTOR>  
CONVECTIVE_EXCANGE_FACTOR::registration; 
ForcedConvection_1D_ElemDataCardFactorsEnum::EnumRegistration<WALL_CONDUCTIVITY_FACTOR>  
WALL_CONDUCTIVITY_FACTOR::registration; 
ForcedConvection_1D_ElemDataCardFactorsEnum::EnumRegistration<FLUID_CONDUCTIVITY_FACTOR>  
FLUID_CONDUCTIVITY_FACTOR::registration; 

