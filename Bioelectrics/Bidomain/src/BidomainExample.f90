!> \file
!> \author Chris Bradley
!> \brief This is an example program to solve a Bidomain equation using OpenCMISS calls.
!>
!> \section LICENSE
!>
!> Version: MPL 1.1/GPL 2.0/LGPL 2.1
!>
!> The contents of this file are subject to the Mozilla Public License
!> Version 1.1 (the "License"); you may not use this file except in
!> compliance with the License. You may obtain a copy of the License at
!> http://www.mozilla.org/MPL/
!>
!> Software distributed under the License is distributed on an "AS IS"
!> basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See the
!> License for the specific language governing rights and limitations
!> under the License.
!>
!> The Original Code is OpenCMISS
!>
!> The Initial Developer of the Original Code is University of Auckland,
!> Auckland, New Zealand and University of Oxford, Oxford, United
!> Kingdom. Portions created by the University of Auckland and University
!> of Oxford are Copyright (C) 2007 by the University of Auckland and
!> the University of Oxford. All Rights Reserved.
!>
!> Contributor(s):
!>
!> Alternatively, the contents of this file may be used under the terms of
!> either the GNU General Public License Version 2 or later (the "GPL"), or
!> the GNU Lesser General Public License Version 2.1 or later (the "LGPL"),
!> in which case the provisions of the GPL or the LGPL are applicable instead
!> of those above. If you wish to allow use of your version of this file only
!> under the terms of either the GPL or the LGPL, and not to allow others to
!> use your version of this file under the terms of the MPL, indicate your
!> decision by deleting the provisions above and replace them with the notice
!> and other provisions required by the GPL or the LGPL. If you do not delete
!> the provisions above, a recipient may use your version of this file under
!> the terms of any one of the MPL, the GPL or the LGPL.
!>

!> \example Bioelectrics/Bidomain/src/BidomainExample.f90
!! Example program to solve a Bidomain equation using openCMISS calls.
!! \par Latest Builds:
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/Bioelectrics/Bidomain/build-intel'>Linux Intel Build</a>
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/Bioelectrics/Bidomain/build-gnu'>Linux GNU Build</a>
!!
!<

!> Main program
PROGRAM BIDOMAINEXAMPLE

  USE OPENCMISS
  USE MPI

#ifdef WIN32
  USE IFQWIN
#endif

  IMPLICIT NONE



  !Test program parameters

  REAL(CMISSDP), PARAMETER :: HEIGHT=1.0_CMISSDP
  REAL(CMISSDP), PARAMETER :: WIDTH=2.0_CMISSDP
  REAL(CMISSDP), PARAMETER :: LENGTH=3.0_CMISSDP

  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystemUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: RegionUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: BasisUserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: GeneratedMeshUserNumber=4
  INTEGER(CMISSIntg), PARAMETER :: MeshUserNumber=5
  INTEGER(CMISSIntg), PARAMETER :: DecompositionUserNumber=6
  INTEGER(CMISSIntg), PARAMETER :: GeometricFieldUserNumber=7
  INTEGER(CMISSIntg), PARAMETER :: ParabolicEquationsSetFieldUserNumber=8
  INTEGER(CMISSIntg), PARAMETER :: EllipticEquationsSetFieldUserNumber=9
  INTEGER(CMISSIntg), PARAMETER :: DependentFieldUserNumber=10
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumber=11
  INTEGER(CMISSIntg), PARAMETER :: CellMLUserNumber=12
  INTEGER(CMISSIntg), PARAMETER :: CellMLModelsFieldUserNumber=13
  INTEGER(CMISSIntg), PARAMETER :: CellMLStateFieldUserNumber=14
  INTEGER(CMISSIntg), PARAMETER :: CellMLIntermediateFieldUserNumber=15
  INTEGER(CMISSIntg), PARAMETER :: CellMLParametersFieldUserNumber=16
  INTEGER(CMISSIntg), PARAMETER :: ParabolicEquationsSetUserNumber=17
  INTEGER(CMISSIntg), PARAMETER :: EllipticEquationsSetUserNumber=18
  INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=19

  !Program types
  
  !Program variables

  INTEGER(CMISSIntg) :: NUMBER_OF_ARGUMENTS,ARGUMENT_LENGTH,STATUS
  INTEGER(CMISSIntg) :: NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS,NUMBER_GLOBAL_Z_ELEMENTS, &
    & INTERPOLATION_TYPE,NUMBER_OF_GAUSS_XI
  CHARACTER(LEN=255) :: COMMAND_ARGUMENT,Filename

  LOGICAL :: EXPORT_FIELD

  INTEGER(CMISSIntg) :: N,CELL_TYPE

  INTEGER(CMISSIntg) :: n98ModelIndex

  INTEGER(CMISSIntg) :: gNacomponent,stimcomponent,node_idx

  REAL(CMISSDP) :: X,Y,DISTANCE,gK1_VALUE,gNa_VALUE

  INTEGER(CMISSIntg), PARAMETER :: NUMBER_OF_ELEMENTS=25

  REAL(CMISSDP), PARAMETER :: STIM_VALUE = 100.0_CMISSDP
  REAL(CMISSDP), PARAMETER :: STIM_STOP = 0.10_CMISSDP
  REAL(CMISSDP), PARAMETER :: TIME_STOP = 1.50_CMISSDP
  REAL(CMISSDP), PARAMETER :: ODE_TIME_STEP = 0.00001_CMISSDP
  REAL(CMISSDP), PARAMETER :: PDE_TIME_STEP = 0.001_CMISSDP
  REAL(CMISSDP), PARAMETER :: CONDUCTIVITY_I = 0.1_CMISSDP
  REAL(CMISSDP), PARAMETER :: CONDUCTIVITY_E = 0.2_CMISSDP

  !CMISS variables

  TYPE(CMISSBasisType) :: Basis
  TYPE(CMISSBoundaryConditionsType) :: BoundaryConditions
  TYPE(CMISSCellMLType) :: CellML
  TYPE(CMISSCellMLEquationsType) :: CellMLEquations
  TYPE(CMISSControlLoopType) :: ControlLoop
  TYPE(CMISSCoordinateSystemType) :: CoordinateSystem,WorldCoordinateSystem
  TYPE(CMISSDecompositionType) :: Decomposition
  TYPE(CMISSEquationsType) :: Equations
  TYPE(CMISSEquationsSetType) :: ParabolicEquationsSet,EllipticEquationsSet
  TYPE(CMISSFieldType) :: GeometricField,ParabolicEquationsSetField,EllipticEquationsSetField,DependentField, &
    & MaterialsField
  TYPE(CMISSFieldType) :: CellMLModelsField,CellMLStateField,CellMLIntermediateField,CellMLParametersField
  TYPE(CMISSFieldsType) :: Fields
  TYPE(CMISSGeneratedMeshType) :: GeneratedMesh  
  TYPE(CMISSMeshType) :: Mesh
  TYPE(CMISSProblemType) :: Problem
  TYPE(CMISSRegionType) :: Region,WorldRegion
  TYPE(CMISSSolverType) :: Solver
  TYPE(CMISSSolverEquationsType) :: SolverEquations

#ifdef WIN32
  !Quickwin type
  LOGICAL :: QUICKWIN_STATUS=.FALSE.
  TYPE(WINDOWCONFIG) :: QUICKWIN_WINDOW_CONFIG
#endif
  
  !Generic CMISS variables
  
  INTEGER(CMISSIntg) :: NumberOfComputationalNodes,ComputationalNodeNumber
  INTEGER(CMISSIntg) :: ParabolicEquationsSetIndex,EllipticEquationsSetIndex,CellMLIndex
  INTEGER(CMISSIntg) :: FirstNodeNumber,LastNodeNumber
  INTEGER(CMISSIntg) :: FirstNodeDomain,LastNodeDomain,NodeDomain
  INTEGER(CMISSIntg) :: Err
    
#ifdef WIN32
  !Initialise QuickWin
  QUICKWIN_WINDOW_CONFIG%TITLE="General Output" !Window title
  QUICKWIN_WINDOW_CONFIG%NUMTEXTROWS=-1 !Max possible number of rows
  QUICKWIN_WINDOW_CONFIG%MODE=QWIN$SCROLLDOWN
  !Set the window parameters
  QUICKWIN_STATUS=SETWINDOWCONFIG(QUICKWIN_WINDOW_CONFIG)
  !If attempt fails set with system estimated values
  IF(.NOT.QUICKWIN_STATUS) QUICKWIN_STATUS=SETWINDOWCONFIG(QUICKWIN_WINDOW_CONFIG)
#endif

  !Intialise OpenCMISS
  CALL CMISSInitialise(WorldCoordinateSystem,WorldRegion,Err)

  CALL CMISSErrorHandlingModeSet(CMISSTrapError,Err)

  !Get the computational nodes information
  CALL CMISSComputationalNumberOfNodesGet(NumberOfComputationalNodes,Err)
  CALL CMISSComputationalNodeNumberGet(ComputationalNodeNumber,Err)
    
  NUMBER_OF_ARGUMENTS = COMMAND_ARGUMENT_COUNT()
  IF(NUMBER_OF_ARGUMENTS >= 4) THEN
    !If we have enough arguments then use the first four for setting up the problem. The subsequent arguments may be used to
    !pass flags to, say, PETSc.
    CALL GET_COMMAND_ARGUMENT(1,COMMAND_ARGUMENT,ARGUMENT_LENGTH,STATUS)
    IF(STATUS>0) CALL HANDLE_ERROR("Error for command argument 1.")
    READ(COMMAND_ARGUMENT(1:ARGUMENT_LENGTH),*) NUMBER_GLOBAL_X_ELEMENTS
    IF(NUMBER_GLOBAL_X_ELEMENTS<=0) CALL HANDLE_ERROR("Invalid number of X elements.")
    CALL GET_COMMAND_ARGUMENT(2,COMMAND_ARGUMENT,ARGUMENT_LENGTH,STATUS)
    IF(STATUS>0) CALL HANDLE_ERROR("Error for command argument 2.")
    READ(COMMAND_ARGUMENT(1:ARGUMENT_LENGTH),*) NUMBER_GLOBAL_Y_ELEMENTS
    IF(NUMBER_GLOBAL_Y_ELEMENTS<=0) CALL HANDLE_ERROR("Invalid number of Y elements.")
    CALL GET_COMMAND_ARGUMENT(3,COMMAND_ARGUMENT,ARGUMENT_LENGTH,STATUS)
    IF(STATUS>0) CALL HANDLE_ERROR("Error for command argument 3.")
    READ(COMMAND_ARGUMENT(1:ARGUMENT_LENGTH),*) NUMBER_GLOBAL_Z_ELEMENTS
    IF(NUMBER_GLOBAL_Z_ELEMENTS<0) CALL HANDLE_ERROR("Invalid number of Z elements.")
    CALL GET_COMMAND_ARGUMENT(4,COMMAND_ARGUMENT,ARGUMENT_LENGTH,STATUS)
    IF(STATUS>0) CALL HANDLE_ERROR("Error for command argument 4.")
    READ(COMMAND_ARGUMENT(1:ARGUMENT_LENGTH),*) INTERPOLATION_TYPE
    IF(INTERPOLATION_TYPE<=0) CALL HANDLE_ERROR("Invalid Interpolation specification.")
  ELSE
    !If there are not enough arguments default the problem specification 
    NUMBER_GLOBAL_X_ELEMENTS=2
    NUMBER_GLOBAL_Y_ELEMENTS=2
    NUMBER_GLOBAL_Z_ELEMENTS=0
    INTERPOLATION_TYPE=1
  ENDIF
  
  NUMBER_GLOBAL_X_ELEMENTS=NUMBER_OF_ELEMENTS
  NUMBER_GLOBAL_Y_ELEMENTS=NUMBER_OF_ELEMENTS
  NUMBER_GLOBAL_Z_ELEMENTS=0
  INTERPOLATION_TYPE=1
  
  !Start the creation of a new RC coordinate system
  CALL CMISSCoordinateSystemTypeInitialise(CoordinateSystem,Err)
  CALL CMISSCoordinateSystemCreateStart(CoordinateSystemUserNumber,CoordinateSystem,Err)
  IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
    !Set the coordinate system to be 2D
    CALL CMISSCoordinateSystemDimensionSet(CoordinateSystem,2,Err)
  ELSE
    !Set the coordinate system to be 3D
    CALL CMISSCoordinateSystemDimensionSet(CoordinateSystem,3,Err)
  ENDIF
  !Finish the creation of the coordinate system
  CALL CMISSCoordinateSystemCreateFinish(CoordinateSystem,Err)

  !Start the creation of the region
  CALL CMISSRegionTypeInitialise(Region,Err)
  CALL CMISSRegionCreateStart(RegionUserNumber,WorldRegion,Region,Err)
  CALL CMISSRegionLabelSet(Region,"BidomainRegion",Err)
  !Set the regions coordinate system to the RC coordinate system that we have created
  CALL CMISSRegionCoordinateSystemSet(Region,CoordinateSystem,Err)
  !Finish the creation of the region
  CALL CMISSRegionCreateFinish(Region,Err)

  !Start the creation of a basis (default is trilinear lagrange)
  CALL CMISSBasisTypeInitialise(Basis,Err)
  CALL CMISSBasisCreateStart(BasisUserNumber,Basis,Err)
  SELECT CASE(INTERPOLATION_TYPE)
  CASE(1,2,3,4)
    CALL CMISSBasisTypeSet(Basis,CMISSBasisLagrangeHermiteTPType,Err)
  CASE(7,8,9)
    CALL CMISSBasisTypeSet(Basis,CMISSBasisSimplexType,Err)
  CASE DEFAULT
    CALL HANDLE_ERROR("Invalid interpolation type.")
  END SELECT
  SELECT CASE(INTERPOLATION_TYPE)
  CASE(1)
    NUMBER_OF_GAUSS_XI=2
  CASE(2)
    NUMBER_OF_GAUSS_XI=3
  CASE(3,4)
    NUMBER_OF_GAUSS_XI=4
  CASE DEFAULT
    NUMBER_OF_GAUSS_XI=0 !Don't set number of Gauss points for tri/tet
  END SELECT
  IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
    !Set the basis to be a bi-interpolation basis
    CALL CMISSBasisNumberOfXiSet(Basis,2,Err)
    CALL CMISSBasisInterpolationXiSet(Basis,[INTERPOLATION_TYPE,INTERPOLATION_TYPE],Err)
    IF(NUMBER_OF_GAUSS_XI>0) THEN
      CALL CMISSBasisQuadratureNumberOfGaussXiSet(Basis,[NUMBER_OF_GAUSS_XI,NUMBER_OF_GAUSS_XI],Err)
    ENDIF
  ELSE
    !Set the basis to be a tri-interpolation basis
    CALL CMISSBasisNumberOfXiSet(Basis,3,Err)
    CALL CMISSBasisInterpolationXiSet(Basis,[INTERPOLATION_TYPE,INTERPOLATION_TYPE,INTERPOLATION_TYPE],Err)
    IF(NUMBER_OF_GAUSS_XI>0) THEN
      CALL CMISSBasisQuadratureNumberOfGaussXiSet(Basis,[NUMBER_OF_GAUSS_XI,NUMBER_OF_GAUSS_XI,NUMBER_OF_GAUSS_XI],Err)
    ENDIF
  ENDIF
  !Finish the creation of the basis
  CALL CMISSBasisCreateFinish(Basis,Err)

  !Start the creation of a generated mesh in the region
  CALL CMISSGeneratedMeshTypeInitialise(GeneratedMesh,Err)
  CALL CMISSGeneratedMeshCreateStart(GeneratedMeshUserNumber,Region,GeneratedMesh,Err)
  !Set up a regular x*y*z mesh
  CALL CMISSGeneratedMeshTypeSet(GeneratedMesh,CMISSGeneratedMeshRegularMeshType,Err)
  !Set the default basis
  CALL CMISSGeneratedMeshBasisSet(GeneratedMesh,Basis,Err)   
  !Define the mesh on the region
  IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
    CALL CMISSGeneratedMeshExtentSet(GeneratedMesh,[WIDTH,HEIGHT],Err)
    CALL CMISSGeneratedMeshNumberOfElementsSet(GeneratedMesh,[NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS],Err)
  ELSE
    CALL CMISSGeneratedMeshExtentSet(GeneratedMesh,[WIDTH,HEIGHT,LENGTH],Err)
    CALL CMISSGeneratedMeshNumberOfElementsSet(GeneratedMesh,[NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS, &
      & NUMBER_GLOBAL_Z_ELEMENTS],Err)
  ENDIF    
  !Finish the creation of a generated mesh in the region
  CALL CMISSMeshTypeInitialise(Mesh,Err)
  CALL CMISSGeneratedMeshCreateFinish(GeneratedMesh,MeshUserNumber,Mesh,Err)
  !Create a decomposition
  CALL CMISSDecompositionTypeInitialise(Decomposition,Err)
  CALL CMISSDecompositionCreateStart(DecompositionUserNumber,Mesh,Decomposition,Err)
  !Set the decomposition to be a general decomposition with the specified number of domains
  CALL CMISSDecompositionTypeSet(Decomposition,CMISSDecompositionCalculatedType,Err)
  CALL CMISSDecompositionNumberOfDomainsSet(Decomposition,NumberOfComputationalNodes,Err)
  !Finish the decomposition
  CALL CMISSDecompositionCreateFinish(Decomposition,Err)
  
  !Start to create a default (geometric) field on the region
  CALL CMISSFieldTypeInitialise(GeometricField,Err)
  CALL CMISSFieldCreateStart(GeometricFieldUserNumber,Region,GeometricField,Err)
  !Set the decomposition to use
  CALL CMISSFieldMeshDecompositionSet(GeometricField,Decomposition,Err)
  !Set the domain to be used by the field components.
  CALL CMISSFieldComponentMeshComponentSet(GeometricField,CMISSFieldUVariableType,1,1,Err)
  CALL CMISSFieldComponentMeshComponentSet(GeometricField,CMISSFieldUVariableType,2,1,Err)
  IF(NUMBER_GLOBAL_Z_ELEMENTS/=0) THEN
    CALL CMISSFieldComponentMeshComponentSet(GeometricField,CMISSFieldUVariableType,3,1,Err)
  ENDIF
  !Finish creating the field
  CALL CMISSFieldCreateFinish(GeometricField,Err)

  !Update the geometric field parameters
  CALL CMISSGeneratedMeshGeometricParametersCalculate(GeometricField,GeneratedMesh,Err)
        
  !Create the first (parabolic) equations set for the bidomain equations
  CALL CMISSEquationsSetTypeInitialise(ParabolicEquationsSet,Err)
  CALL CMISSFieldTypeInitialise(ParabolicEquationsSetField,Err)
  CALL CMISSEquationsSetCreateStart(ParabolicEquationsSetUserNumber,Region,GeometricField,CMISSEquationsSetBioelectricsClass, &
    & CMISSEquationsSetBidomainEquationType,CMISSEquationsSetFirstBidomainSubtype,ParabolicEquationsSetFieldUserNumber, &
    & ParabolicEquationsSetField,ParabolicEquationsSet,Err)
  !Finish creating the equations set
  CALL CMISSEquationsSetCreateFinish(ParabolicEquationsSet,Err)

  !Create the equations set dependent field variables
  CALL CMISSFieldTypeInitialise(DependentField,Err)
  CALL CMISSEquationsSetDependentCreateStart(ParabolicEquationsSet,DependentFieldUserNumber,DependentField,Err)
  !Finish the equations set dependent field variables
  CALL CMISSEquationsSetDependentCreateFinish(ParabolicEquationsSet,Err)
  
  !Create the equations set materials field variables
  CALL CMISSFieldTypeInitialise(MaterialsField,Err)
  CALL CMISSEquationsSetMaterialsCreateStart(ParabolicEquationsSet,MaterialsFieldUserNumber,MaterialsField,Err)
  !Finish the equations set materials field variables
  CALL CMISSEquationsSetMaterialsCreateFinish(ParabolicEquationsSet,Err)
  
  !Create the second (elliptic) equations set for the bidomain equations
  CALL CMISSEquationsSetTypeInitialise(EllipticEquationsSet,Err)
  CALL CMISSEquationsSetCreateStart(EllipticEquationsSetUserNumber,Region,GeometricField,CMISSEquationsSetBioelectricsClass, &
    & CMISSEquationsSetBidomainEquationType,CMISSEquationsSetSecondBidomainSubtype,EllipticEquationsSetFieldUserNumber, &
    & EllipticEquationsSetField,EllipticEquationsSet,Err)
  !Finish creating the equations set
  CALL CMISSEquationsSetCreateFinish(EllipticEquationsSet,Err)

  !Create the equations set dependent field variables using the field from the first (parabolic) equations set
  CALL CMISSEquationsSetDependentCreateStart(EllipticEquationsSet,DependentFieldUserNumber,DependentField,Err)
  !Finish the equations set dependent field variables
  CALL CMISSEquationsSetDependentCreateFinish(EllipticEquationsSet,Err)
  
  !Create the equations set materials field variables using the field from the first (parabolic) equations set
  CALL CMISSEquationsSetMaterialsCreateStart(EllipticEquationsSet,MaterialsFieldUserNumber,MaterialsField,Err)
  !Finish the equations set materials field variables
  CALL CMISSEquationsSetMaterialsCreateFinish(EllipticEquationsSet,Err)
  
  !Set Am
  CALL CMISSFieldComponentValuesInitialise(MaterialsField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,193.6_CMISSDP,Err)
  !Set Cm
  CALL CMISSFieldComponentValuesInitialise(MaterialsField,CMISSFieldUVariableType,CMISSFieldValuesSetType,2,0.014651_CMISSDP,Err)
  !Set conductivity
  CALL CMISSFieldComponentValuesInitialise(MaterialsField,CMISSFieldUVariableType,CMISSFieldValuesSetType,3,CONDUCTIVITY_I,Err)
  CALL CMISSFieldComponentValuesInitialise(MaterialsField,CMISSFieldUVariableType,CMISSFieldValuesSetType,4,CONDUCTIVITY_I,Err)
  IF(NUMBER_GLOBAL_Z_ELEMENTS/=0) THEN
    CALL CMISSFieldComponentValuesInitialise(MaterialsField,CMISSFieldUVariableType,CMISSFieldValuesSetType,5,CONDUCTIVITY_I,Err)
    CALL CMISSFieldComponentValuesInitialise(MaterialsField,CMISSFieldUVariableType,CMISSFieldValuesSetType,6,CONDUCTIVITY_E,Err)
    CALL CMISSFieldComponentValuesInitialise(MaterialsField,CMISSFieldUVariableType,CMISSFieldValuesSetType,7,CONDUCTIVITY_E,Err)
    CALL CMISSFieldComponentValuesInitialise(MaterialsField,CMISSFieldUVariableType,CMISSFieldValuesSetType,8,CONDUCTIVITY_E,Err)
  ELSE
    CALL CMISSFieldComponentValuesInitialise(MaterialsField,CMISSFieldUVariableType,CMISSFieldValuesSetType,5,CONDUCTIVITY_E,Err)
    CALL CMISSFieldComponentValuesInitialise(MaterialsField,CMISSFieldUVariableType,CMISSFieldValuesSetType,6,CONDUCTIVITY_E,Err)
  ENDIF
  
  !Create the CellML environment for the source field
  CALL CMISSCellMLTypeInitialise(CellML,Err)
  CALL CMISSCellMLCreateStart(CellMLUserNumber,Region,CellML,Err)
  !Import a Noble 1998 model from a file
  CALL CMISSCellMLModelImport(CellML,"n98.xml",n98ModelIndex,Err)
  ! Now we have imported all the models we are able to specify which variables from the model we want:
  !   - to set from this side
  CALL CMISSCellMLVariableSetAsKnown(CellML,n98ModelIndex,"fast_sodium_current/g_Na",Err)
  CALL CMISSCellMLVariableSetAsKnown(CellML,n98ModelIndex,"membrane/IStim",Err)
  !   - to get from the CellML side
  CALL CMISSCellMLVariableSetAsWanted(CellML,n98ModelIndex,"membrane/i_K1",Err)
  CALL CMISSCellMLVariableSetAsWanted(CellML,n98ModelIndex,"membrane/i_to",Err)
  CALL CMISSCellMLVariableSetAsWanted(CellML,n98ModelIndex,"membrane/i_K",Err)
  CALL CMISSCellMLVariableSetAsWanted(CellML,n98ModelIndex,"membrane/i_K_ATP",Err)
  CALL CMISSCellMLVariableSetAsWanted(CellML,n98ModelIndex,"membrane/i_Ca_L_K",Err)
  CALL CMISSCellMLVariableSetAsWanted(CellML,n98ModelIndex,"membrane/i_b_K",Err)
  CALL CMISSCellMLVariableSetAsWanted(CellML,n98ModelIndex,"membrane/i_NaK",Err)
  CALL CMISSCellMLVariableSetAsWanted(CellML,n98ModelIndex,"membrane/i_Na",Err)
  CALL CMISSCellMLVariableSetAsWanted(CellML,n98ModelIndex,"membrane/i_b_Na",Err)
  CALL CMISSCellMLVariableSetAsWanted(CellML,n98ModelIndex,"membrane/i_Ca_L_Na",Err)
  CALL CMISSCellMLVariableSetAsWanted(CellML,n98ModelIndex,"membrane/i_NaCa",Err)
  !Finish the CellML environment
  CALL CMISSCellMLCreateFinish(CellML,Err)

  !Start the creation of CellML <--> OpenCMISS field maps
  CALL CMISSCellMLFieldMapsCreateStart(CellML,Err)
  !Now  set up the field variable component <--> CellML model variable mappings.
  !Map Vm
  CALL CMISSCellMLCreateFieldToCellMLMap(CellML,DependentField,CMISSFieldUVariableType,1,CMISSFieldValuesSetType, &
    & n98ModelIndex,"membrane/V",CMISSFieldValuesSetType,Err)
  CALL CMISSCellMLCreateCellMLToFieldMap(CellML,n98ModelIndex,"membrane/V",CMISSFieldValuesSetType, &
    & DependentField,CMISSFieldUVariableType,1,CMISSFieldValuesSetType,Err)
  !Finish the creation of CellML <--> OpenCMISS field maps
  CALL CMISSCellMLFieldMapsCreateFinish(CellML,Err)

  !todo - get Vm initial value.
  CALL CMISSFieldComponentValuesInitialise(DependentField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,-92.5_CMISSDP,Err)

  !Start the creation of the CellML models field
  CALL CMISSFieldTypeInitialise(CellMLModelsField,Err)
  CALL CMISSCellMLModelsFieldCreateStart(CellMLModelsFieldUserNumber,CellML,CellMLModelsField,Err)
  !Finish the creation of the CellML models field
  CALL CMISSCellMLModelsFieldCreateFinish(CellML,Err)
  !Set up the models field
  !DO N=1,(NUMBER_GLOBAL_X_ELEMENTS+1)*(NUMBER_GLOBAL_Y_ELEMENTS+1)*(NUMBER_GLOBAL_Z_ELEMENTS+1)
  !  IF(N < 5) THEN
  !    CELL_TYPE = 1
  !  ELSE
  !    CELL_TYPE = 2
  !  ENDIF
  !  CALL CMISSFieldParameterSetUpdateNode(CellMLModelsField, CMISSFieldUVariableType, CMISSFieldValuesSetType,1,N,1,CELL_TYPE,Err)
  !END DO

  !Start the creation of the CellML state field
  CALL CMISSFieldTypeInitialise(CellMLStateField,Err)
  CALL CMISSCellMLStateFieldCreateStart(CellMLStateFieldUserNumber,CellML,CellMLStateField,Err)
  !Finish the creation of the CellML state field
  CALL CMISSCellMLStateFieldCreateFinish(CellML,Err)

  !Start the creation of the CellML intermediate field
  CALL CMISSFieldTypeInitialise(CellMLIntermediateField,Err)
  CALL CMISSCellMLIntermediateFieldCreateStart(CellMLIntermediateFieldUserNumber,CellML,CellMLIntermediateField,Err)
  !Finish the creation of the CellML intermediate field
  CALL CMISSCellMLIntermediateFieldCreateFinish(CellML,Err)
  
  !Start the creation of CellML parameters field
  CALL CMISSFieldTypeInitialise(CellMLParametersField,Err)
  CALL CMISSCellMLParametersFieldCreateStart(CellMLParametersFieldUserNumber,CellML,CellMLParametersField,Err)
  !Finish the creation of CellML parameters
  CALL CMISSCellMLParametersFieldCreateFinish(CellML,Err)
 
  !Create the parabolic equations set equations
  CALL CMISSEquationsTypeInitialise(Equations,Err)
  CALL CMISSEquationsSetEquationsCreateStart(ParabolicEquationsSet,Equations,Err)
  !Set the equations matrices sparsity type
  CALL CMISSEquationsSparsityTypeSet(Equations,CMISSEquationsSparseMatrices,Err)
  !Set the equations set output
  CALL CMISSEquationsOutputTypeSet(Equations,CMISSEquationsNoOutput,Err)
  !CALL CMISSEquationsOutputTypeSet(Equations,CMISSEquationsTimingOutput,Err)
  !CALL CMISSEquationsOutputTypeSet(Equations,CMISSEquationsMatrixOutput,Err)
  !CALL CMISSEquationsOutputTypeSet(Equations,CMISSEquationsElementMatrixOutput,Err)
  !Finish the equations set equations
  CALL CMISSEquationsSetEquationsCreateFinish(ParabolicEquationsSet,Err)

  !Create the elliptic equations set equations
  CALL CMISSEquationsTypeInitialise(Equations,Err)
  CALL CMISSEquationsSetEquationsCreateStart(EllipticEquationsSet,Equations,Err)
  !Set the equations matrices sparsity type
  CALL CMISSEquationsSparsityTypeSet(Equations,CMISSEquationsSparseMatrices,Err)
  !Set the equations set output
  CALL CMISSEquationsOutputTypeSet(Equations,CMISSEquationsNoOutput,Err)
  !CALL CMISSEquationsOutputTypeSet(Equations,CMISSEquationsTimingOutput,Err)
  !CALL CMISSEquationsOutputTypeSet(Equations,CMISSEquationsMatrixOutput,Err)
  !CALL CMISSEquationsOutputTypeSet(Equations,CMISSEquationsElementMatrixOutput,Err)
  !Finish the equations set equations
  CALL CMISSEquationsSetEquationsCreateFinish(EllipticEquationsSet,Err)

  CALL CMISSCellMLFieldComponentGet(CellML,n98ModelIndex,CMISSCellMLParametersFieldType,"membrane/IStim",stimcomponent,Err)
  !Set the Stimulus at half the bottom nodes
  DO node_idx=1,NUMBER_OF_ELEMENTS/2
    CALL CMISSDecompositionNodeDomainGet(Decomposition,node_idx,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL CMISSFieldParameterSetUpdateNode(CellMLParametersField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,node_idx, &
        & stimcomponent,STIM_VALUE,Err)
    ENDIF
  ENDDO

  !Set up the g_Na gradient
  CALL CMISSCellMLFieldComponentGet(CellML,n98ModelIndex,CMISSCellMLParametersFieldType,"fast_sodium_current/g_Na", &
    & gNacomponent,Err)
  !Loop over the nodes
  DO node_idx=1,(NUMBER_OF_ELEMENTS+1)*(NUMBER_OF_ELEMENTS+1)
    CALL CMISSFieldParameterSetGetNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,node_idx,1, &
      & X,Err)
    CALL CMISSFieldParameterSetGetNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,node_idx,2, &
      & Y,Err)
    DISTANCE=SQRT(X**2+Y**2)/SQRT(2.0_CMISSDP)
    gNa_VALUE=2.0_CMISSDP*(DISTANCE+0.5_CMISSDP)*385.5e-3_CMISSDP
    CALL CMISSFieldParameterSetUpdateNode(CellMLParametersField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,node_idx, &
      & gNacomponent,gNa_VALUE,Err)
  ENDDO
  
  !Start the creation of a problem.
  CALL CMISSProblemTypeInitialise(Problem,Err)
  CALL CMISSProblemCreateStart(ProblemUserNumber,Problem,Err)
  !Set the problem to be a standard Laplace problem
  CALL CMISSProblemSpecificationSet(Problem,CMISSProblemBioelectricsClass,CMISSProblemBidomainEquationType, &
    & CMISSProblemBidomainGudunovSplitSubtype,Err)
  !Finish the creation of a problem.
  CALL CMISSProblemCreateFinish(Problem,Err)

  !Start the creation of the problem control loop
  CALL CMISSProblemControlLoopCreateStart(Problem,Err)
  !Get the control loop
  CALL CMISSControlLoopTypeInitialise(ControlLoop,Err)
  CALL CMISSProblemControlLoopGet(Problem,CMISSControlLoopNode,ControlLoop,Err)
  !Set the times
  CALL CMISSControlLoopTimesSet(ControlLoop,0.0_CMISSDP,STIM_STOP,PDE_TIME_STEP,Err)
  !Set the output type
  !CALL CMISSControlLoopOutputTypeSet(ControlLoop,CMISSControlLoopNoOutput,Err)
  CALL CMISSControlLoopOutputTypeSet(ControlLoop,CMISSControlLoopProgressOutput,Err)
  !CALL CMISSControlLoopOutputTypeSet(ControlLoop,CMISSControlLoopTimingOutput,Err)
  !Finish creating the problem control loop
  CALL CMISSProblemControlLoopCreateFinish(Problem,Err)
 
  !Start the creation of the problem solvers
  CALL CMISSProblemSolversCreateStart(Problem,Err)
  !Get the first (DAE) solver
  CALL CMISSSolverTypeInitialise(Solver,Err)
  CALL CMISSProblemSolverGet(Problem,CMISSControlLoopNode,1,Solver,Err)
  !Set the DAE time step to by 10 us
  CALL CMISSSolverDAETimeStepSet(Solver,ODE_TIME_STEP,Err)
  !CALL CMISSSolverOutputTypeSet(Solver,CMISSSolverNoOutput,Err)
  CALL CMISSSolverOutputTypeSet(Solver,CMISSSolverProgressOutput,Err)
  !CALL CMISSSolverOutputTypeSet(Solver,CMISSSolverTimingOutput,Err)
  !CALL CMISSSolverOutputTypeSet(Solver,CMISSSolverSolverOutput,Err)
  !Get the second (Parabolic) solver
  CALL CMISSSolverTypeInitialise(Solver,Err)
  CALL CMISSProblemSolverGet(Problem,CMISSControlLoopNode,2,Solver,Err)
  !CALL CMISSSolverOutputTypeSet(Solver,CMISSSolverNoOutput,Err)
  CALL CMISSSolverOutputTypeSet(Solver,CMISSSolverProgressOutput,Err)
  !CALL CMISSSolverOutputTypeSet(Solver,CMISSSolverTimingOutput,Err)
  !CALL CMISSSolverOutputTypeSet(Solver,CMISSSolverSolverOutput,Err)
  !CALL CMISSSolverOutputTypeSet(Solver,CMISSSolverSolverMatrixOutput,Err)
  !Get the third (Elliptic) solver
  CALL CMISSSolverTypeInitialise(Solver,Err)
  CALL CMISSProblemSolverGet(Problem,CMISSControlLoopNode,3,Solver,Err)
  !CALL CMISSSolverOutputTypeSet(Solver,CMISSSolverNoOutput,Err)
  CALL CMISSSolverOutputTypeSet(Solver,CMISSSolverProgressOutput,Err)
  !CALL CMISSSolverOutputTypeSet(Solver,CMISSSolverTimingOutput,Err)
  !CALL CMISSSolverOutputTypeSet(Solver,CMISSSolverSolverOutput,Err)
  !CALL CMISSSolverOutputTypeSet(Solver,CMISSSolverSolverMatrixOutput,Err)
  !Finish the creation of the problem solver
  CALL CMISSProblemSolversCreateFinish(Problem,Err)

  !Start the creation of the problem solver CellML equations
  CALL CMISSProblemCellMLEquationsCreateStart(Problem,Err)
  !Get the first solver  
  !Get the CellML equations
  CALL CMISSSolverTypeInitialise(Solver,Err)
  CALL CMISSProblemSolverGet(Problem,CMISSControlLoopNode,1,Solver,Err)
  CALL CMISSCellMLEquationsTypeInitialise(CellMLEquations,Err)
  CALL CMISSSolverCellMLEquationsGet(Solver,CellMLEquations,Err)
  !Add in the CellML environement
  CALL CMISSCellMLEquationsCellMLAdd(CellMLEquations,CellML,CellMLIndex,Err)
  !Finish the creation of the problem solver CellML equations
  CALL CMISSProblemCellMLEquationsCreateFinish(Problem,Err)

  !Start the creation of the problem solver equations
  CALL CMISSProblemSolverEquationsCreateStart(Problem,Err)
  !Get the second (parabolic) solver  
  !Get the solver equations
  CALL CMISSSolverTypeInitialise(Solver,Err)
  CALL CMISSProblemSolverGet(Problem,CMISSControlLoopNode,2,Solver,Err)
  CALL CMISSSolverEquationsTypeInitialise(SolverEquations,Err)
  CALL CMISSSolverSolverEquationsGet(Solver,SolverEquations,Err)
  !Set the solver equations sparsity
  CALL CMISSSolverEquationsSparsityTypeSet(SolverEquations,CMISSSolverEquationsSparseMatrices,Err)
  !CALL CMISSSolverEquationsSparsityTypeSet(SolverEquations,CMISSSolverEquationsFullMatrices,Err)  
  !Add in the equations set
  CALL CMISSSolverEquationsEquationsSetAdd(SolverEquations,ParabolicEquationsSet,ParabolicEquationsSetIndex,Err)
  !Get the third (elliptic) solver  
  !Get the solver equations
  CALL CMISSSolverTypeInitialise(Solver,Err)
  CALL CMISSProblemSolverGet(Problem,CMISSControlLoopNode,3,Solver,Err)
  CALL CMISSSolverEquationsTypeInitialise(SolverEquations,Err)
  CALL CMISSSolverSolverEquationsGet(Solver,SolverEquations,Err)
  !Set the solver equations sparsity
  CALL CMISSSolverEquationsSparsityTypeSet(SolverEquations,CMISSSolverEquationsSparseMatrices,Err)
  !CALL CMISSSolverEquationsSparsityTypeSet(SolverEquations,CMISSSolverEquationsFullMatrices,Err)  
  !Add in the equations set
  CALL CMISSSolverEquationsEquationsSetAdd(SolverEquations,EllipticEquationsSet,EllipticEquationsSetIndex,Err)
  !Finish the creation of the problem solver equations
  CALL CMISSProblemSolverEquationsCreateFinish(Problem,Err)

  !Start the creation of the equations set boundary conditions
  CALL CMISSBoundaryConditionsTypeInitialise(BoundaryConditions,Err)
  CALL CMISSSolverTypeInitialise(Solver,Err)
  CALL CMISSProblemSolverGet(Problem,CMISSControlLoopNode,2,Solver,Err)
  CALL CMISSSolverEquationsTypeInitialise(SolverEquations,Err)
  CALL CMISSSolverSolverEquationsGet(Solver,SolverEquations,Err)
  CALL CMISSSolverEquationsBoundaryConditionsCreateStart(SolverEquations,BoundaryConditions,Err)
  !Set the first node to 0.0 and the last node to 1.0
  FirstNodeNumber=1
  IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
    LastNodeNumber=(NUMBER_GLOBAL_X_ELEMENTS+1)*(NUMBER_GLOBAL_Y_ELEMENTS+1)
  ELSE
    LastNodeNumber=(NUMBER_GLOBAL_X_ELEMENTS+1)*(NUMBER_GLOBAL_Y_ELEMENTS+1)*(NUMBER_GLOBAL_Z_ELEMENTS+1)
  ENDIF
  CALL CMISSDecompositionNodeDomainGet(Decomposition,FirstNodeNumber,1,FirstNodeDomain,Err)
  CALL CMISSDecompositionNodeDomainGet(Decomposition,LastNodeNumber,1,LastNodeDomain,Err)
  IF(FirstNodeDomain==ComputationalNodeNumber) THEN
    !CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField,CMISSFieldUVariableType,1,1,FirstNodeNumber,1, &
    !  & CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  ENDIF
  IF(LastNodeDomain==ComputationalNodeNumber) THEN
    !CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField,CMISSFieldUVariableType,1,1,LastNodeNumber,1, &
    !  & CMISSBoundaryConditionFixed,1.0_CMISSDP,Err)
  ENDIF
  !Finish the creation of the elliptic equations set boundary conditions
  CALL CMISSSolverEquationsBoundaryConditionsCreateFinish(SolverEquations,Err)

  !Set boundary conditions for the elliptic equations
  CALL CMISSBoundaryConditionsTypeInitialise(BoundaryConditions,Err)
  CALL CMISSSolverTypeInitialise(Solver,Err)
  CALL CMISSProblemSolverGet(Problem,CMISSControlLoopNode,3,Solver,Err)
  CALL CMISSSolverEquationsTypeInitialise(SolverEquations,Err)
  CALL CMISSSolverSolverEquationsGet(Solver,SolverEquations,Err)
  CALL CMISSSolverEquationsBoundaryConditionsCreateStart(SolverEquations,BoundaryConditions,Err)
  !Set the first node to 0.0 and the last node to 1.0
  FirstNodeNumber=1
  IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
    LastNodeNumber=(NUMBER_GLOBAL_X_ELEMENTS+1)*(NUMBER_GLOBAL_Y_ELEMENTS+1)
  ELSE
    LastNodeNumber=(NUMBER_GLOBAL_X_ELEMENTS+1)*(NUMBER_GLOBAL_Y_ELEMENTS+1)*(NUMBER_GLOBAL_Z_ELEMENTS+1)
  ENDIF
  CALL CMISSDecompositionNodeDomainGet(Decomposition,FirstNodeNumber,1,FirstNodeDomain,Err)
  CALL CMISSDecompositionNodeDomainGet(Decomposition,LastNodeNumber,1,LastNodeDomain,Err)
  IF(FirstNodeDomain==ComputationalNodeNumber) THEN
    !CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField,CMISSFieldVVariableType,1,1,FirstNodeNumber,1, &
    !  & CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  ENDIF
  IF(LastNodeDomain==ComputationalNodeNumber) THEN
    !CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField,CMISSFieldVVariableType,1,1,LastNodeNumber,1, &
    !  & CMISSBoundaryConditionFixed,1.0_CMISSDP,Err)
  ENDIF
  !Finish the creation of the elliptic equations set boundary conditions
  CALL CMISSSolverEquationsBoundaryConditionsCreateFinish(SolverEquations,Err)

  !Solve the problem for the first STIM_STOP
  CALL CMISSProblemSolve(Problem,Err)

  !Now turn the stimulus off
  !Set the Stimulus at node 1
  DO node_idx=1,NUMBER_OF_ELEMENTS/2
    CALL CMISSDecompositionNodeDomainGet(Decomposition,node_idx,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL CMISSFieldParameterSetUpdateNode(CellMLParametersField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,node_idx, &
        & stimcomponent,0.0_CMISSDP,Err)
    ENDIF
  ENDDO !node_idx

  !Set the time loop from STIM_STOP to TIME_STOP
  CALL CMISSControlLoopTimesSet(ControlLoop,STIM_STOP,TIME_STOP,PDE_TIME_STEP,Err)
  
  !Solve the problem for the next 900 ms
  !CALL CMISSProblemSolve(Problem,Err)
  
  EXPORT_FIELD=.FALSE.
  IF(EXPORT_FIELD) THEN
    CALL CMISSFieldsTypeInitialise(Fields,Err)
    CALL CMISSFieldsTypeCreate(Region,Fields,Err)
    CALL CMISSFieldIONodesExport(Fields,"BidomainExample","FORTRAN",Err)
    CALL CMISSFieldIOElementsExport(Fields,"BidomainExample","FORTRAN",Err)
    CALL CMISSFieldsTypeFinalise(Fields,Err)
  ENDIF
  
  !Finialise CMISS
  CALL CMISSFinalise(Err)

  WRITE(*,'(A)') "Program successfully completed."
  
  STOP
  
CONTAINS

  SUBROUTINE HANDLE_ERROR(ERROR_STRING)

    CHARACTER(LEN=*), INTENT(IN) :: ERROR_STRING

    WRITE(*,'(">>ERROR: ",A)') ERROR_STRING(1:LEN_TRIM(ERROR_STRING))
    STOP

  END SUBROUTINE HANDLE_ERROR
    
END PROGRAM BIDOMAINEXAMPLE
