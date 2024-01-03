!> @brief Global variables definition for thermal hydraulic calculation
!!
!> This module is used to define the variables used in the thermal calculation of the pebble-bed HTGR.
!> To avoid conflicts with other global variables definitions, variables definitions for this module should include the TH identifier as much as possible.
!> @author XUDONGYU
!> @version 1.0
!> @date 
MODULE GlobalTHVariables
    USE GlobalTHConstants !< Global constants definition for thermal hydraulic calculation
    USE THMeshHeader      !< A module for the mesh of thermal hydraulic calculations
    USE THCompositionHeader !< A module for the composition of thermal hydraulic calculations
    IMPLICIT NONE
    PUBLIC
    CHARACTER(LEN=TH_MAX_WORD_LEN) :: CASENAME
    !< Limit of convergence
    REAL(TH_KDUBLE)  :: Th_epsi_gas_tem           = 1.0D-5 !< Relative criterion of convergence for gas temperature
    REAL(TH_KDUBLE)  :: Th_epsi_gas_mass_flow     = 1.0D-2 !< Criterion of convergence for mass flow
    REAL(TH_KDUBLE)  :: Th_epsi_avg_gas_tem       = 2.0D-2 !< Relative criterion of convergence of the avg. gas temperature in the outer iterations between gas temperature and mass flow
    REAL(TH_KDUBLE)  :: The_extropolation_factor  = 0.5D0 !< Extrapolation factor for the mass flow
    REAL(TH_KDUBLE)  :: Th_max_sor = 1.7D0 !< Maximum relaxation factor
    REAL(TH_KDUBLE)  :: Th_min_sor = 0.6D0 !< Minimum relaxation factor
     REAL(TH_KDUBLE)  :: Th_sor = 0.6D0 !< Relaxation factor
    !< Maximum iterations
    INTEGER          :: Th_max_sor_change_time = 100 !< Maximum number of changes of the relaxation factor
    INTEGER          :: Th_max_iter_gas_tem       = 100 !< Maximum number of iterations for gas temperature
    INTEGER          :: Th_max_iter_gas_mass_flow = 500 !< Maximum number of iterations for mass flow
    INTEGER          :: Th_max_avg_gas_tem        = 5 !< Maximum number of outer iterations between gas temperature and mass flow

    TYPE(THMesh)                :: Th_mesh !< the mesh of thermal hydraulic calculations

    REAL(TH_KDUBLE),ALLOCATABLE :: T_fluid(:,:,:) !< The temperature of fluid
    REAL(TH_KDUBLE),ALLOCATABLE :: T_solid(:,:,:) !< The temperature of solid
    REAL(TH_KDUBLE),ALLOCATABLE :: P_fluid(:,:,:) !< The pressure of solid
    REAL(TH_KDUBLE),ALLOCATABLE :: Helium_density(:,:,:) !< The density of fluid
    REAL(TH_KDUBLE),ALLOCATABLE :: Axial_mass_flow(:,:,:) !< The axial interphase friction factor 
    REAL(TH_KDUBLE),ALLOCATABLE :: Radial_mass_flow(:,:,:) !< The radial interphase friction factor 
    REAL(TH_KDUBLE),ALLOCATABLE :: Circumferential_mass_flow(:,:,:) !< The circumferential interphase friction factor 
    ! 这种一个网格内分方向的变量能不能存在一个数组里面，然后增加一个维度用于区分方向
    REAL(TH_KDUBLE),ALLOCATABLE :: Hydraulic_diameter(:) !< The hydraulic diameter
    REAL(TH_KDUBLE),ALLOCATABLE :: Additional_pressure_drop(:) !< The additional pressure drop relative to computed pressure drop over the length of the channel 
    REAL(TH_KDUBLE),ALLOCATABLE :: Source_mass_flow(:) !< Source of mass flow OR Mass flow according to conservation law 
    REAL(TH_KDUBLE),ALLOCATABLE :: Source_temperature_in(:) !< Temperature of inlet gas
    REAL(TH_KDUBLE),ALLOCATABLE :: Xgeo(:) !< xgeo
    REAL(TH_KDUBLE),ALLOCATABLE :: Qtv(:) !< Qtv
    REAL(TH_KDUBLE),ALLOCATABLE :: Ltv(:) !< Ltv
    REAL(TH_KDUBLE),ALLOCATABLE :: Feld(:) !< FELD
    REAL(TH_KDUBLE),ALLOCATABLE :: Coefficient_heat_transition(:) !<  Internal calculation of the coefficient of heat transition
    TYPE(CompositionData),ALLOCATABLE  :: Compositions(:) !< The information of composition 
    TYPE(ChannelData),ALLOCATABLE  :: Channels(:) !< The information of channel 

    REAL(TH_KDUBLE) :: Gas_heat_capacity !< Specific heat capacity of the gas. (J/kg/K) Default value = 5195. (He)
    REAL(TH_KDUBLE) :: Gas_prandtl_number !< Prandtl-constant of the gas. Default value = 0.66
    REAL(TH_KDUBLE) :: Gas_pressure !< Pressure of the gas. (bar)
    REAL(TH_KDUBLE) :: Prandtl_number !< Prandtl-constant of the gas

    INTEGER,ALLOCATABLE :: material_channel2(:) !< The material id of channel type 2
    REAL(TH_KDUBLE),ALLOCATABLE :: volume_channel2(:) !< The volume of channel type 2
    INTEGER,ALLOCATABLE :: upper_channel2(:,:,:) !< The upper boundary of channel type 2
    INTEGER,ALLOCATABLE :: lower_channel2(:,:,:) !< The lower boundary of channel type 2
    
END MODULE GlobalTHVariables
