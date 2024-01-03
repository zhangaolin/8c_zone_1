 !> @brief A module for the composition of thermal hydraulic calculations
!! 
!> @author XUDONGYU
!> @version 1.0
!> @date 
MODULE THCompositionHeader
    USE GlobalTHConstants !< Global constants definition for thermal hydraulic calculation
    IMPLICIT NONE
    !> @class ChannelData
    !> @brief The type of information of Composition to thermal hydraulic calculations
    !> @date 
    TYPE CompositionData
        INTEGER  :: composition_id !< The id of this composition
        INTEGER  :: is_heat_exchange !< Whether there is heat exchange between solid material zone with the coolant
        INTEGER  :: method_cal_heat_capacity !< The method of calculate heat capacity
        INTEGER  :: method_cal_thermal_conductivity !< The method of calculate thermal conductivity
        INTEGER  :: direction_radiation !< The direction of radiation; 0: Radiation in radial direction;1: Radiation in axial direction.
        REAL(TH_KDUBLE)  :: ntvar !< In case of fluid zone (IFTV = 1) provide NTVAR time dependent temperatures on card TX12 
        REAL(TH_KDUBLE)  :: vol_fraction_soild_material !< Volumetric fraction of solid material in this composition. RHO is used for calculation of the heat capacity
        REAL(TH_KDUBLE)  :: heat_capacity !< Heat capacity of the solid material
        REAL(TH_KDUBLE)  :: thermal_conductivity !< Thermal conductivity in solid material zones
        REAL(TH_KDUBLE)  :: lam0
        REAL(TH_KDUBLE)  :: coefficient_heat_radiation !< Coefficient of emission for heat radiation
        REAL(TH_KDUBLE)  :: coefficient_heat_radiation_wall !< Coefficient of emission for heat radiation in outer/lower wall
        REAL(TH_KDUBLE)  :: horizontal_gap !< Horizontal gap
        REAL(TH_KDUBLE)  :: composition_input_temperature !< Fixed temperature of this composition superior to the start-up temperatures of the cards
        REAL(TH_KDUBLE)  :: composition_power_density !< Power density of this composition
        REAL(TH_KDUBLE)  :: coefficient_heat_transfer !< Heat transfer coefficient in fluid zones
        INTEGER  :: is_gas_streaming !< Whether there is gas streaming through the area
        REAL(TH_KDUBLE)  :: diameter_fuel_element !< Diameter of the spherical fuel element.
        INTEGER  :: number_radial_mesh !< Number of radial mesh intervals in the sphere.

    END TYPE CompositionData
    !> @class ChannelData
    !> @brief The type of information of Channel to thermal hydraulic calculations
    !> @date 
    TYPE ChannelData
        INTEGER  :: is_convective_heat_source !< Whether there is convective heat source 
        INTEGER  :: channel_type = 0 !< The type of composition
        REAL(TH_KDUBLE)  :: pressure_begin = -1.0 !< The pressure at beginning of iterations. (bar)
        REAL(TH_KDUBLE)  :: additional_pressure_drop !< The Additional pressure drop relative to computed pressure drop over the length of the channel (only if IFBR = 2)
        REAL(TH_KDUBLE)  :: heat_transition_coefficient !< The coefficient of heat transition
        REAL(TH_KDUBLE)  :: hydraulic_diameter !< The hydraulic diameter 
        REAL(TH_KDUBLE)  :: source_mass_flow !< Source of mass flow
        REAL(TH_KDUBLE)  :: source_temperature_in !< Temperature of inlet gas
        INTEGER  :: is_forced_cooling !< Whether there is forced cooling

    END TYPE ChannelData
    
END MODULE THCompositionHeader