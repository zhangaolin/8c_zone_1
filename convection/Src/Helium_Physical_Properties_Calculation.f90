!> @brief Physical property parameter calculation of helium for thermal hydraulic calculation
!!
!> Physical property parameters are calculated including density, dynamic viscosity, thermal conductivity and a constant value for isobaric specific heat of 5195J/kg·K.
!> These correlations for helium are obtained form <b> H. Petersen. The Properties of Helium: Density, Specific Heats, Viscosity, and Thermal Conductivity at Pressures from 1 to 100 bar and from Room Temperature to about 1800 K.Technical report, Danish Atomic Energy Commission, 1970.</b> The report can been found in \link https://backend.orbit.dtu.dk/ws/portalfiles/portal/52768509/ris_224.pdf \endlink
!> These correlations are valid in the range 1 - P - 100 bar and 273 - Tf - 1800 K, and have been sed by KTA in their analysis of HTRs <b>KTA. Reactor Core Design of High-Temperature Gas-Cooled Reactors Part 1: Calculation of the Material Properties of Helium. Technical Report KTA 3102.1, Nuclear Safety Standards Commission, 1978.<\b>
!> @author XUDONGYU
!> @version 1.0
!> @date 
MODULE HeliumPhysicalPropertiesCalculation
    USE GlobalTHConstants !< Global constants definition for thermal hydraulic calculation

    IMPLICIT NONE

CONTAINS
    !> @brief A subroutine that calculate the helium density
    !! 
    !> @note \f$\rho_h=48.14\frac{P}{T}\left(1+0.4446\frac{P}{T^{1.2}}\right)^{-1}\f$
    !> @note This subroutine has the same functionality as the <b>DICHTE</b> 
    !> @date 
    SUBROUTINE cal_Helium_Density(density,pressure,temperature)
        REAL(TH_KDUBLE), INTENT(INOUT)  :: density     !< The density of helium
        REAL(TH_KDUBLE), INTENT(IN)     :: pressure    !< The pressure of helium
        REAL(TH_KDUBLE), INTENT(IN)     :: temperature !< The temperature of helium
        REAL(TH_KDUBLE)                 :: t_k !< The temperature of helium in Kelvin

        t_k = temperature + TH_CKELVIN
        IF (t_k .LT. TH_REAL_ONE) t_k = TH_REAL_ONE
        density = 48.14*pressure / t_k / (TH_REAL_ONE+0.4446*pressure/(t_k**1.2))
    END SUBROUTINE cal_Helium_Density
    !> @brief A subroutine that calculate the helium viscosity
    !! 
    !> @note \f$\eta_h =3.674\times10_{-7}T^{0.7}\f$
    !> @note This subroutine has the same functionality as the <b>ETHAG</b> 
    !> @date 
    SUBROUTINE cal_Helium_Viscosity(viscosity,temperature)
        REAL(TH_KDUBLE), INTENT(INOUT)  :: viscosity   !< The viscosity of helium
        REAL(TH_KDUBLE), INTENT(IN)     :: temperature !< The temperature of helium
        REAL(TH_KDUBLE)                 :: t_k !< The temperature of helium in Kelvin

        t_k = temperature + TH_CKELVIN
        IF (t_k .LT. TH_REAL_ONE) t_k = TH_REAL_ONE
        viscosity = 3.674D-07*(t_k)**0.7
    END SUBROUTINE cal_Helium_Viscosity
    !> @brief A subroutine that calculate the helium thermal conductivity
    !! 
    !> @note \f$\lambda_h=2.682\times10^{-3}(1+1.123\times10^{-3})T^{0.71(1-2.0\times10^{-4}P)}\f$
    !> @note This subroutine has the same functionality as the <b>XLHE</b> 
    !> @date 
    SUBROUTINE cal_Helium_Conductivity(conductivity,pressure,temperature)
        REAL(TH_KDUBLE), INTENT(INOUT)  :: conductivity   !< The thermal conductivity of helium
        REAL(TH_KDUBLE), INTENT(IN)     :: pressure    !< The pressure of helium
        REAL(TH_KDUBLE), INTENT(IN)     :: temperature !< The temperature of helium
        REAL(TH_KDUBLE)                 :: t_k !< The temperature of helium in Kelvin

        t_k = temperature + TH_CKELVIN
        IF (t_k .LT. TH_REAL_ONE) t_k = TH_REAL_ONE
        conductivity = 0.71*(1.0-2.0D-04*pressure)
        t_k = t_k**conductivity
        conductivity = 2.682D-03*(1+1.123D-03*pressure)*t_k
    END SUBROUTINE cal_Helium_Conductivity

END MODULE HeliumPhysicalPropertiesCalculation