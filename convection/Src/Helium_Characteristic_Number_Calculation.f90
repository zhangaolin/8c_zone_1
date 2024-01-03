!> @brief Characteristic number calculation of helium for thermal hydraulic calculation
!!
!> Characteristic number calculation of helium are calculated including Relo number
!!
!> This module contains not only the calculation of characteristic numbers, but also the calculation of scalar velocities
!> @author XUDONGYU
!> @version 1.0
!> @date 
MODULE HeliumCharacteristicNumberCalculation
    USE GlobalTHConstants !< Global constants definition for thermal hydraulic calculation
    USE THMathematicalCalculation !< A module that defines common mathematical functions
    USE HeliumPhysicalPropertiesCalculation !< Physical property parameter calculation of helium for thermal hydraulic calculation
    IMPLICIT NONE

CONTAINS
    !> @brief A subroutine that calculate the helium Relo number
    !! 
    !> @note This subroutine has the same functionality as the <b>REYN</b> 
    !> @date 
    SUBROUTINE cal_Helium_Re(relo,temperature,density,channel_type,porosity,diameter_fuel_element,tube_diameter,axial_area,radial_area,circumferential_area,axial_mass,radial_mass,circumferential_mass)
        REAL(TH_KDUBLE), INTENT(INOUT)  :: relo     !< The Relo number of this element
        REAL(TH_KDUBLE), INTENT(IN)  :: temperature !< The temperature of helium in this element
        REAL(TH_KDUBLE), INTENT(IN)  :: density !< The density of helium in this element
        INTEGER, INTENT(IN)  :: channel_type !< The type of  this element
        REAL(TH_KDUBLE), INTENT(IN)  :: porosity !< The porosity of  this element
        REAL(TH_KDUBLE), INTENT(IN)  :: diameter_fuel_element !< !< Diameter of the spherical fuel element
        REAL(TH_KDUBLE), INTENT(IN)  :: tube_diameter !< The tube hydraulic diameter of this element
        REAL(TH_KDUBLE), INTENT(IN)  :: axial_area !< The axial area of this element
        REAL(TH_KDUBLE), INTENT(IN)  :: radial_area !< The radial area of this element
        REAL(TH_KDUBLE), INTENT(IN)  :: circumferential_area !< The circumferential area of this element
        REAL(TH_KDUBLE), INTENT(IN)  :: axial_mass !< The axial mass of this element
        REAL(TH_KDUBLE), INTENT(IN)  :: radial_mass !< The radial mass of this element
        REAL(TH_KDUBLE), INTENT(IN)  :: circumferential_mass !< The circumferential mass of this element
        REAL(TH_KDUBLE)  :: viscosity,velocity
        ! 这种分方向的变量，要不要变成一个数组传递进来
        ! I,N,KOM均在函数外部完成计算；TFL-temperature;RHO-Helium_density;IFB-type;EPSIL-porosity;DHYD-tube_diameter;FZQ-axial_area;FRQ-radial_area;MZ-axial_mass;MR-radial_mass
        relo = TH_REAL_ZERO
        viscosity = TH_REAL_ZERO
        velocity = TH_REAL_ZERO
        IF(channel_type .EQ. 0) THEN
            
        ELSE IF (channel_type .EQ. 1) THEN
            CALL cal_Helium_Viscosity(viscosity,temperature)
            CALL cal_Scalar_Velocity(velocity,density,axial_area,radial_area,circumferential_area,axial_mass,radial_mass,circumferential_mass)

            relo = 0.6666*diameter_fuel_element*density*velocity/viscosity/(1-porosity)
        ELSE IF (channel_type .EQ. 2) THEN
            CALL cal_Helium_Viscosity(viscosity,temperature)
            CALL cal_Scalar_Velocity(velocity,density,axial_area,radial_area,circumferential_area,axial_mass,radial_mass,circumferential_mass)

            relo = density*velocity*tube_diameter/porosity/viscosity
        ELSE IF (channel_type .EQ. 5) THEN
            CALL cal_Helium_Viscosity(viscosity,temperature)
            CALL cal_Scalar_Velocity(velocity,density,axial_area,radial_area,circumferential_area,axial_mass,radial_mass,circumferential_mass)

            relo = density*velocity*tube_diameter/porosity/viscosity
        END IF
    END SUBROUTINE cal_Helium_Re
    !> @brief A subroutine that calculate the helium axial Relo number
    !! 
    !> @note This subroutine has the same functionality as the <b>XIRE</b> 
    !> @date 
    SUBROUTINE cal_Helium_AxialRe(relo)
        REAL(TH_KDUBLE), INTENT(INOUT)  :: relo     !< The relo number of this element
        IF (relo .GT. 3.0D4) THEN
            relo = 0.025633189
        ELSE IF (relo .LT. 0.00001) THEN
            relo = 6.4D6
        ELSE
            relo = 64/relo+0.0235*(1-EXP(-4D-4*relo))
        END IF
    END SUBROUTINE cal_Helium_AxialRe
    !> @brief A subroutine that calculate the helium Nussel number in pebble bed
    !! 
    !> @note \f$Nu = 1.27\frac{Pr^{1/3}Re^{0.36}}{\varepsilon^{1.18}}+0.33\frac{Pr^{0.5}Re^{0.86}}{\varepsilon^{1.07}}\f$
    !> @date   
    SUBROUTINE cal_Helium_Nu_Pebble_Bed(nussel,prandtl,relo,porosity)
        REAL(TH_KDUBLE), INTENT(INOUT)  :: nussel     !< The Nussel number of this element
        REAL(TH_KDUBLE), INTENT(IN)  :: prandtl !< The Prandtl number of this element
        REAL(TH_KDUBLE), INTENT(IN)  :: relo !< The Relo number of this element
        REAL(TH_KDUBLE), INTENT(IN)  :: porosity !< The porosity of this element
        IF (relo .GE. 2.0) THEN
            nussel = 1.27*(prandtl**0.33 * relo**0.36) / porosity**1.18
            nussel = nussel + 0.033*(prandtl**0.5 * relo**0.86) / porosity**1.07
        ELSE
            nussel = 4.03
        END IF
    END SUBROUTINE cal_Helium_Nu_Pebble_Bed
    !> @brief A subroutine that calculate the helium thermal conductivity
    !! 
    !> @note \f$Nu = 1.27\frac{Pr^{1/3}Re^{0.36}}{\varepsilon^{1.18}}+0.33\frac{Pr^{0.5}Re^{0.86}}{\varepsilon^{1.07}}\f$
    !> @date   
    SUBROUTINE cal_Helium_Nu_Vertical_Pipes(nussel,prandtl,relo,diameter)
        REAL(TH_KDUBLE), INTENT(INOUT)  :: nussel     !< The Nussel number of this element
        REAL(TH_KDUBLE), INTENT(IN)  :: prandtl !< The Prandtl number of this element
        REAL(TH_KDUBLE), INTENT(IN)  :: relo !< The Relo number of this element
        REAL(TH_KDUBLE), INTENT(IN)  :: diameter !< The diameter of this pipe
        REAL(TH_KDUBLE) :: relo_max = 10000.D0,relo_min = 2300.D0,d_limit = 0.1,nussel_temp 
        IF(relo .GT. relo_max) THEN
            nussel = 0.0214*(relo**0.8 - 100.0)*prandtl**0.4*(1.0+diameter**0.6667)
        ELSE IF (relo .LT. relo_min) THEN
            IF (diameter .LE. d_limit) THEN
                nussel = 3.66**3 + 1.61**3*relo*prandtl*diameter
            ELSE
                nussel = 3.66**3 + 1.61**3*relo*prandtl*diameter
                nussel_temp = 0.664*prandtl**0.3333*(relo*diameter)**0.5
                IF (nussel_temp .GT. nussel) nussel = nussel_temp
            END IF
        ELSE
            IF (diameter .LE. d_limit) THEN
                nussel = 3.66**3 + 1.61**3*relo*prandtl*diameter
                nussel_temp = 0.0214*(relo**0.8 - 100.0)*prandtl**0.4*(1.0+diameter**0.6667)
                CALL linear_Interpolation(nussel,relo,relo_min,relo_max,nussel_temp,nussel)
            ELSE
                nussel = 3.66**3 + 1.61**3*relo*prandtl*diameter
                nussel_temp = 0.664*prandtl**0.3333*(relo*diameter)**0.5
                IF (nussel_temp .GT. nussel) nussel = nussel_temp
                CALL linear_Interpolation(nussel,relo,relo_min,relo_max,nussel_temp,nussel)
            END IF
        END IF
            

    END SUBROUTINE cal_Helium_Nu_Vertical_Pipes
    !> @brief A subroutine that calculate the Scalar Velocity
    !! 
    !> @note This subroutine has the same functionality as the <b>VAU</b> 
    !> @date 
    SUBROUTINE cal_Scalar_Velocity(Velocity,density,axial_area,radial_area,circumferential_area,axial_mass,radial_mass,circumferential_mass)
        REAL(TH_KDUBLE), INTENT(INOUT)  :: Velocity     !< The Scalar Velocity of this element
        REAL(TH_KDUBLE), INTENT(IN)  :: density !< The density of helium in this element
        REAL(TH_KDUBLE), INTENT(IN)  :: axial_area !< The axial area of this element
        REAL(TH_KDUBLE), INTENT(IN)  :: radial_area !< The radial area of this element
        REAL(TH_KDUBLE), INTENT(IN)  :: circumferential_area !< The circumferential area of this element
        REAL(TH_KDUBLE), INTENT(IN)  :: axial_mass !< The axial mass of this element
        REAL(TH_KDUBLE), INTENT(IN)  :: radial_mass !< The radial mass of this element
        REAL(TH_KDUBLE), INTENT(IN)  :: circumferential_mass !< The circumferential mass of this element
        REAL(TH_KDUBLE)  :: temp_axial,temp_radial,temp_circumferential !< Temporary variables
        ! 这种分方向的变量，要不要变成一个数组传递进来
        ! I,N均在函数外部完成计算；RHO-Helium_density;FZQ-axial_area;FRQ-radial_area;MZ-axial_mass;MR-radial_mass
        temp_circumferential = TH_REAL_ZERO
        temp_radial = TH_REAL_ZERO
        temp_axial = TH_REAL_ZERO

        IF (ABS(circumferential_mass)*circumferential_area .GT. TH_REAL_ZERO) THEN !< a*b > 0 -> a>0&b>0
            temp_circumferential = (circumferential_mass/circumferential_area) * (circumferential_mass/circumferential_area)
        END IF
        IF (ABS(radial_mass)*radial_area .GT. TH_REAL_ZERO) THEN !< a*b > 0 -> a>0&b>0
            temp_radial = (radial_mass/radial_area) * (radial_mass/radial_area)
        END IF
        IF (ABS(axial_mass)*axial_area .GT. TH_REAL_ZERO) THEN
            temp_axial = (axial_mass/axial_area) * (axial_mass/axial_area)
        END IF

        Velocity = SQRT(temp_circumferential+temp_radial+temp_axial)/density

    END SUBROUTINE cal_Scalar_Velocity

END MODULE HeliumCharacteristicNumberCalculation