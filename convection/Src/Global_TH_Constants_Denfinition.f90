!> @brief Global constants definition for thermal hydraulic calculation
!!
!> This module is used to define the constants used in the thermal calculation of the pebble-bed HTGR.
!> To avoid conflicts with other global constant definitions, constants definitions for this module should include the TH identifier as much as possible.
!!
!> @author XUDONGYU
!> @version 1.0
!> @date 
MODULE GlobalTHConstants

    IMPLICIT NONE
    PUBLIC
    !< Kind number for integer, real and character
    INTEGER, PARAMETER  :: TH_KINT             = SELECTED_INT_KIND(r=7) !< A high precision integer kind
    INTEGER, PARAMETER  :: TH_KERAL            = 4                      !< A normal percision real kind
    INTEGER, PARAMETER  :: TH_KDUBLE           = 4                      !< A high percision real kind
    INTEGER, PARAMETER  :: TH_MAX_LINE_LEN     = 800                    !< The max length of a line
    INTEGER, PARAMETER  :: TH_MAX_WORD_LEN     = 800                    !< The max length of a word
    INTEGER, PARAMETER  :: TH_MAX_WORDS        = 800                    !< The max number of word per line
    !< Comparision criteria
    !< Mathmatic constants
    INTEGER, PARAMETER  :: TH_INT_ZERO = 0 !< THe zero of a int variable
    REAL(TH_KDUBLE), PARAMETER  :: HALF = 0.5D0 !< The half of a real variable
    REAL(TH_KDUBLE), PARAMETER  :: TH_REAL_ZERO = 0.0D0 !< The zero of a real variable
    REAL(TH_KDUBLE), PARAMETER  :: TH_REAL_ONE = 1.0D0 !< The one of a real variable
    REAL(TH_KDUBLE), PARAMETER  :: TH_REAL_HUNDREAD = 1.0D2 !< The one hundred of a real variable
    REAL(TH_KDUBLE), PARAMETER  :: TH_REAL_MILLION = 1.0D6 !< The one million of a real variable
    REAL(TH_KDUBLE), PARAMETER  :: TH_REAL_INFINITY = 1.0D20 !< The positive infinity of a real variable
    REAL(TH_KDUBLE), PARAMETER  :: TH_ANGLE_CIRCUMFERENCE = 360.0D0 !< The value of angle of circumference
    REAL(TH_KDUBLE), PARAMETER  :: TH_PI = 3.14159D0 !< The value of PI
    !< Physics constants
    REAL(TH_KDUBLE), PARAMETER  :: TH_CKELVIN = 273.D0 !< The temperature transfer between celsius temperature and thermodynamic temperature
    REAL(TH_KDUBLE), PARAMETER  :: TH_GRAVITY = 9.81D0 !< the vaule of gravity

    
    CHARACTER(LEN=TH_MAX_WORD_LEN) :: INPUT_FILE = "INPUT.DAT" !< The input file name
    CHARACTER(LEN=TH_MAX_WORD_LEN) :: OUTPUT_FILE = "OUTPUT.DAT" !< The output file name
    INTEGER :: INUNIT = 101 !< The input file unit number
    INTEGER :: OUTUNIT = 201 !< The output file unit number

    CHARACTER(LEN=TH_MAX_WORD_LEN), PARAMETER :: CHAR_SENTINEL = '----' !< The sentinel character

    
END MODULE GlobalTHConstants
