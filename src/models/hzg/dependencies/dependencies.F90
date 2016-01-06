#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: fabm_hzg_dependencies --- dependency test model
!
! !INTERFACE:
   module fabm_hzg_dependencies
!
! !DESCRIPTION:
!
! This model writes diagnostics of its own dependencies
!
! !USES:
   use fabm_types

   implicit none

!  default: all is private.
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public type_hzg_dependencies
!
! !PRIVATE DATA MEMBERS:
!
! !REVISION HISTORY:!
!  Original author(s): Richard Hofmeister
!
!
! !PUBLIC DERIVED TYPES:
   type,extends(type_base_model) :: type_hzg_dependencies
!     Variable identifiers
      type (type_dependency_id)            :: id_temp,id_par
      type (type_diagnostic_variable_id)   :: id_dtemp,id_dpar

      contains

!     Model procedures
      procedure :: initialize
      procedure :: do

   end type type_hzg_dependencies
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the model
!
! !INTERFACE:
   subroutine initialize(self,configunit)
!
! !DESCRIPTION:
!
! !INPUT PARAMETERS:
   class (type_hzg_dependencies),intent(inout),target  :: self
   integer,                      intent(in)            :: configunit
!
! !REVISION HISTORY:
!  Original author(s): Richard Hofmeister
!
! !LOCAL VARIABLES:

!EOP
!-----------------------------------------------------------------------
!BOC

   ! Register diagnostic variables
   call self%register_diagnostic_variable(self%id_dtemp,'tempdiag','degC', &
         'diagnostic temperature', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_dpar,'pardiag','W/m**2', &
         'diagnostic PAR', output=output_instantaneous)

   ! Register dependencies
   call self%register_dependency(self%id_temp,standard_variables%temperature)
   call self%register_dependency(self%id_par,standard_variables%downwelling_photosynthetic_radiative_flux)

   return

   end subroutine initialize
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Right hand sides of model
!
! !INTERFACE:
   subroutine do(self,_ARGUMENTS_DO_)
!
! !DESCRIPTION:
!
! !INPUT PARAMETERS:
   class (type_hzg_dependencies),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_
!
! !REVISION HISTORY:
!  Original author(s): Richard Hofmeister
!
! !LOCAL VARIABLES:
   real(rk) :: temp, par
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   ! Retrieve current (local) dependencies
   _GET_(self%id_temp,temp)
   _GET_(self%id_par,par)

   ! Export diagnostic variables
   _SET_DIAGNOSTIC_(self%id_dtemp,temp)
   _SET_DIAGNOSTIC_(self%id_dpar ,par)

   ! Leave spatial loops (if any)
   _LOOP_END_

   end subroutine do
!EOC

   end module fabm_hzg_dependencies

