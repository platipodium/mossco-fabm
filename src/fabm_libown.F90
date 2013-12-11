#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: FABM model library [Fortran 2003 models only]
!
! !INTERFACE:
   module fabm_library
!
! !USES:
   ! FABM modules
   use fabm_types

#ifdef _FABM_F2003_
   ! Specific biogeochemical models
   use fabm_hzg_omexdia_p
   use fabm_hzg_maecs
!   use fabm_iow_spm
   ! ADD_NEW_FORTRAN2003_MODEL_HERE - required
#endif

   implicit none
!
!  default: all is private.
   private

#ifdef _FABM_F2003_
   type,extends(type_abstract_model_factory),public :: type_model_factory
      contains
      procedure,nopass :: create => fabm_library_create_model
   end type
#endif
!
! !REVISION HISTORY:!
!  Original author(s): Jorn Bruggeman
!
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Function returning specific biogeochemical model given a model name.
!
! !INTERFACE:
   function fabm_library_create_model(modelname,instancename,parent,configunit) result(model)
!
! !INPUT PARAMETERS:
      character(*),intent(in)           :: modelname,instancename
      integer,     intent(in)           :: configunit
      _CLASS_ (type_model_info),target  :: parent
      _CLASS_ (type_model_info),pointer :: model
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
!EOP
!-----------------------------------------------------------------------
!BOC
      nullify(model)

#ifdef _FABM_F2003_
      select case (modelname)
!         case ('examples_mean');       allocate(type_examples_mean::model)
         case ('hzg_omexdia_p');       allocate(type_hzg_omexdia_p::model)
         case ('hzg_maecs');           allocate(type_hzg_maecs::model)
!         case ('iow_spm');             allocate(type_iow_spm::model)
         ! ADD_NEW_FORTRAN2003_MODEL_HERE - required
         case default
!            if ( modelname(1:4) .eq. 'aed_' ) &
!               model => aed_create_model(configunit,modelname,instancename,parent);
      end select

      if (.not.associated(model)) return

      ! If the model has not been initialized, do so now.
      ! This is the default - simulaneously creating and initializing the model is now deprecated.
      if (.not.associated(model%parent)) call parent%add_child(model,configunit,instancename)
#endif

   end function fabm_library_create_model
!EOC

!-----------------------------------------------------------------------

   end module fabm_library

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
