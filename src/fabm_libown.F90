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
   use fabm_types, only: type_base_model_factory, type_base_model, factory

   ! Specific biogeochemical models
   use fabm_hzg_omexdia_p
   use fabm_hzg_maecs
!   use fabm_iow_spm
   ! ADD_NEW_FORTRAN2003_MODEL_HERE - required

   implicit none
!
!  default: all is private.
   private

   public :: fabm_create_model_factory

   type,extends(type_base_model_factory) :: type_model_factory
      contains
      procedure :: create
   end type

!EOP
!-----------------------------------------------------------------------

   contains

   subroutine fabm_create_model_factory()
      if (.not.associated(factory)) then
         allocate(type_model_factory::factory)

         !call factory%add(aed_model_factory)
         ! Add new additional model factories here

      end if
   end subroutine
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Function returning specific biogeochemical model given a model name.
!
! !INTERFACE:
   subroutine create(self,name,model)
!
! !INPUT PARAMETERS:
      class (type_model_factory),intent(in) :: self
      character(*),              intent(in) :: name
      class (type_base_model),pointer       :: model
!
!EOP
!-----------------------------------------------------------------------
!BOC
      nullify(model)


      select case (name)
!         case ('examples_mean');       allocate(type_examples_mean::model)
         case ('hzg_omexdia_p');       allocate(type_hzg_omexdia_p::model)
         case ('hzg_maecs');           allocate(type_hzg_maecs::model)
!         case ('iow_spm');             allocate(type_iow_spm::model)
         ! ADD_NEW_FORTRAN2003_MODEL_HERE - required
         case default
            call self%type_base_model_factory%create(name,model)
      end select

   end subroutine
!EOC

!-----------------------------------------------------------------------

   end module fabm_library

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
