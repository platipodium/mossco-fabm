module hzg_model_library

   use fabm_types, only: type_base_model_factory,type_base_model

   use hzg_omexdia_p
   use hzg_omexdia_cnp
   use hzg_omexdia_p_mpb
   use hzg_mpb
   use hzg_mpb_cnp
   use hzg_jelly
   use hzg_n2pzdq
   use hzg_medmac
   use hzg_maecs
   use hzg_benthic_pool
   use hzg_Ndepoden
   use fabm_hzg_dependencies
   use hzg_kristineb
   ! Add new HZG models here

   implicit none

   private

   type,extends(type_base_model_factory) :: type_factory
   contains
      procedure :: create
   end type

   type (type_factory),save,target,public :: hzg_model_factory

contains

   subroutine create(self,name,model)

     class (type_factory),intent(in) :: self
     character(*),        intent(in) :: name
     class (type_base_model),pointer :: model

     select case (name)
         case ('omexdia_p'); allocate(type_hzg_omexdia_p::model)
         case ('omexdia_cnp'); allocate(type_hzg_omexdia_cnp::model)
         case ('omexdia_p_mpb'); allocate(type_hzg_omexdia_p_mpb::model)
         case ('mpb'); allocate(type_hzg_mpb::model)
         case ('mpb_cnp'); allocate(type_hzg_mpb_cnp::model)
         case ('jelly'); allocate(type_hzg_jelly::model)
         case ('n2pzdq'); allocate(type_hzg_n2pzdq::model)
         case ('maecs'); allocate(type_hzg_maecs::model)
         case ('medmac'); allocate(type_hzg_medmac::model)
         case ('kristineb'); allocate(type_hzg_kristineb::model)
         case ('Ndepoden'); allocate(type_hzg_Ndepoden::model)
         case ('benthic_pool'); allocate(type_hzg_benthic_pool::model)
         case ('dependencies'); allocate(type_hzg_dependencies::model)
         ! Add case statements for new models here
     end select

   end subroutine create

end module
