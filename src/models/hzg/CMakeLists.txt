add_library(fabm_models_hzg OBJECT
            hzg_model_library.F90
            omexdia_p/omexdia_p.F90
            omexdia_cnp/omexdia_cnp.F90
            omexdia_p_mpb/omexdia_p_mpb.F90
            mpb/mpb.F90
            mpb_cnp/mpb_cnp.F90
            jelly/jelly.F90
            n2pzdq/n2pzdq.F90
            Ndepoden/Ndepoden.F90
            benthic_pool/benthic_pool.F90
            medmac/medmac.F90
            maecs/maecs.F90
            maecs/maecs_types.F90
            #maecs/maecs_do.F90
            maecs/maecs_functions.F90
            maecs/maecs_grazing.F90
            maecs/maecs_primprod.F90
            dependencies/dependencies.F90
            kristineb/kristineb_gen.F90
           )

add_dependencies(fabm_models_hzg fabm_base)

#-----------------------------------------------------------------------
# Copyright (C) 2014 - Karsten Bolding and Jorn Bruggeman (BB)         !
#-----------------------------------------------------------------------
