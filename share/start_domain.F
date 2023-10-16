!WRF:MEDIATION_LAYER:ADT_BARRIER
!

SUBROUTINE start_domain ( grid , allowed_to_read )

   USE module_domain
   USE module_configure

   IMPLICIT NONE

   !  Input Arguments.
   TYPE (domain)          :: grid
   LOGICAL, INTENT(IN)    :: allowed_to_read
   !  Local data.
   INTEGER :: idum1, idum2
   LOGICAL :: already_went_there

   INTERFACE
#if ((EM_CORE == 1) && (DA_CORE != 1))
#    include "start_domain_em.int"
#endif
#if ((EM_CORE == 1) && (WRFPLUS == 1))
#    include "../wrftladj/start_domain_em_tl.int"
#    include "../wrftladj/start_domain_em_ad.int"
#endif
#if (DA_CORE != 1)
      SUBROUTINE calc_ts_locations( grid )
         USE module_domain
         TYPE (domain) :: grid
      END SUBROUTINE calc_ts_locations

      SUBROUTINE calc_track_locations( grid )
         USE module_domain
         TYPE (domain) :: grid
      END SUBROUTINE calc_track_locations
#endif
   END INTERFACE

   CALL set_scalar_indices_from_config ( head_grid%id , idum1, idum2 )
   
   !Test to not reinitialize tslist
   already_went_there = .FALSE.
   IF ( grid%itimestep > 0 ) THEN
      already_went_there = .TRUE.
   END IF

#if ((EM_CORE == 1) && (DA_CORE != 1))

   WRITE(*,*) "CALC_TS_LOCATION BEGINNING OF START_DOMAIN, grid%id = ", grid%id, &
   & ", grid%have_calculated_tslocs", grid%have_calculated_tslocs, &
   & ", grid%itimestep", grid%itimestep, &
   & ", already_went_there", already_went_there

#if ( WRFPLUS == 1 )
   IF ( model_config_rec%dyn_opt .EQ. dyn_em_tl ) THEN
      CALL g_start_domain_em( grid, allowed_to_read  &
!
# include "actual_new_args.inc"
!
                         )
   ELSE IF ( model_config_rec%dyn_opt .EQ. dyn_em_ad ) THEN
      CALL a_start_domain_em( grid, allowed_to_read  &
!
# include "actual_new_args.inc"
!
                         )
   ELSE
#endif

      CALL start_domain_em( grid, allowed_to_read  &
!
# include "actual_new_args.inc"
!
                         )

#if ( WRFPLUS == 1 )
   ENDIF
#endif

#endif

#if (DA_CORE != 1)
   IF ( .not. already_went_there ) THEN
      WRITE(*,*) "CALC_TS_LOCATION, grid%id = ", grid%id, &
      & ", grid%have_calculated_tslocs", grid%have_calculated_tslocs, &
      & ", grid%itimestep", grid%itimestep
      CALL calc_ts_locations( grid )
      CALL calc_track_locations( grid )
   END IF
#endif

END SUBROUTINE start_domain

