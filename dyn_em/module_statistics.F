!WRF:MODEL_LAYER:DYNAMICS
!
MODULE module_statistics

  USE module_model_constants, ONLY : R_v, R_d, T0, g
  USE module_domain, ONLY : domain
  USE module_state_description, ONLY : P_QV

CONTAINS

!-------------------------------------------------------------------------------
  !The calculation of time-average values is basic :
  ! <x>_n = sum_{i=0}^{i=n}(x_i) / n
  ! Then :
  ! <x>_0 = 0.                                      -> computed in zero_statistics
  ! <x>_{n+1} = (n*<x>_n + x_{n+1}) / (n+1)         -> computed in upd_statistics
  !There might exist better methods

  SUBROUTINE zero_statistics(grid, output_stats, its, ite, jts, jte, kts, kte)

    IMPLICIT NONE

    TYPE(domain)             :: grid
    INTEGER, INTENT(IN   )   :: output_stats
    INTEGER, INTENT(IN   )   :: its, ite, jts, jte, kts, kte !tile limit indices

    INTEGER :: i,j,k

    if (output_stats > 0) then
      DO j=jts,jte
      DO i=its,ite
        grid%mut_avg(i,j) = 0.
      ENDDO
      ENDDO
      DO j=jts,jte
      DO k=kts,kte
      DO i=its,ite
        grid%u_avg(i,k,j) = 0.
        grid%v_avg(i,k,j) = 0.
        grid%w_avg(i,k,j) = 0.
        grid%q_avg(i,k,j) = 0.
        grid%t_avg(i,k,j) = 0.
        grid%p_avg(i,k,j) = 0.
        grid%z_avg(i,k,j) = 0.
      END DO
      END DO
      END DO
      if (output_stats > 1) then
        DO j=jts,jte
        DO k=kts,kte
        DO i=its,ite
          grid%u2_avg(i,k,j) = 0.
          grid%uv_avg(i,k,j) = 0.
          grid%uw_avg(i,k,j) = 0.
          grid%uq_avg(i,k,j) = 0.
          grid%ut_avg(i,k,j) = 0.
          grid%v2_avg(i,k,j) = 0.
          grid%vw_avg(i,k,j) = 0.
          grid%vq_avg(i,k,j) = 0.
          grid%vt_avg(i,k,j) = 0.
          grid%w2_avg(i,k,j) = 0.
          grid%wq_avg(i,k,j) = 0.
          grid%wt_avg(i,k,j) = 0.
          grid%q2_avg(i,k,j) = 0.
          grid%qt_avg(i,k,j) = 0.
          grid%t2_avg(i,k,j) = 0.
        END DO
        END DO
        END DO
        if (output_stats > 2) then
          DO j=jts,jte
          DO k=kts,kte
          DO i=its,ite
            grid%u3_avg(i,k,j) = 0.
            grid%v3_avg(i,k,j) = 0.
            grid%w3_avg(i,k,j) = 0.
            grid%q3_avg(i,k,j) = 0.
            grid%t3_avg(i,k,j) = 0.
          END DO
          END DO
          END DO
        endif
      endif
    endif
    
    return
  END SUBROUTINE zero_statistics

  SUBROUTINE upd_statistics(grid, output_stats, output_dry_theta_fluxes, its, ite, jts, jte, kts, kte)
    ! Description
    ! This routine is called after each call of solve_em so after each big timestep
    ! It progressively calculates average values of the prognostic variables in mass grid and higher order
    ! statistics if output_stats is higher
    ! variances and covariances must be calculated in post_process programs, for example :
    !     <u'v'> = <uv> - <u><v>
    ! This routine is called in med_after_solve_io (in share/mediation_integrate.F)
    ! The principle is based on WRFlux calculation but the calculation is different because I experienced
    ! wrong results in nested domains for mean velocities (U_MEAN, UX_MEAN, ...) with WRFlux. I couldn't find
    ! the origin of this error.
    IMPLICIT NONE

    TYPE(domain)             :: grid
    INTEGER, INTENT(IN   )   :: output_stats !0=no stats, 1=mean, 2=1+2nd order moments, 3=2+3rd order moments
    LOGICAL, INTENT(IN   )   :: output_dry_theta_fluxes !.true. = dry theta, .false. = moist theta
    INTEGER, INTENT(IN   )   :: its, ite, jts, jte, kts, kte !tile limit indices

    INTEGER :: i,j,k

    REAL :: local_count, count_p1, one_over_2g, one_over_count_p1,  &
    u,  & !destaggered u
    v,  & !destaggered v
    w,  & !destaggered w
    q,  & !qvapor
    t,  & !dry or moist theta perturbation (T-300K)
    p,  & !pressure perturbation (P-PB)
    z     !z perturbation (Z - PHB/g)

    if (output_stats > 0) then
      ! Note : These two lines are different than in upd_avgflx update because avgflx_count is incremented after
      !        the call of upd_avgflx (in solve_em.F) but stats_count before the call of upd_statistics
      count_p1 = real(grid%stats_count)
      local_count = count_p1 - 1.0

      one_over_2g = 0.5 / g
      one_over_count_p1 = 1.0 / count_p1

      DO j=jts,jte
      DO i=its,ite
        grid%mut_avg(i,j) = (local_count*grid%mut_avg(i,j)  +  grid%mut(i,j) ) * one_over_count_p1
      ENDDO
      ENDDO

      DO j=jts,jte
      DO k=kts,kte
      DO i=its,ite
        q = grid%moist(i,k,j,P_QV)
        w = (grid%w_2(i,k,j)+grid%w_2(i,k+1,j)) * 0.5
        z = (grid%ph_2(i,k,j)+grid%ph_2(i,k+1,j)) * one_over_2g
        p = grid%p(i,k,j)
#if (EM_CORE == 1)
        u = ((grid%u_2(i,k,j) + grid%u_2(i+1,k,j))*grid%cosa(i,j)  &
            -(grid%v_2(i,k,j) + grid%v_2(i+1,k,j))*grid%sina(i,j)) * 0.5
        v = ((grid%v_2(i,k,j) + grid%v_2(i,k,j+1))*grid%cosa(i,j)  &
            +(grid%u_2(i,k,j) + grid%u_2(i,k,j+1))*grid%sina(i,j)) * 0.5
        IF (output_dry_theta_fluxes) THEN
          IF (grid%use_theta_m == 1) THEN
            t = (grid%t_2(i,k,j)+T0) / (1.+R_v/R_d*q) - T0
          ELSE
            t = grid%t_2(i,k,j)
          ENDIF
        ELSE
          IF (grid%use_theta_m == 1) THEN
            t = grid%t_2(i,k,j)
          ELSE
            t = (grid%t_2(i,k,j)+T0) * (1.+R_v/R_d*q) - T0
          ENDIF
        ENDIF
#else
        u = grid%u(i,k,j)
        v = grid%v(i,k,j)
        t = grid%t(i,k,j) !I don't know if it is dry or moist
#endif
        !The calculation of time-average values is basic :
        ! <x>_n = sum_{i=0}^{i=n}(x_i) / n
        ! Then :
        ! <x>_0 = 0.                                      -> computed in zero_statistics
        ! <x>_{n+1} = (n*<x>_n + x_{n+1}) / (n+1)         -> computed in upd_statistics (here)
        grid%u_avg(i,k,j) = (local_count*grid%u_avg(i,k,j)  +  u ) * one_over_count_p1
        grid%v_avg(i,k,j) = (local_count*grid%v_avg(i,k,j)  +  v ) * one_over_count_p1
        grid%w_avg(i,k,j) = (local_count*grid%w_avg(i,k,j)  +  w ) * one_over_count_p1
        grid%q_avg(i,k,j) = (local_count*grid%q_avg(i,k,j)  +  q ) * one_over_count_p1
        grid%t_avg(i,k,j) = (local_count*grid%t_avg(i,k,j)  +  t ) * one_over_count_p1
        grid%p_avg(i,k,j) = (local_count*grid%p_avg(i,k,j)  +  p ) * one_over_count_p1
        grid%z_avg(i,k,j) = (local_count*grid%z_avg(i,k,j)  +  z ) * one_over_count_p1
        if (output_stats > 1) then
          grid%u2_avg(i,k,j) = (local_count*grid%u2_avg(i,k,j)  +  u*u ) * one_over_count_p1
          grid%uv_avg(i,k,j) = (local_count*grid%uv_avg(i,k,j)  +  u*v ) * one_over_count_p1
          grid%uw_avg(i,k,j) = (local_count*grid%uw_avg(i,k,j)  +  u*w ) * one_over_count_p1
          grid%uq_avg(i,k,j) = (local_count*grid%uq_avg(i,k,j)  +  u*q ) * one_over_count_p1
          grid%ut_avg(i,k,j) = (local_count*grid%ut_avg(i,k,j)  +  u*t ) * one_over_count_p1
          grid%v2_avg(i,k,j) = (local_count*grid%v2_avg(i,k,j)  +  v*v ) * one_over_count_p1
          grid%vw_avg(i,k,j) = (local_count*grid%vw_avg(i,k,j)  +  v*w ) * one_over_count_p1
          grid%vq_avg(i,k,j) = (local_count*grid%vq_avg(i,k,j)  +  v*q ) * one_over_count_p1
          grid%vt_avg(i,k,j) = (local_count*grid%vt_avg(i,k,j)  +  v*t ) * one_over_count_p1
          grid%w2_avg(i,k,j) = (local_count*grid%w2_avg(i,k,j)  +  w*w ) * one_over_count_p1
          grid%wq_avg(i,k,j) = (local_count*grid%wq_avg(i,k,j)  +  w*q ) * one_over_count_p1
          grid%wt_avg(i,k,j) = (local_count*grid%wt_avg(i,k,j)  +  w*t ) * one_over_count_p1
          grid%q2_avg(i,k,j) = (local_count*grid%q2_avg(i,k,j)  +  q*q ) * one_over_count_p1
          grid%qt_avg(i,k,j) = (local_count*grid%qt_avg(i,k,j)  +  q*t ) * one_over_count_p1
          grid%t2_avg(i,k,j) = (local_count*grid%t2_avg(i,k,j)  +  t*t ) * one_over_count_p1
          if (output_stats > 2) then
            grid%u3_avg(i,k,j) = (local_count*grid%u3_avg(i,k,j)  +  u*u*u ) * one_over_count_p1
            grid%v3_avg(i,k,j) = (local_count*grid%v3_avg(i,k,j)  +  v*v*v ) * one_over_count_p1
            grid%w3_avg(i,k,j) = (local_count*grid%w3_avg(i,k,j)  +  w*w*w ) * one_over_count_p1
            grid%q3_avg(i,k,j) = (local_count*grid%q3_avg(i,k,j)  +  q*q*q ) * one_over_count_p1
            grid%t3_avg(i,k,j) = (local_count*grid%t3_avg(i,k,j)  +  t*t*t ) * one_over_count_p1
          endif
        endif
      ENDDO
      ENDDO
      ENDDO
    endif
    return
  END SUBROUTINE upd_statistics

end MODULE module_statistics
