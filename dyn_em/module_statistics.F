!WRF:MODEL_LAYER:DYNAMICS
!
MODULE module_statistics

  USE module_model_constants, ONLY : R_v, R_d, T0, g, rcp, p0
  USE module_domain, ONLY : domain
  USE module_state_description, ONLY : P_QV, P_QC, P_QR

CONTAINS

!-------------------------------------------------------------------------------
  !The calculation of time-average values is basic :
  ! <x>_n = \sum_{i=0}^{i=n}(x_i) / n
  ! Then :
  ! <x>_0 = 0.                                      -> computed in zero_statistics
  ! <x>_{n+1} = (n*<x>_n + x_{n+1}) / (n+1)         -> computed in upd_statistics

  SUBROUTINE zero_statistics(grid, output_stats, its, ite, jts, jte, kts, kte)

    IMPLICIT NONE

    TYPE(domain)             :: grid
    INTEGER, INTENT(IN   )   :: output_stats
    INTEGER, INTENT(IN   )   :: its, ite, jts, jte, kts, kte !tile limit indices

    INTEGER :: i,j,k
    if (output_stats > 0) then
      DO j=jts,jte
      DO i=its,ite
        grid%mu_avg(i,j)   = 0.
        grid%psfc_avg(i,j) = 0.
        grid%q2m_avg(i,j)  = 0.
        grid%t2m_avg(i,j)  = 0.
        grid%pt2m_avg(i,j) = 0.
        grid%u10m_avg(i,j) = 0.
        grid%v10m_avg(i,j) = 0.
        grid%tsk_avg(i,j)  = 0.
        grid%ust_avg(i,j)  = 0.
        grid%pblh_avg(i,j) = 0.
      ENDDO
      ENDDO
      DO j=jts,jte
      DO k=kts,kte
      DO i=its,ite
        grid%zp_avg(i,k,j) = 0.
        grid%u_avg(i,k,j) = 0.
        grid%v_avg(i,k,j) = 0.
        grid%w_avg(i,k,j) = 0.
        grid%pt_avg(i,k,j) = 0.
        grid%ptv_avg(i,k,j) = 0.
        grid%qv_avg(i,k,j) = 0.
        grid%q_avg(i,k,j) = 0.
        grid%p_avg(i,k,j) = 0.
        grid%tke_avg(i,k,j) = 0.
        grid%t_avg(i,k,j) = 0.
        grid%rhod_avg(i,k,j) = 0.
        grid%rhov_avg(i,k,j) = 0.
        grid%rho_avg(i,k,j) = 0.
      END DO
      END DO
      END DO
      if (output_stats > 1) then
        DO j=jts,jte
        DO k=kts,kte
        DO i=its,ite
          grid%u2_avg(i,k,j) = 0.
          grid%v2_avg(i,k,j) = 0.
          grid%w2_avg(i,k,j) = 0.
          grid%pt2_avg(i,k,j) = 0.
          grid%ptv2_avg(i,k,j) = 0.
          grid%q2_avg(i,k,j) = 0.
          grid%qv2_avg(i,k,j) = 0.
          grid%zp2_avg(i,k,j) = 0.

          grid%uv_avg(i,k,j) = 0.
          grid%uw_avg(i,k,j) = 0.
          grid%vw_avg(i,k,j) = 0.

          grid%uq_avg(i,k,j) = 0.
          grid%vq_avg(i,k,j) = 0.
          grid%wq_avg(i,k,j) = 0.

          grid%uqv_avg(i,k,j) = 0.
          grid%vqv_avg(i,k,j) = 0.
          grid%wqv_avg(i,k,j) = 0.

          grid%upt_avg(i,k,j) = 0.
          grid%vpt_avg(i,k,j) = 0.
          grid%wpt_avg(i,k,j) = 0.

          grid%uptv_avg(i,k,j) = 0.
          grid%vptv_avg(i,k,j) = 0.
          grid%wptv_avg(i,k,j) = 0.

          grid%up_avg(i,k,j) = 0.
          grid%vp_avg(i,k,j) = 0.
          grid%wp_avg(i,k,j) = 0.
        END DO
        END DO
        END DO
        if (output_stats > 2) then
          DO j=jts,jte
          DO k=kts,kte
          DO i=its,ite
            grid%uuu_avg(i,k,j) = 0.
            grid%vvv_avg(i,k,j) = 0.
            grid%www_avg(i,k,j) = 0.
            grid%uuv_avg(i,k,j) = 0.
            grid%uvv_avg(i,k,j) = 0.
            grid%uuw_avg(i,k,j) = 0.
            grid%uww_avg(i,k,j) = 0.
            grid%vvw_avg(i,k,j) = 0.
            grid%vww_avg(i,k,j) = 0.
            grid%uvw_avg(i,k,j) = 0.
          END DO
          END DO
          END DO
        endif
      endif
    endif
    
    return
  END SUBROUTINE zero_statistics

  SUBROUTINE upd_statistics(grid, output_stats, use_theta_m, its, ite, jts, jte, kts, kte)
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
    INTEGER, INTENT(IN   )   :: use_theta_m  !0=dry theta, 1=moist theta
    INTEGER, INTENT(IN   )   :: its, ite, jts, jte, kts, kte !tile limit indices

    INTEGER :: i,j,k, tke_opt

    REAL :: local_count, count_p1, one_over_2g, one_over_count_p1,  &
    u,    & !destaggered u
    v,    & !destaggered v
    w,    & !destaggered w
    q,    & !total mixing ratio (in fact only rain + cloud + vapor)
    qv,   & !qvapor
    pt,   & !dry theta perturbation (PT-300K)
    ptv,  & !moist theta perturbation (PTV-300K)
    p,    & !pressure perturbation (P-P_BASE)
    p_hyd,& !hydrostatic pressure
    zp,   & !zp perturbation (PHI - PHI_BASE)/g
    t,    & !temperature perturbation (T-300K)
    rhod, & !dry density
    rhov, & !water vapor density
    rho,  & !total density
    tke,  & !tke_pbl in RANS with MYJ scheme, tke_sfs for LES with TKE scheme, 0 otherwhile
    pb,   & !base pressure
    tv      !moist temperature

    if (output_stats > 0) then
      ! Note : These two lines are different than in upd_avgflx update because avgflx_count is incremented after
      !        the call of upd_avgflx (in solve_em.F) but stats_count before the call of upd_statistics
      count_p1 = real(grid%stats_count)
      local_count = count_p1 - 1.0

      one_over_2g = 0.5 / g
      one_over_count_p1 = 1.0 / count_p1

      ! some settings before loop to be sure compiler knows they are invariant in loop
      IF (SIZE(grid%tke_2,1) * SIZE(grid%tke_2,3) > 1) THEN
        tke_opt = 1
      ELSE IF (SIZE(grid%tke_pbl,1) * SIZE(grid%tke_pbl,3) > 1) THEN
        tke_opt = 2
      ELSE 
        tke_opt = 0
      END IF

      DO j=jts,jte
      DO i=its,ite
        grid%mu_avg(i,j)   = (local_count*grid%mu_avg(i,j)    +  grid%mu_2(i,j) ) * one_over_count_p1
        grid%psfc_avg(i,j) = (local_count*grid%psfc_avg(i,j)  +  grid%psfc(i,j) ) * one_over_count_p1
        grid%q2m_avg(i,j)  = (local_count*grid%q2m_avg(i,j)   +  grid%q2(i,j) )   * one_over_count_p1
        grid%t2m_avg(i,j)  = (local_count*grid%t2m_avg(i,j)   +  grid%t2(i,j) )   * one_over_count_p1
        grid%pt2m_avg(i,j) = (local_count*grid%pt2m_avg(i,j)  +  grid%th2(i,j) )  * one_over_count_p1
        grid%u10m_avg(i,j) = (local_count*grid%u10m_avg(i,j)  +  grid%u10(i,j) )  * one_over_count_p1
        grid%v10m_avg(i,j) = (local_count*grid%v10m_avg(i,j)  +  grid%v10(i,j) )  * one_over_count_p1
        grid%tsk_avg(i,j)  = (local_count*grid%tsk_avg(i,j)   +  grid%tsk(i,j) )  * one_over_count_p1
        grid%ust_avg(i,j)  = (local_count*grid%ust_avg(i,j)   +  grid%ust(i,j) )  * one_over_count_p1
        grid%pblh_avg(i,j) = (local_count*grid%pblh_avg(i,j)  +  grid%pblh(i,j) ) * one_over_count_p1
      ENDDO
      ENDDO

      DO j=jts,jte
      DO k=kts,kte
      DO i=its,ite
        qv = grid%moist(i,k,j,P_QV)
        q = qv + grid%moist(i,k,j,P_QC) + grid%moist(i,k,j,P_QR) ! We don't considere ice and snow here
        w = (grid%w_2(i,k,j)+grid%w_2(i,k+1,j)) * 0.5
        zp = (grid%ph_2(i,k,j)+grid%ph_2(i,k+1,j)) * one_over_2g
        p = grid%p(i,k,j)
        p_hyd = grid%p_hyd(i,k,j)
        IF (tke_opt == 1) THEN
          tke = grid%tke_2(i,k,j)
        ELSE IF (tke_opt == 2) THEN
          tke = grid%tke_pbl(i,k,j)
        ELSE 
          tke = 0
        END IF
#if (EM_CORE == 1)
        u = ((grid%u_2(i,k,j) + grid%u_2(i+1,k,j))*grid%cosa(i,j)  &
            -(grid%v_2(i,k,j) + grid%v_2(i+1,k,j))*grid%sina(i,j)) * 0.5
        v = ((grid%v_2(i,k,j) + grid%v_2(i,k,j+1))*grid%cosa(i,j)  &
            +(grid%u_2(i,k,j) + grid%u_2(i,k,j+1))*grid%sina(i,j)) * 0.5
        IF (use_theta_m == 1) THEN
          ptv = grid%t_2(i,k,j)
          pt = (ptv+T0) / (1.+R_v/R_d*qv) - T0
        ELSE
          pt = grid%t_2(i,k,j)
          ptv = (pt+T0) * (1.+R_v/R_d*qv) - T0
        ENDIF
#else
        u = grid%u(i,k,j)
        v = grid%v(i,k,j)
        pt = grid%pt(i,k,j)
        ptv = (pt+T0) * (1.+R_v/R_d*qv) - T0
#endif
        pb = grid%pb(i,k,j)
        t = (pt+T0) * (p0/(p+pb))**(-rcp)
        tv = (ptv+T0) * (p0/(p+pb))**(-rcp)
        rhod = (p+pb)/(R_d*tv)
        rhov = rhod*qv
        rho = rhod*(1.+q)
        t = t - T0
        !The calculation of time-average values is basic :
        ! <x>_n = sum_{i=0}^{i=n}(x_i) / n
        ! Then :
        ! <x>_0 = 0.                                      -> computed in zero_statistics
        ! <x>_{n+1} = (n*<x>_n + x_{n+1}) / (n+1)         -> computed in upd_statistics (here)
        grid%zp_avg(i,k,j)      = (local_count*grid%zp_avg(i,k,j)    +  zp )    * one_over_count_p1
        grid%u_avg(i,k,j)       = (local_count*grid%u_avg(i,k,j)     +  u )     * one_over_count_p1
        grid%v_avg(i,k,j)       = (local_count*grid%v_avg(i,k,j)     +  v )     * one_over_count_p1
        grid%w_avg(i,k,j)       = (local_count*grid%w_avg(i,k,j)     +  w )     * one_over_count_p1
        grid%pt_avg(i,k,j)      = (local_count*grid%pt_avg(i,k,j)    +  pt )    * one_over_count_p1
        grid%ptv_avg(i,k,j)     = (local_count*grid%ptv_avg(i,k,j)   +  ptv )   * one_over_count_p1
        grid%qv_avg(i,k,j)      = (local_count*grid%qv_avg(i,k,j)    +  qv )    * one_over_count_p1
        grid%q_avg(i,k,j)       = (local_count*grid%q_avg(i,k,j)     +  q )     * one_over_count_p1
        grid%p_avg(i,k,j)       = (local_count*grid%p_avg(i,k,j)     +  p )     * one_over_count_p1
        grid%p_hyd_avg(i,k,j)   = (local_count*grid%p_hyd_avg(i,k,j) +  p_hyd ) * one_over_count_p1
        grid%tke_avg(i,k,j)     = (local_count*grid%tke_avg(i,k,j)   +  tke )   * one_over_count_p1
        grid%rho_avg(i,k,j)     = (local_count*grid%rho_avg(i,k,j)   +  rho )   * one_over_count_p1
        grid%rhod_avg(i,k,j)    = (local_count*grid%rhod_avg(i,k,j)  +  rhod )  * one_over_count_p1
        grid%rhov_avg(i,k,j)    = (local_count*grid%rhov_avg(i,k,j)  +  rhov )  * one_over_count_p1
        grid%t_avg(i,k,j)       = (local_count*grid%t_avg(i,k,j)     +  t )     * one_over_count_p1
        if (output_stats > 1) then
          grid%u2_avg(i,k,j)   = (local_count*grid%u2_avg(i,k,j)   +  u*u )     * one_over_count_p1
          grid%v2_avg(i,k,j)   = (local_count*grid%v2_avg(i,k,j)   +  v*v )     * one_over_count_p1
          grid%w2_avg(i,k,j)   = (local_count*grid%w2_avg(i,k,j)   +  w*w )     * one_over_count_p1
          grid%pt2_avg(i,k,j)  = (local_count*grid%pt2_avg(i,k,j)  +  pt*pt )   * one_over_count_p1
          grid%ptv2_avg(i,k,j) = (local_count*grid%ptv2_avg(i,k,j) +  ptv*ptv ) * one_over_count_p1
          grid%q2_avg(i,k,j)   = (local_count*grid%q2_avg(i,k,j)   +  q*q )     * one_over_count_p1
          grid%qv2_avg(i,k,j)  = (local_count*grid%qv2_avg(i,k,j)  +  qv*qv )   * one_over_count_p1
          grid%zp2_avg(i,k,j)   = (local_count*grid%zp2_avg(i,k,j)   +  zp*zp )     * one_over_count_p1

          grid%uv_avg(i,k,j)   = (local_count*grid%uv_avg(i,k,j)   +  u*v )     * one_over_count_p1
          grid%uw_avg(i,k,j)   = (local_count*grid%uw_avg(i,k,j)   +  u*w )     * one_over_count_p1
          grid%vw_avg(i,k,j)   = (local_count*grid%vw_avg(i,k,j)   +  v*w )     * one_over_count_p1

          grid%uq_avg(i,k,j)   = (local_count*grid%uq_avg(i,k,j)   +  u*q )     * one_over_count_p1
          grid%vq_avg(i,k,j)   = (local_count*grid%vq_avg(i,k,j)   +  v*q )     * one_over_count_p1
          grid%wq_avg(i,k,j)   = (local_count*grid%wq_avg(i,k,j)   +  w*q )     * one_over_count_p1

          grid%uqv_avg(i,k,j)  = (local_count*grid%uqv_avg(i,k,j)  +  u*qv )    * one_over_count_p1
          grid%vqv_avg(i,k,j)  = (local_count*grid%vqv_avg(i,k,j)  +  v*qv )    * one_over_count_p1
          grid%wqv_avg(i,k,j)  = (local_count*grid%wqv_avg(i,k,j)  +  w*qv )    * one_over_count_p1

          grid%upt_avg(i,k,j)  = (local_count*grid%upt_avg(i,k,j)  +  u*pt )    * one_over_count_p1
          grid%vpt_avg(i,k,j)  = (local_count*grid%vpt_avg(i,k,j)  +  v*pt )    * one_over_count_p1
          grid%wpt_avg(i,k,j)  = (local_count*grid%wpt_avg(i,k,j)  +  w*pt )    * one_over_count_p1

          grid%uptv_avg(i,k,j) = (local_count*grid%uptv_avg(i,k,j) +  u*ptv )   * one_over_count_p1
          grid%vptv_avg(i,k,j) = (local_count*grid%vptv_avg(i,k,j) +  v*ptv )   * one_over_count_p1
          grid%wptv_avg(i,k,j) = (local_count*grid%wptv_avg(i,k,j) +  w*ptv )   * one_over_count_p1

          grid%up_avg(i,k,j)   = (local_count*grid%up_avg(i,k,j)   +  u*p )     * one_over_count_p1
          grid%vp_avg(i,k,j)   = (local_count*grid%vp_avg(i,k,j)   +  v*p )     * one_over_count_p1
          grid%wp_avg(i,k,j)   = (local_count*grid%wp_avg(i,k,j)   +  w*p )     * one_over_count_p1

          if (output_stats > 2) then
            grid%uuu_avg(i,k,j) = (local_count*grid%uuu_avg(i,k,j)  +  u*u*u ) * one_over_count_p1
            grid%vvv_avg(i,k,j) = (local_count*grid%vvv_avg(i,k,j)  +  v*v*v ) * one_over_count_p1
            grid%www_avg(i,k,j) = (local_count*grid%www_avg(i,k,j)  +  w*w*w ) * one_over_count_p1
            grid%uuv_avg(i,k,j) = (local_count*grid%uuv_avg(i,k,j)  +  u*u*v ) * one_over_count_p1
            grid%uvv_avg(i,k,j) = (local_count*grid%uvv_avg(i,k,j)  +  u*v*v ) * one_over_count_p1
            grid%uuw_avg(i,k,j) = (local_count*grid%uuw_avg(i,k,j)  +  u*u*w ) * one_over_count_p1
            grid%uww_avg(i,k,j) = (local_count*grid%uww_avg(i,k,j)  +  u*w*w ) * one_over_count_p1
            grid%vvw_avg(i,k,j) = (local_count*grid%vvw_avg(i,k,j)  +  v*v*w ) * one_over_count_p1
            grid%vww_avg(i,k,j) = (local_count*grid%vww_avg(i,k,j)  +  v*w*w ) * one_over_count_p1
            grid%uvw_avg(i,k,j) = (local_count*grid%uvw_avg(i,k,j)  +  u*v*w ) * one_over_count_p1
          endif
        endif
      ENDDO
      ENDDO
      ENDDO
    endif
    return
  END SUBROUTINE upd_statistics

end MODULE module_statistics
