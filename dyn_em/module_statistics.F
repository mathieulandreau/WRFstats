!WRF:MODEL_LAYER:DYNAMICS
!
MODULE module_avgflx_em

  USE module_bc
  USE module_model_constants
  USE module_wrf_error
  USE module_state_description

CONTAINS

!-------------------------------------------------------------------------------


  SUBROUTINE zero_statistics(grid, output_stats)

    IMPLICIT NONE

    TYPE(domain)             :: grid
    INTEGER, INTENT(IN   )   :: output_stats

    INTEGER        ::      &
      ids, ide, jds, jde, kds, kde,        &
      ims, ime, jms, jme, kms, kme,        &
      ips, ipe, jps, jpe, kps, kpe,        &
      imsx, imex, jmsx, jmex, kmsx, kmex,  &
      ipsx, ipex, jpsx, jpex, kpsx, kpex,  &
      imsy, imey, jmsy, jmey, kmsy, kmey,  &
      ipsy, ipey, jpsy, jpey, kpsy, kpey

    INTEGER :: i,j,k

    CALL get_ijk_from_grid (  grid ,                               &
                              ids, ide, jds, jde, kds, kde,        &
                              ims, ime, jms, jme, kms, kme,        &
                              ips, ipe, jps, jpe, kps, kpe,        &
                              imsx, imex, jmsx, jmex, kmsx, kmex,  &
                              ipsx, ipex, jpsx, jpex, kpsx, kpex,  &
                              imsy, imey, jmsy, jmey, kmsy, kmey,  &
                              ipsy, ipey, jpsy, jpey, kpsy, kpey )

    if (config_flags%output_stats > 0) then
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
      if (config_flags%output_stats > 1) then
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
        if (config_flags%output_stats > 2) then
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

  SUBROUTINE upd_statistics(grid, output_stats)

    IMPLICIT NONE

    TYPE(domain)             :: grid
    INTEGER, INTENT(IN   )   :: output_stats

    INTEGER        ::      &
      ids, ide, jds, jde, kds, kde,        &
      ims, ime, jms, jme, kms, kme,        &
      ips, ipe, jps, jpe, kps, kpe,        &
      imsx, imex, jmsx, jmex, kmsx, kmex,  &
      ipsx, ipex, jpsx, jpex, kpsx, kpex,  &
      imsy, imey, jmsy, jmey, kmsy, kmey,  &
      ipsy, ipey, jpsy, jpey, kpsy, kpey

    INTEGER :: i,j,k

    REAL :: local_count, count_p1, one_over_2g, one_over_count_p1,  &
    u,  & !destaggered u
    v,  & !destaggered v
    w,  & !destaggered w
    q,  & !qvapor
    t,  & !dry theta perturbation (T-300K)
    p,  & !pressure perturbation (P-PB)
    z   & !z perturbation (Z - PHB/g)

    CALL get_ijk_from_grid (  grid ,                               &
                              ids, ide, jds, jde, kds, kde,        &
                              ims, ime, jms, jme, kms, kme,        &
                              ips, ipe, jps, jpe, kps, kpe,        &
                              imsx, imex, jmsx, jmex, kmsx, kmex,  &
                              ipsx, ipex, jpsx, jpex, kpsx, kpex,  &
                              imsy, imey, jmsy, jmey, kmsy, kmey,  &
                              ipsy, ipey, jpsy, jpey, kpsy, kpey )

    local_count = real(grid%avgflx_count)
    count_p1 = local_count + 1.0

    one_over_2g = 0.5/9.81
    one_over_count_p1 = 1.0 / count_p1

    if (output_stats > 0) then

      DO j=jts,jte
      DO i=its,ite
        grid%mut_avg(i,j) = (local_count*grid%mut_avg(i,j)  +  grid%mut(i,j) ) * one_over_count_p1
      ENDDO
      ENDDO

      DO j=jts,jte
      DO k=kts,kte
      DO i=its,ite
#if (EM_CORE == 1)
        u = ((grid%u_2(i,k,j)*grid%cosa(i,j)-grid%v_2(i,k,j)*grid%sina(i,j)) + &
        (grid%u_2(i+1,k,j)*grid%cosa(i+1,j)-grid%v_2(i+1,k,j)*grid%sina(i+1,j))) * 0.5
        v = ((grid%v_2(i,k,j)*grid%cosa(i,j)+grid%u_2(i,k,j)*grid%sina(i,j)) + &
        (grid%v_2(i,k,j+1)*grid%cosa(i,j+1)+grid%u_2(i,k,j+1)*grid%sina(i,j+1))) * 0.5
        IF (grid%use_theta_m == 1) THEN
          t = (grid%t_2(i,1,j)+T0)/(1.+R_v/R_d*grid%moist(i,1,j,P_QV)) - T0
        ELSE
          t =  grid%t_2(i,1,j)
        END IF
#else
        u = grid%u(i,1,j)
        v = grid%v(i,1,j)
        t = grid%t(i,1,j)
#endif
        q = grid%moist(i,1,j,P_QV)
        w = (grid%w_2(i,k,j)+grid%w_2(i,k+1,j)) * 0.5
        z = (grid%ph_2(i,k,j)+grid%ph_2(i,k+1,j)) * one_over_2g
        p = grid%p(i,k,j)

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
    end if
    return
  END SUBROUTINE upd_statistics

  subroutine stagger_z(var_stag, var,                      &
                       cf1, cf2, cf3, cfn, cfn1, fnm, fnp, &
                       k, kms, kme, kde                    )

    REAL, INTENT(OUT) :: var_stag
    REAL, INTENT(IN) :: cf1, cf2, cf3, cfn, cfn1
    INTEGER, INTENT(IN) ::  k, kms, kme, kde
    REAL, DIMENSION(kms:kme) , INTENT(IN) :: var, fnm, fnp

    if ( k .eq. kde) then
      var_stag = cfn*var(kde-1)+cfn1*var(kde-2)
    elseif ( k .eq. 1) then
      var_stag = cf1*var(1)+cf2*var(2)+cf3*var(3)
    else
      var_stag = fnm(k)*var(k)+ fnp(k)*var(k-1)
    endif

  end subroutine stagger_z

end MODULE module_avgflx_em
