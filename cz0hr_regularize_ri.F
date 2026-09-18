!===================================================================
! SUBROUTINE: cz0hr_regularize_ri
! Purpose: Upstream Curvature-Aware Zero-Offset Hyperbolic Regularization
!          (C-Z0HR) of the gradient Richardson number profile BEFORE it is
!          consumed by the MYNN stability closures.
!          The MYNN Sm/Sh closure algebra remains unmodified.
!===================================================================
  PURE SUBROUTINE cz0hr_regularize_ri(kts, kte, dz, ri)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: kts, kte

    REAL(kind_phys), INTENT(IN), &
         DIMENSION(kts:kte) :: dz

    REAL(kind_phys), INTENT(INOUT), &
         DIMENSION(kts:kte) :: ri

    INTEGER :: k, n_levels

    REAL(kind_phys) :: ri_zz
    REAL(kind_phys) :: c_scale
    REAL(kind_phys) :: x
    REAL(kind_phys) :: phi_eps
    REAL(kind_phys) :: ratio
    REAL(kind_phys) :: dz_up
    REAL(kind_phys) :: dz_dn

    REAL(kind_phys), DIMENSION(kts:kte) :: ri_raw

    n_levels = kte - kts + 1

    ! Require at least one interior stencil.
    IF (n_levels < 3) RETURN

    ri_raw(kts:kte) = ri(kts:kte)

    ! Preserve boundary values.
    ri(kts) = ri_raw(kts)
    ri(kte) = ri_raw(kte)

    DO k = kts + 1, kte - 1

       ! Wall k lies between mass levels k-1 and k.
       dz_dn = dz(k-1)
       dz_up = dz(k)

       !------------------------------------------------------------
       ! Second-order pointwise curvature on a nonuniform grid:
       !
       ! R_zz(k) =
       ! 2 [ dz_dn R(k+1)
       !     -(dz_dn+dz_up) R(k)
       !     + dz_up R(k-1) ]
       ! -----------------------------------------------
       ! dz_dn dz_up (dz_dn+dz_up)
       !------------------------------------------------------------
       ri_zz = 2.0_kind_phys * ( &
                   dz_dn * ri_raw(k+1) &
                 - (dz_dn + dz_up) * ri_raw(k) &
                 + dz_up * ri_raw(k-1) ) &
                 / (dz_dn * dz_up * (dz_dn + dz_up))

       ! Grid-resolution curvature scale.
       ! dz_dn is retained deliberately for the upstream formulation.
       c_scale = 0.5_kind_phys * ABS(ri_zz) * dz_dn**2 &
               + cz0hr_eps_c * cz0hr_ric

       ! C-infinity hyperbolic distance from Ri_c.
       x       = ri_raw(k) - cz0hr_ric
       phi_eps = SQRT(x**2 + cz0hr_eps**2)

       ! Upstream C-Z0HR mapping.
       ratio = (phi_eps + c_scale) / &
               (phi_eps + (1.0_kind_phys + cz0hr_alpha) * c_scale)

       ri(k) = cz0hr_ric + x * ratio

    END DO

  END SUBROUTINE cz0hr_regularize_ri