     
      
      turb_ct(inot) = A(inot,FLOOR(simulated_time)+1_iwp)
   
   END SUBROUTINE wtm_dct_from_matlab
      
   SUBROUTINE wtm_forces_admr


       IMPLICIT NONE

       CHARACTER (LEN=2) ::  turbine_id

       INTEGER(iwp) ::  i, j, k          !< loop indices
       INTEGER(iwp) ::  inot             !< turbine loop index (turbine id)
       INTEGER(iwp) ::  iialpha, iir     !<
       INTEGER(iwp) ::  rseg, rseg_int   !<
       INTEGER(iwp) ::  ring, ring_int   !<
       INTEGER(iwp) ::  ii, jj, kk       !<
    
       REAL(wp)     ::  sin_rot, cos_rot   !<
       REAL(wp)     ::  sin_yaw, cos_yaw   !<
       
       REAL(wp) ::  aa, bb, cc, dd  !< interpolation distances
       REAL(wp) ::  gg              !< interpolation volume var  
       
       REAL(wp) ::  dist_u_3d, dist_v_3d, dist_w_3d  !< smearing distances

       
!
!      Variables for pitch control
       REAL(wp)     ::  torque_max=0.0_wp
       LOGICAL      ::  pitch_sw=.FALSE.

       INTEGER(iwp), DIMENSION(1) :: lct=0
       REAL(wp), DIMENSION(1)     :: rad_d=0.0_wp


       CALL cpu_log( log_point_s(61), 'wtm_forces', 'start' )



!
!--       Set forces to zero for each new time step:
          thrust(:,:,:)         = 0.0_wp
          torque_y(:,:,:)       = 0.0_wp
          torque_z(:,:,:)       = 0.0_wp
          torque_total(:)       = 0.0_wp
          rot_tend_x(:,:,:)     = 0.0_wp
          rot_tend_y(:,:,:)     = 0.0_wp
          rot_tend_z(:,:,:)     = 0.0_wp
          thrust_rotor(:)       = 0.0_wp
!
!--       Loop over number of turbines:
          DO inot = 1, nturbines

             cos_yaw = COS(phi_yaw(inot))
             sin_yaw = SIN(phi_yaw(inot))
!
!--          Loop over rings of each turbine:
             DO ring = 1, nrings(inot)

                thrust_seg(:)   = 0.0_wp
                torque_seg_y(:) = 0.0_wp
                torque_seg_z(:) = 0.0_wp
!
!--             Determine distance between each ring (center) and the hub:
                cur_r = (ring - 0.5_wp) * delta_r(inot)

!
!--             Loop over segments of each ring of each turbine:
                DO rseg = 1, nsegs(ring,inot)
!
!--                !-----------------------------------------------------------!
!--                !-- Determine coordinates of the ring segments            --!
!--                !-----------------------------------------------------------!
!
!--                Determine angle of ring segment towards zero degree angle of
!--                rotor system (at zero degree rotor direction vectors aligned
!--                with y-axis):
                   phi_rotor = rseg * 2.0_wp * pi / nsegs(ring,inot)
                   cos_rot   = COS( phi_rotor )
                   sin_rot   = SIN( phi_rotor )

!--                Now the direction vectors can be determined with respect to
!--                the yaw and tilt angle:
                   re(1) = cos_rot * sin_yaw
                   re(2) = cos_rot * cos_yaw    
                   re(3) = sin_rot

                   rote = MATMUL( rot_coord_trans(inot,:,:), re )
!
!--                Coordinates of the single segments (center points):
                   rbx(ring,rseg) = rcx(inot) + cur_r * rote(1)
                   rby(ring,rseg) = rcy(inot) + cur_r * rote(2)
                   rbz(ring,rseg) = rcz(inot) + cur_r * rote(3)

!--                !-----------------------------------------------------------!
!--                !-- Interpolation of the velocity components from the     --!
!--                !-- cartesian grid point to the coordinates of each ring  --!
!--                !-- segment (follows a method used in the particle model) --!
!--                !-----------------------------------------------------------!

                   u_int(inot,ring,rseg)     = 0.0_wp
                   u_int_1_l(inot,ring,rseg) = 0.0_wp

                   v_int(inot,ring,rseg)     = 0.0_wp
                   v_int_1_l(inot,ring,rseg) = 0.0_wp

                   w_int(inot,ring,rseg)     = 0.0_wp
                   w_int_1_l(inot,ring,rseg) = 0.0_wp

!
!--                Interpolation of the u-component:

                   ii =   rbx(ring,rseg) * ddx
                   jj = ( rby(ring,rseg) - 0.5_wp * dy ) * ddy
                   kk = ( rbz(ring,rseg) - 0.5_wp * dz ) / dz
!
!--                Interpolate only if all required information is available on
!--                the current PE:
                   IF ( ( ii >= nxl )  .AND.  ( ii <= nxr ) )  THEN
                      IF ( ( jj >= nys )  .AND.  ( jj <= nyn ) )  THEN

                         aa = ( ( ii + 1          ) * dx - rbx(ring,rseg) ) *  &
                              ( ( jj + 1 + 0.5_wp ) * dy - rby(ring,rseg) )
                         bb = ( rbx(ring,rseg) - ii * dx )                  *  &
                              ( ( jj + 1 + 0.5_wp ) * dy - rby(ring,rseg) )
                         cc = ( ( ii+1            ) * dx - rbx(ring,rseg) ) *  &
                              ( rby(ring,rseg) - ( jj + 0.5_wp ) * dy )
                         dd = ( rbx(ring,rseg) -              ii * dx )     *  &
                              ( rby(ring,rseg) - ( jj + 0.5_wp ) * dy ) 
                         gg = dx * dy

                         u_int_l = ( aa * u(kk,jj,ii)     +                    &
                                     bb * u(kk,jj,ii+1)   +                    &
                                     cc * u(kk,jj+1,ii)   +                    &
                                     dd * u(kk,jj+1,ii+1)                      &
                                   ) / gg

                         u_int_u = ( aa * u(kk+1,jj,ii)     +                  &
                                     bb * u(kk+1,jj,ii+1)   +                  &
                                     cc * u(kk+1,jj+1,ii)   +                  &
                                     dd * u(kk+1,jj+1,ii+1)                    &
                                   ) / gg

                         u_int_1_l(inot,ring,rseg) = u_int_l          +        &
                                     ( rbz(ring,rseg) - zu(kk) ) / dz *        &
                                     ( u_int_u - u_int_l )

                      ELSE 
                         u_int_1_l(inot,ring,rseg) = 0.0_wp
                      ENDIF
                   ELSE
                      u_int_1_l(inot,ring,rseg) = 0.0_wp
                   ENDIF


!
!--                Interpolation of the v-component:
                   ii = ( rbx(ring,rseg) - 0.5_wp * dx ) * ddx
                   jj =   rby(ring,rseg)                 * ddy
                   kk = ( rbz(ring,rseg) + 0.5_wp * dz ) / dz 
!
!--                Interpolate only if all required information is available on
!--                the current PE:
                   IF ( ( ii >= nxl )  .AND.  ( ii <= nxr ) )  THEN
                      IF ( ( jj >= nys )  .AND.  ( jj <= nyn ) )  THEN

                         aa = ( ( ii + 1 + 0.5_wp ) * dx - rbx(ring,rseg) ) *  &
                              ( ( jj + 1 )          * dy - rby(ring,rseg) )
                         bb = ( rbx(ring,rseg)     - ( ii + 0.5_wp ) * dx ) *  &
                              ( ( jj + 1 ) * dy          - rby(ring,rseg) )
                         cc = ( ( ii + 1 + 0.5_wp ) * dx - rbx(ring,rseg) ) *  &
                              ( rby(ring,rseg)           -        jj * dy )
                         dd = ( rbx(ring,rseg)     - ( ii + 0.5_wp ) * dx ) *  &
                              ( rby(ring,rseg)           -        jj * dy )
                         gg = dx * dy

                         v_int_l = ( aa * v(kk,jj,ii)     +                    &
                                     bb * v(kk,jj,ii+1)   +                    &
                                     cc * v(kk,jj+1,ii)   +                    &
                                     dd * v(kk,jj+1,ii+1)                      &
                                   ) / gg

                         v_int_u = ( aa * v(kk+1,jj,ii)     +                  &
                                     bb * v(kk+1,jj,ii+1)   +                  &
                                     cc * v(kk+1,jj+1,ii)   +                  &
                                     dd * v(kk+1,jj+1,ii+1)                    &
                                  ) / gg

                         v_int_1_l(inot,ring,rseg) = v_int_l +                 &
                                     ( rbz(ring,rseg) - zu(kk) ) / dz *        &
                                     ( v_int_u - v_int_l )

                      ELSE
                         v_int_1_l(inot,ring,rseg) = 0.0_wp
                      ENDIF
                   ELSE
                      v_int_1_l(inot,ring,rseg) = 0.0_wp
                   ENDIF


!
!--                Interpolation of the w-component:
                   ii = ( rbx(ring,rseg) - 0.5_wp * dx ) * ddx
                   jj = ( rby(ring,rseg) - 0.5_wp * dy ) * ddy
                   kk =   rbz(ring,rseg)                 / dz
!
!--                Interpolate only if all required information is available on
!--                the current PE:
                   IF ( ( ii >= nxl )  .AND.  ( ii <= nxr ) )  THEN
                      IF ( ( jj >= nys )  .AND.  ( jj <= nyn ) )  THEN

                         aa = ( ( ii + 1 + 0.5_wp ) * dx - rbx(ring,rseg) ) *  &
                              ( ( jj + 1 + 0.5_wp ) * dy - rby(ring,rseg) )
                         bb = ( rbx(ring,rseg)     - ( ii + 0.5_wp ) * dx ) *  &
                              ( ( jj + 1 + 0.5_wp ) * dy - rby(ring,rseg) )
                         cc = ( ( ii + 1 + 0.5_wp ) * dx - rbx(ring,rseg) ) *  &
                              ( rby(ring,rseg)     - ( jj + 0.5_wp ) * dy )
                         dd = ( rbx(ring,rseg)     - ( ii + 0.5_wp ) * dx ) *  &
                              ( rby(ring,rseg)     - ( jj + 0.5_wp ) * dy )
                         gg = dx * dy

                         w_int_l = ( aa * w(kk,jj,ii)     +                    &
                                     bb * w(kk,jj,ii+1)   +                    &
                                     cc * w(kk,jj+1,ii)   +                    &
                                     dd * w(kk,jj+1,ii+1)                      &
                                   ) / gg

                         w_int_u = ( aa * w(kk+1,jj,ii)     +                  &
                                     bb * w(kk+1,jj,ii+1)   +                  &
                                     cc * w(kk+1,jj+1,ii)   +                  &
                                     dd * w(kk+1,jj+1,ii+1)                    &
                                    ) / gg

                         w_int_1_l(inot,ring,rseg) = w_int_l +                 &
                                     ( rbz(ring,rseg) - zw(kk) ) / dz *        &
                                     ( w_int_u - w_int_l )
                      ELSE
                         w_int_1_l(inot,ring,rseg) = 0.0_wp
                      ENDIF
                   ELSE
                      w_int_1_l(inot,ring,rseg) = 0.0_wp
                   ENDIF

                ENDDO
             ENDDO

          ENDDO

!
!--       Exchange between PEs (information required on each PE):
#if defined( __parallel )
          CALL MPI_ALLREDUCE( u_int_1_l, u_int, nturbines * MAXVAL(nrings) *   &
                              MAXVAL(nsegs), MPI_REAL, MPI_SUM, comm2d, ierr )
          CALL MPI_ALLREDUCE( v_int_1_l, v_int, nturbines * MAXVAL(nrings) *   &
                              MAXVAL(nsegs), MPI_REAL, MPI_SUM, comm2d, ierr )
          CALL MPI_ALLREDUCE( w_int_1_l, w_int, nturbines * MAXVAL(nrings) *   &
                              MAXVAL(nsegs), MPI_REAL, MPI_SUM, comm2d, ierr )
#else
          u_int = u_int_1_l
          v_int = v_int_1_l
          w_int = w_int_1_l
#endif


!
!--       Loop over number of turbines:

          DO inot = 1, nturbines
pit_loop: DO

             IF ( pitch_sw )  THEN
                torque_total(inot) = 0.0_wp
                thrust_rotor(inot) = 0.0_wp
                pitch_add(inot)    = pitch_add(inot) + 0.25_wp
!                 IF ( myid == 0 ) PRINT*, 'Pitch', inot, pitch_add(inot)
             ELSE
                cos_yaw = COS(phi_yaw(inot))
                sin_yaw = SIN(phi_yaw(inot))
                IF ( pitch_control )  THEN
                   pitch_add(inot) = MAX(pitch_add_old(inot) - pitch_rate *    &
                                         dt_3d , 0.0_wp )
                ENDIF
             ENDIF

!
!--          Loop over rings of each turbine:
             DO ring = 1, nrings(inot)
!
!--             Determine distance between each ring (center) and the hub:
                cur_r = (ring - 0.5_wp) * delta_r(inot)
!
!--             Loop over segments of each ring of each turbine:
                DO rseg = 1, nsegs(ring,inot)
!
!--                Determine angle of ring segment towards zero degree angle of
!--                rotor system (at zero degree rotor direction vectors aligned
!--                with y-axis):
                   phi_rotor = rseg * 2.0_wp * pi / nsegs(ring,inot)
                   cos_rot   = COS(phi_rotor)
                   sin_rot   = SIN(phi_rotor)
!
!--                Now the direction vectors can be determined with respect to
!--                the yaw and tilt angle:
                   re(1) = cos_rot * sin_yaw
                   re(2) = cos_rot * cos_yaw
                   re(3) = sin_rot

!                  The current unit vector in azimuthal direction:                         
                   rea(1) = - sin_rot * sin_yaw
                   rea(2) = - sin_rot * cos_yaw
                   rea(3) =   cos_rot

!
!--                To respect the yawing angle for the calculations of
!--                velocities and forces the unit vectors perpendicular to the
!--                rotor area in direction of the positive yaw angle are defined:
                   ren(1) =   cos_yaw
                   ren(2) = - sin_yaw
                   ren(3) = 0.0_wp
!
!--                Multiplication with the coordinate transformation matrix
!--                gives the final unit vector with consideration of the rotor
!--                tilt:
                   rote = MATMUL( rot_coord_trans(inot,:,:), re )
                   rota = MATMUL( rot_coord_trans(inot,:,:), rea )
                   rotn = MATMUL( rot_coord_trans(inot,:,:), ren )
!
!--                Coordinates of the single segments (center points):
                   rbx(ring,rseg) = rcx(inot) + cur_r * rote(1)

                   rby(ring,rseg) = rcy(inot) + cur_r * rote(2)

                   rbz(ring,rseg) = rcz(inot) + cur_r * rote(3)

!
!--                !-----------------------------------------------------------!
!--                !-- Calculation of various angles and relative velocities --!
!--                !-----------------------------------------------------------!
!
!--                In the following the 3D-velocity field is projected its
!--                components perpedicular and parallel to the rotor area
!--                The calculation of forces will be done in the rotor-
!--                coordinates y' and z.
!--                The yaw angle will be reintroduced when the force is applied
!--                on the hydrodynamic equations
!
!--                Projection of the xy-velocities relative to the rotor area
!
!--                Velocity perpendicular to the rotor area:
                   u_rot = u_int(inot,ring,rseg)*rotn(1) +                     &
                   v_int(inot,ring,rseg)*rotn(2) +                             &
                   w_int(inot,ring,rseg)*rotn(3)
!
!--                Projection of the 3D-velocity vector in the azimuthal
!--                direction:
                   vtheta(rseg) = rota(1) * u_int(inot,ring,rseg) +            & 
                                  rota(2) * v_int(inot,ring,rseg) +            &
                                  rota(3) * w_int(inot,ring,rseg)
!
!--                Determination of the angle phi_rel between the rotor plane
!--                and the direction of the flow relative to the rotor:

                   phi_rel(rseg) = ATAN( u_rot /                               &
                                         ( omega_rot(inot) * cur_r -           &
                                           vtheta(rseg) ) )

!
!--                Interpolation of the local pitch angle from tabulated values
!--                to the current radial position:

                   lct=minloc(ABS(cur_r-lrd))
                   rad_d=cur_r-lrd(lct)
                   
                   IF (cur_r == 0.0_wp) THEN
                      alpha_attack(rseg) = 0.0_wp
                   ELSE IF (cur_r >= lrd(size(ard))) THEN
                      alpha_attack(rseg) = ( ard(size(ard)) +                  &
                                             ard(size(ard)-1) ) / 2.0_wp
                   ELSE 
                      alpha_attack(rseg) = ( ard(lct(1)) *  &
                                             ( ( lrd(lct(1)+1) - cur_r ) /     &
                                               ( lrd(lct(1)+1) - lrd(lct(1)) ) &
                                             ) ) + ( ard(lct(1)+1) *           &
                                             ( ( cur_r - lrd(lct(1)) ) /       &
                                               ( lrd(lct(1)+1) - lrd(lct(1)) ) ) )
                   ENDIF

!
!--                In Fortran radian instead of degree is used as unit for all
!--                angles. Therefore, a transformation from angles given in
!--                degree to angles given in radian is necessary here:
                   alpha_attack(rseg) = alpha_attack(rseg) *                   &
                                        ( (2.0_wp*pi) / 360.0_wp )
!
!--                Substraction of the local pitch angle to obtain the local
!--                angle of attack:
                   alpha_attack(rseg) = phi_rel(rseg) - alpha_attack(rseg)
!
!--                Preliminary transformation back from angles given in radian
!--                to angles given in degree:
                   alpha_attack(rseg) = alpha_attack(rseg) *                   &
                                        ( 360.0_wp / (2.0_wp*pi) )
!
!--                Correct with collective pitch angle:
                   alpha_attack = alpha_attack - pitch_add(inot)

!
!--                Determination of the magnitude of the flow velocity relative
!--                to the rotor:
                   vrel(rseg) = SQRT( u_rot**2 +                               &
                                      ( omega_rot(inot) * cur_r -              &
                                        vtheta(rseg) )**2 )

!
!--                !-----------------------------------------------------------!
!--                !-- Interpolation of chord as well as lift and drag       --!
!--                !-- coefficients from tabulated values                    --!
!--                !-----------------------------------------------------------!

!
!--                Interpolation of the chord_length from tabulated values to
!--                the current radial position:

                   IF (cur_r == 0.0_wp) THEN
                      chord(rseg) = 0.0_wp
                   ELSE IF (cur_r >= lrd(size(crd))) THEN
                      chord(rseg) = (crd(size(crd)) + ard(size(crd)-1)) / 2.0_wp
                   ELSE 
                      chord(rseg) = ( crd(lct(1)) *                            &
                            ( ( lrd(lct(1)+1) - cur_r ) /                      &
                              ( lrd(lct(1)+1) - lrd(lct(1)) ) ) ) +            &
                            ( crd(lct(1)+1) *                                  &
                            ( ( cur_r-lrd(lct(1)) ) /                          &
                              ( lrd(lct(1)+1) - lrd(lct(1)) ) ) )
                   ENDIF

!
!--                Determine index of current angle of attack, needed for
!--                finding the appropriate interpolated values of the lift and
!--                drag coefficients (-180.0 degrees = 0, +180.0 degrees = 36000,
!--                so one index every 0.01 degrees):
                   iialpha = CEILING( ( alpha_attack(rseg) + 180.0_wp )        &
                                      * ( 1.0_wp / accu_cl_cd_tab ) )
!
!--                Determine index of current radial position, needed for
!--                finding the appropriate interpolated values of the lift and
!--                drag coefficients (one index every 0.1 m):
                   iir = CEILING( cur_r * 10.0_wp )
!
!--                Read in interpolated values of the lift and drag coefficients
!--                for the current radial position and angle of attack:
                   turb_cl(rseg) = turb_cl_tab(iialpha,iir)
                   turb_cd(rseg) = turb_cd_tab(iialpha,iir)

!
!--                Final transformation back from angles given in degree to
!--                angles given in radian:
                   alpha_attack(rseg) = alpha_attack(rseg) *                   &
                                        ( (2.0_wp*pi) / 360.0_wp )

!
!--                !-----------------------------------------------------!
!--                !-- Calculation of the forces                       --!
!--                !-----------------------------------------------------!

!
!--                Calculate the pre_factor for the thrust and torque forces:

                   pre_factor = 0.5_wp * (vrel(rseg)**2) * 3.0_wp *  &
                                chord(rseg) * delta_r(inot) / nsegs(ring,inot)

!
!--                Calculate the thrust force (x-component of the total force)
!--                for each ring segment:
                   thrust_seg(rseg) = pre_factor *                             &
                                      ( turb_cl(rseg) * COS(phi_rel(rseg)) +   &
                                        turb_cd(rseg) * SIN(phi_rel(rseg)) )

!
!--                Determination of the second of the additional forces acting
!--                on the flow in the azimuthal direction: force vector as basis
!--                for torque (torque itself would be the vector product of the
!--                radius vector and the force vector):
                   torque_seg = pre_factor *                                   &
                                ( turb_cl(rseg) * SIN(phi_rel(rseg)) -         &
                                  turb_cd(rseg) * COS(phi_rel(rseg)) )
!
!--                Decomposition of the force vector into two parts:
!--                One acting along the y-direction and one acting along the
!--                z-direction of the rotor coordinate system:

                   torque_seg_y(rseg) = -torque_seg * sin_rot
                   torque_seg_z(rseg) =  torque_seg * cos_rot

!
!--                Add the segment thrust to the thrust of the whole rotor
                   thrust_rotor(inot) = thrust_rotor(inot) +                   &
                                        thrust_seg(rseg)                   
                   

                   torque_total(inot) = torque_total(inot) + (torque_seg * cur_r)

                ENDDO   !-- end of loop over ring segments

!
!--             Restore the forces into arrays containing all the segments of
!--             each ring:
                thrust_ring(ring,:)   = thrust_seg(:)
                torque_ring_y(ring,:) = torque_seg_y(:)
                torque_ring_z(ring,:) = torque_seg_z(:)


             ENDDO   !-- end of loop over rings


             CALL cpu_log( log_point_s(62), 'wtm_controller', 'start' )

             
             IF ( speed_control )  THEN
!
!--             Calculation of the current generator speed for rotor speed control
             
!                                     
!--             The acceleration of the rotor speed is calculated from 
!--             the force balance of the accelerating torque
!--             and the torque of the rotating rotor and generator
                om_rate = ( torque_total(inot) * air_dens * gear_eff -         &
                            gear_ratio * torque_gen_old(inot) ) /              &
                          ( inertia_rot +                                      & 
                            gear_ratio * gear_ratio * inertia_gen ) * dt_3d

!
!--             The generator speed is given by the product of gear gear_ratio
!--             and rotor speed
                omega_gen(inot) = gear_ratio * ( omega_rot(inot) + om_rate )     
             
             ENDIF
             
             IF ( pitch_control )  THEN

!
!--             If the current generator speed is above rated, the pitch is not
!--             saturated and the change from the last time step is within the 
!--             maximum pitch rate, then the pitch loop is repeated with a pitch
!--             gain
                IF ( (  omega_gen(inot)  > rated_genspeed   )  .AND.           &
                     ( pitch_add(inot) < 25.0_wp ) .AND.                       &
                     ( pitch_add(inot) < pitch_add_old(inot) +                 & 
                       pitch_rate * dt_3d  ) ) THEN 
                   pitch_sw = .TRUE.
!
!--                Go back to beginning of pit_loop                   
                   CYCLE pit_loop
                ENDIF
                
!
!--             The current pitch is saved for the next time step
                pitch_add_old(inot) = pitch_add(inot)
                pitch_sw = .FALSE.
             ENDIF
             EXIT pit_loop             
          ENDDO pit_loop ! Recursive pitch control loop


!
!--          Call the rotor speed controller
             
             IF ( speed_control )  THEN
!
!--             Find processor at i_hub, j_hub             
                IF ( ( nxl <= i_hub(inot) )  .AND.  ( nxr >= i_hub(inot) ) )   &
                   THEN
                   IF ( ( nys <= j_hub(inot) )  .AND.  ( nyn >= j_hub(inot) ) )&
                      THEN
                      CALL wtm_speed_control( inot )
                   ENDIF
                ENDIF
                               
             ENDIF


             CALL cpu_log( log_point_s(62), 'wtm_controller', 'stop' )

             CALL cpu_log( log_point_s(63), 'wtm_smearing', 'start' )


!--          !-----------------------------------------------------------------!
!--          !--                  Regularization kernel                      --!
!--          !-- Smearing of the forces and interpolation to cartesian grid  --!
!--          !-----------------------------------------------------------------!
!
!--          The aerodynamic blade forces need to be distributed smoothly on
!--          several mesh points in order to avoid singular behaviour
!
!--          Summation over sum of weighted forces. The weighting factor
!--          (calculated in user_init) includes information on the distance
!--          between the center of the grid cell and the rotor segment under
!--          consideration
!
!--          To save computing time, apply smearing only for the relevant part
!--          of the model domain:
!
!--
!--          Calculation of the boundaries:
             i_smear(inot) = CEILING( ( rr(inot) * ABS( roty(inot,1) ) +       &
                                        eps_min ) / dx )
             j_smear(inot) = CEILING( ( rr(inot) * ABS( roty(inot,2) ) +       &
                                        eps_min ) / dy )

             DO i = MAX( nxl, i_hub(inot) - i_smear(inot) ),                   &
                    MIN( nxr, i_hub(inot) + i_smear(inot) )
                DO j = MAX( nys, j_hub(inot) - j_smear(inot) ),                &
                        MIN( nyn, j_hub(inot) + j_smear(inot) )
                   DO k = MAX( nzb_u_inner(j,i)+1, k_hub(inot) - k_smear(inot) ), &
                                k_hub(inot) + k_smear(inot)
                      DO ring = 1, nrings(inot)
                         DO rseg = 1, nsegs(ring,inot)
!
!--                         Determine the square of the distance between the
!--                         current grid point and each rotor area segment:
                            dist_u_3d = ( i * dx               - rbx(ring,rseg) )**2 + &
                                        ( j * dy + 0.5_wp * dy - rby(ring,rseg) )**2 + &
                                        ( k * dz - 0.5_wp * dz - rbz(ring,rseg) )**2
                            dist_v_3d = ( i * dx + 0.5_wp * dx - rbx(ring,rseg) )**2 + &
                                        ( j * dy               - rby(ring,rseg) )**2 + &
                                        ( k * dz - 0.5_wp * dz - rbz(ring,rseg) )**2
                            dist_w_3d = ( i * dx + 0.5_wp * dx - rbx(ring,rseg) )**2 + &
                                        ( j * dy + 0.5_wp * dy - rby(ring,rseg) )**2 + &
                                        ( k * dz               - rbz(ring,rseg) )**2

!
!--                         3D-smearing of the forces with a polynomial function
!--                         (much faster than the old Gaussian function), using
!--                         some parameters that have been calculated in user_init.
!--                         The function is only similar to Gaussian function for
!--                         squared distances <= eps_min2:
                            IF ( dist_u_3d <= eps_min2 ) THEN
                            thrust(k,j,i) = thrust(k,j,i) +                    &
                                            thrust_ring(ring,rseg) *           &
                                            ( ( pol_a * dist_u_3d - pol_b ) *  & 
                                             dist_u_3d + 1.0_wp ) * eps_factor
                            ENDIF
                            IF ( dist_v_3d <= eps_min2 ) THEN
                            torque_y(k,j,i) = torque_y(k,j,i) +                &
                                              torque_ring_y(ring,rseg) *       &
                                              ( ( pol_a * dist_v_3d - pol_b ) *&
                                               dist_v_3d + 1.0_wp ) * eps_factor
                            ENDIF
                            IF ( dist_w_3d <= eps_min2 ) THEN
                            torque_z(k,j,i) = torque_z(k,j,i) +                &
                                              torque_ring_z(ring,rseg) *       &
                                              ( ( pol_a * dist_w_3d - pol_b ) *&
                                               dist_w_3d + 1.0_wp ) * eps_factor
                            ENDIF

                         ENDDO  ! End of loop over rseg
                      ENDDO     ! End of loop over ring
              
!
!--                   Rotation of force components:
                      rot_tend_x(k,j,i) = rot_tend_x(k,j,i) +                  &
                                      thrust(k,j,i)*rotx(inot,1) +             &
                                      torque_y(k,j,i)*roty(inot,1) +           &
                                      torque_z(k,j,i)*rotz(inot,1)
                                
                      rot_tend_y(k,j,i) = rot_tend_y(k,j,i) +                  &
                                      thrust(k,j,i)*rotx(inot,2) +             &
                                      torque_y(k,j,i)*roty(inot,2) +           &
                                      torque_z(k,j,i)*rotz(inot,2)
                                
                      rot_tend_z(k,j,i) = rot_tend_z(k,j,i) +                  &
                                      thrust(k,j,i)*rotx(inot,3) +             &
                                      torque_y(k,j,i)*roty(inot,3) +           &
                                      torque_z(k,j,i)*rotz(inot,3)                               

                   ENDDO        ! End of loop over k
                ENDDO           ! End of loop over j
             ENDDO              ! End of loop over i

             CALL cpu_log( log_point_s(63), 'wtm_smearing', 'stop' )          
                   
          ENDDO                  !-- end of loop over turbines

                
          IF ( yaw_control )  THEN
          
          ! TBD: MOVE TO WTM_YAW_CONTROL
!
!--          Allocate arrays for yaw control at first call
!--          Can't be allocated before dt_3d is set
             IF ( start_up )  THEN
                WDLON = NINT( 30.0_wp / dt_3d )  ! 30s running mean array
                ALLOCATE( wd30(1:nturbines,1:WDLON) )
                wd30 = 999.0_wp                  ! Set to dummy value
                ALLOCATE( wd30_l(1:WDLON) )
                
                WDSHO = NINT( 2.0_wp / dt_3d )   ! 2s running mean array
                ALLOCATE( wd2(1:nturbines,1:WDSHO) )
                wd2 = 999.0_wp                   ! Set to dummy value
                ALLOCATE( wd2_l(1:WDSHO) )
                start_up = .FALSE.
             ENDIF          

!
!--          Calculate the inflow wind speed 
!--
!--          Loop over number of turbines:
             DO inot = 1, nturbines
!
!--             Find processor at i_hub, j_hub             
                IF ( ( nxl <= i_hub(inot) )  .AND.  ( nxr >= i_hub(inot) ) )   &
                   THEN
                   IF ( ( nys <= j_hub(inot) )  .AND.  ( nyn >= j_hub(inot) ) )&
                      THEN

                      u_inflow_l(inot) = u(k_hub(inot),j_hub(inot),i_hub(inot))

                      wdir_l(inot) = -1.0_wp * ATAN2(                          &
                         0.5_wp * ( v(k_hub(inot),j_hub(inot),i_hub(inot)+1) + &
                                    v(k_hub(inot),j_hub(inot),i_hub(inot)) ) , &
                         0.5_wp * ( u(k_hub(inot),j_hub(inot)+1,i_hub(inot)) + &
                                    u(k_hub(inot),j_hub(inot),i_hub(inot)) ) )

                      CALL wtm_yawcontrol( inot )

                      phi_yaw_l(inot) = phi_yaw(inot)

                   ENDIF
                ENDIF
                   
             ENDDO                                 !-- end of loop over turbines

!
!--          Transfer of information to the other cpus
#if defined( __parallel )         
             CALL MPI_ALLREDUCE( u_inflow_l, u_inflow, nturbines, MPI_REAL,    &
                                 MPI_SUM, comm2d, ierr )
             CALL MPI_ALLREDUCE( wdir_l, wdir, nturbines, MPI_REAL, MPI_SUM,   &
                                 comm2d, ierr )
             CALL MPI_ALLREDUCE( phi_yaw_l, phi_yaw, nturbines, MPI_REAL,      &
                                 MPI_SUM, comm2d, ierr )
#else
             u_inflow = u_inflow_l
             wdir     = wdir_l
             phi_yaw  = phi_yaw_l
#endif
             DO inot = 1, nturbines
!             
!--             Update rotor orientation                
                CALL wtm_rotate_rotor( inot )

             ENDDO ! End of loop over turbines
                           
          END IF 
          
          IF ( speed_control )  THEN
!
!--          Transfer of information to the other cpus
!              CALL MPI_ALLREDUCE( omega_gen, omega_gen_old, nturbines,        &
!                                  MPI_REAL,MPI_SUM, comm2d, ierr )
#if defined( __parallel )    
             CALL MPI_ALLREDUCE( torque_gen, torque_gen_old, nturbines,        &
                                 MPI_REAL, MPI_SUM, comm2d, ierr )
             CALL MPI_ALLREDUCE( omega_rot_l, omega_rot, nturbines,            &
                                 MPI_REAL, MPI_SUM, comm2d, ierr )
             CALL MPI_ALLREDUCE( omega_gen_f, omega_gen_f_old, nturbines,      &
                                 MPI_REAL, MPI_SUM, comm2d, ierr )
#else
             torque_gen_old  = torque_gen
             omega_rot       = omega_rot_l
             omega_gen_f_old = omega_gen_f
#endif
            
          ENDIF

          DO inot = 1, nturbines

             IF ( myid == 0 ) THEN
                IF ( openfile_turb_mod(400+inot)%opened )  THEN
                   WRITE ( 400+inot, 106 ) simulated_time, omega_rot(inot),    &
                             omega_gen(inot), torque_gen_old(inot),            &
                             torque_total(inot), pitch_add(inot),              &
                             torque_gen_old(inot)*omega_gen(inot)*gen_eff,     &
                             torque_total(inot)*omega_rot(inot)*air_dens,      &
                             thrust_rotor(inot),                               & 
                             wdir(inot)*180.0_wp/pi,                           &
                             (phi_yaw(inot))*180.0_wp/pi                   
                             
                ELSE

                   WRITE ( turbine_id,'(I2.2)')  inot
                   OPEN ( 400+inot, FILE=( 'TURBINE_PARAMETERS'//turbine_id ), &
                                            FORM='FORMATTED' )
                   WRITE ( 400+inot, 105 ) inot
                   WRITE ( 400+inot, 106 ) simulated_time, omega_rot(inot),    &
                             omega_gen(inot), torque_gen_old(inot),            &
                             torque_total(inot), pitch_add(inot),              &
                             torque_gen_old(inot)*omega_gen(inot)*gen_eff,     &
                             torque_total(inot)*omega_rot(inot)*air_dens,      &
                             thrust_rotor(inot),                               & 
                             wdir(inot)*180.0_wp/pi,                           &                   
                             (phi_yaw(inot))*180.0_wp/pi
                ENDIF
             ENDIF

!--          Set open flag
             openfile_turb_mod(400+inot)%opened = .TRUE.
          ENDDO                                    !-- end of loop over turbines

       CALL cpu_log( log_point_s(61), 'wtm_forces', 'stop' )
       
!
!--    Formats
       105 FORMAT ('Turbine control data for turbine ',I2,1X,':'/ &
              &'----------------------------------------'/ &
              &'   Time   RSpeed  GSpeed  ', &
               'GenTorque  AeroTorque  Pitch  Power(Gen)  Power(Rot)  ',       &
               'RotThrust  WDirection  YawOrient')

       106 FORMAT (F9.3,2X,F7.3,2X,F7.2,2X,F9.1,3X,F9.1,1X,F6.2,2X,F10.1,2X,   &
                   F10.1,1X,F9.1,2X,F7.2,1X,F7.2)


    END SUBROUTINE wtm_forces_admr

    
!------------------------------------------------------------------------------!
! Description:
! ------------
!> Yaw controller for the wind turbine model
!------------------------------------------------------------------------------!
    SUBROUTINE wtm_yawcontrol( inot )
    
       USE constants
       USE kinds
                
       IMPLICIT NONE
      
       INTEGER(iwp)             :: inot
       INTEGER(iwp)             :: i_wd_30
       REAL(wp)                 :: missal

       i_wd_30 = 0_iwp

!  
!--    The yaw controller computes a 30s running mean of the wind direction.
!--    If the difference between turbine alignment and wind direction exceeds
!--    5°, the turbine is yawed. The mechanism stops as soon as the 2s-running
!--    mean of the missalignment is smaller than 0.5°.
!--    Attention: If the timestep during the simulation changes significantly
!--    the lengths of the running means change and it does not correspond to
!--    30s/2s anymore.
!--    ! Needs to be modified for these situations !
!--    For wind from the east, the averaging of the wind direction could cause
!--    problems and the yaw controller is probably flawed. -> Routine for
!--    averaging needs to be improved!
!
!--    Check if turbine is not yawing
       IF ( .NOT. doyaw(inot) )  THEN
!
!--       Write current wind direction into array
          wd30_l    = wd30(inot,:)
          wd30_l    = CSHIFT( wd30_l, SHIFT=-1 )
          wd30_l(1) = wdir(inot)
!
!--       Check if array is full ( no more dummies )
          IF ( .NOT. ANY( wd30_l == 999.) ) THEN 

             missal = SUM( wd30_l ) / SIZE( wd30_l ) - phi_yaw(inot)
!
!--          Check if turbine is missaligned by more than max_miss
             IF ( ABS( missal ) > max_miss )  THEN
!
!--             Check in which direction to yaw          
                yawdir(inot) = SIGN( 1.0_wp, missal )
!
!--             Start yawing of turbine
                phi_yaw(inot) = phi_yaw(inot) + yawdir(inot) * yaw_speed * dt_3d
                doyaw(inot) = .TRUE.
                wd30_l = 999.  ! fill with dummies again
             ENDIF
          ENDIF
         
          wd30(inot,:) = wd30_l

!      
!--    If turbine is already yawing:
!--    Initialize 2 s running mean and yaw until the missalignment is smaller
!--    than min_miss

       ELSE
!
!--       Initialize 2 s running mean
          wd2_l = wd2(inot,:)
          wd2_l = CSHIFT( wd2_l, SHIFT = -1 )
          wd2_l(1) = wdir(inot)
!      
!--       Check if array is full ( no more dummies )
          IF ( .NOT. ANY( wd2_l == 999.0_wp ) ) THEN
!
!--          Calculate missalignment of turbine        
             missal = SUM( wd2_l - phi_yaw(inot) ) / SIZE( wd2_l )
!
!--          Check if missalignment is still larger than 0.5 degree and if the
!--          yaw direction is still right
             IF ( ( ABS( missal ) > min_miss )  .AND.                          &
                  ( yawdir(inot) == SIGN( 1.0_wp, missal ) ) )  THEN
!
!--             Continue yawing        
                phi_yaw(inot) = phi_yaw(inot) + yawdir(inot) * yaw_speed * dt_3d
             ELSE
!
!--             Stop yawing        
                doyaw(inot) = .FALSE.
                wd2_l = 999.0_wp ! fill with dummies again
             ENDIF
          ELSE
!
!--          Continue yawing
             phi_yaw(inot) = phi_yaw(inot) + yawdir(inot) * yaw_speed * dt_3d
          ENDIF
      
          wd2(inot,:) = wd2_l
            
       ENDIF
     
    END SUBROUTINE wtm_yawcontrol 


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialization of the speed control
!------------------------------------------------------------------------------!
    SUBROUTINE wtm_init_speed_control


       IMPLICIT NONE

!
!--    If speed control is set, remaining variables and control_parameters for 
!--    the control algorithm are calculated
!
!--    Calculate slope constant for region 15
       slope15   = ( slope2 * min_reg2 * min_reg2 ) / ( min_reg2 - min_reg15 )
!
!--    Calculate upper limit of slipage region
       vs_sysp   = rated_genspeed / 1.1_wp
!
!--    Calculate slope of slipage region
       slope25   = ( rated_power / rated_genspeed ) /                          &
                   ( rated_genspeed - vs_sysp )
!
!--    Calculate lower limit of slipage region
       min_reg25 = ( slope25 - SQRT( slope25 * ( slope25 - 4.0_wp *            &
                                                 slope2 * vs_sysp ) ) ) /      &
                   ( 2.0_wp * slope2 )
!
!--    Frequency for the simple low pass filter
       Fcorner   = 0.25_wp
! 
!--    At the first timestep the torque is set to its maximum to prevent
!--    an overspeeding of the rotor
       torque_gen_old(:) = max_torque_gen  
     
    END SUBROUTINE wtm_init_speed_control


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Simple controller for the regulation of the rotor speed
!------------------------------------------------------------------------------!
    SUBROUTINE wtm_speed_control( inot )


       IMPLICIT NONE

       INTEGER(iwp)             :: inot
       
         

!
!--    The controller is based on the fortran script from Jonkman 
!--    et al. 2009 "Definition of a 5 MW Reference Wind Turbine for
!--    offshore system developement"

!
!--    The generator speed is filtered by a low pass filter 
!--    for the control of the generator torque       
       lp_coeff = EXP( -2.0_wp * 3.14_wp * dt_3d * Fcorner )
       omega_gen_f(inot) = ( 1.0_wp - lp_coeff ) * omega_gen(inot) + lp_coeff *&
                           omega_gen_f_old(inot)

       IF ( omega_gen_f(inot) <= min_reg15 )  THEN
!                        
!--       Region 1: Generator torque is set to zero to accelerate the rotor:
          torque_gen(inot) = 0
       
       ELSEIF ( omega_gen_f(inot) <= min_reg2 )  THEN
!                        
!--       Region 1.5: Generator torque is increasing linearly with rotor speed:
          torque_gen(inot) = slope15 * ( omega_gen_f(inot) - min_reg15 )
                         
       ELSEIF ( omega_gen_f(inot) <= min_reg25 )  THEN
!
!--       Region 2: Generator torque is increased by the square of the generator
!--                 speed to keep the TSR optimal:
          torque_gen(inot) = slope2 * omega_gen_f(inot) * omega_gen_f(inot)
       
       ELSEIF ( omega_gen_f(inot) < rated_genspeed )  THEN
!                        
!--       Region 2.5: Slipage region between 2 and 3:
          torque_gen(inot) = slope25 * ( omega_gen_f(inot) - vs_sysp )
       
       ELSE
!                        
!--       Region 3: Generator torque is antiproportional to the rotor speed to
!--                 keep the power constant:
          torque_gen(inot) = rated_power / omega_gen_f(inot)
       
       ENDIF
!                        
!--    Calculate torque rate and confine with a max
       trq_rate = ( torque_gen(inot) - torque_gen_old(inot) ) / dt_3d
       trq_rate = MIN( MAX( trq_rate, -1.0_wp * max_trq_rate ), max_trq_rate )
!                        
!--    Calculate new gen torque and confine with max torque                         
       torque_gen(inot) = torque_gen_old(inot) + trq_rate * dt_3d
       torque_gen(inot) = MIN( torque_gen(inot), max_torque_gen )                                             
!
!--    Overwrite values for next timestep                        
       omega_rot_l(inot) = omega_gen(inot) / gear_ratio

    
    END SUBROUTINE wtm_speed_control    


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Application of the additional forces generated by the wind turbine on the 
!> flow components (tendency terms)
!> Call for all grid points
!------------------------------------------------------------------------------!
    SUBROUTINE wtm_tendencies( component )

    
       IMPLICIT NONE

       INTEGER(iwp) ::  component   !< prognostic variable (u,v,w)
       INTEGER(iwp) ::  i           !< running index
       INTEGER(iwp) ::  j           !< running index
       INTEGER(iwp) ::  k           !< running index


       SELECT CASE ( component )

       CASE ( 1 )
!
!--       Apply the x-component of the force to the u-component of the flow:
          IF ( simulated_time >= time_turbine_on )  THEN
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   DO  k = nzb_u_inner(j,i)+1, k_hub(1) + k_smear(1)
!
!--                   Calculate the thrust generated by the nacelle and the tower
                      tend_nac_x = 0.5_wp * nac_cd_surf(k,j,i) *               &
                                         SIGN( u(k,j,i)**2 , u(k,j,i) )     
                      tend_tow_x   = 0.5_wp * tow_cd_surf(k,j,i) *             &
                                         SIGN( u(k,j,i)**2 , u(k,j,i) ) 
                                                   
                      tend(k,j,i) = tend(k,j,i) - rot_tend_x(k,j,i)            &
                                  - tend_nac_x - tend_tow_x
                   ENDDO
                ENDDO
             ENDDO
          ENDIF

       CASE ( 2 )
!
!--       Apply the y-component of the force to the v-component of the flow:
          IF ( simulated_time >= time_turbine_on )  THEN
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   DO  k = nzb_v_inner(j,i)+1, k_hub(1) + k_smear(1)
                      tend_nac_y = 0.5_wp * nac_cd_surf(k,j,i) *               &
                                         SIGN( v(k,j,i)**2 , v(k,j,i) )     
                      tend_tow_y   = 0.5_wp * tow_cd_surf(k,j,i) *             &
                                         SIGN( v(k,j,i)**2 , v(k,j,i) )                      
                      tend(k,j,i) = tend(k,j,i) - rot_tend_y(k,j,i)            &
                                  - tend_nac_y - tend_tow_y
                   ENDDO
                ENDDO
             ENDDO
          ENDIF

       CASE ( 3 )
!
!--       Apply the z-component of the force to the w-component of the flow:
          IF ( simulated_time >= time_turbine_on )  THEN
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   DO  k = nzb_w_inner(j,i)+1,  k_hub(1) + k_smear(1)
                      tend(k,j,i) = tend(k,j,i) - rot_tend_z(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF


       CASE DEFAULT

          WRITE( message_string, * ) 'unknown prognostic variable: ', component
          CALL message( 'wtm_tendencies', 'PA04??', 1, 2, 0, 6, 0 ) 

       END SELECT


    END SUBROUTINE wtm_tendencies


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Application of the additional forces generated by the wind turbine on the 
!> flow components (tendency terms)
!> Call for grid point i,j
!------------------------------------------------------------------------------!
    SUBROUTINE wtm_tendencies_ij( i, j, component )


       IMPLICIT NONE

       INTEGER(iwp) ::  component   !< prognostic variable (u,v,w)
       INTEGER(iwp) ::  i           !< running index
       INTEGER(iwp) ::  j           !< running index
       INTEGER(iwp) ::  k           !< running index

       SELECT CASE ( component )

       CASE ( 1 )
!
!--       Apply the x-component of the force to the u-component of the flow:
          IF ( simulated_time >= time_turbine_on )  THEN

             DO  k = nzb_u_inner(j,i)+1,  k_hub(1) + k_smear(1)
!
!--             Calculate the thrust generated by the nacelle and the tower 
                tend_nac_x = 0.5_wp * nac_cd_surf(k,j,i) *                     &
                                   SIGN( u(k,j,i)**2 , u(k,j,i) )     
                tend_tow_x   = 0.5_wp * tow_cd_surf(k,j,i) *                   &
                                   SIGN( u(k,j,i)**2 , u(k,j,i) ) 
                tend(k,j,i) = tend(k,j,i) - rot_tend_x(k,j,i)                  &
                            - tend_nac_x - tend_tow_x
             ENDDO
          ENDIF

       CASE ( 2 )
!
!--       Apply the y-component of the force to the v-component of the flow:
          IF ( simulated_time >= time_turbine_on )  THEN
             DO  k = nzb_v_inner(j,i)+1,  k_hub(1) + k_smear(1)
                tend_nac_y = 0.5_wp * nac_cd_surf(k,j,i) *                     &
                                   SIGN( v(k,j,i)**2 , v(k,j,i) )     
                tend_tow_y   = 0.5_wp * tow_cd_surf(k,j,i) *                   &
                                   SIGN( v(k,j,i)**2 , v(k,j,i) )                      
                tend(k,j,i) = tend(k,j,i) - rot_tend_y(k,j,i)                  &
                            - tend_nac_y - tend_tow_y
             ENDDO
          ENDIF

       CASE ( 3 )
!
!--       Apply the z-component of the force to the w-component of the flow:
          IF ( simulated_time >= time_turbine_on )  THEN
             DO  k = nzb_w_inner(j,i)+1,  k_hub(1) + k_smear(1)
                tend(k,j,i) = tend(k,j,i) - rot_tend_z(k,j,i)
             ENDDO
          ENDIF


       CASE DEFAULT

          WRITE( message_string, * ) 'unknown prognostic variable: ', component
          CALL message( 'wtm_tendencies', 'PA04??', 1, 2, 0, 6, 0 ) 

       END SELECT


    END SUBROUTINE wtm_tendencies_ij

 END MODULE wind_turbine_model_mod
