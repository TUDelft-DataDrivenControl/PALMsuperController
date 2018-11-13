 MODULE wind_turbine_model_mod

    USE arrays_3d,                                                             &
        ONLY:  dd2zu, tend, u, v, w, zu, zw

    USE constants

    USE control_parameters,                                                    &
        ONLY:  dt_3d, dz, message_string, simulated_time,                      &
               current_timestep_number

    USE cpulog,                                                                &
        ONLY:  cpu_log, log_point_s

    USE grid_variables,                                                        &
        ONLY:  ddx, dx, ddy, dy

    USE indices,                                                               &
        ONLY:  nbgp, nx, nxl, nxlg, nxr, nxrg, ny, nyn, nyng, nys, nysg, nz,   &
               nzb, nzb_u_inner, nzb_v_inner, nzb_w_inner, nzt

    USE kinds

    USE pegrid


    IMPLICIT NONE

    PRIVATE

    LOGICAL ::  wind_turbine=.FALSE.   !< switch for use of wind turbine model
!
!-- String for specfic user and controller
    CHARACTER(LEN=100) :: str = '/home/sboersma3/palm/current_version/USER_CODE/2turb/matlab/'

!
!-- Variables specified in the namelist wind_turbine_par

    CHARACTER(LEN=4) ::  turb_mod = 'admr' !< type of turbine model

    INTEGER(iwp) ::  nairfoils = 8   !< number of airfoils of the used turbine model (for ADM-R and ALM)
    INTEGER(iwp) ::  nturbines = 1   !< number of turbines

    LOGICAL ::  pitch_control = .FALSE.   !< switch for use of pitch controller
    LOGICAL ::  speed_control = .FALSE.   !< switch for use of speed controller
    LOGICAL ::  yaw_control   = .FALSE.   !< switch for use of yaw controller

    REAL(wp) ::  segment_length  = 1.0_wp          !< length of the segments, the rotor area is divided into
                                                   !< (in tangential direction, as factor of MIN(dx,dy,dz))
    REAL(wp) ::  segment_width   = 0.5_wp          !< width of the segments, the rotor area is divided into
                                                   !< (in radial direction, as factor of MIN(dx,dy,dz))
    REAL(wp) ::  time_turbine_on = 0.0_wp          !< time at which turbines are started
    REAL(wp) ::  tilt            = 0.0_wp          !< vertical tilt of the rotor [degree] ( positive = backwards )

    REAL(wp), DIMENSION(1:100) ::  dtow             = 0.0_wp  !< tower diameter [m]
    REAL(wp), DIMENSION(1:100) ::  omega_rot        = 0.0_wp  !< inital or constant rotor speed [rad/s]
    REAL(wp), DIMENSION(1:100) ::  phi_yaw          = 0.0_wp  !< yaw angle [degree] ( clockwise, 0 = facing west )
    REAL(wp), DIMENSION(1:100) ::  pitch_add        = 0.0_wp  !< constant pitch angle
    REAL(wp), DIMENSION(1:100) ::  rcx        = 9999999.9_wp  !< position of hub in x-direction
    REAL(wp), DIMENSION(1:100) ::  rcy        = 9999999.9_wp  !< position of hub in y-direction
    REAL(wp), DIMENSION(1:100) ::  rcz        = 9999999.9_wp  !< position of hub in z-direction
    REAL(wp), DIMENSION(1:100) ::  rnac             = 0.0_wp  !< nacelle diameter [m]
    REAL(wp), DIMENSION(1:100) ::  rr              = 63.0_wp  !< rotor radius [m]
    REAL(wp), DIMENSION(1:100) ::  turb_cd_nacelle = 0.85_wp  !< drag coefficient for nacelle
    REAL(wp), DIMENSION(1:100) ::  turb_cd_tower    = 1.2_wp  !< drag coefficient for tower

!
!-- Variables specified in the namelist for speed controller
!-- Default values are from the NREL 5MW research turbine (Jonkman, 2008)

    REAL(wp) ::  rated_power    = 5296610.0_wp    !< rated turbine power [W]
    REAL(wp) ::  gear_ratio     = 97.0_wp         !< Gear ratio from rotor to generator
    REAL(wp) ::  inertia_rot    = 34784179.0_wp   !< Inertia of the rotor [kg/m2]
    REAL(wp) ::  inertia_gen    = 534.116_wp      !< Inertia of the generator [kg/m2]
    REAL(wp) ::  gen_eff        = 0.944_wp        !< Electric efficiency of the generator
    REAL(wp) ::  gear_eff       = 1.0_wp          !< Loss between rotor and generator
    REAL(wp) ::  air_dens       = 1.225_wp        !< Air density to convert to W [kg/m3]
    REAL(wp) ::  rated_genspeed = 121.6805_wp     !< Rated generator speed [rad/s]
    REAL(wp) ::  max_torque_gen = 47402.91_wp     !< Maximum of the generator torque [Nm]
    REAL(wp) ::  slope2         = 2.332287_wp     !< Slope constant for region 2
    REAL(wp) ::  min_reg2       = 91.21091_wp     !< Lower generator speed boundary of region 2 [rad/s]
    REAL(wp) ::  min_reg15      = 70.16224_wp     !< Lower generator speed boundary of region 1.5 [rad/s]
    REAL(wp) ::  max_trq_rate   = 15000.0_wp      !< Max generator torque increase [Nm/s]
    REAL(wp) ::  pitch_rate     = 8.0_wp          !< Max pitch rate [degree/s]


!
!-- Variables specified in the namelist for yaw control

    REAL(wp) ::  yaw_speed = 0.005236_wp   !< speed of the yaw actuator [rad/s]
    REAL(wp) ::  max_miss = 0.08726_wp     !< maximum tolerated yaw missalignment [rad]
    REAL(wp) ::  min_miss = 0.008726_wp    !< minimum yaw missalignment for which the actuator stops [rad]

!
!-- Set flag for output files TURBINE_PARAMETERS
    TYPE file_status
       LOGICAL ::  opened, opened_before
    END TYPE file_status
    
    TYPE(file_status), DIMENSION(500) :: openfile_turb_mod =                   &
                                         file_status(.FALSE.,.FALSE.)

!
!-- Variables for initialization of the turbine model

    INTEGER(iwp) ::  inot         !< turbine loop index (turbine id)
    INTEGER(iwp) ::  nsegs_max    !< maximum number of segments (all turbines, required for allocation of arrays)
    INTEGER(iwp) ::  nrings_max   !< maximum number of rings (all turbines, required for allocation of arrays)
    INTEGER(iwp) ::  ring         !< ring loop index (ring number)
    INTEGER(iwp) ::  rr_int       !<
    INTEGER(iwp) ::  upper_end    !<

    INTEGER(iwp), DIMENSION(1) ::  lct   !<

    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  i_hub     !< index belonging to x-position of the turbine
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  i_smear   !< index defining the area for the smearing of the forces (x-direction)
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  j_hub     !< index belonging to y-position of the turbine
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  j_smear   !< index defining the area for the smearing of the forces (y-direction)
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  k_hub     !< index belonging to hub height
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  k_smear   !< index defining the area for the smearing of the forces (z-direction)
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  nrings    !< number of rings per turbine
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  nsegs_total !< total number of segments per turbine 

    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE ::  nsegs   !< number of segments per ring and turbine

!
!-  parameters for the smearing from the rotor to the cartesian grid   
    REAL(wp) ::  pol_a            !< parameter for the polynomial smearing fct
    REAL(wp) ::  pol_b            !< parameter for the polynomial smearing fct
    REAL(wp) ::  delta_t_factor   !< 
    REAL(wp) ::  eps_factor       !<  
    REAL(wp) ::  eps_min          !< 
    REAL(wp) ::  eps_min2         !< 
    REAL(wp) ::  sqrt_arg         !< 

!
!-- Variables for the calculation of lift and drag coefficients
    REAL(wp), DIMENSION(:), ALLOCATABLE  ::  ard     !< 
    REAL(wp), DIMENSION(:), ALLOCATABLE  ::  crd     !< 
    REAL(wp), DIMENSION(:), ALLOCATABLE  ::  delta_r !< radial segment length
    REAL(wp), DIMENSION(:), ALLOCATABLE  ::  lrd     !<
    
    REAL(wp) ::  accu_cl_cd_tab = 0.1_wp  !< Accuracy of the interpolation of 
                                          !< the lift and drag coeff [deg] 

    REAL(wp), DIMENSION(:,:), ALLOCATABLE :: turb_cd_tab   !< table of the blade drag coefficient
    REAL(wp), DIMENSION(:,:), ALLOCATABLE :: turb_cl_tab   !< table of the blade lift coefficient

    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  nac_cd_surf  !< 3d field of the tower drag coefficient 
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  tow_cd_surf  !< 3d field of the nacelle drag coefficient

!
!-- Variables for the calculation of the forces
     
    REAL(wp) ::  cur_r                       !< 
    REAL(wp) ::  delta_t                     !<  tangential segment length
    REAL(wp) ::  phi_rotor                   !< 
    REAL(wp) ::  pre_factor                  !<  
    REAL(wp) ::  torque_seg                  !< 
    REAL(wp) ::  u_int_l                     !< 
    REAL(wp) ::  u_int_u                     !< 
    REAL(wp) ::  u_rot                       !< 
    REAL(wp) ::  v_int_l                     !< 
    REAL(wp) ::  v_int_u                     !< 
    REAL(wp) ::  w_int_l                     !< 
    REAL(wp) ::  w_int_u                     !<
!
!-  Tendencies from the nacelle and tower thrust
    REAL(wp) ::  tend_nac_x = 0.0_wp  !< 
    REAL(wp) ::  tend_tow_x = 0.0_wp  !< 
    REAL(wp) ::  tend_nac_y = 0.0_wp  !< 
    REAL(wp) ::  tend_tow_y = 0.0_wp  !<

    REAL(wp), DIMENSION(:), ALLOCATABLE ::  alpha_attack !< 
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  chord        !< 
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  omega_gen    !< curr. generator speed
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  phi_rel      !< 
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  pitch_add_old!<
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  torque_total !<
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  thrust_rotor !<
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  turb_cl      !< 
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  turb_cd      !< 
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  vrel         !< 
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  vtheta       !< 

    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  rbx, rby, rbz     !< coordinates of the blade elements
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  rotx, roty, rotz  !< normal vectors to the rotor coordinates

!
!-  Fields for the interpolation of velocities on the rotor grid
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  u_int       !< 
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  u_int_1_l   !< 
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  v_int       !< 
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  v_int_1_l   !< 
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  w_int       !< 
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  w_int_1_l   !< 
    
!
!-  rotor tendencies on the segments 
    REAL(wp), DIMENSION(:), ALLOCATABLE :: thrust_seg   !< 
    REAL(wp), DIMENSION(:), ALLOCATABLE :: torque_seg_y !< 
    REAL(wp), DIMENSION(:), ALLOCATABLE :: torque_seg_z !<    

!
!-  rotor tendencies on the rings
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  thrust_ring       !< 
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  torque_ring_y     !< 
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  torque_ring_z     !<
    
!
!-  rotor tendencies on rotor grids for all turbines
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  thrust      !<
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  torque_y    !< 
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  torque_z    !< 

!
!-  rotor tendencies on coordinate grid
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  rot_tend_x  !< 
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  rot_tend_y  !< 
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  rot_tend_z  !< 
!    
!-  variables for the rotation of the rotor coordinates        
    REAL(wp), DIMENSION(1:100,1:3,1:3) ::  rot_coord_trans  !< matrix for rotation of rotor coordinates
    
    REAL(wp), DIMENSION(1:3) ::  rot_eigen_rad   !< 
    REAL(wp), DIMENSION(1:3) ::  rot_eigen_azi   !< 
    REAL(wp), DIMENSION(1:3) ::  rot_eigen_nor   !< 
    REAL(wp), DIMENSION(1:3) ::  re              !< 
    REAL(wp), DIMENSION(1:3) ::  rea             !< 
    REAL(wp), DIMENSION(1:3) ::  ren             !< 
    REAL(wp), DIMENSION(1:3) ::  rote            !< 
    REAL(wp), DIMENSION(1:3) ::  rota            !< 
    REAL(wp), DIMENSION(1:3) ::  rotn            !< 

!
!-- Fixed variables for the speed controller

    LOGICAL  ::  start_up = .TRUE.   !< 
    
    REAL(wp) ::  Fcorner             !< corner freq for the controller low pass filter
    REAL(wp) ::  min_reg25           !< min region 2.5
    REAL(wp) ::  om_rate             !< rotor speed change
    REAL(wp) ::  slope15             !< slope in region 1.5
    REAL(wp) ::  slope25             !< slope in region 2.5
    REAL(wp) ::  trq_rate            !< torque change
    REAL(wp) ::  vs_sysp             !< 
    REAL(wp) ::  lp_coeff            !< coeff for the controller low pass filter 
    
    REAL(wp), DIMENSION(:), ALLOCATABLE :: omega_gen_old   !< last gen. speed
    REAL(wp), DIMENSION(:), ALLOCATABLE :: omega_gen_f     !< filtered gen. sp
    REAL(wp), DIMENSION(:), ALLOCATABLE :: omega_gen_f_old !< last filt. gen. sp
    REAL(wp), DIMENSION(:), ALLOCATABLE :: torque_gen      !< generator torque
    REAL(wp), DIMENSION(:), ALLOCATABLE :: torque_gen_old  !< last gen. torque

    REAL(wp), DIMENSION(100) :: omega_rot_l = 0.0_wp !< local rot speed [rad/s]
!
!-- Fixed variables for the yaw controller

    REAL(wp), DIMENSION(:)  , ALLOCATABLE ::  yawdir           !< direction to yaw
    REAL(wp), DIMENSION(:)  , ALLOCATABLE ::  phi_yaw_l        !< local (cpu) yaw angle
    REAL(wp), DIMENSION(:)  , ALLOCATABLE ::  wd30_l           !< local (cpu) long running avg of the wd
    REAL(wp), DIMENSION(:)  , ALLOCATABLE ::  wd2_l            !< local (cpu) short running avg of the wd
    REAL(wp), DIMENSION(:)  , ALLOCATABLE ::  wdir             !< wind direction at hub
    REAL(wp), DIMENSION(:)  , ALLOCATABLE ::  u_inflow         !< wind speed at hub
    REAL(wp), DIMENSION(:)  , ALLOCATABLE ::  wdir_l           !< 
    REAL(wp), DIMENSION(:)  , ALLOCATABLE ::  u_inflow_l       !<
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  wd30             !<
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  wd2              !<
    LOGICAL,  DIMENSION(1:100)            ::  doyaw = .FALSE.  !<
    INTEGER(iwp)                          ::  WDLON            !<
    INTEGER(iwp)                          ::  WDSHO            !<

!
!-- ct-curve (adm)
    REAL(wp),DIMENSION(:), ALLOCATABLE ::  vct
    REAL(wp),DIMENSION(:), ALLOCATABLE ::  dct
!
!-- adm    
    REAL(wp),DIMENSION(:), ALLOCATABLE ::  aif
    REAL(wp),DIMENSION(:), ALLOCATABLE ::  sums_u
    REAL(wp),DIMENSION(:), ALLOCATABLE ::  u_rotor
    REAL(wp),DIMENSION(:), ALLOCATABLE ::  u_rotor_old
    REAL(wp),DIMENSION(:), ALLOCATABLE ::  turb_ct
    
    REAL(wp),DIMENSION(:,:,:), ALLOCATABLE ::  turb_area


    REAL(wp),DIMENSION(:), ALLOCATABLE ::  thrust_adm
    REAL(wp),DIMENSION(:), ALLOCATABLE ::  power_rotor 
    REAL(wp),DIMENSION(:), ALLOCATABLE ::  power_adm 
    
    SAVE


    INTERFACE wtm_parin
       MODULE PROCEDURE wtm_parin
    END INTERFACE wtm_parin
    
    INTERFACE wtm_check_parameters
       MODULE PROCEDURE wtm_check_parameters
    END INTERFACE wtm_check_parameters
        
    INTERFACE wtm_init_arrays
       MODULE PROCEDURE wtm_init_arrays
    END INTERFACE wtm_init_arrays

    INTERFACE wtm_init
       MODULE PROCEDURE wtm_init
    END INTERFACE wtm_init
    
    INTERFACE wtm_read_blade_tables
       MODULE PROCEDURE wtm_read_blade_tables
    END INTERFACE wtm_read_blade_tables
    
    !-- adm
    INTERFACE wtm_read_thrust_curve
       MODULE PROCEDURE wtm_read_thrust_curve
    END INTERFACE wtm_read_thrust_curve
           
    INTERFACE wtm_forces
       MODULE PROCEDURE wtm_forces
       MODULE PROCEDURE wtm_yawcontrol
    END INTERFACE wtm_forces
    
    !-- adm
    INTERFACE wtm_forces_adm
       MODULE PROCEDURE wtm_forces_adm
    END INTERFACE
    
    INTERFACE wtm_forces_admr
       MODULE PROCEDURE wtm_forces_admr
    END INTERFACE    
    
    INTERFACE wtm_rotate_rotor
       MODULE PROCEDURE wtm_rotate_rotor
    END INTERFACE wtm_rotate_rotor
    
    INTERFACE wtm_speed_control
       MODULE PROCEDURE wtm_init_speed_control
       MODULE PROCEDURE wtm_speed_control
    END INTERFACE wtm_speed_control
    
    !-- adm
    INTERFACE wtm_surface_rotor
       MODULE PROCEDURE wtm_surface_rotor
    END INTERFACE wtm_surface_rotor

    INTERFACE wtm_tendencies
       MODULE PROCEDURE wtm_tendencies
       MODULE PROCEDURE wtm_tendencies_ij
    END INTERFACE wtm_tendencies
    
    
    PUBLIC wtm_check_parameters, wtm_forces, wtm_init, wtm_init_arrays,        &
           wtm_parin, wtm_tendencies, wtm_tendencies_ij, wind_turbine 

 CONTAINS


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Parin for &wind_turbine_par for wind turbine model
!------------------------------------------------------------------------------!
    SUBROUTINE wtm_parin


       IMPLICIT NONE
       
       INTEGER(iwp) ::  ierrn       !<

       CHARACTER (LEN=80) ::  line  !< dummy string that contains the current line of the parameter file 

       NAMELIST /wind_turbine_par/   air_dens, dtow, gear_eff, gear_ratio,     &
                                  gen_eff, inertia_gen, inertia_rot, max_miss, &
                                  max_torque_gen, max_trq_rate, min_miss,      &
                                  min_reg15, min_reg2, nairfoils, nturbines,   &
                                  omega_rot, phi_yaw, pitch_add, pitch_control,&
                                  rated_genspeed, rated_power, rcx, rcy, rcz,  &
                                  rnac, rr, segment_length, segment_width,     &
                                  slope2, speed_control, tilt, time_turbine_on,&
                                  turb_cd_nacelle, turb_cd_tower,              &
                                  yaw_control, yaw_speed, turb_mod

!
!--    Try to find wind turbine model package
       REWIND ( 11 )
       line = ' '
       DO  WHILE ( INDEX( line, '&wind_turbine_par' ) == 0 )
          READ ( 11, '(A)', END=10 )  line
       ENDDO
       BACKSPACE ( 11 )

!
!--    Read user-defined namelist
       READ ( 11, wind_turbine_par, IOSTAT=ierrn )
       
       IF ( ierrn < 0 )  THEN
          message_string = 'errors in \$wind_turbine_par'
          CALL message( 'wtm_parin', 'PA0???', 1, 2, 0, 6, 0 )
       ENDIF
       
!
!--    Set flag that indicates that the wind turbine model is switched on
       wind_turbine = .TRUE.

 10    CONTINUE   ! TBD Change from continue, mit ierrn machen


    END SUBROUTINE wtm_parin

    SUBROUTINE wtm_check_parameters

    
       IMPLICIT NONE
       
       SELECT CASE ( turb_mod )
          CASE ( 'admr' )    
             IF ( ( .NOT.speed_control ) .AND. pitch_control )  THEN
                message_string = 'pitch_control = .TRUE. requires '//          &
                                'speed_control = .TRUE.'
                CALL message( 'wtm_check_parameters', 'PA0???', 1, 2, 0, 6, 0 )
             ENDIF
            
             IF ( ANY( omega_rot(1:nturbines) <= 0.0 ) )  THEN
                message_string = 'omega_rot <= 0.0, Please set omega_rot to '//&
                                'a value larger than zero'
                CALL message( 'wtm_check_parameters', 'PA0???', 1, 2, 0, 6, 0 )
             ENDIF
             
       END SELECT
       
       
       IF ( ANY( rcx(1:nturbines) == 9999999.9_wp ) .OR.                       &
            ANY( rcy(1:nturbines) == 9999999.9_wp ) .OR.                       &
            ANY( rcz(1:nturbines) == 9999999.9_wp ) )  THEN
          
          message_string = 'rcx, rcy, rcz '                                 // &
                           'have to be given for each turbine.'          
          CALL message( 'wtm_check_parameters', 'PA0???', 1, 2, 0, 6, 0 )          
          
       ENDIF

 
    END SUBROUTINE wtm_check_parameters 
    
                                        
!------------------------------------------------------------------------------!
! Description:
! ------------
!> Allocate wind turbine model arrays
!------------------------------------------------------------------------------!
    SUBROUTINE wtm_init_arrays


       IMPLICIT NONE

       REAL(wp) ::  delta_r_factor   !< 
       REAL(wp) ::  delta_r_init     !< 

!
!--    To be able to allocate arrays with dimension of rotor rings and segments,
!--    the maximum possible numbers of rings and segments have to be calculated:

       ALLOCATE( nrings(1:nturbines) )
       ALLOCATE( delta_r(1:nturbines) )

       nrings(:)  = 0
       delta_r(:) = 0.0_wp

!
!--    Thickness (radial) of each ring and length (tangential) of each segment:
       delta_r_factor = segment_width
       delta_t_factor = segment_length
       delta_r_init   = delta_r_factor * MIN( dx, dy, dz)
       delta_t        = delta_t_factor * MIN( dx, dy, dz)

       DO inot = 1, nturbines
!
!--       Determine number of rings:
          nrings(inot) = NINT( rr(inot) / delta_r_init )

          delta_r(inot) = rr(inot) / nrings(inot)

       ENDDO

       nrings_max = MAXVAL(nrings)

       ALLOCATE( nsegs(1:nrings_max,1:nturbines) )
       ALLOCATE( nsegs_total(1:nturbines) )

       nsegs(:,:)     = 0
       nsegs_total(:) = 0


       DO inot = 1, nturbines
          DO ring = 1, nrings(inot)
!
!--          Determine number of segments for each ring:
             nsegs(ring,inot) = MAX( 8, CEILING( delta_r_factor * pi *         &
                                                 ( 2.0_wp * ring - 1.0_wp ) /  &
                                                 delta_t_factor ) )
          ENDDO
!
!--       Total sum of all rotor segments:
          nsegs_total(inot) = SUM( nsegs(:,inot) )

       ENDDO

!
!--    Maximum number of segments per ring:
       nsegs_max = MAXVAL(nsegs)

!
!--    Allocate 1D arrays (dimension = number of turbines)
       ALLOCATE( i_hub(1:nturbines) )
       ALLOCATE( i_smear(1:nturbines) )
       ALLOCATE( j_hub(1:nturbines) )
       ALLOCATE( j_smear(1:nturbines) )
       ALLOCATE( k_hub(1:nturbines) )
       ALLOCATE( k_smear(1:nturbines) )
       ALLOCATE( torque_total(1:nturbines) )
       ALLOCATE( thrust_rotor(1:nturbines) )

!
!--    Allocation of the 1D arrays for speed pitch_control
       ALLOCATE( omega_gen(1:nturbines) )
       ALLOCATE( omega_gen_old(1:nturbines) )
       ALLOCATE( omega_gen_f(1:nturbines) )
       ALLOCATE( omega_gen_f_old(1:nturbines) )
       ALLOCATE( pitch_add_old(1:nturbines) )
       ALLOCATE( torque_gen(1:nturbines) )
       ALLOCATE( torque_gen_old(1:nturbines) )

!
!--    Allocation of the 1D arrays for yaw control
       ALLOCATE( yawdir(1:nturbines) )
       ALLOCATE( u_inflow(1:nturbines) )
       ALLOCATE( wdir(1:nturbines) )
       ALLOCATE( u_inflow_l(1:nturbines) )
       ALLOCATE( wdir_l(1:nturbines) )
       ALLOCATE( phi_yaw_l(1:nturbines) )
       
!
!--    Allocate 1D arrays (dimension = number of rotor segments)
       ALLOCATE( alpha_attack(1:nsegs_max) )
       ALLOCATE( chord(1:nsegs_max) )
       ALLOCATE( phi_rel(1:nsegs_max) )
       ALLOCATE( thrust_seg(1:nsegs_max) )
       ALLOCATE( torque_seg_y(1:nsegs_max) )
       ALLOCATE( torque_seg_z(1:nsegs_max) )
       ALLOCATE( turb_cd(1:nsegs_max) )
       ALLOCATE( turb_cl(1:nsegs_max) )
       ALLOCATE( vrel(1:nsegs_max) )
       ALLOCATE( vtheta(1:nsegs_max) )

!
!--    Allocate 2D arrays (dimension = number of rotor rings and segments)
       ALLOCATE( rbx(1:nrings_max,1:nsegs_max) )
       ALLOCATE( rby(1:nrings_max,1:nsegs_max) )
       ALLOCATE( rbz(1:nrings_max,1:nsegs_max) )
       ALLOCATE( thrust_ring(1:nrings_max,1:nsegs_max) )
       ALLOCATE( torque_ring_y(1:nrings_max,1:nsegs_max) )
       ALLOCATE( torque_ring_z(1:nrings_max,1:nsegs_max) )

!
!--    Allocate additional 2D arrays
       ALLOCATE( rotx(1:nturbines,1:3) )
       ALLOCATE( roty(1:nturbines,1:3) )
       ALLOCATE( rotz(1:nturbines,1:3) )

!
!--    Allocate 3D arrays (dimension = number of grid points)
       ALLOCATE( nac_cd_surf(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
       ALLOCATE( rot_tend_x(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
       ALLOCATE( rot_tend_y(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
       ALLOCATE( rot_tend_z(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
       ALLOCATE( thrust(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
       ALLOCATE( torque_y(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
       ALLOCATE( torque_z(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
       ALLOCATE( tow_cd_surf(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )

!
!--    Allocate additional 3D arrays
       ALLOCATE( u_int(1:nturbines,1:nrings_max,1:nsegs_max) )
       ALLOCATE( u_int_1_l(1:nturbines,1:nrings_max,1:nsegs_max) )
       ALLOCATE( v_int(1:nturbines,1:nrings_max,1:nsegs_max) )
       ALLOCATE( v_int_1_l(1:nturbines,1:nrings_max,1:nsegs_max) )
       ALLOCATE( w_int(1:nturbines,1:nrings_max,1:nsegs_max) )
       ALLOCATE( w_int_1_l(1:nturbines,1:nrings_max,1:nsegs_max) )

!
!--    All of the arrays are initialized with a value of zero:
       i_hub(:)                 = 0
       i_smear(:)               = 0
       j_hub(:)                 = 0
       j_smear(:)               = 0
       k_hub(:)                 = 0
       k_smear(:)               = 0
       
       torque_total(:)          = 0.0_wp
       thrust_rotor(:)          = 0.0_wp

       omega_gen(:)             = 0.0_wp
       omega_gen_old(:)         = 0.0_wp
       omega_gen_f(:)           = 0.0_wp
       omega_gen_f_old(:)       = 0.0_wp
       pitch_add_old(:)         = 0.0_wp
       torque_gen(:)            = 0.0_wp
       torque_gen_old(:)        = 0.0_wp
       
       yawdir(:)                = 0.0_wp
       wdir(:)                  = 0.0_wp
       u_inflow(:)              = 0.0_wp

!
!--    Allocate 1D arrays (dimension = number of rotor segments)
       alpha_attack(:)          = 0.0_wp
       chord(:)                 = 0.0_wp
       phi_rel(:)               = 0.0_wp
       thrust_seg(:)            = 0.0_wp
       torque_seg_y(:)          = 0.0_wp
       torque_seg_z(:)          = 0.0_wp
       turb_cd(:)               = 0.0_wp
       turb_cl(:)               = 0.0_wp
       vrel(:)                  = 0.0_wp
       vtheta(:)                = 0.0_wp

       rbx(:,:)                 = 0.0_wp
       rby(:,:)                 = 0.0_wp
       rbz(:,:)                 = 0.0_wp
       thrust_ring(:,:)         = 0.0_wp
       torque_ring_y(:,:)       = 0.0_wp
       torque_ring_z(:,:)       = 0.0_wp

       rotx(:,:)                = 0.0_wp
       roty(:,:)                = 0.0_wp
       rotz(:,:)                = 0.0_wp

       nac_cd_surf(:,:,:)       = 0.0_wp
       rot_tend_x(:,:,:)        = 0.0_wp
       rot_tend_y(:,:,:)        = 0.0_wp
       rot_tend_z(:,:,:)        = 0.0_wp
       thrust(:,:,:)            = 0.0_wp
       torque_y(:,:,:)          = 0.0_wp
       torque_z(:,:,:)          = 0.0_wp
       tow_cd_surf(:,:,:)       = 0.0_wp

       u_int(:,:,:)             = 0.0_wp
       u_int_1_l(:,:,:)         = 0.0_wp
       v_int(:,:,:)             = 0.0_wp
       v_int_1_l(:,:,:)         = 0.0_wp
       w_int(:,:,:)             = 0.0_wp
       w_int_1_l(:,:,:)         = 0.0_wp
       
       SELECT CASE ( turb_mod )
         CASE ( 'adm' )
                  
            ALLOCATE( aif(1:nturbines) )
            ALLOCATE( sums_u(1:nturbines) )
            ALLOCATE( u_rotor(1:nturbines) )
            ALLOCATE( u_rotor_old(1:nturbines) )
            ALLOCATE( turb_ct(1:nturbines) )
            
            ALLOCATE( turb_area(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )

            ALLOCATE( thrust_adm(1:nturbines) )
            ALLOCATE( power_rotor(1:nturbines) )
            ALLOCATE( power_adm(1:nturbines) )

            aif(:)                = 0.0_wp
            sums_u(:)             = 0.0_wp
            u_rotor(:)            = 0.0_wp
            u_rotor_old(:)        = 0.0_wp
            turb_ct(:)            = 0.0_wp
            
            turb_area(:,:,:)      = 0.0_wp

            thrust_adm(:)         = 0.0_wp

       END SELECT

    END SUBROUTINE wtm_init_arrays


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialization of the wind turbine model
!------------------------------------------------------------------------------!
    SUBROUTINE wtm_init

    
       IMPLICIT NONE

       INTEGER(iwp) ::  i  !< running index
       INTEGER(iwp) ::  j  !< running index
       INTEGER(iwp) ::  k  !< running index
       
!
!--    Help variables for the smearing function       
       REAL(wp) ::  eps_kernel       !<       
       
!
!--    Help variables for calculation of the tower drag       
       INTEGER(iwp) ::  tower_n      !<
       INTEGER(iwp) ::  tower_s      !<
!
!--    Help variables for the calulaction of the nacelle drag
       INTEGER(iwp) ::  i_ip         !<
       INTEGER(iwp) ::  i_ipg        !<
       
       REAL(wp) ::  yvalue               
       REAL(wp) ::  dy_int           !< 
       REAL(wp) ::  dz_int           !<
       
              
       INTEGER(iwp), DIMENSION(:), ALLOCATABLE :: index_nacb       !< 
       INTEGER(iwp), DIMENSION(:), ALLOCATABLE :: index_nacl       !< 
       INTEGER(iwp), DIMENSION(:), ALLOCATABLE :: index_nacr       !< 
       INTEGER(iwp), DIMENSION(:), ALLOCATABLE :: index_nact       !< 
       
       ALLOCATE( index_nacb(1:nturbines) )
       ALLOCATE( index_nacl(1:nturbines) )
       ALLOCATE( index_nacr(1:nturbines) )
       ALLOCATE( index_nact(1:nturbines) )


       IF ( speed_control)  THEN
       
          CALL wtm_speed_control

       ENDIF

!
!------------------------------------------------------------------------------!
!--    Calculation of parameters for the regularization kernel
!--    (smearing of the forces)
!------------------------------------------------------------------------------!
!
!--    In the following, some of the required parameters for the smearing will
!--    be calculated:

!--    The kernel is set equal to twice the grid spacing which has turned out to
!--    be a reasonable value (see e.g. Troldborg et al. (2013), Wind Energy,
!--    DOI: 10.1002/we.1608):
       eps_kernel = 2.0_wp * dx
!
!--    The zero point (eps_min) of the polynomial function must be the following
!--    if the integral of the polynomial function (for values < eps_min) shall
!--    be equal to the integral of the Gaussian function used before:
       eps_min = ( 105.0_wp / 32.0_wp )**( 1.0_wp / 3.0_wp ) *                 &
                 pi**( 1.0_wp / 6.0_wp ) * eps_kernel
!
!--    Square of eps_min:
       eps_min2 = eps_min**2
!
!--    Parameters in the polynomial function:
       pol_a = 1.0_wp / eps_min**4
       pol_b = 2.0_wp / eps_min**2
!
!--    Normalization factor which is the inverse of the integral of the smearing
!--    function:
       eps_factor = 105.0_wp / ( 32.0_wp * pi * eps_min**3 )
       
!--    Change tilt angle to rad:
       tilt = tilt * pi / 180.0_wp
      
!
!--    Change yaw angle to rad:
       phi_yaw(:) = phi_yaw(:) * pi / 180.0_wp


       DO inot = 1, nturbines
!
!--       Rotate the rotor coordinates in case yaw and tilt are defined
          CALL wtm_rotate_rotor( inot )
          
!
!--       Determine the indices of the hub height
          i_hub(inot) = INT(   rcx(inot)                 / dx )
          j_hub(inot) = INT( ( rcy(inot) + 0.5_wp * dy ) / dy )
          k_hub(inot) = INT( ( rcz(inot) + 0.5_wp * dz ) / dz )

          
          SELECT CASE ( turb_mod )
             CASE ( 'adm' )
!
!--       Determining the boundaries of the rotor surface

                i_smear(inot) = CEILING( ( rr(inot) ) / dx )
                j_smear(inot) = CEILING( ( rr(inot) ) / dy )
                k_smear(inot) = CEILING( ( rr(inot) ) / dz )        
            
             CASE ( 'admr' )
!
!--       Determining the area to which the smearing of the forces is applied.
!--       As smearing now is effectively applied only for distances smaller than
!--       eps_min, the smearing area can be further limited and regarded as a
!--       function of eps_min:
             i_smear(inot) = CEILING( ( rr(inot) + eps_min ) / dx )
             j_smear(inot) = CEILING( ( rr(inot) + eps_min ) / dy ) ! TBD: can be zero?
             k_smear(inot) = CEILING( ( rr(inot) + eps_min ) / dz )
          END SELECT  
       ENDDO

!
!------------------------------------------------------------------------------!
!--    Determine the area within each grid cell that overlaps with the area
!--    of the nacelle and the tower (needed for calculation of the forces)
!------------------------------------------------------------------------------!
!
!--    Note: so far this is only a 2D version, in that the mean flow is
!--    perpendicular to the rotor area.

       
       DO inot = 1, nturbines                     ! loop over number of turbines
!
!--       Determine the grid index (u-grid) that corresponds to the location of
!--       the rotor center (reduces the amount of calculations in the case that
!--       the mean flow is perpendicular to the rotor area):
          i = i_hub(inot)

!           IF (.TRUE. == .FALSE.)  THEN
!           SELECT CASE ( turb_mod )
!              CASE ( 'adm' )
!              PRINT*, 'calc rotor surface' , inot
!                 CALL wtm_surface_rotor(inot)
!              PRINT*, 'end calc rotor surface', myid, inot
!           END SELECT
!           ENDIF      
!
!--       Determine the left and the right edge of the nacelle (corresponding
!--       grid point indices):
          index_nacl(inot) = INT( ( rcy(inot) - rnac(inot) + 0.5_wp * dy ) / dy )
          index_nacr(inot) = INT( ( rcy(inot) + rnac(inot) + 0.5_wp * dy ) / dy )
!
!--       Determine the bottom and the top edge of the nacelle (corresponding
!--       grid point indices).The grid point index has to be increased by 1, as
!--       the first level for the u-component (index 0) is situated below the
!--       surface. All points between z=0 and z=dz/s would already be contained
!--       in grid box 1.
          index_nacb(inot) = INT( ( rcz(inot) - rnac(inot) ) / dz ) + 1
          index_nact(inot) = INT( ( rcz(inot) + rnac(inot) ) / dz ) + 1

!
!--       Determine the indices of the grid boxes containing the left and 
!--       the right boundaries of the tower:
          tower_n = ( rcy(inot) + 0.5_wp * dtow(inot) - 0.5_wp * dy ) / dy
          tower_s = ( rcy(inot) - 0.5_wp * dtow(inot) - 0.5_wp * dy ) / dy

!
!--       Determine the fraction of the grid box area overlapping with the tower
!--       area and multiply it with the drag of the tower:
          IF ( ( nxlg <= i )  .AND.  ( nxrg >= i ) )  THEN

             DO  j = nys, nyn
!
!--             Loop from south to north boundary of tower
                IF ( ( j >= tower_s )  .AND.  ( j <= tower_n ) )  THEN

                   DO  k = nzb, nzt

                      IF ( k == k_hub(inot) )  THEN
                         IF ( tower_n - tower_s >= 1 )  THEN
!
!--                      leftmost and rightmost grid box:
                            IF ( j == tower_s )  THEN
                               tow_cd_surf(k,j,i) = ( rcz(inot) -              &
                                    ( k_hub(inot) * dz - 0.5_wp * dz ) )  *    & ! extension in z-direction
                                  ( ( tower_s + 1.0_wp + 0.5_wp ) * dy    -    &
                                    ( rcy(inot) - 0.5_wp * dtow(inot) ) ) *    & ! extension in y-direction
                                  turb_cd_tower(inot)
                            ELSEIF ( j == tower_n )  THEN
                               tow_cd_surf(k,j,i) = ( rcz(inot)            -   &
                                    ( k_hub(inot) * dz - 0.5_wp * dz ) )  *    & ! extension in z-direction
                                  ( ( rcy(inot) + 0.5_wp * dtow(inot) )   -    &
                                    ( tower_n + 0.5_wp ) * dy )           *    & ! extension in y-direction
                                  turb_cd_tower(inot)
!
!--                         grid boxes inbetween
!--                         (where tow_cd_surf = grid box area):
                            ELSE
                               tow_cd_surf(k,j,i) = ( rcz(inot) -              &
                                    ( k_hub(inot) * dz - 0.5_wp * dz ) )  *    &
                                    dy * turb_cd_tower(inot)
                            ENDIF
!
!--                      tower lies completely within one grid box:
                         ELSE
                            tow_cd_surf(k,j,i) = ( rcz(inot)                 - &
                                       ( k_hub(inot) * dz - 0.5_wp * dz ) ) *  &
                                       dtow(inot) * turb_cd_tower(inot)
                         ENDIF
!
!--                  In case that k is smaller than k_hub the following actions
!--                  are carried out:
                      ELSEIF ( k < k_hub(inot) )  THEN
                      
                         IF ( ( tower_n - tower_s ) >= 1 )  THEN
!
!--                         leftmost and rightmost grid box:
                            IF ( j == tower_s )  THEN                         
                               tow_cd_surf(k,j,i) = dz * (                     &
                                      ( tower_s + 1 + 0.5_wp ) * dy         -  &
                                      ( rcy(inot) - 0.5_wp * dtow(inot) )      &
                                                        ) * turb_cd_tower(inot)
                            ELSEIF ( j == tower_n )  THEN
                               tow_cd_surf(k,j,i) = dz * (                     &
                                      ( rcy(inot) + 0.5_wp * dtow(inot) )   -  &
                                      ( tower_n + 0.5_wp ) * dy                &
                                                         ) * turb_cd_tower(inot)
!
!--                         grid boxes inbetween
!--                         (where tow_cd_surf = grid box area):
                            ELSE
                               tow_cd_surf(k,j,i) = dz * dy * turb_cd_tower(inot)
                            ENDIF
!
!--                         tower lies completely within one grid box:
                         ELSE
                            tow_cd_surf(k,j,i) = dz * dtow(inot) *             &
                                                turb_cd_tower(inot)
                         ENDIF ! end if larger than grid box

                      ENDIF    ! end if k == k_hub

                   ENDDO       ! end loop over k

                ENDIF          ! end if inside north and south boundary of tower

             ENDDO             ! end loop over j

          ENDIF                ! end if hub inside domain + ghostpoints
       
          
          CALL exchange_horiz( tow_cd_surf, nbgp )

       ENDDO   ! end of loop over turbines

       tow_cd_surf   = tow_cd_surf   / ( dx * dy * dz )      ! Normalize tower drag
       nac_cd_surf = nac_cd_surf / ( dx * dy * dz )      ! Normalize nacelle drag

       SELECT CASE ( turb_mod )
         CASE ( 'adm' )
!
!-          Initialize rotor surface and read thrust curve
            CALL wtm_surface_rotor      
            CALL wtm_read_thrust_curve
            turb_area = turb_area / ( dx * dy * dz )
         CASE ( 'admr' )
!
!-          read blade geometry and airfoil tables
            CALL wtm_read_blade_tables
       END SELECT
    
    END SUBROUTINE wtm_init
    
    SUBROUTINE wtm_surface_rotor
    
       INTEGER(iwp) ::  inot !< turbine number
    
       INTEGER(iwp) ::  j  !< running index
       INTEGER(iwp) ::  k  !< running index
       
!--    Help variables for the calculaction of the rotor surface
       INTEGER(iwp) ::  i_ip         !<
       INTEGER(iwp) ::  i_ipg        !<
       
       REAL(wp) ::  yvalue               
       REAL(wp) ::  dy_int           !< 
       REAL(wp) ::  dz_int           !< 
       
       REAL(wp) ::  sqrt_arg_r           !<
       
       REAL(wp), DIMENSION(:,:), ALLOCATABLE :: circle_points  !<
    
       dy_int = dy / 10000.0       
!
!--    Allocation of the array containing information on the intersection points
!--    between rotor disk and the numerical grid:
       upper_end = ( ny + 1 ) * 10000 
       ALLOCATE( circle_points(1:2,1:upper_end) )          
       
       DO inot = 1, nturbines                     ! loop over number of turbines
       
          circle_points(:,:) = 0.0
       
          DO  i_ip = 1, upper_end
              yvalue     = dy_int*(i_ip-0.5_wp) - 0.5_wp * dy
              sqrt_arg_r = rr(inot)**2.0_wp - (yvalue-rcy(inot))**2.0_wp
              IF ( sqrt_arg_r >= 0.0_wp ) THEN
!
!--             bottom intersection point
                circle_points(1,i_ip) = rcz(inot) - SQRT( sqrt_arg_r )
! 
!--             top intersection point
                circle_points(2,i_ip) = rcz(inot) + SQRT( sqrt_arg_r )
             ELSE
                circle_points(:,i_ip) = -111111
             ENDIF
          ENDDO
     
          IF ( ( nxlg <= i_hub(inot) )  .AND.  ( nxrg >= i_hub(inot) ) )  THEN

             DO j = MAX( nys, j_hub(inot) - j_smear(inot) ),                   &
                    MIN( nyn, j_hub(inot) + j_smear(inot) )
                DO k = MAX( nzb_u_inner(j,i_hub(inot))+1,                      & 
                            k_hub(inot) - k_smear(inot) ),                     &
                            k_hub(inot) + k_smear(inot)       


!--                For all other cases Riemann integrals are calculated.
!--                Here, the points on the circle that have been determined
!--                above are used in order to calculate the overlap between the
!--                gridbox and the rotor area (area approached by 10000
!--                rectangulars dz_int * dy_int):
                   DO  i_ipg = 1, 10000
                      dz_int = dz
                      i_ip = j * 10000 + i_ipg
!
!--                   Determine the vertical extension dz_int of the circle
!--                   within the current grid box:
                      IF ( ( circle_points(2,i_ip) < zw(k) ) .AND.             &
                           ( circle_points(2,i_ip) >= zw(k-1) ) ) THEN
                         dz_int = dz_int - ( zw(k) - circle_points(2,i_ip) )
                      ENDIF
                      IF ( ( circle_points(1,i_ip) <= zw(k) ) .AND.            &
                           ( circle_points(1,i_ip) > zw(k-1) ) ) THEN
                         dz_int = dz_int - ( circle_points(1,i_ip) - zw(k-1) )
                      ENDIF
                      IF ( zw(k-1) > circle_points(2,i_ip) ) THEN
                         dz_int = 0.0_wp
                      ENDIF
                      IF ( zw(k) < circle_points(1,i_ip) ) THEN
                         dz_int = 0.0_wp
                      ENDIF 
                         turb_area(k,j,i_hub(inot)) =                          &
                                    turb_area(k,j,i_hub(inot)) + dy_int * dz_int
               
!                         PRINT*, 'myid k,j,turb area', myid, k, j, turb_area(k,j,i_hub(inot))
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
       ENDDO
      
    END SUBROUTINE wtm_surface_rotor


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Read in layout of the rotor blade , the lift and drag tables
!> and the distribution of lift and drag tables along the blade
!------------------------------------------------------------------------------!
!
    SUBROUTINE wtm_read_blade_tables


       IMPLICIT NONE

       INTEGER(iwp) ::  ii   !< running index
       INTEGER(iwp) ::  jj   !< running index
    
       INTEGER(iwp) ::  ierrn       !<
    
       CHARACTER(200) :: chmess     !< Read in string 

       INTEGER(iwp) ::  dlen        !< no. rows of local table 
       INTEGER(iwp) ::  dlenbl      !< no. rows of cd, cl table
       INTEGER(iwp) ::  ialpha      !< table position of current alpha value
       INTEGER(iwp) ::  iialpha     !< 
       INTEGER(iwp) ::  iir         !< 
       INTEGER(iwp) ::  radres      !< radial resolution
       INTEGER(iwp) ::  t1          !< no. of airfoil
       INTEGER(iwp) ::  t2          !< no. of airfoil
       INTEGER(iwp) ::  trow        !< 
       INTEGER(iwp) ::  dlenbl_int  !< no. rows of interpolated cd, cl tables
    
       REAL(wp) :: alpha_attack_i   !<
       REAL(wp) :: weight_a         !< 
       REAL(wp) :: weight_b         !< 

       INTEGER(iwp), DIMENSION(:), ALLOCATABLE :: ttoint1    !<
       INTEGER(iwp), DIMENSION(:), ALLOCATABLE :: ttoint2    !<
    
       REAL(wp), DIMENSION(:), ALLOCATABLE :: turb_cd_sel1   !< 
       REAL(wp), DIMENSION(:), ALLOCATABLE :: turb_cd_sel2   !<
       REAL(wp), DIMENSION(:), ALLOCATABLE :: turb_cl_sel1   !< 
       REAL(wp), DIMENSION(:), ALLOCATABLE :: turb_cl_sel2   !<
       REAL(wp), DIMENSION(:), ALLOCATABLE :: read_cl_cd     !< read in var array
              
       REAL(wp), DIMENSION(:), ALLOCATABLE    :: alpha_attack_tab   !< 
       REAL(wp), DIMENSION(:), ALLOCATABLE    :: trad1              !< 
       REAL(wp), DIMENSION(:), ALLOCATABLE    :: trad2              !<          
       REAL(wp), DIMENSION(:,:), ALLOCATABLE  :: turb_cd_table      !<
       REAL(wp), DIMENSION(:,:), ALLOCATABLE  :: turb_cl_table      !<
                                          
       ALLOCATE ( read_cl_cd(1:2*nairfoils+1) )

!
!--    Read in the distribution of lift and drag tables along the blade, the
!--    layout of the rotor blade and the lift and drag tables:

       OPEN ( 201, FILE='WTM_DATA', STATUS='OLD', FORM='FORMATTED', IOSTAT=ierrn )

       IF ( ierrn /= 0 )  THEN
          message_string = 'file WTM_DATA does not exist'
          CALL message( 'wtm_init', 'PA0???', 1, 2, 0, 6, 0 )
       ENDIF
!
!--    Read distribution table:

       dlen = 0

       READ ( 201, '(3/)' )

       rloop3: DO
          READ ( 201, *, IOSTAT=ierrn ) chmess
          IF ( ierrn < 0  .OR.  chmess == '#'  .OR.  chmess == '')  EXIT rloop3
          dlen = dlen + 1
       ENDDO rloop3

       ALLOCATE( trad1(1:dlen), trad2(1:dlen), ttoint1(1:dlen), ttoint2(1:dlen))

       DO jj = 1,dlen+1
          BACKSPACE ( 201, IOSTAT=ierrn )
       ENDDO

       DO jj = 1,dlen
          READ ( 201, * ) trad1(jj), trad2(jj), ttoint1(jj), ttoint2(jj)
       ENDDO

!
!--    Read layout table:

       dlen = 0 

       READ ( 201, '(3/)')

       rloop1: DO
          READ ( 201, *, IOSTAT=ierrn ) chmess
          IF ( ierrn < 0  .OR.  chmess == '#'  .OR.  chmess == '')  EXIT rloop1
          dlen = dlen + 1
       ENDDO rloop1

       ALLOCATE( lrd(1:dlen), ard(1:dlen), crd(1:dlen) )
       DO jj = 1, dlen+1
          BACKSPACE ( 201, IOSTAT=ierrn )
       ENDDO             
       DO jj = 1, dlen
          READ ( 201, * ) lrd(jj), ard(jj), crd(jj) 
       ENDDO

!
!--    Read tables (turb_cl(alpha),turb_cd(alpha) for the different profiles:

       dlen = 0

       READ ( 201, '(3/)' )

       rloop2: DO
          READ ( 201, *, IOSTAT=ierrn ) chmess
          IF ( ierrn < 0  .OR.  chmess == '#'  .OR.  chmess == '')  EXIT rloop2
          dlen = dlen + 1
       ENDDO rloop2 

       ALLOCATE( alpha_attack_tab(1:dlen), turb_cl_table(1:dlen,1:nairfoils),  &
                 turb_cd_table(1:dlen,1:nairfoils) )

       DO jj = 1,dlen+1
          BACKSPACE ( 201, IOSTAT=ierrn )
       ENDDO 

       DO jj = 1,dlen
          READ ( 201, * ) read_cl_cd
          alpha_attack_tab(jj) = read_cl_cd(1)
          DO ii= 1, nairfoils
             turb_cl_table(jj,ii) = read_cl_cd(ii*2)
             turb_cd_table(jj,ii) = read_cl_cd(ii*2+1)
          ENDDO

       ENDDO

       dlenbl = dlen 

       CLOSE ( 201 )

!
!--    For each possible radial position (resolution: 0.1 m --> 630 values) and
!--    each possible angle of attack (resolution: 0.01 degrees --> 36000 values!)
!--    determine the lift and drag coefficient by interpolating between the
!--    tabulated values of each table (interpolate to current angle of attack)
!--    and between the tables (interpolate to current radial position):

       ALLOCATE( turb_cl_sel1(0:dlenbl) )  
       ALLOCATE( turb_cl_sel2(0:dlenbl) )  
       ALLOCATE( turb_cd_sel1(0:dlenbl) )
       ALLOCATE( turb_cd_sel2(0:dlenbl) )

       radres     = INT( rr(1) * 10.0_wp ) + 1_iwp
       dlenbl_int = INT( 360.0_wp / accu_cl_cd_tab ) + 1_iwp


       ALLOCATE( turb_cl_tab(0:dlenbl_int,1:radres) )
       ALLOCATE( turb_cd_tab(0:dlenbl_int,1:radres) )


       DO iir = 1, radres ! loop over radius

          DO iialpha = 1, dlenbl_int  ! loop over angles

             cur_r = ( iir - 1_iwp ) * 0.1_wp             
             alpha_attack_i = -180.0_wp + REAL( iialpha-1 ) * accu_cl_cd_tab
             ialpha = 1

             DO WHILE ( alpha_attack_i > alpha_attack_tab(ialpha) )
                ialpha = ialpha + 1
             ENDDO
!
!--          Find position in table
             lct = MINLOC( ABS( trad1 - cur_r ) )
!                lct(1) = lct(1)

             IF ( ( trad1(lct(1)) - cur_r ) .GT. 0.0 )  THEN
                lct(1) = lct(1) - 1
             ENDIF

             trow = lct(1)
!
!--          Calculate weights for interpolation
             weight_a = ( trad2(trow) - cur_r ) / ( trad2(trow) - trad1(trow) )
             weight_b = ( cur_r - trad1(trow) ) / ( trad2(trow) - trad1(trow) )
             t1 = ttoint1(trow)
             t2 = ttoint2(trow)

             IF ( t1 .EQ. t2 ) THEN  ! if both are the same, the weights are NaN
                weight_a = 0.5_wp    ! then do interpolate in between same twice
                weight_b = 0.5_wp    ! using 0.5 as weight
             ENDIF

             IF ( t1 == 0 .AND. t2 == 0 ) THEN
                turb_cd_sel1 = 0.0_wp
                turb_cd_sel2 = 0.0_wp
                turb_cl_sel1 = 0.0_wp
                turb_cl_sel2 = 0.0_wp
             ELSE
                turb_cd_sel1 = turb_cd_table(:,t1)
                turb_cd_sel2 = turb_cd_table(:,t2)
                turb_cl_sel1 = turb_cl_table(:,t1)
                turb_cl_sel2 = turb_cl_table(:,t2)
             ENDIF

!
!--          Interpolation of lift and drag coefficiencts on fine grid of radius 
!--          segments and angles of attack

             turb_cl_tab(iialpha,iir) = ( alpha_attack_tab(ialpha) -           &
                                          alpha_attack_i ) /                   &
                                        ( alpha_attack_tab(ialpha) -           &
                                          alpha_attack_tab(ialpha-1) ) *       &
                                        ( weight_a * turb_cl_sel1(ialpha-1) +  &
                                          weight_b * turb_cl_sel2(ialpha-1) ) +&
                                        ( alpha_attack_i             -         &
                                          alpha_attack_tab(ialpha-1) ) /       &
                                        ( alpha_attack_tab(ialpha) -           &
                                          alpha_attack_tab(ialpha-1) ) *       &
                                        ( weight_a * turb_cl_sel1(ialpha) +    &
                                          weight_b * turb_cl_sel2(ialpha) )
             turb_cd_tab(iialpha,iir) = ( alpha_attack_tab(ialpha) -           &
                                          alpha_attack_i ) /                   &
                                        ( alpha_attack_tab(ialpha) -           &
                                          alpha_attack_tab(ialpha-1) ) *       &
                                        ( weight_a * turb_cd_sel1(ialpha-1) +  &
                                          weight_b * turb_cd_sel2(ialpha-1) ) +&
                                        ( alpha_attack_i             -         &
                                          alpha_attack_tab(ialpha-1) ) /       &
                                        ( alpha_attack_tab(ialpha) -           &
                                          alpha_attack_tab(ialpha-1) ) *       &
                                        ( weight_a * turb_cd_sel1(ialpha) +    &
                                          weight_b * turb_cd_sel2(ialpha) )
   
          ENDDO   ! end loop over angles of attack

       ENDDO   ! end loop over radius
    
    END SUBROUTINE wtm_read_blade_tables

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Read in layout of the rotor blade , the lift and drag tables
!> and the distribution of lift and drag tables along the blade
!------------------------------------------------------------------------------!
!      

   SUBROUTINE wtm_read_thrust_curve
   
   
      INTEGER(iwp) ::  jj   !< running index
      INTEGER(iwp) ::  ierrn       !<
      INTEGER(iwp) ::  dlen        !< no. rows of local table 
      
      CHARACTER(200) :: chmess     !< Read in string 

      OPEN ( 201, FILE='WTM_DATA', STATUS='OLD', FORM='FORMATTED', IOSTAT=ierrn )

      IF ( ierrn /= 0 )  THEN
         message_string = 'file WTM_DATA does not exist'
         CALL message( 'wtm_init', 'PA0???', 1, 2, 0, 6, 0 )
      ENDIF 

      dlen = 0 ! Bei keiner Headerzeile

!     TBD: Header auslesen
!     TBD: Fehlermeldung wenn falscher Header ( auch in admr )
      
      rloop4: DO
         READ ( 201, *, IOSTAT=ierrn ) chmess
         IF ( ierrn < 0  .OR.  chmess == '#'  .OR.  chmess == '')  EXIT rloop4
         dlen = dlen + 1
      ENDDO rloop4      
      
! 10     READ( 201, '( X )', END=11 )
!       dlen = dlen + 1
!       GOTO 10
! 11     CONTINUE

      ALLOCATE( vct(1:dlen), dct(1:dlen) )
      
      DO jj = 1, dlen+1
         BACKSPACE ( 201, IOSTAT=ierrn )
      ENDDO
         
      !READ( 201, '( X )' )

      DO jj = 1, dlen
         READ ( 201, * ) vct(jj), dct(jj)              
      ENDDO
         
      IF ( myid == 0 ) THEN
         DO jj = 1, dlen
            !PRINT *, '1', vct(jj), dct(jj)
         ENDDO  
      ENDIF
	  
      CLOSE( 201 )
   
   END SUBROUTINE
    
!------------------------------------------------------------------------------!
! Description:
! ------------
!> The projection matrix for the coordinate system of the rotor disc in respect
!> to the yaw and tilt angle of the rotor is calculated
!------------------------------------------------------------------------------!
    SUBROUTINE wtm_rotate_rotor( inot )


       IMPLICIT NONE

       INTEGER(iwp) :: inot
!
!--    Calculation of the rotation matrix for the application of the tilt to
!--    the rotors
       rot_eigen_rad(1) = SIN( phi_yaw(inot) )    ! x-component of the radial eigenvector
       rot_eigen_rad(2) = COS( phi_yaw(inot) )    ! y-component of the radial eigenvector 
       rot_eigen_rad(3) = 0.0_wp                  ! z-component of the radial eigenvector

       rot_eigen_azi(1) = 0.0_wp                  ! x-component of the azimuth eigenvector
       rot_eigen_azi(2) = 0.0_wp                  ! y-component of the azimuth eigenvector
       rot_eigen_azi(3) = 1.0_wp                  ! z-component of the azimuth eigenvector

       rot_eigen_nor(1) =  COS( phi_yaw(inot) )   ! x-component of the normal eigenvector
       rot_eigen_nor(2) = -SIN( phi_yaw(inot) )   ! y-component of the normal eigenvector
       rot_eigen_nor(3) = 0.0_wp                  ! z-component of the normal eigenvector
    
!
!--    Calculation of the coordinate transformation matrix to apply a tilt to
!--    the rotor. If tilt = 0, rot_coord_trans is a unit matrix.

       rot_coord_trans(inot,1,1) = rot_eigen_rad(1)**2                   *     &
                                   ( 1.0_wp - COS( tilt ) ) + COS( tilt ) 
       rot_coord_trans(inot,1,2) = rot_eigen_rad(1) * rot_eigen_rad(2)   *     &
                                   ( 1.0_wp - COS( tilt ) )              -     &
                                   rot_eigen_rad(3) * SIN( tilt )
       rot_coord_trans(inot,1,3) = rot_eigen_rad(1) * rot_eigen_rad(3)   *     &
                                   ( 1.0_wp - COS( tilt ) )              +     &
                                   rot_eigen_rad(2) * SIN( tilt )
       rot_coord_trans(inot,2,1) = rot_eigen_rad(2) * rot_eigen_rad(1)   *     &
                                   ( 1.0_wp - COS( tilt ) )              +     &
                                   rot_eigen_rad(3) * SIN( tilt )
       rot_coord_trans(inot,2,2) = rot_eigen_rad(2)**2                   *     &
                                   ( 1.0_wp - COS( tilt ) ) + COS( tilt ) 
       rot_coord_trans(inot,2,3) = rot_eigen_rad(2) * rot_eigen_rad(3)   *     &
                                   ( 1.0_wp - COS( tilt ) )              -     &
                                   rot_eigen_rad(1) * SIN( tilt )
       rot_coord_trans(inot,3,1) = rot_eigen_rad(3) * rot_eigen_rad(1)   *     &
                                   ( 1.0_wp - COS( tilt ) )              -     &
                                   rot_eigen_rad(2) * SIN( tilt )
       rot_coord_trans(inot,3,2) = rot_eigen_rad(3) * rot_eigen_rad(2)   *     &
                                   ( 1.0_wp - COS( tilt ) )              +     &
                                   rot_eigen_rad(1) * SIN( tilt )
       rot_coord_trans(inot,3,3) = rot_eigen_rad(3)**2                   *     &
                                   ( 1.0_wp - COS( tilt ) ) + COS( tilt )

!
!--    Vectors for the Transformation of forces from the rotor's spheric
!--    coordinate system to the cartesian coordinate system
       rotx(inot,:) = MATMUL( rot_coord_trans(inot,:,:), rot_eigen_nor )
       roty(inot,:) = MATMUL( rot_coord_trans(inot,:,:), rot_eigen_rad )
       rotz(inot,:) = MATMUL( rot_coord_trans(inot,:,:), rot_eigen_azi )
    
    END SUBROUTINE wtm_rotate_rotor


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculation of the forces generated by the wind turbine
!------------------------------------------------------------------------------!
   SUBROUTINE wtm_forces
   
      IF ( simulated_time >= time_turbine_on ) THEN
    
         SELECT CASE ( turb_mod )

            CASE ( 'adm' )
                
                CALL wtm_forces_adm
            
            CASE ( 'admr' )
                
                CALL wtm_forces_admr
                
         END SELECT        
    
      
         !IF ( simulated_time >= 1000 .AND. simulated_time <= 1100) THEN ! (SB)
         !  u = u + 0.004*u  
         !ENDIF

         !IF ( simulated_time >= 2000 .AND. simulated_time <= 2100) THEN ! (SB)
         !  u = u - 0.005*u  
         !ENDIF
         

         !IF ( simulated_time >= 3000 .AND. simulated_time <= 3100) THEN ! (SB)
         !  v = v - 0.01  
         !ENDIF


      ENDIF
   
   END SUBROUTINE wtm_forces
   
   
   SUBROUTINE wtm_forces_adm

      IMPLICIT NONE

      CHARACTER (LEN=2) ::  turbine_id
      
      LOGICAL      ::  uopt_exists = .FALSE. !< 

      INTEGER(iwp) ::  i, j, k          !< loop indices
      INTEGER(iwp) ::  inot             !< turbine loop index (turbine id)
      INTEGER(iwp), DIMENSION(1) ::  lct              !< 
      
      REAL(wp), DIMENSION(1:100) ::  ufree
      REAL(wp), DIMENSION(1)     ::  ufree_d
      REAL(wp)                   ::  u2_ind      !< induction velocity
      
      thrust_rotor = 0.0_wp            !< reset rotor thrust for every time step
      
!
!-- wait until control signals are present (SB)
    IF ( MOD(current_timestep_number,2) == 0 ) THEN ! read control signal at t={0,1,2,...,N} 
      

      DO WHILE ( .NOT.uopt_exists )  ! wait for text files with control signals (in matlab, first make ui.txt, and then uflag.txt)      
          INQUIRE(FILE=TRIM(str)//"uflag.txt", EXIST=uopt_exists)   
      ENDDO                       
    
    ENDIF

!
!--          Evaluation of the rotor disk-averaged velocity and the velocity at
!--          the position of the hub. By default, the rotor disk-averaged
!--          velocity is used as reference velocity in the calculation of the
!--          rotor thrust:
!
!--          Loop over number of turbines contained in the current subdomain:

      DO inot = 1, nturbines

         sums_u(inot) = 0.0_wp

!
!--             Only those points that are really found on the disk area are
!--             regarded for the calculation of the disk averaged velocity:

         IF ( ( nxl <= i_hub(inot) )  .AND.  ( nxr >= i_hub(inot) ) )  THEN
         
            DO j = MAX( nys, j_hub(inot) - j_smear(inot) ),                &
                    MIN( nyn, j_hub(inot) + j_smear(inot) )
                DO k = MAX( nzb_u_inner(j,i_hub(inot))+1, k_hub(inot) - k_smear(inot) ), &
                            k_hub(inot) + k_smear(inot)
    !
    !--               ONLY FOR FLOW FROM X-DIRECTION AND SMALL INCLINATIONS 
    !--             
                sums_u(inot) = sums_u(inot) +                                   &
                               u(k,j,i_hub(inot))**3 * turb_area(k,j,i_hub(inot))
                                    
                ENDDO
            ENDDO
         
         ENDIF
         
      ENDDO         !-- end loop over number of turbines
      
      ! u_rotor is at sample k-1  (SB) 
      DO inot = 1, nturbines
         IF ( ( nxl <= i_hub(inot) )  .AND.  ( nxr >= i_hub(inot) ) )   &
           THEN
           IF ( ( nys <= j_hub(inot) )  .AND.  ( nyn >= j_hub(inot) ) )&
              THEN 
                u_rotor_old(inot) = u_rotor(inot)
           ENDIF
         ENDIF 
      ENDDO  
      
      CALL MPI_ALLREDUCE( sums_u, u_rotor, nturbines, MPI_REAL, &
                          MPI_SUM, comm2d, ierr )
      

      DO inot = 1,nturbines
        u_rotor(inot) = u_rotor(inot) * dx * dy * dz /                      &
                         ( rr(inot) * rr(inot) * pi )
        u_rotor(inot) = u_rotor(inot)**(1./3.)
      ENDDO
                          
      ! u_rotor is at sample k (SB)

      !DO inot = 1, nturbines
      !   IF ( ( nxl <= i_hub(inot) )  .AND.  ( nxr >= i_hub(inot) ) )   &
      !     THEN
      !     IF ( ( nys <= j_hub(inot) )  .AND.  ( nyn >= j_hub(inot) ) )&
      !        THEN 
      !          PRINT*, u_rotor(inot), u_rotor_old(inot)
      !     ENDIF
      !   ENDIF 
      !ENDDO 


!
!-- Calculation of the axial induction factor from the thrust table
      DO inot = 1, nturbines
         ufree(inot) = u_rotor(inot) / ( 1 - aif(inot) )
         
         !PRINT *, '1', aif(inot)

         lct = MINLOC( ABS( ufree(inot) - vct ) )
         ufree_d = ufree(inot) - vct(lct)

         IF ( ufree_d(1) > 0 )  THEN
            turb_ct(inot)  = ( 1-ufree_d(1) ) * dct ( lct(1) ) +               &
                              ufree_d(1) * dct( lct(1)+1 )
         ELSEIF ( ufree_d(1) < 0 )  THEN
            turb_ct(inot) = (1+ufree_d(1)) * dct( lct(1) ) +                   & 
                              ABS( ufree_d(1) ) * dct( lct(1)-1 ) 
         ELSEIF ( ufree_d(1) == 0 )  THEN
            turb_ct(inot) = dct( lct(1) )
         ENDIF

!
!-- Add time series to turb_ct (SB)          
		 
         !CALL wtm_ct_from_matlab( inot )
!
!-- Read control signal from external file (SB)
                
                ! read control signals from text files
                WRITE ( turbine_id,'(I2.2)')  inot

                OPEN( 50, FILE=TRIM(str)//"u" &
                    //TRIM(turbine_id)//".txt",FORM="FORMATTED",STATUS="OLD",ACTION="READ" )  

                READ(50,*) turb_ct(inot)
                !turb_ct(inot) = 0.55

                CLOSE( 50 )


!
!--      Determination of the axial induction factor 
         aif(inot) = -1.0_wp *                                                 & 
                     SQRT( 0.25_wp - ( turb_ct(inot) / 4.0_wp ) ) +            &
                     0.5_wp         

      ENDDO             !-- end loop over number of turbines

!
!-- Calculation of thrust forces
!-- Loop over number of turbines    

      thrust_rotor(:) = 0.0_wp
      power_rotor(:)  = 0.0_wp
    
      DO inot = 1, nturbines
         DO i = MAX( nxl, i_hub(inot) - i_smear(inot) ),                       &
               MIN( nxr, i_hub(inot) + i_smear(inot) )
            DO j = MAX( nys, j_hub(inot) - j_smear(inot) ),                    &
                  MIN( nyn, j_hub(inot) + j_smear(inot) )
               DO k = MAX( nzb_u_inner(j,i)+1, k_hub(inot) - k_smear(inot) ),  &
                           k_hub(inot) + k_smear(inot)
!
!--               Calculate and apply the thrust force to the flow
!--               only for positive flow velocities to avoid a
!--               negative feedback and finally the crash of the
!--               simulation due to extremly high (negative)
!--               velocities:

!
!--               Calculate the thrust force of the rotor disk:
!
!--               Options for the induction velocity
!
!--               (1) u_free (Constant over the rotor surface)
!--                   equivalent to constant thrust over rotor surface
!                  u2_ind = ( u_free(inot) *                                 &
!                           cos(phi_yaw(inot)) * cos(tilt) )**2               
!--               (2) u local corrected by the induction factor
!--                   equivalent to constant ct over rotor surface
                  u2_ind = ( u(k,j,i) / ( 1 - aif(inot) ) *                    &
                            cos(phi_yaw(inot)) * cos(tilt) )**2
!
!--               Exclude negative velocities
                  u2_ind = MAX( u2_ind, 0.0_wp )
!
!--               Calculate tendencies
                  rot_tend_x(k,j,i) = 0.5_wp * turb_ct(inot) *                 &
                           u2_ind *                                            &
                           turb_area(k,j,i) *                                  &
                           cos(phi_yaw(inot)) * cos(tilt)
               
                  rot_tend_y(k,j,i) = 0.5_wp * turb_ct(inot) *                 &
                           u2_ind *                                            &
                           turb_area(k,j,i) *                                  &
                           (-1.0_wp) * sin(phi_yaw(inot))
                           
                  rot_tend_z(k,j,i) = 0.5_wp * turb_ct(inot) *                 &
                           u2_ind *                                        &
                           turb_area(k,j,i) *                                  &
                           cos(phi_yaw(inot)) * (-1.0_wp) * sin(tilt)
!
!--               Sum up the rotor thrust over the whole rotor disk
!--               (for time-series output of total rotor thrust):
                  thrust_rotor(inot) = thrust_rotor(inot) +                    & 
                                       rot_tend_x(k,j,i) * dx * dy * dz 
               
               ENDDO
            ENDDO
         ENDDO
            power_rotor(inot)  = power_rotor( inot ) +                         &
                                 thrust_rotor( inot ) * u_rotor( inot )

      ENDDO                       !-- end loop over number of turbines


!
!-- Write measurement to external file (SB)
!-- Find processor at i_hub, j_hub  
    IF ( MOD(current_timestep_number,2) == 1 ) THEN ! write out measurement at t={0.5,1.5,2.5,...,N.5} 
      
     IF ( uopt_exists ) THEN ! delete the uflag.txt
       OPEN( UNIT=50, FILE =TRIM(str)//"uflag.txt",FORM="FORMATTED" )
       CLOSE( 50 , STATUS="DELETE")
       uopt_exists = .FALSE.
     ENDIF

      DO inot = 1, nturbines             

        IF ( ( nxl <= i_hub(inot) )  .AND.  ( nxr >= i_hub(inot) ) )   &
           THEN
           IF ( ( nys <= j_hub(inot) )  .AND.  ( nyn >= j_hub(inot) ) )&
              THEN
                
                ! power_rotor(inot) is at sample k and power_adm(inot) at samlpe k-1 before calling MPI_ALLREDUCE
                !PRINT*, '1', simulated_time, power_rotor(inot), power_adm(inot)
                !PRINT*, '2', simulated_time, thrust_rotor(inot), thrust_adm(inot) 

                WRITE ( turbine_id,'(I2.2)')  inot

                OPEN( 50, FILE= TRIM(str)//"y" &
                    //TRIM(turbine_id)//".txt",FORM="FORMATTED",STATUS="REPLACE",ACTION="WRITE" )  

                WRITE( 50, *) simulated_time-dt_3d, inot, thrust_adm(inot), power_adm(inot), u_rotor_old(inot) 
                CLOSE( 50 )

           ENDIF
        ENDIF
      ENDDO
    ENDIF


      CALL MPI_ALLREDUCE( thrust_rotor, thrust_adm, nturbines, MPI_REAL, &
                          MPI_SUM, comm2d, ierr )
      CALL MPI_ALLREDUCE( power_rotor, power_adm, nturbines, MPI_REAL, &
                          MPI_SUM, comm2d, ierr )

!-- Printing out the outputs of the adm model     

      DO inot = 1, nturbines
         
         IF (myid == 0) THEN
            IF ( openfile_turb_mod(400+inot)%opened ) THEN
               WRITE (400+inot, 106) simulated_time, u_rotor(inot), ufree(inot), &
                              turb_ct(inot), aif(inot), phi_yaw(inot)*180.0_wp/pi,  &
                              thrust_adm(inot), power_adm(inot)                                   
                                  
            ELSE

               WRITE (turbine_id, '(I2.2)') inot
               OPEN ( 400+inot, FILE=( 'TURBINE_PARAMETERS'//turbine_id ), &
                                        FORM='FORMATTED' )
               WRITE (400+inot, 105) inot
               WRITE (400+inot, 106) simulated_time, u_rotor(inot), ufree(inot), &
                              turb_ct(inot), aif(inot), phi_yaw(inot)*180.0_wp/pi,  &
                              thrust_adm(inot), power_adm(inot)                                
                                 
            ENDIF
          ENDIF

       openfile_turb_mod(400+inot)%opened = .TRUE.
       ENDDO

       105 FORMAT ('Actuator disc model data for turbine ', I2,1X,':'/ &
               &'--------------------------------------------'/ &
               &'    Time      Ur        Uinf  ', &
                '    Ct_adm    a        Yaw(deg)', &
                '  Thrust     Power  '   )

       106 FORMAT (F9.3,1X,F9.3,1X,F9.3,1X,F9.3,1X,F9.3,1X,F9.3,1X,F9.1,2X,F9.1,2X,F9.1)


   END SUBROUTINE wtm_forces_adm

! Change thrust coefficient in time with time series defined in Matlab (SB)
   
   SUBROUTINE wtm_ct_from_matlab ( inot )
      
      IMPLICIT NONE
      
      INTEGER(iwp)    :: inot              !< turbine loop index (turbine id)

      ! A=array with excitation signal from Matlab. A(u1(t1), u2(t1),...,uN(t1),u1(t2),....)
      ! A = [u1(t1) u1(t2) ... u1(tN) ; u2(t1) u2(t2) ... u2(tN)] 
      ! This A matrix contains the input signals and is generated with generate_new_wt_model.m 
      ! Run generate_new_wt_model.m and copy the new generated wind_turbine_model_mod.f90 into USER_CODE 


       REAL(wp), DIMENSION(6,1) :: CT

       CT = RESHAPE((/0.89,0.89,0.89,0.89,0.89,0.89/),SHAPE(CT))
     
      
       turb_ct(inot) = CT(inot,1_iwp)

   END SUBROUTINE wtm_ct_from_matlab
      
   SUBROUTINE wtm_forces_admr


       IMPLICIT NONE

       CHARACTER (LEN=2) ::  turbine_id
       
       LOGICAL      ::  uopt_exists = .FALSE. !< 

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
!-- wait until control signals are present (SB)
    IF ( MOD(current_timestep_number,2) == 0 ) THEN ! read control signal at t={0,h,2h,...,Nh} 
      
      DO WHILE ( .NOT.uopt_exists )  ! wait for text files with control signals (in matlab, first make ui.txt, and then uflag.txt)      
          INQUIRE(FILE=TRIM(str)//"uflag.txt", EXIST=uopt_exists)   
      ENDDO                     
    
    ENDIF

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

                      !CALL wtm_yawcontrol( inot )

                      !-- Add time series (SB)                   
                      CALL wtm_yaw_from_matlab( inot )

                      phi_yaw_l(inot) = phi_yaw(inot)

                   ENDIF
                ENDIF
!
!--             Write measurement to external file (SB)
                IF ( MOD(current_timestep_number,2) == 1 ) THEN ! write out measurement at t={h,3h,5h,...,N} 
      
                  IF ( uopt_exists ) THEN ! delete the uflag.txt
                    OPEN( UNIT=50, FILE =TRIM(str)//"uflag.txt",FORM="FORMATTED" )
                    CLOSE( 50 , STATUS="DELETE")
                    uopt_exists = .FALSE.
                  ENDIF
               
                  IF ( inot==1 ) THEN  
                    CALL write_variables  
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

!
!-- Change thrust coefficient in time with time series defined in Matlab (SB)
   
    SUBROUTINE wtm_yaw_from_matlab ( inot )
      
      IMPLICIT NONE

      CHARACTER (LEN=2) :: turbine_id
      INTEGER(iwp)      :: inot              !< turbine loop index (turbine id)
   
!
!-- Read control signal from external file (SB)
                
       WRITE ( turbine_id,'(I2.2)')  inot

       OPEN( 50, FILE=TRIM(str)//"u" &
           //TRIM(turbine_id)//".txt",FORM="FORMATTED",STATUS="OLD",ACTION="READ" )  

       READ(50,*) phi_yaw(inot)

       CLOSE( 50 )
!

       !PRINT*, inot , phi_yaw(inot)
       
   END SUBROUTINE wtm_yaw_from_matlab

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
!--    5, the turbine is yawed. The mechanism stops as soon as the 2s-running
!--    mean of the missalignment is smaller than 0.5.
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
    
! Description:
! ------------
!> Masked data output in .txt for first masks defined in _pd3 (SB). 
!------------------------------------------------------------------------------!
 SUBROUTINE write_variables

    USE arrays_3d,                                                             &
        ONLY:  u,v
        
    USE control_parameters,                                                    &
        ONLY:  domask, domask_no, mask_i, mask_j, mask_k,                      &
               mask_size, mask_size_l, mask_start_l, mid,                      &
               current_timestep_number, simulated_time
        
    USE kinds
    
    USE pegrid

    IMPLICIT NONE
    
    INTEGER(iwp) ::  av          !< 
    INTEGER(iwp) ::  ngp         !< 
    INTEGER(iwp) ::  i           !< 
    INTEGER(iwp) ::  if          !< 
    INTEGER(iwp) ::  j           !< 
    INTEGER(iwp) ::  k           !< 
    INTEGER(iwp) ::  n           !< 
    INTEGER(iwp) ::  sender      !< 
    INTEGER(iwp) ::  ind(6)      !< 
    INTEGER(iwp) ::  Np = 2.0

    CHARACTER(LEN=9)  ::  filename_out          !< (SB)

    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  local_pf        !<
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  total_pf        !<
    REAL(wp), DIMENSION(:,:,:), POINTER     ::  to_be_resorted  !<

       
!
! Write measurement after every Np samples to folder

      mid  = 1_iwp
      av   = 0_iwp

!-- Return, if nothing to output
      IF ( domask_no(mid,av) == 0 )  RETURN
!
!-- Allocate total and local output arrays.
      IF ( myid == 0 )  THEN
         ALLOCATE( total_pf(mask_size(mid,1),mask_size(mid,2),mask_size(mid,3)) )
      ENDIF

      ALLOCATE( local_pf(mask_size_l(mid,1),mask_size_l(mid,2), &
                       mask_size_l(mid,3)) )

      IF ( myid == 0  )  THEN
        OPEN( 50, FILE=TRIM(str)//"time.txt" & 
          ,FORM="FORMATTED",STATUS="REPLACE",ACTION="WRITE" )                  
        WRITE( 50, "(F9.3)", ADVANCE='NO') simulated_time 
        CLOSE( 50 )
      ENDIF
!
!-- Loop over all variables to be written.
      if = 1

      DO  WHILE ( domask(mid,av,if)(1:1) /= ' ' )
!
!--    Reallocate local_pf on PE 0 since its shape changes during MPI exchange
         IF ( myid == 0  .AND.  if > 1 )  THEN
            DEALLOCATE( local_pf )
            ALLOCATE( local_pf(mask_size_l(mid,1),mask_size_l(mid,2), &
                             mask_size_l(mid,3)) )
         ENDIF
!
!--    Store the variable chosen.
         SELECT CASE ( TRIM( domask(mid,av,if) ) )
          
            CASE ( 'u' )
               to_be_resorted => u                          
               filename_out = "u"
          
            CASE ( 'v' )
               to_be_resorted => v
               filename_out = "v"
           
            CASE DEFAULT

               PRINT*, "Mask is not handled."

         END SELECT

!
!--    Resort the array to be output
         DO  i = 1, mask_size_l(mid,1)
           DO  j = 1, mask_size_l(mid,2)
             DO  k = 1, mask_size_l(mid,3)
               local_pf(i,j,k) =  to_be_resorted(mask_k(mid,k), &
                                   mask_j(mid,j),mask_i(mid,i))
             ENDDO
           ENDDO
         ENDDO

!--      PE0 receives partial arrays from all processors of the respective mask
!--      and outputs them. Here a barrier has to be set, because otherwise 
!--      "-MPI- FATAL: Remote protocol queue full" may occur.
         CALL MPI_BARRIER( comm2d, ierr )

         ngp = mask_size_l(mid,1) * mask_size_l(mid,2) * mask_size_l(mid,3)
         IF ( myid == 0 )  THEN
!
!--        Local array can be relocated directly.
           total_pf( &
            mask_start_l(mid,1):mask_start_l(mid,1)+mask_size_l(mid,1)-1, &
            mask_start_l(mid,2):mask_start_l(mid,2)+mask_size_l(mid,2)-1, &
            mask_start_l(mid,3):mask_start_l(mid,3)+mask_size_l(mid,3)-1 ) &
               = local_pf
!
!--         Receive data from all other PEs.
            DO  n = 1, numprocs-1
!
!--             Receive index limits first, then array.
!--             Index limits are received in arbitrary order from the PEs.
                CALL MPI_RECV( ind(1), 6, MPI_INTEGER, MPI_ANY_SOURCE, 0,  &
                     comm2d, status, ierr )
!
!--             Not all PEs have data for the mask
                IF ( ind(1) /= -9999 )  THEN
                   ngp = ( ind(2)-ind(1)+1 ) * (ind(4)-ind(3)+1 ) *  &
                         ( ind(6)-ind(5)+1 )
                   sender = status(MPI_SOURCE)
                   DEALLOCATE( local_pf )
                   ALLOCATE(local_pf(ind(1):ind(2),ind(3):ind(4),ind(5):ind(6)))
                   CALL MPI_RECV( local_pf(ind(1),ind(3),ind(5)), ngp,  &
                        MPI_REAL, sender, 1, comm2d, status, ierr )
                   total_pf(ind(1):ind(2),ind(3):ind(4),ind(5):ind(6)) &
                        = local_pf
                ENDIF
            ENDDO
                                  
            OPEN( 50, FILE=TRIM(str)//TRIM(filename_out)// & 
                "k.txt",FORM="FORMATTED",STATUS="REPLACE",ACTION="WRITE" )
       
            DO i=1,mask_size(mid,1),1 
              DO j=1,mask_size(mid,2),1
               ! The z-index of the hub in total_pf will be corrected by the Matlab program
                      WRITE( 50, "(F7.3)", ADVANCE="NO") total_pf(i,j,5)
              ENDDO
              WRITE( 50, *) ' '
            ENDDO
            !WRITE( 50, *) '#' !> to end the file and Matlab knows it can continue
            CLOSE( 50 ) 
            
          ELSE
!
!--          If at least part of the mask resides on the PE, send the index
!--          limits for the target array, otherwise send -9999 to PE0.
             IF ( mask_size_l(mid,1) > 0 .AND.  mask_size_l(mid,2) > 0 .AND. &
                  mask_size_l(mid,3) > 0  ) &
                  THEN
                ind(1) = mask_start_l(mid,1)
                ind(2) = mask_start_l(mid,1) + mask_size_l(mid,1) - 1
                ind(3) = mask_start_l(mid,2)
                ind(4) = mask_start_l(mid,2) + mask_size_l(mid,2) - 1
                ind(5) = mask_start_l(mid,3)
                ind(6) = mask_start_l(mid,3) + mask_size_l(mid,3) - 1
             ELSE
                ind(1) = -9999; ind(2) = -9999
                ind(3) = -9999; ind(4) = -9999
                ind(5) = -9999; ind(6) = -9999
             ENDIF
             CALL MPI_SEND( ind(1), 6, MPI_INTEGER, 0, 0, comm2d, ierr )
!
!--          If applicable, send data to PE0.
             IF ( ind(1) /= -9999 )  THEN
                CALL MPI_SEND( local_pf(1,1,1), ngp, MPI_REAL, 0, 1, comm2d, &
                     ierr )
             ENDIF
          ENDIF
!
!--       A barrier has to be set, because otherwise some PEs may proceed too
!--       fast so that PE0 may receive wrong data on tag 0.
          CALL MPI_BARRIER( comm2d, ierr )

        if = if + 1

      ENDDO ! loop over the masks
!
!-- Deallocate temporary arrays.
      DEALLOCATE( local_pf )
      IF ( myid == 0 )  THEN
         DEALLOCATE( total_pf )
      ENDIF

 
 END SUBROUTINE write_variables

 END MODULE wind_turbine_model_mod
