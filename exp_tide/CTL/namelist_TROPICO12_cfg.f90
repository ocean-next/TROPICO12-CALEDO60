!!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!! NEMO/OCE  Configuration namelist : overwrite default values defined in SHARED/namelist_ref
!!
!!    TROPICO12 (AGRIF MOTHER domain)
!!
!!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!! NEMO/OCE  :  1 - Domain & run manager (namrun, namcfg, namdom, namtsd, namcrs, namc1d, namc1d_uvd)
!! namelists    2 - Surface boundary (namsbc, namsbc_flx, namsbc_blk, namsbc_cpl,
!!                                    namsbc_sas, namtra_qsr, namsbc_rnf,
!!                                    namisf, namsbc_apr, 
!!                                    namsbc_ssr, namsbc_wave, namberg)
!!              3 - lateral boundary (namlbc, namagrif, nambdy, nambdy_tide)
!!              4 - top/bot boundary (namdrg, namdrg_top, namdrg_bot, nambbc, nambbl)
!!              5 - Tracer           (nameos, namtra_adv, namtra_ldf, namtra_eiv, namtra_dmp)
!!              6 - dynamics         (namdyn_adv, namdyn_vor, namdyn_hpg, namdyn_spg, namdyn_ldf)
!!              7 - Vertical physics (namzdf, namzdf_ric, namzdf_tke, namzdf_gls, namzdf_iwm)
!!              8 - diagnostics      (namnc4, namtrd, namspr, namflo, namhsb)
!!              9 - Obs & Assim      (namobs, nam_asminc)
!!             10 - miscellaneous    (nammpp, namctl, namsto)
!!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

!!======================================================================
!!              ***  Domain & Run management namelists  ***           !!
!!                                                                    !!
!!   namrun       parameters of the run
!!   namdom       space and time domain
!!   namcfg       parameters of the configuration                       (default: user defined GYRE)
!!   namwad       Wetting and drying                                    (default: OFF)
!!   namtsd       data: temperature & salinity                          (default: OFF)
!!   namcrs       coarsened grid (for outputs and/or TOP)               (ln_crs =T)
!!   namc1d       1D configuration options                              ("key_c1d")
!!   namc1d_dyndmp 1D newtonian damping applied on currents             ("key_c1d")
!!   namc1d_uvd   1D data (currents)                                    ("key_c1d")
!!======================================================================
!
!-----------------------------------------------------------------------
&namrun        !   parameters of the run
!-----------------------------------------------------------------------
   cn_exp      =  'TROPICO12_NST-TRPC12NT0'  !  experience name
   nn_it000    =    105265   !  first time step
   nn_itend    =    109728   !  last  time step (std 5840)
   nn_date0    =    20140101   !  date at nit_0000 (format yyyymmdd) used if ln_rstart=F or (ln_rstart=T and nn_rstctl=0 or 1)
   nn_leapy    =       1   !  Leap year calendar (1) or not (0)
   ln_rstart   =  .true.   !  start from rest (F) or from a restart file (T)
      ln_1st_euler = .false.  !  =T force a start with forward time step (ln_rstart=T)
      nn_rstctl   = 2      !  restart control ==> activated only if ln_rstart=T
      !                          !    = 0 nn_date0 read in namelist ; nn_it000 : read in namelist
      !                          !    = 1 nn_date0 read in namelist ; nn_it000 : check consistancy between namelist and restart
      !                          !    = 2 nn_date0 read in restart  ; nn_it000 : check consistancy between namelist and restart
      cn_ocerst_in    = 'TROPICO12_NST-TRPC12NT0_00105264_restart_oce'   !  suffix of ocean restart name (input)
      cn_ocerst_indir = '/scratch/cnt0024/ige2071/brodeau/NEMO/TROPICO12/TROPICO12_NST-TRPC12NT0-R/00105264'         !  directory from which to read input ocean restarts
      cn_ocerst_out   = 'restart_oce'   !  suffix of ocean restart name (output)
      cn_ocerst_outdir = '/scratch/cnt0024/ige2071/brodeau/NEMO/TROPICO12/TROPICO12_NST-TRPC12NT0-R/00109728'         !  directory in which to write output ocean restarts
   nn_istate   =       0   !  output the initial state (1) or not (0)
   ln_rst_list = .false.   !  output restarts at list of times using nn_stocklist (T) or at set frequency with nn_stock (F)
   nn_stock    = 4464  ! 1year @ 5400 s  frequency of creation of a restart file (modulo referenced to 1)
   nn_write    = 4464   ! 1year @ 5400 s   frequency of write in the output file   (modulo referenced to nn_it000)
   ln_mskland  = .false.   !  mask land points in NetCDF outputs
   ln_clobber  = .true.    !  clobber (overwrite) an existing file
/
!-----------------------------------------------------------------------
&namdom        !   time and space domain
!-----------------------------------------------------------------------
   rn_Dt       = 600.     !  time step for the dynamics and tracer
   ln_meshmask = .false.   !  =T create a mesh file
/
!-----------------------------------------------------------------------
&namcfg        !   parameters of the configuration                      (default: use namusr_def in namelist_cfg)
!-----------------------------------------------------------------------
   ln_read_cfg = .true.    !  (=T) read the domain configuration file
      cn_domcfg = 'domain_cfg'  ! domain configuration filename
/
!-----------------------------------------------------------------------
&namtsd        !    Temperature & Salinity Data  (init/dmp)             (default: OFF)
!-----------------------------------------------------------------------
   !                       ! =T  read T-S fields for:
   ln_tsd_init = .false.    !  ocean initialisation
   ln_tsd_dmp  = .false.         !  T-S restoring   (see namtra_dmp)

   cn_dir      = './'      !  root directory for the T-S data location
   !___________!_________________________!___________________!___________!_____________!________!___________!__________________!__________!_______________!
   !           !  file name              ! frequency (hours) ! variable  ! time interp.!  clim  ! 'yearly'/ ! weights filename ! rotation ! land/sea mask !
   !           !                         !  (if <0  months)  !   name    !   (logical) !  (T/F) ! 'monthly' !                  ! pairing  !    filename   !
   sn_tem = 'votemper_GLORYS2V4-TROPICO12_L125_201201', -12.   ,'votemper' ,   .false.   , .true. , 'yearly'  ,    ''    ,    ''    ,    ''
   sn_sal = 'vosaline_GLORYS2V4-TROPICO12_L125_201201', -12.   ,'vosaline' ,   .false.   , .true. , 'yearly'  ,    ''    ,    ''    ,    ''
/
!-----------------------------------------------------------------------
&namwad        !   Wetting and Drying (WaD)                             (default: OFF)
!-----------------------------------------------------------------------
/
!-----------------------------------------------------------------------
&namcrs        !   coarsened grid (for outputs and/or TOP)              (ln_crs =T)
!-----------------------------------------------------------------------
/
!-----------------------------------------------------------------------
&namc1d        !   1D configuration options                             ("key_c1d" default: PAPA station)
!-----------------------------------------------------------------------
/
!-----------------------------------------------------------------------
&namc1d_dyndmp !   U & V newtonian damping                              ("key_c1d" default: OFF)
!-----------------------------------------------------------------------
/
!-----------------------------------------------------------------------
&namc1d_uvd    !   data: U & V currents                                 ("key_c1d" default: OFF)
!-----------------------------------------------------------------------
/

!!======================================================================
!!            ***  Surface Boundary Condition namelists  ***          !!
!!                                                                    !!
!!   namsbc          surface boundary condition manager                 (default: NO selection)
!!   namsbc_flx      flux               formulation                     (ln_flx     =T)
!!   namsbc_blk      Bulk formulae formulation                          (ln_blk     =T)
!!   namsbc_cpl      CouPLed            formulation                     ("key_oasis3" )
!!   namsbc_sas      Stand-Alone Surface module                         (SAS_SRC  only)
!!   namsbc_iif      Ice-IF: use observed ice cover                     (nn_ice = 1   )
!!   namtra_qsr      penetrative solar radiation                        (ln_traqsr  =T)
!!   namsbc_ssr      sea surface restoring term (for T and/or S)        (ln_ssr     =T)
!!   namsbc_rnf      river runoffs                                      (ln_rnf     =T)
!!   namsbc_apr      Atmospheric Pressure                               (ln_apr_dyn =T)
!!   namsbc_wave     external fields from wave model                    (ln_wave    =T)
!!   namberg         iceberg floats                                     (ln_icebergs=T)
!!   namsbc_fwb      freshwater-budget adjustment                       (nn_fwb > 0)
!!======================================================================
!
!-----------------------------------------------------------------------
&namsbc        !   Surface Boundary Condition manager                   (default: NO selection)
!-----------------------------------------------------------------------
   nn_fsbc     = 6   !  frequency of SBC module call
!!!   nn_fsbc     =  1        !#lolo: if not skin param  Frequency? of SBC module call
                           !     (also = the frequency of sea-ice & iceberg model call)
                     ! Type of air-sea fluxes
   ln_blk      = .true.    !  Bulk formulation                          (T => fill namsbc_blk )
                     ! Sea-ice :
   nn_ice      = 0         ! NOTHING !!!
                     ! Misc. options of sbc :
   ln_traqsr   = .true.    !  Light penetration in the ocean            (T => fill namtra_qsr)
   ln_dm2dc    = .false.   !#lolo  daily mean to diurnal cycle on short wave
   ln_ssr      = .false.   !#lolo  Sea Surface Restoring on T and/or S       (T => fill namsbc_ssr)
   nn_fwb      = 0         !  FreshWater Budget: =0 unchecked
      !                    !     =1 global mean of e-p-r set to zero at each time step
      !                    !     =2 annual global mean of e-p-r set to zero
   ln_rnf      = .true.    !  runoffs                                   (T => fill namsbc_rnf)
   ln_apr_dyn  = .true.    !  Patm gradient added in ocean & ice Eqs.   (T => fill namsbc_apr )
   ln_wave     = .false.   !  Activate coupling with wave  (T => fill namsbc_wave)
   nn_lsm      = 0         !  =0 land/sea mask for input fields is not applied (keep empty land/sea mask filename field) ,
                           !  =1:n number of iterations of land/sea mask application for input fields (fill land/sea mask filename field)
/
!-----------------------------------------------------------------------
&namsbc_flx    !   surface boundary condition : flux formulation        (ln_flx =T)
!-----------------------------------------------------------------------
/
!-----------------------------------------------------------------------
&namsbc_blk    !   namsbc_blk  generic Bulk formula          (ln_blk =T)
!-----------------------------------------------------------------------
   !                    !  bulk algorithm :
   ln_NCAR      = .false.    ! "NCAR"      algorithm   (Large and Yeager 2008)
   ln_COARE_3p0 = .false.    ! "COARE 3.0" algorithm   (Fairall et al. 2003)
   ln_COARE_3p6 = .false.    ! "COARE 3.6" algorithm   (Edson et al. 2013)
   ln_ECMWF     = .true.     ! "ECMWF"     algorithm   (IFS cycle 45r1)
   ln_ANDREAS   = .false.    ! "ANDREAS"   algorithm   (Andreas et al. 2015)
      rn_zqt       =  2.     !  Air temperature & humidity reference height (m)
      rn_zu        = 10.     !  Wind vector reference height (m)
      nn_iter_algo =  5      !  Number of iterations in bulk param. algo ("stable ABL + weak wind" requires more)
      ln_skin_cs   = .false. !  use the cool-skin parameterization  => use at least nn_iter_algo > 10
      ln_skin_wl   = .false. !  use the warm-layer parameterization => use at least nn_iter_algo > 10
   !
   ln_humi_sph = .false. !  humidity "sn_humi" is specific humidity  [kg/kg]
   ln_humi_dpt = .true.  !  humidity "sn_humi" is dew-point temperature [K]
   ln_humi_rlh = .false. !  humidity "sn_humi" is relative humidity     [%]
   ln_tair_pot = .false. !  air temperature read in "sn_tair" is already POTENTIAL TEMPERATURE, NOT ABSOLUTE (ECMWF => ln_tair_pot=.false.)
   !!   
   !
   cn_dir      = './FATM/' !  root directory for the bulk data location
   !___________!_________________________!___________________!___________!_____________!________!___________!______________________________________!__________!_______________!
   !           !  file name              ! frequency (hours) ! variable  ! time interp.!  clim  ! 'yearly'/ !       weights filename               ! rotation ! land/sea mask !
   !           !                         !  (if <0  months)  !   name    !   (logical) !  (T/F) ! 'monthly' !                                      ! pairing  !    filename   !
   sn_wndi =  'u10_ERA5_tropico-box_surface_1h', 1., 'u10',     .true.   , .false. , 'yearly'  , 'weight_bicubic_ERA5_tropico-box_605x341-TROPICO12.nc'  , 'Uwnd'   , ''
   sn_wndj =  'v10_ERA5_tropico-box_surface_1h', 1., 'v10',     .true.   , .false. , 'yearly'  , 'weight_bicubic_ERA5_tropico-box_605x341-TROPICO12.nc'  , 'Vwnd'   , ''
   sn_qsr  = 'ssrd_ERA5_tropico-box_surface_1h', 1.,'ssrd',     .false.  , .false. , 'yearly'  , 'weight_bilinear_ERA5_tropico-box_605x341-TROPICO12.nc' , ''       , ''
   sn_qlw  = 'strd_ERA5_tropico-box_surface_1h', 1.,'strd',     .false.  , .false. , 'yearly'  , 'weight_bilinear_ERA5_tropico-box_605x341-TROPICO12.nc' , ''       , ''
   sn_tair =  't2m_ERA5_tropico-box_surface_1h', 1., 't2m',     .true.   , .false. , 'yearly'  , 'weight_bilinear_ERA5_tropico-box_605x341-TROPICO12.nc' , ''       , ''
   sn_humi =  'd2m_ERA5_tropico-box_surface_1h', 1., 'd2m',     .true.   , .false. , 'yearly'  , 'weight_bilinear_ERA5_tropico-box_605x341-TROPICO12.nc' , ''       , ''
   sn_prec =   'tp_ERA5_tropico-box_surface_1h', 1.,  'tp',     .false.  , .false. , 'yearly'  , 'weight_bilinear_ERA5_tropico-box_605x341-TROPICO12.nc' , ''       , ''
   sn_snow =   'tp_ERA5_tropico-box_surface_1h', 1.,  'tp',     .false.  , .false. , 'yearly'  , 'weight_bilinear_ERA5_tropico-box_605x341-TROPICO12.nc' , ''       , ''
   sn_slp  =  'msl_ERA5_tropico-box_surface_1h', 1., 'msl',     .true.   , .false. , 'yearly'  , 'weight_bilinear_ERA5_tropico-box_605x341-TROPICO12.nc' , ''       , ''
/
!-----------------------------------------------------------------------
&namsbc_abl    !   Atmospheric Boundary Layer formulation           (ln_abl = T)
!-----------------------------------------------------------------------
/
!-----------------------------------------------------------------------
&namsbc_sas    !   Stand-Alone Surface module: ocean data               (SAS_SRC  only)
!-----------------------------------------------------------------------
/
!-----------------------------------------------------------------------
&namsbc_iif    !   Ice-IF : use observed ice cover                      (nn_ice = 1)
!-----------------------------------------------------------------------
/
!-----------------------------------------------------------------------
&namtra_qsr    !   penetrative solar radiation                          (ln_traqsr =T)
!-----------------------------------------------------------------------
   !                       ! type of penetration                        (default: NO selection)
   ln_qsr_rgb  = .true.       !#lolo  RGB light penetration (Red-Green-Blue)
   !
   nn_chldta   =      1       !  RGB : Chl data (=1) or cst value (=0)

   cn_dir      = './'      !  root directory for the chlorophyl data location
   !___________!_________________________!___________________!___________!_____________!________!___________!__________________!__________!_______________!
   !           !  file name              ! frequency (hours) ! variable  ! time interp.!  clim  ! 'yearly'/ ! weights filename ! rotation ! land/sea mask !
   !           !                         !  (if <0  months)  !   name    !   (logical) !  (T/F) ! 'monthly' !                  ! pairing  !    filename   !
   sn_chl      = 'chla_ESACCI-OC-TROPICO12_clim',        -1.     , 'chla'    ,   .TRUE.    , .TRUE. , 'yearly'  , ''               , ''       , ''
/
!-----------------------------------------------------------------------
&namsbc_ssr    !   surface boundary condition : sea surface restoring   (ln_ssr =T)
!-----------------------------------------------------------------------
/
!-----------------------------------------------------------------------
&namsbc_rnf    !   runoffs                                              (ln_rnf =T)
!-----------------------------------------------------------------------
   ln_rnf_mouth = .true.    !  specific treatment at rivers mouths
      rn_hrnf     =  15.e0    !  depth over which enhanced vertical mixing is used    (ln_rnf_mouth=T)
      rn_avt_rnf  =   2.e-3   !  value of the additional vertical mixing coef. [m2/s] (ln_rnf_mouth=T)
   rn_rfact    =   1.e0    !  multiplicative factor for runoff

   cn_dir      = './'      !  root directory for the runoff data location
   !___________!_________________________!___________________!___________!_____________!________!___________!__________________!__________!_______________!
   !           !  file name              ! frequency (hours) ! variable  ! time interp.!  clim  ! 'yearly'/ ! weights filename ! rotation ! land/sea mask !
   !           !                         !  (if <0  months)  !   name    !   (logical) !  (T/F) ! 'monthly' !                  ! pairing  !    filename   !
   sn_rnf = 'runoff_v2.4',                 -1. , 'sorunoff', .true.  , .true. , 'yearly'  , ''       , ''       , ''
   sn_cnf = 'runoff_v2.4',                  0. , 'socoefr', .false.  , .true. , 'yearly'  , ''       , ''       , ''
   sn_s_rnf    = 'runoffs'            ,    24.   , 'rosaline',   .true.     , .true. , 'yearly'  , ''       , ''       , ''
   sn_t_rnf    = 'runoffs'            ,    24.   , 'rotemper',   .true.     , .true. , 'yearly'  , ''       , ''       , ''
   sn_dep_rnf  = 'runoffs'            ,     0.   , 'rodepth' ,   .false.    , .true. , 'yearly'  , ''       , ''       , ''
/
!-----------------------------------------------------------------------
&namsbc_apr    !   Atmospheric pressure used as ocean forcing           (ln_apr_dyn =T)
!-----------------------------------------------------------------------
   rn_pref     = 101000.   !  reference atmospheric pressure   [N/m2]/
   ln_ref_apr  = .false.   !  ref. pressure: global mean Patm (T) or a constant (F)
   ln_apr_obc  = .false.   !  inverse barometer added to OBC ssh data

   cn_dir = './FATM/'        !  root directory for the Patm data location
   !___________!_________________________!___________________!___________!_____________!________!___________!__________________!______$
   !           !  file name              ! frequency (hours) ! variable  ! time interp.!  clim  ! 'yearly'/ ! weights filename ! rotat$
   !           !                         !  (if <0  months)  !   name    !   (logical) !  (T/F) ! 'monthly' !                  ! pairi$
   sn_apr  =  'msl_ERA5_tropico-box_surface_1h', 1., 'msl',     .true.   , .false. , 'yearly'  , 'weight_bilinear_ERA5_tropico-box_605x341-TROPICO12.nc' , ''       , ''
/
!-----------------------------------------------------------------------
&namisf       !  Top boundary layer (ISF)                               (default: OFF)
!-----------------------------------------------------------------------
/
!-----------------------------------------------------------------------
&namsbc_wave   ! External fields from wave model                        (ln_wave=T)
!-----------------------------------------------------------------------
/
!-----------------------------------------------------------------------
&namberg       !   iceberg parameters                                   (default: OFF)
!-----------------------------------------------------------------------
/
!-----------------------------------------------------------------------
&namsbc_fwb    !   freshwater-budget adjustment                         (nn_fwb > 0)
!-----------------------------------------------------------------------
/

!!======================================================================
!!               ***  Lateral boundary condition  ***                 !!
!!                                                                    !!
!!   namlbc        lateral momentum boundary condition                  (default: NO selection)
!!   namagrif      agrif nested grid   (read by child model only)       ("key_agrif")
!!   nam_tide      Tidal forcing                                        (default: OFF)
!!   nambdy        Unstructured open boundaries                         (default: OFF)
!!   nambdy_dta    Unstructured open boundaries - external data         (see  nambdy)
!!   nambdy_tide   tidal forcing at open boundaries                     (default: OFF)
!!======================================================================
!
!-----------------------------------------------------------------------
&namlbc        !   lateral momentum boundary condition                  (default: NO selection)
!-----------------------------------------------------------------------
   !                       !  free slip  !   partial slip  !   no slip   ! strong slip
   rn_shlat    =  0.       !  shlat = 0  !  0 < shlat < 2  !  shlat = 2  !  2 < shlat
   ln_vorlat   = .false.   !#lolo?  consistency of vorticity boundary condition with analytical Eqs.
/
!-----------------------------------------------------------------------
&namagrif      !  AGRIF zoom                                            ("key_agrif")
!-----------------------------------------------------------------------
/
!-----------------------------------------------------------------------
&nam_tide      !   tide parameters                                      (default: OFF)
!-----------------------------------------------------------------------
   ln_tide     = .true.       ! Activate tides
      nn_tide_var   = 1          !  Variant of tidal parameter set and tide-potential computation
      !                          !     (1: default; 0: compatibility with previous versions)
      ln_tide_dia   = .true.     !  Enable tidal diagnostic output
      ln_tide_pot   = .true.                !  use tidal potential forcing
         rn_tide_gamma = 0.7                   ! Tidal tilt factor
         ln_scal_load  = .true.                ! Use scalar approximation for
            rn_scal_load = 0.094               !     load potential
         ln_read_load  = .false.               ! Or read load potential from file
            cn_tide_load = 'tide_LOAD_grid_T.nc'  ! filename for load potential
            !
      ln_tide_ramp  = .false.               !  Use linear ramp for tides at startup
         rn_tide_ramp_dt = 0.               !  ramp duration in days
      sn_tide_cnames(1) = 'M2'              !  name of constituent - all tidal components must be set in namelist_cfg
      sn_tide_cnames(2) = 'S2'
      sn_tide_cnames(3) = 'N2'
      sn_tide_cnames(4) = 'K1'
      sn_tide_cnames(5) = 'O1'
/
!-----------------------------------------------------------------------
&nambdy        !  unstructured open boundaries                          (default: OFF)
!-----------------------------------------------------------------------
   ln_bdy         = .true.    !  Use unstructured open boundaries
   nb_bdy         = 5         !  number of open boundary sets
   !
   ln_coords_file = .false.,.false.,.false.,.false.,.false.   !  =T : read bdy coordinates from file
      cn_coords_file = 'none','none','none','none','none'
   !
   ln_mask_file   = .false.   !  =T : read mask from file
      cn_mask_file   = 'none'
   !
   cn_dyn2d    = 'flather','flather','flather','flather','flather'
   nn_dyn2d_dta   =  3,3,3,3,3        !  = 0, bdy data are equal to the initial state
   !                          !  = 1, bdy data are read in 'bdydata   .nc' files
   !                          !  = 2, use tidal harmonic forcing data from files
   !                          !  = 3, use external data AND tidal harmonic forcing
   cn_dyn3d      =  'frs','frs','frs','frs','frs'     !
   nn_dyn3d_dta  =  1,1,1,1,1         !  = 0, bdy data are equal to the initial state
   !                          !  = 1, bdy data are read in 'bdydata   .nc' files
   cn_tra        =  'frs','frs','frs','frs','frs'     !
   nn_tra_dta    =  1,1,1,1,1         !  = 0, bdy data are equal to the initial state
   !                          !  = 1, bdy data are read in 'bdydata   .nc' files
   !
   ln_tra_dmp    =.false.,.false.,.false.,.false.,.false.     !  open boudaries conditions for tracers
   ln_dyn3d_dmp  =.false.,.false.,.false.,.false.,.false.     !  open boundary condition for baroclinic velocities
   rn_time_dmp   =  1.,1.,1.,1.,1.        !  Damping time scale in days
   rn_time_dmp_out = 1.,1.,1.,1.,1.       !  Outflow damping time scale
   !                          !  = 1, bdy data are read in 'bdydata   .nc' files
   cn_ice        =  'none','none','none','none','none'
   nn_ice_dta    =  0,0,0,0,0         !  = 0, bdy data are equal to the initial state
   !                          !  = 1, bdy data are read in 'bdydata   .nc' files
   nn_rimwidth   =  10,10,10,10,10   !  width of the relaxation zone
   ln_vol        = .false.    !  total volume correction (see nn_volctl parameter)
   nn_volctl     =  1         !  = 0, the total water flux across open boundaries is zero
/
!
!-----------------------------------------------------------------------
&nambdy_index  !  structured open boundaries definition
!-----------------------------------------------------------------------
   ctypebdy ='S'      !! S1        ! Open boundary type (W,E,S or N)
   nbdyind  =  2                   ! indice of velocity row or column
   nbdybeg  = 332                  ! indice of segment start
   nbdyend  = 1716                 ! indice of segment end
/
!-----------------------------------------------------------------------
&nambdy_index  !  structured open boundaries definition
!-----------------------------------------------------------------------
   ctypebdy ='S'      !! S2        ! Open boundary type (W,E,S or N)
   nbdyind  = 57                   ! indice of velocity row or column
   nbdybeg  = 69                   ! indice of segment start
   nbdyend  = 349                  ! indice of segment end
/
!-----------------------------------------------------------------------
&nambdy_index  !  structured open boundaries definition
!!-----------------------------------------------------------------------
   ctypebdy ='W'      !! W1        ! Open boundary type (W,E,S or N)
   nbdyind  = 58                   ! indice of velocity row or column
   nbdybeg  = 91                   ! indice of segment start
   nbdyend  = 139                  ! indice of segment end
/
!-----------------------------------------------------------------------
&nambdy_index  !  structured open boundaries definition
!!-----------------------------------------------------------------------
   ctypebdy ='W'      !! W2        ! Open boundary type (W,E,S or N)
   nbdyind  =  2                   ! indice of velocity row or column
   nbdybeg  = 487                  ! indice of segment start
   nbdyend  = 925                  ! indice of segment end
/
!-----------------------------------------------------------------------
&nambdy_index  !  structured open boundaries definition
!-----------------------------------------------------------------------
   ctypebdy ='N'      !! N1        ! Open boundary type (W,E,S or N)
   nbdyind  = 924                  ! indice of velocity row or column
   nbdybeg  =   2                  ! indice of segment start
   nbdyend  = 1342                 ! indice of segment end
/
!!
!-----------------------------------------------------------------------
&nambdy_dta    ! S1 ! open boundaries - external data                       (see nam_bdy)
!-----------------------------------------------------------------------
   ln_zinterp  = .false.
   ln_full_vel = .true.
   cn_dir      = './BDY/'
   bn_ssh = 'sossheig_GLORYS2V4-BDY_t_S1_TROPICO12_L125_tr21', -1., 'sossheig',    .true.   , .false.,  'yearly'  ,    ''            ,   ''     ,     ''
   bn_u2d = 'none'                          ,         24.    , 'vobtcrtx',    .true.   , .false.,  'yearly'  ,    ''            ,   ''     ,     ''
   bn_v2d = 'none'                          ,         24.    , 'vobtcrty',    .true.   , .false.,  'yearly'  ,    ''            ,   ''     ,     ''
   bn_u3d = 'vozocrtx_GLORYS2V4-BDY_u_S1_TROPICO12_L125_tr21', -1., 'vozocrtx',    .true.   , .false.,  'yearly'  ,    ''            ,   ''     ,     ''
   bn_v3d = 'vomecrty_GLORYS2V4-BDY_v_S1_TROPICO12_L125_tr21', -1., 'vomecrty',    .true.   , .false.,  'yearly'  ,    ''            ,   ''     ,     ''
   bn_tem = 'votemper_GLORYS2V4-BDY_t_S1_TROPICO12_L125_tr21', -1., 'votemper',    .true.   , .false.,  'yearly'  ,    ''            ,   ''     ,     ''
   bn_sal = 'vosaline_GLORYS2V4-BDY_t_S1_TROPICO12_L125_tr21', -1., 'vosaline',    .true.   , .false.,  'yearly'  ,    ''            ,   ''     ,     ''
/
!-----------------------------------------------------------------------
&nambdy_dta    ! S2 ! open boundaries - external data                       (see nam_bdy)
!-----------------------------------------------------------------------
   ln_zinterp  = .false.
   ln_full_vel = .true.
   cn_dir      = './BDY/'
   bn_ssh = 'sossheig_GLORYS2V4-BDY_t_S2_TROPICO12_L125_tr21', -1., 'sossheig',    .true.   , .false.,  'yearly'  ,    ''            ,   ''     ,     ''
   bn_u2d = 'none'                          ,         24.    , 'vobtcrtx',    .true.   , .false.,  'yearly'  ,    ''            ,   ''     ,     ''
   bn_v2d = 'none'                          ,         24.    , 'vobtcrty',    .true.   , .false.,  'yearly'  ,    ''            ,   ''     ,     ''
   bn_u3d = 'vozocrtx_GLORYS2V4-BDY_u_S2_TROPICO12_L125_tr21', -1., 'vozocrtx',    .true.   , .false.,  'yearly'  ,    ''            ,   ''     ,     ''
   bn_v3d = 'vomecrty_GLORYS2V4-BDY_v_S2_TROPICO12_L125_tr21', -1., 'vomecrty',    .true.   , .false.,  'yearly'  ,    ''            ,   ''     ,     ''
   bn_tem = 'votemper_GLORYS2V4-BDY_t_S2_TROPICO12_L125_tr21', -1., 'votemper',    .true.   , .false.,  'yearly'  ,    ''            ,   ''     ,     ''
   bn_sal = 'vosaline_GLORYS2V4-BDY_t_S2_TROPICO12_L125_tr21', -1., 'vosaline',    .true.   , .false.,  'yearly'  ,    ''            ,   ''     ,     ''
/
!-----------------------------------------------------------------------
&nambdy_dta    ! W1 ! open boundaries - external data                       (see nam_bdy)
!-----------------------------------------------------------------------
   ln_zinterp  = .false.
   ln_full_vel = .true.
   cn_dir      = './BDY/'
   bn_ssh = 'sossheig_GLORYS2V4-BDY_t_W1_TROPICO12_L125_tr21', -1., 'sossheig',    .true.   , .false.,  'yearly'  ,    ''            ,   ''     ,     ''
   bn_u2d = 'none'                          ,         24.    , 'vobtcrtx',    .true.   , .false.,  'yearly'  ,    ''            ,   ''     ,     ''
   bn_v2d = 'none'                          ,         24.    , 'vobtcrty',    .true.   , .false.,  'yearly'  ,    ''            ,   ''     ,     ''
   bn_u3d = 'vozocrtx_GLORYS2V4-BDY_u_W1_TROPICO12_L125_tr21', -1., 'vozocrtx',    .true.   , .false.,  'yearly'  ,    ''            ,   ''     ,     ''
   bn_v3d = 'vomecrty_GLORYS2V4-BDY_v_W1_TROPICO12_L125_tr21', -1., 'vomecrty',    .true.   , .false.,  'yearly'  ,    ''            ,   ''     ,     ''
   bn_tem = 'votemper_GLORYS2V4-BDY_t_W1_TROPICO12_L125_tr21', -1., 'votemper',    .true.   , .false.,  'yearly'  ,    ''            ,   ''     ,     ''
   bn_sal = 'vosaline_GLORYS2V4-BDY_t_W1_TROPICO12_L125_tr21', -1., 'vosaline',    .true.   , .false.,  'yearly'  ,    ''            ,   ''     ,     ''
/
!-----------------------------------------------------------------------
&nambdy_dta    ! W2 ! open boundaries - external data                       (see nam_bdy)
!-----------------------------------------------------------------------
   ln_zinterp  = .false.
   ln_full_vel = .true.
   cn_dir      = './BDY/'
   bn_ssh = 'sossheig_GLORYS2V4-BDY_t_W2_TROPICO12_L125_tr21', -1., 'sossheig',    .true.   , .false.,  'yearly'  ,    ''            ,   ''     ,     ''
   bn_u2d = 'none'                          ,         24.    , 'vobtcrtx',    .true.   , .false.,  'yearly'  ,    ''            ,   ''     ,     ''
   bn_v2d = 'none'                          ,         24.    , 'vobtcrty',    .true.   , .false.,  'yearly'  ,    ''            ,   ''     ,     ''
   bn_u3d = 'vozocrtx_GLORYS2V4-BDY_u_W2_TROPICO12_L125_tr21', -1., 'vozocrtx',    .true.   , .false.,  'yearly'  ,    ''            ,   ''     ,     ''
   bn_v3d = 'vomecrty_GLORYS2V4-BDY_v_W2_TROPICO12_L125_tr21', -1., 'vomecrty',    .true.   , .false.,  'yearly'  ,    ''            ,   ''     ,     ''
   bn_tem = 'votemper_GLORYS2V4-BDY_t_W2_TROPICO12_L125_tr21', -1., 'votemper',    .true.   , .false.,  'yearly'  ,    ''            ,   ''     ,     ''
   bn_sal = 'vosaline_GLORYS2V4-BDY_t_W2_TROPICO12_L125_tr21', -1., 'vosaline',    .true.   , .false.,  'yearly'  ,    ''            ,   ''     ,     ''
/
!-----------------------------------------------------------------------
&nambdy_dta    ! N1 ! open boundaries - external data                       (see nam_bdy)
!-----------------------------------------------------------------------
   ln_zinterp  = .false.
   ln_full_vel = .true.
   cn_dir      = './BDY/'
   bn_ssh = 'sossheig_GLORYS2V4-BDY_t_N1_TROPICO12_L125_tr21', -1., 'sossheig',    .true.   , .false.,  'yearly'  ,    ''            ,   ''     ,     ''
   bn_u2d = 'none'                          ,         24.    , 'vobtcrtx',    .true.   , .false.,  'yearly'  ,    ''            ,   ''     ,     ''
   bn_v2d = 'none'                          ,         24.    , 'vobtcrty',    .true.   , .false.,  'yearly'  ,    ''            ,   ''     ,     ''
   bn_u3d = 'vozocrtx_GLORYS2V4-BDY_u_N1_TROPICO12_L125_tr21', -1., 'vozocrtx',    .true.   , .false.,  'yearly'  ,    ''            ,   ''     ,     ''
   bn_v3d = 'vomecrty_GLORYS2V4-BDY_v_N1_TROPICO12_L125_tr21', -1., 'vomecrty',    .true.   , .false.,  'yearly'  ,    ''            ,   ''     ,     ''
   bn_tem = 'votemper_GLORYS2V4-BDY_t_N1_TROPICO12_L125_tr21', -1., 'votemper',    .true.   , .false.,  'yearly'  ,    ''            ,   ''     ,     ''
   bn_sal = 'vosaline_GLORYS2V4-BDY_t_N1_TROPICO12_L125_tr21', -1., 'vosaline',    .true.   , .false.,  'yearly'  ,    ''            ,   ''     ,     ''
/
!-----------------------------------------------------------------------
&nambdy_tide ! S1  !  tidal forcing at open boundaries                      (default: OFF)
!-----------------------------------------------------------------------
   filtide          = './BDY/bdytide_FES-TROPICO12_S1_'   !  file name root of tidal forcing files
   ln_bdytide_2ddta = .false.                   !
/
!-----------------------------------------------------------------------
&nambdy_tide ! S2  !  tidal forcing at open boundaries                      (default: OFF)
!-----------------------------------------------------------------------
   filtide          = './BDY/bdytide_FES-TROPICO12_S2_'   !  file name root of tidal forcing files
   ln_bdytide_2ddta = .false.                   !
/
!-----------------------------------------------------------------------
&nambdy_tide ! W1  !  tidal forcing at open boundaries                      (default: OFF)
!-----------------------------------------------------------------------
   filtide          = './BDY/bdytide_FES-TROPICO12_W1_'   !  file name root of tidal forcing files
   ln_bdytide_2ddta = .false.                   !
/
!-----------------------------------------------------------------------
&nambdy_tide ! W2  !  tidal forcing at open boundaries                      (default: OFF)
!-----------------------------------------------------------------------
   filtide          = './BDY/bdytide_FES-TROPICO12_W2_'   !  file name root of tidal forcing files
   ln_bdytide_2ddta = .false.                   !
/
!-----------------------------------------------------------------------
&nambdy_tide ! N1  !  tidal forcing at open boundaries                      (default: OFF)
!-----------------------------------------------------------------------
   filtide          = './BDY/bdytide_FES-TROPICO12_N1_'   !  file name root of tidal forcing files
   ln_bdytide_2ddta = .false.                   !
/

!!======================================================================
!!                ***  Top/Bottom boundary condition  ***             !!
!!                                                                    !!
!!   namdrg        top/bottom drag coefficient                          (default: NO selection)
!!   namdrg_top    top    friction                                      (ln_drg_OFF=F & ln_isfcav=T)
!!   namdrg_bot    bottom friction                                      (ln_drg_OFF=F)
!!   nambbc        bottom temperature boundary condition                (default: OFF)
!!   nambbl        bottom boundary layer scheme                         (default: OFF)
!!======================================================================
!
!-----------------------------------------------------------------------
&namdrg        !   top/bottom drag coefficient                          (default: NO selection)
!-----------------------------------------------------------------------
   ln_drg_OFF  = .false.   !  free-slip       : Cd = 0                  (F => fill namdrg_bot
   ln_lin      = .false.   !      linear  drag: Cd = Cd0 Uc0                   &   namdrg_top)
   ln_non_lin  = .false.   !#lolo  non-linear  drag: Cd = Cd0 |U|
   ln_loglayer = .true.   !#lolo  logarithmic drag: Cd = vkarmn/log(z/z0) |U|
   !
   ln_drgimp   = .true.    !  implicit top/bottom friction flag
      ln_drgice_imp = .false. ! implicit ice-ocean drag
/
!-----------------------------------------------------------------------
&namdrg_top    !   TOP friction                                         (ln_drg_OFF =F & ln_isfcav=T)
!-----------------------------------------------------------------------
/
!-----------------------------------------------------------------------
&namdrg_bot    !   BOTTOM friction                                      (ln_drg_OFF =F)
!-----------------------------------------------------------------------
   rn_Cd0      =  3.e-3   !#lolo  drag coefficient [-]
   rn_Uc0      =  0.4      !  ref. velocity [m/s] (linear drag=Cd0*Uc0)
   rn_Cdmax    =  0.1      !  drag value maximum [-] (logarithmic drag)
   rn_ke0      =  2.5e-3       !#lolo: MUST BE 0 because tidal motion in???  background kinetic energy  [m2/s2] (non-linear cases)
   !!   rn_z0       =  3.e-3    !#lolo: Seems WAY too small!!!  roughness [m] (ln_loglayer=T)
   rn_z0       =  1.e-2    !#lolo: Better!!!  roughness [m] (ln_loglayer=T)
   ln_boost    = .false.   !#lolo  =T regional boost of Cd0 ; =F constant
      rn_boost =  50.         !  local boost factor  [-]
/
!-----------------------------------------------------------------------
&nambbc        !   bottom temperature boundary condition                (default: OFF)
!-----------------------------------------------------------------------
/
!-----------------------------------------------------------------------
&nambbl        !   bottom boundary layer scheme                         (default: OFF)
!-----------------------------------------------------------------------
   ln_trabbl   = .false.      !#lolo  Bottom Boundary Layer parameterisation flag
/

!!======================================================================
!!                        Tracer (T-S) namelists                      !!
!!                                                                    !!
!!   nameos        equation of state                                    (default: NO selection)
!!   namtra_adv    advection scheme                                     (default: NO selection)
!!   namtra_ldf    lateral diffusion scheme                             (default: NO selection)
!!   namtra_mle    mixed layer eddy param. (Fox-Kemper param.)          (default: OFF)
!!   namtra_eiv    eddy induced velocity param.                         (default: OFF)
!!   namtra_dmp    T & S newtonian damping                              (default: OFF)
!!======================================================================
!
!-----------------------------------------------------------------------
&nameos        !   ocean Equation Of Seawater                           (default: NO selection)
!-----------------------------------------------------------------------
   ln_teos10   = .true.          !  = Use TEOS-10
/
!-----------------------------------------------------------------------
&namtra_adv    !   advection scheme for tracer                          (default: NO selection)
!-----------------------------------------------------------------------
   ln_traadv_fct = .true.  !  FCT scheme
      nn_fct_h   =  4            !#lolo  =2/4, horizontal 2nd / 4th order 
      nn_fct_v   =  4            !#lolo  =2/4, vertical   2nd / COMPACT 4th order 
/
!-----------------------------------------------------------------------
&namtra_ldf    !   lateral diffusion scheme for tracers                 (default: NO selection)
!-----------------------------------------------------------------------
   !                       !  Operator type:
   ln_traldf_lap   =  .true.   !    laplacian operator
   !                       !  Direction of action:
   ln_traldf_lev   = .false.   !  iso-level
   ln_traldf_hor   = .false.   !#lolo  horizontal  (geopotential)
   ln_traldf_iso   = .true.   !#lolo  iso-neutral (standard operator)
   ln_traldf_triad = .false.   !  iso-neutral (triad    operator)
   !
   !                       !  Coefficients:
   nn_aht_ijk_t    = 20        !  space/time variation of eddy coefficient:
      !                             !   =-20 (=-30)    read in eddy_diffusivity_2D.nc (..._3D.nc) file
      !                             !   =  0           constant 
      !                             !   = 10 F(k)      =ldf_c1d 
      !                             !   = 20 F(i,j)    =ldf_c2d 
      !                             !   = 21 F(i,j,t)  =Treguier et al. JPO 1997 formulation
      !                             !   = 30 F(i,j,k)  =ldf_c2d * ldf_c1d
      !                             !   = 31 F(i,j,k,t)=F(local velocity and grid-spacing)
      !                        !  time invariant coefficients:  aht0 = 1/2  Ud*Ld   (lap case) 
      !                             !                           or   = 1/12 Ud*Ld^3 (blp case)
      rn_Ud        = 0.02           !#lolo  lateral diffusive velocity [m/s] (nn_aht_ijk_t= 0, 10, 20, 30)
/
!-----------------------------------------------------------------------
&namtra_mle    !   mixed layer eddy parametrisation (Fox-Kemper)       (default: OFF)
!-----------------------------------------------------------------------
   ln_mle      = .false.  !#lolo (T) use the Mixed Layer Eddy (MLE) parameterisation
/
!-----------------------------------------------------------------------
&namtra_eiv    !   eddy induced velocity param.                         (default: OFF)
!-----------------------------------------------------------------------
   ln_ldfeiv   = .false.   ! use eddy induced velocity parameterization
/
!-----------------------------------------------------------------------
&namtra_dmp    !   tracer: T & S newtonian damping                      (default: OFF)
!-----------------------------------------------------------------------
   ln_tradmp   =  .false.  !  add a damping term (using resto.nc coef.)
/

!!======================================================================
!!                      ***  Dynamics namelists  ***                  !!
!!                                                                    !!
!!   nam_vvl       vertical coordinate options                          (default: z-star)
!!   namdyn_adv    formulation of the momentum advection                (default: NO selection)
!!   namdyn_vor    advection scheme                                     (default: NO selection)
!!   namdyn_hpg    hydrostatic pressure gradient                        (default: NO selection)
!!   namdyn_spg    surface pressure gradient                            (default: NO selection)
!!   namdyn_ldf    lateral diffusion scheme                             (default: NO selection)
!!   namdta_dyn    offline TOP: dynamics read in files                  (OFF_SRC only)
!!======================================================================
!
!-----------------------------------------------------------------------
&nam_vvl       !   vertical coordinate options                          (default: z-star)
!-----------------------------------------------------------------------
   ln_vvl_zstar  = .true.           !  z-star vertical coordinate
/
!-----------------------------------------------------------------------
&namdyn_adv    !   formulation of the momentum advection                (default: NO selection)
!-----------------------------------------------------------------------
   ln_dynadv_ubs = .true.   !#lolo!  flux form - 3rd order UBS      scheme
/
!-----------------------------------------------------------------------
&namdyn_vor    !   Vorticity / Coriolis scheme                          (default: NO selection)
!-----------------------------------------------------------------------
   ln_dynvor_een = .true.  !  energy & enstrophy scheme
   nn_e3f_typ = 0          !  type of e3f (EEN, ENE, ENS, MIX only)  =0  e3f = mi(mj(e3t))/4
   !                       !                                         =1  e3f = mi(mj(e3t))/mi(mj( tmask))
/
!-----------------------------------------------------------------------
&namdyn_hpg    !   Hydrostatic pressure gradient option                 (default: NO selection)
!-----------------------------------------------------------------------
   ln_hpg_sco  = .true.   !  s-coordinate (standard jacobian formulation)
/
!-----------------------------------------------------------------------
&namdyn_spg    !   surface pressure gradient                            (default: NO selection)
!-----------------------------------------------------------------------
   ln_dynspg_exp  = .false.   ! explicit free surface
   ln_dynspg_ts   = .true.    ! split-explicit free surface
      ln_bt_fw      = .true.     !#lolo Forward integration of barotropic Eqs.
      ln_bt_av      = .true.     ! Time filtering of barotropic variables
         nn_bt_flt     = 2          !#lolo Time filter choice  = 0 None
         !                          !                     = 1 Boxcar over   nn_e sub-steps
         !                          !                     = 2 Boxcar over 2*nn_e  "    "
      ln_bt_auto    = .true.     ! Number of sub-step defined from:
         rn_bt_cmax   =  0.7        ! =T : the Maximum Courant Number allowed
         nn_e         = 30          ! =F : the number of sub-step in rn_Dt seconds
      rn_bt_alpha   = 0.         ! Temporal diffusion parameter (if ln_bt_av=F)
/
!-----------------------------------------------------------------------
&namdyn_ldf    !   lateral diffusion on momentum                        (default: NO selection)
!-----------------------------------------------------------------------
   ln_dynldf_OFF = .true.     !#lolo: because UBS for dyn!  No operator (i.e. no explicit diffusion)
/

!!======================================================================
!!                     vertical physics namelists                     !!
!!                                                                    !!
!!    namzdf        vertical physics manager                            (default: NO selection)
!!    namzdf_ric    richardson number vertical mixing                   (ln_zdfric=T)
!!    namzdf_tke    TKE vertical mixing                                 (ln_zdftke=T)
!!    namzdf_gls    GLS vertical mixing                                 (ln_zdfgls=T)
!!    namzdf_osm    OSM vertical diffusion                              (ln_zdfosm=T)
!!    namzdf_iwm    tidal mixing parameterization                       (ln_zdfiwm=T)
!!======================================================================
!
!-----------------------------------------------------------------------
&namzdf        !   vertical physics manager                             (default: NO selection)
!-----------------------------------------------------------------------
   !                       ! adaptive-implicit vertical advection
   ln_zad_Aimp = .true.      !#lolo  Courant number dependent scheme (Shchepetkin 2015)
   !
   !                       ! type of vertical closure (required)
   ln_zdfcst   = .false.      !  constant mixing
   ln_zdftke   = .true.       !  Turbulent Kinetic Energy closure       (T =>   fill namzdf_tke)
   !
   !                       ! convection
   ln_zdfevd   = .true.       !  Enhanced Vertical Diffusion scheme
      nn_evdm  =    0            !  evd apply on tracer (=0) or on tracer and momentum (=1)
      rn_evd   =  100.           !  evd mixing coefficient [m2/s]
   !                       !  Coefficients
!!JMM eORCA12.L75-GJM2020: because param de Casimir !!!
   rn_avm0     =   1.e-5      !#lolo?  vertical eddy viscosity   [m2/s]       (background Kz if ln_zdfcst=F)
   rn_avt0     =   1.e-6      !#lolo?  vertical eddy diffusivity [m2/s]       (background Kz if ln_zdfcst=F)
   nn_avb      =    0         !  profile for background avt & avm (=1) or not (=0)
   nn_havtb    =    1         !  horizontal shape for avtb (=1) or not (=0)
/
!-----------------------------------------------------------------------
&namzdf_ric    !   richardson number dependent vertical diffusion       (ln_zdfric =T)
!-----------------------------------------------------------------------
/
!-----------------------------------------------------------------------
&namzdf_tke    !   turbulent eddy kinetic dependent vertical diffusion  (ln_zdftke =T)
!-----------------------------------------------------------------------
   rn_ediff    =   0.1     !  coef. for vertical eddy coef. (avt=rn_ediff*mxl*sqrt(e) )
   rn_ediss    =   0.7     !  coef. of the Kolmogoroff dissipation
   rn_ebb      =  67.83    !  coef. of the surface input of tke (=67.83 suggested when ln_mxl0=T)
   rn_emin     =   1.e-6   !  minimum value of tke [m2/s2]
   rn_emin0    =   1.e-4   !  surface minimum value of tke [m2/s2]
   rn_bshear   =   1.e-20  ! background shear (>0) currently a numerical threshold (do not change it)
   nn_pdl      =   1       !  Prandtl number function of richarson number (=1, avt=pdl(Ri)*avm) or not (=0, avt=avm)
   nn_mxl      =   3       !  mixing length: = 0 bounded by the distance to surface and bottom
   !                       !                 = 1 bounded by the local vertical scale factor
   !                       !                 = 2 first vertical derivative of mixing length bounded by 1
   !                       !                 = 3 as =2 with distinct dissipative an mixing length scale
   ln_mxl0     = .true.    !  surface mixing length scale = F(wind stress) (T) or not (F)
      nn_mxlice    = 2        ! type of scaling under sea-ice
      !                       !    = 0 no scaling under sea-ice
      !                       !    = 1 scaling with constant sea-ice thickness
      !                       !    = 2 scaling with mean sea-ice thickness ( only with SI3 sea-ice model )
      !                       !    = 3 scaling with maximum sea-ice thickness
      rn_mxlice   = 10.       ! max constant ice thickness value when scaling under sea-ice ( nn_mxlice=1)
   rn_mxl0     =   0.04    !  surface  buoyancy lenght scale minimum value
   ln_mxhsw    = .false.   !  surface mixing length scale = F(wave height)
   ln_lc       = .true.    !  Langmuir cell parameterisation (Axell 2002)
      rn_lc       =   0.15    !  coef. associated to Langmuir cells
   nn_etau     =   1       !  penetration of tke below the mixed layer (ML) due to NIWs
   !                          !        = 0 none ; = 1 add a tke source below the ML
   !                          !        = 2 add a tke source just at the base of the ML
   !                          !        = 3 as = 1 applied on HF part of the stress           (ln_cpl=T)
      rn_efr      =   0.05    !  fraction of surface tke value which penetrates below the ML (nn_etau=1 or 2)
      nn_htau     =   1       !  type of exponential decrease of tke penetration below the ML
      !                       !        = 0  constant 10 m length scale
      !                       !        = 1  0.5m at the equator to 30m poleward of 40 degrees
   nn_eice     =   1       !  attenutaion of langmuir & surface wave breaking under ice
   !                       !           = 0 no impact of ice cover on langmuir & surface wave breaking
   !                       !           = 1 weigthed by 1-TANH(10*fr_i)
   !                       !           = 2 weighted by 1-fr_i
   !                       !           = 3 weighted by 1-MIN(1,4*fr_i)   
   nn_bc_surf   =     1    !  surface condition (0/1=Dir/Neum) ! Only applicable for wave coupling (ln_cplwave=1)
   nn_bc_bot    =     1    !  bottom condition (0/1=Dir/Neum) ! Only applicable for wave coupling (ln_cplwave=1)
/
!-----------------------------------------------------------------------
&namzdf_gls    !   GLS vertical diffusion                               (ln_zdfgls =T)
!-----------------------------------------------------------------------
/
!-----------------------------------------------------------------------
&namzdf_osm    !   OSM vertical diffusion                               (ln_zdfosm =T)
!-----------------------------------------------------------------------
/
!-----------------------------------------------------------------------
&namosm_mle    !   mixed layer eddy parametrisation (Fox-Kemper)       (default: OFF)
!-----------------------------------------------------------------------
/
!-----------------------------------------------------------------------
&namzdf_mfc     !   Mass Flux Convection
!-----------------------------------------------------------------------
/
!-----------------------------------------------------------------------
&namzdf_iwm    !    internal wave-driven mixing parameterization        (ln_zdfiwm =T)
!-----------------------------------------------------------------------
/

!!======================================================================
!!                  ***  Diagnostics namelists  ***                   !!
!!                                                                    !!
!!   namtrd       dynamics and/or tracer trends                         (default: OFF)
!!   namhsb       Heat and salt budgets                                 (default: OFF)
!!   namdiu       Cool skin and warm layer models                       (default: OFF)
!!   namdiu       Cool skin and warm layer models                       (default: OFF)
!!   namflo       float parameters                                      (default: OFF)
!!   nam_diadct   transports through some sections                      (default: OFF)
!!   nam_dia25h   25h Mean Output                                       (default: OFF)
!!   namnc4       netcdf4 chunking and compression settings             ("key_netcdf4")
!!======================================================================
!
!-----------------------------------------------------------------------
&namtrd        !   trend diagnostics                                    (default: OFF)
!-----------------------------------------------------------------------
/
!-----------------------------------------------------------------------
&namhsb        !  Heat and salt budgets                                 (default: OFF)
!-----------------------------------------------------------------------
/
!-----------------------------------------------------------------------
&namdiu        !   Cool skin and warm layer models                      (default: OFF)
!-----------------------------------------------------------------------
/
!-----------------------------------------------------------------------
&namflo        !   float parameters                                     (default: OFF)
!-----------------------------------------------------------------------
/
!-----------------------------------------------------------------------
&nam_diadct    !   transports through some sections                     (default: OFF)
!-----------------------------------------------------------------------
/
!-----------------------------------------------------------------------
&nam_dia25h    !  25h Mean Output                                       (default: OFF)
!-----------------------------------------------------------------------
/
!-----------------------------------------------------------------------
&namnc4        !   netcdf4 chunking and compression settings            ("key_netcdf4")
!-----------------------------------------------------------------------
/

!!======================================================================
!!               ***  Observation & Assimilation  ***                 !!
!!                                                                    !!
!!   namobs       observation and model comparison                      (default: OFF)
!!   nam_asminc   assimilation increments                               ('key_asminc')
!!======================================================================
!
!-----------------------------------------------------------------------
&namobs        !  observation usage switch                              (default: OFF)
!-----------------------------------------------------------------------
/
!-----------------------------------------------------------------------
&nam_asminc    !   assimilation increments                              ('key_asminc')
!-----------------------------------------------------------------------
/

!!======================================================================
!!                  ***  Miscellaneous namelists  ***                 !!
!!                                                                    !!
!!   nammpp            Massively Parallel Processing
!!   namctl            Control prints                                   (default: OFF)
!!   namsto            Stochastic parametrization of EOS                (default: OFF)
!!======================================================================
!
!-----------------------------------------------------------------------
&nammpp        !   Massively Parallel Processing
!-----------------------------------------------------------------------
   ln_listonly =  .false.   !  do nothing else than listing the best domain decompositions (with land domains suppression)
   !                        !  if T: the largest number of cores tested is defined by max(mppsize, jpni*jpnj)
   ln_nnogather =  .false.
   jpni        =   0       !  number of processors following i (set automatically if < 1), see also ln_listonly = T
   jpnj        =   0       !  number of processors following j (set automatically if < 1), see also ln_listonly = T
   nn_hls      =   1       !  halo width (applies to both rows and columns)
/
!-----------------------------------------------------------------------
&namctl        !   Control prints                                       (default: OFF)
!-----------------------------------------------------------------------
sn_cfctl%l_oceout  = .false.
/
!-----------------------------------------------------------------------
&namsto        ! Stochastic parametrization of EOS                      (default: OFF)
!-----------------------------------------------------------------------
/
