#undef DEBUG
#include <misc.h>

module cldcond

!---------------------------------------------------------------------------------
! Purpose:
!
! Provides the CAM interface to the prognostic cloud water and ice parameterization
!
! Author: Byron Boville  Sept 04, 2002
!
!---------------------------------------------------------------------------------

  use shr_kind_mod,  only: r8=>shr_kind_r8
  use ppgrid,        only: pcols, pver, pverp
  use physconst,     only: gravit, latvap, latice
  use abortutils,    only: endrun
  use perf_mod,      only: t_startf, t_stopf  !_EXTERNAL
!  use ref_pres,       only: top_lev=>trop_cloud_top_lev 
!++sxj
  use pmgrid,        only: masterproc,iam    
!--sxj
  implicit none
  private
  save

  public :: cldcond_register, cldcond_init_cnst, cldcond_implements_cnst
  public :: cldcond_init, cldcond_tend
  public :: cldcond_zmconv_detrain
  public :: cldcond_sediment
!++sxj
#ifdef MG08
  public ::  stratiform_tend
  !public :: cldcond_init, stratiform_tend
!#else
!  public :: cldcond_init, cldcond_tend
!  public :: cldcond_zmconv_detrain
!  public :: cldcond_sediment
#endif
!--sxj

! Private module data

#ifdef SPP_DEEPCONV
  !++wy
  integer    :: sppnum1_idx = 0
  integer    :: cape1_idx   = 0
  !--wy
#endif

#ifdef MICROP-MG
  integer, parameter :: ncnst=4                      ! number of constituents
  character(len=8), dimension(ncnst), parameter :: & ! constituent names
     cnst_names = (/'CLDLIQ', 'CLDICE','NUMLIQ','NUMICE'/)

  integer :: &
       ixcldice,     &! cloud ice water index
       ixcldliq,     &! cloud liquid water index
       ixnumliq     ,&! cloud liquid number index
       ixnumice     ,&! cloud ice water index
!	   ixzzzice  	  ! cloud ice ZZZ clu
  ! Physics buffer indices
  integer  ::  qcwat_idx
  integer  ::  lcwat_idx
  integer  ::  iccwat_idx
  integer  ::  nlwat_idx
  integer  ::  niwat_idx
  integer  ::  cc_t_idx
  integer  ::  cc_qv_idx
  integer  ::  cc_ql_idx
  integer  ::  cc_qi_idx
  integer  ::  cc_nl_idx
  integer  ::  cc_ni_idx
  integer  ::  cc_qlst_idx
  integer  ::  tcwat_idx
  integer  ::  cld_idx
  integer  ::  cldo_idx
  integer  ::  ast_idx
  integer  ::  aist_idx
  integer  ::  alst_idx
  integer  ::  qist_idx
  integer  ::  qlst_idx
  integer  ::  concld_idx
  integer  ::  rhdfda_idx
  integer  ::  rhu00_idx
  integer  ::  rel2_idx
  integer  ::  rei2_idx
  integer  ::  concldql_idx
  integer  ::  fice_idx
  integer  ::  sh_frac_idx
  integer  ::  dp_frac_idx
  integer  ::  qini_idx
  integer  ::  cldliqini_idx
  integer  ::  cldiceini_idx
  integer  ::  tini_idx
  integer  ::  qme_idx
  integer  ::  prain_idx
  integer  ::  nevapr_idx
  integer  ::  wsedl_idx
  integer  ::  rei_idx
  integer  ::  rel_idx
  integer  ::  rel_fn_idx
  integer  ::  dei_idx
  integer  ::  mu_idx
  integer  ::  lambdac_idx
  integer  ::  iciwp_idx
  integer  ::  iclwp_idx
  integer  ::  deiconv_idx
  integer  ::  muconv_idx
  integer  ::  lambdaconv_idx
  integer  ::  iciwpst_idx
  integer  ::  iclwpst_idx
  integer  ::  iciwpconv_idx
  integer  ::  iclwpconv_idx
  integer  ::  des_idx
  integer  ::  icswp_idx
  integer  ::  cldfsnow_idx
  integer  ::  ls_mrprc_idx
  integer  ::  ls_mrsnw_idx
!---
! shallow convection
!
   integer    ::     icwmrsh_idx    = 0
   integer    ::      rprdsh_idx    = 0
   integer    ::     rprdtot_idx    = 0
   integer    ::      cldtop_idx    = 0
   integer    ::      cldbot_idx    = 0
   integer    ::        cush_idx    = 0
   integer    :: nevapr_shcu_idx    = 0
   integer    ::       shfrc_idx    = 0

!---
! deep convection
!
   integer    ::  icwmrdp_idx    = 0
   integer    ::  rprddp_idx     = 0
   integer    ::  nevapr_dpcu_idx= 0
!
! Vertical diffusion
!
  integer :: tke_idx   ! TKE and eddy diffusivity indices for fields in the physics buffer
  integer :: kvh_idx
  integer :: kvm_idx
  integer :: turbtype_idx
  integer :: smaw_idx
  integer :: tauresx_idx
  integer :: tauresy_idx
  integer :: wgustd_index
!
! conv_water
!
  integer :: sh_cldliq1_idx
  integer :: sh_cldice1_idx

#else
!-----------------------------------------------------------

#ifdef MG08
  integer, parameter :: ncnst=4                     ! number of constituents
  character(len=8), dimension(ncnst), parameter :: & ! constituent names
  cnst_names = (/'CLDLIQ', 'CLDICE','NUMLIQ','NUMICE'/)
  integer ::       &
     ixnumliq,     & ! cloud liquid number index
     ixnumice,     & ! cloud ice water index
!	 ixzzzice,     & ! cloud ice ZZZ index clu
     kvh_idx
#else
  integer, parameter :: ncnst=2                      ! number of constituents
  character(len=8), dimension(ncnst), parameter :: & ! constituent names
  cnst_names = (/'CLDLIQ', 'CLDICE'/)
#endif

  integer :: &
       ixcldice,     &! cloud ice water index
       ixcldliq,     &! cloud liquid water index
       qcwat_idx,    &! qcwat index in physics buffer
       tcwat_idx,    &! tcwat index in physics buffer
       cld_idx,      &! cld index in physics buffer
       lcwat_idx      ! lcwat index in physics buffer
#endif

#ifdef BYANG20
   integer    ::       dt_shcu_idx     ! T tendency due to shallow convection
   integer    ::       dq_shcu_idx     ! Q Tendency due to shallow convection
   integer    ::       org_dpcu_idx    ! Memory of convection
   integer    ::       prec_dpcu_idx   ! Deep convective precipitation rate
#endif

contains

!===============================================================================
  subroutine cldcond_sediment(state, ptend, dtime, &
                             ls_reffrain, ls_reffsnow, ls_flxprc,  ls_flxsnw, & 
#ifdef UWMT
       cloud, icefrac, landfrac, ocnfrac, prec, snow, landm, snowh, prsed, wsedl)
#else
       cloud, icefrac, landfrac, ocnfrac, prec, snow, landm, snowh, prsed)
#endif
!
! Interface to sedimentation of cloud liquid and ice particles
!
! NOTE: initial implementation by B.A. Boville from earlier code by P.J. Rasch
!
! B. A. Boville, Sept 20, 2002
!
!-----------------------------------------------------------------------
    use physics_types,    only: physics_state, physics_ptend
    use pkg_cld_sediment, only: cld_sediment_vel, cld_sediment_tend
    use history,          only: outfld

! Arguments
    type(physics_state), intent(in)    :: state   ! state variables
    type(physics_ptend), intent(inout) :: ptend   ! package tendencies

    real(r8), intent(in)  :: cloud(pcols,pver)    ! cloud fraction
    real(r8), intent(in)  :: icefrac (pcols)      ! sea ice fraction (fraction)
    real(r8), intent(in)  :: landfrac(pcols)      ! land fraction (fraction)
    real(r8), intent(in)  :: ocnfrac (pcols)      ! ocean fraction (fraction)
    real(r8), intent(in)  :: dtime                ! timestep

    real(r8), intent(out) :: prec(pcols)          ! surface flux of total cloud water
    real(r8), intent(out) :: snow(pcols)          ! surface flux of cloud ice
    real(r8), intent(in) :: landm(pcols)          ! land fraction ramped over water
    real(r8), intent(in) :: snowh(pcols)         ! Snow depth over land, water equivalent (m)
    real(r8), optional, intent(out) :: prsed(pcols,pver)          ! flux of total cloud water

!------
! wtw

   real(r8), intent(out) :: ls_flxprc(pcols,pverp)    ! stratiform interface gbm flux_cloud_rain+snow (kg m^-2 s^-1) ^M
   real(r8), intent(out) :: ls_flxsnw(pcols,pverp)    ! stratiform interface gbm flux_cloud_snow (kg m^-2 s^-1)

   real(r8), intent(out) :: ls_reffrain(pcols,pver)    ! rain effective drop radius (microns)
   real(r8), intent(out) :: ls_reffsnow(pcols,pver)    ! snow effective drop size (microns)

   real(r8)  :: ls_flxrain(pcols,pverp) 

! Local variables
    integer  :: i,k                               ! loop indexes
    integer  :: ncol                              ! number of atmospheric columns in chunk
    integer  :: lchnk                             ! chunk index
    real(r8) :: rain(pcols)                       ! surface flux of cloud liquid
    real(r8) :: pvliq(pcols,pver+1)               ! vertical velocity of cloud liquid drops (Pa/s)
    real(r8) :: pvice(pcols,pver+1)               ! vertical velocity of cloud ice particles (Pa/s)
#ifdef UWMT
    real(r8), intent(out) :: wsedl(pcols,pver)    ! Sedimentation velocity of liquid stratus cloud droplet [ m/s ]
#endif
!-----------------------------------------------------------------------
    ncol = state%ncol

    ptend%name         = 'pcwsediment'
    ptend%ls           = .TRUE.
    ptend%lq(1)        = .TRUE.
    ptend%lq(ixcldice) = .TRUE.
    ptend%lq(ixcldliq) = .TRUE.

    call cld_sediment_vel (ncol,                                    &
         icefrac, landfrac, ocnfrac, state%pmid, state%pdel, state%t, &
         cloud, state%q(:,:,ixcldliq), state%q(:,:,ixcldice), pvliq, pvice, landm, snowh , &
         ls_reffrain, ls_reffsnow)

#ifdef UWMT
         wsedl(:ncol,:pver) = pvliq(:ncol,:pver)/gravit/(state%pmid(:ncol,:pver)/(287.15_r8*state%t(:ncol,:pver)))
#endif

    call cld_sediment_tend (ncol, dtime ,                                       &
         state%pint, state%pmid           , state%pdel            , state%t     , &
         cloud  ,    state%q(:,:,ixcldliq), state%q(:,:,ixcldice) , pvliq, pvice, &
         ptend%q(:,:,ixcldliq), ptend%q(:,:,ixcldice), ptend%q(:,:,1), ptend%s  , &
         rain   , snow,  prsed  , &
         ls_flxrain, ls_flxsnw )                   ! wtw

   ls_flxprc(:,:) = ls_flxrain(:,:) + ls_flxsnw(:,:)   ! wtw 

! convert rain and snow from kg/m2 to m/s
    snow(:ncol) = snow(:ncol)*1.0e-3
    rain(:ncol) = rain(:ncol)*1.0e-3
! compute total precip (m/s)
    prec(:ncol) = rain(:ncol) + snow(:ncol)

! record history variables
    lchnk = state%lchnk
    call outfld('DQSED'   ,ptend%q(:,:,1)       , pcols,lchnk)
    call outfld('DISED'   ,ptend%q(:,:,ixcldice), pcols,lchnk)
    call outfld('DLSED'   ,ptend%q(:,:,ixcldliq), pcols,lchnk)
    call outfld('HSED'    ,ptend%s              , pcols,lchnk)
    call outfld('PRECSED' ,prec                 , pcols,lchnk)
    call outfld('SNOWSED' ,snow                 , pcols,lchnk)
    call outfld('RAINSED' ,rain                 , pcols,lchnk)

    return
  end subroutine cldcond_sediment

!===============================================================================
  subroutine cldcond_zmconv_detrain(dlf, cld, state, ptend)
!
! Partition the detrained condensed water from the ZM convection scheme.
!
! The ZM scheme does not have an ice phase, so the detrained water is paritioned
! soley between cloud liquid water and the environment. The ice/liquid partitioning 
! happens in cldcond_tend.
!
! NOTE: initial implementation by B.A. Boville just moves code here from TPHYSBC.
!
! B. A. Boville, Sept 09, 2002
!
!-----------------------------------------------------------------------
    use physics_types, only: physics_state, physics_ptend
    use history,       only: outfld

    implicit none

! Arguments
    type(physics_state), intent(in  )  :: state   ! state variables
    type(physics_ptend), intent(inout) :: ptend   ! package tendencies

    real(r8), intent(in) :: dlf(pcols,pver)       ! detrained water from ZM
    real(r8), intent(in) :: cld(pcols,pver)       ! cloud fraction

! Local variables
    integer :: i,k                            ! loop indexes
!-----------------------------------------------------------------------

    ptend%name         = 'pcwdetrain'
!!$    ptend%ls           = .TRUE.
!!$    ptend%lq(1)        = .TRUE.
!!$    ptend%lq(ixcldice) = .TRUE.
    ptend%lq(ixcldliq) = .TRUE.
!
! Put all of the detraining cloud water from convection into the large scale cloud.
! It all goes in liquid for the moment.
    do k = 1,pver
       do i = 1,state%ncol
!!$          ptend%q(i,k,1)        = dlf(i,k) * (1.-cld(i,k))
!!$          ptend%s(i,k)          =-dlf(i,k) * (1.-cld(i,k))*latvap
!!$          ptend%q(i,k,ixcldice) = 0.
!!$          ptend%q(i,k,ixcldliq) = dlf(i,k) * cld(i,k)
          ptend%q(i,k,ixcldliq) = dlf(i,k)
       end do
    end do
    call outfld('ZMDLF' ,dlf  , pcols,state%lchnk)

    return
  end subroutine cldcond_zmconv_detrain

!===============================================================================
  subroutine cldcond_register()
!
! Register the constituents (cloud liquid and cloud ice) and the fields
! in the physics buffer.
! 
!-----------------------------------------------------------------------
    use constituents, only: cnst_add, advected, nonadvec, ppcnst
    use physconst,    only: mwdry, cpair
    use phys_buffer,  only: pbuf_times, pbuf_add

    implicit none

!    logical, parameter :: cldw_adv=.false.  ! true => cloud water is treated as advected tracer
    logical, parameter :: cldw_adv=.true.  ! true => cloud water is treated as advected tracer

    integer flag
!-----------------------------------------------------------------------
#ifdef SPP_DEEPCONV
!++wy
   call pbuf_add('SPPNUM1', 'global', 1, 1, 1, sppnum1_idx)
   call pbuf_add('CAPE1',   'global', 1, 1, 1, cape1_idx)
!--wy
#endif

!-------------
#ifdef MICROP-MG

! Register cloud water and determine index (either advected or non-adv).
    if (cldw_adv) then
       flag = advected
    else 
       flag = nonadvec
    endif
    call cnst_add(cnst_names(1), flag, mwdry, cpair, 0._r8, ixcldliq, &
         longname='Grid box averaged liquid condensate amount')
    call cnst_add(cnst_names(2), flag, mwdry, cpair, 0._r8, ixcldice, &
         longname='Grid box averaged ice condensate amount')
    
    call cnst_add(cnst_names(3), flag, mwdry, cpair, 0._r8, ixnumliq, &
          longname='Grid box averaged cloud liquid number')
    call cnst_add(cnst_names(4), flag, mwdry, cpair, 0._r8, ixnumice, &
          longname='Grid box averaged cloud ice number')
!++clu
!	call cnst_add(cnst_names(5), flag, mwdry, cpair, 0._r8, ixzzzice, &
!          longname='Grid box averaged cloud ice ZZZ')
!--clu
  ! Request physics buffer space for fields that persist across timesteps.

    call pbuf_add('QCWAT',   'global',  1, pver, pbuf_times,   qcwat_idx)
    call pbuf_add('LCWAT',   'global',  1, pver, pbuf_times,   lcwat_idx)
    call pbuf_add('ICCWAT',  'global',  1, pver, pbuf_times,  iccwat_idx)
    call pbuf_add('NLWAT',   'global',  1, pver, pbuf_times,   nlwat_idx)
    call pbuf_add('NIWAT',   'global',  1, pver, pbuf_times,   niwat_idx)
    call pbuf_add('CC_T',    'global',  1, pver, pbuf_times,    cc_t_idx)
    call pbuf_add('CC_qv',   'global',  1, pver, pbuf_times,   cc_qv_idx)
    call pbuf_add('CC_ql',   'global',  1, pver, pbuf_times,   cc_ql_idx)
    call pbuf_add('CC_qi',   'global',  1, pver, pbuf_times,   cc_qi_idx)
    call pbuf_add('CC_nl',   'global',  1, pver, pbuf_times,   cc_nl_idx)
    call pbuf_add('CC_ni',   'global',  1, pver, pbuf_times,   cc_ni_idx)
    call pbuf_add('CC_qlst', 'global',  1, pver, pbuf_times, cc_qlst_idx)
    call pbuf_add('TCWAT',   'global',  1, pver, pbuf_times,   tcwat_idx)
    call pbuf_add('CLD',     'global',  1, pver, pbuf_times,     cld_idx)
    call pbuf_add('CLDO',    'global',  1, pver, pbuf_times,    cldo_idx)
    call pbuf_add('AST',     'global',  1, pver, pbuf_times,     ast_idx)
    call pbuf_add('AIST',    'global',  1, pver, pbuf_times,    aist_idx)
    call pbuf_add('ALST',    'global',  1, pver, pbuf_times,    alst_idx)
    call pbuf_add('QIST',    'global',  1, pver, pbuf_times,    qist_idx)
    call pbuf_add('QLST',    'global',  1, pver, pbuf_times,    qlst_idx)
    call pbuf_add('CONCLD',  'global',  1, pver, pbuf_times,  concld_idx)
    call pbuf_add('RHDFDA',  'global',  1, pver, pbuf_times,  rhdfda_idx)
    call pbuf_add('RHU00',   'global',  1, pver, pbuf_times,   rhu00_idx)
    call pbuf_add('REL2',    'global',  1, pver, pbuf_times,    rel2_idx)
    call pbuf_add('REI2',    'global',  1, pver, pbuf_times,    rei2_idx)

  ! Physics buffer variables for convective cloud properties.

    call pbuf_add('CONCLDQL',   'physpkg', 1, pver, 1, concldql_idx)
    call pbuf_add('FICE',       'physpkg', 1, pver, 1, fice_idx)
    call pbuf_add('SH_FRAC',    'physpkg', 1, pver, 1, sh_frac_idx)
    call pbuf_add('DP_FRAC',    'physpkg', 1, pver, 1, dp_frac_idx)

    call pbuf_add('QINI'      , 'physpkg', 1,pver, 1, qini_idx)
    call pbuf_add('CLDLIQINI' , 'physpkg', 1,pver, 1, cldliqini_idx)
    call pbuf_add('CLDICEINI' , 'physpkg', 1,pver, 1, cldiceini_idx)
    call pbuf_add('TINI'      , 'physpkg', 1,pver, 1, tini_idx)

    call pbuf_add('QME',        'physpkg', 1, pver, 1, qme_idx)
    call pbuf_add('PRAIN' ,     'physpkg', 1, pver, 1, prain_idx)
    call pbuf_add('NEVAPR' ,    'physpkg', 1, pver, 1, nevapr_idx)

    call pbuf_add('WSEDL',      'physpkg', 1, pver, 1, wsedl_idx)

    call pbuf_add('REI',        'physpkg', 1, pver, 1, rei_idx)
    call pbuf_add('REL',        'physpkg', 1, pver, 1, rel_idx)

    call pbuf_add('REL_FN',     'physpkg', 1, pver, 1, rel_fn_idx)          ! REL at fixed number for indirect rad forcing

    call pbuf_add('DEI',        'physpkg', 1, pver, 1, dei_idx)          ! Mitchell ice effective diameter for radiation
    call pbuf_add('MU',         'physpkg', 1, pver, 1, mu_idx)          ! Size distribution shape parameter for radiation
    call pbuf_add('LAMBDAC',    'physpkg', 1, pver, 1, lambdac_idx)          ! Size distribution shape parameter for radiation
    call pbuf_add('ICIWP',      'physpkg', 1, pver, 1, iciwp_idx)          ! In cloud ice water path for radiation
    call pbuf_add('ICLWP',      'physpkg', 1, pver, 1, iclwp_idx)          ! In cloud liquid water path for radiation

    call pbuf_add('DEICONV',    'physpkg', 1, pver, 1, deiconv_idx)          ! Convective ice effective diameter for radiation
    call pbuf_add('MUCONV',     'physpkg', 1, pver, 1, muconv_idx)          ! Convective size distribution shape parameter for radiation
    call pbuf_add('LAMBDACONV', 'physpkg', 1, pver, 1, lambdaconv_idx)          ! Convective size distribution shape parameter for radiation
    call pbuf_add('ICIWPST',    'physpkg', 1, pver, 1, iciwpst_idx)          ! Stratiform only in cloud ice water path for radiation
    call pbuf_add('ICLWPST',    'physpkg', 1, pver, 1, iclwpst_idx)          ! Stratiform in cloud liquid water path for radiation
    call pbuf_add('ICIWPCONV',  'physpkg', 1, pver, 1, iciwpconv_idx)          ! Convective only in cloud ice water path for radiation
    call pbuf_add('ICLWPCONV',  'physpkg', 1, pver, 1, iclwpconv_idx)          ! Convective in cloud liquid water path for radiation

    call pbuf_add('DES',        'physpkg', 1, pver, 1, des_idx)          ! Snow effective diameter for radiation
    call pbuf_add('ICSWP',      'physpkg', 1, pver, 1, icswp_idx)          ! In cloud snow water path for radiation
    call pbuf_add('CLDFSNOW',   'physpkg', 1, pver ,pbuf_times, cldfsnow_idx) ! Cloud fraction for liquid drops + snow

    call pbuf_add('LS_MRPRC',    'physpkg', 1, pver, 1, ls_mrprc_idx)
    call pbuf_add('LS_MRSNW',    'physpkg', 1, pver, 1, ls_mrsnw_idx)
!
! shallow convection
!
    call pbuf_add( 'ICWMRSH'  , 'physpkg' ,  1,  pver,  1,       icwmrsh_idx )
    call pbuf_add( 'RPRDSH'   , 'physpkg' ,  1,  pver,  1,        rprdsh_idx )
    call pbuf_add( 'RPRDTOT'  , 'physpkg' ,  1,  pver,  1,       rprdtot_idx )
    call pbuf_add( 'CLDTOP'   , 'physpkg' ,  1,  1,     1,        cldtop_idx )
    call pbuf_add( 'CLDBOT'   , 'physpkg' ,  1,  1,     1,        cldbot_idx )
    call pbuf_add( 'cush'     , 'global'  ,  1,  1,     pbuf_times, cush_idx )
    call pbuf_add( 'NEVAPR_SHCU', 'physpkg', 1,  pver,  1,   nevapr_shcu_idx )

    call pbuf_add( 'shfrc'  ,  'physpkg' ,  1,  pver,  1,  shfrc_idx )
!
! deep convection
!
    call pbuf_add('ICWMRDP' , 'physpkg', 1,pver,      1,         icwmrdp_idx)
    call pbuf_add('RPRDDP' , 'physpkg', 1,pver,       1,          rprddp_idx)
    call pbuf_add('NEVAPR_DPCU' , 'physpkg', 1,pver,      1, nevapr_dpcu_idx)

!
!  Vertical diffusion
!

    call pbuf_add( 'TKE',      'global',  1,  pverp,  pbuf_times,  tke_idx )
    call pbuf_add( 'KVH',      'global',  1,  pverp,  pbuf_times,  kvh_idx )
    call pbuf_add( 'KVM',      'global',  1,  pverp,  pbuf_times,  kvm_idx )
    call pbuf_add( 'TURBTYPE', 'global',  1,  pverp,  pbuf_times,  turbtype_idx )
    call pbuf_add( 'SMAW',     'global',  1,  pverp,  pbuf_times,  smaw_idx )
    call pbuf_add( 'TAURESX',  'global',  1,  1,      pbuf_times,  tauresx_idx )
    call pbuf_add( 'TAURESY',  'global',  1,  1,      pbuf_times,  tauresy_idx )
    call pbuf_add( 'WGUSTD',   'global',  1,  1,      1,           wgustd_index )
!-------
! conv_water_register
!
! these calls were already done in convect_shallow...so here I add the same fields to the physics buffer with a "1" at the end

    call pbuf_add('SH_CLDLIQ1', 'physpkg', 1, pver,  1, sh_cldliq1_idx)  ! shallow gbm cloud liquid water (kg/kg)
    call pbuf_add('SH_CLDICE1', 'physpkg', 1, pver,  1, sh_cldice1_idx)  ! shallow gbm cloud ice water (kg/kg)

!----------
! used in original version
!
    call pbuf_add('FRACS',   'global', 1,pver,ppcnst   , flag     ) ! wzz, 2008.8.28
    call pbuf_add('CSTCLD' , 'physpkg', 1,pver,      1, flag)   ! wzz, 2008.2.28
    call pbuf_add('CONPRE' , 'physpkg', 1,pver,      1, flag)   ! wzz, 2008.2.28
    call pbuf_add('CSTPRE' , 'physpkg', 1,pver,      1, flag)   ! wzz, 2008.2.28

#else
!-------------------------------------------
! Register cloud water and determine index (either advected or non-adv).
    if (cldw_adv) then
       flag = advected
    else
       flag = nonadvec
    endif
    call cnst_add(cnst_names(1), flag, mwdry, cpair, 0._r8, ixcldliq, &
         longname='Grid box averaged liquid condensate amount')
    call cnst_add(cnst_names(2), flag, mwdry, cpair, 0._r8, ixcldice, &
         longname='Grid box averaged ice condensate amount')

!++sxj
    if ( masterproc) write(6,*) "SXJ cldcond.F90 cldcond_register ncnst",ncnst
#ifdef MG08
    call cnst_add(cnst_names(3), flag, mwdry, cpair, 0._r8, ixnumliq, &
          longname='Grid box averaged cloud liquid number')
    call cnst_add(cnst_names(4), flag, mwdry, cpair, 0._r8, ixnumice, &
          longname='Grid box averaged cloud ice number')
!++clu
!	call cnst_add(cnst_names(5), flag, mwdry, cpair, 0._r8, ixzzzice, &
!         longname='Grid box averaged cloud ice ZZZ')
!--clu
!    call pbuf_add('REI',    'physpkg', 1, pver, 1, flag)
!    call pbuf_add('REL',    'physpkg', 1, pver, 1, flag)
    call pbuf_add('REL_FN', 'physpkg', 1, pver, 1, flag) ! REL at fixed number for indirect rad forcing
    call pbuf_add( 'KVH',   'global',  1, pverp,  pbuf_times,  kvh_idx )
    if ( masterproc) write(6,*) "SXJ cldcond.F90 cldcond_register ixnumliq,ixnumice",ixnumliq,ixnumice
#endif
!--sxj
    call pbuf_add('REI',    'physpkg', 1, pver, 1, flag)
    call pbuf_add('REL',    'physpkg', 1, pver, 1, flag)

! Request physics buffer space for fields that persist across timesteps.
    call pbuf_add('QCWAT', 'global', 1,pver,pbuf_times, qcwat_idx)
    call pbuf_add('TCWAT', 'global', 1,pver,pbuf_times, tcwat_idx)
    call pbuf_add('CLD',   'global', 1,pver,pbuf_times, cld_idx)
    call pbuf_add('LCWAT', 'global', 1,pver,pbuf_times, lcwat_idx)
    call pbuf_add('FRACS', 'global', 1,pver,ppcnst    , flag     ) ! wzz, 2008.8.28

    call pbuf_add('QINI' , 'physpkg', 1,pver,      1, flag)
    call pbuf_add('TINI' , 'physpkg', 1,pver,      1, flag)
    call pbuf_add('CONCLD' , 'physpkg', 1,pver,      1, flag)   ! wzz, 2008.2.28
    call pbuf_add('CSTCLD' , 'physpkg', 1,pver,      1, flag)   ! wzz, 2008.2.28
    call pbuf_add('CONPRE' , 'physpkg', 1,pver,      1, flag)   ! wzz, 2008.2.28
    call pbuf_add('CSTPRE' , 'physpkg', 1,pver,      1, flag)   ! wzz, 2008.2.28

#ifdef UWMT
    call pbuf_add('WSEDL', 'physpkg', 1, pver, 1, flag)
    call pbuf_add('QRL'  , 'global',  1, pver, 1, flag) ! longwave  radiative heating rate
    call pbuf_add('LHFLX', 'global',  1, 1, pbuf_times, flag)
    call pbuf_add('SHFLX', 'global',  1, 1, pbuf_times, flag)
    call pbuf_add('TAUX', 'global',  1, 1, pbuf_times, flag)
    call pbuf_add('TAUY', 'global',  1, 1, pbuf_times, flag)
    call pbuf_add('QFLX', 'global',  1, 1, pbuf_times, flag)
    call pbuf_add('LHFLX_RES', 'global',  1, 1, pbuf_times, flag)
    call pbuf_add('SHFLX_RES', 'global',  1, 1, pbuf_times, flag)
    call pbuf_add('TAUX_RES', 'global',  1, 1, pbuf_times, flag)
    call pbuf_add('TAUY_RES', 'global',  1, 1, pbuf_times, flag)
    call pbuf_add('QFLX_RES', 'global',  1, 1, pbuf_times, flag)
#endif

#ifdef BYANG20
    call pbuf_add( 'DT_SHCU',      'global',  1,  pver,  1,  dt_shcu_idx )
    call pbuf_add( 'DQ_SHCU',      'global',  1,  pver,  1,  dq_shcu_idx )
    call pbuf_add( 'ORG_DPCU',     'global',  1,     1,  1,  org_dpcu_idx )
    call pbuf_add( 'PREC_DPCU',    'global',  1,     1,  1,  prec_dpcu_idx )
#endif

#endif

!------------------------------------
  end subroutine cldcond_register

!===============================================================================

  function cldcond_implements_cnst(name)
!----------------------------------------------------------------------- 
! 
! Purpose: return true if specified constituent is implemented by this package
! 
! Author: B. Eaton
! 
!-----------------------------------------------------------------------
     implicit none
!-----------------------------Arguments---------------------------------

     character(len=*), intent(in) :: name   ! constituent name
     logical :: cldcond_implements_cnst     ! return value
!---------------------------Local workspace-----------------------------
     integer :: m
!-----------------------------------------------------------------------

     cldcond_implements_cnst = .false.
     do m = 1, ncnst
        if (name == cnst_names(m)) then
           cldcond_implements_cnst = .true.
           return
        end if
     end do
  end function cldcond_implements_cnst

!===============================================================================
  subroutine cldcond_init_cnst(name, q)
!
! Initialize the cloud water mixing ratios (liquid and ice), if they are
! not read from the initial file
! 
!-----------------------------------------------------------------------
    use pmgrid,        only: plon, plev, plat

    implicit none

! Arguments
    character(len=*), intent(in)  :: name                ! constituent name
    real(r8),         intent(out) :: q(plon,plev,plat)   ! mass mixing ratio
!-----------------------------------------------------------------------

    if ( name == 'CLDLIQ' ) then
       q = 0.0
       return
    else if ( name == 'CLDICE' ) then
       q = 0.0
       return
    end if


#ifdef MICROP-MG
    if ( name == 'NUMLIQ' ) then
       q = 0.0_r8
       return
    else if ( name == 'NUMICE' ) then
       q = 0.0_r8
       return
    endif
#endif

!++sxj
#ifdef MG08
    if ( name == 'NUMLIQ' ) then
       q = 0.0_r8
       return
    else if ( name == 'NUMICE' ) then
       q = 0.0_r8
       return
    endif
#endif
!--sxj

  end subroutine cldcond_init_cnst

!===============================================================================
  subroutine cldcond_init()
!
! Initialize the cloud water parameterization
! 
!-----------------------------------------------------------------------
    use cldwat,        only: inimc
!++sxj
#ifdef MG08
   use cldwat2m,      only: inimmc
!#else
!   use cldwat,        only: inimc
#endif
!--sxj
    use history,       only: addfld, add_default, phys_decomp
    use physconst,     only: tmelt, rh2o, rhodair

    implicit none
!-----------------------------------------------------------------------

! initialization routine for prognostic cloud water
    call inimc (tmelt, rhodair/1000.0, gravit, rh2o )
!++sxj
#ifdef MG08
    if ( masterproc) write(6,*) "SXJ cldcond.F90 cldcond_init define MG08"
    call inimmc 
!#else
!   call inimc (tmelt, rhodair/1000.0, gravit, rh2o )
#endif
!--sxj

! register history variables
    call addfld ('FWAUT   ','fraction',pver, 'A','Relative importance of liquid autoconversion' ,phys_decomp)
    call addfld ('FSAUT   ','fraction',pver, 'A','Relative importance of ice autoconversion'    ,phys_decomp)
    call addfld ('FRACW   ','fraction',pver, 'A','Relative  importance of rain accreting liquid',phys_decomp)
    call addfld ('FSACW   ','fraction',pver, 'A','Relative  importance of snow accreting liquid',phys_decomp)
    call addfld ('FSACI   ','fraction',pver, 'A','Relative  importance of snow accreting ice'   ,phys_decomp)
    call addfld ('CME     ','kg/kg/s ',pver, 'A','Rate of cond-evap within the cloud'           ,phys_decomp)
    call addfld ('ZMDLF   ','kg/kg/s ',pver, 'A','Detrained liquid water from ZM convection'    ,phys_decomp)
    call addfld ('PRODPREC','kg/kg/s ',pver, 'A','Rate of conversion of condensate to precip'   ,phys_decomp)
    call addfld ('EVAPPREC','kg/kg/s ',pver, 'A','Rate of evaporation of falling precip'        ,phys_decomp)
    call addfld ('EVAPSNOW','kg/kg/s ',pver, 'A','Rate of evaporation of falling snow'          ,phys_decomp)
    call addfld ('HPROGCLD','W/kg'    ,pver, 'A','Heating from prognostic clouds'               ,phys_decomp)
    call addfld ('HEVAP   ','W/kg'    ,pver, 'A','Heating from evaporation of falling precip'   ,phys_decomp)
    call addfld ('HMELT   ','W/kg'    ,pver, 'A','Heating from snow melt'                       ,phys_decomp)
    call addfld ('HREPART ','W/kg'    ,pver, 'A','Heating from cloud ice/liquid repartitioning' ,phys_decomp)
    call addfld ('FICE    ','fraction',pver, 'A','Fractional ice content within cloud'          ,phys_decomp)
    call addfld ('ICWMR   ','kg/kg   ',pver, 'A','Prognostic in-cloud water mixing ratio'       ,phys_decomp)
    call addfld ('ICIMR   ','kg/kg   ',pver, 'A','Prognostic in-cloud ice mixing ratio'         ,phys_decomp)
    call addfld ('PCSNOW  ','m/s'     ,1   , 'A','Snow fall from prognostic clouds'             ,phys_decomp)

    call addfld ('DQSED   ','kg/kg/s' ,pver, 'A','Water vapor tendency from cloud sedimentation',phys_decomp)
    call addfld ('DLSED   ','kg/kg/s' ,pver, 'A','Cloud liquid tendency from sedimentation'     ,phys_decomp)
    call addfld ('DISED   ','kg/kg/s' ,pver, 'A','Cloud ice tendency from sedimentation'        ,phys_decomp)
    call addfld ('HSED    ','W/kg'    ,pver, 'A','Heating from cloud sediment evaporation'      ,phys_decomp)
    call addfld ('SNOWSED ','m/s'     ,1   , 'A','Snow from cloud ice sedimentation'            ,phys_decomp)
    call addfld ('RAINSED ','m/s'     ,1   , 'A','Rain from cloud liquid sedimentation'         ,phys_decomp)
    call addfld ('PRECSED ','m/s'     ,1   , 'A','Precipitation from cloud sedimentation'       ,phys_decomp)
	
    call add_default ('FICE    ', 1, ' ')
!  These fields removed 10/30/2003 per CRB decision.
!    call add_default ('CME     ', 1, ' ')
!    call add_default ('ZMDLF   ', 1, ' ')
!    call add_default ('PRODPREC', 1, ' ')
!    call add_default ('EVAPPREC', 1, ' ')
!    call add_default ('EVAPSNOW', 1, ' ')
!    call add_default ('HPROGCLD', 1, ' ')
!    call add_default ('HEVAP   ', 1, ' ')
!    call add_default ('HMELT   ', 1, ' ')
!    call add_default ('HREPART ', 1, ' ')
!    call add_default ('ICWMR   ', 1, ' ')
!    call add_default ('ICIMR   ', 1, ' ')
!    call add_default ('PCSNOW  ', 1, ' ')
!    call add_default ('DQSED   ', 1, ' ')
!    call add_default ('DLSED   ', 1, ' ')
!    call add_default ('DISED   ', 1, ' ')
!    call add_default ('HSED    ', 1, ' ')
!    call add_default ('SNOWSED ', 1, ' ')
!    call add_default ('RAINSED ', 1, ' ')
!    call add_default ('PRECSED ', 1, ' ')


!++sxj
#ifdef MG08

   !call addfld ('HCME    ','W/kg'    ,pver, 'A','Heating from cond-evap within the cloud'      ,phys_decomp) 

    call addfld ('aer_act ','cm-3     ',pver, 'A','aerosol activate number                  '    ,phys_decomp)
    call add_default ('aer_act ', 1, ' ')
    call addfld ('aer_data  ','kg/kg',pver, 'A','read in aerosol data for test '       ,phys_decomp)
    call add_default('aer_data  ', 1, ' ')  !

    call addfld ('CMEICE  ','kg/kg/s ',pver, 'A','Rate of cond-evap of ice within the cloud'    ,phys_decomp)
    call addfld ('CMELIQ  ','kg/kg/s ',pver, 'A','Rate of cond-evap of liq within the cloud'    ,phys_decomp)
    call addfld ('ICE2PR  ','kg/kg/s ',pver, 'A','Rate of conversion of ice to precip'          ,phys_decomp)
    call addfld ('LIQ2PR  ','kg/kg/s ',pver, 'A','Rate of conversion of liq to precip'          ,phys_decomp)
    call addfld ('HFREEZ  ','W/kg'    ,pver, 'A','Heating rate due to freezing of precip'       ,phys_decomp)
    call addfld ('REPARTICE','kg/kg/s',pver, 'A','Cloud ice tendency from cloud ice/liquid repartitioning' ,phys_decomp)
    call addfld ('REPARTLIQ','kg/kg/s',pver, 'A','Cloud liq tendency from cloud ice/liquid repartitioning' ,phys_decomp)

    call addfld ('MPDT   ','K/s     ',pver, 'A','T tendency - Morrison microphysics',phys_decomp)
    call addfld ('MPDQ   ','kg/kg/s     ',pver, 'A','Q tendency - Morrison microphysics',phys_decomp)
    call addfld ('ICWNC   ','m-3     ',pver, 'A','Prognostic in-cloud water number conc'       ,phys_decomp)
    call addfld ('ICINC   ','m-3     ',pver, 'A','Prognostic in-cloud ice number conc'         ,phys_decomp)
    call addfld ('EFFLIQ  ','Micron  ',pver, 'A','Prognostic droplet effective radius'       ,phys_decomp)
    call addfld ('EFFLIQ_IND','Micron  ',pver, 'A','Prognostic droplet effective radius (indirect effect)'       ,phys_decomp)
    call addfld ('EFFICE  ','Micron  ',pver, 'A','Prognostic ice effeictive radius'         ,phys_decomp)
    call addfld ('WSUB  ','m/s     ',pver, 'A','Diagnostic sub-grid vertical velocity'         ,phys_decomp)
    call addfld ('CDNUMC  ','#/m2    ',1, 'A','Vertically-integrated droplet concentration'         ,phys_decomp)
    call add_default('ICWNC   ', 1, ' ')  
    call add_default('ICINC   ', 1, ' ')  
    call add_default('EFFLIQ  ', 1, ' ')  
    call add_default('EFFICE  ', 1, ' ')  
    call add_default('CDNUMC  ', 1, ' ')  
    call add_default('WSUB    ', 1, ' ')  
    call addfld ('CCN1    ','#/cm3   ',pver, 'A','CCN concentration at S=0.02%',phys_decomp)
    call addfld ('CCN2    ','#/cm3   ',pver, 'A','CCN concentration at S=0.05%',phys_decomp)
    call addfld ('CCN3    ','#/cm3   ',pver, 'A','CCN concentration at S=0.1%',phys_decomp)
    call add_default('CCN3    ', 1, ' ')  
    call addfld ('CCN4    ','#/cm3   ',pver, 'A','CCN concentration at S=0.2%',phys_decomp)
    call addfld ('CCN5    ','#/cm3   ',pver, 'A','CCN concentration at S=0.5%',phys_decomp)
    call addfld ('CCN6    ','#/cm3   ',pver, 'A','CCN concentration at S=1.0%',phys_decomp)
    ! diagnostic precip
    call addfld ('QRAIN   ','kg/kg   ',pver, 'A','Diagnostic grid-mean rain mixing ratio'         ,phys_decomp)
    call addfld ('QSNOW   ','kg/kg   ',pver, 'A','Diagnostic grid-mean snow mixing ratio'         ,phys_decomp)
    call addfld ('NRAIN   ','m-3     ',pver, 'A','Diagnostic grid-mean rain number conc'         ,phys_decomp)
    call addfld ('NSNOW   ','m-3     ',pver, 'A','Diagnostic grid-mean snow number conc'         ,phys_decomp)
    ! averaging for cloud particle number and size
    call addfld ('AWNC   ','m-3     ',pver, 'A','Average cloud water number conc'         ,phys_decomp)
    call addfld ('AWNI   ','m-3     ',pver, 'A','Average cloud ice number conc'         ,phys_decomp)
    call addfld ('AREL  ','Micron  ',pver, 'A','Average droplet effective radius'       ,phys_decomp)
    call addfld ('AREI  ','Micron  ',pver, 'A','Average ice effective radius'       ,phys_decomp)
    call add_default('AWNC   ', 1, ' ')  
    call add_default('AWNI   ', 1, ' ')  
    call add_default('AREL   ', 1, ' ')  
    call add_default('AREI   ', 1, ' ')  
	
	! frequency arrays for above
    call addfld ('FREQL  ','fraction  ',pver, 'A','Fractional occurance of liquid'       ,phys_decomp)
    call addfld ('FREQI  ','fraction  ',pver, 'A','Fractional occurance of ice'       ,phys_decomp)
    call add_default('FREQL   ', 1, ' ')  
    call add_default('FREQI   ', 1, ' ')  

    ! Average rain and snow mixing ratio (Q), number (N) and diameter (D), with frequency
    call addfld ('AQRAIN   ','kg/kg   ',pver, 'A','Average rain mixing ratio'         ,phys_decomp)
    call addfld ('AQSNOW   ','kg/kg   ',pver, 'A','Average snow mixing ratio'         ,phys_decomp)
    call addfld ('ANRAIN   ','m-3     ',pver, 'A','Average rain number conc'         ,phys_decomp)
    call addfld ('ANSNOW   ','m-3     ',pver, 'A','Average snow number conc'         ,phys_decomp)
    call addfld ('ADRAIN   ','Micron  ',pver, 'A','Average rain effective Diameter'         ,phys_decomp)
    call addfld ('ADSNOW   ','Micron  ',pver, 'A','Average snow effective Diameter'         ,phys_decomp)
    call addfld ('FREQR  ','fraction  ',pver, 'A','Fractional occurance of rain'       ,phys_decomp)
    call addfld ('FREQS  ','fraction  ',pver, 'A','Fractional occurance of snow'       ,phys_decomp)
    ! Average cloud top particle size and number (liq, ice) and frequency
    call addfld ('ACTREL  ','Micron  ',1, 'A','Average Cloud Top droplet effective radius'       ,phys_decomp)
    call addfld ('ACTREI  ','Micron  ',1, 'A','Average Cloud Top ice effective radius'       ,phys_decomp)
    call addfld ('ACTNL  ','1/m3  ',1, 'A','Average Cloud Top droplet number'       ,phys_decomp)
    call addfld ('ACTNI  ','1/m3  ',1, 'A','Average Cloud Top ice number'       ,phys_decomp)
    call addfld ('FCTL  ','fraction  ',1, 'A','Fractional occurance of cloud top liquid'       ,phys_decomp)
    call addfld ('FCTI  ','fraction  ',1, 'A','Fractional occurance of cloud top ice'       ,phys_decomp)
    call add_default('ACTREL  ', 1, ' ')  
    call add_default('ACTREI  ', 1, ' ')  
    call add_default('ACTNL   ', 1, ' ')  
    call add_default('ACTNI   ', 1, ' ')  

    !korolev
   call addfld ('uzstar'  ,'m/s     ',pver, 'A','korolev Uzmax',phys_decomp)
   call addfld ('uzzero'  ,'m/s     ',pver, 'A','korolev Uzmin',phys_decomp)
   call addfld ('sigmaw'  ,'m/s     ',pver, 'A','sigmaw',phys_decomp)
   call addfld ('Frac1_BF','fraction',pver, 'A','korolev frac1 U>Uzmax     ',phys_decomp)
   call addfld ('Frac2_BF','fraction',pver, 'A','korolev frac2 BF          ',phys_decomp)
   call addfld ('Frac3_BF','fraction',pver, 'A','korolev frac3 U<Uzmin     ',phys_decomp)
   call addfld ('BFstate' ,'state   ',pver, 'A','235<T<273 and liqice exist',phys_decomp)
   !call add_default('uzstar'  ,1,' ')
   !call add_default('uzzero'  ,1,' ')
   !call add_default('sigmaw'  ,1,' ')
   !call add_default('Frac1_BF',1,' ')
   !call add_default('Frac2_BF',1,' ')
   !call add_default('Frac3_BF',1,' ')
   !call add_default('BFstate' ,1,' ')

   call addfld ('niicBF'  ,'#/m3   ',pver, 'A','in-cloud ice number concentration     ',phys_decomp)
   call addfld ('ncicBF'  ,'#/m3   ',pver, 'A','in-cloud liquid number concentration  ',phys_decomp)
   call addfld ('MixNIN'  ,'#/m3   ',pver, 'A','IN number concentration ',phys_decomp)
   call addfld ('ice_rn'  ,'um/cm3 ',pver, 'A','ice mean radius * number concentration',phys_decomp)
   call addfld ('wtr_rn'  ,'um/cm3 ',pver, 'A','wtr mean radius * number concentration',phys_decomp)
   !call add_default('niicBF'  ,1,' ')
   !call add_default('ncicBF'  ,1,' ')
   !call add_default('ice_rn'  ,1,' ')
   !call add_default('wtr_rn'  ,1,' ')
   !call add_default('MixNIN'  ,1,' ')

   call addfld ('sta1','frequency',pver, 'A','gamil ice and droplet both grow      ',phys_decomp)
   call addfld ('sta2','frequency',pver, 'A','gamil BF ice grow and droplet shrink ',phys_decomp)
   call addfld ('sta3','frequency',pver, 'A','gamil droplet grow and ice shrink    ',phys_decomp)
   call addfld ('sta4','frequency',pver, 'A','gamil ice and droplet both shrink    ',phys_decomp)
   call addfld ('icetend','kg/kg/s',pver, 'A','ice growth tendency',phys_decomp)
   call addfld ('wtrtend','kg/kg/s',pver, 'A','wtr growth tendency',phys_decomp)
   call addfld ('Mixcmei','kg/kg/s',pver, 'A','ice growth tendency',phys_decomp)
   call addfld ('Mixcmel','kg/kg/s',pver, 'A','ice growth tendency',phys_decomp)
   call addfld ('Mixberg','kg/kg/s',pver, 'A','ice growth tendency',phys_decomp)
   call addfld ('qiBF','kg/kg',pver, 'A','ice mass mixing ratio ',phys_decomp)
   call addfld ('qcBF','kg/kg',pver, 'A','wtr mass mixing ratio ',phys_decomp)
   call addfld ('FiceBF','fraction',pver, 'A','ice Fraction in mixed phase cloud ',phys_decomp)
   !call add_default('sta1',1,' ')
   !call add_default('sta2',1,' ')
   !call add_default('sta3',1,' ')
   !call add_default('sta4',1,' ')
   !call add_default('icetend',1,' ')
   !call add_default('wtrtend',1,' ')
   !call add_default('qiBF',1,' ')
   !call add_default('qcBF',1,' ')
   !call add_default('FiceBF',1,' ')
   !call add_default('Mixcmei',1,' ')
   !call add_default('Mixcmel',1,' ')
   !call add_default('Mixberg',1,' ')

   !  ice nucleatin cirrus
  call addfld ('INhet    ','#/m3   ',pver, 'A','Number Concentation  contribtuion from het nucleation in ice cloud',phys_decomp)
  call addfld ('INhom    ','#/m3   ',pver, 'A','Number Concentation  contribtuion from hom nucleation in ice cloud',phys_decomp)
  call addfld ('INice    ','#/m3   ',pver, 'A','Number Concentation  contribtuion from bot hom and het in ice cloud',phys_decomp)
  call addfld ('INFrehom  ','frequency',pver, 'A','hom IN frequency ice cloud',phys_decomp)
  call addfld ('INFrehet  ','frequency',pver, 'A','het IN frequency ice cloud',phys_decomp)
  call addfld ('INFreIN   ','frequency',pver, 'A','IN frequency ice cloud',phys_decomp)
  !call add_default ('INhet  ', 1, ' ')
  !call add_default ('INhom  ', 1, ' ')
  !call add_default ('INice  ', 1, ' ')
  !call add_default ('INFrehom', 1, ' ')
  !call add_default ('INFrehet', 1, ' ')
  !call add_default ('INFreIN ', 1, ' ')
!liuym+++
  call addfld ('pre_out','kg/kg/s',pver, 'A','qc tendency due to',phys_decomp)
  call addfld ('prds_out','kg/kg/s',pver, 'A','qc tendency due to',phys_decomp)
  call addfld ('cmel_out','kg/kg/s',pver, 'A','qc tendency due to',phys_decomp)
  call addfld ('cmei_out','kg/kg/s',pver, 'A','qc tendency due to',phys_decomp)
  call addfld ('bergs_out','kg/kg/s',pver, 'A','qc tendency due to',phys_decomp)
  call addfld ('berg_out','kg/kg/s',pver, 'A','qc tendency due to',phys_decomp)

  call addfld ('npccn_out','#/kg/s',pver, 'A','number conc tendency due to',phys_decomp)
  call addfld ('nsubc_out','#/kg/s',pver, 'A','number conc tendency due to',phys_decomp)
  call addfld ('nnuccd_out','#/kg/s',pver, 'A','number conc tendency due to',phys_decomp)
  call addfld ('nsubi_out','#/kg/s',pver, 'A','number conc tendency due to',phys_decomp)
  call addfld ('nsubs_out','#/kg/s',pver, 'A','number conc tendency due to',phys_decomp)
  call addfld ('nsubr_out','#/kg/s',pver, 'A','number conc tendency due to',phys_decomp)
  call addfld ('prc_out','kg/kg/s',pver, 'A','qc tendency due to autoconversion of cloud droplets',phys_decomp)
  call addfld ('nprc_out','#/kg/s',pver, 'A','number conc tendency due to  autoconversion of cloud droplets',phys_decomp)
  call addfld ('nprc1_out','#/kg/s',pver, 'A','number qr tendency due to autoconversion of cloud droplets',phys_decomp)
   call addfld ('prci_out','kg/kg/s',pver, 'A','mixing rat tendency due to autoconversion of cloud ice to snow',phys_decomp)
   call addfld ('nprci_out','#/kg/s',pver, 'A','number conc tendency due to autoconversion of cloud ice to snow',phys_decomp)
   call addfld ('mnuccc_out','kg/kg/s',pver, 'A','mixing ratio tendency due to freezing of cloud water',phys_decomp)
  call addfld ('nnuccc_out','#/kg/s',pver, 'A','number conc tendency due to freezing of cloud water',phys_decomp)
  call addfld ('nsagg_out','#/kg/s',pver, 'A','ns tendency due to self-aggregation of snow',phys_decomp)
  call addfld ('psacws_out','kg/kg/s',pver, 'A','mixing rat tendency due to collection of droplets by snow',phys_decomp)
  call addfld ('npsacws_out','#/kg/s',pver, 'A','number conc tendency due to collection of droplets by snow',phys_decomp)
  call addfld ('pracs_out','kg/kg/s',pver, 'A','mixing rat tendency due to collection of rain by snow',phys_decomp)
  call addfld ('npracs_out','#/kg/s',pver, 'A','number conc tendency due to collection of rain by snow',phys_decomp)
  call addfld ('mnuccr_out','kg/kg/s',pver, 'A','mixing rat tendency due to freezing of rain',phys_decomp)
  call addfld ('nnuccr_out','#/kg/s',pver, 'A','number conc tendency due to freezing of rain',phys_decomp)
  call addfld ('pra_out','kg/kg/s',pver, 'A','mixing rat tendnency due to accretion of droplets by rain',phys_decomp)
  call addfld ('npra_out','#/kg/s',pver, 'A','nc tendnency due to accretion of droplets by rain',phys_decomp)
  call addfld ('nragg_out','#/kg/s',pver, 'A','nr tendency due to self-collection of rain',phys_decomp)
  call addfld ('prai_out','kg/kg/s',pver, 'A','mixing rat tendency due to accretion of cloud ice by snow',phys_decomp)
  call addfld ('nprai_out','#/kg/s',pver, 'A','number conc tendency due to accretion of cloud ice by snow',phys_decomp)
!liuym---
!++clu
  call addfld ('ZICE'     , 'm6/kg   ', pver, 'A', 'in-cloud ZZZ of ice'                                     ,phys_decomp)
  call addfld ('QICE'     , 'kg/kg   ', pver, 'A', 'in-cloud q of ice'                                       ,phys_decomp)
  call addfld ('NICE'     , '1/kg    ', pver, 'A', 'in-cloud n of ice'                                       ,phys_decomp)
  call addfld ('MIUICE'   , '-'       , pver, 'A', 'miu of ice'                                              ,phys_decomp)
  call addfld ('FREZICE'  , '-'       , pver, 'A', 'occur frequency ice'                                     ,phys_decomp)

  ! call add_default ('ZICE'  , 1, ' ')
  ! call add_default ('NICE'  , 1, ' ')
  ! call add_default ('QICE'  , 1, ' ')
  ! call add_default ('MIUICE', 1, ' ')
  ! call add_default ('FREZICE', 1, ' ')
  
  call addfld ('Mufallqi' , 'm/s '    , pver, 'A', 'mass-weighted ice fall velocity'                         ,phys_decomp)
  call addfld ('Mufallni' , 'm/s '    , pver, 'A', 'number-weighted ice fall velocity'                       ,phys_decomp)
  call addfld ('Mufallpro', '-'       , pver, 'A', 'fall velocity probability'                               ,phys_decomp)
  call addfld ('MuR3_fall', '-'       , pver, 'A', 'ice 3 moment radius (fallspeed)'                         ,phys_decomp)

  call addfld ('MuR3_before_fall', '-', pver, 'A', 'ice 3 moment radius (fallspeed) before fall'             ,phys_decomp)
  call addfld ('Mubeforefallpro', '-' , pver, 'A', 'before fall velocity probability'                        ,phys_decomp)
  call addfld ('MuR3_after_fall', '-' , pver, 'A', 'ice 3 moment radius (fallspeed) after fall'              ,phys_decomp)
  call addfld ('Muafterfallpro', '-'  , pver, 'A', 'after ice 3 moment radius (fallspeed)'                   ,phys_decomp)
  call addfld ('Mufqprci' , '-'       , pver, 'A', 'autoconversion of cloud ice to snow;  mass ratio'        ,phys_decomp)
  call addfld ('Mufnprci' , '-'       , pver, 'A', 'autoconversion of cloud ice to snow;  number ratio'      ,phys_decomp)
  call addfld ('hxepsi'  , '1/s'      , pver, 'A', 'bergeron process,in-cloud 1/sat relaxation timescale for ice' ,phys_decomp)
  call addfld ('hxdvi'  , 'mm'        , pver, 'A', 'volume-mean diameter' 									 ,phys_decomp)
  ! call add_default ('Mufallqi', 1, ' ')
  ! call add_default ('Mufallni', 1, ' ')
  ! call add_default ('Mufallpro', 1, ' ')
  ! call add_default ('MuR3_fall', 1, ' ')

  ! call add_default ('MuR3_before_fall', 1, ' ')
  ! call add_default ('Mubeforefallpro', 1, ' ')
  ! call add_default ('MuR3_after_fall', 1, ' ')
  ! call add_default ('Muafterfallpro', 1, ' ')
  
  ! call add_default ('Mufqprci', 1, ' ')
  ! call add_default ('Mufnprci', 1, ' ')
  ! call add_default ('hxepsi'  , 1, ' ')
!  call addfld ('ICECLDF'  , 'fraction'  ,pver, 'A','Ice cloud fraction',phys_decomp )
!  call addfld ('LIQCLDF'  , 'fraction'  ,pver, 'A','Liq cloud fraction',phys_decomp )
!  call addfld ('RHCLOUD' , 'fraction'  ,pver, 'A','CLU cloud fraction',phys_decomp )  
!--clu
#endif
!--sxj

    return
  end subroutine cldcond_init

!===============================================================================
!!$  subroutine cldcond_tend(state, ptend, dt, pbuf)
  subroutine cldcond_tend(state, ptend, dt, &
       tcwato, qcwato, lcwato, precip, snow, icefrac, rhdfda, rhu00, cldn, evapprec, prodprec, cme, snowh)
!
! Compute the tendency of cloud water mixing ratios (liquid and ice) from
! microphysical processes and sedimenation of cloud particles.
! 
! Author: Byron Boville  Sept 04, 2002
!  modified pjr: march 2003 to repartition snow/rain
!
!-----------------------------------------------------------------------

    use physics_types, only: physics_state, physics_ptend
!!$    use phys_buffer,   only: pbuf_size_max, pbuf_fld
    use history,       only: outfld
    use cldwat,        only: pcond, cldwat_fice
    use physconst,     only: tmelt
!wtw
    use wu_cdnc,       only: get_cdnc
!wtw
    implicit none

! Arguments
    type(physics_state), intent(inout) :: state          ! state variables
    type(physics_ptend), intent(inout) :: ptend          ! package tendencies
    real(r8),            intent(in)    :: dt             ! timestep
!!$    type(pbuf_fld), intent(inout), dimension(pbuf_size_max) :: pbuf  ! physics buffer
    real(r8), intent(in)  :: tcwato(pcols,pver)        !cloud water old temperature
    real(r8), intent(in)  :: qcwato(pcols,pver)        ! cloud water old q
    real(r8), intent(in)  :: lcwato(pcols,pver)        ! cloud liquid water old q
    real(r8), intent(in)  :: icefrac(pcols)            ! sea ice fraction (fraction)
    real(r8), intent(in)  :: rhdfda(pcols,pver)        ! dRh/dcloud, old 
    real(r8), intent(in)  :: rhu00 (pcols,pver)        ! Rh threshold for cloud, old
    real(r8), intent(in)  :: cldn(pcols,pver)          !new cloud fraction
    real(r8), intent(in) :: snowh(pcols)         ! Snow depth over land, water equivalent (m)

    real(r8), intent(out) :: precip(pcols)             ! sfc flux of precip (m/s)
    real(r8), intent(out) :: snow (pcols)              ! sfc flux of snow   (m/s)
    real(r8), intent(out) :: evapprec(pcols,pver)          ! local evaporation of precipitation
    real(r8), intent(out) :: prodprec(pcols,pver)          ! local production of precipitation
    real(r8), intent(out) :: cme     (pcols,pver)          ! local condensation - evaporation of cloud water

! Local variables
    integer :: lchnk                          ! chunk identifier
    integer :: ncol                           ! number of atmospheric columns in chunk
    integer :: i,k                            ! loop indexes
!!$    real(r8), pointer, dimension(:)   :: buffld1  ! physics buffer field1
!!$    real(r8), pointer, dimension(:,:) :: buffld2  ! physics buffer field2
    real(r8) :: rdt                          ! 1./dt
    real(r8) :: qtend (pcols,pver)            ! moisture tendencies
    real(r8) :: ttend (pcols,pver)            ! temperature tendencies
    real(r8) :: ltend (pcols,pver)            ! cloud liquid water tendencies
!    real(r8) :: cme     (pcols,pver)          ! local condensation - evaporation of cloud water
    real(r8) :: evapheat(pcols,pver)          ! heating rate due to evaporation of precip
!    real(r8) :: evapprec(pcols,pver)          ! local evaporation of precipitation
    real(r8) :: evapsnow(pcols,pver)          ! local evaporation of snow
    real(r8) :: prfzheat(pcols,pver)          ! heating rate due to freezing of precip (W/kg)
    real(r8) :: meltheat(pcols,pver)          ! heating rate due to phase change of precip
    real(r8) :: prodsnow(pcols,pver)          ! local production of snow
    real(r8) :: totcw   (pcols,pver)          ! total cloud water mixing ratio
    real(r8) :: fice    (pcols,pver)          ! Fractional ice content within cloud
    real(r8) :: fsnow   (pcols,pver)          ! Fractional snow production
    real(r8) :: repartht(pcols,pver)          ! heating rate due to phase repartition of input precip
    real(r8) :: icimr(pcols,pver)             ! in cloud ice mixing ratio
    real(r8) :: icwmr(pcols,pver)             ! in cloud water mixing ratio
    real(r8) fwaut(pcols,pver)              
    real(r8) fsaut(pcols,pver)              
    real(r8) fracw(pcols,pver)              
    real(r8) fsacw(pcols,pver)              
    real(r8) fsaci(pcols,pver)              
    real(r8) ice2pr(pcols,pver)   ! rate of conversion of ice to precip
    real(r8) liq2pr(pcols,pver)   ! rate of conversion of liquid to precip
    real(r8) liq2snow(pcols,pver)   ! rate of conversion of liquid to snow
    real(r8) hs1, qv1, ql1, qi1, qs1, qr1, fice2, pr1, w1, w2, w3, fliq, res(pcols,pver)
    real(r8) temp(pcols), w4, wl, wv, wi, wlf, wvf, wif, qif, qlf, qvf, res2
!wtw
    real(r8) :: rel(pcols,pver)      ! Liquid effective drop size (microns)
    real(r8) :: cdnc(pcols,pver)     ! Liquid cloud droplet number concentration (1/cm3)
    real(r8) :: rei(pcols,pver)      ! Ice effective drop size (microns)

!-----------------------------------------------------------------------

! Set output flags
    ptend%name         = 'cldwat'
    ptend%ls           = .true.
    ptend%lq(1)        = .true.
    ptend%lq(ixcldice) = .true.
    ptend%lq(ixcldliq) = .true.

! Initialize chunk id and size
    lchnk = state%lchnk
    ncol  = state%ncol

! associate local pointers with fields in the physics buffer
!!$    buffld1 => pbuf(ixbuffld1)%fld_ptr(1,1:pcols,1,     lchnk,1)
!!$    buffld2 => pbuf(ixbuffld2)%fld_ptr(1,1:pcols,1:pver,lchnk,1)

! Define fractional amount of cloud condensate in ice phase
    call cldwat_fice(ncol, state%t, fice, fsnow)

! compute total cloud water
    totcw(:ncol,:pver) = state%q(:ncol,:pver,ixcldice) + state%q(:ncol,:pver,ixcldliq)

! save input cloud ice
    repartht(:ncol,:pver) = state%q(:ncol,:pver,ixcldice)

! Repartition ice and liquid
    state%q(:ncol,:pver,ixcldice) = totcw(:ncol,:pver) * fice(:ncol,:pver)
    state%q(:ncol,:pver,ixcldliq) = totcw(:ncol,:pver) * (1.0_r8 - fice(:ncol,:pver))

! Determine heating from change in cloud ice
    repartht(:ncol,:pver) = latice/dt * (state%q(:ncol,:pver,ixcldice) - repartht(:ncol,:pver))

! calculate the tendencies for moisture, temperature and cloud fraction
    rdt = 1./dt
    qtend(:ncol,:pver) = (state%q(:ncol,:pver,1) - qcwato(:ncol,:pver))*rdt
    ttend(:ncol,:pver) = (state%t(:ncol,:pver)   - tcwato(:ncol,:pver))*rdt
    ltend(:ncol,:pver) = (totcw  (:ncol,:pver)   - lcwato(:ncol,:pver))*rdt

! call microphysics package to calculate tendencies
    call t_startf('pcond')

!wtw
    call get_cdnc (lchnk   ,ncol    , state,   cldn,  rel, rei, cdnc )
!
    call pcond (lchnk,   ncol, &
         cdnc,                 &   ! wutw
         state%t  , ttend      , state%q(1,1,1), qtend   , state%omega, &
         totcw    , state%pmid , state%pdel , cldn       , fice       , fsnow, &
         cme      , prodprec   , prodsnow, evapprec   , evapsnow   , evapheat   , prfzheat, &
         meltheat , precip     , snow       , dt         , fwaut      , &
         fsaut    , fracw      , fsacw      , fsaci      , ltend      , &
         rhdfda   , rhu00      , icefrac    , state%zi   , ice2pr, liq2pr, liq2snow, snowh)
    call t_stopf('pcond')

! make it interactive
    do k = 1,pver
       do i = 1,ncol
          ptend%s(i,k)          = cme(i,k) * (latvap + latice*fice(i,k)) &
               + evapheat(i,k) + prfzheat(i,k) + meltheat(i,k) + repartht(i,k)

          ptend%q(i,k,1)        =-cme(i,k) + evapprec(i,k)

          ptend%q(i,k,ixcldice) =cme(i,k)*fice(i,k)      - ice2pr(i,k)
          ptend%q(i,k,ixcldliq) =cme(i,k)*(1.-fice(i,k)) - liq2pr(i,k)
       end do
    end do

! sanity checks to be removed after 5 year run
    res(:ncol,:) = state%q(:ncol,:,1) + ptend%q(:ncol,:,1)*dt
    if (any(res(:ncol,:) < 0)) then
       do k = 1,pver
          do i = 1,ncol
             if (res(i,k).lt.0.) then
                write (6,*) ' predicted neg vapor ', i,k,lchnk, res(i,k)
                write (6,*) ' oldq, cme(i,k)*dt, evapprec(i,k)*dt ', state%q(i,k,1), &
                   -cme(i,k)*dt, evapprec(i,k)*dt, ptend%q(i,k,1)*dt
                call endrun('CLDCOND_TEND')
             endif
          end do
       end do
    endif
    res(:ncol,:) = state%q(:ncol,:,ixcldice) + ptend%q(:ncol,:,ixcldice)*dt
    if (any(res(:ncol,:) < 0)) then
       do k = 1,pver
          do i = 1,ncol
             if (res(i,k).lt.0.) then
                write (6,*) ' cldcond_tend: predicted neg ice ', i,k,lchnk, res(i,k), fice(i,k)
                write (6,*) ' oldice, cme(i,k)*dt*(fice(i,k)), ice2pr*dt, ptend_ice*dt', &
                   state%q(i,k,ixcldice), &
                   cme(i,k)*dt*(fice(i,k)), ice2pr(i,k)*dt, ptend%q(i,k,ixcldice)*dt
                write (6,*) ' oldcond, cme(i,k)*dt, prodprec*dt, ptend_cond*dt', &
                   state%q(i,k,ixcldice)+state%q(i,k,ixcldliq), &
                   cme(i,k)*dt, prodprec(i,k)*dt, &
                   (ptend%q(i,k,ixcldice)+ptend%q(i,k,ixcldliq))*dt
                call endrun('CLDCOND_TEND')
             endif
          end do
       end do
    endif
! end sanity check

#ifdef DEBUG
    do i = 1,ncol
       pr1 = 0
       hs1 = 0
       qv1 = 0
       ql1 = 0
       qi1 = 0
       qs1 = 0
       qr1 = 0
       w1 = 0
       wl = 0
       wv = 0
       wi = 0
       wlf = 0
       wvf = 0
       wif = 0
       do k = 1,pver

          if (lchnk.eq.248.and.i.eq.12) then

             write (6,*) 
             write (6,*) ' input state, t, q, l, i ', k, state%t(i,k), state%q(i,k,1), state%q(i,k,ixcldliq),  state%q(i,k,ixcldice)
             write (6,*) ' rain, snow, total from components before accumulation ', qr1, qs1, qr1+qs1
             write (6,*) ' total precip before accumulation                      ', k, pr1

             wv = wv + state%q(i,k,1       )*state%pdel(i,k)/gravit
             wl = wl + state%q(i,k,ixcldliq)*state%pdel(i,k)/gravit
             wi = wi + state%q(i,k,ixcldice)*state%pdel(i,k)/gravit

             qvf = state%q(i,k,1) + ptend%q(i,k,1)*dt
             qlf = state%q(i,k,ixcldliq) + ptend%q(i,k,ixcldliq)*dt
             qif = state%q(i,k,ixcldice) + ptend%q(i,k,ixcldice)*dt

             if (qvf.lt.0.) then
                write (6,*) ' qvf is negative *******', qvf
             endif
             if (qlf.lt.0.) then
                write (6,*) ' qlf is negative *******', qlf
             endif
             if (qif.lt.0.) then
                write (6,*) ' qif is negative *******', qif
             endif
             write (6,*) ' qvf, qlf, qif ', qvf, qlf, qif

             wvf = wvf + qvf*state%pdel(i,k)/gravit
             wlf = wlf + qlf*state%pdel(i,k)/gravit
             wif = wif + qif*state%pdel(i,k)/gravit

             hs1 = hs1 + ptend%s(i,k)*state%pdel(i,k)/gravit
             pr1 = pr1 + state%pdel(i,k)/gravit*(prodprec(i,k)-evapprec(i,k))
             qv1 = qv1 - (cme(i,k)-evapprec(i,k))*state%pdel(i,k)/gravit    ! vdot
             w1  = w1  + (cme(i,k)-prodprec(i,k))*state%pdel(i,k)/gravit    ! cdot
             qi1 = qi1 + ((cme(i,k))*fice(i,k)        -ice2pr(i,k) )*state%pdel(i,k)/gravit   ! idot
             ql1 = ql1 + ((cme(i,k))*(1._r8-fice(i,k))-liq2pr(i,k) )*state%pdel(i,k)/gravit   ! ldot

             qr1 = qr1 &
                  + ( liq2pr(i,k)-liq2snow(i,k)   &     ! production of rain
                   -(evapprec(i,k)-evapsnow(i,k)) &     ! rain evaporation
                    )*state%pdel(i,k)/gravit
             qs1 = qs1 &
                  + ( ice2pr(i,k) + liq2snow(i,k) &     ! production of snow.Note last term has phase change
                     -evapsnow(i,k)               &     ! snow evaporation
                    )*state%pdel(i,k)/gravit

             if (state%t(i,k).gt.tmelt) then
                qr1 = qr1 + qs1
                qs1 = 0.
             endif
             write (6,*) ' rain, snow, total after accumulation ', qr1, qs1, qr1+qs1
             write (6,*) ' total precip after accumulation      ', k, pr1
             write (6,*)
             write (6,*) ' layer prodprec, evapprec, pdel ', prodprec(i,k), evapprec(i,k), state%pdel(i,k)
             write (6,*) ' layer prodsnow, ice2pr+liq2snow ', prodsnow(i,k), ice2pr(i,k)+liq2snow(i,k)
             write (6,*) ' layer prodprec-prodsnow, liq2pr-liq2snow ', prodprec(i,k)-prodsnow(i,k), liq2pr(i,k)-liq2snow(i,k)
             write (6,*) ' layer evapsnow, evaprain ', k, evapsnow(i,k), evapprec(i,k)-evapsnow(i,k)
             write (6,*) ' layer ice2pr, liq2pr, liq2snow ', ice2pr(i,k), liq2pr(i,k), liq2snow(i,k)
             write (6,*) ' layer ice2pr+liq2pr, prodprec ', ice2pr(i,k)+liq2pr(i,k), prodprec(i,k)
             write (6,*)
             write (6,*) ' qv1 vapor removed from col after accum  (vdot)   ', k, qv1
             write (6,*) ' - (precip produced - vapor removed) after accum  ', k, -pr1-qv1
             write (6,*) ' condensate produce after accum                   ', k, w1
             write (6,*) ' liq+ice tends accum                              ', k, ql1+qi1
             write (6,*) ' change in total water after accum                ', k, qv1+ql1+qi1
             write (6,*) ' imbalance in colum after accum                   ', k, qs1+qr1+qv1+ql1+qi1
             write (6,*) ' fice at this lev ', fice(i,k)
             write (6,*)

             res2 = abs((qs1+qr1+qv1+ql1+qi1)/max(abs(qv1),abs(ql1),abs(qi1),abs(qs1),abs(qr1),1.e-36))
             write (6,*) ' relative residual in column method 1             ', k, res2

             write (6,*) ' relative residual in column method 2             ', k, abs((qs1+qr1+qv1+ql1+qi1)/max(abs(qv1+ql1+qi1),1.e-36))
!            if (abs((qs1+qr1+qv1+ql1+qi1)/(qs1+qr1+1.e-36)).gt.1.e-14) then
             if (res2.gt.1.e-14) then
                call endrun ('CLDCOND_TEND')
             endif

!             w3  = cme(i,k) * (latvap + latice*fice(i,k)) &
!               + evapheat(i,k) + prfzheat(i,k) + meltheat(i,k)

             res2 = qs1+qr1-pr1
             w4 = max(abs(qs1),abs(qr1),abs(pr1)) 
             if (w4.gt.0.)  then
                if (res2/w4.gt.1.e-14) then
                   write (6,*) ' imbalance in precips calculated two ways '
                   write (6,*) ' res2/w4, pr1, qr1, qs1, qr1+qs1 ', &
                        res2/w4, pr1, qr1, qs1, qr1+qs1
!                   call endrun()
                endif
             endif
             if (k.eq.pver) then
                write (6,*) ' pcond returned precip, rain and snow rates ', precip(i), precip(i)-snow(i), snow(i)
                write (6,*) ' I calculate ', pr1, qr1, qs1
!               call endrun
                write (6,*) ' byrons water check ', wv+wl+wi-pr1*dt, wvf+wlf+wif
             endif
          write (6,*)
          endif
       end do
    end do
#endif

#ifdef DEBUG
    if (.true.) then
    do i = 1,ncol
       if (snow(i).gt.0.01/8.64e4.and.state%t(i,pver).gt.273.16) then
          write (6,*) ' cldcond: snow, temp, ', i, lchnk, &
               snow(i), state%t(i,pver)
          write (6,*) ' t ', state%t(i,:)
          write (6,*) ' fsaut ', fsaut(i,:)
          write (6,*) ' fsacw ', fsacw(i,:)
          write (6,*) ' fsaci ', fsaci(i,:)
          write (6,*) ' meltheat ', meltheat(i,:)
          call endrun ('CLDCOND_TEND')
       endif

       if (snow(i)*8.64e4.lt.-1.e-5) then
          write (6,*) ' neg snow ', snow(i)*8.64e4
          write (6,*) ' cldcond: snow, temp, ', i, lchnk, &
               snow(i), state%t(i,pver)
          write (6,*) ' t ', state%t(i,:)
          write (6,*) ' fsaut ', fsaut(i,:)
          write (6,*) ' fsacw ', fsacw(i,:)
          write (6,*) ' fsaci ', fsaci(i,:)
          write (6,*) ' meltheat ', meltheat(i,:)
          call endrun ('CLDCOND_TEND')
       endif
    end do
    endif
#endif

! Compute in cloud ice and liquid mixing ratios
    do k=1,pver
       do i = 1,ncol
          icimr(i,k) = (state%q(i,k,ixcldice) + dt*ptend%q(i,k,ixcldice)) / max(0.01_r8,cldn(i,k))
          icwmr(i,k) = (state%q(i,k,ixcldliq) + dt*ptend%q(i,k,ixcldliq)) / max(0.01_r8,cldn(i,k))
       end do
    end do

! convert precipitation from kg/m2 to m/s
    snow  (:ncol) = snow  (:ncol)/1000.
    precip(:ncol) = precip(:ncol)/1000.

! record history variables
    call outfld('FWAUT',fwaut,    pcols,lchnk)
    call outfld('FSAUT',fsaut,    pcols,lchnk)
    call outfld('FRACW',fracw,    pcols,lchnk)
    call outfld('FSACW',fsacw,    pcols,lchnk)
    call outfld('FSACI',fsaci,    pcols,lchnk)
    call outfld('ICIMR',icimr,    pcols,lchnk)
    call outfld('ICWMR',icwmr,    pcols,lchnk)

    call outfld('PCSNOW'  ,snow    , pcols,lchnk)
    call outfld('FICE'    ,fice    , pcols,lchnk)
    call outfld('CME'     ,cme     , pcols,lchnk)
    call outfld('PRODPREC',prodprec, pcols,lchnk)
    call outfld('EVAPPREC',evapprec, pcols,lchnk)
    call outfld('EVAPSNOW',evapsnow, pcols,lchnk)
    call outfld('HPROGCLD',ptend%s , pcols,lchnk)
    call outfld('HEVAP   ',evapheat, pcols,lchnk)
    call outfld('HMELT'   ,meltheat, pcols,lchnk)
    call outfld('HREPART' ,repartht, pcols,lchnk)

! update boundary quantities
!!$    ptend%hflx_srf = 0.
!!$    ptend%hflx_top = 0.
!!$    ptend%cflx_srf = 0.
!!$    ptend%cflx_top = 0.
  end subroutine cldcond_tend

!++sxj
#ifdef MG08

! Interface to sedimentation, detrain, cloud fraction and  microphysics subroutines
!
   subroutine stratiform_tend(state, ptend_all, dtime, &
       icefrac, landfrac, ocnfrac, &
       snowh,    dlf, &
       ls_flxprc,   ls_flxsnw,  ls_reffrain, ls_reffsnow , &       ! wtw     
       cmfmc,  cmfmc2,   &
       effliq, effice,    &                          ! wtw
       tdeepcld, qdeepcld,deepcld ,tshallowcld, qshallowcld, shallowcld  ,& 
       cld_env  ,t_env   ,q_env     , &
       cnt, cnb, lcl, cstcld, &  ! sxj  add deepcld ,shallowcld,cnt, cnb
       cld, concld, cldst, prain, nevapr, qme,  &  ! sxj
       ts,      sst,      zdu,  &  !sxj
       prec_sed, snow_sed, prec_pcw, snow_pcw, &   ! sxj  prec_str snow_str
       tcwato, qcwato, lcwato )  !sxj
       !pbuf )  ! sxj get pbuf from phys_buffer module
    use shr_kind_mod,     only: r8 => shr_kind_r8
    use ppgrid
    use cloud_fraction,   only: cldfrc  ! sxj bcc
    use physics_types,    only: physics_state, physics_ptend, physics_tend
    use physics_types,    only: physics_ptend_init, physics_update
    use physics_types,    only: physics_ptend_sum,  physics_state_copy, physics_tend_init
    use history,          only: outfld
    use phys_buffer,      only: pbuf,pbuf_size_max, pbuf_fld, pbuf_old_tim_idx, pbuf_get_fld_idx
    use constituents,     only: cnst_get_ind
    use cldwat,           only: pcond, cldwat_fice
    use cldwat2m,         only: mmicro_pcond
    use prescribed_aerosols,    only: naer_all,naer !naer_all=14,naer=12
#ifdef MAC_SP
    use prescribed_aerosols,      only: get_aerosol,  get_int_scales  ! sxj 
    use prescribed_aerosols_cdnc, only: get_aerosol_cdnc
#else
    use prescribed_aerosols,      only: get_aerosol,  get_int_scales  ! sxj
#endif

    use time_manager,     only: is_first_step
    implicit none
!   Parameters
    real(r8) pnot                  ! reference pressure
    parameter (pnot = 1.e5_r8)
!   Input arguments
    type(physics_state), intent(in)    :: state     ! state variables
    type(physics_ptend), intent(out) :: ptend_all   ! package tendencies
!    type(pbuf_fld),      intent(inout), dimension(pbuf_size_max) :: pbuf       ! sxj get pbuf from phys_buffer module
    real(r8), intent(in)  :: dtime                ! timestep
    real(r8), intent(in)  :: icefrac (pcols)      ! sea ice fraction (fraction)
    real(r8), intent(in)  :: landfrac(pcols)      ! land fraction (fraction)
    real(r8), intent(in)  :: ocnfrac (pcols)      ! ocean fraction (fraction)
    real(r8), intent(in) :: snowh(pcols)          ! Snow depth over land, water equivalent (m)
    real(r8), intent(in) :: dlf(pcols,pver)       ! detrained water from ZM
    !real(r8), intent(in) :: rliq(pcols)           ! sxj  vertical integral of liquid not yet in q(ixcldliq) 
    real(r8), intent(in) :: cmfmc(pcols,pverp)    ! convective mass flux--m sub c
    real(r8), intent(in) :: cmfmc2(pcols,pverp)    ! shallow convective mass flux--m sub c
    real(r8), intent(inout) :: deepcld(pcols,pver)    ! sxj 
    real(r8), intent(inout) :: shallowcld(pcols,pver)  ! sxj 

    real(r8), intent(inout) :: qdeepcld(pcols,pver)      ! water vapor inside deep convective cld
    real(r8), intent(inout) :: qshallowcld(pcols,pver)   ! water vapor inside shallow convective cld
    real(r8), intent(inout) :: tdeepcld(pcols,pver)
    real(r8), intent(inout) :: tshallowcld(pcols,pver)

    real(r8), intent(out) :: q_env(pcols,pver)
    real(r8), intent(out) :: t_env(pcols,pver)
    real(r8), intent(out) :: cld_env(pcols,pver)

    real(r8), intent(in) :: cnt(pcols)            ! sxj top level of convection
    real(r8), intent(in) :: cnb(pcols)            ! sxj  bottom level of convection
    integer , intent(inout) :: lcl(pcols)            ! sxj 
    real(r8), intent(in) :: ts(pcols)             ! surface temperature
    real(r8), intent(in) :: sst(pcols)            ! sea surface temperature  !! used for calculate cloud fraction
    real(r8), intent(in)  :: tcwato(pcols,pver)        !!! sxj cloud water old temperature  
    real(r8), intent(in)  :: qcwato(pcols,pver)        !!! sxj cloud water old q
    real(r8), intent(in)  :: lcwato(pcols,pver)        !!! sxj cloud liquid water old q  
    real(r8), intent(in) :: zdu(pcols,pver)       ! detrainment rate from deep convection
    !real(r8), intent(out)  :: prec_str(pcols)  ! sxj [Total] sfc flux of precip from stratiform (m/s)
    !real(r8), intent(out)  :: snow_str(pcols)  ! sxj [Total] sfc flux of snow from stratiform   (m/s)
    real(r8), intent(out)  :: prec_sed(pcols)  ! surface flux of total cloud water from sedimentation
    real(r8), intent(out)  :: snow_sed(pcols) ! surface flux of cloud ice from sedimentation
    real(r8), intent(out)  :: prec_pcw(pcols) ! sfc flux of precip from microphysics(m/s)
    real(r8), intent(out)  :: snow_pcw(pcols) ! sfc flux of snow from microphysics (m/s)
    real(r8), intent(inout)  :: cstcld(pcols,pver)   ! sxj  ,    wzz for layered cloud cover
    real(r8), intent(inout)  ::    cld(pcols,pver)   ! sxj  
    real(r8), intent(inout)  :: concld(pcols,pver)   ! sxj
    real(r8), intent(inout)  ::  cldst(pcols,pver)   ! sxj 
    real(r8), intent(out)  ::   qme(pcols,pver)    ! sxj
    real(r8), intent(out)  :: prain(pcols,pver)    ! sxj
    real(r8), intent(out)  :: nevapr(pcols,pver)   ! sxj
!--------------
! wtw

    real(r8), intent(out) :: ls_flxprc(pcols,pver+1)  ! grid-box average (rain +snow) flux (kg m^-2 s^-1)
    real(r8), intent(out) :: ls_flxsnw(pcols,pver+1)  ! grid-box average snow flux (kg m^-2 s^-1)
    real(r8), intent(out) :: ls_reffrain(pcols,pver)    ! rain effective drop radius (microns)
    real(r8), intent(out) :: ls_reffsnow(pcols,pver)    ! snow effective drop size (microns)
!--------------
!liuym+++
        real(r8) ::       prc_out(pcols,pver)
        real(r8) ::       nprc_out(pcols,pver)
        real(r8) ::       nprc1_out(pcols,pver)
        real(r8) ::       prci_out(pcols,pver)
        real(r8) ::       nprci_out(pcols,pver)
        real(r8) ::       mnuccc_out(pcols,pver)
        real(r8) ::       nnuccc_out(pcols,pver)
        real(r8) ::       nsagg_out(pcols,pver)
        real(r8) ::       psacws_out(pcols,pver)
        real(r8) ::       npsacws_out(pcols,pver)
        real(r8) ::       pracs_out(pcols,pver)
        real(r8) ::       npracs_out(pcols,pver)
        real(r8) ::       mnuccr_out(pcols,pver)
        real(r8) ::       nnuccr_out(pcols,pver)
        real(r8) ::       nnuccd_out(pcols,pver)
        real(r8) ::       pra_out(pcols,pver)
        real(r8) ::       npra_out(pcols,pver)
        real(r8) ::       nragg_out(pcols,pver)
        real(r8) ::       prai_out(pcols,pver)
        real(r8) ::       nprai_out(pcols,pver)
!liuym---

    ! Local variables
    type(physics_state)  :: state1       ! local copy of the state variable
    type(physics_tend ) :: tend          ! Physics tendencies (empty, needed for physics_update call)
    type(physics_ptend)  :: ptend_loc    ! package tendencies
    integer i,k,m
    integer :: lchnk                  ! chunk identifier
    integer :: ncol                   ! number of atmospheric columns
    ! physics buffer fields
    integer itim, ifld
    !real(r8), pointer, dimension(:,:) :: qcwat   ! cloud water old q
    !real(r8), pointer, dimension(:,:) :: tcwat   ! cloud water old temperature
    !real(r8), pointer, dimension(:,:) :: lcwat   ! cloud liquid water old q
    !real(r8), pointer, dimension(:,:) ::  cld    ! cloud fraction
    !real(r8), pointer, dimension(:,:) ::  concld ! convective cloud fraction
    !real(r8), pointer, dimension(:,:) :: qme
    !real(r8), pointer, dimension(:,:) :: prain
    !real(r8), pointer, dimension(:,:) :: nevapr
    real(r8), pointer, dimension(:,:) :: rel     ! liquid effective drop radius (microns)
    real(r8), pointer, dimension(:,:) :: rei     ! ice effective drop size (microns)
    real(r8) :: shfrc(pcols,pver)                ! cloud fraction from shallow convection scheme
    real(r8), pointer, dimension(:,:) :: rel_fn  ! ice effective drop size at fixed number (indirect effect) (microns)
    real(r8), pointer, dimension(:,:) :: kkvh    ! heat flux for cldwat  !
   

    !! local variables for stratiform_sediment
    real(r8) :: rain(pcols)                       ! surface flux of cloud liquid
    real(r8) :: pvliq(pcols,pver+1)               ! vertical velocity of cloud liquid drops (Pa/s)
    real(r8) :: pvice(pcols,pver+1)               ! vertical velocity of cloud ice particles (Pa/s)
    real(r8) :: prec_str(pcols)  ! sxj [Total] sfc flux of precip from stratiform (m/s)
    real(r8) :: snow_str(pcols)  ! sxj [Total] sfc flux of snow from stratiform   (m/s)
    !! local variables for cldfrc
    !real(r8)  cldst(pcols,pver)     ! cloud fraction
    real(r8)  clc(pcols)            ! column convective cloud amount
    real(r8) rhdfda(pcols,pver)     ! d_RH/d_cloud_fraction    ====wlin
    real(r8) rhu00(pcols,pver)      ! RH limit, U00             ====wlin
    real(r8) relhum(pcols,pver)       ! RH, output to determine drh/da
    real(r8) rhu002(pcols,pver)       ! same as rhu00 but for perturbed rh
    real(r8) cld2(pcols,pver)         ! same as cld but for perturbed rh
    real(r8) concld2(pcols,pver)      ! same as concld but for perturbed rh
    real(r8) cldst2(pcols,pver)       ! same as cldst but for perturbed rh
    real(r8) relhum2(pcols,pver)      ! RH after  perturbation
    real(r8) :: pmid(pcols,pver)      ! midpoint pressures
    real(r8) :: t(pcols,pver)         ! temperature
    real(r8) :: q(pcols,pver)         ! specific humidity
    real(r8) :: omga(pcols,pver)      ! vertical pressure velocity
    real(r8) :: phis(pcols)           ! surface geopotential
    real(r8) :: pdel(pcols,pver)      ! pressure depth of layer
    real(r8) :: ps(pcols)             ! surface pressure
    real(r8) :: pint(pcols,pver+1)     ! interface pressures
    real(r8) :: zm(pcols,pver)        !  sxj 
!	real(r8) :: icecldf(pcols,pver) 	!clu
!	real(r8) :: liqcldf(pcols,pver)		!clu
!	real(r8) :: rhcloud(pcols,pver)	!clu
    !! local variables for microphysics
    real(r8) :: rdtime                        ! 1./dtime
    real(r8) :: qtend (pcols,pver)            ! moisture tendencies
    real(r8) :: ttend (pcols,pver)            ! temperature tendencies
    real(r8) :: ltend (pcols,pver)            ! cloud liquid water tendencies
    real(r8) :: evapheat(pcols,pver)          ! heating rate due to evaporation of precip
    real(r8) :: evapsnow(pcols,pver)          ! local evaporation of snow
    real(r8) :: prfzheat(pcols,pver)          ! heating rate due to freezing of precip (W/kg)
    real(r8) :: meltheat(pcols,pver)          ! heating rate due to phase change of precip
    real(r8) :: cmeheat (pcols,pver)          ! heating rate due to phase change of precip
    real(r8) :: prodsnow(pcols,pver)          ! local production of snow
    real(r8) :: totcw   (pcols,pver)          ! total cloud water mixing ratio
    real(r8) :: fice    (pcols,pver)          ! Fractional ice content within cloud
    real(r8) :: fsnow   (pcols,pver)          ! Fractional snow production
    real(r8) :: repartht(pcols,pver)          ! heating rate due to phase repartition of input precip
    real(r8) :: icimr(pcols,pver)             ! in cloud ice mixing ratio
    real(r8) :: icwmr(pcols,pver)             ! in cloud water mixing ratio
    real(r8) fwaut(pcols,pver)
    real(r8) fsaut(pcols,pver)
    real(r8) fracw(pcols,pver)
    real(r8) fsacw(pcols,pver)
    real(r8) fsaci(pcols,pver)
    real(r8) cmeice(pcols,pver)   ! Rate of cond-evap of ice within the cloud
    real(r8) cmeliq(pcols,pver)   ! Rate of cond-evap of liq within the cloud
    real(r8) ice2pr(pcols,pver)   ! rate of conversion of ice to precip
    real(r8) liq2pr(pcols,pver)   ! rate of conversion of liquid to precip
    real(r8) liq2snow(pcols,pver) ! rate of conversion of liquid to snow
    real(r8) temp(pcols)
    real(r8) res(pcols,pver)
    !! variables for morrison microphysics
    real(r8) :: dum1,dum2
    real(r8) :: qc(pcols,pver)
    real(r8) :: qi(pcols,pver)
    real(r8) :: nc(pcols,pver)
    real(r8) :: ni(pcols,pver)
	real(r8) :: zi(pcols,pver)		!+clu
    real(r8) :: icinc(pcols,pver)             ! in cloud ice number conc
    real(r8) :: cdnumc(pcols)                 ! vertically-integrated droplet concentration
    real(r8) :: icwnc(pcols,pver)             ! in cloud water number conc

!-----------------------
! wtw
!
    real(r8), intent(out) :: effliq(pcols,pver)            ! in cloud liq eff rad
    real(r8), intent(out) :: effice(pcols,pver)            ! in cloud ice eff rad
!------------------
    real(r8) :: effliq_fn(pcols,pver)         ! in cloud liq eff rad at fixed number concentration
    real(r8) :: wsub(pcols,pver)              ! sub-grid vertical velocity (m/s)
    !! output from mmicro_pcond
    real(r8) :: tlat(pcols,pver)
    real(r8) :: qvlat(pcols,pver)
    real(r8) :: qcten(pcols,pver)
    real(r8) :: qiten(pcols,pver)
    real(r8) :: ncten(pcols,pver)
    real(r8) :: niten(pcols,pver)
    real(r8) :: effc(pcols,pver)
    real(r8) :: effc_fn(pcols,pver)     ! liquid effective radius at fixed number (for indirect calc)
    real(r8) :: effi(pcols,pver)
!++clu
	real(r8) :: allcld_ice(pcols,pver)                 ! All-cloud cloud ice
	real(r8) :: allcld_liq(pcols,pver)                 ! All-cloud liquid
	
	real(r8) :: zin(pcols,pver)
	real(r8) :: ziten(pcols,pver)
	real(r8) :: zitend(pcols,pver)

   real(r8) :: icizzz(pcols,pver)                      ! In cloud ice zzz, m6/kg
   real(r8) :: miuice(pcols,pver)                      ! miu of ice [-]  
   real(r8) :: icizzz_out(pcols,pver)                 !  
   real(r8) :: icinc_out(pcols,pver)                 !    
   real(r8) :: miuice_out(pcols,pver)                 ! 
   real(r8) :: frezzzmiuice_out(pcols,pver)                 !  occur frequency  1- 0-

	real(r8) :: Mufallni(pcols,pver)  ! number-weighted ice fall velocity, m/s 
	real(r8) :: Mufallqi(pcols,pver)  ! mass-weighted ice fall velocity, m/s
	real(r8) :: Mufallpro(pcols,pver) !
	real(r8) :: MuR3_fall(pcols,pver) !

	real(r8) :: MuR3_before_fall(pcols,pver)     
	real(r8) :: Mubeforefallpro(pcols,pver)     
	real(r8) :: MuR3_after_fall(pcols,pver)     
	real(r8) :: Muafterfallpro(pcols,pver)
	
	real(r8) :: Mufqprci(pcols,pver)  ! autoconversion of cloud ice to snow;  mass ratio, -
	real(r8) :: Mufnprci(pcols,pver)  ! autoconversion of cloud ice to snow;  number ratio, -
	real(r8) :: hxepsi(pcols,pver)
	real(r8) :: hxdvi(pcols,pver)
!--clu
    real(r8) :: prect(pcols)
    real(r8) :: preci(pcols)
    !! averaging arrays for effective radius and number....
    real(r8) :: efiout(pcols,pver)
    real(r8) :: efcout(pcols,pver)
    real(r8) :: ncout(pcols,pver)
    real(r8) :: niout(pcols,pver)
    real(r8) :: freqi(pcols,pver)
    real(r8) :: freql(pcols,pver)
    !! average cloud top radius & number
    real(r8) :: ctrel(pcols)
    real(r8) :: ctrei(pcols)
    real(r8) :: ctnl(pcols)
    real(r8) :: ctni(pcols)
    real(r8) :: fcti(pcols)
    real(r8) :: fctl(pcols)
!++sxj
    real(r8) mass_to_mmr(pcols,pver) ! conversion of layer mass to mass mixing ratio
    real(r8) aer_mmr(pcols,pver,naer_all) ! aerosol mass mixing ratio
    real(r8) aer_mass(pcols,pver,naer_all) ! aerosol mass  
    real(r8) scales(naer_all)  
!--sxj    
	!! WATER TRACERS
    integer itrip			            ! counter of water tracer triplets
    integer ixwtvap, ixwtliq, ixwtice   ! constituent indicies for vap, liq, ice of triplet
    real(r8) dum(pcols,pver)	        ! dummy argument
    !real(r8) :: wtprec_sed(pcols,pwspc) ! surface flux of total cloud water from sedimentation
    !real(r8) :: wtsnow_sed(pcols,pwspc) ! surface flux of cloud ice from sedimentation
    !real(r8) :: wtrain_sed(pcols,pwspc) ! surface flux of cloud liquid


!
	lchnk = state%lchnk
    ncol  = state%ncol


    call physics_state_copy(state,state1)   ! copy state to local state1.
    call physics_ptend_init(ptend_loc)  ! initialize local ptend type
    call physics_ptend_init(ptend_all)  ! initialize output ptend type
    call physics_tend_init(tend)        ! tend here is just a null place holder

!! Associate pointers with physics buffer fields
    itim = pbuf_old_tim_idx()
    !ifld = pbuf_get_fld_idx('QCWAT')
    !qcwat => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)
    !ifld = pbuf_get_fld_idx('TCWAT')
    !tcwat => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)
    !ifld = pbuf_get_fld_idx('LCWAT')
    !lcwat => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)

    !ifld = pbuf_get_fld_idx('CLD')
    !cld => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)
    !ifld = pbuf_get_fld_idx('CONCLD')
    !concld => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)

    !ifld = pbuf_get_fld_idx('QME')
    !qme  => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)
    !ifld = pbuf_get_fld_idx('PRAIN')
    !prain  => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)
    !ifld = pbuf_get_fld_idx('NEVAPR')
    !nevapr  => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)

    ifld = pbuf_get_fld_idx('REL')
    rel  => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)
    ifld = pbuf_get_fld_idx('REI')
    rei  => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)

    ifld = pbuf_get_fld_idx('REL_FN')
    rel_fn  => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)
    ifld = pbuf_get_fld_idx('KVH')
    kkvh => pbuf(ifld)%fld_ptr(1,1:pcols,1:pverp,lchnk, 1)
    !! if first timestep, initialize heatflux....in pbuf at all time levels.
    if (is_first_step()) then
       kkvh(:,:)= 0.0_r8
    else
       !if ( masterproc) write(6,*) "SXJ MG.F90  kkvh(2,:)",kkvh(2,:)
    endif

!++sxj  Note that, In RK, state is the state after detrainment. 
    totcw(:ncol,:pver) = state%q(:ncol,:pver,ixcldice) + state%q(:ncol,:pver,ixcldliq)
    rdtime = 1._r8/dtime 
    qtend(:ncol,:pver) = (state%q(:ncol,:pver,1) - qcwato(:ncol,:pver))*rdtime
    ttend(:ncol,:pver) = (state%t(:ncol,:pver)   - tcwato(:ncol,:pver))*rdtime
    ltend(:ncol,:pver) = (totcw  (:ncol,:pver)   - lcwato(:ncol,:pver))*rdtime
!--sxj 

	! ++detrain ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Put the detraining cloud water from convection into the cloud and environment.
    ptend_loc%name  = 'pcwdetrain'
	ptend_loc%lq(ixcldliq) = .TRUE.
    ptend_loc%lq(ixcldice) = .TRUE.
    ptend_loc%lq(ixnumliq) = .TRUE.
    ptend_loc%lq(ixnumice) = .TRUE.
!	ptend_loc%lq(ixzzzice) = .TRUE. 		!clu
    ptend_loc%ls           = .TRUE.
    ! Put all of the detraining cloud water from convection into the large scale cloud.
    ! put detraining cloud water into liq and ice based on temperature partition
    do k = 1,pver
       do i = 1,state1%ncol
          if (state1%t(i,k) > 268.15_r8) then
             dum1=0.0_r8
          else if (state1%t(i,k) < 238.15_r8) then
             dum1=1.0_r8
			 if (state1%t(i,k) < 50._r8) then   ! sxj add
			    write(*,*) "SXJ state1%t(i,k) < 50._r8"
				call endrun
			 endif
          else
             dum1=(268.15_r8-state1%t(i,k))/30._r8
          end if
          ptend_loc%q(i,k,ixcldliq) = dlf(i,k)*(1._r8-dum1)
          ptend_loc%q(i,k,ixcldice) = dlf(i,k)*dum1
          ! calculate detrainment of Nc  ! 'dif'------detrain water ---sxj---
          ! assume detrained cloud water has mean volume radius of 8 micron
          dum2 = dlf(i,k)*(1._r8-dum1)
   ptend_loc%q(i,k,ixnumliq) = 3._r8*dum2/(4._r8*3.14_r8*5.12e-16_r8* 997._r8) !cntl run rd=8um
!liuym   ptend_loc%q(i,k,ixnumliq) = 3._r8*dum2/(4._r8*3.14_r8*4.096e-15_r8* 997._r8)    !rd=16um
!liuym   ptend_loc%q(i,k,ixnumliq) = 3._r8*dum2/(4._r8*3.14_r8*3.2768e-14_r8* 997._r8)    !rd=32um
!liuym   ptend_loc%q(i,k,ixnumliq) = 3._r8*dum2/(4._r8*3.14_r8*2.62144e-13_r8* 997._r8)    !rd=64um
!!liuym  ptend_loc%q(i,k,ixnumliq) = 3._r8*dum2/(4._r8*3.14_r8*2.16e-16_r8* 997._r8)        !rd=6um
!!liuym  ptend_loc%q(i,k,ixnumliq) = 3._r8*dum2/(4._r8*3.14_r8*5.12e-16_r8* 997._r8)*0.296  ! rd=12um
!!liuym  ptend_loc%q(i,k,ixnumliq) = 3._r8*dum2/(4._r8*3.14_r8*5.12e-16_r8* 997._r8)*0.187  ! rd=14um
!!liuym  ptend_loc%q(i,k,ixnumliq) = 3._r8*dum2/(4._r8*3.14_r8*5.12e-16_r8* 997._r8)*0.125  ! cntl,rd=16um
!!liuym  ptend_loc%q(i,k,ixnumliq) = 3._r8*dum2/(4._r8*3.14_r8*5.12e-16_r8* 997._r8)*0.064 ! rd=20um
          dum2 = dlf(i,k)*dum1
          if (state1%t(i,k) < 233.15_r8) then !!! maybe wrong   233.15_r8 should accord with the T 238.15 above ---sxj---
 ptend_loc%q(i,k,ixnumice) = 3._r8*dum2/(4._r8*3.14_r8*1.563e-14_r8* 500._r8)    !cntl,rdice=25um
!liuym ptend_loc%q(i,k,ixnumice) = 3._r8*dum2/(4._r8*3.14_r8*1.0648e-14_r8* 500._r8)       !rdice=22um
!liuym   ptend_loc%q(i,k,ixnumice) = 3._r8*dum2/(4._r8*3.14_r8*1.728e-15_r8* 500._r8)       !rdice=12um
!liuym   ptend_loc%q(i,k,ixnumice) = 3._r8*dum2/(4._r8*3.14_r8*8.0e-15_r8* 500._r8)       !rdice=20um
!liuym   ptend_loc%q(i,k,ixnumice) = 3._r8*dum2/(4._r8*3.14_r8*1.25e-13_r8* 500._r8)       !rdice=50um
!liuym   ptend_loc%q(i,k,ixnumice) = 3._r8*dum2/(4._r8*3.14_r8*4.21875e-13_r8* 500._r8)       !rdice=75um
!!liuym  ptend_loc%q(i,k,ixnumice) = 3._r8*dum2/(4._r8*3.14_r8*1.563e-14_r8* 500._r8)*1.3 ! control run rdice=29um
!!liuym ptend_loc%q(i,k,ixnumice) = 3._r8*dum2/(4._r8*3.14_r8*1.563e-14_r8* 500._r8)*1.7   !rdice=27um
!!liuym ptend_loc%q(i,k,ixnumice) = 3._r8*dum2/(4._r8*3.14_r8*1.563e-14_r8* 500._r8)*2.0   !rdice=25um
!!   ptend_loc%q(i,k,ixnumice) = 3._r8*dum2/(4._r8*3.14_r8*1.563e-14_r8* 500._r8)*1.7   !  ! sxj ssume detrained ice crystal has mean volume radius of 30 micron 
!++clu
!   ptend_loc%q(i,k,ixzzzice) = 3._r8*dum2/(4._r8*3.14_r8*25.e-6_r8**3*500._r8) * 50.e-6_r8**6+ & ! Deep    Convection
!							   3._r8*dum2/(4._r8*3.14_r8*50.e-6_r8**3*500._r8) * 100.e-6_r8**6   ! Shallow Convection
!--clu
          endif
          ! account for latent heat release during freezing
          ptend_loc%s(i,k) = dlf(i,k)*dum1*latice
       end do
    end do
    call outfld('ZMDLF' ,dlf  , pcols,state1%lchnk)
    ! add tendency from this process to tend from other processes here
    call physics_ptend_sum(ptend_loc, ptend_all, state)
    call physics_update(state1, tend, ptend_loc, dtime)
    call physics_ptend_init(ptend_loc)

!++++ cldfrc ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! cloud fraction after transport and convection,
	! derive the relationship between rh and cld from
    ! the employed c7loud scheme

    pmid(:ncol,:pver) = state1%pmid(:ncol,:pver)
    pint(:ncol,:pverp)= state1%pint(:ncol,:pverp)
    t(:ncol,:pver) = state1%t(:ncol,:pver)
    q(:ncol,:pver) = state1%q(:ncol,:pver,1)
    omga(:ncol,:pver) = state1%omega(:ncol,:pver)
    pdel(:ncol,:pver) = state1%pdel(:ncol,:pver)
    zm(:ncol,:pver) =  state1%zm(:ncol,:pver)  ! sxj 
    ps(:ncol)   = state1%pint(:ncol,pverp)
    phis(:ncol) = state1%phis(:ncol)

    call t_startf("cldfrc")

    qc(:ncol,:pver) = state1%q(:ncol,:pver,ixcldliq)
    qi(:ncol,:pver) = state1%q(:ncol,:pver,ixcldice)
  
    call cldfrc(lchnk,   ncol,    tdeepcld, qdeepcld, deepcld,  tshallowcld, &
               qshallowcld, shallowcld,    &
               pmid, pint, dtime ,dlf,     t,        q,     qi,   qc, omga, phis, &
               cnt,     cnb,     lcl   , cld,    clc,     pdel,   &
               cmfmc,   cmfmc2,  landfrac,snowh,   concld,  cldst,    &
               cld_env  ,t_env   ,q_env     , &
               ts,      sst, ps,      zdu,     ocnfrac, rhu00, &
               relhum,  0, zm, cstcld=cstcld  )


    call cldfrc(lchnk,   ncol,   tdeepcld, qdeepcld, deepcld,  tshallowcld, &
               qshallowcld, shallowcld,    &
               pmid, pint, dtime ,dlf,     t,        q,     qi,   qc, omga, phis, &
               cnt,     cnb,     lcl,  cld2,   clc,     pdel,   &
               cmfmc,   cmfmc2   ,landfrac,snowh,   concld2, cldst2,   &
               cld_env  ,t_env   ,q_env     , &
               ts,      sst, ps,        zdu,   ocnfrac, rhu002, &
               relhum2, 1, zm  )

    call t_stopf("cldfrc")
    ! cldfrc does not define layer cloud for model layer at k=1
    ! so set rhu00(k=1)=2.0 to not calculate cme for this layer
    rhu00(:ncol,1) = 2.0_r8
    ! Add following to estimate rhdfda
    do k=1,pver
    do i=1,ncol
       if (relhum(i,k) < rhu00(i,k) ) then
          rhdfda(i,k)=0.0_r8
       else if (relhum(i,k) >= 1.0_r8 ) then
          rhdfda(i,k)=0.0_r8
       else
          if((cld2(i,k) - cld(i,k) ) < 1.e-4_r8 ) then
             rhdfda(i,k) = 0.01_r8*relhum(i,k)*1.e+4_r8   !instead of 0.0
          else
             rhdfda(i,k)=0.01_r8*relhum(i,k)/(cld2(i,k)-cld(i,k))
          endif
       endif
    enddo
    enddo


!+ mp +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   ! cloud water and ice parameterizations
   call t_startf('stratiform_microphys')
   ! Initialize chunk id and size
   lchnk = state1%lchnk
   ncol  = state1%ncol
   rdtime = 1._r8/dtime
   ! Define fractional amount of cloud condensate in ice phase
   call cldwat_fice(ncol, state1%t, fice, fsnow)
   !tmax_fice  = tmelt - 10._r8       ! max temperature for cloud ice formation
   !tmin_fice  = tmax_fice - 30._r8   ! min temperature for cloud ice formation
   !tmax_fsnow = tmelt                ! max temperature for transition to convective snow
   !tmin_fsnow = tmelt-5._r8          ! min temperature for transition to convective snow
   ptend_loc%name         = 'cldwat2m_tend'  
   ptend_loc%ls           = .true.
   ptend_loc%lq(1)        = .true.
   ptend_loc%lq(ixcldice) = .true.
   ptend_loc%lq(ixcldliq) = .true.
   ptend_loc%lq(ixnumliq) = .true.
   ptend_loc%lq(ixnumice) = .true.
!   ptend_loc%lq(ixzzzice) = .true.  	!+clu
   qc(:ncol,:pver) = state1%q(:ncol,:pver,ixcldliq)
   qi(:ncol,:pver) = state1%q(:ncol,:pver,ixcldice)
   nc(:ncol,:pver) = state1%q(:ncol,:pver,ixnumliq)
   ni(:ncol,:pver) = state1%q(:ncol,:pver,ixnumice)
!   zi(:ncol,:pver) = state1%q(:ncol,:pver,ixzzzice)		!+clu

!++sxj===== read aerosol data   ==========================================
     aer_mmr(:,:,:) = 0.
     aer_mass(:,:,:) = 0.

     call get_int_scales(scales)
#ifdef MAC_SP
     call get_aerosol_cdnc(lchnk, state%pint, aer_mass, scales, state)   ! sxj
#else
     call get_aerosol(lchnk, state%pint, aer_mass, scales, state)   ! sxj
#endif

     mass_to_mmr(:ncol,:)=gravit*state1%rpdel(:ncol,:)
	 do m=1,naer_all
        !aer_mmr(:ncol,:,m)=aer_mass(:ncol,:,m)*mass_to_mmr(:ncol,:)
		aer_mmr(:ncol,:,m)=aer_mass(:ncol,:,m)
     enddo
	 call outfld('aer_data',aer_mmr(1:pcols,1:pver,1), pcols,lchnk)    !!sxj 

	 !if ( masterproc) write(6,*) "SXJ stratiform_tend aer_mmr(5,:,1) ",aer_mmr(5,:,1) 

!--sxj=====================================================================
			  
!      call t_startf('mmicro_pcond')
	  
!++clu  init output vars 
   ziten = 0._r8
!   miuice= 0._r8
   icizzz = 0._r8
   icizzz_out = 0._r8
   icinc_out = 0._r8
   miuice_out = 0._r8
   frezzzmiuice_out = 0._r8
!--clu
	   call t_startf('mmicro_pcond')
		
      call mmicro_pcond (lchnk,   ncol,dtime,&
          state1%t  , ttend      , state1%q(1,1,1), qtend   , ltend,  qc,qi , &
          ls_flxprc,  ls_flxsnw,                                              &    ! wtw
	  nc,ni, state1%pmid , state1%pdel , cld, aer_mmr, rhdfda, rhu00, fice, &  
          tlat,qvlat,qcten,qiten,ncten,niten,effc,effc_fn,effi, &
          prect,preci,kkvh,wsub,zi,ziten,miuice,nevapr,evapsnow,prain,prodsnow,qme,landfrac,snowh,omga ,&  ! sxj add omega For BF process
		  Mufallni,Mufallqi,MuR3_fall,Mufallpro,MuR3_before_fall,Mubeforefallpro,MuR3_after_fall,Muafterfallpro,  &  !!clu_hx
		  Mufqprci,Mufnprci,hxepsi,hxdvi)
!clu+ zi, ziten, miuice ,Mu* and hxepsi  
      call t_stopf('mmicro_pcond')

      do i = 1,ncol
         do k = 1,pver
            ptend_loc%s(i,k)          =tlat(i,k)
            ptend_loc%q(i,k,1)        =qvlat(i,k)
            ptend_loc%q(i,k,ixcldliq) =qcten(i,k)
            ptend_loc%q(i,k,ixcldice) =qiten(i,k)
            ptend_loc%q(i,k,ixnumliq) =ncten(i,k)
            ptend_loc%q(i,k,ixnumice) =niten(i,k)
!			ptend_loc%q(i,k,ixzzzice) =ziten(i,k)
         end do
      end do


      prec_pcw(:ncol) = prect(:ncol)
      snow_pcw(:ncol) = preci(:ncol)
      prec_sed(:ncol) =0._r8              ! has been considered in prect
      snow_sed(:ncol) =0._r8              ! has been considered in preci
      !prec_str(:ncol) = prec_pcw(:ncol)+prec_sed(:ncol)-rliq(:ncol)  ! sxj     rliq has been considered before stratiform_tend
  	  prec_str(:ncol) = prec_pcw(:ncol)+prec_sed(:ncol)             ! sxj   
      snow_str(:ncol) = snow_pcw(:ncol)+snow_sed(:ncol)
	  

      call physics_ptend_sum(ptend_loc,ptend_all, state)     !!!


   ! Set the name of the final package tendencies. Note that this
   ! is a special case in physics_update, so a change here must be
   ! matched there.
      ptend_all%name = 'stratiform'    ! sxj ???
   ! used below


      call physics_update (state1, tend, ptend_loc, dtime)


      call physics_ptend_init(ptend_loc)  ! sxj  ptend_loc has been reset in  physics_update


!      accumulate prec and snow
      ! Save off q and t after cloud water
      do k=1,pver
        !qcwat(:ncol,k) = state1%q(:ncol,k,1)
        !tcwat(:ncol,k) = state1%t(:ncol,k)
        !lcwat(:ncol,k) = state1%q(:ncol,k,ixcldliq) + state1%q(:ncol,k,ixcldice)
         rel(:ncol,k) = effc(:ncol,k)
         rel_fn(:ncol,k) = effc_fn(:ncol,k)
         rei(:ncol,k) = effi(:ncol,k)
      end do

      ls_reffrain(:,:) = effc(:,:)    ! wtw
      ls_reffsnow(:,:) = effi(:,:)    ! wtw

      ! Compute in cloud ice and liquid mixing ratios (output only)
      do k=1,pver
         do i = 1,ncol
            icimr(i,k) = (state1%q(i,k,ixcldice) + dtime*ptend_loc%q(i,k,ixcldice)) /&
               max(0.01_r8,cld(i,k))
            icwmr(i,k) = (state1%q(i,k,ixcldliq) + dtime*ptend_loc%q(i,k,ixcldliq)) /&
               max(0.01_r8,cld(i,k))
            icinc(i,k) = (state1%q(i,k,ixnumice) + dtime*ptend_loc%q(i,k,ixnumice)) /&
               max(0.01_r8,cld(i,k))*state1%pmid(i,k)/(287.15_r8*state1%t(i,k))
            icwnc(i,k) = (state1%q(i,k,ixnumliq) + dtime*ptend_loc%q(i,k,ixnumliq)) /&
               max(0.01_r8,cld(i,k))*state1%pmid(i,k)/(287.15_r8*state1%t(i,k))
!++clu
!			icizzz(i,k) = (state1%q(i,k,ixzzzice)+dtime*ptend_loc%q(i,k,ixzzzice))  /&
!			   max(0.01_r8,cld(i,k))*state1%pmid(i,k)/(287.15_r8*state1%t(i,k))
!--clu			   
            effliq(i,k) = effc(i,k)
            effliq_fn(i,k) = effc_fn(i,k)
            effice(i,k) = effi(i,k)
         end do
      end do


   !--------------------- OUTPUT FIELDS FOR HISTORY ---------------------

      ! Column droplet concentration
      do i = 1,ncol
         cdnumc(i)=0._r8
         do k=1,pver
            cdnumc(i)=cdnumc(i)+state1%q(i,k,ixnumliq)*state1%pdel(i,k)/9.816_r8
         end do
      end do

      ! Averaging for new output fields
      efcout(:,:)=0._r8
      efiout(:,:)=0._r8
      ncout(:,:)=0._r8
      niout(:,:)=0._r8
      freql(:,:)=0._r8
      freqi(:,:)=0._r8
      do i = 1,ncol
         do k=1,pver
            if (cld(i,k).gt.0.01.and.icwmr(i,k).gt.1.e-4_r8) then
               efcout(i,k)=effc(i,k)
               ncout(i,k)=icwnc(i,k)
               freql(i,k)=1._r8
            endif
            if (cld(i,k).gt.0.01.and.icimr(i,k).gt.1.e-6_r8) then
               efiout(i,k)=effi(i,k)
               niout(i,k)=icinc(i,k)
               freqi(i,k)=1._r8
!++clu
			   frezzzmiuice_out(i,k)= 1.0_r8
			   icizzz_out(i,k) = icizzz(i,k)
!               miuice_out(i,k) = miuice(i,k)
               icinc_out(i,k)  = icinc(i,k)
!--clu
            endif
         end do
      end do

      !Add averaged out fields.
      call outfld('AREL',efcout,    pcols,lchnk)
      call outfld('AREI',efiout,    pcols,lchnk)
      call outfld('AWNC',ncout,    pcols,lchnk)
      call outfld('AWNI',niout,    pcols,lchnk)
      call outfld('FREQL',freql,    pcols,lchnk)
      call outfld('FREQI',freqi,    pcols,lchnk)

      !Cloud top effective radius and number....
      fcti(:)=0._r8
      fctl(:)=0._r8
      ctrel(:)=0._r8
      ctrei(:)=0._r8
      ctnl(:)=0._r8
      ctni(:)=0._r8
      do i = 1,ncol
         do k=1,pver
            if (cld(i,k).gt.0.01.and.icwmr(i,k).gt.1.e-7_r8) then
               ctrel(i)=effc(i,k)
               ctnl(i)=icwnc(i,k)
               fctl(i)=1._r8
               exit
            endif

            if (cld(i,k).gt.0.01.and.icimr(i,k).gt.1.e-7_r8) then
               ctrei(i)=effi(i,k)
               ctni(i)=icinc(i,k)
               fcti(i)=1._r8
               exit
            endif

         enddo
      enddo

      call outfld('ACTREL',ctrel,    pcols,lchnk)
      call outfld('ACTREI',ctrei,    pcols,lchnk)
      call outfld('ACTNL',ctnl,    pcols,lchnk)
      call outfld('ACTNI',ctni,    pcols,lchnk)
      call outfld('FCTL',fctl,    pcols,lchnk)
      call outfld('FCTI',fcti,    pcols,lchnk)

      ! microphysics variables to output fields
      call outfld('MPDT',tlat,    pcols,lchnk)
      call outfld('MPDQ',qvlat,    pcols,lchnk)
      call outfld('ICINC',icinc,    pcols,lchnk)
      call outfld('ICWNC',icwnc,    pcols,lchnk)
      call outfld('EFFLIQ',effliq,    pcols,lchnk)
      call outfld('EFFLIQ_IND',effliq_fn,    pcols,lchnk)
      call outfld('EFFICE',effice,    pcols,lchnk)
      call outfld('WSUB',wsub,    pcols,lchnk)
      call outfld('CDNUMC',cdnumc,    pcols,lchnk)

     call outfld('ICIMR',icimr,    pcols,lchnk)
     call outfld('ICWMR',icwmr,    pcols,lchnk)
     call outfld('CME'  ,qme     , pcols,lchnk)
     call outfld('PRODPREC',prain, pcols,lchnk)
     call outfld('EVAPPREC',nevapr, pcols,lchnk)
     !call outfld('EVAPSNOW',evapsnow, pcols,lchnk)

!++clu
	call outfld( 'FREZICE',frezzzmiuice_out,   pcols, lchnk )
!	call outfld( 'MIUICE', miuice_out,   pcols, lchnk )
	call outfld( 'MIUICE', miuice,   pcols, lchnk )
	call outfld( 'ZICE', icizzz_out,   pcols, lchnk )
	call outfld( 'NICE', icinc_out,   pcols, lchnk )
!	call outfld( 'QICE', icimrst_out,   pcols, lchnk )
	
	call outfld('Mufallqi',    Mufallqi,     pcols, lchnk )
	call outfld('Mufallni',    Mufallni,     pcols, lchnk )
	call outfld('Mufallpro',   Mufallpro,    pcols, lchnk )
	call outfld('MuR3_fall',   MuR3_fall,    pcols, lchnk )

	call outfld('MuR3_before_fall',  MuR3_before_fall,   pcols, lchnk )
	call outfld('Mubeforefallpro',   Mubeforefallpro,    pcols, lchnk )
	call outfld('MuR3_after_fall',   MuR3_after_fall,    pcols, lchnk )
	call outfld('Muafterfallpro',   Muafterfallpro,     pcols, lchnk )

	call outfld('Mufqprci',    Mufqprci,     pcols, lchnk )
	call outfld('Mufnprci',    Mufnprci,     pcols, lchnk )
	call outfld('hxepsi',      hxepsi,       pcols, lchnk )
	call outfld('hxdvi',        hxdvi,       pcols, lchnk )
	
!	call outfld('ICECLDF',      icecldf,       pcols, lchnk )
!	call outfld('LIQCLDF',      liqcldf,       pcols, lchnk )
!	call outfld('RHCLOUD',     rhcloud,       pcols, lchnk )
!--clu
      call t_stopf('stratiform_microphys')

	  
	  !if (masterproc) write(6,*) "SXJ stratiform_tend done " 

   endsubroutine stratiform_tend

#endif
!--sxj

end module cldcond

