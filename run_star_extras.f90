! ***********************************************************************
!
!   Copyright (C) 2011  The MESA Team
!
!   this file is part of mesa.
!
!   mesa is free software; you can redistribute it and/or modify
!   it under the terms of the gnu general library public license as published
!   by the free software foundation; either version 2 of the license, or
!   (at your option) any later version.
!
!   mesa is distributed in the hope that it will be useful,
!   but without any warranty; without even the implied warranty of
!   merchantability or fitness for a particular purpose.  see the
!   gnu library general public license for more details.
!
!   you should have received a copy of the gnu library general public license
!   along with this software; if not, write to the free software
!   foundation, inc., 59 temple place, suite 330, boston, ma 02111-1307 usa
!
! ***********************************************************************

module run_star_extras

   use custom_eos

   use star_lib
   use star_def
   use const_def
   use math_lib

   implicit none

contains

   subroutine extras_controls(id, ierr)
      integer, intent(in) :: id
      integer, intent(out) :: ierr
      type (star_info), pointer :: s
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return

      ! this is the place to set any procedure pointers you want to change
      ! e.g., other_wind, other_mixing, other_energy  (see star_data.inc)


      ! the extras functions in this file will not be called
      ! unless you set their function pointers as done below.
      ! otherwise we use a null_ version which does nothing (except warn).

      s% extras_startup => extras_startup
      s% extras_start_step => extras_start_step
      s% extras_check_model => extras_check_model
      s% extras_finish_step => extras_finish_step
      s% extras_after_evolve => extras_after_evolve
      s% how_many_extra_history_columns => how_many_extra_history_columns
      s% data_for_extra_history_columns => data_for_extra_history_columns
      s% how_many_extra_profile_columns => how_many_extra_profile_columns
      s% data_for_extra_profile_columns => data_for_extra_profile_columns

      s% how_many_extra_history_header_items => how_many_extra_history_header_items
      s% data_for_extra_history_header_items => data_for_extra_history_header_items
      s% how_many_extra_profile_header_items => how_many_extra_profile_header_items
      s% data_for_extra_profile_header_items => data_for_extra_profile_header_items

      ! EoS
      
      ! initialize the EoS
      call eos_init_custom_eos()
      write(*,*) "extras_controls: eos initialized"

      s% eos_rq % other_eos_frac => my_other_eos_frac
      s% eos_rq % other_eos_component => my_other_eos_component

   end subroutine extras_controls


   subroutine extras_startup(id, restart, ierr)
      integer, intent(in) :: id
      logical, intent(in) :: restart
      integer, intent(out) :: ierr
      type (star_info), pointer :: s
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
   end subroutine extras_startup


   integer function extras_start_step(id)
      integer, intent(in) :: id
      integer :: ierr
      type (star_info), pointer :: s
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      extras_start_step = 0
   end function extras_start_step


   ! returns either keep_going, retry, or terminate.
   integer function extras_check_model(id)
      integer, intent(in) :: id
      integer :: ierr
      type (star_info), pointer :: s
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      extras_check_model = keep_going
      if (.false. .and. s% star_mass_h1 < 0.35d0) then
         ! stop when star hydrogen mass drops to specified level
         extras_check_model = terminate
         write(*, *) 'have reached desired hydrogen mass'
         return
      end if


      ! if you want to check multiple conditions, it can be useful
      ! to set a different termination code depending on which
      ! condition was triggered.  MESA provides 9 customizeable
      ! termination codes, named t_xtra1 .. t_xtra9.  You can
      ! customize the messages that will be printed upon exit by
      ! setting the corresponding termination_code_str value.
      ! termination_code_str(t_xtra1) = 'my termination condition'

      ! by default, indicate where (in the code) MESA terminated
      if (extras_check_model == terminate) s% termination_code = t_extras_check_model
   end function extras_check_model


   integer function how_many_extra_history_columns(id)
      integer, intent(in) :: id
      integer :: ierr
      type (star_info), pointer :: s
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      how_many_extra_history_columns = 0
   end function how_many_extra_history_columns


   subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
      integer, intent(in) :: id, n
      character (len=maxlen_history_column_name) :: names(n)
      real(dp) :: vals(n)
      integer, intent(out) :: ierr
      type (star_info), pointer :: s
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return

      ! note: do NOT add the extras names to history_columns.list
      ! the history_columns.list is only for the built-in history column options.
      ! it must not include the new column names you are adding here.


   end subroutine data_for_extra_history_columns


   integer function how_many_extra_profile_columns(id)
      integer, intent(in) :: id
      integer :: ierr
      type (star_info), pointer :: s
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      how_many_extra_profile_columns = 0
   end function how_many_extra_profile_columns


   subroutine data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
      integer, intent(in) :: id, n, nz
      character (len=maxlen_profile_column_name) :: names(n)
      real(dp) :: vals(nz,n)
      integer, intent(out) :: ierr
      type (star_info), pointer :: s
      integer :: k
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return

      ! note: do NOT add the extra names to profile_columns.list
      ! the profile_columns.list is only for the built-in profile column options.
      ! it must not include the new column names you are adding here.

      ! here is an example for adding a profile column
      !if (n /= 1) stop 'data_for_extra_profile_columns'
      !names(1) = 'beta'
      !do k = 1, nz
      !   vals(k,1) = s% Pgas(k)/s% P(k)
      !end do

   end subroutine data_for_extra_profile_columns


   integer function how_many_extra_history_header_items(id)
      integer, intent(in) :: id
      integer :: ierr
      type (star_info), pointer :: s
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      how_many_extra_history_header_items = 0
   end function how_many_extra_history_header_items


   subroutine data_for_extra_history_header_items(id, n, names, vals, ierr)
      integer, intent(in) :: id, n
      character (len=maxlen_history_column_name) :: names(n)
      real(dp) :: vals(n)
      type(star_info), pointer :: s
      integer, intent(out) :: ierr
      ierr = 0
      call star_ptr(id,s,ierr)
      if(ierr/=0) return

      ! here is an example for adding an extra history header item
      ! also set how_many_extra_history_header_items
      ! names(1) = 'mixing_length_alpha'
      ! vals(1) = s% mixing_length_alpha

   end subroutine data_for_extra_history_header_items


   integer function how_many_extra_profile_header_items(id)
      integer, intent(in) :: id
      integer :: ierr
      type (star_info), pointer :: s
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      how_many_extra_profile_header_items = 0
   end function how_many_extra_profile_header_items


   subroutine data_for_extra_profile_header_items(id, n, names, vals, ierr)
      integer, intent(in) :: id, n
      character (len=maxlen_profile_column_name) :: names(n)
      real(dp) :: vals(n)
      type(star_info), pointer :: s
      integer, intent(out) :: ierr
      ierr = 0
      call star_ptr(id,s,ierr)
      if(ierr/=0) return

      ! here is an example for adding an extra profile header item
      ! also set how_many_extra_profile_header_items
      ! names(1) = 'mixing_length_alpha'
      ! vals(1) = s% mixing_length_alpha

   end subroutine data_for_extra_profile_header_items


   ! returns either keep_going or terminate.
   ! note: cannot request retry; extras_check_model can do that.
   integer function extras_finish_step(id)
      integer, intent(in) :: id
      integer :: ierr
      type (star_info), pointer :: s
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      extras_finish_step = keep_going

      ! to save a profile,
      ! s% need_to_save_profiles_now = .true.
      ! to update the star log,
      ! s% need_to_update_history_now = .true.

      ! see extras_check_model for information about custom termination codes
      ! by default, indicate where (in the code) MESA terminated
      if (extras_finish_step == terminate) s% termination_code = t_extras_finish_step
   end function extras_finish_step


   subroutine extras_after_evolve(id, ierr)

      integer, intent(in) :: id
      integer, intent(out) :: ierr
      type (star_info), pointer :: s

      ierr = 0
      call star_ptr(id, s, ierr)

      if (ierr /= 0) return

   end subroutine extras_after_evolve

   ! Equation of state functions

   subroutine my_other_eos_component( &
      handle, &
      species, chem_id, net_iso, xa, &
      Rho, log10Rho, T, log10T, &
      res, d_dlnRho_const_T, d_dlnT_const_Rho, d_dxa_const_TRho, ierr)

      ! INPUT
      use chem_def, only: num_chem_isos
      use chem_lib, only: basic_composition_info
      use eos_def, only: eos_handles


      integer, intent(in) :: handle ! eos handle; from star, pass s% eos_handle

      integer, intent(in) :: species ! number of species
      integer, pointer :: chem_id(:) ! maps species to chem id
      ! index from 1 to species
      ! value is between 1 and num_chem_isos
      integer, pointer :: net_iso(:) ! maps chem id to species number
      ! index from 1 to num_chem_isos (defined in chem_def)
      ! value is 0 if the iso is not in the current net
      ! else is value between 1 and number of species in current net
      real(dp), intent(in) :: xa(:) ! mass fractions

      real(dp), intent(in) :: Rho, log10Rho ! the density
      real(dp), intent(in) :: T, log10T ! the temperature

      ! OUTPUT

      real(dp), intent(inout) :: res(:) ! (num_eos_basic_results)
      ! partial derivatives of the basic results wrt lnd and lnT
      real(dp), intent(inout) :: d_dlnRho_const_T(:) ! (num_eos_basic_results)
      ! d_dlnRho(i) = d(res(i))/dlnd|T
      real(dp), intent(inout) :: d_dlnT_const_Rho(:) ! (num_eos_basic_results)
      ! d_dlnT(i) = d(res(i))/dlnT|Rho
      real(dp), intent(inout) :: d_dxa_const_TRho(:,:) ! (num_eos_basic_results, species)
      ! d_dxa(i,j) = d(res(i))/dxa(j)|T,Rho

      integer, intent(out) :: ierr ! 0 means AOK.

      ! local variables
      real(dp) :: X, Y, Z, abar, zbar, z2bar, z53bar, ye, mass_correction, sumx
      integer :: which_eosdt

      ! for accessing extra controls
      type (EoS_General_Info), pointer :: rq

      ! get the composition info
      call basic_composition_info( &
         species, chem_id, xa, X, Y, Z, &
         abar, zbar, z2bar, z53bar, ye, mass_correction, sumx)


      rq => eos_handles(handle)


      if (.not. (rq % eos_integer_ctrl(1) == 1)) then
         write(*,*) "my_other_eos_component: eos_integer_ctrl(1) is not set correctly."
         write(*,*) "Currently, only eos_integer_ctrl(1) = 1 (user-defined eos) is supported."
         ierr = 1
         return
      end if

      which_eosdt = rq % eos_integer_ctrl(1)

      call get1_for_eosdt( &
         handle, which_eosdt, Z, X, abar, zbar, &
         species, chem_id, net_iso, xa, &
         log10Rho, log10T, mass_correction, &
         res, d_dlnRho_const_T, d_dlnT_const_Rho, d_dxa_const_TRho, &
         ierr)

      ! zero phase information
      res(i_phase:i_latent_ddlnRho) = 0d0
      d_dlnT_const_Rho(i_phase:i_latent_ddlnRho) = 0d0
      d_dlnRho_const_T(i_phase:i_latent_ddlnRho) = 0d0

      ! zero all components
      res(i_frac:i_frac+num_eos_frac_results-1) = 0.0
      d_dlnRho_const_T(i_frac:i_frac+num_eos_frac_results-1) = 0.0
      d_dlnT_const_Rho(i_frac:i_frac+num_eos_frac_results-1) = 0.0

   end subroutine my_other_eos_component

   subroutine my_other_eos_frac( &
      handle, &
      species, chem_id, net_iso, xa, &
      Rho, log10Rho, T, log10T, &
      frac, dfrac_dlogRho, dfrac_dlogT, ierr)

      ! INPUT
      use chem_def, only: num_chem_isos
      use auto_diff

      integer, intent(in) :: handle ! eos handle; from star, pass s% eos_handle

      integer, intent(in) :: species ! number of species
      integer, pointer :: chem_id(:) ! maps species to chem id
      ! index from 1 to species
      ! value is between 1 and num_chem_isos
      integer, pointer :: net_iso(:) ! maps chem id to species number
      ! index from 1 to num_chem_isos (defined in chem_def)
      ! value is 0 if the iso is not in the current net
      ! else is value between 1 and number of species in current net
      real(dp), intent(in) :: xa(:) ! mass fractions

      real(dp), intent(in) :: Rho, log10Rho ! the density
      real(dp), intent(in) :: T, log10T ! the temperature

      ! OUTPUT
      ! this routine must provide a fraction (in [0,1]) of the 'other' eos to use
      ! the remaining fraction (1-frac) will be provided by the standard MESA eos
      real(dp), intent(out) :: frac ! fraction of other_eos to use
      real(dp), intent(out) :: dfrac_dlogRho ! its partial derivative at constant T
      real(dp), intent(out) :: dfrac_dlogT   ! its partial derivative at constant Rho

      integer, intent(out) :: ierr ! 0 means AOK.

      ! local variables

      ! for accessing extra controls
      type (EoS_General_Info), pointer :: rq
      integer :: which_eosdt
      
      ! for accessing the EoS tables
      type (EosDT_XZ_Info), pointer :: ep
      real(dp) :: log10Q
      logical :: flag

      ! for blending
      type(auto_diff_real_2var_order1) :: logT_auto, logRho_auto, logQ_auto
      type(auto_diff_real_2var_order1) :: blend, blend_logT, blend_logRho, blend_logQ
      real(dp) :: logT_min_for_any_custom_eos, logT_min_for_all_custom_eos, logT_max_for_any_custom_eos, logT_max_for_all_custom_eos
      real(dp) :: logQ_min_for_any_custom_eos, logQ_min_for_all_custom_eos, logQ_max_for_any_custom_eos, logQ_max_for_all_custom_eos
      real(dp) :: logRho_min_for_any_custom_eos, logRho_min_for_all_custom_eos, logRho_max_for_any_custom_eos, logRho_max_for_all_custom_eos

      include 'formats'

      ierr = 0
      
      ! logRho is val1
      logRho_auto% val = log10Rho
      logRho_auto% d1val1 = 1d0
      logRho_auto% d1val2 = 0d0

      ! logT is val2
      logT_auto% val = log10T
      logT_auto% d1val1 = 0d0
      logT_auto% d1val2 = 1d0

      logQ_auto = logRho_auto - 2d0*logT_auto + 12d0

      ! for blend variables 1 is the custom EoS, 0 is other
      ! (this is the opposite of the final frac)

      ! 1. which eos to use
      rq => eos_handles(handle)
      which_eosdt = rq % eos_integer_ctrl(1) ! see my_other_eos_component for details
      
      ! Note that this code assumes that all iX, iZ of an EoS have the same
      ! logT and logQ bounds. Hence, we load (iX, iZ) = (1, 1) here.
      call load_single_eosDT_table_by_id(which_eosdt, ep, 1, 1, ierr)
      if (ierr /= 0) return

      ! These values need to be adjusted according to your EoS
      logT_min_for_all_custom_eos = ep% logT_min + 0.11d0
      logT_min_for_any_custom_eos = ep% logT_min + 0.10d0
      logT_max_for_all_custom_eos = ep% logT_max - 0.4d0
      logT_max_for_any_custom_eos = ep% logT_max - 0.3d0
      
      logQ_min_for_all_custom_eos = ep% logQ_min
      logQ_min_for_any_custom_eos = ep% logQ_min + 0.1d0
      logQ_max_for_all_custom_eos = ep% logQ_max - 0.1d0
      logQ_max_for_any_custom_eos = ep% logQ_max

      logRho_min_for_all_custom_eos = ep% logQ_min + 2d0 * ep% logT_min - 12d0
      logRho_min_for_any_custom_eos = logRho_min_for_all_custom_eos + 0.1d0
      logRho_max_for_all_custom_eos = ep% logQ_max + 2d0 * ep% logT_max - 12d0
      logRho_max_for_any_custom_eos = logRho_max_for_all_custom_eos - 0.1d0

      ! logT blend
      if (logT_auto < logT_min_for_any_custom_eos) then
         blend_logT = 0d0
      else if (logT_auto < logT_min_for_all_custom_eos) then
         blend_logT = (logT_auto - logT_min_for_any_custom_eos) / (logT_min_for_all_custom_eos - logT_min_for_any_custom_eos)
      else if (logT_auto < logT_max_for_all_custom_eos) then
         blend_logT = 1d0
      else if (logT_auto < logT_max_for_any_custom_eos) then
         blend_logT = (logT_auto - logT_max_for_any_custom_eos) / (logT_max_for_all_custom_eos - logT_max_for_any_custom_eos)
      else
         blend_logT = 0
      end if

      ! logRho blend
      if (logRho_auto < logRho_min_for_any_custom_eos) then
         blend_logRho = 0d0
      else if (logRho_auto < logRho_min_for_all_custom_eos) then
         blend_logRho = (logRho_auto - logRho_min_for_any_custom_eos) / (logRho_min_for_all_custom_eos - logRho_min_for_any_custom_eos)
      else if (logRho_auto < logRho_max_for_all_custom_eos) then
         blend_logRho = 1d0
      else if (logRho_auto < logRho_max_for_any_custom_eos) then
         blend_logRho = (logRho_auto - logRho_max_for_any_custom_eos) / (logRho_max_for_all_custom_eos - logRho_max_for_any_custom_eos)
      else
         blend_logRho = 0
      end if

      ! logQ blend
      if (logQ_auto < logQ_min_for_any_custom_eos) then
         blend_logQ = 0d0
      else if (logQ_auto < logQ_min_for_all_custom_eos) then
         blend_logQ = (logQ_auto - logQ_min_for_any_custom_eos) / (logQ_min_for_all_custom_eos - logQ_min_for_any_custom_eos)
      else if (logQ_auto < logQ_max_for_all_custom_eos) then
         blend_logQ = 1d0
      else if (logQ_auto < logQ_max_for_any_custom_eos) then
         blend_logQ = (logQ_auto - logQ_max_for_any_custom_eos) / (logQ_max_for_all_custom_eos - logQ_max_for_any_custom_eos)
      else
         blend_logQ = 0
      end if

      ! combine blends
      blend = blend_logRho * blend_logT * blend_logQ

      frac = blend% val
      dfrac_dlogRho = blend% d1val1
      dfrac_dlogT = blend% d1val2

   end subroutine my_other_eos_frac


end module run_star_extras

