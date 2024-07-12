module custom_eos
   ! This module allows to load a custom equation of state (EOS) into MESA.

   use interp_2d_lib_db
   use star_lib
   use star_def
   use const_def
   use math_lib
   use eos_def

   use utils_lib, only: is_bad, mesa_error, mv, switch_str

   implicit none

   logical, parameter :: return_ierr_beyond_table_bounds = .true.
   ! from eos_helm_eval
   logical, parameter :: stop_for_is_bad = .false.

   ! the file EOS data
   integer, parameter :: jlogPgas = 1
   integer, parameter :: jlogE = 2
   integer, parameter :: jlogS = 3
   integer, parameter :: jchiRho = 4
   integer, parameter :: jchiT = 5
   integer, parameter :: jCp = 6
   integer, parameter :: jCv = 7
   integer, parameter :: jdE_dRho = 8
   integer, parameter :: jdS_dT = 9
   integer, parameter :: jdS_dRho = 10
   integer, parameter :: jmu = 11
   integer, parameter :: jlogfree_e = 12
   integer, parameter :: jgamma1 = 13
   integer, parameter :: jgamma3 = 14
   integer, parameter :: jgrad_ad = 15
   integer, parameter :: jeta = 16
   integer, parameter :: num_eos_file_vals = 16

   integer, parameter :: file_max_num_logQs = 1000

   ! values for a user-supplied EoS
   integer, parameter :: num_my_eos_Zs = 11
   integer, parameter :: num_my_eos_Xs = 11

   ! test if the default EoS works
   logical, dimension(num_my_eos_Xs, num_my_eos_Zs) :: my_eos_XZ_loaded

   ! the spacing of the qeos_h2o grid is the same; hence, we can use the qeos_sio2 parameters

   ! user-supplied EoS
   type (DT_XZ_Info), target :: my_eos_XZ_struct

   ! test if the default EoS works
   type (EosDT_XZ_Info), dimension(num_my_eos_Xs, num_my_eos_Zs), target :: &
      my_eos_XZ_data

   ! convinience id for a user-supplied EoS
   integer, parameter :: eosdt_my_eos = 1

   ! custom controls
   logical, parameter :: my_use_cache_for_eos = .true.



contains

   subroutine eos_init_custom_eos
      ! for looping
      integer :: i, j

      ! test if the default EoS works
      type (DT_XZ_Info), pointer :: my_eos_XZ_ptr

      real(dp) :: my_eos_spacing(1:11) = &
         (/ 0.00d0, 0.1d0, 0.2d0, 0.3d0, 0.4d0, 0.5d0, 0.6d0, &
         0.7d0, 0.8d0, 0.9d0, 1.0d0 /)

      ! tracking if a table was already loaded
      my_eos_XZ_loaded(:,:) = .false.

      ! user-supplied EoS

      my_eos_XZ_ptr => my_eos_XZ_struct

      ! #Z values in the grid
      my_eos_XZ_ptr % nZs = num_my_eos_Zs

      ! Z values of the EoS grid
      my_eos_XZ_ptr % Zs(1: my_eos_XZ_ptr % nZs) = my_eos_spacing

      my_eos_XZ_ptr % nXs_for_Z(1: my_eos_XZ_ptr % nZs) = &
         (/11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1/)

      do i = 1, my_eos_XZ_ptr % nZs
         ! X values of the EoS XZ grid
         j = my_eos_XZ_ptr % nXs_for_Z(i)
         my_eos_XZ_ptr % Xs_for_Z(1: j, i) = my_eos_spacing(1: j)
      end do

   end subroutine eos_init_custom_eos

   subroutine request_user_to_reinstall
      write(*,'(A)')
      write(*,'(A)')
      write(*,'(A)')
      write(*,'(A)')
      write(*,'(A)')
      write(*,'(A)')
      write(*,*) 'NOTICE: you need to install a new verion of the eos data.'
      write(*,*) 'Please update by removing the directory mesa/data/eosDT_data,'
      write(*,*) 'and rerunning the mesa ./install script.'
      write(*,'(A)')
      write(*,'(A)')
      call mesa_error(__FILE__,__LINE__)
   end subroutine request_user_to_reinstall


   subroutine check_for_error_in_eosDT_data(ierr, fname)
      integer, intent(in) :: ierr
      character (len=*) :: fname
      if (ierr == 0) return
      write(*,*) 'load eos tables ' // trim(fname)
      write(*,'(A)')
      write(*,'(A)')
      write(*,'(A)')
      write(*,'(A)')
      write(*,'(A)')
      write(*,'(a)') 'FATAL ERROR: missing or bad eos data.'
      write(*,'(a)') 'Please update by removing the directories ' &
         // 'mesa/data/eos*_data' &
         // ' and rerunning the mesa ./install script.'
      write(*,'(A)')
      call mesa_error(__FILE__,__LINE__)
   end subroutine check_for_error_in_eosDT_data


   subroutine load_single_eosDT_table_by_id( &
      which_eosdt, ep, ix, iz, ierr)
      use utils_lib
      integer, intent(in) :: which_eosdt
      type (EosDT_XZ_Info), pointer :: ep
      integer,intent(in) :: iz, ix
      integer, intent(out) :: ierr

      if (which_eosdt == eosdt_my_eos) then
         ep => my_eos_XZ_data(ix,iz)
         if (my_eos_XZ_loaded(ix,iz)) return
      else
         ierr = -1
         return
      end if

!$OMP CRITICAL(eosDT_load)
      if (which_eosdt == eosdt_my_eos) then
         if (.not. my_eos_XZ_loaded(ix,iz)) call do_read
      end if
!$OMP END CRITICAL(eosDT_load)

   contains

      subroutine do_read
         call read_one(ix,iz,ierr)
         if (ierr /= 0) return

         if (which_eosdt == eosdt_my_eos) then
            my_eos_XZ_loaded(ix,iz) = .true.
         else
            ierr = -1
            return
         end if
      end subroutine do_read

      subroutine read_one(ix,iz,ierr)
         integer, intent(in) :: ix, iz
         integer, intent(out) :: ierr
         character (len=256) :: fname, cache_filename, temp_cache_filename
         integer :: iounit1, iounit2
         real(dp) :: X, Z
         type (DT_xz_Info), pointer :: xz
         include 'formats'
         iounit1 = alloc_iounit(ierr); if (ierr /= 0) return
         iounit2 = alloc_iounit(ierr); if (ierr /= 0) return

         if (which_eosdt == eosdt_my_eos) then
            xz => my_eos_XZ_struct
         else
            write(*,*) 'unknown which_eosdt supplied'
            ierr = -1
            return
         end if

         call Get_eosDT_Table_Filenames(which_eosdt, xz, &
            ix, iz, fname, cache_filename, temp_cache_filename)
         call Load1_eosDT_Table(which_eosdt, ep, xz, &
            ix, iz, fname, cache_filename, temp_cache_filename, iounit1, iounit2, my_use_cache_for_eos, ierr)
         if (ierr /= 0) then
            write(*,*) 'Load1_eosDT_Table ierr', ierr, ix, iz, X, Z
         end if
         call free_iounit(iounit2)
         call free_iounit(iounit1)
      end subroutine read_one

   end subroutine load_single_eosDT_table_by_id

   subroutine Get_eosDT_Table_Filenames(which_eosdt, xz, &
      ix, iz, fname, cache_filename, temp_cache_filename)
      integer, intent(in) :: which_eosdt
      type (DT_xz_Info), pointer :: xz
      integer, intent(in) :: ix, iz

      character (len=10) :: data_dir = './src/data'
      character (len=*), intent(out) :: fname, cache_filename, temp_cache_filename
      character (len=256) :: Zstr, Xstr, suffix, data_dir_name, data_prefix
      real(dp) :: X, Z

      ! get the (Z, X) value of interest
      Z = xz% Zs(iz)
      X = xz% Xs_for_Z(ix,iz)

      ! create percentage format, i.e., 0.50 -> 50
      call setstr(Z,Zstr)
      call setstr(X,Xstr)

      if (which_eosdt == eosdt_my_eos) then
         ! specify the directory and prefix for the user-supplied EoS
         data_dir_name = '/my_eosDT/'
         data_prefix = 'my_eosDT_'
      else
         write(*,*) 'unknown eosdt'
         return
         call mesa_error(__FILE__,__LINE__)
      end if

      fname = trim(data_dir) //  &
         trim(data_dir_name) // trim(data_prefix) // &
         trim(Zstr) // 'z' // trim(Xstr) // 'x' // '.data'

      cache_filename = trim(eosDT_cache_dir) //  &
         '/' // trim(data_prefix) // &
         trim(Zstr) // 'z' // trim(Xstr) // 'x' // '.bin'

      temp_cache_filename = trim(eosDT_temp_cache_dir) //  &
         '/' // trim(data_prefix) // &
         trim(Zstr) // 'z' // trim(Xstr) // 'x' // '.bin'

   contains

      subroutine setstr(v,str)
         real(dp), intent(in) :: v
         character (len=*) :: str
         if (v > 0.99999d0) then
            str = '100'
         else if (v > 0.09999d0) then
            write(str, '(i2)') floor(100d0 * v + 0.5d0)
         else
            write(str, '(a,i1)') '0', floor(100d0 * v + 0.5d0)
         end if
      end subroutine setstr

   end subroutine Get_eosDT_Table_Filenames

   subroutine Load1_eosDT_Table(which_eosdt, ep, xz, &
      ix, iz, filename, cache_filename, temp_cache_filename, &
      io_unit, cache_io_unit, use_cache, info)
      integer, intent(in) :: which_eosdt
      type (EosDT_XZ_Info), pointer :: ep
      type (DT_xz_Info), pointer :: xz
      integer, intent(in) :: ix, iz
      character (*), intent(in) :: filename, cache_filename, temp_cache_filename
      integer, intent(in) :: io_unit, cache_io_unit
      logical, intent(in) :: use_cache
      integer, intent(out) :: info

      real(dp) :: X, Z, logQ, logT, X_in, Z_in
      integer :: j, i, k, iQ, ios, status
      character (len=1000) :: message
      real(dp), parameter :: tiny = 1d-10
      real(dp), pointer :: tbl(:,:,:,:) ! => ep% tbl1
      real(dp), pointer :: tbl2_1(:), tbl2(:,:,:)
      real(dp), target :: vec_ary(50)
      real(dp), pointer :: vec(:)
      integer :: n

      include 'formats'

      info = 0
      vec => vec_ary
      Z = xz% Zs(iz)
      X = xz% Xs_for_Z(ix,iz)

      write(message,*) 'open ', trim(filename)

      open(UNIT=io_unit, FILE=trim(filename), ACTION='READ', STATUS='OLD', IOSTAT=ios)
      call check_for_error_in_eosDT_data(ios, filename)

      read(io_unit,*,iostat=info)
      if (info /= 0) return

      read(io_unit,'(a)',iostat=info) message
      if (info == 0) call str_to_vector(message, vec, n, info)
      if (info /= 0 .or. n < 11) then
         write(*,'(a)') 'failed while reading ' // trim(filename)
         close(io_unit)
         info = -1
         return
      end if
      ep% version = int(vec(1))
      X_in = vec(2)
      Z_in = vec(3)
      ep% num_logTs = int(vec(4))
      ep% logT_min = vec(5)
      ep% logT_max = vec(6)
      ep% del_logT = vec(7)
      ep% num_logQs = int(vec(8))
      ep% logQ_min = vec(9)
      ep% logQ_max = vec(10)
      ep% del_logQ = vec(11)

      read(io_unit,*,iostat=info)
      if (info /= 0) return

      if (abs(X-X_in) > tiny .or. abs(Z-Z_in) > tiny) then
         write(*,*) 'bad header info in ' // trim(filename)
         info = -1
         close(io_unit)
         if (abs(X-X_in) > tiny) then
            write(*,'(a50,l1)') 'abs(X-X_in) > tiny', abs(X-X_in) > tiny
         end if
         if (abs(Z-Z_in) > tiny) then
            write(*,'(a50,l1)') 'abs(Z-Z_in) > tiny', abs(Z-Z_in) > tiny
         end if
         write(*,'(A)')
         call request_user_to_reinstall
         return
      end if

      if (show_allocations) write(*,2) 'Load1_eosDT_Table ep% tbl1', &
         sz_per_eos_point*nv*ep% num_logQs*ep% num_logTs + ep% num_logQs + ep% num_logTs
      allocate(ep% tbl1(sz_per_eos_point*nv*ep% num_logQs*ep% num_logTs), &
         ep% logQs(ep% num_logQs), ep% logTs(ep% num_logTs),   &
         STAT=info)
      if (info /= 0) then
         write(*,*)  "Info: ",info
         call mesa_error(__FILE__,__LINE__, "Allocation in Load1_eosDT_Table failed, you're likely out of memory")
      end if

      tbl(1:sz_per_eos_point,1:nv,1:ep% num_logQs,1:ep% num_logTs) =>  &
         ep% tbl1(1:sz_per_eos_point*nv*ep% num_logQs*ep% num_logTs)

      ep% logQs(1) = ep% logQ_min
      do i = 2, ep% num_logQs-1
         ep% logQs(i) = ep% logQs(i-1) + ep% del_logQ
      end do
      ep% logQs(ep% num_logQs) = ep% logQ_max

      ep% logTs(1) = ep% logT_min
      do i = 2, ep% num_logTs-1
         ep% logTs(i) = ep% logTs(i-1) + ep% del_logT
      end do
      ep% logTs(ep% num_logTs) = ep% logT_max

      if (use_cache) then
         call Read_EoS_Cache(X, Z, ep, cache_filename, cache_io_unit, ios)
         if (ios == 0) then
            close(io_unit)
            return
         end if
      end if

      status = 0
      allocate(tbl2_1(num_eos_file_vals*ep% num_logQs*ep% num_logTs), STAT=status)
      if (status .ne. 0) then
         info = -1
         return
      end if

      tbl2(1:num_eos_file_vals,1:ep% num_logQs,1:ep% num_logTs) =>  &
         tbl2_1(1:num_eos_file_vals*ep% num_logQs*ep% num_logTs)

      do iQ=1,ep% num_logQs

         read(io_unit,*,iostat=info)
         if (failed('skip line')) return

         read(io_unit,'(a)',iostat=info) message
         if (info == 0) call str_to_double(message, vec(1), info)
         if (failed('read logQ')) return
         logQ = vec(1)

         read(io_unit,*,iostat=info)
         if (failed('skip line')) return

         read(io_unit,*,iostat=info)
         if (failed('skip line')) return

         do i=1,ep% num_logTs

            read(io_unit,'(a)',iostat=info) message
            if (failed('read line')) then
               write(*,'(a)') trim(message)
               write(*,*) trim(filename)
               write(*,*) 'iQ, i', iQ, i
               write(*,*) 'logQ', logQ
               write(*,*) 'bad input line?'
               call mesa_error(__FILE__,__LINE__)
            end if

            call str_to_vector(message, vec, n, info)
            if (info /= 0 .or. n < 1+num_eos_file_vals) then
               write(*,'(a)') trim(message)
               write(*,*) trim(filename)
               write(*,*) 'iQ, i', iQ, i
               write(*,*) 'logQ', logQ
               write(*,*) 'bad input line?'
               call mesa_error(__FILE__,__LINE__)
            end if
            logT = vec(1)
            do j=1,num_eos_file_vals
               tbl2(j,iQ,i) = vec(1+j)
            end do

         enddo

         if(iQ == ep% num_logQs) exit
         read(io_unit,*,iostat=info)
         if (failed('skip line')) return
         read(io_unit,*,iostat=info)
         if (failed('skip line')) return

      end do

      close(io_unit)

      call Make_XEoS_Interpolation_Data(ep, tbl2_1, info)
      deallocate(tbl2_1)
      if (failed('Make_XEoS_Interpolation_Data')) return

      call Check_XEoS_Interpolation_Data(ep)

      if (.not. use_cache) return

      open(unit=cache_io_unit, &
         file=trim(switch_str(temp_cache_filename, cache_filename, use_mesa_temp_cache)), &
         iostat=ios,action='write', form='unformatted')

      if (ios == 0) then
         write(*,'(a)') 'write ' // trim(cache_filename)
         write(cache_io_unit)  &
            X_in, Z_in, ep% num_logTs, ep% logT_min, ep% logT_max, ep% del_logT,  &
            ep% num_logQs, ep% logQ_min, ep% logQ_max, ep% del_logQ, ep% version
         write(cache_io_unit)  &
            ep% tbl1(&
            1:sz_per_eos_point*nv*ep% num_logQs*ep% num_logTs)
         close(cache_io_unit)
         if(use_mesa_temp_cache) call mv(temp_cache_filename, cache_filename,.true.)
      end if


   contains

      subroutine Check_XEoS_Interpolation_Data(ep)
         use utils_lib,only:is_bad
         type (EosDT_XZ_Info), pointer :: ep
         ! for logT > 6.8 and logRho < -10, splines can get bogus higher order terms
         ! replace NaN's and Infinities with 0
         integer :: i, j, iQ, jtemp
         do i = 1, sz_per_eos_point
            do j = 1, nv
               do iQ = 1, ep% num_logQs
                  do jtemp = 1, ep% num_logTs
                     if (is_bad(tbl(i,j,iQ,jtemp))) then
                        tbl(i,j,iQ,jtemp) = 0
                     end if
                  end do
               end do
            end do
         end do
      end subroutine Check_XEoS_Interpolation_Data

      logical function failed(str)
         character (len=*), intent(in) :: str
         failed = (info /= 0)
         if (failed) write(*,*) 'Load1_eosDT_Table failed: ' // trim(str)
      end function failed


   end subroutine Load1_eosDT_Table


   subroutine Make_XEoS_Interpolation_Data(ep, tbl2_1, info)
      use interp_2d_lib_db
      use const_def, only: crad, ln10

      type (EosDT_XZ_Info), pointer :: ep
      real(dp), pointer :: tbl2_1(:) ! =(num_eos_file_vals, ep% num_logQs, ep% num_logTs)
      integer, intent(out) :: info

      real(dp) :: logQs(ep% num_logQs)              ! x vector, strict ascending
      real(dp) :: logTs(ep% num_logTs)                    ! y vector, strict ascending
      real(dp) :: Ts(ep% num_logTs)
      real(dp), allocatable, target :: f1_ary(:) ! data & spline coefficients
      real(dp), pointer :: f1(:), f(:,:,:), ep_tbl(:,:,:,:), tbl2(:,:,:)
      integer :: ibcxmin                   ! bc flag for x=xmin
      real(dp) :: bcxmin(ep% num_logTs)    ! bc data vs. y at x=xmin
      integer :: ibcxmax                   ! bc flag for x=xmax
      real(dp) :: bcxmax(ep% num_logTs)     ! bc data vs. y at x=xmax
      integer :: ibcymin                   ! bc flag for y=ymin
      real(dp) :: bcymin(ep% num_logQs)   ! bc data vs. x at y=ymin
      integer :: ibcymax                   ! bc flag for y=ymax
      real(dp) :: bcymax(ep% num_logQs)   ! bc data vs. x at y=ymax
      integer :: ili_logQs    ! =1: logRho grid is "nearly" equally spaced
      integer :: ili_logTs      ! =1: logT grid is "nearly" equally spaced
      integer :: ier            ! =0 on exit if there is no error.
      real(dp) :: logQ, Rho, logRho, T, P, Cv, chiRho, chiT, logT, logT0, logT1, logQ0, logQ1
      real(dp) :: gamma3, gamma1, grad_ad, Prad, E, S
      integer :: iQ, jtemp, ilogT, ilogQ
      real(dp) :: fval(num_eos_file_vals), df_dx(num_eos_file_vals), df_dy(num_eos_file_vals)

      real(dp) :: x, y, dlnT, energy, lnE, entropy, lnS, Pgas, lnPgas, dlogT, &
         dlnPgas_dlnd, dlnE_dlnd, dlnS_dlnd, dlnPgas_dlnT, dlnE_dlnT, dlnS_dlnT

      integer :: v, vlist(3), var, i, j, num_logQs, num_logTs, ii, jj
      character (len=256) :: message

      include 'formats'

      info = 0

      ! just use "not a knot" bc's at edges of tables
      ibcxmin = 0; bcxmin(:) = 0
      ibcxmax = 0; bcxmax(:) = 0
      ibcymin = 0; bcymin(:) = 0
      ibcymax = 0; bcymax(:) = 0

      num_logQs = ep% num_logQs
      num_logTs = ep% num_logTs

      ep_tbl(1:sz_per_eos_point,1:nv,1:num_logQs,1:num_logTs) =>  &
         ep% tbl1(1:sz_per_eos_point*nv*num_logQs*num_logTs)

      tbl2(1:num_eos_file_vals,1:num_logQs,1:num_logTs) =>  &
         tbl2_1(1:num_eos_file_vals*num_logQs*num_logTs)

      allocate(f1_ary(sz_per_eos_point * ep% num_logQs * ep% num_logTs))

      f1 => f1_ary
      f(1:sz_per_eos_point,1:num_logQs,1:num_logTs) => &
         f1_ary(1:sz_per_eos_point*num_logQs*num_logTs)

      do iQ = 1, ep% num_logQs
         logQs(iQ) = ep% logQ_min + (iQ-1) * ep% del_logQ
      end do

      do jtemp = 1, ep% num_logTs
         logTs(jtemp) = ep% logT_min + (jtemp-1) * ep% del_logT
      end do

      ! copy file eos variables to internal eos interpolation tables
      do j=1,num_logTs
         do i=1,num_logQs
            ep_tbl(1,i_lnPgas,i,j) = tbl2(jlogPgas,i,j)*ln10
            ep_tbl(1,i_lnE,i,j) = tbl2(jlogE,i,j)*ln10
            ep_tbl(1,i_lnS,i,j) = tbl2(jlogS,i,j)*ln10
            ep_tbl(1,i_grad_ad,i,j) = tbl2(jgrad_ad,i,j)
            ep_tbl(1,i_chiRho,i,j) = tbl2(jchiRho,i,j)
            ep_tbl(1,i_chiT,i,j) = tbl2(jchiT,i,j)
            ep_tbl(1,i_Cp,i,j) = tbl2(jCp,i,j)
            ep_tbl(1,i_Cv,i,j) = tbl2(jCv,i,j)
            ep_tbl(1,i_dE_dRho,i,j) = tbl2(jdE_dRho,i,j)
            ep_tbl(1,i_dS_dT,i,j) = tbl2(jdS_dT,i,j)
            ep_tbl(1,i_dS_dRho,i,j) = tbl2(jdS_dRho,i,j)
            ep_tbl(1,i_mu,i,j) = tbl2(jmu,i,j)
            ep_tbl(1,i_lnfree_e,i,j) = max(-4d0,tbl2(jlogfree_e,i,j))*ln10
            ! to protect against non-monotonic interpolation caused by extreme values
            ep_tbl(1,i_gamma1,i,j) = tbl2(jgamma1,i,j)
            ep_tbl(1,i_gamma3,i,j) = tbl2(jgamma3,i,j)
            ep_tbl(1,i_eta,i,j) = tbl2(jeta,i,j)
         end do
      end do

      ! create tables for bicubic spline interpolation
      do v = 1, nv
         do i=1,ep% num_logQs
            do j=1,ep% num_logTs
               f(1,i,j) = ep_tbl(1,v,i,j)
            end do
         end do
         call interp_mkbicub_db( &
            logQs,ep% num_logQs,logTs,ep% num_logTs,f1,ep% num_logQs, &
            ibcxmin,bcxmin,ibcxmax,bcxmax, &
            ibcymin,bcymin,ibcymax,bcymax, &
            ili_logQs,ili_logTs,ier)
         if (ier /= 0) then
            write(*,*) 'interp_mkbicub_db error happened for eos_value', v
            info = 3
            return
         end if
         do i=1,ep% num_logQs
            do j=1,ep% num_logTs
               ep_tbl(2,v,i,j) = f(2,i,j)
               ep_tbl(3,v,i,j) = f(3,i,j)
               ep_tbl(4,v,i,j) = f(4,i,j)
            end do
         end do
      end do


   end subroutine Make_XEoS_Interpolation_Data


   subroutine Read_EoS_Cache(X, Z, ep, cache_filename, io_unit, ios)
      real(dp), intent(in) :: X, Z
      type (EosDT_XZ_Info), pointer :: ep
      character (*), intent(in) :: cache_filename
      integer, intent(in) :: io_unit ! use this for file access
      integer, intent(out) :: ios

      real(dp) :: X_in, Z_in, logT_min_in, logT_max_in, del_logT_in,  &
         logQ_min_in, logQ_max_in, del_logQ_in
      integer :: num_logQs_in, num_logTs_in, version_in, i, j
      real(dp), parameter :: tiny = 1d-10

      include 'formats'

      ios = 0
      open(unit=io_unit,file=trim(cache_filename),action='read', &
         status='old',iostat=ios,form='unformatted')
      if (ios /= 0) return

      read(io_unit, iostat=ios)  &
         X_in, Z_in, num_logTs_in, logT_min_in, logT_max_in, del_logT_in,  &
         num_logQs_in, logQ_min_in, logQ_max_in, del_logQ_in, version_in
      if (ios /= 0) return

      if (ep% version /= version_in) then
         ios = 1
         write(*,*) 'read cache failed for version_in'
      end if
      if (ep% num_logQs /= num_logQs_in) then
         ios = 1
         write(*,*) 'read cache failed for ep% num_logQs'
      end if
      if (ep% num_logTs /= num_logTs_in) then
         ios = 1
         write(*,*) 'read cache failed for ep% num_logTs'
      end if
      if (abs(X-X_in) > tiny) then
         ios = 1
         write(*,*) 'read cache failed for X_in'
      end if
      if (abs(Z-Z_in) > tiny) then
         ios = 1
         write(*,*) 'read cache failed for Z_in'
      end if
      if (abs(ep% logT_min-logT_min_in) > tiny) then
         ios = 1
         write(*,*) 'read cache failed for eos_logT_min'
      end if
      if (abs(ep% logT_max-logT_max_in) > tiny) then
         ios = 1
         write(*,*) 'read cache failed for eos_logT_max'
      end if
      if (abs(ep% del_logT-del_logT_in) > tiny) then
         ios = 1
         write(*,*) 'read cache failed for eos_del_logT'
      end if
      if (abs(ep% logQ_min-logQ_min_in) > tiny) then
         ios = 1
         write(*,*) 'read cache failed for eos_logQ_min'
      end if
      if (abs(ep% logQ_max-logQ_max_in) > tiny) then
         ios = 1
         write(*,*) 'read cache failed for eos_logQ_max'
      end if
      if (abs(ep% del_logQ-del_logQ_in) > tiny) then
         ios = 1
         write(*,*) 'read cache failed for eos_del_logQ'
      end if

      if (ios /= 0) then
         close(io_unit); return
      end if

      read(io_unit, iostat=ios)  &
         ep% tbl1( &
         1:sz_per_eos_point*nv*ep% num_logQs*ep% num_logTs)
      if (ios /= 0) then
         close(io_unit); return
      end if

      close(io_unit)

   end subroutine Read_EoS_Cache

   subroutine Locate_logQ(ep, logQ, iQ, logQ0, logQ1, ierr)
      type (EosDT_xz_Info), pointer :: ep
      real(dp), intent(inout) :: logQ
      integer, intent(out) :: iQ
      real(dp), intent(out) :: logQ0, logQ1
      integer, intent(out) :: ierr
      ierr = 0
      iQ = int((logQ - ep% logQ_min)/ep% del_logQ + 1d-4) + 1
      if (iQ < 1 .or. iQ >= ep% num_logQs) then
         if (iQ < 1) then
            iQ = 1
            logQ0 = ep% logQ_min
            logQ1 = logQ0 + ep% del_logQ
            logQ = logQ0
            if (return_ierr_beyond_table_bounds) ierr = -1
         else
            iQ = ep% num_logQs-1
            logQ0 = ep% logQ_min + (iQ-1)*ep% del_logQ
            logQ1 = logQ0 + ep% del_logQ
            logQ = logQ1
            if (return_ierr_beyond_table_bounds) ierr = -1
         end if
      else
         logQ0 = ep% logQ_min + (iQ-1)*ep% del_logQ
         logQ1 = logQ0 + ep% del_logQ
      end if
      
      if (ierr /= 0) write(*,*) "beyond table bounds: ", logQ, logQ0, logQ1

   end subroutine Locate_logQ

   subroutine Locate_logT(ep, logT, iT, logT0, logT1, ierr)

      type (EosDT_xz_Info), pointer :: ep
      real(dp), intent(inout) :: logT
      integer, intent(out) :: iT
      real(dp), intent(out) :: logT0, logT1
      integer, intent(out) :: ierr
      ierr = 0
      iT = int((logT - ep% logT_min)/ep% del_logT + 1d-4) + 1
      if (iT < 1 .or. iT >= ep% num_logTs) then
         if (iT < 1) then
            iT = 1
            logT0 = ep% logT_min
            logT1 = logT0 + ep% del_logT
            logT = logT0
            if (return_ierr_beyond_table_bounds) ierr = -1
         else
            iT = ep% num_logTs-1
            logT0 = ep% logT_min + (iT-1)*ep% del_logT
            logT1 = logT0 + ep% del_logT
            logT = logT1
            if (return_ierr_beyond_table_bounds) ierr = -1
         end if
      else
         logT0 = ep% logT_min + (iT-1)*ep% del_logT
         logT1 = logT0 + ep% del_logT
      end if
   end subroutine Locate_logT

   subroutine Do_EoS_Interpolations( &
      nvlo, nvhi, n, nx, x, ny, y, fin1, i, j, &
      x0, xget, x1, y0, yget, y1, &
      fval, df_dx, df_dy, ierr)
      integer, intent(in) :: nvlo, nvhi, n, nx, ny
      real(dp), intent(in) :: x(:) ! (nx)
      real(dp), intent(in) :: y(:) ! (ny)
      real(dp), intent(in), pointer :: fin1(:) ! =(4,n,nx,ny)
      integer, intent(in) :: i, j           ! target cell in f
      real(dp), intent(in) :: x0, xget, x1      ! x0 <= xget <= x1;  x0 = xs(i), x1 = xs(i+1)
      real(dp), intent(in) :: y0, yget, y1      ! y0 <= yget <= y1;  y0 = ys(j), y1 = ys(j+1)
      real(dp), intent(inout), dimension(nv) :: fval, df_dx, df_dy
      integer, intent(out) :: ierr

      real(dp) :: xp, xpi, xp2, xpi2, ax, axbar, bx, bxbar, cx, cxi, hx2, cxd, cxdi, hx, hxi
      real(dp) :: yp, ypi, yp2, ypi2, ay, aybar, by, bybar, cy, cyi, hy2, cyd, cydi, hy, hyi
      real(dp) :: sixth_hx2, sixth_hy2, z36th_hx2_hy2
      real(dp) :: sixth_hx, sixth_hxi_hy2, z36th_hx_hy2
      real(dp) :: sixth_hx2_hyi, sixth_hy, z36th_hx2_hy
      integer :: k, ip1, jp1
      real(dp), pointer :: fin(:,:,:,:)

      include 'formats'

      ierr = 0

      fin(1:sz_per_eos_point,1:n,1:nx,1:ny) => &
         fin1(1:sz_per_eos_point*n*nx*ny)

      hx=x1-x0
      hxi=1d0/hx
      hx2=hx*hx

      xp=(xget-x0)*hxi

      xpi=1d0-xp
      xp2=xp*xp
      xpi2=xpi*xpi

      ax=xp2*(3d0-2d0*xp)
      axbar=1d0-ax

      bx=-xp2*xpi
      bxbar=xpi2*xp

      cx=xp*(xp2-1d0)
      cxi=xpi*(xpi2-1d0)
      cxd=3d0*xp2-1d0
      cxdi=-3d0*xpi2+1d0

      hy=y1-y0
      hyi=1d0/hy
      hy2=hy*hy

      yp=(yget-y0)*hyi

      ypi=1d0-yp
      yp2=yp*yp
      ypi2=ypi*ypi

      ay=yp2*(3d0-2d0*yp)
      aybar=1d0-ay

      by=-yp2*ypi
      bybar=ypi2*yp

      cy=yp*(yp2-1d0)
      cyi=ypi*(ypi2-1d0)
      cyd=3d0*yp2-1d0
      cydi=-3d0*ypi2+1d0

      sixth_hx2 = one_sixth*hx2
      sixth_hy2 = one_sixth*hy2
      z36th_hx2_hy2 = sixth_hx2*sixth_hy2

      sixth_hx = one_sixth*hx
      sixth_hxi_hy2 = sixth_hy2*hxi
      z36th_hx_hy2 = sixth_hx*sixth_hy2

      sixth_hx2_hyi = sixth_hx2*hyi
      sixth_hy = one_sixth*hy
      z36th_hx2_hy = sixth_hx2*sixth_hy

      ip1 = i+1
      jp1 = j+1

      !$omp simd
      do k = nvlo, nvhi
         ! bicubic spline interpolation

         ! f(1,i,j) = f(x(i),y(j))
         ! f(2,i,j) = d2f/dx2(x(i),y(j))
         ! f(3,i,j) = d2f/dy2(x(i),y(j))
         ! f(4,i,j) = d4f/dx2dy2(x(i),y(j))

         fval(k) = &
            xpi*( &
            ypi*fin(1,k,i,j)  +yp*fin(1,k,i,jp1)) &
            +xp*(ypi*fin(1,k,ip1,j)+yp*fin(1,k,ip1,jp1)) &
            +sixth_hx2*( &
            cxi*(ypi*fin(2,k,i,j) +yp*fin(2,k,i,jp1))+ &
            cx*(ypi*fin(2,k,ip1,j)+yp*fin(2,k,ip1,jp1))) &
            +sixth_hy2*( &
            xpi*(cyi*fin(3,k,i,j) +cy*fin(3,k,i,jp1))+ &
            xp*(cyi*fin(3,k,ip1,j)+cy*fin(3,k,ip1,jp1))) &
            +z36th_hx2_hy2*( &
            cxi*(cyi*fin(4,k,i,j) +cy*fin(4,k,i,jp1))+ &
            cx*(cyi*fin(4,k,ip1,j)+cy*fin(4,k,ip1,jp1)))

         ! derivatives of bicubic splines
         df_dx(k) = &
            hxi*( &
            -(ypi*fin(1,k,i,j)  +yp*fin(1,k,i,jp1)) &
            +(ypi*fin(1,k,ip1,j)+yp*fin(1,k,ip1,jp1))) &
            +sixth_hx*( &
            cxdi*(ypi*fin(2,k,i,j) +yp*fin(2,k,i,jp1))+ &
            cxd*(ypi*fin(2,k,ip1,j)+yp*fin(2,k,ip1,jp1))) &
            +sixth_hxi_hy2*( &
            -(cyi*fin(3,k,i,j)  +cy*fin(3,k,i,jp1)) &
            +(cyi*fin(3,k,ip1,j)+cy*fin(3,k,ip1,jp1))) &
            +z36th_hx_hy2*( &
            cxdi*(cyi*fin(4,k,i,j) +cy*fin(4,k,i,jp1))+ &
            cxd*(cyi*fin(4,k,ip1,j)+cy*fin(4,k,ip1,jp1)))

         df_dy(k) = &
            hyi*( &
            xpi*(-fin(1,k,i,j) +fin(1,k,i,jp1))+ &
            xp*(-fin(1,k,ip1,j)+fin(1,k,ip1,jp1))) &
            +sixth_hx2_hyi*( &
            cxi*(-fin(2,k,i,j) +fin(2,k,i,jp1))+ &
            cx*(-fin(2,k,ip1,j)+fin(2,k,ip1,jp1))) &
            +sixth_hy*( &
            xpi*(cydi*fin(3,k,i,j) +cyd*fin(3,k,i,jp1))+ &
            xp*(cydi*fin(3,k,ip1,j)+cyd*fin(3,k,ip1,jp1))) &
            +z36th_hx2_hy*( &
            cxi*(cydi*fin(4,k,i,j) +cyd*fin(4,k,i,jp1))+ &
            cx*(cydi*fin(4,k,ip1,j)+cyd*fin(4,k,ip1,jp1)))

      end do

   end subroutine Do_EoS_Interpolations

   subroutine Get1_eosdt_XTable_Results( &
      which_eosdt, ix, iz, logRho_in, logT_in, &
      res, d_dlnd, d_dlnT, ierr)

      integer, intent(in) :: which_eosdt
      integer, intent(in) :: ix, iz
      real(dp), intent(in) :: logRho_in, logT_in
      real(dp), intent(inout), dimension(nv) :: res, d_dlnd, d_dlnT
      integer, intent(out) :: ierr

      real(dp), parameter :: ln10sq = ln10*ln10
      real(dp) :: &
         fval(nv), df_dx(nv), df_dy(nv), &
         df_dlnd(nv), df_dlnT(nv), &
         energy, entropy, P, Pgas, Prad, x, y, &
         dx_dlnd, dx_dlnT, dy_dlnd, dy_dlnT, &
         chiT, chiRho, Cv, gamma1, numeric, &
         dS_dlnd, dS_dlnT, dE_dlnd, dE_dlnT, &
         dPgas_dlnd, dPgas_dlnT, dP_dlnd, dP_dlnT
      real(dp) :: logQ0, logQ1, logT0, logT1, logRho0, logRho1
      integer :: iQ, jtemp, k, j, irho
      type (EosDT_xz_Info), pointer :: ep
      logical, parameter :: show = .false.
      real(dp) :: logRho, logT, logQ

      include 'formats'

      logRho = logRho_in
      logT = logT_in
      logQ = logRho - 2*logT + 12

      ierr = 0
      call load_single_eosDT_table_by_id(which_eosdt, ep, ix, iz, ierr)
      if (ierr /= 0) return

      call Locate_logQ(ep, logQ, iQ, logQ0, logQ1, ierr)
      if (ierr /= 0) then
         write(*,1) 'eosDT failed in Locate_logQ', logQ
         return
      end if

      call Locate_logT(ep, logT, jtemp, logT0, logT1, ierr)
      if (ierr /= 0) then
         write(*,1) 'eosDT failed in Locate_logT', logT
         return
      end if

      call Do_EoS_Interpolations( &
         1, nv, nv, ep% num_logQs, ep% logQs, ep% num_logTs, ep% logTs, &
         ep% tbl1, iQ, jtemp, logQ0, logQ, logQ1, logT0, logT, logT1, &
         fval, df_dx, df_dy, ierr)
      if (ierr /= 0) then
         write(*,1) 'failed in Do_EoS_Interpolations'
         return
      end if

      if (is_bad(fval(i_lnS))) then
         ierr = -1
         if (.not. stop_for_is_bad) return
         write(*,1) 'fval(i_lnS), logRho, logT', fval(i_lnS), logRho, logT
         call mesa_error(__FILE__,__LINE__,'after Do_Interp_with_2nd_derivs')
      end if

      res(i_lnPgas) = fval(i_lnPgas)
      res(i_lnE) = fval(i_lnE)
      res(i_lnS) = fval(i_lnS)

      if (is_bad(res(i_lnS))) then
         ierr = -1
         if (.not. stop_for_is_bad) return
         write(*,1) 'res(i_lnS), logRho, logT', res(i_lnS), logRho, logT
         call mesa_error(__FILE__,__LINE__,'after interpolation')
      end if

      if (is_bad(res(i_lnS)) .or. res(i_lnS) > ln10*100) then
         ierr = -1
         if (.not. stop_for_is_bad) return
         write(*,1) 'res(i_lnS), logRho, logT', res(i_lnS), logRho, logT
         call mesa_error(__FILE__,__LINE__,'after interpolation')
      end if

      res(i_grad_ad) = fval(i_grad_ad)
      
      res(i_chiRho) = fval(i_chiRho)
      res(i_chiT) = fval(i_chiT)

      res(i_Cp) = fval(i_Cp)
      res(i_Cv) = fval(i_Cv)

      res(i_dE_dRho) = fval(i_dE_dRho)
      res(i_dS_dT) = fval(i_dS_dT)
      res(i_dS_dRho) = fval(i_dS_dRho)

      res(i_mu) = fval(i_mu)
      res(i_lnfree_e) = fval(i_lnfree_e)
      res(i_gamma1) = fval(i_gamma1)
      res(i_gamma3) = fval(i_gamma3)
      res(i_eta) = fval(i_eta)

      ! convert df_dx and df_dy to df_dlogRho_c_T and df_dlogT_c_Rho

      ! df_dx is df_dlogQ at const T
      ! df_dy is df_dlogT_c_Rho at const Q
      ! logQ = logRho - 2*logT + 12

      ! f = f(logQ(logRho,logT),logT)
      ! df/dlogRho|T = df/dlogQ|T * dlogQ/dlogRho|T = df_dx
      ! df/dlogT|Rho = df/dlogT|Q + df/dlogQ|T * dlogQ/dlogT|Rho = df_dy - 2*df_dx

      do k=1,nv
         df_dlnd(k) = df_dx(k)/ln10
         df_dlnT(k) = df_dy(k)/ln10 - 2d0*df_dlnd(k)
      end do

      d_dlnd(i_lnPgas) = df_dlnd(i_lnPgas)
      d_dlnd(i_lnE) = df_dlnd(i_lnE)
      d_dlnd(i_lnS) = df_dlnd(i_lnS)
      d_dlnd(i_grad_ad) = df_dlnd(i_grad_ad)
      d_dlnd(i_chiRho) = df_dlnd(i_chiRho)
      d_dlnd(i_chiT) = df_dlnd(i_chiT)

      d_dlnd(i_Cp) = df_dlnd(i_Cp)
      d_dlnd(i_Cv) = df_dlnd(i_Cv)
      d_dlnd(i_dE_dRho) = df_dlnd(i_dE_dRho)
      d_dlnd(i_dS_dT) = df_dlnd(i_dS_dT)
      d_dlnd(i_dS_dRho) = df_dlnd(i_dS_dRho)
      d_dlnd(i_mu) = df_dlnd(i_mu)
      d_dlnd(i_lnfree_e) = df_dlnd(i_lnfree_e)
      d_dlnd(i_gamma1) = df_dlnd(i_gamma1)
      d_dlnd(i_gamma3) = df_dlnd(i_gamma3)
      d_dlnd(i_eta) = df_dlnd(i_eta)

      d_dlnT(i_lnPgas) = df_dlnT(i_lnPgas)
      d_dlnT(i_lnE) = df_dlnT(i_lnE)
      d_dlnT(i_lnS) = df_dlnT(i_lnS)
      d_dlnT(i_grad_ad) = df_dlnT(i_grad_ad)
      d_dlnT(i_chiRho) = df_dlnT(i_chiRho)
      d_dlnT(i_chiT) = df_dlnT(i_chiT)
      d_dlnT(i_Cp) = df_dlnT(i_Cp)
      d_dlnT(i_Cv) = df_dlnT(i_Cv)
      d_dlnT(i_dE_dRho) = df_dlnT(i_dE_dRho)
      d_dlnT(i_dS_dT) = df_dlnT(i_dS_dT)
      d_dlnT(i_dS_dRho) = df_dlnT(i_dS_dRho)
      d_dlnT(i_mu) = df_dlnT(i_mu)
      d_dlnT(i_lnfree_e) = df_dlnT(i_lnfree_e)
      d_dlnT(i_gamma1) = df_dlnT(i_gamma1)
      d_dlnT(i_gamma3) = df_dlnT(i_gamma3)
      d_dlnT(i_eta) = df_dlnT(i_eta)

      if (is_bad(d_dlnd(i_lnS)) .or. is_bad(d_dlnT(i_lnS))) then
         ierr = -1
         if (.not. stop_for_is_bad) return
         write(*,1) 'fval(i_lnS)', fval(i_lnS)
         write(*,1) 'd_dlnd(i_lnS)', d_dlnd(i_lnS)
         write(*,1) 'd_dlnT(i_lnS)', d_dlnT(i_lnS)
         call mesa_error(__FILE__,__LINE__,'Get1_eosdt_XTable_Results')
      end if

   end subroutine Get1_eosdt_XTable_Results

   subroutine Get1_eosdt_for_X( &
      rq, which_eosdt, xz, iz, X, logRho, logT, &
      res, dlnd, dlnT, d_dX, ierr)
      type (EoS_General_Info), pointer :: rq
      integer, intent(in) :: which_eosdt
      type (DT_xz_Info), pointer :: xz
      integer, intent(in) :: iz ! the index in eos_Zs
      real(dp), intent(in) :: X, logRho, logT
      real(dp), intent(inout), dimension(nv) :: &
         res, dlnd, dlnT, d_dX
      integer, intent(out) :: ierr

      real(dp), dimension(nv, 4) :: &
         res_zx, dlnd_zx, dlnT_zx
      real(dp) :: dX, dX1, dX2, dX3, c(4), dcdX(4), denom, delX, coef, dcoef_dX, alfa, beta, dalfa_dX, dbeta_dX, tiny
      character (len=256) :: message
      integer :: ix, ix_lo, ix_hi, j, num_Xs
      logical :: what_we_use_is_equal_spaced

      include 'formats'

      ierr = 0
      tiny = rq% tiny_fuzz

      num_Xs = xz% nXs_for_Z(iz)

      if (xz% Xs_for_Z(1,iz) /= 0d0) then
         write(*, *) 'error: Get1_eosdt_for_X assumes xz% nXs_for_Z(1) == 0'
         call mesa_error(__FILE__,__LINE__)
      end if

      if (X < tiny .or. num_Xs == 1) then
         call Get1_eosdt_XTable_Results( &
            which_eosdt, 1, iz, logRho, logT, &
            res, dlnd, dlnT, ierr)
         d_dX = 0
         return
      end if

      if (X >= xz% Xs_for_Z(num_Xs,iz)) then

         call Get1_eosdt_XTable_Results( &
            which_eosdt, num_Xs, iz, logRho, logT, &
            res, dlnd, dlnT, ierr)
         d_dX = 0

         if (is_bad(res(i_lnS))) then
            ierr = -1
            if (.not. stop_for_is_bad) return
            write(*,1) 'res(i_lnS), logRho, logT', res(i_lnS), logRho, logT
            call mesa_error(__FILE__,__LINE__,'Get1_eosdt_for_X num_Xs')
         end if

         return
      end if

      if (rq% eosDT_use_linear_interp_for_X .or. num_Xs == 2) then
         call do_linear
         return
      end if

      ix_hi = -1
      if (X <= xz% Xs_for_Z(2,iz)) then
         ix_lo = 1; ix_hi = 3
      else if (X >= xz% Xs_for_Z(num_Xs-1,iz)) then
         ix_lo = num_Xs-2; ix_hi = num_Xs
      else
         do ix = 3, num_Xs-1
            if (X <= xz% Xs_for_Z(ix,iz)) then
               ix_lo = ix-2; ix_hi = ix+1; exit
            end if
         end do
      end if

      if (ix_hi < 0) then
         write(*, *) 'X', X
         write(*, *) 'ix_lo', ix_lo
         write(*, *) 'ix_hi', ix_hi
         write(*, *) 'error: Get1_eosdt_for_X logic bug'
         call mesa_error(__FILE__,__LINE__)
      end if

      what_we_use_is_equal_spaced = .true.
      dX1 = xz% Xs_for_Z(ix_lo+1,iz)-xz% Xs_for_Z(ix_lo,iz)
      dX2 = xz% Xs_for_Z(ix_lo+2,iz)-xz% Xs_for_Z(ix_lo+1,iz)
      if (ix_hi-ix_lo==2) then ! check that the 3 table X's are equal spaced
         if (abs(dX1 - dX2) > tiny) what_we_use_is_equal_spaced = .false.
      else ! check that the 4 table X's are equal spaced
         dX3 = xz% Xs_for_Z(ix_hi,iz)-xz% Xs_for_Z(ix_lo+2,iz)
         if (abs(dX1 - dX2) > tiny .or. abs(dX2 - dX3) > tiny) &
            what_we_use_is_equal_spaced = .false.
      end if

      if (.not. what_we_use_is_equal_spaced) then
         call do_linear
         if (is_bad(d_dX(1))) then
            call mesa_error(__FILE__,__LINE__,'Get1_eosdt_for_X bad d_dX; linear')
         end if
         return
      end if

      do ix=ix_lo, ix_hi
         j = ix-ix_lo+1
         call Get1_eosdt_XTable_Results( &
            which_eosdt, ix, iz, logRho, logT, &
            res_zx(:, j), dlnd_zx(:, j), dlnT_zx(:, j), &
            ierr)
         if (ierr /= 0) return
      end do

! zero these for now
      res_zx(i_phase:i_latent_ddlnRho,:) = 0d0
      res_zx(i_frac:i_frac+num_eos_frac_results-1,:) = 0d0

      dlnd_zx(i_phase:i_latent_ddlnRho,:) = 0d0
      dlnd_zx(i_frac:i_frac+num_eos_frac_results-1,:) = 0d0

      dlnT_zx(i_phase:i_latent_ddlnRho,:) = 0d0
      dlnT_zx(i_frac:i_frac+num_eos_frac_results-1,:) = 0d0


      delX = X - xz% Xs_for_Z(ix_lo,iz)
      dX = dX1

      if (ix_hi-ix_lo==2) then

         denom = 2*dX*dX
         c(1) = (2*dX*dX - 3*dX*delX + delX*delX)/denom
         c(2) = 2*(2*dX-delX)*delX/denom
         c(3) = delX*(delX-dX)/denom
         res(:) = c(1)*res_zx(:, 1) + c(2)*res_zx(:, 2) + c(3)*res_zx(:, 3)

         dlnd(:) = &
            c(1)*dlnd_zx(:,1) + &
            c(2)*dlnd_zx(:,2) + &
            c(3)*dlnd_zx(:,3)
         dlnT(:) = &
            c(1)*dlnT_zx(:,1) + &
            c(2)*dlnT_zx(:,2) + &
            c(3)*dlnT_zx(:,3)

         dcdx(1) = (-3*dX + 2*delX)/denom
         dcdx(2) = 2*(2*dX-2*delX)/denom
         dcdx(3) = (2*delX-dX)/denom

         d_dX(:) = &
            dcdX(1)*res_zx(:,1) + &
            dcdX(2)*res_zx(:,2) + &
            dcdX(3)*res_zx(:,3)

         if (is_bad(d_dX(1))) then
            call mesa_error(__FILE__,__LINE__,'Get1_eosdt_for_X bad d_dX; 3')
         end if

      else

         coef = (X-xz% Xs_for_Z(ix_lo+1,iz))/dX
         ! coef = fractional location of X between 2nd and 3rd X's for fit.
         ! coef is weight for the quadratic based on points 2, 3, 4 of fit.
         ! (1-coef) is weight for quadratic based on points 1, 2, 3 of fit.
         coef = min(1d0,max(0d0,coef))
         c(1) = -coef*(coef-1)*(coef-1)/2
         c(2) = (2 - coef*coef*(5 - 3*coef))/2
         c(3) = coef*(1 + coef*(4 - 3*coef))/2
         c(4) = coef*coef*(coef-1)/2
         res(:) = c(1)*res_zx(:, 1) + &
            (c(2)*res_zx(:, 2) + &
            (c(3)*res_zx(:, 3) + &
            c(4)*res_zx(:, 4)))

         dlnd(:) = &
            c(1)*dlnd_zx(:, 1) + &
            (c(2)*dlnd_zx(:, 2) + &
            (c(3)*dlnd_zx(:, 3) + &
            c(4)*dlnd_zx(:, 4)))
         dlnT(:) = &
            c(1)*dlnT_zx(:, 1) + &
            (c(2)*dlnT_zx(:, 2) + &
            (c(3)*dlnT_zx(:, 3) + &
            c(4)*dlnT_zx(:, 4)))

         dcoef_dX = 1d0/dX
         dcdX = 0
         dcdX(1) = -(3*coef*coef-4*coef+1)/2*dcoef_dX
         dcdX(2) = (9*coef*coef-10*coef)/2*dcoef_dX
         dcdX(3) = -(9*coef*coef-8*coef-1)/2*dcoef_dX
         dcdX(4) = coef*(3*coef-2)/2*dcoef_dX

         d_dX(:) = &
            dcdX(1)*res_zx(:,1) + &
            dcdX(2)*res_zx(:,2) + &
            dcdX(3)*res_zx(:,3) + &
            dcdX(4)*res_zx(:,4)

         if (is_bad(d_dX(1))) then
            call mesa_error(__FILE__,__LINE__,'Get1_eosdt_for_X bad d_dX; 4')
         end if

      end if

   contains

      subroutine do_linear

         do ix = 2, num_Xs
            if (xz% Xs_for_Z(ix,iz) >= X) exit
         end do

         j = 1
         call Get1_eosdt_XTable_Results( &
            which_eosdt, ix-1, iz, logRho, logT, &
            res_zx(:,j), dlnd_zx(:,j), dlnT_zx(:,j), &
            ierr)
         if (ierr /= 0) then
            if (.not. stop_for_is_bad) return
            call mesa_error(__FILE__,__LINE__,'Get1_eosdt_for_X')
         end if

         j = 2
         call Get1_eosdt_XTable_Results( &
            which_eosdt, ix, iz, logRho, logT, &
            res_zx(:,j), dlnd_zx(:,j), dlnT_zx(:,j), &
            ierr)
         if (ierr /= 0) then
            if (.not. stop_for_is_bad) return
            call mesa_error(__FILE__,__LINE__,'Get1_eosdt_for_X')
         end if

         ! zero these for now
         res_zx(i_phase:i_latent_ddlnRho,:) = 0d0
         res_zx(i_frac:i_frac+num_eos_frac_results-1,:) = 0d0

         dlnd_zx(i_phase:i_latent_ddlnRho,:) = 0d0
         dlnd_zx(i_frac:i_frac+num_eos_frac_results-1,:) = 0d0

         dlnT_zx(i_phase:i_latent_ddlnRho,:) = 0d0
         dlnT_zx(i_frac:i_frac+num_eos_frac_results-1,:) = 0d0


         alfa = (X - xz% Xs_for_Z(ix,iz))/(xz% Xs_for_Z(ix-1,iz) - xz% Xs_for_Z(ix,iz))
         beta = 1d0 - alfa

         dalfa_dX = 1d0 / (xz% Xs_for_Z(ix-1,iz) - xz% Xs_for_Z(ix,iz))
         dbeta_dX = -dalfa_dX

         do j=1,nv

            res(j) = alfa*res_zx(j,1) + beta*res_zx(j,2)

            dlnd(j) = &
               alfa*dlnd_zx(j,1) + beta*dlnd_zx(j,2)

            dlnT(j) = &
               alfa*dlnT_zx(j,1) + beta*dlnT_zx(j,2)

            d_dX(j) = &
               dalfa_dX*res_zx(j,1) + dbeta_dX*res_zx(j,2)

         end do

      end subroutine do_linear

   end subroutine Get1_eosdt_for_X

   subroutine Get1_eosdt_Results( & ! blend in Z
      rq, which_eosdt, xz, Z, X, logRho, logT, &
      res, dlnd, dlnT, dX, dZ, ierr)
      use chem_def
      type (EoS_General_Info), pointer :: rq
      integer, intent(in) :: which_eosdt
      type (DT_xz_Info), pointer :: xz
      real(dp), intent(in) :: Z, X, logRho, logT
      real(dp), intent(inout), dimension(nv) :: res, dlnd, dlnT, dX, dZ
      integer, intent(out) :: ierr

      real(dp), dimension(nv, 2) :: res_zx, dlnd_zx, dlnT_zx, dX_zx
      real(dp) :: denom, c(2), dcdZ(2), tiny

      integer :: iz, j, ci

      ! alternative interpolation
      integer :: iz_lo, iz_hi
      real(dp) :: dZ_alt, dZ1, dZ2, dZ3, coef, dcoef_dZ, delZ
      real(dp) :: c_alt(4), dcdZ_alt(4)
      logical :: what_we_use_is_equal_spaced
      real(dp), dimension(nv, 4) :: &
         res_zx_alt, dlnd_zx_alt, dlnT_zx_alt, dX_zx_alt

      include 'formats'

      ierr = 0
      tiny = rq% tiny_fuzz

      if (xz% nZs < 3) then
         write(*, *) 'error: Get1_eosdt_Results assumes nZs >= 3'
         call mesa_error(__FILE__,__LINE__)
      end if

      if (xz% Zs(1) /= 0) then
         write(*, *) 'error: Get1_eosdt_Results assumes eos_Zs(1) == 0'
         call mesa_error(__FILE__,__LINE__)
      end if

      if (abs(xz% Zs(1) - 2*xz% Zs(2) + xz% Zs(3)) > tiny) then
         write(*, *) 'error: Get1_eosdt_Results assumes equal spaced Zs(1:3)'
         call mesa_error(__FILE__,__LINE__)
      end if

      if (Z <= max(1d-20,xz% Zs(1))) then
         call Get1_eosdt_for_X( &
            rq, which_eosdt, xz, 1, X, &
            logRho, logT, &
            res, dlnd, dlnT, dX, ierr)
         dZ = 0
         return
      end if

! zero these for now
      res_zx(i_phase:i_latent_ddlnRho,:) = 0d0
      res_zx(i_frac:i_frac+num_eos_frac_results-1,:) = 0d0

      dlnd_zx(i_phase:i_latent_ddlnRho,:) = 0d0
      dlnd_zx(i_frac:i_frac+num_eos_frac_results-1,:) = 0d0

      dlnT_zx(i_phase:i_latent_ddlnRho,:) = 0d0
      dlnT_zx(i_frac:i_frac+num_eos_frac_results-1,:) = 0d0

      if (Z >= xz% Zs(xz% nZs)) then
         call Get1_eosdt_for_X( &
            rq, which_eosdt, xz, xz% nZs, X, &
            logRho, logT, &
            res, dlnd, dlnT, dX, ierr)
         dZ = 0
         return
      end if

      ! method 1: linear interpolation in Z
      if (rq % eos_integer_ctrl(2) == 1) then
         do iz = 2, xz% nZs
            if (Z < xz% Zs(iz)) then
               call do_interp2(iz-1,iz,ierr)
               if (ierr /= 0) return
               exit
            end if
         end do

         do j=1,nv

            res(j) = c(1)*res_zx(j,1) + c(2)*res_zx(j,2)

            dlnd(j) = &
               c(1)*dlnd_zx(j,1) + c(2)*dlnd_zx(j,2)

            dlnT(j) = &
               c(1)*dlnT_zx(j,1) + c(2)*dlnT_zx(j,2)

            dX(j) = &
               c(1)*dX_zx(j,1) + c(2)*dX_zx(j,2)

            dZ(j) = &
               dcdZ(1)*res_zx(j,1) + dcdZ(2)*res_zx(j,2)

         end do

      ! cubic (?) interpolation
      ! essentially the way MESA is doing it for X.
      else if (rq % eos_integer_ctrl(2) == 2) then
         
         iz_hi = -1
         if (Z <= xz % Zs(2)) then
            iz_lo = 1; iz_hi = 3
         else if (Z >= xz % Zs(xz % nZs-1)) then
            iz_lo = xz % nZs-2; iz_hi = xz % nZs
         else
            do iz = 3, xz % nZs-1
               if (Z <= xz % Zs(iz)) then
                  iz_lo = iz-2; iz_hi = iz+1; exit
               end if
            end do
         end if

         if (iz_hi < 0) then
            write(*,*) 'Z', Z
            write(*,*) 'iz_lo', iz_lo
            write(*,*) 'iz_hi', iz_hi
            write(*,*) 'error: Get1_eosdt_Results logic bug'
            call mesa_error(__FILE__,__LINE__)
         end if

         what_we_use_is_equal_spaced = .true.
         dZ1 = xz % Zs(iz_lo+1)-xz % Zs(iz_lo)
         dZ2 = xz % Zs(iz_lo+2)-xz % Zs(iz_lo+1)
         if (iz_hi-iz_lo==2) then ! check that the 3 table Z's are equal spaced
            if (abs(dZ1 - dZ2) > tiny) what_we_use_is_equal_spaced = .false.
         else ! check that the 4 table Z's are equal spaced
            dZ3 = xz % Zs(iz_hi)-xz % Zs(iz_lo+2)
            if (abs(dZ1 - dZ2) > tiny .or. abs(dZ2 - dZ3) > tiny) &
               what_we_use_is_equal_spaced = .false.
         end if

         if (.not. what_we_use_is_equal_spaced) then
            write(*,*) 'error: Get1_eosdt_Results EoS in Z not equally spaced.'
            call mesa_error(__FILE__,__LINE__)
         end if

         do iz=iz_lo, iz_hi
            j = iz-iz_lo+1
            call Get1_eosdt_for_X( &
               rq, which_eosdt, xz, iz, X, &
               logRho, logT, &
               res_zx_alt(:,j), dlnd_zx_alt(:,j), dlnT_zx_alt(:,j), dX_zx_alt(:,j), &
               ierr)
            if (ierr /= 0) return
         end do

         delZ = Z - xz % Zs(iz_lo)
         dZ_alt = dZ1

         if (iz_hi - iz_lo == 2) then
            denom = 2*dZ_alt*dZ_alt
            c_alt(1) = (2*dZ_alt*dZ_alt - 3*dZ_alt*delZ + delZ*delZ)/denom
            c_alt(2) = 2*(2*dZ_alt-delZ)*delZ/denom
            c_alt(3) = delZ*(delZ-dZ_alt)/denom
            res(:) = c_alt(1)*res_zx_alt(:, 1) + c_alt(2)*res_zx_alt(:, 2) + c_alt(3)*res_zx_alt(:, 3)

            dlnd(:) = &
               c_alt(1)*dlnd_zx_alt(:,1) + &
               c_alt(2)*dlnd_zx_alt(:,2) + &
               c_alt(3)*dlnd_zx_alt(:,3)
            
            dlnT(:) = &
               c_alt(1)*dlnT_zx_alt(:,1) + &
               c_alt(2)*dlnT_zx_alt(:,2) + &
               c_alt(3)*dlnT_zx_alt(:,3)
            
            dX(:) = &
               c_alt(1)*dX_zx_alt(:,1) + &
               c_alt(2)*dX_zx_alt(:,2) + &
               c_alt(3)*dX_zx_alt(:,3)
            
            dcdZ_alt(1) = (-3*dZ_alt + 2*delZ)/denom
            dcdZ_alt(2) = 2*(2*dZ_alt-2*delZ)/denom
            dcdZ_alt(3) = (2*delZ-dZ_alt)/denom

            dZ(:) = &
               dcdZ_alt(1)*res_zx_alt(:,1) + &
               dcdZ_alt(2)*res_zx_alt(:,2) + &
               dcdZ_alt(3)*res_zx_alt(:,3)

            if (is_bad(dZ(1))) then
               call mesa_error(__FILE__,__LINE__,'Get1_eosdt_Results bad dZ; 3')
            end if

         else

            coef = (Z - xz % Zs(iz_lo+1))/dZ_alt
            coef = min(1d0,max(0d0,coef))
            c_alt(1) = -coef*(coef-1)*(coef-1)/2
            c_alt(2) = (2 - coef*coef*(5 - 3*coef))/2
            c_alt(3) = coef*(1 + coef*(4 - 3*coef))/2
            c_alt(4) = coef*coef*(coef-1)/2
            res(:) = c_alt(1)*res_zx_alt(:, 1) + &
               (c_alt(2)*res_zx_alt(:, 2) + &
               (c_alt(3)*res_zx_alt(:, 3) + &
               c_alt(4)*res_zx_alt(:, 4)))

            dlnd(:) = &
               c_alt(1)*dlnd_zx_alt(:, 1) + &
               (c_alt(2)*dlnd_zx_alt(:, 2) + &
               (c_alt(3)*dlnd_zx_alt(:, 3) + &
               c_alt(4)*dlnd_zx_alt(:, 4)))
            
            dlnT(:) = &
               c_alt(1)*dlnT_zx_alt(:, 1) + &
               (c_alt(2)*dlnT_zx_alt(:, 2) + &
               (c_alt(3)*dlnT_zx_alt(:, 3) + &
               c_alt(4)*dlnT_zx_alt(:, 4)))

            dX(:) = &
               c_alt(1)*dX_zx_alt(:, 1) + &
               (c_alt(2)*dX_zx_alt(:, 2) + &
               (c_alt(3)*dX_zx_alt(:, 3) + &
               c_alt(4)*dX_zx_alt(:, 4)))
            
            dcoef_dZ = 1d0/dZ_alt
            dcdZ_alt = 0
            dcdZ_alt(1) = -(3*coef*coef-4*coef+1)/2*dcoef_dZ
            dcdZ_alt(2) = (9*coef*coef-10*coef)/2*dcoef_dZ
            dcdZ_alt(3) = -(9*coef*coef-8*coef-1)/2*dcoef_dZ
            dcdZ_alt(4) = coef*(3*coef-2)/2*dcoef_dZ

            dZ(:) = &
               dcdZ_alt(1)*res_zx_alt(:,1) + &
               dcdZ_alt(2)*res_zx_alt(:,2) + &
               dcdZ_alt(3)*res_zx_alt(:,3) + &
               dcdZ_alt(4)*res_zx_alt(:,4)

            if (is_bad(dZ(1))) then
               call mesa_error(__FILE__,__LINE__,'Get1_eosdt_Results bad dZ; 4')
            end if

         endif
      
      else
         write(*, *) 'error: Get1_eosdt_Results assumes eos_integer_ctrl(2) == 1 or 2.'
         call mesa_error(__FILE__,__LINE__)
      end if

   contains

      subroutine do_interp2(iz1, iz2, ierr)
         integer, intent(in) :: iz1, iz2
         integer, intent(out) :: ierr
         real(dp) :: Z1, Z2
         include 'formats'
         ierr = 0
         Z1 = xz% Zs(iz1)
         Z2 = xz% Zs(iz2)
         c(2) = (Z - Z1) / (Z2 - Z1)
         c(1) = 1d0 - c(2)
         dcdZ(2) = 1d0/(Z2 - Z1)
         dcdZ(1) = -dcdZ(2)
         call Get1_eosdt_for_X( &
            rq, which_eosdt, xz, iz1, X, logRho, logT, &
            res_zx(:,1), dlnd_zx(:,1), dlnT_zx(:,1), dX_zx(:,1), &
            ierr)
         if (ierr /= 0) return
         call Get1_eosdt_for_X( &
            rq, which_eosdt, xz, iz2, X, logRho, logT, &
            res_zx(:,2), dlnd_zx(:,2), dlnT_zx(:,2), dX_zx(:,2), &
            ierr)
         if (ierr /= 0) return
      end subroutine do_interp2

   end subroutine Get1_eosdt_Results

   subroutine get1_for_eosdt( &
      handle, which_eosdt, Z, X, abar, zbar, &
      species, chem_id, net_iso, xa, &
      logRho, logT, remaining_fraction, &
      res, d_dlnd, d_dlnT, d_dxa, &
      ierr)
      use chem_def, only: chem_isos
      integer, intent(in) :: handle
      integer, intent(in) :: which_eosdt
      real(dp), intent(in) :: &
         Z, X, abar, zbar, remaining_fraction
      integer, intent(in) :: species
      integer, pointer :: chem_id(:), net_iso(:)
      real(dp), intent(in) :: xa(:)
      real(dp), intent(in) :: logRho, logT
      real(dp), intent(inout), dimension(nv) :: &
         res, d_dlnd, d_dlnT
      real(dp), intent(inout), dimension(nv, species) :: d_dxa
      real(dp), dimension(nv) :: d_dX, d_dZ
      integer, intent(out) :: ierr
      type (EoS_General_Info), pointer :: rq
      type (DT_xz_Info), pointer :: xz
      integer :: i
      rq => eos_handles(handle)

      if (which_eosdt == eosdt_my_eos) then
         xz => my_eos_XZ_struct
      else
         write(*,*) 'unknown which_eosdt supplied'
         ierr = -1
         return
      end if

      call Get1_eosdt_Results( &
         rq, which_eosdt, xz, Z, X, logRho, logT, &
         res, d_dlnd, d_dlnT, d_dX, d_dZ, ierr)

      do i=1,species
         select case(chem_isos% Z(chem_id(i))) ! charge
          case (1) ! X
            d_dxa(:,i) = d_dX
          case (2) ! Y
            d_dxa(:,i) = 0
          case default ! Z
            d_dxa(:,i) = d_dZ
         end select
      end do

   end subroutine get1_for_eosdt

end module custom_eos

