# `mesa_custom_EoS`
`mesa_custom_EoS` is an extension to `MESA` to implement user-defined equations of state. It uses the `run_star_extras` functionality.

## Implementing your own EoS

> [!NOTE]
> Prerequisites: The EoS you want to implement must be formatted like the other `MESA` EoS tables.
### Using the placerholder id
`mesa_custom_EoS` comes with a placeholder id for user-defined EoS. To use it, you need to do the following:
1. Add your EoS tables to the `src/data` folder of your work directory.
2. Either rename the EoS folder inside `src/data` to `my_eosDT` and the EoS tables to `my_eosDT_<Z>z<X>x.data`, where `<Z>` and `<X>` are the $Z$ and $X$ values of the EoS table in percentage (e.g., `my_eosDT_50z20x.data`), or change `data_dir_name` and `data_file_name` in `custom_eos.f90` (line 356 and 357) to the name of your EoS folder and the name of your EoS tables, respectively.
3. Run `./mk` to compile the code and make sure that `eos_integer_ctrl(1)` is set to `5` in your `inlist`.
### Creating your own EoS entry
To implement your own EoS, you need to do the 
following:
1. Add the folder that contains your EoS tables to the `src/data` folder of your work directory.
Then inside `custom_eos.f90`:
1. Define the two parameters `num_<eos_name>_Zs` and `num_<eos_name>_Xs` that contain the number of $X$ and $Z$ values in your EoS tables. For example, if your EoS tables range from $X = 0$ to $X = 1$ in steps of $0.1$ and from $Z = 0$ to $Z = 1$ in steps of $0.1$, you would define:
```Fortran
integer, parameter :: num_<eos_name>_Zs = 11
integer, parameter :: num_<eos_name>_Xs = 11
```
1. Define the array `<eos_name>_XZ_loaded` of boolean parameters with dimension `(eos_<eos_name>_Xs, eos_<eos_name>_Zs)`. This array will be used to check if the EoS tables have been loaded. For example:
```Fortran
logical, dimension(num_<eos_name>_Xs, num_<eos_name>_Zs) :: <eos_name>_XZ_loaded
```
1. Define a `DT_XZ_Info` type target that contains the information about the $(X,Z)$ grid of your EoS tables. Sticking with the example above, you would define:
```Fortran
type (DT_XZ_Info), target :: <eos_name>_XZ_struct
```
1. Define another target of type `EosDT_XZ_Info` of the dimensions `(num_<eos_name>_Zs, num_<eos_name>_Xs)`. This object will store the EoS table data. Again, for the example above:
```Fortran
type (EosDT_XZ_Info), dimension(num_<eos_name>_Xs, num_<eos_name>_Zs), target :: qeos_<eos_name>_data
```
1. Define an id for your EoS with the name `eosdt_<eos_name>`. This is the integer that you will use in the `inlist` to select your EoS. For example:
```Fortran
integer, parameter :: eosdt_<eos_name> = 4
```
1. Inside the `eos_init_custom_eos subroutine`, add a `pointer` of type `DT_XZ_Info` called `<eos_name>_ptr`. It will be used to populate `<eos_name>_XZ_struct`. In addition, set `<eos_name>_XZ_loaded` to `.false.`. The `DT_XZ_Info` type has the following properties: `nZs` stores the number of $Z$ values in your EoS tables, `Zs` is a 1D array that contains these $Z$ values, `nXs_for_Z` is a 1D array that contains the number of $X$ values for each $Z$ value and `Xs_for_Z` is a 2D array that contains the $X$ values for each $Z$ value. Check out the definition of `qeos_sio2` inside the subroutine for explicit numerical values.

For the example from above, we can make use of the already defined `qeos_sio2_XZ_struct` variable:
```Fortran

! pointer to the EoS XZ structure
type (DT_XZ_Info), pointer :: qeos_sio2_XZ_ptr, qeos_h2o_XZ_ptr, cms_qeos_h2o_XZ_ptr, <eos_name>_XZ_ptr

! #Z values in the EoS XZ structure
real(dp) :: qeos_sio2_spacing(1:11) = &
	(/ 0.00d0, 0.1d0, 0.2d0, 0.3d0, 0.4d0, 0.5d0, 0.6d0, &
	0.7d0, 0.8d0, 0.9d0, 1.0d0 /)

! tracking if a table was already loaded
qeos_sio2_XZ_loaded(:,:)=.false.
qeos_h2o_XZ_loaded(:,:)=.false.
cms_qeos_h2o_XZ_loaded(:,:)=.false.
<eos_name>_XZ_loaded(:,:)=.false.

! < ... code populating the other "struct" variables ... >

!> <eos_name>
<eos_name>_XZ_ptr => <eos_name>_XZ_struct
<eos_name>_ptr % nZs = qeos_sio2_XZ_ptr % nZs
<eos_name>_ptr % Zs(1: cms_qeos_h2o_XZ_ptr % nZs) = qeos_sio2_XZ_ptr % Zs(1: <eos_name>_ptr % nZs)
<eos_name>_ptr % nXs_for_Z(1: cms_qeos_h2o_XZ_ptr % nZs) = qeos_sio2_XZ_ptr % nXs_for_Z(1: qeos_sio2_XZ_ptr % nZs)
<eos_name>_XZ_ptr % Xs_for_Z(1: <eos_name>_XZ_ptr % nXs_for_Z(1), 1: <eos_name>_XZ_ptr % nZs) = &
qeos_sio2_XZ_ptr % Xs_for_Z(1: qeos_sio2_XZ_ptr % nXs_for_Z(1), 1: qeos_sio2_XZ_ptr % nZs)
```
8. Now, go through the rest of the code and anywhere the code checks for `which_eosdt`, add an `else if` statement for your EoS. For example:
```Fortran
if (which_eosdt == eosdt_qeos_sio2) then
    xz => qeos_sio2_XZ_struct
else if (which_eosdt == eosdt_qeos_h2o) then
    xz => qeos_h2o_XZ_struct
else if (which_eosdt == eosdt_cms_qeos_h2o) then
    xz => cms_qeos_h2o_XZ_struct
else if (which_eosdt == eosdt_<eos_name>) then
	xz => <eos_name>_XZ_struct
else
    write(*,*) 'unknown which_eosdt supplied'
    ierr = -1
    return
end if
```
9. Lastely, run `./mk` to compile the code and change the `eos_integer_ctrl(1)` to the integer you defined in step 6, i.e. `eosdt_<eos_name>`.