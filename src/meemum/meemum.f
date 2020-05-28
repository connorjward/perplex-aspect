      module calcs

        use, intrinsic :: iso_c_binding
        use, intrinsic :: iso_fortran_env

        implicit none

        include "perplex_parameters.h"

      contains

        ! part 1 of wrapper for meemm (meemum.f)
        ! read the input file (only needs to be called once)
        subroutine init(filename) bind(c)
          character(c_char), dimension(*), intent(in) :: filename

          integer :: i

          character*100 prject,tfname
          common/ cst228 /prject,tfname

          integer iam
          common/ cst4 /iam

          ! ----- PREP WORK -----

          ! prevent input1 crashing by setting project name
          do i = 1, 100
          if (filename(i) == c_null_char) then
            exit
          end if
          prject(i:i) = filename(i)
          end do

          iam = 2
          call vrsion (6)
          call iniprp_wrapper
        end subroutine

        subroutine minimize(p, T, X) bind(c)
          real(c_double), intent(in), value :: T, p
          type(c_ptr), intent(in) :: X
          !type(c_size_t), intent(in), value :: X_len

          real(c_double), dimension(:), pointer :: f_X

          ! source: rlib.f
          integer icomp,istct,iphct,icp
          common/ cst6 /icomp,istct,iphct,icp

          integer :: i

          ! associate C array with Fortran pointer
          !call c_f_pointer(X, f_X, [k5])
          !call c_f_pointer(X, f_X, [icp])
          call c_f_pointer(X, f_X, [10])

          print *, cblk
          print *, f_X

          call abort
          ! set compositions
          !cblk = transfer(f_X, 1.0, k5)
          !cblk = transfer(f_X, 1.0, icp)
          !cblk = f_X
          do i = 1, icp
            cblk(i) = f_X(i)
          end do

          print *, cblk

          ! perform the minimization
          call meemm_wrapper(p, T)
        end subroutine

        !> part 2 of wrapper for meemm (meemum.f)
        ! called at each iteration
        ! pretty much imitates meemm
        subroutine meemm_wrapper(p, T)
          real(c_double), value :: p, T

          integer i

          logical bad

          double precision num

          integer iwt
          common/ cst209 /iwt

          double precision atwt
          common/ cst45 /atwt(k0) 

          double precision v,tr,pr,r,ps
          common/ cst5  /v(l2),tr,pr,r,ps

          integer io3,io4,io9
          common / cst41 /io3,io4,io9

          ! -----------------------------------------

          ! set potential (P, T) values
          ! maybe put in above subroutine
          v(1) = p
          v(2) = T

          ! convert to moles if needed
          if (iwt.eq.1) then 
            do i = 1, jbulk
            cblk(i) = cblk(i)/atwt(i)
            end do 
          end if

          call meemum (bad)

          if (.not.bad) then
            call calpr0 (6)

            if (io3.eq.0) call calpr0 (n3)
         end if 
       end subroutine

       subroutine iniprp_wrapper()
         ! wrapper for iniprp (resub.f)

         ! ----- VARIABLES -----
         ! resub.f : subroutine iniprp
         logical first, err 


         ! ----- WRAPPER -----
         first = .true.

         ! *** prevent prompt for project name ***
         outprt = .true.

         call input1 (first,err)
         call input2 (first)
         call setau1 
         call input9 (first)
         if (lopt(50)) call outsei
         call initlp
       end subroutine

      end module
