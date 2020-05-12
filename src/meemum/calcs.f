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

        !> part 2 of wrapper for meemm (meemum.f)
        ! called at each iteration
        ! pretty much imitates meemm
        subroutine minimize(T, p) bind(c)
          real(c_double), value :: T, p

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
          v(1) = p
          v(2) = T

          ! convert to moles if needed
          if (iwt.eq.1) then 
            do i = 1, jbulk
            cblk(i) = cblk(i)/atwt(i)
            end do 
          end if

          call meemum_wrapper
          !call meemum

          if (.not.bad) then
            call calpr0 (6)

            if (io3.eq.0) call calpr0 (n3)
         end if 
       end subroutine

       ! PRIVATE FUNCTIONS

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

       ! required to do what?
       subroutine meemum_wrapper()
         ! wrapper for meemum (resub.f)

         ! ----- VARIABLES -----

         integer i, idead

         logical nodata, bad

         integer itri(4),jtri(4),ijpt

         double precision wt(3), cum

         double precision v,tr,pr,r,ps
         common/ cst5  /v(l2),tr,pr,r,ps

         integer ipot,jv,iv
         common / cst24 /ipot,jv(l2),iv(l2)

         character*8 vname,xname
         common/ csta2  /xname(k5),vname(l2)

         character*5 cname
         common/ csta4 /cname(k5)

         double precision a,b,c
         common/ cst313 /a(k5,k1),b(k5),c(k1)

         logical gflu,aflu,fluid,shear,lflu,volume,rxn
         common/ cxt20 /gflu,aflu,fluid(k5),shear,lflu,volume,rxn

         integer npt,jdv
         logical fulrnk
         double precision cptot,ctotal
         common/ cst78 /cptot(k19),ctotal,jdv(k19),npt,fulrnk

         integer icomp,istct,iphct,icp
         common/ cst6  /icomp,istct,iphct,icp
         ! ----- PREP WORK -----

         ! ----- WRAPPER -----
         ! normalize vector
         ctotal = 0d0
         do i = 1, icp
         ctotal = ctotal + cblk(i)
         end do 

         do i = 1, icp
         b(i) = cblk(i)/ctotal
         end do
         call incdp0

         call lpopt0 (idead)

         if (idead.eq.0) then
         call getloc (itri,jtri,ijpt,wt,nodata)

         bad = .false.

        else 

        bad = .true.

      end if
      end subroutine


      end module
