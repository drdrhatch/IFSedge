program skim

        use file_utils

        implicit none

        real (kind=8), parameter :: pi=3.141592653589793
        complex (kind=8), parameter :: zi=(0.0,1.0)

        integer (kind=8) :: nkz, nky, nkperp
        real (kind=8) :: dkz, dky, dkperp
        integer (kind=8) :: ikz, iky, ikperp
        real (kind=8) :: ky0, kz0, kperp0
        real (kind=8) :: kperp1, ky1, kz1
        real (kind=8), dimension(:), allocatable :: ky_grid, kz_grid, kperp_grid

        integer (kind=8) :: nfprime, ntprime_i, ntprime_e
        integer (kind=8) :: nomd
        real (kind=8) :: dfprime, dtprime_i, dtprime_e
        real (kind=8) :: domd
        integer (kind=8) :: ifprime, itprime_i, itprime_e
        integer (kind=8) ::  iomd
        real (kind=8) :: tprime0_i, fprime0, tprime0_e
        real (kind=8) :: omd0
        real (kind=8), dimension(:), allocatable :: fprime_grid, tprime_i_grid
        real (kind=8), dimension(:), allocatable :: tprime_e_grid
        real (kind=8), dimension(:), allocatable :: omd_grid
        real (kind=8) :: tprime_i, tprime_e

        real (kind=8) :: vy_1, vy_2, vz_1, vz_2
        real (kind=8) :: fprime, tprime, omd, ky, kz, kperp
        real (kind=8) :: vi, int_tol,sec_tol, fr, na_e, na_z, na_i
        real (kind=8) :: omde, eps, nu0, coll_c, f_t, t_e, theta_c, theta_v
        real (kind=8) :: omde_avg, lambda_debye2
        real (kind=8) :: Ti, Zeff, Z, mu_e, mu_z, fr_te
        complex (kind=8) :: seed1, seed2, seed3
        real (kind=8) :: om1_re, om1_im
        integer (kind=8) :: sec_nsteps, sec_istep
        integer (kind=8) :: int_nsteps, int_istep

        integer(kind=8), dimension(:),allocatable :: ikperp_ref, ikz_ref, iky_ref

        complex (kind=8) :: omega 
        complex (kind=8) :: root, root_ikperp_ref, root_ikzx_ref, root_ikyzx_ref
        complex (kind=8) :: root_ikperp_1, root_ikzx_1, root_ikyzx_1
        complex (kind=8) :: root_tmp_fl 
        integer :: gamma_unit=101, Dmixing_unit=102
        integer :: out_unit=103
        real(kind=8) :: cut_buffer, seed_fr, lower_bnd
        
        call init_file_utils
        call read_input_file
        call init_grids
        call init_output_files

        do iomd = 1, nomd
           omd = omd_grid(iomd)
           print*, 'omd'
           print*, omd
           do ifprime = 1, nfprime
              fprime = fprime_grid(ifprime)
              print*,'omn'
              print*, fprime
              do itprime_i = 1, ntprime_i
                 tprime_i = tprime_i_grid(itprime_i)
                 print*,'omt_i'
                 print*, tprime_i
                    do itprime_e = 1, ntprime_e
 !                      tprime_e = tprime_e_grid(itprime_e)
                       tprime_e = tprime_i
                       print*,'omt_e'
                       print*,tprime_e
                       call ky_scan
                    end do
              end do
           end do
        end do

!        call close_output_file (gamma_unit)
!        call close_output_file (Dmixing_unit)
        call close_output_file (out_unit)
        call finish_file_utils

contains

 subroutine ky_scan 

        implicit none

        complex (kind=8), dimension(:), allocatable :: omega_ky_scan
        real (kind=8) :: kz_gam,kz_Dmix,Dmixing_kz_max
        complex (kind=8) :: omega_kz_max

        allocate(omega_ky_scan(nky))
        omega_ky_scan = 0. -9999.*zi
        iky_ref = minloc(abs(ky_grid-ky1))
        !do iky = 1, nky
        !   ky = ky_grid(iky)
        !   call kz_scan(kz_gam,omega_kz_max,kz_Dmix,Dmixing_kz_max)
        !   omega_ky_scan(iky) = omega_kz_max
        !   if (iky==iky_ref(1).and.aimag(omega_kz_max)>aimag(root_ikyzx_ref)) &
        !      root_ikyzx_ref = omega_kz_max
        !end do

        do iky = iky_ref(1), nky
           ky = ky_grid(iky)
           call kz_scan(kz_gam,omega_kz_max,kz_Dmix,Dmixing_kz_max)
           if (aimag(omega_kz_max)>aimag(omega_ky_scan(iky))) &
             omega_ky_scan(iky) = omega_kz_max
           !if (.not.aimag(omega_ky_scan(iky))==-9999.0) &
           !  write (gamma_unit,&
           !  '(10e12.4)') tprime_i,fprime,tprime_e,omd, &
           !  kperp,ky,kz,omega_ky_scan(iky)
        end do

        do iky = iky_ref(1),1,-1
           ky = ky_grid(iky)
           call kz_scan(kz_gam,omega_kz_max,kz_Dmix,Dmixing_kz_max)
           if (aimag(omega_kz_max)>aimag(omega_ky_scan(iky))) &
             omega_ky_scan(iky) = omega_kz_max
           !if (.not.aimag(omega_ky_scan(iky))==-9999.0) &
           !  write (gamma_unit,&
           !  '(10e12.4)') tprime_i,fprime,tprime_e,omd, &
           !  kperp,ky,kz,omega_ky_scan(iky)
        end do

        root_ikyzx_ref = omega_ky_scan(iky_ref(1))
        root_ikyzx_1 = omega_ky_scan(1)

 end subroutine

 subroutine kz_scan (kz_maxgam,omega_maxgam,kz_maxDmix,Dmixing_max)

        implicit none

        real (kind=8),intent(out) :: kz_maxgam, kz_maxDmix, Dmixing_max
        complex (kind=8),intent(out) :: omega_maxgam
        integer (kind=8), dimension(:), allocatable :: ikz_maxgamma
        integer (kind=8), dimension(:), allocatable :: ikz_maxDmixing
        real (kind=8), dimension(:), allocatable :: Dmixing_kz_scan
        complex (kind=8), dimension(:), allocatable :: omega_kz_scan
        real (kind=8) :: kperp_gam,kperp_Dmix,Dmixing_kperp_max
        complex (kind=8) :: omega_kperp_max

        allocate(Dmixing_kz_scan(nkz),omega_kz_scan(nkz))
        allocate(ikz_maxgamma(1),ikz_maxDmixing(1))
        Dmixing_kz_scan = -9999.0
        omega_kz_scan = 0. -9999.*zi
        ikz_ref = minloc(abs(kz_grid-kz1))
        !do ikz = 1, nkz
        !   kz = kz_grid(ikz)
        !   print*, "ky, kz"
        !   print*, ky, kz
        !   call kperp_scan(kperp_gam,omega_kperp_max,kperp_Dmix,Dmixing_kperp_max)
        !   omega_kz_scan(ikz) = omega_kperp_max
        !   Dmixing_kz_scan(ikz) = Dmixing_kperp_max
        !   if (ikz==ikz_ref(1).and.aimag(omega_kperp_max)>aimag(root_ikzx_ref)) &
        !      root_ikzx_ref = omega_kperp_max
        !end do

        do ikz = ikz_ref(1), nkz
           kz = kz_grid(ikz)
           print*, "ky, kz"
           print*, ky, kz
           call kperp_scan(kperp_gam,omega_kperp_max,kperp_Dmix,Dmixing_kperp_max)
           if (aimag(omega_kperp_max)>aimag(omega_kz_scan(ikz))) &
              omega_kz_scan(ikz) = omega_kperp_max
           if (Dmixing_kperp_max>Dmixing_kz_scan(ikz)) &
              Dmixing_kz_scan(ikz) = Dmixing_kperp_max
          if (.not.aimag(omega_kz_scan(ikz))==-9999.0) &
            write (gamma_unit,&
            '(9e12.4)') tprime_i,fprime,tprime_e,omd*ky, &
            kperp,ky,kz,omega_kz_scan(ikz)
        end do
        write (gamma_unit,*)
!        write (Dmixing_unit,*)

        do ikz = ikz_ref(1),1,-1
           kz = kz_grid(ikz)
           print*, "ky, kz"
           print*, ky, kz
           call kperp_scan(kperp_gam,omega_kperp_max,kperp_Dmix,Dmixing_kperp_max)
           if (aimag(omega_kperp_max)>aimag(omega_kz_scan(ikz))) &
              omega_kz_scan(ikz) = omega_kperp_max
           if (Dmixing_kperp_max>Dmixing_kz_scan(ikz)) &
              Dmixing_kz_scan(ikz) = Dmixing_kperp_max
          if (.not.aimag(omega_kz_scan(ikz))==-9999.0) &
            write (gamma_unit,&
            '(9e12.4)') tprime_i,fprime,tprime_e,omd*ky, &
            kperp,ky,kz,omega_kz_scan(ikz)
        end do
        write (gamma_unit,*)
!        write (Dmixing_unit,*)
        
        ikz_maxgamma = maxloc(aimag(omega_kz_scan))
        ikz_maxDmixing = maxloc(Dmixing_kz_scan)
        kz_maxgam = kz_grid(ikz_maxgamma(1))
        kz_maxDmix = kz_grid(ikz_maxDmixing(1))
        omega_maxgam = omega_kz_scan(ikz_maxgamma(1))
        Dmixing_max = Dmixing_kz_scan(ikz_maxDmixing(1))

        root_ikzx_ref = omega_kz_scan(ikz_ref(1))
        root_ikzx_1 = omega_kz_scan(1)

 end subroutine

 subroutine kperp_scan (kperp_maxgam,omega_maxgam,kperp_maxDmix,Dmixing_max)

        implicit none

        real (kind=8), intent(out) :: kperp_maxgam, kperp_maxDmix, Dmixing_max
        complex (kind=8), intent(out) :: omega_maxgam
        integer (kind=8), dimension(:), allocatable :: ikperp_maxgamma
        integer (kind=8), dimension(:), allocatable :: ikperp_maxDmixing
        real (kind=8), dimension(:), allocatable :: Dmixing_kperp_scan
        complex (kind=8), dimension(:), allocatable :: omega_kperp_scan
        complex (kind=8) :: root_tmp,om_linr,om_linr_a,om_quadr,om_quadr_a
        complex (kind=8) :: root_fl
        real (kind=8) :: Dmixing_tmp


          allocate(Dmixing_kperp_scan(nkperp), omega_kperp_scan(nkperp))
          allocate(ikperp_maxgamma(1),ikperp_maxDmixing(1))
          Dmixing_kperp_scan = -9999.0
          omega_kperp_scan = 0. -9999.*zi
          ikperp_ref = minloc(abs(kperp_grid-kperp1))
          do ikperp = 1, nkperp
             kperp = kperp_grid(ikperp)*ky
             print*, "kperp"
             print*, kperp
             call init_roots
             print*, 'init_roots'
             print*, seed1, seed2, seed3
             call rootfinder_muller(seed2,seed1,seed3,root_tmp)

             !call recur_linear_soln(root_tmp)
             !call recur_quadratic_soln(root_tmp)
             if (sec_istep==sec_nsteps.or.isnan(real(root_tmp))) then
              print*, 'No root is found.'
              root = seed1
              omega_kperp_scan(ikperp) = 0.-9999.*zi
             !else if (abs(aimag(root_tmp))<lower_bnd) then
             ! print*, 'No root is found.'
             ! root = seed1
             ! omega_kperp_scan(ikperp) = 0. -9999.*zi
             else if (aimag(root_tmp)==aimag(root_tmp)-1..or.&
                      isnan(aimag(root_tmp))) then
              print*, 'No root is found.'
              root = seed1
              omega_kperp_scan(ikperp) = 0.-9999.*zi
             else
              print*, "Root is found:"
              print*, root_tmp    
              omega_kperp_scan(ikperp) = root_tmp
              root = root_tmp
             end if
             if (ikperp==ikperp_ref(1).and.aimag(omega_kperp_scan(ikperp))>aimag(root_ikperp_ref)) &
              root_ikperp_ref = omega_kperp_scan(ikperp)
              !!call quadratic_solution(om_quadr,om_quadr_a)
              !call linear_solution(om_linr,om_linr_a)
              !call init_roots
              !seed1 = root_tmp_fl
              !!seed1 = om_quadr
              !!seed2 = seed1*(1.-seed_fr)
              !!seed3 = seed1*(1.+seed_fr)
              !!print*, 'init_roots analytic'
              !!print*, seed1, seed2, seed3
              !!call analytic_rf(seed2,seed1,seed3,root_tmp_fl)
              !!root_fl = root_tmp_fl

              if (.not.aimag(omega_kperp_scan(ikz))==-9999.0) &
              write (out_unit, '(9e12.4)') tprime_i,fprime,tprime_e,&
                             ky,kz,omega_kperp_scan(ikperp),kperp,omd*ky!,&
                             !root_fl,&
                             !om_quadr,om_quadr_a
          end do
          do ikperp = ikperp_ref(1), 1, -1
             kperp = kperp_grid(ikperp)*ky
             print*, "kperp"
             print*, kperp
             call init_roots
             print*, 'init_roots'
             print*, seed1, seed2, seed3
             call rootfinder_muller(seed1,seed2,seed3,root_tmp)

             !call recur_linear_soln(root_tmp)
             !call recur_quadratic_soln(root_tmp)
             if (sec_istep==sec_nsteps.or.isnan(real(root_tmp))) then
              print*, 'No root is found.'
              root = seed1
              omega_kperp_scan(ikperp) = 0.-9999.*zi
             !else if (abs(aimag(root_tmp))<lower_bnd) then
             ! print*, 'No root is found.'
            ! root = seed1
             ! omega_kperp_scan(ikperp) = 0. -9999.*zi
             else if (aimag(root_tmp)==aimag(root_tmp)-1..or.&
                      isnan(aimag(root_tmp))) then
              print*, 'No root is found.'
              root = seed1
              omega_kperp_scan(ikperp) = 0.-9999.*zi
             else
              print*, "Root is found:"
              print*, root_tmp    
              root = root_tmp
              if (aimag(root_tmp)>aimag(omega_kperp_scan(ikperp))) &
                omega_kperp_scan(ikperp) = root_tmp
              !!call quadratic_solution(om_quadr,om_quadr_a)
              !call linear_solution(om_linr,om_linr_a)
              !!call init_roots
              !!print*, 'init_roots analytic'
              !!print*, seed1, seed2, seed3
              !!call analytic_rf(seed2,seed1,seed3,root_tmp_fl)
              if (.not.aimag(omega_kperp_scan(ikz))==-9999.0) &
              write (out_unit, '(9e12.4)') tprime_i,fprime,tprime_e,&
                             ky,kz,omega_kperp_scan(ikperp),kperp,omd*ky!!,&
                             !!root_fl,&
                             !!om_quadr,om_quadr_a
             end if
          end do
          write (out_unit,*) 
          do ikperp = ikperp_ref(1), nkperp
             kperp = kperp_grid(ikperp)*ky
             print*, "kperp"
             print*, kperp
             call init_roots
             print*, 'init_roots'
             print*, seed1, seed2, seed3
             call rootfinder_muller(seed2,seed1,seed3,root_tmp)
!
 !            !call recur_linear_soln(root_tmp)
 !            !call recur_quadratic_soln(root_tmp)
             if (sec_istep==sec_nsteps.or.isnan(real(root_tmp))) then
              print*, 'No root is found.'
              root = seed1
              omega_kperp_scan(ikperp) = 0.-9999.*zi
             !else if (abs(aimag(root_tmp))<lower_bnd) then
             ! print*, 'No root is found.'
             ! root = seed1
             ! omega_kperp_scan(ikperp) = 0.-9999.*zi
             else if (aimag(root_tmp)==aimag(root_tmp)-1..or.&
                      isnan(aimag(root_tmp))) then
              print*, 'No root is found.'
              root = seed1
              omega_kperp_scan(ikperp) = 0.-9999.*zi
             else
              print*, "Root is found:"
              print*, root_tmp    
              root = root_tmp
              if (aimag(root_tmp)>aimag(omega_kperp_scan(ikperp))) &
                omega_kperp_scan(ikperp) = root_tmp
              !!call quadratic_solution(om_quadr,om_quadr_a)
              !call linear_solution(om_linr,om_linr_a)
              !!call init_roots
              !!print*, 'init_roots analytic'
              !!print*, seed1, seed2, seed3
              !!call analytic_rf(seed2,seed1,seed3,root_tmp_fl)
              if (.not.aimag(omega_kperp_scan(ikz))==-9999.0) &
              write (out_unit, '(9e12.4)') tprime_i,fprime,tprime_e,&
                             ky,kz,omega_kperp_scan(ikperp),kperp,omd!!,&
                             !!root_tmp_fl,&
                             !!om_quadr,om_quadr_a
             end if
          end do
          write (out_unit, *)

          ikperp_maxgamma = maxloc(aimag(omega_kperp_scan))
          Dmixing_kperp_scan = aimag(omega_kperp_scan)/(kperp_grid**2+ky**2)
          ikperp_maxDmixing = maxloc(Dmixing_kperp_scan)

          kperp_maxgam = kperp_grid(ikperp_maxgamma(1))
          omega_maxgam = omega_kperp_scan(ikperp_maxgamma(1))
          kperp_maxDmix = kperp_grid(ikperp_maxDmixing(1))
          Dmixing_max = Dmixing_kperp_scan(ikperp_maxDmixing(1))

          if (.not.aimag(omega_maxgam)==-9999.0) then  
              root_ikperp_ref = omega_kperp_scan(ikperp_ref(1))
              root_ikperp_1 = omega_kperp_scan(1)
          end if
          deallocate(Dmixing_kperp_scan,omega_kperp_scan)
          deallocate(ikperp_maxgamma,ikperp_maxDmixing)

  end subroutine

subroutine quadratic_solution(omega_quadr, omega_quadr_a)

        implicit none
        
        complex(kind=8), intent(out) :: omega_quadr, omega_quadr_a
        complex(kind=8) :: omega1, omega2
        complex (kind=8) :: term1, term2, term3, term4
        complex (kind=8) :: coeff_a, coeff_b, coeff_c
        integer(kind=8) :: t,s
        
        s = 2   
        t = 1 
        call vy_integral_simpson(omega,vy_1,vy_2,term1,s,t)
        t = 2
        call vy_integral_simpson(omega,vy_1,vy_2,term2,s,t)
        t = 3 
        call vy_integral_simpson(omega,vy_1,vy_2,term3,s,t)

        term4 = -Ti - 1. + term1

        omega1 = (-term2+sqrt(term2**2-4.*term4*term3))/2./term4
        omega2 = (-term2-sqrt(term2**2-4.*term4*term3))/2./term4

        if (aimag(omega1)>aimag(omega2)) then
           omega_quadr = omega1
        else
           omega_quadr = omega2
        end if

        coeff_a = 1.+Ti-(1.-ky**2+ky**4/2.)
        coeff_b = fprime*ky*(1.-ky**2+ky**4/2.) + &
                  tprime_i*ky*(-ky**2+ky**4) - &
                  omd*ky*(2.-3.*ky**2+2.*ky**4)
        coeff_c = fprime*ky*omd*ky*(2.-3.*ky**2+2.*ky**4) + &
                  tprime_i*ky*omd*ky*(2.-6.*ky**2+6.*ky**4)

        omega1 = (-coeff_b+sqrt(coeff_b**2-4.*coeff_a*coeff_c))/2./coeff_a
        omega2 = (-coeff_b-sqrt(coeff_b**2-4.*coeff_a*coeff_c))/2./coeff_a

        if (aimag(omega1)>aimag(omega2)) then
           omega_quadr_a = omega1
        else
           omega_quadr_a = omega2
        end if

end subroutine

subroutine recur_quadratic_soln(omega_quadr_recur)

        implicit none

        complex (kind=8), intent(out) :: omega_quadr_recur
        complex (kind=8) :: this_int_te
        complex (kind=8) :: this_omega, next_omega
        complex (kind=8) :: omega1, omega2
        complex (kind=8) :: coeff_a, coeff_b, coeff_c
        complex (kind=8) :: term1, term2, term3, term4
        integer (kind=8) :: s,t

!        coeff_a = 2.-(1.-ky**2+ky**4/2.)
!        coeff_b = fprime_i*ky*(1.-ky**2+ky**4/2.) + &
!                  tprime_i*ky*(-ky**2+ky**4) - &
!                  omd/A*ky*(2.-3.*ky**2+2.*ky**4)
!        coeff_c = fprime_i*ky*omd/A*ky*(2.-3.*ky**2+2.*ky**4) + &
!                  tprime_i*ky*omd/A*ky*(2.-6.*ky**2+6.*ky**4)

!        omega1 = (-coeff_b+sqrt(coeff_b**2-4.*coeff_a*coeff_c))/2./coeff_a
!        omega2 = (-coeff_b-sqrt(coeff_b**2-4.*coeff_a*coeff_c))/2./coeff_a

!        if (aimag(omega1)>aimag(omega2)) then
!           this_omega = omega1
!        else
!           this_omega = omega2
!        end if

        this_omega = 0.
        print*,'this_omega'
        print*, this_omega

        s = 2   
        t = 1 
        call vy_integral_simpson(this_omega,vy_1,vy_2,term1,s,t)
        t = 2
        call vy_integral_simpson(this_omega,vy_1,vy_2,term2,s,t)
        t = 3 
        call vy_integral_simpson(this_omega,vy_1,vy_2,term3,s,t)

        term4 = -Ti - 1. + term1

        omega1 = (-term2+sqrt(term2**2-4.*term4*term3))/2./term4
        omega2 = (-term2-sqrt(term2**2-4.*term4*term3))/2./term4

        if (aimag(omega1)>aimag(omega2)) then
           next_omega = omega1
        else
           next_omega = omega2
        end if
        print*,'next_omega'
        print*,next_omega

        if (abs(this_omega-next_omega)>1.E-08) then
                this_omega = next_omega
                print*, 'omega_quad_recur'
                print*, this_omega

                s = 4
                call vy_integral(this_omega,this_int_te,s)

                s = 2   
                t = 1 
                call vy_integral_simpson(this_omega,vy_1,vy_2,term1,s,t)
                t = 2
                call vy_integral_simpson(this_omega,vy_1,vy_2,term2,s,t)
                t = 3 
                call vy_integral_simpson(this_omega,vy_1,vy_2,term3,s,t)

                term4 = -Ti - 1. + term1 + t_e*this_int_te

                omega1 = (-term2+sqrt(term2**2-4.*term4*term3))/2./term4
                omega2 = (-term2-sqrt(term2**2-4.*term4*term3))/2./term4

                if (aimag(omega1)>aimag(omega2)) then
                   next_omega = omega1
                else
                   next_omega = omega2
                end if

        end if

        omega_quadr_recur = next_omega
        
end subroutine

subroutine linear_solution(omega_linr, omega_linr_a)

        implicit none
        
        complex(kind=8), intent(out) :: omega_linr, omega_linr_a
        complex(kind=8) :: omega1, omega2
        complex (kind=8) :: term1, term2, term3
        integer(kind=8) :: t,s
        
        s = 2   
        t = 1 
        call vy_integral_simpson(omega,vy_1,vy_2,term1,s,t)
        t = 2
        call vy_integral_simpson(omega,vy_1,vy_2,term2,s,t)

        term3 = Ti + 1. - term1

        omega1 = term2/term3

        omega_linr = omega1

        omega_linr_a = -(fprime*ky*(1-ky**2+ky**4/2.)+&
                       tprime_i*ky*(-ky**2+ky**4))/(Ti+ky**2-ky**4/2.)

end subroutine

subroutine recur_linear_soln(omega_linr_recur)

        implicit none

        complex (kind=8), intent(out) :: omega_linr_recur
        complex (kind=8) :: this_int_te
        complex (kind=8) :: this_omega, next_omega
        integer (kind=8) :: s

        this_omega = 0.
        next_omega = -(fprime*ky*(1-ky**2+ky**4/2.)+&
                tprime_i*ky*(-ky**2+ky**4))/(Ti+ky**2-ky**4/2.)

        if (abs(this_omega-next_omega)>sec_tol) then
                this_omega = next_omega
                print*, 'omega_linear_recur'
                print*, this_omega
                s = 4
                call vy_integral(this_omega,this_int_te,s)
                next_omega = -(fprime*ky*(1-ky**2+ky**4/2.)+&
                             tprime_i*ky*(-ky**2+ky**4))/&
                             (Ti+ky**2-ky**4/2.-t_e*this_int_te)
        end if

        omega_linr_recur = next_omega
        
end subroutine

 subroutine init_roots
 
        implicit none

        complex(kind=8) :: omega_tmp

        !if ((.not.t_e==0.).and.ky<0.5) then
        !   if (omd==1.) then
        !      call recur_quadratic_soln(omega_tmp)
        !   else
        !      call recur_linear_soln(omega_tmp)
        !   end if
        !   seed1 = omega_tmp
        !else if (iomd==1.and.iky==1.and.ikz==1.and.ikperp==1) then
        if (iomd==1.and.itprime_i==1.and.ifprime==1.and.iky==iky_ref(1).and.ikz==ikz_ref(1).and.ikperp==1) then
           seed1 = om1_re + om1_im*zi
        !else if (itprime_i.and.ifprime==1.and.iky==1.and.ikz==1.and.ikperp==1.and.aimag(root_ikyzx_1)>aimag(root)) then
        !   seed1 = root_ikyzx_1
        !else if (itprime_i.and.ifprime==1.and.ikz==1.and.ikperp==1.and.aimag(root_ikzx_1)>aimag(root)) then
        !   seed1 = root_ikzx_1
        !else if (itprime_i.and.ifprime==1.and.ikperp==1.and.aimag(root_ikperp_1)>aimag(root)) then
        !   seed1 = root_ikperp_1
        !else if (itprime_i.and.ifprime==1.and.iky==iky_ref(1).and.ikz==ikz_ref(1).and.ikperp==ikperp_ref(1).and.aimag(root_ikyzx_ref)>aimag(root)) then
        !   seed1 = root_ikyzx_ref
        !else if (itprime_i.and.ifprime==1.and.ikz==ikz_ref(1).and.ikperp==ikperp_ref(1).and.aimag(root_ikzx_ref)>aimag(root)) then
        !   seed1 = root_ikzx_ref
        !else if (itprime_i.and.ifprime==1.and.ikperp==ikperp_ref(1).and.aimag(root_ikperp_ref)>aimag(root)) then
        !   seed1 = root_ikperp_ref
        else 
           seed1 = root
        end if

        if (aimag(seed1)>0.) then
           seed2 = seed1*(1.0-seed_fr)
           seed3 = seed1*(1.0+seed_fr)
        else
           seed2 = seed1*(1.0+seed_fr)
           seed3 = seed1*(1.0-seed_fr)
        end if

 end subroutine

 subroutine rootfinder_secant(sd1,sd2,rt)

        implicit none

        complex (kind=8), intent(in) :: sd1,sd2
        complex (kind=8), intent(out) :: rt
        complex (kind=8) :: x_2,f_2,x_1,f_1
        complex (kind=8) :: x_tmp, f_tmp

        x_1 = sd1
        call dispersion_relation(x_1,f_1)
        x_2 = sd2
        call dispersion_relation(x_2,f_2)
        print*, 'f_1, f_2'
        print*, f_1, f_2
        x_tmp = x_1
        sec_istep = 0
        do while (abs(f_1)>sec_tol .and. abs(f_2)>sec_tol .and. sec_istep<sec_nsteps) 
            x_tmp = x_1 - f_1 * ((x_1-x_2)/(f_1-f_2))
                call dispersion_relation(x_tmp,f_tmp)
                if (isnan(aimag(f_tmp)).or.isnan(real(f_tmp))) then
                   sec_istep = sec_nsteps
                   exit
                end if
                print*, 'root finder'
                print*, sec_istep, x_tmp, f_tmp
                f_2 = f_1
                f_1 = f_tmp
                x_2 = x_1
                x_1 = x_tmp
                sec_istep = sec_istep +1
        end do
        rt = x_tmp

  end subroutine

 subroutine rootfinder_muller(sd0,sd1,sd2,rt)

        implicit none

        complex (kind=8), intent(in) :: sd0, sd1, sd2
        complex (kind=8), intent(out) :: rt
        complex (kind=8) :: x_0,x_1,x_2,f_0,f_1,f_2
        complex (kind=8) :: h_1,h_2,g_1,g_2
        complex (kind=8) :: coeff_b,h_tmp_p,h_tmp_m
        complex (kind=8) :: x_tmp,f_tmp,h_tmp,g_tmp

        x_0 = sd0
        call dispersion_relation(x_0,f_0)
        x_1 = sd1
        call dispersion_relation(x_1,f_1)
        x_2 = sd2
        call dispersion_relation(x_2,f_2)
        h_1 = x_1 - x_0
        h_2 = x_2 - x_1
        g_1 = (f_1 - f_0)/h_1
        !g_2 = (f_2 - f_1)/h_2
        sec_istep = 0
        x_tmp = x_1
        do while (abs(f_1)>sec_tol.and.abs(f_2)>sec_tol.and.sec_istep<sec_nsteps)
           g_2 = (f_2 - f_1)/h_2
           g_tmp = (g_2 - g_1)/(h_2 + h_1)
           coeff_b = g_2 + h_2*g_tmp
           h_tmp_p = (-2.*f_2)/(coeff_b+sqrt(coeff_b**2-4.*g_tmp*f_2))
           h_tmp_m = (-2.*f_2)/(coeff_b-sqrt(coeff_b**2-4.*g_tmp*f_2))
           if (abs(h_tmp_p)>abs(h_tmp_m)) then
              h_tmp = h_tmp_m
           else
              h_tmp = h_tmp_p
           end if
           x_tmp = x_2 + h_tmp
           call dispersion_relation(x_tmp,f_tmp)
           print*, 'root finder'
           print*, sec_istep, x_tmp, f_tmp
           if (isnan(aimag(f_tmp)).or.isnan(real(f_tmp))) then
              sec_istep = sec_nsteps
              exit
           end if
           x_0 = x_1
           x_1 = x_2
           x_2 = x_tmp
           f_0 = f_1
           f_1 = f_2
           f_2 = f_tmp
           h_1 = h_2
           h_2 = h_tmp
           g_1 = g_2
           sec_istep = sec_istep + 1
        end do
        rt = x_tmp

 end subroutine

 subroutine analytic_rf(sd0,sd1,sd2,rt)

        implicit none

        complex (kind=8), intent(in) :: sd0, sd1, sd2
        complex (kind=8), intent(out) :: rt
        complex (kind=8) :: x_0,x_1,x_2,f_0,f_1,f_2
        complex (kind=8) :: h_1,h_2,g_1,g_2
        complex (kind=8) :: coeff_b,h_tmp_p,h_tmp_m
        complex (kind=8) :: x_tmp,f_tmp,h_tmp,g_tmp

        x_0 = sd0
        call analytic_dr(x_0,f_0)
        x_1 = sd1
        call analytic_dr(x_1,f_1)
        x_2 = sd2
        call analytic_dr(x_2,f_2)
        h_1 = x_1 - x_0
        h_2 = x_2 - x_1
        g_1 = (f_1 - f_0)/h_1
        !g_2 = (f_2 - f_1)/h_2
        sec_istep = 0
        x_tmp = x_1
        do while (abs(f_1)>sec_tol.and.abs(f_2)>sec_tol.and.sec_istep<sec_nsteps)
           g_2 = (f_2 - f_1)/h_2
           g_tmp = (g_2 - g_1)/(h_2 + h_1)
           coeff_b = g_2 + h_2*g_tmp
           h_tmp_p = (-2.*f_2)/(coeff_b+sqrt(coeff_b**2-4.*g_tmp*f_2))
           h_tmp_m = (-2.*f_2)/(coeff_b-sqrt(coeff_b**2-4.*g_tmp*f_2))
           if (abs(h_tmp_p)>abs(h_tmp_m)) then
              h_tmp = h_tmp_m
           else
              h_tmp = h_tmp_p
           end if
           x_tmp = x_2 + h_tmp
           call analytic_dr(x_tmp,f_tmp)
           print*, 'analytic root finder'
           print*, sec_istep, x_tmp, f_tmp
           if (isnan(aimag(f_tmp)).or.isnan(real(f_tmp))) then
              sec_istep = sec_nsteps
              exit
           end if
           x_0 = x_1
           x_1 = x_2
           x_2 = x_tmp
           f_0 = f_1
           f_1 = f_2
           f_2 = f_tmp
           h_1 = h_2
           h_2 = h_tmp
           g_1 = g_2
           sec_istep = sec_istep + 1
        end do
        rt = x_tmp

 end subroutine

subroutine dispersion_relation(omega, rhs)
        
        implicit none

        complex(kind=8), intent(in) :: omega
        complex (kind=8), intent(out) :: rhs
        complex (kind=8) :: integral_e, integral_i, integral_z
        complex (kind=8) :: integral_te, integral_te_s_nu_a, integral_te_l_nu_a
        complex (kind=8) :: integral_te_l_nu, integral_te_s_nu,integral_te_l_nu2
        integer(kind=8) :: s

        if (t_e==0.) then
            integral_te = 0.
        else
            s = 4
            call vy_integral(omega,integral_te,s)
            !call theta_integral_simpson(omega,integral_te,s,t=0)
            print*, 'integral_te'
            print*, integral_te
!            s = 5
!            call theta_integral_simpson(omega,integral_te_l_nu,s,t=0)
!            !call vy_integral(omega,integral_te_l_nu,s)
!            s = 7
!            !call vy_integral(omega,integral_te_l_nu2,s)
!            call theta_integral_simpson(omega,integral_te_l_nu2,s,t=0)
!            s = 6
!            !call vy_integral(omega,integral_te_s_nu,s)
!            call theta_integral_simpson(omega,integral_te_s_nu,s,t=0)

!            integral_te_s_nu_a = f_t*(1.-(-theta*fprime_i*ky)/omega)
!            integral_te_l_nu_a = f_t*(-zi)*eps/3./nu0*&
!                (8.*sqrt(2./pi)*(omega-(-theta*fprime_i*ky))-12.*sqrt(2./pi)&
!                *(-theta*tprime_e*ky))
            
!            write (Dmixing_unit, '(15e12.4)') nu0,omega, &
!                integral_te,integral_te_s_nu,integral_te_s_nu_a,&
!                integral_te_l_nu, integral_te_l_nu2,integral_te_l_nu_a
        end if
        if (na_e==0.) then
            integral_e = 0.
        else
            s = 1
            call vy_integral(omega,integral_e,s)
            print*, 'integral_e'
            print*, integral_e
        end if
        if (na_i==0.) then
            integral_i = 0.
        else
            s = 2
            call vy_integral(omega,integral_i,s)
            print*, 'integral_i'
            print*, integral_i
        end if
        if (na_z==0.) then
            integral_z = 0.
        else
            s = 3
            call vy_integral(omega,integral_z,s)
        end if

        rhs = 1.0 + Zeff/Ti - na_e*integral_e - &
              na_i/Ti*(Z-Zeff)/(Z-1.0)*integral_i - &
              na_z/Ti*(Zeff-1.0)*Z/(Z-1.0)*integral_z - &
              t_e*integral_te + kperp**2*lambda_debye2

end subroutine

subroutine analytic_dr(omega, rhs)
        
        implicit none

        complex(kind=8), intent(in) :: omega
        complex (kind=8), intent(out) :: rhs
        complex (kind=8) :: term1, term2, term4, term5
        complex (kind=8) :: integral_i, integral_e, integral_z
        integer(kind=8) :: s,t 

        s = 2   
        t = 1 
        call vy_integral_simpson(omega,vy_1,vy_2,term1,s,t)
        t = 2
        call vy_integral_simpson(omega,vy_1,vy_2,term2,s,t)
        t = 4 
        call vy_integral_simpson(omega,vy_1,vy_2,term4,s,t)
        t = 5 
        call vy_integral_simpson(omega,vy_1,vy_2,term5,s,t)

        integral_i = term1 + term2/omega + term4/omega**2 + term5/omega**3


        s = 3   
        t = 1 
        call vy_integral_simpson(omega,vy_1,vy_2,term1,s,t)
        t = 2
        call vy_integral_simpson(omega,vy_1,vy_2,term2,s,t)
        t = 4 
        call vy_integral_simpson(omega,vy_1,vy_2,term4,s,t)
        t = 5 
        call vy_integral_simpson(omega,vy_1,vy_2,term5,s,t)

        integral_z =  term1 + term2/omega + term4/omega**2 + term5/omega**3

        rhs = 1. + Zeff/Ti - na_i/Ti*(Z-Zeff)/(Z-1.)*integral_i - &
              na_z/Ti*(Zeff-1.)*Z/(Z-1.)*integral_z 

end subroutine

subroutine vy_integral(omega,integral,s)

        implicit none

        integer(kind=8), intent(in) :: s
        complex(kind=8), intent(in) :: omega
        complex (kind=8), intent(out) :: integral
        complex (kind=8) :: integral1, integral2, integral3, integral4
        real (kind=8) :: discrim1, discrim2
        real (kind=8) :: vy_d, vy_l, vy_r, vy_b
        
        vy_b = cut_buffer
        integral1 = 0.
        integral2 = 0.
        integral3 = 0.
        integral4 = 0.

        if (.not.omd==0.) then
          if (s==2) then
           discrim1 = real(kz**2-4.*omd*ky*(omd*ky/2.*vy_1**2-omega))
           discrim2 = real(kz**2-4.*omd*ky*(omd*ky/2.*vy_2**2-omega))
           vy_d = sqrt(real(omega)*2./omd/ky+kz**2/omd**2/ky**2/2.)
          else if (s==1) then
           !discrim1 = 1.
           !discrim2 = 1.
           !vy_d = 0.
           discrim1 = real((kz/sqrt(Ti*mu_e))**2-&
             4.*omd*ky/(-Ti)*(-omd*ky/(-Ti)/2.*vy_1**2-omega))
           discrim2 = real((kz/sqrt(Ti*mu_e))**2-&
             4.*omd*ky/(-Ti)*(-omd*ky/(-Ti)/2.*vy_2**2-omega))
           vy_d = sqrt(2./omd/ky*(-Ti)*(real(omega)+&
             (kz/sqrt(Ti*mu_e))**2/omd/ky*(-Ti)/4.))
          else if (s==3) then
           discrim1 = 1.
           discrim2 = 1.
           vy_d = 0.
           !discrim1 = real((sqrt(1./mu_z)*kz)**2-&
           !  (1./Z)*4.*kv/A*(-1./Z*kv/2./A*vy_1**2-omega))
           !discrim2 = real((sqrt(1./mu_z)*kz)**2-&
           !  (1./Z)*4.*kv/A*(-1./Z*kv/2./A*vy_2**2-omega))
           !vy_d = sqrt(2.*A/kv*Z*(real(omega)+&
           !  (sqrt(1./mu_z)*kz)**2*A/kv*Z/4.))
          else if (s==4.or.s==5.or.s==6.or.s==7) then
           discrim1 = 1.
           discrim2 = 1.
           vy_d = 0.
          end if
          if ((.not.aimag(omega)>0.).and.discrim1*discrim2<0.) then
           vy_l = vy_d - vy_b
           vy_r = vy_d + vy_b
           if (vy_l>vy_1 .and. vy_r<vy_2) then
              call vy_integral_simpson(omega,vy_1,vy_l,integral1,s,t=0)
              call vy_integral_midpoint(omega,vy_l,vy_d,integral2,s)
              call vy_integral_midpoint(omega,vy_d,vy_r,integral3,s)
              call vy_integral_simpson(omega,vy_r,vy_2,integral4,s,t=0)
              integral = integral1+integral2+integral3+integral4 
           else if (vy_l<vy_1 .and. vy_r<vy_2) then
              integral1 = 0.
              call vy_integral_midpoint(omega,vy_1,vy_d,integral2,s)
              call vy_integral_midpoint(omega,vy_d,vy_r,integral3,s)
              call vy_integral_simpson(omega,vy_r,vy_2,integral4,s,t=0)
              integral = integral1+integral2+integral3+integral4 
           else if (vy_l>vy_1 .and. vy_r>vy_2) then
              call vy_integral_simpson(omega,vy_1,vy_l,integral1,s,t=0)
              call vy_integral_midpoint(omega,vy_l,vy_d,integral2,s)
              call vy_integral_midpoint(omega,vy_d,vy_2,integral3,s)
              integral4 = 0.
              integral = integral1+integral2+integral3+integral4 
           end if
          else !s==4
           call vy_integral_simpson(omega,vy_1,vy_2,integral1,s,t=0) 
           integral2 = 0.
           integral3 = 0.
           integral4 = 0.
           integral = integral1+integral2+integral3+integral4 
           !print*,'vy_integral, omega, s'
           !print*, integral, omega, s
          end if
        end if

        !print*, 'integral'
        !print*, omega,integral,s

end subroutine

subroutine vy_integral_midpoint(omega,vy_a,vy_b,integral,s)

        implicit none

        integer(kind=8), intent(in) :: s
        integer(kind=8) :: t
        real(kind=8),intent(in) :: vy_a, vy_b
        complex(kind=8), intent(in) :: omega
        complex (kind=8), intent(out) :: integral
        real(kind=8) :: vy_m, dvy, vy_s, vy_f
        complex(kind=8) :: f_m

        integral = 0.
        dvy = 0.01*cut_buffer
        vy_s = vy_a
        vy_f = vy_s+dvy
        do while (vy_s<vy_b)
                vy_m = (vy_s+vy_f)/2.
                call vz_integral(vy_m,f_m,s,omega,t=0)
                integral = integral + f_m*dvy
                vy_s = vy_f
                vy_f = vy_f + dvy
                if (vy_f>vy_b) then
                    vy_f = vy_b
                    dvy = vy_f-vy_s
                end if
        end do

        !print*, 'vy_integral_midpoint'
        !print*, omega,integral,s,t

end subroutine

subroutine theta_integral_simpson(omega,integral,s,t)

        implicit none

        integer(kind=8), intent(in) :: s,t
        complex(kind=8), intent(in) :: omega
        complex (kind=8), intent(out) :: integral
        real(kind=8) :: x_left, x_center, x_right, dx
        real(kind=8) :: theta_a, theta_b
        complex (kind=8) :: f_left, f_center, f_right
        complex (kind=8) :: i_trapezoid, i_simpson
        integer :: theta_istep

        theta_istep = 0.
        theta_a = theta_c
        theta_b = pi - theta_c
        dx = 0.01
        x_left = theta_a
        x_right = x_left + dx
        theta_v = x_left
        call vy_integral_simpson(omega,vy_1,vy_2,f_left,s,t)
        integral = 0.0
        do while (x_left<theta_b.and.theta_istep<int_nsteps)
                x_center = 0.5*(x_left+x_right)
                theta_v = x_center
                call vy_integral_simpson(omega,vy_1,vy_2,f_center,s,t)
                theta_v = x_right
                call vy_integral_simpson(omega,vy_1,vy_2,f_right,s,t)
                i_trapezoid = 0.5*(f_left+2.0*f_center+f_right)*0.5*dx
                i_simpson = dx*(f_left+4.0*f_center+f_right)/6.
                do while (abs(i_trapezoid-i_simpson)>int_tol)
                        dx=fr*dx*(abs(i_trapezoid-i_simpson)/int_tol)**(-1.0/3.0)
                        if (dx>0.2) dx=0.2
                        x_right = x_left + dx
                        if (x_right>theta_b) then
                                x_right = theta_b
                                dx = x_right - x_left
                        end if
                        x_center = 0.5*(x_left+x_right)
                        theta_v = x_center
                        call vy_integral_simpson(omega,vy_1,vy_2,f_center,s,t)
                        theta_v = x_right
                        call vy_integral_simpson(omega,vy_1,vy_2,f_right,s,t)
                        i_trapezoid = 0.5*(f_left+2.0*f_center+f_right)*0.5*dx
                        i_simpson = dx*(f_left+4.0*f_center+f_right)/6.
                end do
                integral = integral + i_simpson
                x_left = x_right
                dx = fr*dx*(abs(i_trapezoid-i_simpson)/int_tol)**(-1.0/3.0)
                if (dx>0.2) dx=0.2
                x_right = x_left + dx
                if (x_right>theta_b) then
                        x_right = theta_b
                        dx = x_right - x_left
                end if
                f_left = f_right
                theta_istep = theta_istep+1
        end do
        if (theta_istep==int_nsteps) then
           int_istep = int_nsteps
           integral = 0.
        end if

  end subroutine

subroutine vy_integral_simpson(omega,vy_a,vy_b,integral,s,t)

        implicit none

        integer(kind=8), intent(in) :: s,t
        real(kind=8),intent(in) :: vy_a, vy_b
        complex(kind=8), intent(in) :: omega
        complex (kind=8), intent(out) :: integral
        real(kind=8) :: x_left, x_center, x_right, dx
        complex (kind=8) :: f_left, f_center, f_right
        complex (kind=8) :: i_trapezoid, i_simpson
        integer :: vy_istep

        vy_istep = 0.
        dx = 0.1
        x_left = vy_a
        x_right = x_left + dx
        call vz_integral(x_left,f_left,s,omega,t)
        integral = 0.0

        do while (x_left<vy_b.and.vy_istep<int_nsteps)
                x_center = 0.5*(x_left+x_right)
                call vz_integral(x_center,f_center,s,omega,t)
                call vz_integral(x_right,f_right,s,omega,t)
                i_trapezoid = 0.5*(f_left+2.0*f_center+f_right)*0.5*dx
                i_simpson = dx*(f_left+4.0*f_center+f_right)/6.
                do while (abs(i_trapezoid-i_simpson)>int_tol)
                        dx=fr*dx*(abs(i_trapezoid-i_simpson)/int_tol)**(-1.0/3.0)
                        if (dx>0.2) dx=0.2
                        x_right = x_left + dx
                        if (x_right>vy_b) then
                                x_right = vy_b
                                dx = x_right - x_left
                        end if
                        x_center = 0.5*(x_left+x_right)
                        call vz_integral(x_center,f_center,s,omega,t)
                        call vz_integral(x_right,f_right,s,omega,t)
                        i_trapezoid = 0.5*(f_left+2.0*f_center+f_right)*0.5*dx
                        i_simpson = dx*(f_left+4.0*f_center+f_right)/6.
                end do
                integral = integral + i_simpson
                x_left = x_right
                dx = fr*dx*(abs(i_trapezoid-i_simpson)/int_tol)**(-1.0/3.0)
                if (dx>0.2) dx=0.2
                x_right = x_left + dx
                if (x_right>vy_b) then
                        x_right = vy_b
                        dx = x_right - x_left
                end if
                f_left = f_right
                vy_istep = vy_istep+1
        end do
        if (vy_istep==int_nsteps) then
           int_istep = int_nsteps
           integral = 0.
        end if

        !print*, 'vy_integral_simpson'
        !print*, omega,integral

  end subroutine

  subroutine vz_integral(vy,integral,s,omega,t)

        implicit none

        integer(kind=8), intent(in) :: s,t
        real(kind=8), intent(in) :: vy
        complex(kind=8), intent(in) :: omega
        real(kind=8) :: vz_a, vz_b
        complex (kind=8), intent(out) :: integral
        real(kind=8) :: x_left, x_center, x_right, dx
        complex (kind=8) :: f_left, f_center, f_right
        complex (kind=8) :: i_trapezoid, i_simpson
        complex (kind=8) :: res_v1, res_v2
        complex (kind=8) :: shift_vz
        real(kind=8) :: rho,bj,dj,fj
        integer(kind=8) :: n
        integer(kind=8) :: vz_istep
        real(kind=8) :: vz_te, vy_te

        if (t==0.) then
           if (s==2) then
              call vz_integral_residue(vy,omega,s,res_v1,res_v2,&
                shift_vz)
           else if (s==1) then
              call vz_integral_residue(vy,omega,s,res_v1,res_v2,&
                shift_vz)
           else if (s==3) then
              call vz_integral_residue(vy,omega,s,res_v1,res_v2,&
                shift_vz)
           else if (s==4.or.s==5.or.s==6.or.s==6) then
              res_v1 = 0.
              res_v2 = 0.
              shift_vz = 0.
           end if
        else
           res_v1 = 0.
           res_v2 = 0.
           shift_vz = 0.
        end if

        if ((.not.s==4).and.(.not.s==5).and.(.not.s==6).and.(.not.s==7)) then
           vz_istep = 0
           vz_a = vz_1
           vz_b = vz_2
           dx = 0.1
           x_left = vz_a
           x_right = x_left + dx
           f_left = integrand(x_left,vy,omega,shift_vz,s,t)
           integral = 0.0

           do while (x_left<vz_b.and.vz_istep<int_nsteps)
                x_center = 0.5*(x_left+x_right)
                f_center = integrand(x_center,vy,omega,shift_vz,s,t)
                f_right = integrand(x_right,vy,omega,shift_vz,s,t)
                i_trapezoid = 0.5*(f_left+2.0*f_center+f_right)*0.5*dx
                i_simpson = dx*(f_left+4.0*f_center+f_right)/6.
                do while (abs(i_trapezoid-i_simpson)>int_tol)
                        dx=fr*dx*(abs(i_trapezoid-i_simpson)/int_tol)**(-1.0/3.0)
                        if (dx>0.2) dx=0.2
                        x_right = x_left + dx
                        if (x_right>vz_b) then
                                x_right = vz_b
                                dx = x_right - x_left
                        end if
                        x_center = 0.5*(x_left+x_right)
                        f_center = integrand(x_center,vy,omega,shift_vz,s,t)
                        f_right = integrand(x_right,vy,omega,shift_vz,s,t)
                        i_trapezoid = 0.5*(f_left+2.0*f_center+f_right)*0.5*dx
                        i_simpson = dx*(f_left+4.0*f_center+f_right)/6.
                end do
                integral = integral + i_simpson
                x_left = x_right
                dx = fr*dx*(abs(i_trapezoid-i_simpson)/int_tol)**(-1.0/3.0)
                if (dx>0.2) dx=0.2
                x_right = x_left + dx
                if (x_right>vz_b) then 
                        x_right = vz_b
                        dx = x_right - x_left
                end if
                f_left = f_right
                vz_istep = vz_istep+1
           end do
        else
           ! trapped electrons vz=v*cos(theta),vy=v*sin(theta)
           !vz_te = vy*cos(theta_v)
           !vy_te = vy*sin(theta_v)
           vy_te = vy
           vz_te = 0.
           integral = integrand(vz_te,vy_te,omega,shift_vz,s,t)
        end if

        if (vz_istep==int_nsteps) then
           int_istep = int_nsteps
           integral = 0.
        else
           n=0
           if (s==2) then
              !rho = sqrt(gyy*ky**2+kx**2)*vy
              rho = kperp*vy
              call bjndd (n, rho, bj, dj, fj)
              if (isnan(bj)) bj = 1.
              integral = &
                 bj**2*sqrt(1.0/pi/2.0)*vy*exp(-vy**2/2.)*(integral+res_v1+res_v2)
           else if (s==1) then
              !rho = sqrt(theta*mu_e)*sqrt(gyy*ky**2+kx**2)*vy
              rho = sqrt(mu_e/Ti)*kperp*vy
              call bjndd (n, rho, bj, dj, fj)
              if (isnan(bj)) bj = 1.
              integral = &
                 bj**2*sqrt(1.0/pi/2.0)*vy*exp(-vy**2/2.)*(integral+res_v1+res_v2)
           else if (s==3) then
              !rho = sqrt(mu_z)/Z*sqrt(gyy*ky**2+kx**2)*vy
              rho = sqrt(mu_z)/Z*kperp*vy
              call bjndd (n, rho, bj, dj, fj)
              if (isnan(bj)) bj = 1.
              integral = &
                 bj**2*sqrt(1.0/pi/2.0)*vy*exp(-vy**2/2.)*(integral+res_v1+res_v2)
           else if (s==4.or.s==5.or.s==6.or.s==7) then
              !rho = sqrt(theta*mu_e)*sqrt(ky**2+kx**2)*vy
              !call bjndd (n, rho, bj, dj, fj)
              !if (isnan(bj)) bj = 1.
              integral = &
                 f_t*sqrt(2.0/pi)*vy**2*exp(-vy**2/2.)*integral
                 !sin(theta_v)*sqrt(1./pi/2.0)*vy**2*exp(-vy**2/2.)*integral
           end if
        end if

  end subroutine

  subroutine vz_integral_residue(vy,omega,s,res_v1,res_v2,shift_vz)

        implicit none

        real(kind=8), intent(in) :: vy
        complex(kind=8), intent(in) :: omega
        integer(kind=8), intent(in) :: s
        complex (kind=8), intent(out) :: res_v1, res_v2, shift_vz
        real(kind=8) :: real_discrim
        real(kind=8) :: v0,coeff_a,coeff_b,coeff_d
        complex(kind=8) :: dv,v1,v2
        complex(kind=8) :: coeff_c, coeff_e, coeff_f, coeff_g, tmp_f, tmp_g
        complex(kind=8) :: tmp_res1, tmp_res2

        if (s==1.or.s==4) then
           tprime = tprime_e
        else
           tprime = tprime_i
        end if
        if (s==2) then
          if (omd==0.) then
           if (kz==0.) then
             v1 = 100.
             v2 = 100.
             tmp_res1 = 0.
             tmp_res2 = 0.
           else
             v1 = omega/kz
             v2 = 0.
             tmp_res1 = 2.*pi*zi*exp(-v1**2/2.)*&
                        (-1./kz*(omega-fprime*ky-&
                        (-3./2.+v1**2/2.+vy**2/2.)*&
                        tprime*ky))
             tmp_res2 = 0.
           end if
          else
           coeff_a = omd*ky
           coeff_b = kz
           coeff_c = omd*ky/2.*vy**2-omega
           coeff_d = 0.5*tprime*ky
           coeff_e = (-3./2.+vy**2/2.)*tprime*ky+fprime*ky-omega
           v0 = -coeff_b/2./coeff_a
           dv = sqrt(coeff_b**2-4.*coeff_a*coeff_c)/coeff_a/2.
           real_discrim = real(coeff_b**2-4.*coeff_a*coeff_c)
           if (real_discrim<0..and.aimag(omega)<0.) then
              v1 = v0 - dv
              v2 = v0 + dv
           else
              v1 = v0 + dv
              v2 = v0 - dv
           end if
           tmp_res1 = 2*pi*zi*exp(-v1**2/2.)*(coeff_d*&
                      v1**2+coeff_e)/coeff_a/(v1-v2)
           tmp_res2 = 2*pi*zi*exp(-v2**2/2.)*(coeff_d*&
                      v2**2+coeff_e)/coeff_a/(v2-v1)
          end if
        else if (s==1) then
          if (omd==0.) then
           if (kz==0.) then
             v1 = 100.
             v2 = 100.
             tmp_res1 = 0.
             tmp_res2 = 0.
           else
             v1 = omega/kz/sqrt(1./Ti/mu_e)
             v2 = 0.
             tmp_res1 = 2.*pi*zi*exp(-v1**2/2.)*&
                        (-1./kz/sqrt(1./Ti/mu_e)*&
                        (omega+1./Ti*fprime*ky+&
                        (-3./2.+v1**2/2.+vy**2/2.)*&
                        1./Ti*tprime*ky))
             tmp_res2 = 0.
           end if
          else
           coeff_a = -1./Ti*omd*ky
           coeff_b = sqrt(1./Ti/mu_e)*kz
           coeff_c = -1./Ti*omd/2.*ky*vy**2-omega
           coeff_d = -1./Ti*0.5*tprime*ky
           coeff_e = -1./Ti*((-3./2.+vy**2/2.)*tprime*ky+fprime*ky)-omega
           v0 = -coeff_b/2./coeff_a
           dv = sqrt(coeff_b**2-4.*coeff_a*coeff_c)/coeff_a
           real_discrim = real(coeff_b**2-4.*coeff_a*coeff_c)
           if (real_discrim<0..and.aimag(omega)<0.) then
              v1 = v0-0.5*dv
              v2 = v0+0.5*dv
           else
              v1 = v0+0.5*dv
              v2 = v0-0.5*dv
           end if
           tmp_res1 = 2*pi*zi*exp(-v1**2/2.)*(coeff_d*&
                      v1**2+coeff_e)/coeff_a/(v1-v2)
           tmp_res2 = 2*pi*zi*exp(-v2**2/2.)*(coeff_d*&
                      v2**2+coeff_e)/coeff_a/(v2-v1)
          end if
        else if (s==3) then
          if (omd==0.) then
           if (kz==0.) then
             v1 = 100.
             v2 = 100.
             tmp_res1 = 0.
             tmp_res2 = 0.
           else
             v1 = omega/kz/sqrt(1./mu_z)
             v2 = 0.
             tmp_res1 = 2.*pi*zi*exp(-v1**2/2.)*&
                        (-1./kz/sqrt(1./mu_z)*&
                        (omega-1./Z*fprime*ky-&
                        (-3./2.+v1**2/2.+vy**2/2.)*&
                        1./Z*tprime*ky))
             tmp_res2 = 0.
           end if
          else
           coeff_a = 1./Z*omd*ky
           coeff_b = sqrt(1./mu_z)*kz
           coeff_c = 1./Z*omd/2.*ky*vy**2-omega
           coeff_d = 1./Z*0.5*tprime*ky
           coeff_e = 1./Z*((-3./2.+vy**2/2.)*tprime*ky+fprime*ky)-omega
           v0 = -coeff_b/2./coeff_a
           dv = sqrt(coeff_b**2-4.*coeff_a*coeff_c)/coeff_a
           real_discrim = real(coeff_b**2-4.*coeff_a*coeff_c)
           if (real_discrim<0..and.aimag(omega)<0.) then
              v1 = v0-0.5*dv
              v2 = v0+0.5*dv
           else
              v1 = v0+0.5*dv
              v2 = v0-0.5*dv
           end if
           tmp_res1 = 2*pi*zi*exp(-v1**2/2.)*(coeff_d*&
                      v1**2+coeff_e)/coeff_a/(v1-v2)
           tmp_res2 = 2*pi*zi*exp(-v2**2/2.)*(coeff_d*&
                      v2**2+coeff_e)/coeff_a/(v2-v1)
          end if
        end if

        !print*, 'roots v1,v2'
        !print*,v1,v2
        if (aimag(v1) > 0.) then
           shift_vz = 0.
           res_v1 = 0.
           res_v2 = 0.
        else if (aimag(v1) > -vi) then
           shift_vz = -zi*2.*vi
           res_v1 = 0.
           res_v2 = -tmp_res2
        else
           shift_vz = 0.
           res_v1 = tmp_res1
           res_v2 = -tmp_res2
        end if
        !print*, 'residues'
        !print*,res_v1,res_v2

  end subroutine

  function integrand(vz_prime,vy,omega,shift_vz,s,t)

        real(kind=8) :: vz_prime,vy
        complex(kind=8) :: shift_vz,vz
        complex(kind=8) :: omega,integrand,int_tmp
        complex(kind=8) :: omega_star_n, omega_star_t, omega_d
        complex(kind=8) :: omega_parallel, coll_te
        integer(kind=8) :: s,t

        if (s==1.or.s==4.or.s==5.or.s==6.or.s==7) then
           tprime = tprime_e
        else
           tprime = tprime_i
        end if

        vz = vz_prime + shift_vz

        if (s==2) then
                omega_star_n = fprime*ky
                omega_star_t = 0.5*(-3.0+vy**2+&
                        vz**2)*tprime*ky
                omega_d = omd*ky*(0.5*vy**2+vz**2)
                omega_parallel = kz*vz
        else if (s==1) then
                omega_star_n = -1./Ti*fprime*ky
                omega_star_t = -1./Ti*0.5*(-3.0+vy**2+&
                        vz**2)*tprime*ky
                omega_d = -1./Ti*omd*ky*(0.5*vy**2+vz**2)
                omega_parallel = sqrt(1./Ti/mu_e)*kz*vz
        else if (s==3) then
                omega_star_n = 1.0/Z*fprime*ky
                omega_star_t = 1.0/Z*0.5*(-3.0+vy**2+&
                        vz**2)*tprime*ky
                omega_d = 1.0/Z*omd*ky*(0.5*vy**2+vz**2)
                omega_parallel = sqrt(1.0/mu_z)*kz*vz
        else if (s==4.or.s==5.or.s==6.or.s==7) then
                !omega_star_n = -fprime*ky/Ti
                !omega_star_t = -(-3./2.+vy**2/2.+&
                !        vz**2/2.)*tprime*ky/Ti
                !omega_d = -omde/Ti*omd*ky*(vy**2/2.+vz**2)*0.53
                !coll_te = zi*2.**(3./2.)*nu0/eps/(sqrt(vy**2+vz**2))**3
                omega_star_n = -fprime*ky/Ti
                omega_star_t = -(-3./2.+vy**2/2.)&
                               *tprime*ky/Ti
                omega_d = -omde/Ti*omd*ky*vy**2*omde_avg
                coll_te = coll_c*zi*2.**(3./2.)*nu0*&
                          (1.+vy/(2.+vy))&
                          /eps/vy**3
                
        end if
        
        if (t==0.and.(.not.s==4).and.(.not.s==5).and.(.not.s==6).and.(.not.s==7)) then
                int_tmp = exp(-vz**2/2.)*&
                  (omega-omega_star_n-omega_star_t)/&
                  (omega-omega_parallel-omega_d)
        else if (t==0.and.s==4) then
                !int_tmp = exp(-vz**2/2.)*&
                !  (omega-omega_star_n-omega_star_t)/&
                int_tmp = (omega-omega_star_n-omega_star_t)/&
                  (omega-omega_d+coll_te)
        else if (t==0.and.s==5) then
                int_tmp = (omega-omega_star_n-omega_star_t)/&
                  coll_te*(1.-omega/coll_te)
        else if (t==0.and.s==7) then
                int_tmp = (omega-omega_star_n-omega_star_t)/&
                  coll_te
        else if (t==0.and.s==6) then
                int_tmp = (omega-omega_star_n-omega_star_t)/&
                  omega
        else if (t==1) then
                int_tmp = exp(-vz**2/2.)
        else if (t==2) then
                int_tmp = exp(-vz**2/2.)*&
                   (omega_d-omega_star_n-omega_star_t)
        else if (t==3) then
                int_tmp = exp(-vz**2/2.)*&
                   (-omega_d*(omega_star_n+omega_star_t))
        else if (t==4) then
                int_tmp = exp(-vz**2/2.)*&
                   (omega_parallel**2-omega_d*(omega_star_n+omega_star_t))
        else if (t==5) then
                int_tmp = exp(-vz**2/2.)*&
                   ((-omega_star_n-omega_star_t)*omega_parallel**2)
        end if

        integrand = int_tmp

  end function

  subroutine read_input_file

        use file_utils, only: input_unit_exist, input_unit

        implicit none

        integer(kind=8) :: in_file
        logical :: exist

        namelist / parameters / sec_nsteps,int_nsteps,nkz,dkz,kz0,kz1, &
                                nky,dky,ky0,ky1,dfprime,nfprime, &
                                fprime0,dtprime_i,ntprime_i, &
                                tprime0_i,dtprime_e,ntprime_e, &
                                tprime0_e, nomd,domd,omd0, &
                                vy_1,vy_2,vz_1,vz_2,om1_re,om1_im, &
                                eps,nu0,coll_c,vi,int_tol,sec_tol,fr, &
                                Ti,Zeff,Z,mu_e,mu_z,na_e,na_z,na_i, &
                                nkperp, dkperp, kperp0, kperp1, cut_buffer, &
                                seed_fr, int_nsteps, lower_bnd, t_e, omde, &
                                kperp, fr_te, lambda_debye2
        
        sec_nsteps = 10
        int_nsteps = 1000

        fr_te = 1.
        lambda_debye2 = 0.

        nkz = 1
        dkz = 0.1
        kz0 = 0.0
        kz1 = 0.0

        nky = 1
        dky = 0.1
        ky0 = 0.0
        ky1 = 0.0
        kperp = 1.

        dfprime = 0.5
        nfprime = 1.0
        fprime0 = 0.0

        ntprime_i = 1
        dtprime_i = 1.0
        tprime0_i = 0.0

        ntprime_e = 1
        dtprime_e = 1.0
        tprime0_e = 0.0

        nomd = 1
        domd = 1.0
        omd0 = 0.0

        vy_1 = 0.0
        vy_2 = 0.0
        vz_1 = 8.0
        vz_2 = 8.0

        om1_re = 1.0
        om1_im = 0.1

        nu0 = 0.05
        coll_c = 1.
        eps = 0.2
        omde = 1.
        vi  = 1.0
        sec_tol = 1.0E-05
        int_tol = 1.0E-05
        fr = 0.5

        Ti = 1.0
        Zeff = 1.65
        Z = 5
        mu_e = 5.45D-04
        mu_z = 10.8
        na_e = 0.0
        na_z = 1.0
        na_i = 1.0
        t_e = 1.0
        
        nkperp = 1
        dkperp = 0.5
        kperp0 = 0.0
        kperp1 = 10.0
        
        cut_buffer = 0.001
        seed_fr = 0.95
        lower_bnd = 1E-16
        
    in_file = input_unit_exist ("parameters", exist)
    if (exist) read (unit=input_unit("parameters"), nml=parameters)

    if (t_e==0.) then
        f_t = 0.
        theta_c = pi/2.
        omde_avg = 1.
    else
        f_t = sqrt(2.*eps/(1.+eps))*fr_te
        theta_c = acos(f_t)
        omde_avg = (3./4.*(pi-2.*theta_c)+&
                   sin(2.*(pi-theta_c))/8.-&
                   sin(2.*theta_c)/8.)/(pi-2.*theta_c)
    end if

    if (t_e ==1.) then
        print*,'fraction of trapped electrons'
        print*, f_t
        print*,'critical theta'
        print*, theta_c
        print*,'average sin^2(theta)/2.+cos^2(theta)'
        print*, omde_avg
    end if

    root_tmp_fl = om1_re + om1_im*zi
  end subroutine read_input_file

subroutine init_grids

        implicit none

        allocate(ky_grid(nky), kz_grid(nkz), kperp_grid(nkperp))

        do ikz = 1, nkz
                kz_grid(ikz) = ikz*dkz+kz0
        end do
        do iky = 1, nky
                ky_grid(iky) = iky*dky+ky0
        end do
        do ikperp = 1, nkperp
                kperp_grid(ikperp) = ikperp*dkperp+kperp0
        end do

        allocate(fprime_grid(nfprime), tprime_i_grid(ntprime_i))
        allocate(tprime_e_grid(ntprime_e))
        allocate(omd_grid(nomd))

        do ifprime = 1, nfprime
                fprime_grid(ifprime) = ifprime*dfprime+fprime0
        end do 
        do itprime_i = 1, ntprime_i
                tprime_i_grid(itprime_i) = itprime_i*dtprime_i+tprime0_i
        end do 

        do itprime_e = 1, ntprime_e
                tprime_e_grid(itprime_e) = itprime_e*dtprime_e+tprime0_e
        end do 
        do iomd = 1, nomd
                omd_grid(iomd) = iomd*domd+omd0
        end do

        allocate(ikperp_ref(1), ikz_ref(1), iky_ref(1))


end subroutine

subroutine init_output_files

        call open_output_file(gamma_unit,'.datgam')
        write (gamma_unit,'(9a12)') "tprime_i","fprime_i","tprime_e",&
                "omd","kperp","ky","kz","omega","gamma"
!        !call open_output_file(Dmixing_unit,'.datDmixing')
!        !write (Dmixing_unit,'(9a12)') "tprime_i","fprime_i","tprime_e","fprime_e",&
!        !        "omd","kx","ky","kz","Dmixing"
!        call open_output_file(Dmixing_unit,'.te')
!        write (Dmixing_unit,'(15a12)') "nu0","omega","","full_int","",&
!                "s_nu_int","","s_nu_int_a","","l_nu_int","","l_nu_int2","","l_nu_int_a",""
        call open_output_file(out_unit,'.dat')
        !write (out_unit,'(10a12)') "tprime_i","fprime_i","tprime_e","fprime_e",&
        !        "omd","kperp","ky","kz","omega","gamma"
        write (out_unit,'(9a12)') "tprime_i","fprime","tprime_e",&
                "ky","kz","omega","gamma","kperp","omd"!!,"omlina","gmlina",&
                !!"omquad","gmquad","omquada","gmquada" 
end subroutine

end program
