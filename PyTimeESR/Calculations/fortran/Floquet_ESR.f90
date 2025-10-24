program Floquet_ESR

! Long-time solution for CW ESR
! based on the theory by Galvez et al
! Phys. Rev. B Physical Review B 104 (24), 245435 (2021)

!
! gnu licence 3.0 (c) J. Reina Galvez & N. Lorente
!
!  Version 1.0 Including Lamb shift with cutoff
!  January 2023

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!This is the main file and calls in:
Use declarations !all variables, constants etc are defined here
Use OpenFiles !all units of files in open statements are here
!Use feed_back_on
Use io !(input/output) the reading and writing routines are here
Use timer
Use H_QD !(Hamiltonian Quantum Dot) contains the Hamiltonian and solution
Use QME_F !(Quantum Master Equation) contains rates and The Matrix
Use Transport_F ! computation of electron current
Use SpinTools
Use Matrix_Coeff
use output
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                          input run's values                            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     call clock ('STEP 1:: Reading INPUT and Hamiltonian for the QME', 1)

     call reading_input ( gamma_R_0, A_R, gamma_L_0, A_L, phi, orb,&
       NF,NCF,GammaC,Cutoff,redimension,Nd, frequency,&
       bias_R, bias_L, Spin_polarization_R, Spin_polarization_L, Temperature, &
       Electrode, transport_exponentR, transport_exponentL, B_R, B_L, p_maxR, p_maxL, write_populations, &
       write_coherences, spinflo,Ef,FermiP, VDC)

!      call clock ('Finished reading INPUT for the QME ', 2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                    solve the spin+orbital system                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      call clock ('STEP 2:: Diagonalizing Hamiltonian ', 1)

    call Hamiltonian (N,lambda,Delta,Eigenvalues,H,Ss,spinX,spinY,spinZ,&
    &spin2_T,runs,Nm,hx,hy,hz,Ndim,&
    &redimension,orb,bias_R, bias_L,Cutoff)

     call clock ('STEP 2:: Finished diagonalizing Hamiltonian and reading', 2)

    Nmatrix=Ndim*Ndim*NF ! total dimension matrix A

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                      ALLOCATES                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     allocate (G (Ndim, Ndim, Ndim, Ndim, NF, 2))
     allocate (GC (Ndim, Ndim, Ndim, Ndim, NF, 2))
     allocate (GA (Ndim, Ndim, Ndim, Ndim, NF), GCA (Ndim, Ndim, Ndim, Ndim, NF))
     allocate (A (Nmatrix, Nmatrix),A_old(Nmatrix,Nmatrix))
     allocate (B (Nmatrix))
     allocate (Rho (Ndim, Ndim, NF))
     allocate (rho2(Ndim,Ndim))
     allocate (WORK (Nmatrix))
     allocate (RWORK (Nmatrix))
     allocate (SWORK (Nmatrix*(Nmatrix+1)))
     allocate (IPIV (Nmatrix))
     allocate (curr(NF))

! loop on driving freq
!uencies. Parallelized with coarrays

! we open the outfile here too
        k=0


! add more floquet numbers (\pm 2) if the driving is large enough to make them important

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                    Solve QME in Floquet basis                          !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   call clock ('STEP 3:: Calculating rates off resonance freq and feedback loop', 1)
! the off resonance frequency should be too large

! First evaluate rates in the H_QD basis and Floquet indices
!      call clock ('STEP 3:: Computing rates ', 1)
  
  !     right electrode
  call rates (Ndim, orb, frequency, gamma_R_0, lambda, Spin_polarization_R,&
            NCF, p_maxR, A_R, phi, B_R, GammaC, bias_R, Delta, Cutoff, Temperature,& 
            N_int, transport_exponentR, GA, GCA)
  G(:,:,:,:,:,1)  = GA(:,:,:,:,:)
  GC(:,:,:,:,:,1) = GCA(:,:,:,:,:)
  
  !     left electrode
  call rates (Ndim, orb, frequency, gamma_L_0, lambda, Spin_polarization_L,&
    NCF, p_maxL, A_L, phi, B_L, GammaC, bias_L, Delta, Cutoff, Temperature,& 
    N_int, transport_exponentL, GA, GCA)
  G(:,:,:,:,:,2)  = GA(:,:,:,:,:)
  GC(:,:,:,:,:,2) = GCA(:,:,:,:,:)
  
  call coeff_matrix (Ndim, frequency, NF, NCF, G, Nmatrix, Rho)

! Compute the DC electron current
  call Current (Ndim, NF, NCF, Rho, GC(:,:,:,:,:,2-Electrode), curr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                    output the results                                  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Current
      
!      call clock ('STEP 6:: Writing output ', 1)
            
    write(*,*) ''
    write(*,*) '************************************************************************************************'
    write(*,*) 'Calculation done and output written!!!'
    write(*,*) ''
    write(*,*) 'Every output is divided by Floquet numbers 0,-1,1 (n indicates the minus sign)'
    write(*,*) ''
    
    call check_populations (Rho, NCF, Ndim, bias_R, bias_L, frequency)
    call write_rate_out (G, NCF, Ndim)
    call write_cur_out (curr, NCF)
    
    if (write_populations) then
        call write_pop_out (Rho, NCF, Ndim)
    endif

    if (write_coherences) then
        call write_coh_out (Rho, NCF)
    endif

    if (spinflo) then
        call write_spin_out (Rho, NCF, Nm, Ndim, Ss, spinX, spinY, spinZ, spin2_T, H, hx, hy, hz,&
                              Sx, Sy, Sz, Sh, spin2_ave)
    endif

    write(*,*) 'No zero real and imag left and right rates for Floquet numbers -2,-1,0,1,2 written in rates_floquet.dat'
    write(*,*) 'while rates_floquet0.dat contains only the Floquet zero'
    write(*,*) ''
    write(*,*) 'DC current wrote in Current_0.dat'

     call clock ('Final STEP:: Everything done', 2)

end program Floquet_ESR
