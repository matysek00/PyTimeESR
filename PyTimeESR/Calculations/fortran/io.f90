module io
Use declarations
Use OpenFiles
CONTAINS
   subroutine reading_input ( gamma_R_0, A_R, gamma_L_0, A_L,phi,orb, &
       NF,NCF,GammaC,Cutoff,redimension,Nd, frequency,&
       bias_R, bias_L, Spin_polarization_R, Spin_polarization_L, Temperature, &
       Electrode, transport_exponentR, transport_exponentL, B_R, B_L, p_maxR, p_maxL, &
       write_populations, write_coherences, spinflo, Ef, FermiP, VDC)
   implicit none
   integer :: Ninterval, Nfreq, Ntime, n, i, Electrode, FermiP, p_maxR, p_maxL
   real (q) :: Spin_polarization_R, Spin_polarization_L, Temperature
   real (q) :: gamma_R_0,gamma_L_0,phi,Cutoff, frequency
   real (q) :: bias_R, bias_L, step_freq, B_L, B_R
   real (q) :: transport_exponentR, transport_exponentL
   real(q), intent(out) :: Ef,VDC
   logical :: write_populations, write_coherences, redimension,presence2,spinflo
   integer :: Nd,NF,NCF,orb
   complex (qc) :: A_R,A_L,gammaC

   open (unit_input, file='Floquet.input', status='old')

       read (unit_input, *)  ! separator
       read (unit_input, *)  ! separator
       read (unit_input, *)  ! separator
       read (unit_input, *) NCF ! Maximum Floquet index, number of Floquets: 
       read (unit_input, *) frequency ! number of points in frequency swept
       read (unit_input, *)  ! separator
       read (unit_input, *) gamma_R_0 ! in meV : gamma_R_0= 2*pi*W_R_0*W_R_0*rho
       read (unit_input, *) gamma_L_0 ! in meV : gamma_L_0= 2*pi*W_R_0*W_L_0*rho
       read (unit_input, *) A_R ! real and imag  driving: gamma_R_1/gamma_R_0
       read (unit_input, *) A_L ! real and imag driving: gamma_L_1/gamma_L_0
       read (unit_input, *) phi ! phase shift in rad - azimuthal angle of the bloch sphere
       read (unit_input, *) Cutoff ! For convergence in Lambshift (meV)
       read (unit_input, *) gammaC ! Broadening of Green's function (meV), Real for SOMO, Imag for SUMO
       read (unit_input, *) N_int ! Number of points in I11 and I21 (see Manual)
       read (unit_input, *) orb ! set number of orbital, 1 will not use any Diego Lehman coeff
       read (unit_input, *)  ! separator
       read (unit_input, *) FermiP ! 0 Ef closer to muR if gammaR>>gammaL, 1 EF to muL if gammaR>>gammaL
       ! 2 Ef at muR if gammaR>>gammaL and or at muL if gammaR<<gammaL (simple Ef definition)
       read (unit_input, *) bias_R ! mV right electrode
       read (unit_input, *) bias_L ! mV left electrode
       read (unit_input, *) Temperature
       read (unit_input, *) Spin_polarization_R !  between -1 and 1
       read (unit_input, *) Spin_polarization_L !  between -1 and 1
       read (unit_input, *) Electrode !  0 is left and 1 is right
       read (unit_input, *) transport_exponentR !gamma(E) = gamma_0 e^(c*E/2) 
       read (unit_input, *) transport_exponentL !gamma(E) = gamma_0 e^(c*E/2) 
       read (unit_input, *)  ! separator
       read (unit_input, *) B_R 
       read (unit_input, *) B_L
       read (unit_input, *) p_maxR
       read (unit_input, *) p_maxL
       read (unit_input, *)  ! separator
       read (unit_input, *) write_populations ! .true. will write the populations only
       read (unit_input, *) write_coherences ! .true. will write the full density matrix
       read (unit_input, *) spinflo ! .true. will write the spin components
       read (unit_input, *)  ! separator
       read (unit_input, *) redimension !if .true. reduce states to open plus a few closed channels
       read (unit_input, *) Nd  !new dimension
       close (unit_input)

       NF = 2*NCF +1
       
       
       if (Cutoff<500*abs(gammaC))then
       print*,''
       print*,'WARNING'
       print*,'Low Cutoff given the broadening of the impurity Greens function'
       print*,'Increase Cutoff (Cutoff > 500*gammaC) to ensure convergence in the wide band limit'
       print*,''
       endif
       
! Then gamma_R_1/gamma_R_0 is exactly "the driving"
       print *, ' '
       print *, ' The driving is:'
       print *, '   Left electrode  =', 100*A_L,'   Right electrode =', 100*A_R, 'in %'
       print *, ' '
       print *, ' The applied bias is:'
       print *, '   Left electrode  =', bias_L,'   Right electrode =', bias_R, ' in mV'
       print *, ' '
       
       Ef=(bias_L*gamma_L_0/(gamma_L_0+gamma_R_0)+bias_R*gamma_R_0/(gamma_L_0+gamma_R_0))*(1-FermiP)+&
       &(bias_L*gamma_R_0/(gamma_L_0+gamma_R_0)+bias_R*gamma_L_0/(gamma_L_0+gamma_R_0))*(FermiP)

       if ((gamma_R_0.gt.gamma_L_0).and.(FermiP.eq.2)) then
              Ef= bias_R
       elseif ((gamma_R_0.gt.gamma_L_0).and.(FermiP.eq.2)) then
              Ef=bias_L
       elseif (FermiP.eq.3) then
              Ef=0
       endif

       print *, ' The Fermi energy is then',Ef
       print *, ' VDC = mu_L-mu_R=', bias_L - bias_R, 'Applied to the tip'
       print *, ' '
       print *, ' The couplings are:'
       print *, '   Left coupling  =', gamma_L_0,'   Right coupling =', gamma_R_0, 'in meV'
       print *, ' '
       print *,'Broadening (meV), Cutoff (meV), Temperature (K), Polarization L,R:'
       print *, gammaC,Cutoff,Temperature,Spin_polarization_L,Spin_polarization_R 

       call atomic_units (GammaC, Cutoff,&
         gamma_R_0, gamma_L_0,Temperature, bias_R, &
         bias_L, Ef, frequency, B_L, B_R, &
         transport_exponentR, transport_exponentL) ! Freq transformed from frequency to radial freq

       VDC = bias_L - bias_R

       return 
    end subroutine reading_input

!
! change into atomic units using the parameters of declarations.f90
!

   subroutine  atomic_units (GammaC, Cutoff,&
         gamma_R_0, gamma_L_0,Temperature, bias_R, &
         bias_L, Ef, frequency, B_L, B_R, &
         transport_exponentR, transport_exponentL) ! Freq transformed from frequency to radial freq
       implicit none
    integer, parameter :: q = SELECTED_REAL_KIND(10)
    integer, parameter :: qc = SELECTED_REAL_KIND(10)
    real (q), parameter :: pi_d = 3.14159265358979323846_q
    real (q), parameter :: Hartree = 27211.386245988_q ! meV
    real (q), parameter :: BohrMagneton=5.7883818060E-5_q !eV/T
    real (q), parameter :: GHz = 4.135667696786E-3_q ! meV 
    real (q), parameter :: time_unit = 2.418884326505E-8_q ! nanoseconds
    real (q), parameter :: pA = .6623618237510E10_q !a.u. to pA 
   real (q) :: Temperature, Ef, frequency
   real (q) :: gamma_R_0, gamma_L_0!, gamma_L_1, gamma_R_1
   real (q) :: bias_R, bias_L, B_L, B_R
   real (q) :: transport_exponentR, transport_exponentL
   real (q) :: tolerance, Cutoff
   logical :: runs
   complex (qc) :: gammaC

       tolerance = tolerance/time_unit
       gamma_R_0 = gamma_R_0/Hartree
       gamma_L_0 = gamma_L_0/Hartree
       gammaC = gammaC/Hartree
       Cutoff=Cutoff/Hartree
       bias_R = bias_R/Hartree
       bias_L = bias_L/Hartree
       Temperature = Temperature * 25.852_q / (Hartree*300._q)
       frequency = frequency*2.*pi_d*time_unit
       Ef=Ef/Hartree
       B_L = B_L/Hartree
       B_R = B_R/Hartree
       transport_exponentR = transport_exponentR*Hartree
       transport_exponentL = transport_exponentL*Hartree

   end subroutine atomic_units


end module io

