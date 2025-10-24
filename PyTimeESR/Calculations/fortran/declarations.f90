module declarations

  implicit none

! PARAMETERS

!TODO: some repetions some unused variables
  integer, parameter :: kk = SELECTED_INT_KIND (10)
  integer, parameter :: q = SELECTED_REAL_KIND(10)
  integer, parameter :: qs = SELECTED_REAL_KIND(5)
  integer, parameter :: qc = SELECTED_REAL_KIND(10)
  real (q), parameter :: pi_d = 3.14159265358979323846_q
  real (q), parameter :: sqr2=1.41421356237309504880_q
  real (q), parameter :: sqr3=1.73205080756887729352_q
  real (q), parameter :: sqr5=2.23606797749978969640_q
  real (q), parameter :: sqr6=2.44948974278317809819_q
  real (q), parameter :: sqr7=2.64575131106459059050_q
  real (q), parameter :: sqr8=2.82842712474619009760_q
  real (q), parameter :: Hartree = 27211.386245988_q ! meV
  real (q), parameter :: BohrMagneton=5.7883818060E-5_q !eV/T
  real (q), parameter :: GHz = 4.135667696786E-3_q ! meV 
  real (q), parameter :: time_unit = 2.418884326505E-8_q ! nanoseconds
  real (q), parameter :: pA = .6623618237510E10_q !a.u. to pA 
  complex (qc), parameter :: zero=(0._q,0._q), ui = (0._q,1._q)
  complex (qc), parameter ::  one = (1._q, 0._q)
  complex (qc), parameter ::  ur = (1._q, 0._q)

! numbers
  integer :: NF, INFO, ITER, LDA, LDB, LDX, Nmatrix, NRHS
  integer :: Nd, i_m, Nplot, Ndim, N_int, Ndim_old, p_maxR, p_maxL, pmax
  integer :: i,ii, N, Nm, Np, i_, j, j1, j2, l, i_omega, k, FermiP, p, p_ind
  integer :: j3, j4, i1, i2, i3, i4, u, v, i_sigma, i2p, n_index, nfour
  integer :: Electrode, NCF,orb,i_feed
  real (q) :: eps_QD, U_Hubbard, p_mod,phi,Ef,sum_rule
  real (q) :: px, py, pz, pxx, pyy, pzz, suma, omega, WW, gau
  real (q) :: bias_R, bias_L, bias, Spin_polarization_R, Spin_polarization_L, Spin_polarization
  real (q) :: transport_exponentR, transport_exponentL, transport_exponent
  real (q) :: Temperature, gamma_R_0, gamma_R_1, gamma_L_0, gamma_L_1, Cutoff, gamma_0
  real (q) :: VDC, Bdrive, B_L, B_R, frequency
  complex (qc) :: A_L, A_R, spin2_ave, GammaC, Adrive
  complex (qc) :: g_dn, g_up, bessel_contribution, ubessel_contribution
  real (q) :: e, eFermi, step_e, esq, gaushift, WWsq, gausian, f, uf, te
  real (q) :: imG, rGammaC, iGammaC, rGammaCsq, iGammaCsq
  real (q) :: D

! arrays
  integer, allocatable :: IPIV (:)
  integer, allocatable :: mol1(:), mol2(:), N_in (:), N_block (:)
  real (q), allocatable :: RWORK (:),curr(:)
  real (q), allocatable :: W (:), Jexch (:,:), nx(:), ny(:), nz(:)
  real (q), allocatable :: hx (:), hy (:), hz (:)
  real (q), allocatable :: gx (:), gy (:), gz (:)
  real (q), allocatable :: B20 (:), B22 (:), B40 (:), B44 (:)
  real (q), allocatable :: Spin (:), Eigenvalues (:), Delta (:,:)
  real (q), allocatable :: Amplitude_seq (:,:), Freq_seq (:,:)
  real (q), allocatable :: Phase_seq (:) 
  complex (qc), allocatable :: WORK (:), SWORK (:), A_old(:,:)
  complex (qc), allocatable :: A_fast_left(:,:), A_fast_right(:,:)
  complex (qc), allocatable :: Ss(:,:,:,:), Sn (:,:,:), H_el(:,:)
  complex (qc), allocatable :: H (:,:), SprodS (:,:,:), SprodS2 (:,:,:)
  complex (qc), allocatable :: Identity (:,:), A (:,:), B (:), Rho2 (:,:), XX(:)
  complex (qc), allocatable :: Sp (:,:), Sm (:,:)
  complex (qc), allocatable :: Sp2 (:,:), Sm2 (:,:)
  complex (qc), allocatable :: Sp4 (:,:), Sm4 (:,:)
  complex (qc), allocatable :: Sx(:), Sy(:), Sz(:), Sh(:), S2(:)
  complex (qc), allocatable :: spinX (:,:,:), spinY(:,:,:), spinZ(:,:,:)
  complex (qc), allocatable :: Sx_u (:,:,:), Sy_u (:,:,:), Sz_u (:,:,:)
  complex (qc), allocatable :: spin2_T(:,:), spin2(:,:,:)
  complex (qc), allocatable :: lambda (:,:,:), Lvluj(:,:), Ljulv(:,:)
  complex (qc), allocatable :: G (:,:,:,:,:,:), GC (:,:,:,:,:,:), GA (:,:,:,:,:), GCA (:,:,:,:,:)
  complex (qc), allocatable :: rho (:,:,:)
  complex (qc), allocatable :: fermiR_a(:,:), fermiL_a(:,:), fermi(:)
  complex (qc), allocatable :: ufermiR_a(:,:), ufermiL_a(:,:), ufermi(:)
  complex (qc), allocatable :: Kbess(:), Jbess(:)
  complex (qc), allocatable :: QME_Tensor_A(:,:,:,:,:,:), QME_Tensor_R(:,:,:,:,:,:), QME_Tensor_L(:,:,:,:,:,:)
  complex (qc), allocatable :: QME_Tensor_D(:,:,:,:,:,:), QME_Tensor(:,:,:,:,:,:)

! logical
  logical :: presence,runs,redimension,faster,presence2
  logical :: write_populations, write_coherences, spinflo

! character
   character ( len = 100 ) :: Name_output, output_file, output_fourier, output_ESR, UPLO, filename


end module declarations
