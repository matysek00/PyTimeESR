module H_QD
Use declarations
Use OpenFiles
Use algebra
Use KrondProd
implicit none

CONTAINS

   subroutine create_basis(Nm, eps_QD, U_Hubbard, Spin, N, H_el, N_in, N_block)
    integer, intent(in) :: Nm 
    real(q), intent(in) :: Spin(Nm)
    real(q), intent(inout) :: eps_QD, U_Hubbard
    
    integer, intent(out) ::  N, N_in(Nm), N_block(Nm)
    real(q), intent(out) :: H_el(4,4)

    integer :: i_m


        H_el = zero
        H_el(1,1) = eps_QD; H_el(2,2) = eps_QD ! energy of the one electron state
        H_el(3,3)=0; H_el(4,4) = 2*(eps_QD)+U_Hubbard

! the basis here is a tensorial product of the first state times the second times...
! each state is given by a spin matrix
! the spin matrix is called Ss (:,:,:,:)
!                  the first entry Ss (i,:,:,:) refers to the site
!                  the second entry Ss (:,i,:,:) refers to the x, y, z component of the spin
!                  the third and fourth entry coincide with the Hamiltonian dimensions Hamiltonian (:,:) 

! Total dimension (Hamiltonian entries) :: N
! each spin is of dimension N_in
! up to a spin i_m the total dimension is the product of previous N_in, stored in N_block
          N_in(1)=4  !int(2*Spin(1)+1)+1+1=4
          N_block (1) = N_in (1)
! Test to remove after debugging
!         print *, 'N_in',1, N_in(1)
!         print *, 'N_block',1, N_block(1)
       do i_m = 2, Nm
          N_in (i_m)=int(2*Spin(i_m)+1)
          N_block (i_m) = N_block (i_m-1) * N_in (i_m)
!         print *, 'N_in',i_m, N_in(i_m)
!         print *, 'N_block',i_m, N_block(i_m)
       enddo
          N = N_block (Nm) ! full dimension
          
          print*,''
          print*, 'Key energies for transport:'
          print*, 'Ionization energy:',eps_QD*Hartree
          print*, 'Coulomb repulsion:',U_Hubbard*Hartree
          print*, 'On site energy plus Ionization:',(U_Hubbard+eps_QD)*Hartree
          print*,''
    
   end subroutine create_basis


   subroutine Hamiltonian (Nm,Np,N,N_in, N_block,&
    &H_el, hloc, spinaxis, Spin, B20, B22, B40, B44,&
    &Jexch, mol1, mol2,& 
    & lambda, Delta, H, W, Ss)!&
!    &,spin_first)
   implicit none
   !F2PY INTENT(HIDE) :: Nm
   !F2PY INTENT(HIDE) :: Np
   !F2PY INTENT(IN)   :: N
   !F2PY INTENT(IN)   :: N_in(Nm)
   !F2PY INTENT(IN)   :: N_block(Nm)
   !F2PY INTENT(IN)   :: H_el(4,4)
   !F2PY INTENT(IN)   :: hloc(Nm, 3)
   !F2PY INTENT(IN)   :: spinaxis(Nm, 3)
   !F2PY INTENT(IN)   :: Spin(Nm)
   !F2PY INTENT(IN)   :: B20(Nm)
   !F2PY INTENT(IN)   :: B22(Nm)
   !F2PY INTENT(IN)   :: B40(Nm)
   !F2PY INTENT(IN)   :: B44(Nm)
   !F2PY INTENT(IN)   :: Jexch(Np, 3)
   !F2PY INTENT(IN)   :: mol1(Np)
   !F2PY INTENT(IN)   :: mol2(Np)
   !F2PY INTENT(OUT)  :: lambda (N,N,2)
   !F2PY INTENT(OUT)  :: Delta
   !F2PY INTENT(OUT)  :: H(N,N)
   !F2PY INTENT(OUT)  :: W(N)
   !F2PY INTENT(OUT)  :: SS (Nm,4,N,N)

   integer, intent(in) :: Nm, Np, N
   integer, intent(in) :: mol1(Np), mol2(Np), N_in(Nm), N_block(Nm)

   real(q), intent(inout) :: hloc(Nm,3), B20(Nm), B22(Nm), B40(Nm), B44(Nm), Jexch(Np, 3), spinaxis(Nm,3) 
   real(q), intent (in) ::  Spin(Nm), H_el(4,4)

   real (q), intent (out) ::  W (N), Delta (N,N)
   complex (qc), intent (out) ::  H(N,N), lambda (N,N,2), Ss (Nm,4,N,N)
   
!
! Create many body configurations for spin excitation
! or ESR dynamics
!
! Brute force diagonalization of a model Hamiltonian
! that contains
!         spins, anisotropic exchange, local magnetic fields, Stephen Operators,
!         an exchange interaction between electron site and first spin
!         finite Hubbard U
!
!
! gnu licence 3.0 (c) J. Reina Galvez & N. Lorente
!
!
! in this code, sites or spins are considered to be molecules, this explains the notation


   ! finally, information on the electronic level:
  ! TODO: move this check to python 
     if (int(2*Spin(1)+1) /= 2) then
        write (*, *) ' '
        write (*,*) 'ERROR: '
        write (*,*) ' You are using a spin different from 1/2 for the transport electron!'
        write (*,*) ' Stop. '
        write (*,*) '  '
        stop
     endif

      
! MEMORY
      
      allocate (Identity (N,N))
      allocate (Sn (3, N, N)) ! rotated spin matrix
      allocate (SProdS (3, N,N)) ! The following are operated spins 
      allocate (SProdS2 (3, N,N)) 
      allocate (Sp (N,N))
      allocate (Sm (N,N))
      allocate (Sp2 (N,N))
      allocate (Sm2 (N,N))
      allocate (Sp4 (N,N))
      allocate (Sm4 (N,N))

      
! initializations
     H = zero
     Ss = zero
     Sn = zero
     SProdS = zero
     SProdS2 = zero
     Identity = zero
        do i =1,N
           Identity (i,i) = ur
        enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                                                     !
!    We generate the tensorial product of spin matrices (aka the basis set):                          !
!                                                                                                     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      call Kronecker_Product (N, Nm, N_in, N_block, Ss, Sx_u, Sy_u, Sz_u)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                                                     !
!    We generate the Hamiltonian                                                                      !
!                                                                                                     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     STEPS:
!
!     ZERO: The electronic contribution to the Hamiltonian
!           we use the same algo as for the spin Ss

       do i2 = 1, N_in (1)
        do i3 = 1, N_in (1)

            do i4 = 1, N/N_block(1)

       H (i4 + (i2-1)*(N/N_block(1)),  &
       &  i4 + (i3-1)*(N/N_block(1)) )  =  &
       &  H_el(i2,i3);

            enddo
         enddo
        enddo


!     FIRST: local magnetic fields

      lmol: do i_ = 1, Nm ! Loop on molecules or sites
   
      H (:,:) = H (:,:)+Ss (i_, 1, :,:)*hloc(i_,1)+Ss (i_, 2, :,:)*hloc(i_,2)+Ss (i_, 3, :,:)*hloc(i_,3)

!     SECOND: local Stephen operators
!             We need to compute spin matrices up to 4th power.

      ! We define new directions following the preferential axis given in the input file
      ! in this way we can take the internal "z" axis to be always along the preferential axis
      ! First normalize vector input:
          p_mod=sqrt(spinaxis(i_,1)**2+spinaxis(i_,2)**2+spinaxis(i_,3)**2)
          spinaxis(i_,:) = spinaxis(i_,:)/p_mod
      ! we make vectors px,py,pz,pxx,... perpendicular to nx,ny,nz
       if (spinaxis(i_,3) /= 0) then
          if (spinaxis(i_,2) /= 0 ) then
          px=-spinaxis(i_,2); py=spinaxis(i_,1); pz=0
          p_mod=sqrt(px**2+py**2); px=px/p_mod; py=py/p_mod
          pxx=spinaxis(i_,1); pyy=spinaxis(i_,2);pzz=-(spinaxis(i_,1)**2+spinaxis(i_,2)**2)/spinaxis(i_,3);
          p_mod=sqrt(pxx**2+pyy**2+pzz**2); pxx=pxx/p_mod; pyy=pyy/p_mod; pzz=pzz/p_mod
          else if (spinaxis(i_,1) /=0) then
          px=-spinaxis(i_,2); py=spinaxis(i_,1); pz=0
          p_mod=sqrt(px**2+py**2); px=px/p_mod; py=py/p_mod
          pxx=spinaxis(i_,1); pyy=spinaxis(i_,2);pzz=-(spinaxis(i_,1)**2+spinaxis(i_,2)**2)/spinaxis(i_,3);
          p_mod=sqrt(pxx**2+pyy**2+pzz**2); pxx=pxx/p_mod; pyy=pyy/p_mod; pzz=pzz/p_mod
          else
          pxx = 1.0; pyy = 0.0; pzz=0.0; px = 0.0; py = 1.0; pz=0.0
          endif
       else 
          px=-spinaxis(i_,2); py=spinaxis(i_,1); pz=0
          p_mod=sqrt(px**2+py**2); px=px/p_mod; py=py/p_mod
          pxx=0; pyy=0;pzz=1._q
       endif


      ! rotated spins, inside the site loop on i_

       Sn (3,:,:)=Ss(i_,1,:,:)*spinaxis(i_,1)+Ss(i_,2,:,:)*spinaxis(i_,2)+Ss(i_,3,:,:)*spinaxis(i_,3)
       Sn (2,:,:)=Ss(i_,1,:,:)*px+Ss(i_,2,:,:)*py+Ss(i_,3,:,:)*pz
       Sn (1,:,:)=Ss(i_,1,:,:)*pxx+Ss(i_,2,:,:)*pyy+Ss(i_,3,:,:)*pzz

      ! We perform matrix multiplication and use the definitions of Stephens Operators
      ! choosing a few ones up to 4th order -See our paper in J. Phys. Chem. A 124, 2318  (2020)
      ! We add over all molecules the contribution of each local anisotropy, this seems
      ! to work, at least for not extremely coupled spins


        call  MatSquareSpin (N, Sn, SprodS)
        call  MatSquareSpin (N, SprodS, SprodS2)

        H (:,:) = H (:,:)+ B20(i_)*SprodS(3,:,:)  !Longitudinal anisotropy !3 being include in B20 -> 3.*B20(i_)*SprodS(3,:,:) 
        H (:,:) = H (:,:)+ B22(i_)*(SprodS(1,:,:)-SprodS(2,:,:)) ! Transversal anisotropy
!         H (:,:) = H (:,:)+ B22(i_)*(SprodS(1,:,:)+SprodS(2,:,:)) ! Transversal anisotropy
        H (:,:) = H (:,:)+ B40(i_)*35.*SprodS2(3,:,:) ! Fourth order: B40 ain't pretty
        H (:,:) = H (:,:)- B40(i_)*30.*(Spin(i_)*(Spin(i_)+1._q)-2*Spin(i_))*SprodS(3,:,:)
        H (:,:) = H (:,:)+ B40(i_)*(3*(Spin(i_)*(Spin(i_)+1))**2-6*(Spin(i_)*(Spin(i_)+1)))*Identity (:,:)

        ! change to circular spins: S+ and S-
        Sp (:,:) = Sn (1, :,:)+ui*Sn (2, :,:)
        Sm (:,:) = Sn (1, :,:)-ui*Sn (2, :,:)

        call MatSpSm (N, Sp, Sm, Sp2, Sm2)
        call MatSpSm (N, Sp2, Sm2, Sp4, Sm4)

        H = H + B44(i_)*0.5*(Sp4+Sm4) ! Fourth order: B44, sort of fourth-order longitudinal anisotropy
!         H = H + B22(i_)*0.5*(Sp2+Sm2) ! Transversal anisotropy
!         H = H - B44(i_)*0.5*ui*(Sp4-Sm4) ! Fourth order: B44, sort of fourth-order longitudinal anisotropy
   
      enddo lmol
! Test to remove after debugging
!     do i = 1, N
!         write (*, *) (j,H(i,j), j=1, N)
!     enddo


!     THIRD: anisotropic intermolecular exchange interactions
!     we keep the original axis, not the one of the anisotropy

      lpair: do i = 1, Np ! loop on pairs

!     if (i ==1) then
!     print *, 'computing first pair!!!'
!     else if (i == 2) then
!     print *, 'computing second pair!!!'
!     else
!     print *, 'computing', i,'th pair!!!'
!     endif

      call MatProdSpin (N, i, mol1, mol2, Ss, SProdS) ! mol1 and mol2 contain the indices of the paired molecules

! Test: remove after debugging
!     print *, 'First Spin'
!     do j1 =1,N
!      write (*, '(16g14.4)') (j, Real(Ss(1,3,j1,j)), j=1, 8)
!     enddo
!     print *, 'Second Spin'
!     do j1 =1,N
!      write (*, '(16g14.4)') (j, Real(Ss(2,3,j1,j)), j=1, 8)
!     enddo
!     print *, 'Product:'
!     do j1 =1,N
!      write (*, '(16g14.4)') (j, Real(SProdS(3,j1,j)), j=1, 8)
!     enddo

      H (:,:) = H (:,:)+Jexch (i, 1)*SProdS (1,:,:)+Jexch (i, 2)*SProdS (2,:,:)+Jexch (i, 3)*SProdS (3,:,:)

      enddo lpair
! Test: remove after debugging
!      print *, 'The Hamiltonian is:'
!    do i = 1,N
!      write (*,'(16g14.4)')  (j,Real(H(i,j)), aimag(H(i,j)), j=1,N)
!    enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                                                     !
!        DIAGONALIZE                                                                                  !
!                                                                                                     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!    print *, ' '
!    print *, 'Beginning of diagonalization'
!    print *, ' '
!          print *, 'size H:', size(H,1), N
!     do i =1, N
!     write(126,*) i,(dble(H(i,j)*hartree),aimag(H(i,j)*hartree), j= 1, N)
!     enddo
      call DIAGONALIZE (N,H,W)
      
      do i = 1, N
      do j = 1, N
        Delta (i,j) = W(i)-W(j)
      enddo
      enddo

      ! Calculation of lambda
      
      lambda = zero

      do i = 1, N
      do j = 1, N
      do i_sigma= 1, 2 !1 spin down as always in this code

! contribution from |0Xsigma|

      i2 = 3 ! 3 is |0>
      i3 = i_sigma ! 1 is down and 2 is up as corresponds to the basis set

      do i4 = 1, N/N_block(1)

       lambda (i,j,i_sigma) = lambda (i,j,i_sigma) +  &
         conjg(H (i4 + (i2-1)*(N/N_block(1)), i)) * H (i4 + (i3-1)*(N/N_block(1)), j)

      enddo

! contribution from |\bar{sigm}aX4| where 4 is the doubly occupied state (singlet)

      i2 = (-1)**(i_sigma+1)+i_sigma! if i_sigma=1 then this yields 2
                            ! if i_sigma=2 then this yields 1
      i3 = 4 ! 4 is the doubly occupied state

      do i4 = 1, N/N_block(1)

       lambda (i,j,i_sigma) = lambda (i,j,i_sigma) +  &
         conjg(H (i4 + (i2-1)*(N/N_block(1)), i)) * H (i4 + (i3-1)*(N/N_block(1)), j)

      enddo


    enddo
! test lambda
!      write (123,*) i,j, 1, lambda (i,j,1)
!      write (123,*) i,j, 2, lambda (i,j,2)

    enddo
    enddo
    
   print *, 'The Hamiltonian calculation is DONE!'
   write (*,*) '**************************************************'

  end subroutine Hamiltonian
  

  subroutine spin_states(N, Nplot, Nm, H, Ss, spin_plot, spin2_T)
    implicit none
    !F2PY INTENT(HIDE) :: N
    !F2PY INTENT(HIDE) :: Nm
    !F2PY INTENT(IN) :: Nplot
    !F2PY INTENT(IN) :: H
    !F2PY INTENT(IN) :: SS
    !F2PY INTENT(IN) :: spin_plot(Nplot, Nplot,Nm,3)
    !F2PY INTENT(IN) :: spin2_T(Nplot,Nplot)
    integer, intent(in) :: N, Nm
    integer, intent(inout) :: Nplot
    complex (qc), intent(in) :: Ss(Nm, 4, N,N), H(N, N)
    complex (qc), intent(out) :: spin_plot(Nplot, Nplot, Nm, 3), spin2_T(Nplot, Nplot)


  if (Nplot > N) then
    Nplot = N
  endif

  spin_plot = zero
  
  do j = 1, Nm
  do i=1, Nplot
  do ii=1, Nplot
  do j1=1, N
  do j2=1, N
    spin_plot(i,ii,j,:) = spin_plot(i,ii,j,:)+conjg(H(j1,i))*Ss (j, :3, j1, j2)*H(J2, ii)
  enddo
  enddo
  enddo
  enddo 
  enddo

    spin2_T = zero  ! spin square
  
    do i = 1, Nplot
    do ii= 1, Nplot
    do j = 1, Nm
    do j4 = 1, Nm
    do j1=1, N
    do j2=1, N
    do j3=1,N
      spin2_T(i,ii)=spin2_T(i,ii)+conjg(H(j1,i))*(Ss (j, 1, j1, j3)*Ss (j4, 1, j3, j2)+  &
        &        Ss (j, 2, j1, j3)*Ss (j4, 2, j3, j2)+ Ss (j, 3, j1, j3)*Ss (j4, 3, j3, j2) &
        &        )*H(J2, ii)
    enddo
    enddo
    enddo
    enddo
    enddo
    enddo        
    enddo

  end subroutine spin_states

  subroutine WriteEigen (Ndim, H, W)
    implicit none
    !F2PY INTENT(HIDE) :: Ndim
    !F2PY INTENT(IN) :: H
    !F2PY INTENT(IN) :: W
    integer, intent(in) :: Ndim
    complex(qc), intent(in) :: H(Ndim, Ndim), W(Ndim)
   print *, 'States and energies written in: Eigenvalues and eigenvectors .dat'
!    print *, ' '
! We go ahead and write the ouput 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                                                     !
!        Write output in a useful way                                                                 !
!                                                                                                     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    print *, 'Eigenvalues in meV are written in Hamiltonian_Eigen.dat'
!    print *, ' '

        print *, 'These are the eigenvalues (GHz)::'
        write (*, '(8g14.4)') ((W(i))*Hartree/GHz, i=1,N)
        print *, 'Eigenvalues (meV)::'
        write (*, '(8g14.4)') ((W(i))*Hartree, i=1,N)
        print *, ' '

       open (unit_spinoutput, file='Eigenstates.dat')
       open (unit_spinoutput+1, file='Eigenvalues.dat')
       open (unit_spinoutput+2, file='Matrix_transition.dat')

       do i =1, Ndim
          write (unit_spinoutput,*) i,(real(H(j,i)),aimag(H(j,i)), j= 1, Ndim) ! IF REDIMENSION THEN EXTRA ZEROS WILL APPEAR FROM THE OLD BASIS CUT
          write (unit_spinoutput+1,*) i,hartree*W(i),hartree*W(i)/GHz,hartree*(W(i)-W(1)),hartree*(W(i)-W(1))/GHz
       enddo

        write(unit_spinoutput+2,*)
        write(unit_spinoutput+2,*) 'Energy transition between the different Hamiltonian states (in matrix form)::'
        write(unit_spinoutput+2,*)
        write(unit_spinoutput+2,*) '', (j, j=1,Ndim)
        do i=1,Ndim
            write(unit_spinoutput+2,*) i,(hartree*(W(i)-W(j)),j=1,Ndim)
        end do

        write(unit_spinoutput+2,*) ''
        write(unit_spinoutput+2,*) 'Energy transition between the different Hamiltonian states (in matrix form)::'
        write(unit_spinoutput+2,*)
        write(unit_spinoutput+2,*) '', (j, j=1,Ndim)        
        do i=1,Ndim
            write(unit_spinoutput+2,*) i,((W(i)-W(j))*Hartree/GHz,j=1,Ndim)
        end do
       close (unit_spinoutput)
       close (unit_spinoutput+1)
       close (unit_spinoutput+2)
  end subroutine WriteEigen

! Plotting the spin distribution for the first Nplot states 

subroutine WriteSpinDist (N, Nplot, Nm, W, spin_plot, spin2_T)
        implicit none
        !F2PY INTENT(HIDE) :: N
        !F2PY INTENT(HIDE) :: Nplot
        !F2PY INTENT(HIDE) :: Nm
        !F2PY INTENT(IN) :: W(N)
        !F2PY INTENT(IN) :: spin_plot(Nplot, Nplot, Nm, 3)
        !F2PY INTENT(IN) :: spin2_T(Nplot, Nplot)

        integer, intent(in) :: N, Nm
        integer, intent(inout) :: Nplot
        real(q), intent(in) ::  W(N)
        complex (qc), intent(in) :: spin_plot(Nplot, Nplot, Nm, 3), spin2_T(Nplot, Nplot)

    if (Nplot > N) then
      Nplot = N
    endif
    
    open (unit_spinoutput, file='Spin_distribution.dat')
    
    do i=1, Nplot
       write (unit_spinoutput,*) 'State=',i
       write (unit_spinoutput,*) '#  Site, Sx, Sy, Sz'
    enddo

    do i=1,Nplot
        write (unit_spinoutput,*) i, (real(spin_plot(i,i,j,1)), real(spin_plot(i,i,j,2)), real(spin_plot(i,i,j,3)),j=1,Nm)
    enddo
      
    print *, 'Spins written in file:  Spin_distribution.dat'
    print *, ' '        

    close (unit_spinoutput)

    open (unit_spinoutput, file='Full_energies_and_spin.dat') ! costy writing if total spin is large
    write (unit_spinoutput,*) ' State ', ' Excitation Energy (GHz) ', ' (meV) ', ' Spin^2 ', ' Spin '

    do i=1, Nplot
    write (unit_spinoutput,*) i, (W(i))*Hartree/GHz, (W(i))*Hartree, &
        real(spin2_T(i,i)), 0.5*(sqrt(1+4*real(spin2_T(i,i)))-1)
    enddo
    
!     do i=1, Nplot
!         spin_first=real(spinX(i,1))**2+ real(spinY(i,1))**2 +real(spinZ(i,1))**2
!     enddo
    
    
    close (unit_spinoutput)
  end subroutine WriteSpinDist
    
!
! Cutting off dimensions for sizeable calculations
! 
     subroutine redimensioning (Ndim, Delta, bias, Cutoff)
     implicit none
     integer, intent(inout) :: Ndim
     integer :: i, l
     real (q), intent(in) :: bias, Cutoff
     real (q), dimension (:,:), intent(in) :: Delta

       l = 0
        
     do i = 1, Ndim

       if (Delta (i,1) <= bias+10*Cutoff) then
           l = l+1
       endif

     enddo

    write (*,*) 'Redimensioning from',Ndim,'to'
    Ndim = l
    write(*,*) 'New dimension',l,'based on energy differences.'
    write (*,*) ' '


     return
     end subroutine redimensioning


end module H_QD


! TODO: move to PYTHON 
!!!! Begin Printing 

!
! Save it, unformatted so it does not need to be recalculated if it 
! exists on the running folder
!
!      open (unit_spinoutput, file=Name_output, form='unformatted')
!      write (unit_spinoutput) H, W
!      close (unit_spinoutput)
!
!!!! End Printing 

!!! Begin Redimension 
!    Ndim=N  
!    p=0
!    if ((redimension).and.(Nd.gt.N/2).and.(Nd<N)) then
!        
!!     write (*,*) ' '
!    write (*,*) 'Redimensioning from',N,'to New dimension',Nd,'from the input file.'
!     Ndim = Nd
!!     write (*,*) 'Extra zeros will appear in the eigenvectors due to the old basis being cut'
!!
!    p=1 ! it will indicate that a redimensionalition was done so we do not do another one
!    write (*,*) ' '
!    else
!    write (*,*) ' '
!    write(*,*) 'Input dimension',Nd,' below min dimension',N/2,'above or equal max dimension',Ndim
!    write(*,*) 'or redimension is set to be .FALSE.. Redimensioning based on energies differences'
!!     redimension=.false.
!    write (*,*) ' '
!   endif
!    
!    ! redimension of the problem if the input redimension does not apply:
!    ! Nd is a number larger than N or more states are disconnected
!    
!    if (abs(bias_L).lt.abs(bias_R))then
!        bias=abs(bias_R)
!    else
!        bias=abs(bias_L)
!    endif
!
!    if ((redimension).and.(p.eq.0)) then
!        call redimensioning (Ndim, Delta, bias, Cutoff)
!    elseif (p.eq.0) then
!        Ndim=N
!        write (*,*) 'Full dimension is taken',N
!        write (*,*) ' '
!    endif
!!! End Redimensions 