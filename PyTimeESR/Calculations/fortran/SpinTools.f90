module SpinTools
Use OpenFiles
Use declarations, only: q, qc, zero
CONTAINS
    subroutine SpinFloquet(Nm, Ndim, Ss, spinX, spinY, spinZ,spin2_T, H, rho, hx, hy, hz,&
    &Sx,Sy,Sz,Sh,spin2_ave)
     implicit none
     integer :: Nm, Ndim, i,ii, j, m, Ndim_old
     real (q), intent (in) :: hx(:), hy(:), hz(:)
     real (q) :: mod_h
     complex (qc), intent(out), allocatable :: Sx(:),Sy(:),Sz(:),Sh(:)
     complex (qc), intent(out) :: spin2_ave
     complex (qc), intent (in) :: Ss (:,:,:,:), H (:,:), rho(:,:)
     complex (qc), intent (in) :: spinX(:,:,:),spinY(:,:,:),spinZ(:,:,:),spin2_T(:,:)

     
     
    Ndim_old=size(H(1,:))
    spin2_ave=zero
     
     do i = 1, Ndim
     do ii = 1, Ndim
        spin2_ave=spin2_ave+spin2_T(i,ii)*rho(ii,i)
!         write(*,*)spin2_ave 
     enddo
     enddo
     
     
        allocate (Sx(Nm),Sy(Nm),Sz(Nm),Sh(Nm))
        Sx = zero; Sy=zero; Sz=zero

     
    do m=1,Nm ! loop over the molecules
      ! Trace over density matrix rho
     do i = 1, Ndim
     do ii = 1, Ndim
      Sx (m) = Sx (m) + spinX (i,ii,m)*rho(ii,i)
      Sy (m) = Sy (m) + spinY (i,ii,m)*rho(ii,i)
      Sz (m) = Sz (m) + spinZ (i,ii,m)*rho(ii,i)
     enddo
     enddo
     
!     do m=1,Nm ! loop over the molecules
! write(222,*)MATMUL(spinX(:,:,m)*rho(:,:)),MATMUL(spinY(:,:,m)*rho(:,:)),MATMUL(spinZ(:,:,m)*rho(:,:))
!      enddo
!      enddo
    
! Projection of spin along the local magnetic field axis
      mod_h = sqrt(hx (m)**2 + hy (m)**2 + hz (m)**2)
      if (mod_h==0) then
      print *, ' '
      print *, 'WARNING!! Magnetic field on spin number:', m
      print *, 'WARNING!! is ZERO! We cannot project along that axis!'
      print *, 'WARNING!! We set the projection to axis z'
      print *, ' '
      Sh (m) = Sz (m)
      else
      Sh (m) = Sx (m)*hx (m) + Sy (m)*hy (m) + Sz (m)*hz (m) 
      Sh (m) = Sh (m) / mod_h
      endif
      
     enddo
!         deallocate(Sx,Sy,Sz,Sh,S2)
    
      return

      end subroutine SpinFloquet
end module SpinTools
