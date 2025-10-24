module output
use declarations
use OpenFiles !all units of files in open statements are here
use SpinTools
CONTAINS

    subroutine write_rate_out(G, NCF, Ndim)
    implicit none
    ! Declare all arguments and local variables
    complex(q), intent(in) :: G(:,:,:,:,:,:)
    integer, intent(in) :: NCF, Ndim

    integer :: l, j, u, v, pn
    complex(qc) :: G_temp(2*NCF+1, 2)

    open (unit_rates, file='rates_floquet.dat')
    open (unit_rates+1, file='rates_floquet0.dat')

    write (unit_rates, *) 'l j u v n Re(GL_n)(meV) Im(GL_n)(meV) Re(GR_n)(meV) Im(GR_n)(meV)'
    write (unit_rates+1, *) 'l j u v Re(GL_0)(meV) Im(GL_0)(meV) Re(GR_0)(meV) Im(GR_0)(meV)'
    do l=1,Ndim
    do j=1,Ndim
    do u=1,Ndim
    do v=1,Ndim
        G_temp = G(l,j,u,v,:,:)*hartree ! convert to Hartree

        if ((G_temp(NCF+1,1)+G_temp(NCF+1,2)) == (0.0_q,0.0_q)) then
            !           if the rates are zero, we do not write them
            cycle
        endif
        
!       writing the zeroth component explicitely as in the loop it would be written twice
        write (unit_rates,*) l,j,u,v,0,& 
                dble(G_temp(NCF+1,2)), dimag(G_temp((NCF+1),2)),&
                dble(G_temp(NCF+1,1)), dimag(G_temp((NCF+1),1))

        pnloop: do pn = 1, NCF
            if ((G_temp(pn+NCF+1,1)+G_temp(pn+NCF+1,2)).eq.zero) then
!           write until reaching a vanishing fourier component. 
            exit pnloop
            endif
        
            write (unit_rates,*) l,j,u,v,pn,& 
                dble(G_temp(pn+NCF+1,2)), dimag(G_temp(pn+(NCF+1),2)),&
                dble(G_temp(pn+NCF+1,1)), dimag(G_temp(pn+(NCF+1),1))

            write (unit_rates,*) l,j,u,v,-pn,& 
                dble(G_temp(-pn+NCF+1,2)), dimag(G_temp(-pn+(NCF+1),2)),&
                dble(G_temp(-pn+NCF+1,1)), dimag(G_temp(-pn+(NCF+1),1))

        enddo pnloop
            

!       simplified version of writing the rates, only 0th component  
        write (unit_rates+1,*) l,j,u,v,& 
            dble(G_temp(NCF+1,2)), dimag(G_temp(NCF+1,2)),&
            dble(G_temp(NCF+1,1)), dimag(G_temp(NCF+1,1))
        
    enddo
    enddo
    enddo
    enddo

    close(unit_rates)
    close(unit_rates+1)

    return
end subroutine write_rate_out


subroutine check_populations(Rho, NCF, Ndim, bias_R, bias_L, frequency)
    implicit none
    real(q), intent(in) :: bias_R, bias_L, frequency
    complex(qc), intent(in) :: Rho(:,:,:)
    integer, intent(in) :: NCF, Ndim
    integer :: l
    real(q) :: sum_rule

    open (unit_error, file='strange_populations.dat')
    do l=1,Ndim

        sum_rule  = sum_rule  + dble (Rho (l,l,NCF))
!       TODO: This is not working I need to fix it
        if((dble(Rho(l,l,NCF))<-1.E-8).or.(dble(Rho(l,l,NCF))>1.000001_q).or.(dabs(dimag(Rho(l,l,NCF)))>1E-6)&
                &.or.(sum_rule>1.000001_q))then

            write(unit_error,*) Rho(l,l,NCF),l,sum_rule,Ndim,bias_R*hartree,bias_L*hartree,frequency/(2*pi_d*time_unit)
            write(*,*) 'WEIRD RESULTS! -> CHECK strange_populations.dat but we continue...'
        end if
    enddo
    close (unit_error)
    
    return
end subroutine check_populations


subroutine write_cur_out(current, NCF)
    implicit none
    real(q), intent(in) :: current(:)
    integer, intent(in) :: NCF
    integer :: pn

    open (unit_curr, file='Current_0.dat')
    write (unit_curr, *) 'n Current(pA)'
    write (unit_curr, *) 0, current(NCF+1)*pA
    do pn = 1, NCF
        if ((current(pn+NCF+1)).eq.zero) then
!           write until reaching a vanishing fourier component. 
            exit
        endif
        write (unit_curr, *) pn, current(pn+NCF+1)*pA
        write (unit_curr, *) -pn, current(-pn+NCF+1)*pA
    enddo
    close (unit_curr)

    return
end subroutine write_cur_out


subroutine write_pop_out(Rho, NCF, Ndim)
    implicit none
    integer, intent(in) :: NCF, Ndim
    complex (qc), intent(in) :: Rho(:,:,:)

    integer :: l, pn 
    
    
    write(*,*) 'Populations are written in POPULATIONS.dat'
    write(*,*) ''
    open  (unit_pop, file='POPULATIONS.dat')

    write (unit_pop,*) 'n Rho(1,1,n) Rho(2,2,n) ... Rho(Ndim,Ndim,n)'
    write (unit_pop,*) 0, (dble(Rho (l,l,NCF+1)), l=1, Ndim)
    
    do pn = 0, NCF
        write (unit_pop,*) pn, (dble(Rho (l,l,pn+NCF+1)), l=1, Ndim)
        write (unit_pop,*) -pn, (dble(Rho (l,l,-pn+NCF+1)), l=1, Ndim)
    enddo
    close (unit_pop)

    return
end subroutine write_pop_out


subroutine write_coh_out(Rho, NCF)
    implicit none
    integer, intent(in) :: NCF
    complex(qc), intent(in) :: Rho(:,:,:)
    integer :: pn
    
    write(*,*) 'Coherences are written in COHERENCES_'
    write(*,*) ''
    
    open  (unit_coh, file='COHERENCES.dat')
    write (unit_coh,*) 'u v pn Re(Rho(u,v,pn)) Im(Rho(u,v,pn))'
    
    do u = 1, Ndim
    do v = 1, Ndim
        write (unit_coh,*) u,v,0, dble(Rho (u,v,NCF+1)), dimag(Rho (u,v,NCF+1))

        pnloop: do pn = 1, NCF    
        if ((Rho (u,v,pn+NCF+1)).eq.zero) then
            exit pnloop
        endif
        write (unit_coh,*) u,v,pn, dble(Rho (u,v,pn+NCF+1)), dimag(Rho (u,v,pn+NCF+1))
        write (unit_coh,*) u,v,-pn, dble(Rho (u,v,-pn+NCF+1)), dimag(Rho (u,v,-pn+NCF+1))
        enddo pnloop

    enddo
    enddo
    
    close (unit_coh)
    return
end subroutine write_coh_out


subroutine write_spin_out (Rho, NCF, Nm, Ndim, Ss, spinX, spinY, spinZ, spin2_T, H, hx, hy, hz,&
                           Sx, Sy, Sz, Sh, spin2_ave)
    implicit none
    integer, intent(in) :: NCF, Nm, Ndim
    real(q), intent(in) ::  hx(:), hy(:), hz(:)
    complex(q), intent(in) :: Ss(:,:,:,:), spinX(:,:,:), spinY(:,:,:), spinZ(:,:,:), spin2_T(:,:)
    complex(q), intent(in) :: Rho(:,:,:), H(:,:)
    complex(q), intent(inout) :: spin2_ave
    complex(q), intent(inout), allocatable :: Sx(:), Sy(:), Sz(:), Sh(:)
    
    ! Declare local variables

    integer ::  pn, l
    
    write(*,*) 'Spin Sx,Sy,Sz,Sh per site are written in SpinFloquet_'
    write(*,*) ''
    open (unit_floq, file='SpinFloquet.dat')
    write (unit_floq,*) 'n ((Sx(1) Sy(1) Sz(1) Sh(1)) (Sx(2) Sy(2) Sz(2)) ...) spin2_ave real(sqrt(1+4*(spin2_ave))-1)*0.5'
    
    do pn = -NCF, NCF
               
        call SpinFloquet (Nm, Ndim, Ss, spinX, spinY, spinZ, spin2_T, H, Rho(:,:,pn+NCF+1), hx, hy, hz,&
                  Sx,Sy,Sz,Sh,spin2_ave)
        
        write (unit_floq,*) pn,&
            (real(Sx(l)), real(Sy(l)), real(Sz(l)), real(Sh(l)),  l=1, Nm),&
            real(spin2_ave),real(sqrt(1+4*(spin2_ave))-1)*0.5
        
    enddo
    close (unit_floq)

    return
end subroutine write_spin_out

end module output