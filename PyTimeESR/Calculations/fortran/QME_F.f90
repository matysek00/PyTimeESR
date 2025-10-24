module QME_F
Use declarations
CONTAINS
!
! Rate computation
!
   subroutine rates (Ndim, orb, frequency, gamma_0, lambda, Spin_polarization,&
        NF, p_max, Adrive, Bdrive, GammaC, bias, Delta, Cutoff, Temperature,& 
        N_int, transport_exponent, GA, GCA)
     implicit none
        !F2PY INTENT(HIDE) :: Ndim
        !F2PY INTENT(HIDE) :: orb
        !F2PY COMPLEX :: lambda(Ndim, Ndim, orb)
        !F2PY COMPLEX :: Delta(Ndim, Ndim)
        !F2PY INTENT(OUT) :: GA(Ndim, Ndim, Ndim, Ndim, NF)
        !F2PY INTENT(OUT) :: GCA(Ndim, Ndim, Ndim, Ndim, NF)
     integer :: Ndim, NF, p_max, orb, N_int
     complex (qc), intent (in), dimension(Ndim,Ndim,orb) :: lambda 
     real (q), intent (in), dimension(Ndim, Ndim) :: Delta
     real (q), intent (in) :: gamma_0, Temperature, frequency, Spin_polarization
     real (q), intent (in) :: Cutoff, bias, Bdrive, transport_exponent
     complex (qc), intent (in) :: GammaC, Adrive

     !complex (qc), intent (out), dimension(Ndim, Ndim, Ndim, Ndim, 2*NCF+1) :: GA, GCA 
     complex (qc), intent (out) :: GA (Ndim,Ndim,Ndim,Ndim,NF), GCA (Ndim,Ndim,Ndim,Ndim,NF)
!       G_Alpha and GC_Alpha 
     complex (qc), dimension(2*p_max-1) :: fermi, ufermi     
     complex (q), dimension(2*p_max-1) :: Kbess
     complex (qc) :: g_up, g_dn
     complex (qc), dimension(Ndim, Ndim) :: Lvluj, Ljulv
     complex (qc) :: bessel_contribution, ubessel_contribution

     integer :: j, u, nfour, n_index
        
        g_up = 0.5*gamma_0*(1+Spin_polarization)
        g_dn = 0.5*gamma_0*(1-Spin_polarization)

!       Calculate Contribution of Bessel functions
!       K(p) = J(p) + .5*A*(J(p-1)+ J(p+1))
        print *, Bdrive/frequency, p_max
        Kbess = Bessel_K(Bdrive/frequency, Adrive, p_max)

!       loop on states:
        level_j: do j=1,Ndim
        level_u: do u=1,Ndim

!       Skip lambda is zero
        if((lambda(u,j,1).eq.zero).and.(lambda(u,j,2).eq.zero).and.&
                &((lambda(j,u,1).eq.zero).and.(lambda(j,u,2).eq.zero)))then
                GA (:,:, j,u,:) = zero
                GCA (:,:, j,u,:) = zero
                cycle 
                ! TODO: should this include other spin channels?
!               Let the record show that this cycle avoided a 400 line if statement 
!               when it was originally added.
        end if

!       Green's functions integrals involving occupation factors
        call ExtendedFermiIntegral(Delta(j,u), frequency, bias, p_max-1, Temperature, Cutoff, &
                        GammaC, N_int, transport_exponent, fermi, ufermi)

        call Orbital_overlaps(lambda, j, u, orb, g_up, g_dn, Ndim, Lvluj, Ljulv)
        fourier_component: do n_index = 1, NF 
                nfour = n_index - (NF-1)/2 - 1

!               contribution of bessel functions
!               bessel_cont  = sum_p K*_{p-n} K_p  I(p)                              
!               ubessel_cont = sum_p K*_p K_{p+n}  uI(p) 
                call compute_bessel_contribution(Kbess, fermi, ufermi, p_max, nfour, &
                              bessel_contribution, ubessel_contribution)
                
                ! TODO: ADD PHASE
                GA  (:,:,j,u,n_index) = 0.5*(Lvluj*bessel_contribution + Ljulv*ubessel_contribution)
                GCA (:,:,j,u,n_index) = 0.5*(Lvluj*bessel_contribution - Ljulv*ubessel_contribution)
                
        enddo fourier_component
        enddo level_u
        enddo level_j
      return
      
   end subroutine rates 

!
!       Fermi occupation function with limitors to avoid underflow/overflows
!
        function FermiDist (e) result (eFermi)
        implicit none
        real (q) :: eFermi
        real (q), intent(in) :: e
        
        if ( e > 0._q ) then
                eFermi = exp(-e)/(exp(-e)+1._q)
             else
               eFermi = 1._q/(exp(e)+1._q)
        endif
      end function FermiDist

      subroutine ExtendedFermiIntegral ( D, frequency, V, p_max, T, Cutoff, GammaC, N, transport_exponent,  fermiA, ufermiA)
        implicit none
        real (q), intent(in) :: D, V, T, Cutoff, frequency, transport_exponent
        complex (qc), intent(in) :: GammaC
        real (q) :: e, step_e
        integer :: i, N, p, p_max, p_ind
        complex (qc), dimension(2*p_max+1):: fermiA, ufermiA
        real (q), dimension(N) :: f, te
        complex (qc) :: denom, udenom
        ! I11_p = i/pi \int_{-Cutoff}^{Cutoff} dE t(E)*f(E,V)/(E-D+p*frequency+ui\Gamma_C)
        ! I21_p = -i/pi \int_{-Cutoff}^{Cutoff} dE t(E)*(1-f(E,V))/(E+D+p*frequency-ui\Gamma_C)
        ! f(E,V) = \frac{1}{\exp(\beta (E-V)) + 1} Fermi distribution with V as fermi level
        ! t(E) = exp(cE) WKB aproximation to first order 
        
        ! TODO: the normal integrals could be included in this by setting p_max = 0
        
        ! calculate all fermi contributions they are the same for all integrals
        step_e = 2*Cutoff/(N-1)/T ! rescaling with to units T = 1 
        e= -(Cutoff + V)/T
        
        fstep: do i = 1, N
             e = e + step_e
             f(i) = FermiDist (e) 
             te(i) = exp (transport_exponent*e*T/2)
        enddo fstep

        step_e = 2*Cutoff/(N-1) ! now in atomic units
        
        ploop: do p = -p_max, p_max
             p_ind = p+p_max+1
             
             ! The constant part of the denominator 
             denom =  - D + p*frequency + ui*GammaC
             udenom = D + p*frequency - ui*GammaC
   
             e = - Cutoff
             fermiA(p_ind) = 0.5*f(1)/(e+denom)
             ufermiA(p_ind) = 0.5*(1-f(1))/(e+udenom)   
             
             istep: do i = 2, N-1
                  e= e + step_e
                  fermiA(p_ind)=fermiA(p_ind) + te(i)*f(i)/(e+denom)
                  ufermiA(p_ind)=ufermiA(p_ind) + te(i)*(1-f(i))/(e+udenom)
             enddo istep
             
             e = Cutoff
             fermiA(p_ind) = fermiA(p_ind) + 0.5*f(N)/(e+denom)
             ufermiA(p_ind) = ufermiA(p_ind) + 0.5*(1-f(N))/(e+udenom)
   
             fermiA(p_ind) = step_e*ui*fermiA(p_ind)/pi_d
             ufermiA(p_ind) = -step_e*ui*ufermiA(p_ind)/pi_d
   
        enddo ploop
        return
        end subroutine ExtendedFermiIntegral

      
        function Bessel_K(z, Adrive, p_max) result (Kbess)
        implicit none
        real (q), intent (in) :: z
        complex(qc), intent(in) :: Adrive
        integer, intent (in) :: p_max
        complex (qc), dimension(2*p_max-1) :: Kbess
        complex(qc), dimension(2*p_max+1) :: Jbess
        integer :: p
   
                Jbess(p_max+1:) = Bessel_JN(0, p_max, z)
                print *, z, Jbess(p_max+1)
               

                ! J(-p) = (-1)**p J(p) for p = 0,1,2,...
                negative_bessel : do p = 0, p_max -1
                        Jbess(p+1) = ((-1)**(p_max-p))*Jbess(2*p_max+1-p)
                enddo negative_bessel
                
                
                print *, '\sum J_p = ',  sum(Jbess)
                if (sum(Jbess) == 0) then
                        print *, 1 
                        stop 'Error: p_max is to large. Use a smaller value'
                print *, 'if the sum is far bellow 1 use a larger p_max'
                endif

                ! K(p) = J(p) + .5*A*(J(p-1)+ J(p+1))
                Kbess = Jbess(2:2*p_max) + 0.5 * Adrive * (Jbess(1:2*p_max-1) + Jbess(3:2+p_max+1))
                
        return
        end function Bessel_K

        subroutine Orbital_overlaps (lambda, j, u, orb, g_up, g_dn, Ndim, overlapvluj, overlapjulv)
        implicit none
        !F2PY INTENT(HIDE) :: Ndim
        !F2PY INTENT(HIDE) :: orb
        !F2PY INTENT(IN) :: lambda(Ndim, Ndim, orb)
        !F2PY INTENT(IN) :: g_up
        !F2PY INTENT(IN) :: g_dn
        !F2PY INTENT(OUT) :: overlapvluj(Ndim, Ndim)
        !F2PY INTENT(OUT) :: overlapjulv(Ndim, Ndim)
        
        integer, intent(in):: Ndim, j, u, orb
        complex (qc), dimension(:, :, :), intent (in) :: lambda 
        complex (qc), intent (in) :: g_up, g_dn
        complex (qc), dimension(Ndim, Ndim), intent (out) :: overlapvluj, overlapjulv
        integer :: v, l, lorb 
                
                overlapvluj = zero
                overlapjulv = zero
                
                level_v: do v=1,Ndim
                level_l: do l=1, Ndim
                        
                        orbital: do lorb = 1, orb, 2 
                        overlapvluj(v,l) = lambda(v,l,lorb)*conjg(lambda(u,j,lorb))*g_up+&
                                lambda (v,l,lorb+1)*conjg(lambda(u,j,lorb+1))*g_dn
                        overlapjulv(v,l) = lambda(j,u,lorb)*conjg(lambda(l,v,lorb))*g_up+&
                                lambda (j,u,lorb+1)*conjg(lambda(l,v,lorb+1))*g_dn
                        enddo orbital

                enddo level_l
                enddo level_v
   
        return
        end subroutine Orbital_overlaps


        subroutine compute_bessel_contribution(K, fermi, ufermi, p_max, nfour, result_bessel, result_ubessel)
        integer, intent(in) :: p_max, nfour
        complex (qc), intent(in), dimension(2*p_max-1) :: K
        complex (qc), intent(in), dimension(2*p_max-1) :: fermi, ufermi
        complex (qc), intent(out) :: result_bessel, result_ubessel
   
             integer :: p
             result_bessel = 0.0_qc 
             result_ubessel = 0.0_qc
             
             ! sum_p K*_{p-n} K_p  I(p)
             ! sum_p K*_p K_{p+n}  uI(p)
             ! assumes that for |p|>p_max K_p = 0 
             ! for n<0 goes from p=1-p_max to p=p_max-1-n
             ! for n>0 goes from p=1-p_max to p=p_max-1
             ! thus p, p+n, p-n are all in the range 1-p_max to p_max-1
             bessel: do p = max(1, 1-nfour), min(2*p_max-1, 2*p_max-nfour-1)
                  result_bessel  = result_bessel  + conjg(K(p)) * K(p+nfour) * fermi(p+nfour)
                  result_ubessel = result_ubessel + conjg(K(p)) * K(p+nfour) * ufermi(p)
             enddo bessel
             
             return
        end subroutine compute_bessel_contribution

end module QME_F