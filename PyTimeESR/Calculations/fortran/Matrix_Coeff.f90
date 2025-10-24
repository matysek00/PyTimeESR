module Matrix_Coeff
    Use declarations
    CONTAINS

subroutine coef_matrix_transport (Ndim, NF, NCF, GA, QME_Tensor_A)
    implicit none
    integer, intent (in) :: Ndim, NF, NCF
    complex (qc), intent (in) :: GA (Ndim, Ndim, Ndim, Ndim, NF)

    !complex (qc), intent(out), allocatable :: QME_Tensor_A (:,:,:,:,:,:)
    complex (qc), intent(out) :: QME_Tensor_A (Ndim, Ndim, NF, Ndim, Ndim, NF)
    integer :: l, j, v, u, n, m 
            
    !allocate (QME_Tensor_A (Ndim, Ndim, NF, Ndim, Ndim, NF))
    !allocate (QME_Tensor_A (NF, Ndim, Ndim, NF, Ndim, Ndim))
    QME_Tensor_A = zero
    
    fourier_n: do n= 1, NF
    fourier_m: do m = max(1,n-NCF),min(NF,n+NCF)
    level_l: do l=1, Ndim
    level_j: do j=1, Ndim
    level_v: do v=1, Ndim
    level_u: do u=1, Ndim
!       Gamma_{j,v,v,u:n-m} rho_{l,u,m}
        QME_Tensor_A (j,l,n, u,l,m) = QME_Tensor_A (j,l,n, u,l,m) + ui*GA (j,v,v,u,n-m+NCF+1)

!       Gamma_{l,v,v,u:m-n} rho_{u,j,m}
        QME_Tensor_A (j,l,n, j,u,m) = QME_Tensor_A (j,l,n, j,u,m) + ui*CONJG(GA (l,v,v,u,m-n+NCF+1))

!       Gamma_{u,j,l,v:m-n} rho_{v,u,m}
        QME_Tensor_A (j,l,n, u,v,m) = QME_Tensor_A (j,l,n, u,v,m) - ui*CONJG(GA (u,j,l,v,m-n+NCF+1))
           
!       Gamma_{v,l,j,u:n-m} rho_{v,u,m}
        QME_Tensor_A (j,l,n, u,v,m) = QME_Tensor_A (j,l,n, u,v,m) - ui*(GA (v,l,j,u,n-m+NCF+1))
            
    enddo level_u
    enddo level_v
    enddo level_j
    enddo level_l
    enddo fourier_m
    enddo fourier_n
    
    return 
end subroutine coef_matrix_transport


subroutine coeff_matrix (Ndim, frequency, NF, G, Rho)
    implicit none

    !F2PY INTENT(HIDE) :: Ndim 
    !F2PY INTENT(HIDE) :: Ndim 
    !F2PY INTENT(HIDE) :: NF
    !F2PY COMPLEX :: G(Ndim, Ndim, Ndim, Ndim, NF, 2)
    !F2PY INTENT(OUT) :: Rho(Ndim, Ndim, NF)

    integer,intent (in) :: Ndim, NF
    real (q), intent(in) :: frequency
    complex (qc), intent (in) :: G (Ndim, Ndim, Ndim, Ndim, NF, 2)
    complex (qc), intent (out) :: Rho (Ndim, Ndim, NF)
    
    !complex (qc), allocatable :: QME_Tensor_R (:,:,:,:,:,:), QME_Tensor_L (:,:,:,:,:,:)
    complex (qc) :: QME_Tensor_R (Ndim,Ndim,NF,Ndim,Ndim,NF), QME_Tensor_L (Ndim,Ndim,NF,Ndim,Ndim,NF)
    complex (qc) :: QME_Tensor (Ndim,Ndim,NF,Ndim,Ndim,NF)
    complex (qc), allocatable :: A (:,:), B(:)

    integer :: l, j, n, i2, Nmatrix, NCF
    integer, allocatable:: IPIV(:)

!   Solve The Quantum Master Equation writen as
!       rho_{l,j,n} A_{l,j,n u,v,m}  = delta(v,1) * delta(u,1) * delta (n,0)
!       B = delta(v,1) * delta(u,1) * delta (n,0)
!       B is the right hand side of the equation 
    
    NCF = (NF-1)/2
    Nmatrix=Ndim*Ndim*NF

    !allocate (QME_Tensor (Ndim,Ndim,NF,Ndim,Ndim,NF))
    allocate (IPIV(Nmatrix))
    allocate (A(Nmatrix, Nmatrix), B(Nmatrix))
    
!   In dot zeroth-order dynamics 
    QME_Tensor = zero
    do n=1,NF
    do l = 1, Ndim
    do j = 1, Ndim
        QME_Tensor (j,l,n, j,l,n) = Delta (l,j) + frequency*(n-NCF-1)
    enddo
    enddo
    enddo

!   Right and left electrode hopping
    call coef_matrix_transport (Ndim, NF, NCF, G(:,:,:,:,:,1), QME_Tensor_R) 
    call coef_matrix_transport (Ndim, NF, NCF, G(:,:,:,:,:,2), QME_Tensor_L) 
    
    QME_Tensor = QME_Tensor + QME_Tensor_L + QME_Tensor_R      

!   I don't understand why we do this
    do n=1,NF
        QME_Tensor(1,1,n, :,:,:) = zero
        do l = 1, Ndim
            QME_Tensor(1,1,n, l,l,n) = one!        t0(:)=t0(:)/time_unit; t1(:)=t1(:)/time_unit
            !        t_initial=t_initial/time_unit; t_final=t_final/time_unit
        enddo  
    enddo
    
    A = reshape(QME_Tensor, (/Nmatrix, Nmatrix/))
    i2 = 0
    
!   X definition entails detailed balance
    B = zero
    B (1+ndim*ndim*(NCF)) = one
    
!   Solve the equations A*X = B. X retuned in B 
!   A  is n x m matrix
!   B, X are  l x k matrix
!   zgsev (n, l, A, m, IPIV, B, l, INFO)
    call zgesv(Nmatrix, 1, A, Nmatrix, IPIV, B, Nmatrix, INFO)
    
    if (INFO > 0) then
        print *, 'The linear problem is singular and cannot be solved, check input: dimension and bias. INFO=', INFO
!       i, U(i,i) computed in DOUBLE PRECISION is
!       exactly zero.  The factorization has been completed,
!       but the factor U is exactly singular, so the solution
!       could not be computed.
    else if (INFO < 0) then
        print *, 'if INFO = -i, the i-th argument had an illegal value. INFO=', INFO
    endif

    Rho = reshape(B, (/Ndim, Ndim, NF/))
    do n = 1, NF
!   Make the order of indicies consistent with the documentation
        Rho (:,:, n) = transpose(Rho (:,:, n))
    enddo
    
    !deallocate (QME_Tensor_R, QME_Tensor_L)
    deallocate (A, B)
    return
end subroutine coeff_matrix 
end module Matrix_Coeff