! 2D Hubbard Model - Unified Basis Generator & Exact Diagonalization
! Compile with LAPACK: gfortran hubbard2d.f90 -llapack -o hubbard2d

program hubbard_2d_unified
    implicit none
    
    integer :: lx, ly, ns, nup, ndn, ne, ts
    integer :: i, j, k, l, inf, l1, idx
    integer :: count_up, count_dn, u_idx, d_idx
    integer, allocatable :: valid_up(:), valid_dn(:)
    
    integer :: x, y, n_neighbor, num_neighbors
    integer :: c1_up, a1_up, c1_dn, a1_dn, phase
    real*8 :: u, t, u1
    real*8, allocatable :: hmatrix(:,:), eig(:), work(:), uterm(:), temp_h(:,:)
    integer, allocatable :: sup(:), sdn(:)
    integer :: neighbors(4)
    integer :: temp
    real :: start_time, end_time

    ! Start timing
    call cpu_time(start_time)

    ! Requested Parameters: 2 Sites (2x1 grid), 2 Electrons (1 Up, 1 Down)
    lx = 2
    ly = 4
    nup = 4
    ndn = 4
    ns = lx * ly
    ne = nup + ndn
    
    ! Maximum sizes for sectors
    allocate(valid_up(2**ns), valid_dn(2**ns))

    count_up = 0
    count_dn = 0

    ! Generate Spin-Up Sector
    do i = 0, (2**ns) - 1
        if (count_bits(i, ns) == nup) then
            count_up = count_up + 1
            valid_up(count_up) = i
        endif
    enddo

    ! Generate Spin-Down Sector
    do i = 0, (2**ns) - 1
        if (count_bits(i, ns) == ndn) then
            count_dn = count_dn + 1
            valid_dn(count_dn) = i
        endif
    enddo

    ts = count_up * count_dn
    write(*,*) "Grid: ", lx, "x", ly, " | Sites =", ns, " | Electrons =", ne, " | States =", ts
    
    ! Set T and U
    u = 1.d0
    t = 1.d0
    write(*,*) "Parameters: U = ", u, ", t = ", t, " (Open Boundary Conditions)"
    
    ! The original code used t = -t internally for standard physical convention
    t = -t
    
    allocate(sup(ts), sdn(ts), hmatrix(ts,ts), eig(ts))
    allocate(work(ts*(3+ts/2)), uterm(ts), temp_h(ts,ts))

    ! Tensor product: Pair every down-state with every up-state in memory
    idx = 1
    do d_idx = 1, count_dn
        do u_idx = 1, count_up
            sdn(idx) = valid_dn(d_idx)
            sup(idx) = valid_up(u_idx)
            idx = idx + 1
        enddo
    enddo

    write(*,*) "---------------processing-------------"
    hmatrix = 0.d0
    uterm = 0.d0

    ! ------------------------------------------------------------------
    ! Calculating U term (Diagonal)
    ! ------------------------------------------------------------------
    do i = 1, ts
        temp = 0
        do j = 0, ns - 1
            if (btest(sup(i), j) .and. btest(sdn(i), j)) then
                temp = temp + 1
            endif
        enddo
        hmatrix(i, i) = temp + hmatrix(i, i)
        uterm(i) = temp + uterm(i)
    enddo

    ! ------------------------------------------------------------------
    ! Calculating t term (Off-Diagonal 2D Hopping)
    ! ------------------------------------------------------------------
    do k = 1, ts
        do i = 0, ns - 1
            x = mod(i, lx)
            y = i / lx
            
            ! Identify valid 2D Open Boundary neighbors for site i
            num_neighbors = 0
            if (x + 1 < lx) then 
                num_neighbors = num_neighbors + 1
                neighbors(num_neighbors) = i + 1      ! Right
            endif
            if (x - 1 >= 0) then 
                num_neighbors = num_neighbors + 1
                neighbors(num_neighbors) = i - 1      ! Left
            endif
            if (y + 1 < ly) then 
                num_neighbors = num_neighbors + 1
                neighbors(num_neighbors) = i + lx      ! Top
            endif
            if (y - 1 >= 0) then 
                num_neighbors = num_neighbors + 1
                neighbors(num_neighbors) = i - lx      ! Bottom
            endif

            ! --------- Up States Hopping -------------------    
            if (.not. btest(sup(k), i)) then    
                c1_up = ibset(sup(k), i)        
                do n_neighbor = 1, num_neighbors
                    j = neighbors(n_neighbor)
                    a1_up = ibclr(c1_up, j)
                    
                    phase = get_phase(sup(k), i, j)
                    
                    do l = 1, ts
                        if (l /= k) then
                            if (a1_up == sup(l) .and. sdn(l) == sdn(k)) then
                                hmatrix(k, l) = t * dble(phase)                
                            endif
                        endif 
                    enddo
                enddo
            endif

            ! --------- Down States Hopping -----------------
            if (.not. btest(sdn(k), i)) then    
                c1_dn = ibset(sdn(k), i)
                do n_neighbor = 1, num_neighbors
                    j = neighbors(n_neighbor)
                    a1_dn = ibclr(c1_dn, j)
                    
                    phase = get_phase(sdn(k), i, j)
                    
                    do l = 1, ts
                        if (l /= k) then
                            if (a1_dn == sdn(l) .and. sup(l) == sup(k)) then
                                hmatrix(k, l) = t * dble(phase)    
                            endif
                        endif
                    enddo
                enddo
            endif
        enddo !i
    enddo !k

    temp_h = hmatrix
    l1 = ts * (3 + ts / 2)
    hmatrix = 0.d0
    u1 = u    

    hmatrix = temp_h
    do i = 1, ts
        hmatrix(i, i) = u1 * uterm(i)
    enddo

    !write(*,*) "The Hamiltonian matrix is:"
    !do j = 1, ts
    !    do k = 1, ts
    !        write(*, "(F8.2)", advance='no') hmatrix(k, j)
    !    enddo
    !    write(*,*) "" 
    !enddo
        
    ! Diagonalizing the matrix
    write(*,*) "Diagonalizing... U =", u1
    call dsyev('V', 'U', ts, hmatrix, ts, eig, work, l1, inf)
    
    write(*,*) "INFO=", inf
    write(*,*) "Eigenvalues:"
    do j = 1, ts    
        write(*, "(F10.4)") eig(j)
    enddo
        
    call cpu_time(end_time)
    write(*, "(A, F8.4, A)") "Execution Time: ", end_time - start_time, " seconds."

contains

    ! Helper function to act as a population count operator
    integer function count_bits(state, total_sites)
        integer, intent(in) :: state, total_sites
        integer :: b_idx, temp_c
        temp_c = 0
        do b_idx = 0, total_sites - 1
            if (btest(state, b_idx)) then
                temp_c = temp_c + 1
            endif
        enddo
        count_bits = temp_c
    end function count_bits

    ! Calculates the Fermionic Jordan-Wigner Phase Sign
    integer function get_phase(state, site1, site2)
        integer, intent(in) :: state, site1, site2
        integer :: start_s, end_s, b_idx, count_b
        
        start_s = min(site1, site2)
        end_s = max(site1, site2)
        count_b = 0
        
        ! Count occupied states strictly between the two hopping sites
        do b_idx = start_s + 1, end_s - 1
            if (btest(state, b_idx)) count_b = count_b + 1
        enddo
        
        if (mod(count_b, 2) /= 0) then
            get_phase = -1
        else
            get_phase = 1
        endif
    end function get_phase

end program hubbard_2d_unified

! =====================================================================
! EXPECTED OUTPUT (Hardcoded for t=1, U=1, 2 Sites, 2 Electrons (1 Up, 1 Dn)):
! =====================================================================
! Grid:  2 x 1 | Sites = 2 | Electrons = 2 | States = 4
! Parameters: U =  1.0000000000000000 , t =  1.0000000000000000 (Open Boundary Conditions)
! ---------------processing-------------
! The Hamiltonian matrix is:
!     1.00   -1.00   -1.00    0.00
!    -1.00    0.00    0.00   -1.00
!    -1.00    0.00    0.00   -1.00
!     0.00   -1.00   -1.00    1.00
! Diagonalizing... U =  1.0000000000000000
! INFO=           0
! Eigenvalues:
!    -1.5616
!     0.0000
!     1.0000
!     2.5616
! Execution Time:   0.0015 seconds.
! =====================================================================
