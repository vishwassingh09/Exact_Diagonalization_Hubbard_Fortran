! 2D Hubbard Model - Matrix-Free Lanczos Diagonalization (ARPACK)
! Compile: gfortran -O3 hubbard2d_lanczos.f90 -larpack -lopenblas -o hubbard_run

module hubbard_data
    implicit none
    integer :: lx, ly, ns, nup, ndn, ne, ts
    integer :: count_up, count_dn
    integer, allocatable :: valid_up(:), valid_dn(:)
    real*8, allocatable :: uterm(:)
    real*8 :: u, t
contains
    ! Binary search to find the state index rapidly
    integer function find_state(state_array, n, target)
        integer, intent(in) :: state_array(:), n, target
        integer :: left, right, mid
        left = 1
        right = n
        find_state = -1
        do while (left <= right)
            mid = (left + right) / 2
            if (state_array(mid) == target) then
                find_state = mid
                return
            else if (state_array(mid) < target) then
                left = mid + 1
            else
                right = mid - 1
            endif
        enddo
    end function find_state

    ! Population count
    integer function count_bits(state, total_sites)
        integer, intent(in) :: state, total_sites
        integer :: b_idx, temp_c
        temp_c = 0
        do b_idx = 0, total_sites - 1
            if (btest(state, b_idx)) temp_c = temp_c + 1
        enddo
        count_bits = temp_c
    end function count_bits

    ! Jordan-Wigner Phase Sign
    integer function get_phase(state, site1, site2)
        integer, intent(in) :: state, site1, site2
        integer :: start_s, end_s, b_idx, count_b
        start_s = min(site1, site2)
        end_s = max(site1, site2)
        count_b = 0
        do b_idx = start_s + 1, end_s - 1
            if (btest(state, b_idx)) count_b = count_b + 1
        enddo
        if (mod(count_b, 2) /= 0) then
            get_phase = -1
        else
            get_phase = 1
        endif
    end function get_phase

    ! Matrix-vector multiplication: Y = H * X
    subroutine compute_Hx(x_vec, y_vec)
        real*8, intent(in) :: x_vec(ts)
        real*8, intent(out) :: y_vec(ts)
        integer :: k, i, j, x_coord, y_coord, n_neighbor
        integer :: u_idx, d_idx, l, phase
        integer :: c1_up, a1_up, c1_dn, a1_dn, new_u_idx, new_d_idx
        integer :: num_neighbors, neighbors(4)

        y_vec = 0.d0

        !$OMP PARALLEL DO PRIVATE(k, u_idx, d_idx, i, x_coord, y_coord, num_neighbors, neighbors, &
        !$OMP c1_up, a1_up, new_u_idx, phase, l, c1_dn, a1_dn, new_d_idx, n_neighbor, j)
        do k = 1, ts
            ! Apply diagonal U term
            y_vec(k) = y_vec(k) + u * uterm(k) * x_vec(k)

            ! Decode global index 'k' into independent up/down indices
            d_idx = (k - 1) / count_up + 1
            u_idx = mod(k - 1, count_up) + 1

            ! Apply Hopping Terms
            do i = 0, ns - 1
                x_coord = mod(i, lx)
                y_coord = i / lx
                num_neighbors = 0
                if (x_coord + 1 < lx) then 
                    num_neighbors = num_neighbors + 1; neighbors(num_neighbors) = i + 1
                endif
                if (x_coord - 1 >= 0) then 
                    num_neighbors = num_neighbors + 1; neighbors(num_neighbors) = i - 1
                endif
                if (y_coord + 1 < ly) then 
                    num_neighbors = num_neighbors + 1; neighbors(num_neighbors) = i + lx
                endif
                if (y_coord - 1 >= 0) then 
                    num_neighbors = num_neighbors + 1; neighbors(num_neighbors) = i - lx
                endif

                ! Up Spin Hopping
                if (.not. btest(valid_up(u_idx), i)) then
                    c1_up = ibset(valid_up(u_idx), i)
                    do n_neighbor = 1, num_neighbors
                        j = neighbors(n_neighbor)
                        if (btest(c1_up, j)) then
                            a1_up = ibclr(c1_up, j)
                            new_u_idx = find_state(valid_up, count_up, a1_up)
                            if (new_u_idx > 0) then
                                phase = get_phase(valid_up(u_idx), i, j)
                                l = (d_idx - 1) * count_up + new_u_idx
                                y_vec(k) = y_vec(k) + t * dble(phase) * x_vec(l)
                            endif
                        endif
                    enddo
                endif

                ! Down Spin Hopping
                if (.not. btest(valid_dn(d_idx), i)) then
                    c1_dn = ibset(valid_dn(d_idx), i)
                    do n_neighbor = 1, num_neighbors
                        j = neighbors(n_neighbor)
                        if (btest(c1_dn, j)) then
                            a1_dn = ibclr(c1_dn, j)
                            new_d_idx = find_state(valid_dn, count_dn, a1_dn)
                            if (new_d_idx > 0) then
                                phase = get_phase(valid_dn(d_idx), i, j)
                                l = (new_d_idx - 1) * count_up + u_idx
                                y_vec(k) = y_vec(k) + t * dble(phase) * x_vec(l)
                            endif
                        endif
                    enddo
                endif
            enddo
        enddo
        !$OMP END PARALLEL DO
    end subroutine compute_Hx
end module hubbard_data

program hubbard_2d_lanczos
    use hubbard_data
    implicit none

    integer :: i, j, temp
    real :: start_time, end_time

    ! ARPACK Variables
    integer :: ido, nev, ncv, info, lworkl, ierr
    integer :: iparam(11), ipntr(11)
    logical, allocatable :: select_vec(:)
    real*8, allocatable :: resid(:), v(:,:), workd(:), workl(:), d(:,:)
    real*8 :: tol, sigma

    call cpu_time(start_time)

    ! Model Parameters
    lx = 2
    ly = 2
    nup = 2
    ndn = 2
    u = 1.d0
    t = -1.d0
    ns = lx * ly
    ne = nup + ndn

    allocate(valid_up(2**ns), valid_dn(2**ns))
    count_up = 0; count_dn = 0

    do i = 0, (2**ns) - 1
        if (count_bits(i, ns) == nup) then
            count_up = count_up + 1; valid_up(count_up) = i
        endif
    enddo
    do i = 0, (2**ns) - 1
        if (count_bits(i, ns) == ndn) then
            count_dn = count_dn + 1; valid_dn(count_dn) = i
        endif
    enddo

    ts = count_up * count_dn
    write(*,*) "Grid: ", lx, "x", ly, " | Sites =", ns, " | Electrons =", ne, " | States =", ts
    write(*,*) "Parameters: U = ", u, ", t = ", -t

    ! Precompute U interaction term
    allocate(uterm(ts))
    uterm = 0.d0
    do i = 1, ts
        temp = 0
        do j = 0, ns - 1
            if (btest(valid_up(mod(i-1, count_up)+1), j) .and. btest(valid_dn((i-1)/count_up+1), j)) temp = temp + 1
        enddo
        uterm(i) = dble(temp)
    enddo

    ! Initialize ARPACK
    nev = 1          ! Number of eigenvalues requested (Ground State)
    ncv = 10         ! Number of Lanczos vectors generated at each iteration
    tol = 0.d0       ! Machine precision
    lworkl = ncv**2 + 8*ncv
    ido = 0
    info = 0

    allocate(resid(ts), v(ts, ncv), workd(3*ts), workl(lworkl))
    allocate(select_vec(ncv), d(ncv, 2))

    iparam(1) = 1    ! Exact shifts
    iparam(3) = 300  ! Maximum iterations
    iparam(7) = 1    ! Standard eigenvalue problem

    write(*,*) "Starting Matrix-Free Lanczos Diagonalization..."

    ! ARPACK Reverse Communication Loop
    10 continue
    call dsaupd(ido, 'I', ts, 'SA', nev, tol, resid, ncv, v, ts, iparam, ipntr, workd, workl, lworkl, info)

    if (ido == 1 .or. ido == -1) then
        ! ARPACK requests Y = H * X
        call compute_Hx(workd(ipntr(1):), workd(ipntr(2):))
        goto 10
    endif

    if (info < 0) then
        write(*,*) "Error with dsaupd, info = ", info
    else
        ! Extract the eigenvalue
        call dseupd(.true., 'A', select_vec, d, v, ts, sigma, 'I', ts, 'SA', nev, tol, resid, ncv, v, ts, iparam, ipntr, workd, workl, lworkl, ierr)
        
        write(*,*) "----------------------------------------"
        write(*,*) "Ground State Energy: ", d(1, 1)
        write(*,*) "----------------------------------------"
    endif

    call cpu_time(end_time)
    write(*, "(A, F8.4, A)") "Execution Time: ", end_time - start_time, " seconds."
end program hubbard_2d_lanczos
