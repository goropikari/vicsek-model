! To compile, type following command.
! gfortran vicsek.f95 gettime.c -o vicsek && ./vicsek
program main
implicit none
character filename*128
integer:: i,j, iter
integer:: iflag=0, gettime
real:: rand
integer:: num_particle, steps
parameter(num_particle=300)
real*8:: pi=atan(1.0)*4.0
real*8:: x(num_particle), y(num_particle), vx(num_particle), vy(num_particle), theta(num_particle)
real*8:: mean_ang(num_particle)=0, D(num_particle, num_particle), randarray(num_particle)
real*8:: L, vc, dt, eta, r

call system('mkdir -p data')

!L:             system size
!num_particle:  the number of particles
!vc:            constant speed
!dt:            time step
!steps:         the number of run steps
!eta:           Order
!r:             interaction radius
steps=500
L=5.0; vc=1.0; dt=0.1; eta=0.1; r=0.5
! initial condition
i=rand(gettime())
do i=1,num_particle
    theta(i) = unirand(-pi,pi)
    x(i) = unirand(0.0d0,L)
    y(i) = unirand(0.0d0,L)
    vx(i) = vc * cos(theta(i))
    vy(i) = vc * sin(theta(i))
end do

! start vicsek model
do iter=1,steps
    call filewrite(x,y,vx,vy,iter,num_particle) 
    call ppdist(x,y,D,num_particle,L) 
    call cal_mean_ang(D, mean_ang, theta, num_particle, r) 

    
    !i=rand(gettime())
    do i=1,num_particle
        randarray(i) = unirand(-pi,pi)
    end do
    theta = mean_ang + eta * randarray
    vx = vc * cos(theta); vy = vc * sin(theta)
    x = x + vx * dt; y = y + vy * dt

    ! periodic boundary condition
    x = dmod((x + L), L); y = dmod((y + L), L)
end do


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! function, subroutine
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
    ! uniform random number
    function unirand(a,b)
        real*8:: unirand, a, b
!       real:: rand
!       integer:: i
        unirand = a + (b-a)*rand(iflag)
    end function unirand

    ! pairwise distance
    function dist(x1, x2, y1, y2, L)
        real*8 dist, x1, x2, y1, y2, L, dx, dy
        dx = dabs(x1 - x2)
        if (dx > L-dx) then
            dx = L-dx
        end if

        dy = dabs(y1 - y2)
        if (dy > L-dy) then
            dy = L-dy
        end if  
        dist = sqrt(dx**2 + dy**2)
    end function dist

    subroutine ppdist(x,y,D,num_particle,L)
        integer:: num_particle, i, j
        real*8:: x(num_particle), y(num_particle), D(num_particle,num_particle), L

        do i=1,num_particle
            do j=1,num_particle
                D(i,j) = dist(x(i),x(j), y(i), y(j), L)
            end do
        end do
    return
    end subroutine ppdist

    ! mean angle
    subroutine cal_mean_ang(D, mean_ang, theta, N, r)
        integer:: N, near_count, i, j
        real*8:: D(N,N), mean_ang(N), theta(N), list_ang(N), mean_sin, mean_cos, r

        do i=1,N
            near_count = 0
            do j=1,N
                if (D(i,j) .le. r) then
                    list_ang(near_count+1) = theta(j)
                    near_count = near_count + 1
                end if
            end do

            mean_sin=0.0d0; mean_cos=0.0d0
            do j=1,near_count
                mean_sin = mean_sin + sin(list_ang(j))
                mean_cos = mean_cos + cos(list_ang(j))
            end do
            mean_sin = mean_sin / dble(near_count)
            mean_cos = mean_cos / dble(near_count)

            mean_ang(i) = atan2(mean_sin, mean_cos)
        end do
    end subroutine cal_mean_ang

    ! file write
    subroutine filewrite(x,y,vx,vy,iter,N)
        integer:: iter, N
        real*8:: x(N), y(N), vx(N), vy(N)
        write (filename, '("./data/",i5.5, ".csv")') iter-1
            open(17,file=filename, status='replace')
                do i=1,N
                    write (17,*) x(i), ",", y(i), ",", vx(i), ",", vy(i)
                end do
            close(17)
    end subroutine filewrite
end program main
