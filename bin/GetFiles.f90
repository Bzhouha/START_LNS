program read_write_data 
    use iso_fortran_env, only: REAL64
    implicit none
    integer :: i,j,k,l
    integer :: in,jn,kn,ln
    real(REAL64), dimension(:, :, :), allocatable :: xx, yy, zz
    real(REAL64), dimension(:, :, :, :), allocatable :: qq
    real(REAL64), dimension(5) :: tmp

    write(*,*) "输入流场大小:"
    read(*,*) in,jn,kn

    ln=5

    call random_seed()

    allocate(xx(in,jn,kn),yy(in,jn,kn),zz(in,jn,kn))
    allocate(qq(5,in,jn,kn))
    
    do k=1,kn 
        do j=1,jn 
            do i=1,in 
                call random_number(qq(:,i,j,k))
            enddo
        enddo
    enddo

    do k=1,kn
        do j=1,jn 
            do i=1,in
                xx(i,j,k) = i
                yy(i,j,k) = j
                zz(i,j,k) = k
            enddo
        enddo
    enddo
    select case (kn)
    case(1)
        open(23,file="..//files//in//grid.dat",action='write',status='replace',form='unformatted')
        write(23) "x,y"
        write(23) in,jn
        write(23) xx,yy
        close(23)
        ln=5
        open(23,file="..//files//in//flow.dat",action='write',status='replace',form='unformatted')
        write(23) "rho,u,v,w,T"
        write(23) in,jn,ln
        write(23) qq
        close(23)
    case default
        open(21,file="..//files//in//grid.dat",action='write',status='replace',form='unformatted')
        write(21) "x,y,z"
        write(21) in,jn,kn
        write(21) xx,yy,zz
        close(21)
        ln=5
        open(22,file="..//files//in//flow.dat",action='write',status='replace',form='unformatted')
        write(22) "rho,u,v,w,T"
        write(22) in,jn,kn,ln
        write(22) qq
        close(22)
    end select

end program read_write_data
