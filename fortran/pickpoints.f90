module pickpoints
    implicit none

    contains

    subroutine pick_points(z,f,nmin,M,zp,fp)
        implicit none
        complex(kind=16),intent(in) :: z(:),f(:)            ! Input points and values 
        integer,intent(in) :: nmin                                     ! lowest index to use from the input points. nmin < 0, mirror values are added.
        integer,intent(in) :: M                                        ! number of points to pick.
        complex(kind=16),allocatable :: zp(:),fp(:)         ! Picked points and values to use in the analytical continuation
        integer :: nadd
       
        ! deallocate zp and fp if they are already allocated
        if(allocated(zp)) deallocate(zp)
        if(allocated(fp)) deallocate(fp)
        allocate(zp(M))
        allocate(fp(M))

        if (nmin >= 0) then
            zp = z(nmin+1:nmin+M)
            fp = f(nmin+1:nmin+M)
        else
            nadd = -nmin
            zp(1:nadd) = conjg(z(nadd:1:-1)) 
            zp(nadd+1:) = z(1:M-nadd)
            fp(1:nadd) = conjg(f(nadd:1:-1)) 
            fp(nadd+1:) = f(1:M-nadd)
        endif
    end subroutine

end module
