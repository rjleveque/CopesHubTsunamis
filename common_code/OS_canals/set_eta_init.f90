
subroutine set_eta_init(mbc,mx,my,xlow,ylow,dx,dy,t,veta)

    ! This routine is called only if variable_eta_init = .true.
    ! The input is a single grid patch at time t, as specified by
    !     mbc,mx,my,xlow,ylow,dx,dy,t
    ! The output is an array 
    !     veta(1-mbc:mx+mbc,1-mbc:my+mbc)
    ! with the desired initial surface elevation at all points on this patch.

    ! This is called by qinit and also by filpatch and filval when refining.

    ! This default version sets veta to the value sea_level everywhere and
    ! then adjusts this based on any dtopo deformation that has occured
    ! up to this time t, as determined by interpolation from the 
    ! (possibly time-dependent) ! dtopo files.


    use topo_module
    use geoclaw_module, only: sea_level

    implicit none

    ! Arguments
    integer, intent(in) :: mbc,mx,my
    real(kind=8), intent(in) :: xlow,ylow,dx,dy,t
    real(kind=8), intent(inout) :: veta(1-mbc:mx+mbc,1-mbc:my+mbc)
    
    ! Local
    integer :: i,j,i1,i2,j1,j2,idtopo,jdtopo,kdtopo, &
               index0_dtopowork,ij,m,k
    real(kind=8) :: x,y
    real(kind=8) :: xofyl, xofyu

    ! parameters defining ruled rectangle containing lake / canals:
    integer, parameter :: mbdry = 8
    real(kind=8), parameter,dimension(mbdry) :: &
        ys = (/46.9472,46.9519,46.952,46.9752,46.9914,47.005648,47.012,47.0206/)
    real(kind=8), parameter,dimension(mbdry) :: &
        xl = (/-124.141,-124.168,-124.168,-124.168,-124.168,-124.168,-124.168,-124.168/)
    real(kind=8), parameter,dimension(mbdry) :: &
        xu = (/-124.1323,-124.1323,-124.1323,-124.1418,-124.1418,-124.1483,-124.1483,-124.1537/)
    
    veta = sea_level  ! initialize to sea_level, update below at some (i,j)

    do i=1,mx
        x = xlow + (i-0.5d0)*dx
        do j=1,my
            y = ylow + (j-0.5d0)*dy
            do k=1,mbdry-1
                if ((y >= ys(k)) .and. (y < ys(k+1))) then
                    xofyl = xl(k) + (xl(k+1)-xl(k)) * &
                            (y-ys(k))/((ys(k+1)-ys(k)))
                    xofyu = xu(k) + (xu(k+1)-xu(k)) * &
                            (y-ys(k))/((ys(k+1)-ys(k)))
                    if ((x > xofyl) .and. (x < xofyu)) then
                        veta(i,j) = 2.d0  ! initial level of canals
                        endif
                    exit
                    endif
                enddo
            enddo
        enddo
                    


    if (num_dtopo == 0) return

    do m=1,num_dtopo
        if (t < t0dtopo(m)) cycle  ! this dtopo isn't moving yet
        
        ! compute indices in patch that overlap with this dtopo,
        ! cycling to the next dtopo array if we discover there's no overlap:
        
        i1 = max(int(floor((xlowdtopo(m) - xlow)/dx)), 1)
        if (i1 >= mx) cycle
        
        i2 = min(int(floor((xhidtopo(m) - xlow)/dx)), mx)
        if (i2 < 1) cycle

        j1 = max(int(floor((ylowdtopo(m) - ylow)/dy)), 1)
        if (j1 >= my) cycle
        
        j2 = min(int(floor((yhidtopo(m) - ylow)/dy)), my)
        if (j2 < 1) cycle
    
        ! There is some overlap of dtopo with this patch
        ! Next figure out index into time-dependent dtopo based on t:

        if (mtdtopo(m) == 1) then
            ! Special case: instantaneous displacement at one instant in time
            kdtopo = 1
          else
            kdtopo = int(floor((t-t0dtopo(m))/dtdtopo(m)))+1
            kdtopo = min(kdtopo,mtdtopo(m))
            kdtopo = max(kdtopo,1)
          endif

        index0_dtopowork = i0dtopo(m) + (kdtopo-1)*mxdtopo(m)*mydtopo(m)
        !write(6,*) '+++ index0_dtopowork = ',index0_dtopowork

        ! Adjust eta_init by dtopo on part of patch that overlaps dtopo.
        ! Code below assumes dtopo is smooth enough that we can just 
        ! evaluate at one point in space and time, not doing interpolation.
        ! This could be improved as in topo_update.f, but probably not
        ! necessary (??).

        do i=i1,i2
            x = xlow + (i-0.5d0)*dx
            do j=j1,j2
                y = ylow + (j-0.5d0)*dy
                idtopo = int(floor((x-xlowdtopo(m))/dxdtopo(m))) + 1
                idtopo = max(1, min(mxdtopo(m)-1, idtopo))
                jdtopo = int(floor((yhidtopo(m)-y)/dydtopo(m))) + 1
                jdtopo = max(1, min(mydtopo(m)-1, jdtopo))
                ij = index0_dtopowork + (jdtopo-1)*mxdtopo(m) + idtopo-1

                veta(i,j) = veta(i,j) + dtopowork(ij)

                enddo
            enddo

        enddo
    return

end subroutine set_eta_init
