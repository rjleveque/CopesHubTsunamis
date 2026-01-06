program computeB

    use amr_module, only: dbugunit, parmunit, outunit, inunit, matunit
    use amr_module, only: mxnest, rinfinity, iinfinity
    use amr_module, only: xupper, yupper, xlower, ylower
    use amr_module, only: hxposs, hyposs, intratx, intraty, kratio
    use amr_module, only: cfl, cflv1, cflmax, evol

    use amr_module, only: checkpt_style, checkpt_interval, tchk, nchkpt
    use amr_module, only: rstfile

    use amr_module, only: max1d, maxvar, maxlv

    use amr_module, only: method, mthlim, use_fwaves, numgrids
    use amr_module, only: nghost, mwaves, mcapa, auxtype
    use amr_module, only: tol, tolsp, flag_richardson, flag_gradient

    use amr_module, only: nghost, mthbc
    use amr_module, only: xperdom, yperdom, spheredom

    use amr_module, only: nstop, nout, iout, tfinal, tout, output_style
    use amr_module, only: output_format, printout, verbosity_regrid
    use amr_module, only: output_q_components, output_aux_components
    use amr_module, only: output_aux_onlyonce, matlabu

    use amr_module, only: lfine, lentot, iregridcount, avenumgrids
    use amr_module, only: tvoll, tvollCPU, rvoll, rvol, mstart, possk, ibuff
    use amr_module, only: timeRegridding,timeUpdating, timeValout
    use amr_module, only: timeBound,timeStepgrid, timeFlagger,timeBufnst
    use amr_module, only: timeBoundCPU,timeStepGridCPU,timeRegriddingCPU
    use amr_module, only: timeValoutCPU,timeTick,timeTickCPU
    use amr_module, only: kcheck, iorder, lendim, lenmax, memsize

    use amr_module, only: dprint, eprint, edebug, gprint, nprint, pprint
    use amr_module, only: rprint, sprint, tprint, uprint

    use amr_module, only: t0, tstart_thisrun
    use amr_module, only: tick_clock_start, tick_cpu_start


    ! Data modules
    use geoclaw_module, only: set_geo
    use topo_module, only: read_topo_settings, read_dtopo_settings
    use qinit_module, only: set_qinit
    use refinement_module, only: set_refinement
    use storm_module, only: set_storm
    use friction_module, only: setup_variable_friction
    use gauges_module, only: set_gauges, num_gauges, gauges
    use regions_module, only: set_regions
    use fgout_module, only: set_fgout
    use fgmax_module, only: set_fgmax, FG_num_fgrids
    use multilayer_module, only: set_multilayer
    use adjoint_module, only: read_adjoint_data

    use geoclaw_module, only: coordinate_system, earth_radius, deg2rad


    implicit none

    ! Local variables
    integer :: i, iaux, mw, level
    integer :: ndim, nvar, naux, mcapa1, mindim, dimensional_split
    integer :: nstart, nsteps, nv1, nx, ny, lentotsave, num_gauge_SAVE
    integer :: omp_get_max_threads, maxthreads
    real(kind=8) :: time, ratmet, cut, dtinit, dt_max
    logical :: vtime, rest, output_t0    

    ! Timing variables
    integer(kind=8) :: clock_start, clock_finish, clock_rate, ttotal, count_max
    real(kind=8) :: ttotalcpu
    integer(kind=8) :: tick_clock_finish
    integer, parameter :: timing_unit = 48
    character(len=512) :: timing_line, timing_substr
    character(len=*), parameter :: timing_base_name = "timing."
    character(len=*), parameter :: timing_header_format =                      &
                                                  "(' wall time (', i2,')," // &
                                                  " CPU time (', i2,'), "   // &
                                                  "cells updated (', i2,'),')"

    integer :: i1,j1,j,ii,jj,iindex,jindex
    real(kind=8) :: x,y,xm,xp,ym,yp,topo_integral,dx,dy,aa,xcenter,ycenter
    real(kind=8) :: xcent,ycent,xoff,yoff,mod_dry_tolerance,Bint,sea_level
    real(kind=8) :: B(0:2,0:2), xedge(0:2), yedge(0:2)
    logical :: use_Bint

    ! Common block variables
    real(kind=8) :: dxmin, dymin

    common /comfine/ dxmin,dymin

    character(len=364) :: format_string
    character(len=*), parameter :: clawfile = 'claw.data'
    character(len=*), parameter :: amrfile = 'amr.data'
    character(len=*), parameter :: outfile = 'fort.amr'
    character(len=*), parameter :: dbugfile = 'fort.debug'
    character(len=*), parameter :: matfile = 'fort.nplot'
    character(len=*), parameter :: parmfile = 'fort.parameters'

    ! Open parameter and debug files
    open(dbugunit,file=dbugfile,status='unknown',form='formatted')
    open(parmunit,file=parmfile,status='unknown',form='formatted')

    maxthreads = 1    !! default, if no openmp

    ! Open AMRClaw primary parameter file
    call opendatafile(inunit,clawfile)

    ! Number of space dimensions, not really a parameter but we read it in and
    ! check to make sure everyone is on the same page. 
    read(inunit,"(i1)") ndim  
    if (ndim /= 2) then
        print *,'Error ***   ndim = 2 is required,  ndim = ',ndim
        print *,'*** Are you sure input has been converted'
        print *,'*** to Clawpack 5.x form?'
        stop
    endif
          
    ! Domain variables
    read(inunit,*) xlower, ylower
    read(inunit,*) xupper, yupper
    ! that's all we need from claw.data, plus these things:
    coordinate_system = 2
    nvar = 3
    naux = 3
    rest = .false.
    earth_radius = 6367500.0d0
    mod_dry_tolerance = 1d-4
    sea_level = 0.d0

    call read_topo_settings(rest)     ! specifies topography (bathymetry) files
    call set_fgout(rest,nvar)         ! Fixed grid settings
    call set_gauges(rest, nvar, naux) ! Set gauge output
    !write(6,*) '+++ x,y: ',gauges(1)%x,gauges(1)%y
    call set_fgmax()

    ! new code to compute B around each gauge:

    ! eventually need to read desired resolution dx,dy
    dx = 1.d0 / (3.d0*3600.d0)
    dy = dx

    do ii=1,num_gauges
        x = gauges(ii)%x
        y = gauges(ii)%y
        !write(6,*) '+++ x,y: ',gauges(ii)%x,gauges(ii)%y
        
        i1 = int(floor((x-xlower)/dx)) - 1
        j1 = int(floor((y-ylower)/dy)) - 1
        !write(6,*) '+++x,y,lower, i1,j1: ',x,y,xlower,ylower,i1,j1

        do i=0,2
            xm = xlower + (i1+i)*dx
            xp = xm + dx
            xedge(i) = xm
            do j=0,2
                ym = ylower + (j1+j)*dy
                yp = ym + dy
                yedge(j) = ym
                !write(6,661) i,j,xm,xp,ym,yp,-999.d0
                call cellgridintegrate(topo_integral,xm,xp,ym,yp)
                B(i,j) = topo_integral / (dx * dy)
                if (coordinate_system == 2) then
                    aa = deg2rad * earth_radius**2                      &
                            * (sin(yp * deg2rad) - sin(ym * deg2rad)) / dy
                    B(i,j) = B(i,j) / aa
                endif
 661            format('i=',i1,' j=',i1,' xm=',f12.7,' xp=',f12.7,' ym='f12.7, &
                       ' yp=',f12.7,/,'        B=',f10.3)
                !write(6,661) i,j,xm,xp,ym,yp,B(i,j)
            enddo

        enddo

        write(6,*) '=============================================================================='
 662    format('Gauge ',i6,'  At x=',f12.7, '        y=',f12.7)
        write(6,662) gauges(ii)%gauge_num, x, y
 666    format('          xcenter =',f12.7, ' ycenter =',f12.7)
        xcenter = 0.5d0*(xedge(1)+xedge(2))
        ycenter = 0.5d0*(yedge(1)+yedge(2))
        write(6,666) xcenter, ycenter

        ! compute interpolated B value:
        iindex = int(0.5d0 + (x - xedge(0))/dx) - 1  ! shifted since first cell ind=0
        jindex = int(0.5d0 + (y - yedge(0))/dy) - 1
        xcent = xedge(0) + (iindex + 0.5d0)*dx
        ycent = yedge(0) + (jindex + 0.5d0)*dy
        xoff = (x - xcent) / dx
        yoff = (y - ycent) / dy
692     format('+++ iindex = ',i1, ' jindex = ',i1, '  xoff = ',f6.3, '  yoff = ', f6.3)
        !write(6,692) iindex,jindex,xoff,yoff
        Bint = (1.d0 - xoff) * (1.d0 - yoff) * B(iindex,  jindex)         &
                        + xoff*(1.d0 - yoff) * B(iindex+1,jindex)         &
                      + (1.d0 - xoff) * yoff * B(iindex,  jindex+1)       &
                               + xoff * yoff * B(iindex+1,jindex+1)

        use_Bint = (B(iindex  ,jindex  ) < sea_level - mod_dry_tolerance) .and. &
                   (B(iindex+1,jindex  ) < sea_level - mod_dry_tolerance) .and. &
                   (B(iindex  ,jindex+1) < sea_level - mod_dry_tolerance) .and. &
                   (B(iindex+1,jindex+1) < sea_level - mod_dry_tolerance)

 668    format('Computed cell average B in cell = ', f10.3, '  interpolated B = ', f10.3)
        write(6,668) B(1,1), Bint
        if (use_Bint) then
 690        format('With sea_level = ',f10.3, '   will use interpolated B = ', f10.3)
            write(6,690) sea_level, Bint
        else
 691        format('With sea_level = ',f10.3, '   will use cell average B = ', f10.3)
            write(6,691) sea_level, B(1,1)
        endif
    

        ! print grid:
        write(6,*) ''
        write(6,*) '               |                    |                   |               |'
        !write(6,*) '---------------+--------------------+-------------------+---------------+--'
 663    format(f12.7,'    |                    |                   |               |')
 667    format(f12.7,'  --+--------------------+-------------------+---------------+--')
 664    format('                |', f13.3, '       |', f13.3, '      |', f13.3, '  |')
!665    format(4f20.7)
 665    format(f21.7,f23.7,f18.7,f17.7)
        write(6,667) yedge(2) + dy
        do j=2,0,-1
            write(6,664) (B(i,j), i=0,2)
            write(6,667) yedge(j)
        enddo
        write(6,*) '               |                    |                   |               |'
        write(6,665) (xedge(i), i=0,2), xedge(2)+dx
        write(6,*) '=============================================================================='
    
   
    enddo
        


end program computeB
