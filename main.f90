! SMAC法
! 角柱流れ

program fortran

implicit none

! 変数
integer,parameter::nxc = 200
integer,parameter::nyc = 100
integer,parameter::nxd = nxc + 1
integer,parameter::nyd = nyc + 1
real(8),parameter::lx = 0.2d0
real(8),parameter::ly = 0.1d0
real(8),parameter::dx = 0.001d0
real(8),parameter::dy = 0.001d0
real(8),parameter::uin = 0.01d0
real(8),parameter::lt = 50d0

real(8),parameter::Re = 1000d0 !=uin*lx/visc
real(8),parameter::Co = 0.1d0 !クーラン数

real(8), parameter::dens  = 1.0e+3 !密度[kg/m^3]
real(8), parameter::c     = 4000d0 !比熱[j/(kg*K)]
real(8), parameter::ramda_cavity = 0.6d0 !熱伝導率[w/(m*K)]
real(8), parameter::ramda_mold   = 300d0 !
real(8), parameter::heat_transfer_coefficient = 837d0 !熱伝達係数
real(8), parameter::volume = dx*dy
real(8), parameter::surface = dx

real(8), parameter::Liquid_initial_temp = 0d0
real(8), parameter::Prizm_initial_temp = 100d0

real(8)::dt
real(8)::visc,arufa,err_total
real(8)::start_time,end_time
real(8)::time_Step_limit1,time_Step_limit2
integer::i,j,time,nt,outputstep
integer::ic,jc
real(8)::x,y
real(8),dimension(0:nxd+1,0:nyc+1)::u
real(8),dimension(0:nxc+1,0:nyd+1)::v
real(8),dimension(0:nxd+1,0:nyc+1)::u_aux
real(8),dimension(0:nxc+1,0:nyd+1)::v_aux
real(8),dimension(0:nxd+1,0:nyc+1)::flux_u
real(8),dimension(0:nxc+1,0:nyd+1)::flux_v
real(8),dimension(0:nxc+1,0:nyc+1)::p
real(8),dimension(0:nxc+1,0:nyc+1)::phi
real(8),dimension(0:nxc+1,0:nyc+1)::dp
real(8),dimension(nxc    ,nyc    )::theta

real(8),dimension(0:nxc+1,0:nyc+1)::temp
real(8),dimension(0:nxc+1,0:nyc+1)::temp_old
real(8),dimension(0:nxc+1,0:nyc+1)::temp_conv

integer,parameter::prizm_left   = 90
integer,parameter::prizm_right  = 110
integer,parameter::prizm_bottom = 40
integer,parameter::prizm_top    = 60

! 初期設定
dt = Co*dx/uin
nt = nint(lt/dt)
time_Step_limit1 = dens*c*volume/(ramda_cavity*surface/dx*4d0)
time_Step_limit2 = dens*c*volume/((ramda_cavity*surface/dx*3d0) &
                + (surface/(1d0/heat_transfer_coefficient + dx/(2d0*ramda_cavity) + dx/(2d0*ramda_mold))))
print *, time_Step_limit1,time_Step_limit2

! visc = 1d0/1e+6
visc = uin*lx/Re !0.01*0.2/1000
err_total = 1d0/1e+12
arufa = 1.99d0

u(:,:) = uin
v(:,:) = 0d0
u_aux(:,:) = 0d0
v_aux(:,:) = 0d0
p(:,:) = 0d0
phi(:,:) = 0d0
theta(:,:) = 0d0
temp(:,:) = Liquid_initial_temp
flux_u(:,:) = 0d0
flux_v(:,:) = 0d0

u(0,:) = uin*2d0 - u(1,:)
u(:,0) = uin*2d0 - u(:,1)
u(:,nyc+1) = uin*2d0 - u(:,nyc)

v(0,:) = 0d0
v(:,0) = 0d0
v(:,nyd) = 0d0    

temp(prizm_left:prizm_right,prizm_bottom:prizm_top) = Prizm_initial_temp

temp(10:20,10:20) = Prizm_initial_temp

call setBoundary(u,v,uin,prizm_left,prizm_right,prizm_top,prizm_bottom)

! メインルーチン
! outputstep = 500
outputstep = 50

call cpu_time(start_time)
write(*,*) dt,nt
! do time = 1,nt
do time = 1,5000
    write(*,*) time
    call computeAuxiallyVelocity(u_aux,v_aux,u,v,p,dx,dy,dt,visc,uin, &
                        prizm_left,prizm_right,prizm_top,prizm_bottom)
    call computeDivergenceAuxiallyVelocity(theta,u_aux,v_aux,dx,dy,dt, &
                        prizm_left,prizm_right,prizm_top,prizm_bottom)
    call computePressurePoisson(p,phi,dp,u,v,dx,dy,theta,err_total,arufa,dens,visc,dt, &
                        prizm_left,prizm_right,prizm_top,prizm_bottom)
    call computeVelocity(u,v,u_aux,v_aux,dt,dx,dy,dens,phi,uin, &
                        prizm_left,prizm_right,prizm_top,prizm_bottom)

    call computeTemp(temp,temp_old,temp_conv,flux_u,flux_v,u,v)

    if(mod(time,outputstep)==0)then
        call output(u,v,p,phi,time,temp)
    elseif(time > 4900 .and. time < 5000)then
        call output(u,v,p,phi,time,temp)
    endif
enddo
call cpu_time(end_time)
write(*,*) "elapsed time = ",end_time - start_time

open(unit=30,file="elapsed.txt")
write(30,*) end_time - start_time
close(30)

contains 

subroutine output(u,v,p,phi,time,temp)
    implicit none
    real(8),intent(in)::u(0:,0:)
    real(8),intent(in)::v(0:,0:)
    real(8),intent(in)::p(0:,0:)
    real(8),intent(in)::phi(0:,0:)
    integer,intent(in)::time
    real(8),intent(in)::temp(0:,0:)

    integer::i,j
    real(8)::uout,vout
    character*128 filename1,filename2,filename3

    write(filename1,'("pres/pres",i5.5,".txt")') time
    write(filename2,'("phi/phi",i5.5,".txt")') time
    write(filename3,'("temp/temp",i5.5,".txt")') time
    open(unit=10,file=filename1)
    open(unit=20,file=filename2)
    open(unit=30,file=filename3)
    do j=1,nyc
    do i=1,nxc
        uout = (u(i,j)+u(i+1,j))/2d0
        vout = (v(i,j)+v(i,j+1))/2d0
        if(prizm_left <= i .and. i<=prizm_right .and. prizm_bottom<=j .and. j<=prizm_top)then
            uout = 0d0
            vout = 0d0
        elseif(i==prizm_right.and.prizm_bottom<=j .and. j<=prizm_top)then
            uout = 0d0
        elseif(prizm_left <= i .and. i<=prizm_right.and.j==prizm_top)then
            vout = 0d0
        endif
        write(10,'(2i8,3f23.14)') i,j,p(i,j),uout,vout
        write(20,'(2i8,3f23.14)') i,j,phi(i,j),uout,vout
        write(30,'(2i8,3f23.14)') i,j,temp(i,j),uout,vout
    enddo
    enddo
    close(10)
    close(20)
    close(30)

end subroutine

subroutine setBoundary(u,v,uin,prizm_left,prizm_right,prizm_top,prizm_bottom)
    implicit none
    real(8),intent(inout)::u(0:,0:)
    real(8),intent(inout)::v(0:,0:)
    real(8),intent(in) ::uin
    integer,intent(in) ::prizm_left,prizm_right,prizm_top,prizm_bottom

    u(0,:) = uin*2d0 - u(1,:)
    u(:,0) = uin*2d0 - u(:,1)
    u(:,nyc+1) = uin*2d0 - u(:,nyc)
    
    v(0,:) = 0d0
    v(:,1) = 0d0
    v(:,nyd) = 0d0    
    
    ! 流出境界条件
    u(nxd+1,:) = u(nxd,:)
    v(nxc+1,:) = v(nxc,:)
    
    ! 角柱表面
    u(prizm_left +1:prizm_right,prizm_bottom            ) = -u(prizm_left+1 :prizm_right,prizm_bottom-1          )
    u(prizm_left +1:prizm_right,prizm_top               ) = -u(prizm_left+1 :prizm_right,prizm_top+1             )
    u(prizm_left +1            ,prizm_bottom  :prizm_top) = -u(prizm_left-1             ,prizm_bottom  :prizm_top)
    u(prizm_right              ,prizm_bottom  :prizm_top) = -u(prizm_right+2            ,prizm_bottom  :prizm_top)

    v(prizm_left              ,prizm_bottom+1:prizm_top) = -v(prizm_left-1            ,prizm_bottom+1:prizm_top  )
    v(prizm_right             ,prizm_bottom+1:prizm_top) = -v(prizm_right+1           ,prizm_bottom+1:prizm_top  )
    v(prizm_left:prizm_right  ,prizm_bottom+1          ) = -v(prizm_left:prizm_right  ,prizm_bottom-1            )
    v(prizm_left:prizm_right  ,prizm_top               ) = -v(prizm_left:prizm_right  ,               prizm_top+2)
    
    ! 角柱と流体の境界
    u(prizm_left            ,prizm_bottom:prizm_top) = 0d0
    u(prizm_right+1         ,prizm_bottom:prizm_top) = 0d0
    v(prizm_left:prizm_right,prizm_bottom          ) = 0d0
    v(prizm_left:prizm_right,prizm_top+1           ) = 0d0
    
    ! 角柱内部
    u(prizm_left+2:prizm_right-1,prizm_bottom+1:prizm_top-1) = 0d0
    v(prizm_left+1:prizm_right-1,prizm_bottom+2:prizm_top-1) = 0d0

end subroutine

subroutine computeAuxiallyVelocity(u_aux,v_aux,u,v,p,dx,dy,dt,visc,uin, &
                        prizm_left,prizm_right,prizm_top,prizm_bottom)
    implicit none
    real(8),intent(inout)::u_aux(0:,0:)
    real(8),intent(inout)::v_aux(0:,0:)
    real(8),intent(in   )::u(0:,0:)
    real(8),intent(in   )::v(0:,0:)
    real(8),intent(in   )::p(0:,0:)
    real(8),intent(in) ::dx,dy,dt,visc,uin
    integer,intent(in) ::prizm_left,prizm_right,prizm_top,prizm_bottom

    do jc = 1,nyc
    do i  = 1,nxd
        j = jc
        ic = i
        if(prizm_left <= i .and. i<=prizm_right .and. prizm_bottom<=jc .and. jc<=prizm_top)cycle
        if(i==prizm_right+1.and.prizm_bottom<=jc.and.jc<=prizm_top)cycle !角柱の右辺上
        u_aux(i,jc) = u(i,jc) -dt*( (u(i  ,jc ) + u(i-1 ,jc ))/2d0*(u(i  ,jc  ) - u(i-1,jc  ))/dx + &
                                    (u(i+1,jc ) + u(i   ,jc ))/2d0*(u(i+1,jc  ) - u(i  ,jc  ))/dx + &
                                    (v(ic ,j+1) + v(ic-1,j+1))/2d0*(u(i  ,jc+1) - u(i  ,jc  ))/dy + &
                                    (v(ic ,j  ) + v(ic-1,j  ))/2d0*(u(i  ,jc  ) - u(i  ,jc-1))/dy)/2d0 + &
                              dt*visc*((u(i+1,jc  ) - 2d0*u(i,jc) + u(i-1,jc  ))/dx**2 + &
                                       (u(i  ,jc+1) - 2d0*u(i,jc) + u(i  ,jc-1))/dy**2) - &
                                dt*(p(ic,jc) - p(ic-1,jc))/dens/dx 
    enddo
    enddo

    do j  = 1,nyd
    do ic = 1,nxc
        jc = j
        i  = ic
        if(prizm_left <= ic .and. ic<=prizm_right .and. prizm_bottom<=j .and. j<=prizm_top)cycle
        if(prizm_left<=ic .and. ic<=prizm_right .and. j==prizm_top+1)cycle !角柱の上辺上
        v_aux(ic,j) = v(ic,j) -dt*( (u(i  ,jc ) + u(i  ,jc-1))/2d0*(v(ic  ,j  ) - v(ic-1,j  ))/dx + &
                                    (u(i+1,jc ) + u(i+1,jc-1))/2d0*(v(ic+1,j  ) - v(ic  ,j  ))/dx + &
                                    (v(ic ,j+1) + v(ic ,j   ))/2d0*(v(ic  ,j+1) - v(ic  ,j  ))/dy + &
                                    (v(ic ,j  ) + v(ic ,j-1 ))/2d0*(v(ic  ,j  ) - v(ic  ,j-1))/dy)/2d0 + &
                              dt*visc*((v(ic+1,j  ) - 2d0*v(ic,j) + v(ic-1,j  ))/dx**2 + &
                                       (v(ic  ,j+1) - 2d0*v(ic,j) + v(ic  ,j-1))/dy**2) - &
                                dt*(p(ic,jc) - p(ic,jc-1))/dens/dy
    enddo
    enddo

    call setBoundary(u_aux,v_aux,uin,prizm_left,prizm_right,prizm_top,prizm_bottom)

end subroutine

subroutine computeDivergenceAuxiallyVelocity(theta,u_aux,v_aux,dx,dy,dt,prizm_left,prizm_right,prizm_top,prizm_bottom)
    implicit none
    real(8),intent(inout)::theta(:,:)
    real(8),intent(in   )::u_aux(0:,0:)
    real(8),intent(in   )::v_aux(0:,0:)

    real(8),intent(in)::dt,dx,dy
    integer,intent(in)::prizm_left,prizm_right,prizm_top,prizm_bottom

    do jc = 1,nyc
    do ic = 1,nxc
        j = jc
        i = ic
        ! if(prizm_left<=ic.and.ic<=prizm_right .and. prizm_bottom<=jc.and.jc<=prizm_top)cycle
        theta(ic,jc) = (u_aux(i+1,jc) - u_aux(i,jc))/dx + (v_aux(ic,j+1) - v_aux(ic,j))/dy
    enddo
    enddo

end subroutine


subroutine computePressurePoisson(p,phi,dp,u,v,dx,dy,theta,err_total,arufa,dens,visc,dt, &
                                  prizm_left,prizm_right,prizm_top,prizm_bottom)
    implicit none
    real(8),intent(inout)::p(0:,0:)
    real(8),intent(inout)::phi(0:,0:)
    real(8),intent(inout)::dp(0:,0:)
    real(8),intent(in)::u(0:,0:)
    real(8),intent(in)::v(0:,0:)
    real(8),intent(in)::theta(:,:)

    real(8),intent(in)::dx,dy,err_total,arufa,dens,visc,dt
    integer,intent(in)::prizm_left,prizm_right,prizm_top,prizm_bottom
    real(8)::err,err_phi,err_dp,err_p
    integer count

    dp(:,:) = 0d0
    phi(0:,0:) = 0d0
    err = 1d0
    count = 0

    do while(err > err_total)
        ! count = count + 1
        err_p = 0d0
        err_dp = 0d0
        err_phi = 0d0
        do jc=1,nyc
        do ic=1,nxc
            j = jc
            i = ic
            if(prizm_left <= ic .and. ic<=prizm_right .and. prizm_bottom<=jc .and. jc<=prizm_top)cycle
            dp(ic,jc) =( (dy**2d0*(phi(ic+1,jc  ) + phi(ic-1,jc  )) + &
                          dx**2d0*(phi(ic  ,jc+1) + phi(ic  ,jc-1)) - &
                          dx**2d0*dy**2d0*theta(ic,jc)*dens/dt)/(2d0*(dx**2d0+dy**2d0))) - phi(ic,jc)
            phi(ic,jc) = phi(ic,jc) + arufa*dp(ic,jc)
            err_dp   = err_dp  + abs( dp(ic,jc))
            err_phi  = err_phi + abs(phi(ic,jc))
        enddo
        enddo
        phi(prizm_left            ,prizm_bottom:prizm_top) = phi(prizm_left-1          ,prizm_bottom:prizm_top)
        phi(prizm_right           ,prizm_bottom:prizm_top) = phi(prizm_right+1         ,prizm_bottom:prizm_top)
        phi(prizm_left:prizm_right,prizm_bottom          ) = phi(prizm_left:prizm_right,prizm_bottom-1        )
        phi(prizm_left:prizm_right,prizm_top             ) = phi(prizm_left:prizm_right,prizm_top   +1        )

        ! 境界条件　圧力勾配０
        phi(nxc+1,:) = phi(nxc,:)

        if(err_phi < 1.0d-20)then
            err_phi = 1d0
        endif
        err = err_dp/err_phi

        ! write(*,*) err
    enddo
    ! p(0:,0:) = p(0:,0:) + phi(0:,0:)
    p(:, :) = p(:, :) + phi(:, :) !修正ポイント2
    
end subroutine

subroutine computeVelocity(u,v,u_aux,v_aux,dt,dx,dy,dens,phi,uin,prizm_left,prizm_right,prizm_top,prizm_bottom)
    implicit none
    real(8),intent(inout)::u(0:,0:)
    real(8),intent(inout)::v(0:,0:)
    real(8),intent(inout)::phi(0:,0:)
    real(8),intent(in   )::u_aux(0:,0:)
    real(8),intent(in   )::v_aux(0:,0:)

    integer,intent(in)::prizm_left,prizm_right,prizm_top,prizm_bottom
    real(8),intent(in)::dt,dx,dy,dens,uin

    do jc=1,nyc
    do i =1,nxd
        j  = jc
        ic = i
        if(prizm_left <= i .and. i<=prizm_right .and. prizm_bottom<=jc .and. jc<=prizm_top)cycle
        if(prizm_bottom <= jc .and. jc <= prizm_top .and. i == prizm_right+1)cycle
        u(i,jc) = u_aux(i,jc) - dt*(phi(ic,jc)-phi(ic-1,jc))/dx/dens
    enddo
    enddo

    do j=1,nyd
    do ic=1,nxc
        jc = j
        i  = ic
        if(prizm_left <= ic .and. ic<=prizm_right .and. prizm_bottom<=j .and. j<=prizm_top)cycle
        if(prizm_left <= ic .and. ic <= prizm_right .and. j == prizm_top+1)cycle
        v(ic,j) = v_aux(ic,j) - dt*(phi(ic,jc)-phi(ic,jc-1))/dy/dens
    enddo
    enddo
    
    call setBoundary(u,v,uin,prizm_left,prizm_right,prizm_top,prizm_bottom)

end subroutine

subroutine computeTemp(temp,temp_old,temp_conv,flux_u,flux_v,u,v)
    real(8),intent(inout)::temp(0:,0:)
    real(8),intent(inout)::temp_old(0:,0:)
    real(8),intent(inout)::temp_conv(0:,0:)
    real(8),intent(inout)::flux_u(0:,0:)
    real(8),intent(inout)::flux_v(0:,0:)
    real(8),intent(in   )::u(0:,0:)
    real(8),intent(in   )::v(0:,0:)

    integer::i,j,id,jd
    real(8)::conv_u1,conv_u2,conv_v1,conv_v2
    real(8)::temp_w
    real(8)::u_diff,v_diff
    real(8)::res_e,res_w,res_n,res_s
    real(8)::max_temp,min_temp
    real(8)::test1,test2

    real(8)::heat_conv

    temp_old(:,:) = temp(:,:)
    temp_conv(:,:) = temp(:,:)

    !最大値と最小値の設定
    ! max_temp = maxval(temp)
    ! min_temp = minval(temp)
    max_temp = Liquid_initial_temp
    min_temp = Prizm_initial_temp
    do j=1,nyc
    do i=1,nxc
        if(prizm_left <= i .and. i<=prizm_right .and. prizm_bottom<=j .and. j<=prizm_top)cycle
        if(temp(i,j) > max_temp) max_temp = temp(i,j)
        if(temp(i,j) < min_temp) min_temp = temp(i,j)
    enddo
    enddo

    print *, max_temp,min_temp

    ! 移流の影響を計算
    do j =1,nyc
    do id=1,nxd
        i = id
        heat_conv = (temp_old(i-1,j) - temp_old(i,j))*surface*dt*u(id,j)
        if(u(i,j) > 0d0)then
            temp(i  ,j) = temp(i  ,j) + heat_conv
            temp(i-1,j) = temp(i-1,j) - heat_conv            
        else
            temp(i  ,j) = temp(i  ,j) - heat_conv
            temp(i-1,j) = temp(i-1,j) + heat_conv            
        endif
        ! if(i==111.and.j==40)then
        !     print *, temp(i,j),temp_old(i,j),heat_conv
        ! endif

        ! Flux 6.25式（数値流体解析の基礎）
        flux_u(i,j) = u(i,j)*(temp(i,j)+temp(i-1,j))/2d0 &
                    - (visc/dx + abs(u(i,j))/2d0)*(temp(i,j)-temp(i-1,j))
    enddo
    enddo
    do jd=1,nyd
    do i =1,nxc
        j = jd
        heat_conv = (temp_old(i,j-1) - temp_old(i,j))*surface*dt*v(i,jd)
        if(v(i,j) > 0d0)then
            temp(i,j  ) = temp(i  ,j  ) + heat_conv
            temp(i,j-1) = temp(i  ,j-1) - heat_conv            
        else
            temp(i,j  ) = temp(i  ,j  ) - heat_conv
            temp(i,j-1) = temp(i  ,j-1) + heat_conv            
        endif
        ! if(i==111.and.j==39)then
        !     print *, temp(i,j),heat_conv
        ! endif

        ! Flux 6.25式（数値流体解析の基礎）
        flux_v(i,j) = v(i,j)*(temp(i,j)+temp(i,j-1))/2d0 &
                    - (visc/dy + abs(v(i,j))/2d0)*(temp(i,j)-temp(i,j-1))
    enddo
    enddo
    

    do j=1,nyc
    do i=1,nxc
        if(prizm_left <= i .and. i<=prizm_right .and. prizm_bottom<=j .and. j<=prizm_top)cycle
        if(i==prizm_left-1 .and. prizm_bottom<=j .and. j<=prizm_top)then !右側が角柱
            ! eastが東、角柱があるのはeast
            res_e = 1d0/heat_transfer_coefficient + dx/(2d0*ramda_cavity) + dx/(2d0*ramda_mold)
            res_w = dx/ramda_cavity
            res_n = dy/ramda_cavity
            res_s = dy/ramda_cavity
        elseif(i==prizm_right+1 .and. prizm_bottom<=j .and. j<=prizm_top)then !左側が角柱
            ! eastが東、角柱があるのはwest
            res_e = dx/ramda_cavity
            res_w = 1d0/heat_transfer_coefficient + dx/(2d0*ramda_cavity) + dx/(2d0*ramda_mold)
            res_n = dy/ramda_cavity
            res_s = dy/ramda_cavity
        elseif(j==prizm_bottom-1 .and. prizm_left <= i .and. i<=prizm_right)then !角柱があるのは上側
            ! 北側がnorth、角柱があるのは北側
            res_e = dx/ramda_cavity
            res_w = dx/ramda_cavity
            res_n = 1d0/heat_transfer_coefficient + dy/(2d0*ramda_cavity) + dy/(2d0*ramda_mold)
            res_s = dy/ramda_cavity
        elseif(j==prizm_top+1 .and. prizm_left <= i .and. i<=prizm_right)then
            res_e = dx/ramda_cavity
            res_w = dx/ramda_cavity
            res_n = dy/ramda_cavity
            res_s = 1d0/heat_transfer_coefficient + dy/(2d0*ramda_cavity) + dy/(2d0*ramda_mold)
        else
            res_e = dx/ramda_cavity
            res_w = dx/ramda_cavity
            res_n = dy/ramda_cavity
            res_s = dy/ramda_cavity
        endif
        ! TopCASTのソース(熱伝導計算)
        ! temp w = temp(i,j) + temp_conv(i,j)
        ! if(i==111.and.j==39)then
        !     print *, "temp",temp(i,j),temp_old(i,j)
        ! endif
        temp(i,j) = temp(i,j) + dt/(dens*c*volume)*surface &
                                  * ( (temp_old(i+1,j  ) - temp_old(i,j))/res_e & !右（東）
                                    + (temp_old(i-1,j  ) - temp_old(i,j))/res_w & !左（西）
                                    + (temp_old(i  ,j+1) - temp_old(i,j))/res_n &
                                    + (temp_old(i  ,j-1) - temp_old(i,j))/res_s) !修正ポイント1
        test1 = dt/(dens*c*volume)*surface &
        * ( (temp_old(i+1,j  ) - temp_old(i,j))/res_e & !右（東）
          + (temp_old(i-1,j  ) - temp_old(i,j))/res_w & !左（西）
          + (temp_old(i  ,j+1) - temp_old(i,j))/res_n &
          + (temp_old(i  ,j-1) - temp_old(i,j))/res_s) !修正ポイント1
        ! フラックス分を追加
        temp(i,j) = temp(i,j) + dt/(dens*c*volume) &
                    * ( (-flux_u(i+1,j)+flux_u(i,j))/dx + (-flux_v(i,j+1)+flux_v(i,j))/dy )
        test2 = dt/(dens*c*volume)*surface &
        * ( (-flux_u(i+1,j)+flux_u(i,j))/dx + (-flux_v(i,j+1)+flux_v(i,j))/dy )
        ! TopCASTのソース参考（考え方は不明・・・）
        ! if(temp(i,j) > Prizm_initial_temp)then
        !     temp(i,j) = Prizm_initial_temp
        ! elseif(temp(i,j) < Liquid_initial_temp)then
        !     temp(i,j) = Liquid_initial_temp
        ! endif
        if(i==111.and.j==39)then
            print *, temp(i,j),temp_old(i,j)
            print *, (-flux_u(i+1,j)+flux_u(i,j))/dx,(-flux_v(i,j+1)+flux_v(i,j))/dy
            print *, test1,test2
            print *, flux_u(i,j),u(i,j),1d0/(dens*c*volume)
        endif
        ! ! TopCASTの式
        ! temp(i,j) = temp_old(i,j) - dt/volume*surface*(conv_u1 + conv_u2 + conv_v1 + conv_v2) &
        !                           + dt/(dens*c*volume)*surface &
        !                           * ( (temp_old(i+1,j  ) - temp_old(i  ,j  ))/res_e &
        !                             + (temp_old(i  ,j  ) - temp_old(i-1,j  ))/res_w &
        !                             + (temp_old(i  ,j+1) - temp_old(i  ,j  ))/res_n &
        !                             + (temp_old(i  ,j  ) - temp_old(i  ,j-1))/res_s)
        ! ! P160 3.149式
        ! temp(i,j) = temp_old(i,j) + dt*kappa/(c*dens)*( &
        !     (temp_old(i+1,j  ) - 2d0*temp_old(i,j) + temp_old(i-1,j  ))/(dx**2) &
        !   + (temp_old(i  ,j+1) - 2d0*temp_old(i,j) + temp_old(i  ,j-1))/(dy**2) ) &
        !   ! 移流
        !   - dt*(conv_u1 + conv_u2 + conv_v1 + conv_v2) &
        !   - dt*temp_old(i,j)*(u_diff/dx + v_diff/dy)

        ! test = dt/(dens*c*volume)*surface &
        ! * ( (temp_old(i+1,j  ) - temp_old(i  ,j  ))/res_e &
        !   + (temp_old(i  ,j  ) - temp_old(i-1,j  ))/res_w &
        !   + (temp_old(i  ,j+1) - temp_old(i  ,j  ))/res_n &
        !   + (temp_old(i  ,j  ) - temp_old(i  ,j-1))/res_s)
        ! if(i==111.and.j==39)then
        !     print *, temp(i,j)
        ! endif
        ! if(i==109.and.j==39)then
        !     print *, temp(i-2,j),temp(i,j),temp(i+2,j),temp(i+4,j)
        ! endif
    enddo
    enddo
    print *, temp(109-2,39),temp(109,39),temp(109+1,39),temp(109+2,39)
    print *, temp(91+2,39) ,temp(91,39) ,temp(91-1,39) ,temp(91-2,39)

    temp(0    ,:    ) = temp(1  ,:  )
    temp(nxc+1,:    ) = temp(nxc,:  )
    temp(:    ,0    ) = temp(:  ,1  )
    temp(:    ,nyc+1) = temp(:  ,nyc)

        

end subroutine

end program