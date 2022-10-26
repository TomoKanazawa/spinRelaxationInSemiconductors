module field
  use config_mod
  implicit none

  !This module gives carrier density and solves poisson's equation
  
  contains

  subroutine charge !This subroutine gives carrier distribution
    implicit none
  
  
    real(double) x,x1,xb,y,y1,yb
    integer i,j,n
  
    cden=0d0
  
    !---( セル電荷雲法 )----p.119参照
  
    do n=1,inum
      if (ipt(n)==9) cycle
      x=pt(n,5)/dx                                      !粒子のx座標
      y=pt(n,6)/dy                                      !粒子のy座標
      i=max(1,min(int(x)+1,mx))                         !粒子に最も近い格子番号(x成分)
      j=max(1,min(int(y)+1,my))                         !粒子に最も近い格子番号(y成分)
      xb=real(i-1)                                      !上記の格子のx座標
      yb=real(j-1)                                      !上記の格子のy座標
      x1=1d0-(x-xb)                                     !粒子から遠いほうの格子までの距離(x成分)
      y1=1d0-(y-yb)                                     !粒子から遠いほうの格子までの距離(y成分)
      cden(i,j)=cden(i,j)+x1*y1                         !(i,j)の格子に属する電荷雲
      cden(i+1,j)=cden(i+1,j)+(1d0-x1)*y1               !(i+1,j)の格子に属する電荷雲
      cden(i,j+1)=cden(i,j+1)+x1*(1d0-y1)               !(i,j+1)の格子に属する電荷雲
      cden(i+1,j+1)=cden(i+1,j+1)+(1d0-x1)*(1d0-y1)     !(i+1,j+1)の格子に属する電荷雲
    enddo
  
    cden=cden*epp/dx/dy                                 !格子内の電荷密度
    cden(1,:)   = cden(1,:)*2d0                         !端は２倍
    cden(mx1,:) = cden(mx1,:)*2d0                       !端は２倍
    cden(:,1)   = cden(:,1)*2d0                         !端は２倍
    cden(:,my1)   = cden(:,my1)*2d0                     !端は２倍
  
    return
  end subroutine charge

      
      
  subroutine poisson  !This subroutine solves poisson equation
    implicit none

    integer :: i,j,loop

    real(double) :: phi(mx1,my1), const(mx1,my1),error,accel,allow,before

    const = -cden * q / eps
    phi = 0d0     !境界条件(境界で phi=0 )の役割も兼ねる
    error = 0d0   !反復にようるphiの値の変化を格納(解からのズレ)
    allow = 1d-4  !ズレの許容値
    accel = 1.9d0 !SOR法の収束加速パラメータ

    do loop = 1, MAX_LOOP
      do i = 2, mx1-1
        do j = 2, my1-1
          before = phi(i,j)
          phi(i,j) = (1d0 - accel) * phi(i,j) &
                          + accel * (phi(i+1,j) + phi(i-1,j) + phi(i,j+1) + phi(i,j-1) &
                          - const(i,j) * dx**2) / 4d0
          if (error < abs(phi(i,j) - before)) then
            error = abs(phi(i,j) - before)
          end if
        end do
      end do

      if (error <= allow) then
        do i = 1, mx1 !境界条件
          fy(i,1)   = 0d0
          fy(i,my1) = 0d0
        end do

        do i = 1, my1 !境界条件
          fx(1,i)   = 0d0
          fx(mx1,i) = 0d0
        end do

        do i = 2, mx1-1
          do j = 2, my1-1 !f(x) = - dU/dx
            fx(i,j) = (phi(i-1,j) - phi(i+1,j)) / (2d0 * dx)
            fy(i,j) = (phi(i,j-1) - phi(i,j+1)) / (2d0 * dy)
          end do
        end do

        return
      end if
    end do

    write(*,*) "Error : the loop in the poisson subroutine didnt end!", error

    stop
  end subroutine poisson



    
    



  
end module field