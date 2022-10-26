module progress
  use scatter
  implicit none

  !This module gives time development of particles

  contains

  subroutine emcd  !粒子のドリフトと散乱をΔtだけ進める関数
    implicit none
      
    integer n,loop
      
    tdt = t + dt

    do n=1,inum
      if (ipt(n) == 9) cycle

      t1 = t
      t_s = pt(n,4)
      t_r = pt(n,8)
      t_dp = pt(n,9)
      t_ey = pt(n,10)
      t_bap = pt(n,11)


      do loop = 1, MAX_LOOP
        if (t_s>tdt .and. t_r>tdt .and. t_dp>tdt .and. t_ey>tdt .and. t_bap>tdt) exit

        if (t_r <= tdt) then
          ipt(n) = 9
          exit

        else if (t_s <= t_dp .and. t_s <= t_ey .and. t_s <= t_bap) then
          tau = t_s - t1
          if (tau < 0) then
            write(*,*) "Error: tau < 0 at the ts loop in emcd!",n,loop,t_s,t1
            close(10)
            close(11)
            stop
          endif
          call drift(pt(n,:), ipt(n))
          call scat(pt(n,:), ipt(n))

        else if (t_dp < t_s .and. t_dp <= t_ey .and. t_dp <= t_bap) then
          tau = t_dp - t1
          if (tau < 0) then
            write(*,*) "Error: tau < 0 at the dp loop in emcd!",n,loop,t_dp,t1
            close(10)
            close(11)
            stop
          endif
          call drift(pt(n,:), ipt(n))
          call spin_relaxation_dp(pt(n,:),ipt(n))
        
        else if (t_ey < t_s .and. t_ey < t_dp .and. t_ey <= t_bap) then
          tau = t_ey - t1
          if (tau < 0) then
            write(*,*) "Error: tau < 0 at the ey loop in emcd!",n,loop,t_ey,t1
            close(10)
            close(11)
            stop
          endif
          call drift(pt(n,:), ipt(n))
          call spin_relaxation_ey(pt(n,:),ipt(n))
          
        else if (t_bap < t_s .and. t_bap < t_dp .and. t_bap < t_ey) then
          tau = t_bap - t1
          if (tau < 0) then
            write(*,*) "Error: tau < 0 at the bap loop in emcd!",n,loop,t_bap,t1
            close(10)
            close(11)
            stop
          endif
          call drift(pt(n,:), ipt(n))
          call spin_relaxation_bap(pt(n,:))

        else 
          write(*,*) 'Error : something is wrong in the emcd loop!',n,loop,t_s,t_dp,t_ey,t_bap
          close(10)
          close(11)
          stop
        end if
        
        if (loop == MAX_LOOP) then
          write(*,*) "Error: progress loop didnt end!", n, t_s,t_r,t_dp, t_ey, t_bap
          close(10)
          close(11)
          stop
        endif
      enddo

      tau = tdt - t1
      call drift(pt(n,:), ipt(n))
      pt(n,4) = t_s
      pt(n,9) = t_dp
      pt(n,10) = t_ey
      pt(n,11) = t_bap
    enddo

    return
  end subroutine emcd


 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



  subroutine drift(p, iv)
    implicit none
    
    real(double), intent(inout) :: p(11)
    integer, intent(inout) :: iv
    real(double) sq,dkx,dky,hmt,sk,qht
    integer i,j
    
    if (iv==9) return
     
    i   = max(1,min(int(p(5)/dx+1.5d0),mx1))
    j   = max(1,min(int(p(6)/dy+1.5d0),my1))
    qht = q / h * tau
    dkx = -qht*fx(i,j)
    dky = -qht*fy(i,j)
    hmt = hm(iv)*tau
    sk  = p(1)*p(1)+p(2)*p(2)+p(3)*p(3)
    sq  = sqrt(1d0 + 4d0 * af(iv) * hhm(iv) * sk)
    p(5)  = p(5)+hmt*(p(1)+0.5d0*dkx)/sq
    p(6)  = p(6)+hmt*(p(2)+0.5d0*dky)/sq
    p(1)  = p(1)+dkx
    p(2)  = p(2)+dky
    
    call surf(p)
    return
  end subroutine drift



  subroutine spin_relaxation_bap(pt)
    implicit none
    real(double) ::pt(11),pt_before,hden

    pt_before = pt(7)

    pt(7) = -pt(7)

    bapef(jt) = abs(pt(7) - pt_before)       !y方向のスピンの変化量を蓄積

    hden = inum * epp / (xmax * ymax * wd)        !ホール濃度（=キャリア濃度）

    t1 = t_bap
    t_bap = t1 - h**2 / (param_bap * hden * lattice**6 * sqrt(am(iv))**3) * log(rnd())

    return
  end subroutine spin_relaxation_bap



  subroutine spin_relaxation_ey(pt,ipt)
    implicit none

    real(double) ::pt(11),pt_before
    integer :: ipt

    pt_before = pt(7)

    pt(7) = -pt(7)

    eyef(jt)  = abs(pt(7) - pt_before)          !y方向のスピンの変化量を蓄積

    t1 = t_ey
    t_ey = t1 - 9d0/8d0 / (del*q/(bgap*q + del*q))**2 &
          / (1d0-am(ipt)/am0)**2 / (gstate_en*q * kb * temp) * (bgap*q)**2 * t_p &
          * log(rnd())                             !次の散乱時刻

    return
  end subroutine spin_relaxation_ey


    
  subroutine spin_relaxation_dp(pt,ipt)
    implicit none

    real(double) :: pt(11), pt_before
    integer :: ipt

    pt_before = pt(7)
    alpha_dp = 4d0*eta_dp/sqrt(3d0-eta_dp) * am(iv)/am0
    gamma_dp = alpha_dp * h**3 / sqrt(2d0 * am(iv)**3 * bgap*q)

    pt(7) = -pt(7)

    dpef(jt)  = abs(pt(7) - pt_before)                             !y方向のスピンの変化量を蓄積

    t1 = t_dp
    t_dp = t1 - h**8 / (16d0 * kb*temp * am(ipt)**3 * (gamma_dp * gstate_en*q)**2 * t_p ) &
              * log(rnd())                                            !次にDP効果が起きる時刻

    return
  end subroutine spin_relaxation_dp



  subroutine surf(p)  !境界処理 : ワープ
    implicit none
  
    real(double), intent(inout) :: p(11)
  
    if (p(5) < 0) then
      p(5) = xmax + p(5)
    elseif (p(5) >= xmax) then
      p(5) = p(5) - xmax
    endif

    if (p(6) <0) then
      p(6) = ymax + p(6)
    elseif (p(6) >= ymax) then
      p(6) = p(6) - ymax
    endif

    return
  end subroutine surf



  subroutine renew
    implicit none

    integer :: n,loop,nmax
  
    nmax = inum
  
    do n= 1, nmax
      do loop = 1, MAX_LOOP
        if (ipt(inum) /= 9) exit
        inum=inum-1
        if (loop == MAX_LOOP) then
          write(*,*) "Error: renew loop 1 didnt end!"
          close(10)
          close(11)
          stop
        endif
      enddo

      if (n >= inum) exit
      
      if (ipt(n)==9) then
        pt(n,:)=pt(inum,:)
        ipt(n)=ipt(inum)
        ipt(inum)=9
        inum=inum-1
      endif
  
      do loop = 1, MAX_LOOP
        if (ipt(inum) /= 9) exit
        inum=inum-1
        if (loop==MAX_LOOP) then
          write(*,*) "Error: renew loop 2 didnt end!"
          close(10)
          close(11)
          stop
        endif
      enddo

      if (n >= inum) exit
    enddo

    return
  end subroutine renew

end module progress