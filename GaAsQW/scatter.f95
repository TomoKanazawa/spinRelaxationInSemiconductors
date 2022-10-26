module scatter
  use initialize_cnp_mod
  implicit none

  !This module gives scatting rate as a arrangement which is a function of energy

  real(double) kx,ky,kz,x
  real(double) skx, sky, skz, sk, ki, sq, kf, r1, r2, kh, eih, k12, e12

  integer :: test

  contains

  subroutine gaas_scat_rate  !様々なキャリアネルギーに応じた散乱レートをあらかじめ配列に格納
    implicit none

    real(double) :: absk(2)

    af(1)  = (1d0 - am(1)/am0)**2 / (eg(1)*q)            !非放物線係数 (1.5)
    af(2)  = (1d0 - am(2)/am0)**2 / (eg_l*q)
  
    do ie=1,iemax
      ei = de * real(ie)                                 ![eV]
      absk(1) = smh(1)*sqrt(ei*q * (1 + af(1)*(ei*q)))   !波数の大きさP.112(4.8)
      absk(2) = smh(2)*sqrt(ei*q * (1 + af(2)*(ei*q)))
  
  !---( Gamma 谷の電子の散乱レート )---
  !---( 有極性光学フォノン散乱 )---
      ef = ei - hwo
      if (ef>0d0) then
        gammai = ei*q * (1d0 + af(1)*ei*q)
        gammaf = ef*q * (1d0 + af(1)*ef*q)
        la = (2d0 * (1d0 + af(1)*ei*q) * (1d0 + af(1)*ef*q) + af(1) * (gammai + gammaf))**2
        lb = -2d0 * af(1) * sqrt(gammai * gammaf) * (4d0 * (1d0 + af(1)*ei*q) * (1d0 + af(1)*ef*q) &
           + af(1) * (gammai + gammaf) )
        lc = 4d0 * (1d0 + af(1)*ei*q) * (1d0 + af(1)*ef*q) * (1d0 + 2d0*af(1)*ei*q) * (1d0 + 2d0*af(1)*ef*q)
        f0 = (la * log(abs((sqrt(gammai) + sqrt(gammaf)) / (sqrt(gammai) - sqrt(gammaf)))) + lb) / lc !(2.138)
        swk(1,1,ie) = poe * absk(1) / (ei*q * (1d0 + af(1)*ei*q)) &
                    * (1d0 + 2d0*af(1)*ef*q) * f0!散乱レート(2.137)
      else
        swk(1,1,ie)=0d0
      endif
  
      ef = ei + hwo
      gammai = ei*q * (1d0 + af(1)*ei*q)
      gammaf = ef*q * (1d0 + af(1)*ef*q)
      la = (2d0 * (1d0 + af(1)*ei*q) * (1d0 + af(1)*ef*q) + af(1) * (gammai + gammaf))**2
      lb = -2d0 * af(1) * sqrt(gammai * gammaf) * (4d0 * (1d0 + af(1)*ei*q) * (1d0 + af(1)*ef*q) &
         + af(1) * (gammai + gammaf) )
      lc = 4d0 * (1d0 + af(1)*ei*q) * (1d0 + af(1)*ef*q) * (1d0 + 2d0*af(1)*ei*q) * (1d0 + 2d0*af(1)*ef*q)
      f0 = (la * log(abs((sqrt(gammai) + sqrt(gammaf)) / (sqrt(gammai) - sqrt(gammaf)))) + lb) / lc !(2.138)
      swk(1,2,ie) = swk(1,1,ie) + poa * absk(1) / (ei*q * (1d0 + af(1)*ei*q)) &
                  * (1d0 + 2d0*af(1)*ef*q) * f0

  !---( バンド間フォノン散乱, from Gamma to L)---
  !---(非有極性光学フォノン散乱)---
      ef = ei - hwij + ec(1) - ec(2)
      if (ef > 0d0) then
        fij = ((1d0 + af(1)*ei*q) * (1d0 + af(2)*ef*q)) / ((1d0 + 2d0*af(1)*ei*q) * (1d0 + 2d0*af(2)*ef*q))
        swk(1,3,ie) = swk(1,2,ie) + z2*ope &
                    * sqrt(2d0 * am(2))**3 * sqrt(ef*q * (1d0 + af(2)*ef*q)) / (4d0 * pi**2 * h**2) &
                    * (1d0 + 2d0*af(2)*ef*q) * fij  !P.61(2.135)
      else
        swk(1,3,ie)=swk(1,2,ie)
      endif
  
      ef=ei+hwij+ec(1)-ec(2)
      if (ef>0d0) then
        fij = ((1d0 + af(1)*ei*q) * (1d0 + af(2)*ef*q)) / ((1d0 + 2d0*af(1)*ei*q) * (1d0 + 2d0*af(2)*ef*q))
        swk(1,4,ie) = swk(1,3,ie) + z2*opa &
                    * sqrt(2d0 * am(2))**3 * sqrt(ef*q * (1d0 + af(2)*ef*q)) / (4d0 * pi**2 * h**2) &
                    * (1d0 + 2d0*af(2)*ef*q) * fij  !P.61(2.135) !P.61(2.135)
      else
        swk(1,4,ie)=swk(1,3,ie)
      endif

  !---( 音響フォノン散乱 )---
      swk(1,5,ie)= swk(1,4,ie) + aco * sqrt(2d0 * am(1))**3 * sqrt(ei*q * (1d0 + af(1)*ei*q)) / (4d0 * pi**2 * h**3)&
                 * ((1d0 + af(1)*ei*q)**2 + (af(1)*ei*q)**2/3d0) &
                 / (1d0 + 2d0*af(1)*ei*q) !P.41(2.64)

  !---( 不純物散乱 )---
      sef = sqrt( ei*q * (1d0 + af(1) * ei*q) )
      ak = smh(1) * sef
      qq = qd**2 * (4d0 * ak**2 + qd**2)
      swk(1,6,ie) = swk(1,5,ie) + bimp * (sqrt(2d0 * am(1)) / h)**3 * sqrt(ei*q) / (4d0 * pi**2) / qq  !(2.43)
  
  !---( プラズモン散乱 )---
      qmax = sqrt(q**2 * cden0 / (eps * kb*temp))
      qmin = absk(1) * abs(1d0 - sqrt(1d0 + h*wp(1)/(ei*q)))
      swk(1,7,ie) = swk(1,6,ie) + q**2 * wp(1) * absk(1) / (8d0*pi * eps * (ei*q)) &
                  * 1.0d0/(exp(h*wp(1) / kb*temp)-1.0d0) * log(qmax/qmin)
      
      ef=ei-hwo
      if (ef>0d0) then
        qmax = sqrt(q**2 * cden0 / (eps * kb*temp))
        qmin = absk(1) * abs(1d0 - sqrt(1d0 - h*wp(1)/(ei*q)))
        swk(1,8,ie) = swk(1,7,ie) + q**2 * wp(1) * absk(1) / (8d0*pi * eps * (ei*q)) &
                  * (1.0d0/(exp(h*wp(1) / kb*temp)-1.0d0) + 1d0) * log(qmax/qmin)
      else
        swk(1,8,ie) = swk(1,7,ie)
      end if

  !---( L 谷の電子の散乱レート )---
  !---( 有極性光学フォノン散乱 )---
      ef=ei-hwo
      if (ef>0d0) then
        gammai = ei*q * (1d0 + af(2)*ei*q)
        gammaf = ef*q * (1d0 + af(2)*ef*q)
        la = (2d0 * (1d0 + af(2)*ei*q) * (1d0 + af(2)*ef*q) + af(2) * (gammai + gammaf))**2
        lb = -2d0 * af(2) * sqrt(gammai * gammaf) * (4d0 * (1d0 + af(2)*ei*q) * (1d0 + af(2)*ef*q) &
           + af(2) * (gammai + gammaf) )
        lc = 4d0 * (1d0 + af(2)*ei*q) * (1d0 + af(2)*ef*q) * (1d0 + 2d0*af(2)*ei*q) * (1d0 + 2d0*af(2)*ef*q)
        f0 = (la * log(abs((sqrt(gammai) + sqrt(gammaf)) / (sqrt(gammai) - sqrt(gammaf)))) + lb) / lc !(2.138)
        swk(2,1,ie) = poe * absk(2) / (ei*q * (1d0 + af(2)*ei*q)) &
                    * (1d0 + 2d0*af(2)*ef*q) * f0!散乱レート(2.137)
      else
        swk(2,1,ie)=0d0
      endif
  
      ef = ei + hwo
      gammai = ei*q * (1d0 + af(2)*ei*q)
      gammaf = ef*q * (1d0 + af(2)*ef*q)
      la = (2d0 * (1d0 + af(2)*ei*q) * (1d0 + af(2)*ef*q) + af(2) * (gammai + gammaf))**2
      lb = -2d0 * af(2) * sqrt(gammai * gammaf) * (4d0 * (1d0 + af(2)*ei*q) * (1d0 + af(2)*ef*q) &
         + af(2) * (gammai + gammaf) )
      lc = 4d0 * (1d0 + af(2)*ei*q) * (1d0 + af(2)*ef*q) * (1d0 + 2d0*af(2)*ei*q) * (1d0 + 2d0*af(2)*ef*q)
      f0 = (la * log(abs((sqrt(gammai) + sqrt(gammaf)) / (sqrt(gammai) - sqrt(gammaf)))) + lb) / lc !(2.138)
      swk(2,2,ie) = swk(2,1,ie) + poa * absk(2) / (ei*q * (1d0 + af(2)*ei*q)) &
                  * (1d0 + 2d0*af(2)*ef*q) * f0!散乱レート(2.137)

  !---( バンド間フォノン散乱, from L to L )---
  !---(非有極性光学フォノン散乱)---
      ef = ei-hwe
      if (ef > 0d0) then
        fij = ((1d0 + af(2)*ei*q) * (1d0 + af(2)*ef*q)) / ((1d0 + 2d0*af(2)*ei*q) * (1d0 + 2d0*af(2)*ef*q))
        swk(2,3,ie) = swk(2,2,ie) + (z2-1d0)*ope &
                    * sqrt(2d0 * am(2))**3 * sqrt(ef*q * (1d0 + af(2)*ef*q)) / (4d0 * pi**2 * h**2) &
                    * (1d0 + 2d0*af(2)*ef*q) * fij  !P.61(2.135)
      else
        swk(2,3,ie)=swk(2,2,ie)
      endif
  
      ef=ei+hwe
      fij = ((1d0 + af(2)*ei*q) * (1d0 + af(2)*ef*q)) / ((1d0 + 2d0*af(2)*ei*q) * (1d0 + 2d0*af(2)*ef*q))
      swk(2,4,ie) = swk(2,3,ie) + (z2-1d0)*opa &
                  * sqrt(2d0 * am(2))**3 * sqrt(ef*q * (1d0 + af(2)*ef*q)) / (4d0 * pi**2 * h**2) &
                  * (1d0 + 2d0*af(2)*ef*q) * fij  !P.61(2.135)

  !---( バンド間フォノン散乱, from L to Gamma )---
  !---(非有極性光学フォノン散乱)---
      ef=ei-hwij+ec(2)-ec(1)
      if (ef>0d0) then
        fij = ((1d0 + af(2)*ei*q) * (1d0 + af(1)*ef*q)) / ((1d0 + 2d0*af(2)*ei*q) * (1d0 + 2d0*af(1)*ef*q))
        swk(2,5,ie) = swk(2,4,ie) + ope &
                    * sqrt(2d0 * am(1))**3 * sqrt(ef*q * (1d0 + af(1)*ef*q)) / (4d0 * pi**2 * h**2) &
                    * (1d0 + 2d0*af(1)*ef*q) * fij 
      else
        swk(2,5,ie)=swk(2,4,ie)
      endif
  
      ef=ei+hwij+ec(2)-ec(1)
      if (ef>0d0) then
        fij = ((1d0 + af(2)*ei*q) * (1d0 + af(1)*ef*q)) / ((1d0 + 2d0*af(2)*ei*q) * (1d0 + 2d0*af(1)*ef*q))
        swk(2,6,ie) = swk(2,5,ie) + opa &
                    * sqrt(2d0 * am(1))**3 * sqrt(ef*q * (1d0 + af(1)*ef*q)) / (4d0 * pi**2 * h**2) &
                    * (1d0 + 2d0*af(1)*ef*q) * fij 
      else
        swk(2,6,ie)=swk(2,5,ie)
      endif

  !---( 音響フォノン散乱 )---
      swk(2,7,ie)= swk(1,6,ie) + aco * sqrt(2d0 * am(2))**3 * sqrt(ei*q * (1d0 + af(2)*ei*q)) / (4d0 * pi**2 * h**3)&
                 * ((1d0 + af(2)*ei*q)**2 + (af(2)*ei*q)**2/3d0) &
                 / (1d0 + 2d0*af(2)*ei*q) 

  !---( 不純物散乱 )---
      sef = sqrt(ei*q * (1d0 + af(2) * ei*q) )
      ak = smh(2) * sef
      qq = qd**2 * (4d0 * ak**2 + qd**2)
      swk(2,8,ie) = swk(2,7,ie) + bimp * (sqrt(2d0 * am(2)) / h)**3 * sqrt(ei*q) / (4d0 * pi**2) / qq  !(2.43)

  !---( プラズモン散乱 )---
      qmax = sqrt(q**2 * cden0 / (eps * kb*temp))
      qmin = absk(2) * abs(1d0 - sqrt(1d0 + h*wp(2)/(ei*q)))
      swk(2,9,ie) = swk(2,8,ie) + q**2 * wp(2) * absk(2) / (8d0*pi * eps * (ei*q)) &
                  * 1.0d0/(exp(h*wp(2) / kb*temp)-1.0d0) * log(qmax/qmin)
      
      ef=ei-hwo
      if (ef>0d0) then
        qmax = sqrt(q**2 * cden0 / (eps * kb*temp))
        qmin = absk(2) * abs(1d0 - sqrt(1d0 - h*wp(2)/(ei*q)))
        swk(2,10,ie) = swk(2,9,ie) + q**2 * wp(2) * absk(2) / (8d0*pi * eps * (ei*q)) &
                  * (1.0d0/(exp(h*wp(2) / kb*temp)-1.0d0) + 1d0) * log(qmax/qmin)
      else
        swk(2,10,ie) = swk(2,9,ie)
      end if
    enddo
  
  !---( 散乱レートの総和の計算 )---

    test = 10

    write(*,'(e10.3,",",e10.3,",",e10.3,",",e10.3,",",e10.3)')&
         swk(1,2,test), swk(1,4,test)-swk(1,2,test), swk(1,5,test)-swk(1,4,test),&
         swk(1,6,test)-swk(1,5,test), swk(1,8,test)-swk(1,6,test)

    write(*,'(e10.3,",",e10.3,",",e10.3,",",e10.3,",",e10.3)')&
         swk(2,2,test), swk(2,6,test)-swk(2,2,test), swk(2,7,test)-swk(2,6,test),&
         swk(2,8,test)-swk(2,7,test), swk(2,10,test)-swk(2,8,test)

    gamma=swk(1,8,1)
    do ie=1,iemax
      if (swk(1,8,ie)>gamma) then 
        gamma=swk(1,8,ie)
      end if
      if (swk(2,10,ie)>gamma) then 
        gamma=swk(2,10,ie)
      end if
    enddo
    do i=1,8
      do ie=1,iemax
        swk(1,i,ie)=swk(1,i,ie)/gamma
      enddo
    enddo
    do i=1,10
      do ie=1,iemax
        swk(2,i,ie)=swk(2,i,ie)/gamma
      enddo
    enddo
    return
  end subroutine gaas_scat_rate



  subroutine scat(pt, iv) !キャリアエネルギーに応じた散乱機構の選択
    real(double), intent(inout) :: pt(11)
    integer, intent(inout) :: iv

    kx = pt(1)
    ky = pt(2)
    kz = pt(3)
    x  = pt(5)

    skx = kx**2
    sky = ky**2
    skz = kz**2
    sk  = skx + sky + skz
    if (sk==0d0) return
    ki = sqrt(sk)                                  !キャリアの波数の絶対値
    sq = sqrt(1d0 + 4d0 * af(iv) * hhm(iv) * sk)
    ei = (sq - 1d0) / (2d0 * af(iv)) / q           !キャリアのエネルギー [eV]
    ie = max(1,min(int(ei/de)+1,iemax))            !整数化したエネルギーieへの置き換え

    if (iv==1) then
      r1=rnd()

      if (r1<=swk(iv,1,ie)) then
        ef=ei-hwo
        if (ef<=0d0) return
        call fin20(kx,ky,kz)
      elseif (r1<=swk(iv,2,ie)) then
        ef=ei+hwo
        call fin20(kx,ky,kz)
      elseif (r1<=swk(iv,3,ie)) then
        ef=ei-hwij+ec(1)-ec(2)
        if (ef<=0d0) return
        iv=2
        call fin40(kx,ky,kz,iv)
      elseif (r1<=swk(iv,4,ie)) then
        ef=ei+hwij+ec(1)-ec(2)
        if (ef<=0d0) return
        iv=2
        call fin40(kx,ky,kz,iv)
      elseif (r1<=swk(iv,5,ie)) then
        ef=ei
        kf=ki
        call fin40(kx,ky,kz,iv)
      elseif (r1<=swk(iv,6,ie)) then
        ef=ei
        r2=rnd()
        cb=1d0-2d0*r2/(1d0+(1d0-r2)*4d0*sk/qd**2)
        kf=ki
        call fin30(kx,ky,kz)
      elseif (r1<=swk(iv,7,ie)) then
        ef=ei+hwo
        call fin20(kx,ky,kz)
      elseif (r1<=swk(iv,8,ie)) then
        ef=ei+hwo
        call fin20(kx,ky,kz)
      endif
    else
      r1=rnd()

      if (r1<=swk(iv,1,ie)) then
        ef=ei-hwo
        if (ef<=0d0) return
        call fin20(kx,ky,kz)
      elseif (r1<=swk(2,2,ie)) then
        ef=ei+hwo
        call fin20(kx,ky,kz)
      elseif (r1<=swk(2,3,ie)) then
        ef=ei-hwe
        if (ef<=0d0) return
        call fin40(kx,ky,kz,iv)
      elseif (r1<=swk(2,4,ie)) then
        ef=ei+hwe
        call fin40(kx,ky,kz,iv)
      elseif (r1<=swk(2,5,ie)) then
        ef=ei-hwij+ec(2)-ec(1)
        if (ef<=0d0) return
        iv=1
        call fin40(kx,ky,kz,iv)
      elseif (r1<=swk(2,6,ie)) then
        ef=ei+hwij+ec(2)-ec(1)
        if (ef<=0d0) return
        iv=1
        call fin40(kx,ky,kz,iv)
      elseif (r1<=swk(2,7,ie)) then
        ef=ei
        kf=ki
        call fin40(kx,ky,kz,iv)
      elseif (r1<=swk(2,8,ie)) then
        ef=ei
        r2=rnd()
        cb=1d0-2d0*r2/(1d0+(1d0-r2)*4d0*sk/qd**2)
        kf=ki
        call fin30(kx,ky,kz)
      elseif (r1<=swk(iv,9,ie)) then
        ef=ei+hwo
        call fin20(kx,ky,kz)
      elseif (r1<=swk(iv,10,ie)) then
        ef=ei+hwo
        call fin20(kx,ky,kz)
      endif
    endif
    pt(1) = kx
    pt(2) = ky
    pt(3) = kz

    t1 = t_s
    t_s = t1 - log(rnd()) / gamma

    return
  end subroutine scat



  subroutine fin20(kx,ky,kz)
    real(double), intent(inout) :: kx,ky,kz
    real(double) sb,fai,cf,sf,skk,a11,a12,a13,a21,a22,a23,a32,a33,x1,x2,x3,f
    !---( 散乱後の状態決定 )---
    f   = 2d0 * sqrt(ei * ef) / (sqrt(ei) - sqrt(ef))**2
    if (f<=0d0) return
    cb  = (1d0 + f - (1d0 + 2d0 * f)**rnd()) / f

    sb  = sqrt(1d0-cb*cb)
    fai = 2d0*pi*rnd()
    cf  = cos(fai)
    sf  = sin(fai)
    skk = sqrt(skx+sky)
    a11 = ky/skk
    a12 = kx*kz/(skk*ki)
    a13 = kx/ki
    a21 =-kx/skk
    a22 = ky*kz/(skk*ki)
    a23 = ky/ki
    a32 =-skk/ki
    a33 = kz/ki
    x1  = kf*sb*cf
    x2  = kf*sb*sf
    x3  = kf*cb
    kx  = a11*x1+a12*x2+a13*x3
    ky  = a21*x1+a22*x2+a23*x3
    kz  =        a32*x2+a33*x3
    return
  end subroutine fin20



  subroutine fin30(kx,ky,kz)  !不純物散乱の散乱後の状態決定
    real(double), intent(inout) :: kx,ky,kz
    real(double) sb,fai,cf,sf,skk,a11,a12,a13,a21,a22,a23,a32,a33,x1,x2,x3
    !---( 散乱後の状態決定 )---
    sb  = sqrt(1d0-cb*cb)
    fai = 2d0*pi*rnd()
    cf  = cos(fai)
    sf  = sin(fai)
    skk = sqrt(skx+sky)
    a11 = ky/skk
    a12 = kx*kz/(skk*ki)
    a13 = kx/ki
    a21 =-kx/skk
    a22 = ky*kz/(skk*ki)
    a23 = ky/ki
    a32 =-skk/ki
    a33 = kz/ki
    x1  = kf*sb*cf
    x2  = kf*sb*sf
    x3  = kf*cb
    kx  = a11*x1+a12*x2+a13*x3
    ky  = a21*x1+a22*x2+a23*x3
    kz  =        a32*x2+a33*x3
    return
  end subroutine fin30



  subroutine fin40(kx,ky,kz,iv)
    integer, intent(in) :: iv
    real(double), intent(inout) :: kx,ky,kz
    real(double) cs,sn,fai
    kf  = smh(iv)*sqrt(ef*q*(1d0+af(iv)*ef*q))
    cs  = 1d0-2d0*rnd()
    sn  = sqrt(1d0-cs*cs)
    fai = 2d0*pi*rnd()
    kx  = kf*cs
    ky  = kf*sn*cos(fai)
    kz  = kf*sn*sin(fai)
    return
  end subroutine fin40

end module scatter
