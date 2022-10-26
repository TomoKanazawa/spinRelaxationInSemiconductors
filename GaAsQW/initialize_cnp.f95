module initialize_cnp_mod
  use shooting_method_mod
  implicit none

  !This module initializes carriers and particles
  !Note that the particles are representatives of the carriers for reducing caliculation load

  real(double) cb,cf,fai,sb,sf
  integer j,n,iv,m,npij

  contains

  subroutine initialize_carrier  !キャリア濃度を計算
    implicit none
  
    real(double) :: en
    integer :: loop,n
    en = 0d0
    n = 0
    efel  = ec(1) - bgap/2                                      !ガンマ点のフェルミエネルギー
  
    puls_en = pw / hz &                                         !パルスのエネルギー [J]
              * (1 - ((1-n1)/(1+n1))**2 ) &                     !AlGaAsと真空の境界での表面反射（１回）
              * (1 - ((n1-n2)/(n1+n2))**2)**7                   !AlGaAsとGaAsの境界での表面反射（７回）
    sim_en = puls_en / (pi * lay_r**2) * (xmax*ymax)            !シミュレーションでは小さい範囲のみを考える [J]
  
    do loop = 1, MAX_LOOP
      en  = en + gstate_en*q + bgap*q
      n = n +1
      if (en >= sim_en) exit
      if(n >= pmax) then
        write(*,*) 'Error ; the exited particle number exceeds pmax!'
        close(10)
        close(11)
        stop
      end if
      if (loop == MAX_LOOP) then
        write(*,*) 'Error : the loop in initialize_carrier didnt end!'
        close(10)
        close(11)
        stop
      end if
    end do
  
    ex_cden = n / (xmax * ymax * wd)                             !励起されるキャリア数密度
  
  
    eq_cden = 2d0 * sqrt(2d0*pi * am(1) * kb*temp / h**2)**3 &
              * exp(-bgap*q/(2d0 * kb*temp))                     !熱平衡時のキャリア密度
  
    cden0 = ex_cden + eq_cden                                    !キャリア（電子）密度の初期値
    exr = ex_cden / cden0                                        !励起により発生したキャリアの割合
    eqr = eq_cden / cden0                                        !熱平衡時に伝導体にすでにあるキャリアの割合
  
  
    epp = cden0 * dx*dy*wd / real(ncon)                          !粒子1個当たりのキャリアの個数
    qd   = sqrt(q**2 * cden0 / (eps * kb*temp))                  !不純物散乱のためのパラメータ
    bimp   = 2d0 * pi * cden0 * q**4 / (h * eps**2)              !不純物散乱のためのパラメータ
    wp(1)     = sqrt(q**2 * cden0 / (eps * am(1)))               !プラズマ周波数
    wp(2)     = sqrt(q**2 * cden0 / (eps * am(2)))
  
    return
  end subroutine initialize_carrier


    
  subroutine initialize_pt  !キャリアの代表点である粒子の初期設定
    implicit none
    
    real(double) :: rnd_value
    call initialize_random  !乱数生成機構の初期設定
      
    pt = 0d0
    
    !gamma_dp = 11.0d0 * q * 1d-30  !Refernce : Spin-orbit coupling in bulk GaAs
        
    n=0
    do i=1,mx1
      do j=1,my1
        npij = int(cden0*dx*dy*wd/epp+0.5d0)                    !格子中の粒子の数(四捨五入)
        if ((i==1).or.(i==mx1)) then
          npij=npij/2
        end if
        if ((j==1).or.(j==my1)) then
          npij=npij/2
        end if
        if (npij==0) cycle
          do m=1,npij
            n = n + 1
            if (n > pmax) then
              write(*,*) 'Error : number of particles exceeds', pmax
              close(10)
              close(11)
              stop                                              !粒子数が足りないとプログラムを終了する
            endif
      
            iv  = 1                                             !バンドを識別する添え字
            ei  = - kb*temp * log(rnd()) / q                    !キャリアのエネルギー[eV](4.4)
            ak  = smh(iv)*sqrt(ei*q * (1d0 + af(iv) * ei*q))    !キャリアの波数kの大きさ(4.5)
            cb  = 1d0-2d0*rnd()                                 !波数ベクトルの天頂角のcos
            sb  = sqrt(1d0-cb*cb)                               !波数ベクトルの天頂角のsin
            fai = 2d0*pi*rnd()                                  !波数ベクトルの方位角
            sf  = sin(fai)
            cf  = cos(fai)

            pt(n,1) = ak*sb*cf                                 !波数ベクトルのx成分
            pt(n,2) = ak*sb*sf                                 !波数ベクトルのy成分
            pt(n,3) = ak*cb                                    !波数ベクトルのz成分
            !pt(n,1) = ak                                       !運動量緩和時間を調べる時以外はコメントアウト
            !pt(n,2) = 0d0
            !pt(n,3) = 0d0

            pt(n,4) = -log(rnd()) / gamma                       !キャリアの自由時間(3.5)

            pt(n,5) = dx*(real(i-1)+rnd()-0.5d0)                !キャリアのx座標
            pt(n,6) = dy*(real(j-1)+rnd()-0.5d0)                !キャリアのy座標

            pt(n,7) = h/2.0d0                                   !電子のスピン角運動量(z軸方向)    

            rnd_value = rnd()
            if( rnd_value - eqr >= 0 ) then
              pt(n,8) = -tr * log((rnd_value - eqr) / exr)      !再結合時間 [s]
            else
              pt(n,8) = 1000                                    !熱平衡時のキャリア密度相当のキャリアは再結合しない
            end if 

            eta_dp = del*q / (bgap*q + del*q)
            alpha_dp = 4d0*eta_dp/sqrt(3d0-eta_dp) * am(iv)/am0
            gamma_dp = alpha_dp * h**3 / sqrt(2d0 * am(iv)**3 * bgap*q)
            pt(n,9) = -h**8 &
                      / (16d0 * kb*temp * am(iv)**3 * (gamma_dp * gstate_en*q)**2 * t_p ) &
                      * log(rnd())                              !DP効果のスピン緩和時間

            pt(n,10) = -9d0/8d0 / (del*q/(bgap*q + del*q))**2 &
                       / (1d0-am(iv)/am0)**2 / (gstate_en*q * kb*temp) * (bgap*q)**2 * t_p &
                       * log(rnd())                             !EY効果のスピン緩和時間

            pt(n,11) = - h**2 / (param_bap * cden0 * lattice**6 * sqrt(am(iv))**3) &
                       * log(rnd())                             !BAP効果のスピン緩和時間 
                       
            if (i==1)   pt(n,5)=dx*0.5d0*rnd()
            if (j==1)   pt(n,6)=dy*0.5d0*rnd()
            if (i==mx1) pt(n,5)=xmax-dx*0.5d0*rnd()
            if (j==my1) pt(n,6)=ymax-dy*0.5d0*rnd()
            ipt(n)  = iv
          enddo
      enddo
    enddo
    inum=n
    do n=inum+1,pmax
      ipt(n)=9
    enddo
      
    return
  end subroutine initialize_pt

end module initialize_cnp_mod