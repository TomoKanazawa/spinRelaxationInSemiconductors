program main
    use progress
    use field
    implicit none

    ! Monte Carlo Simulation of spin relaxation along z axis 
    ! in AlGaAs/GaAs quantum well
    !
    !               x axis
    !           _______________
    !          |               |                  x axis (y axis)
    !          |               |              _______________
    !  y sxis  | semiconductor |      z axis |_______________|
    !          |               |       
    !          |_______________|
    !


    open(unit=10,file='output_individual')
    open(unit=11,file='output_average')

    call config
    call shooting_method
    call initialize_carrier

    allocate(ipt(pmax),pt(pmax,11),cden(mx1,my1),an(mx1,my1), &
    swk(2,10,iemax),dpef(jtl),eyef(jtl),bapef(jtl),fx(mx1,my1),fy(mx1,my1))

    fx      = 0d0  !電界のX成分
    fy      = 0d0  !電界のy成分
    an      = 0d0  !格子中のキャリア数
    cden    = 0d0  !キャリア濃度密度
    spin    = 0d0
    dpef    = 0d0  !dt秒間のスピン緩和の内のDP効果の寄与の割合
    eyef    = 0d0
    bapef   = 0d0
    dpav    = 0d0
    eyav    = 0d0
    bapav   = 0d0
    dpall   = 0d0 !時刻tまでの全スピン緩和の内のDP効果の寄与の割合
    eyall   = 0d0
    bapall  = 0d0
    momentum= 0d0
    energy  = 0d0
    t       = 0d0


    call gaas_scat_rate
    call initialize_pt  !ドープ濃度に応じたキャリアの振り分け

    write(*,*) "1 =", bimp/2d0 * (sqrt(2d0 * am(1)) / h)**3 * sqrt(0.001d0*q) / (4d0 * pi**2)
    write(*,*) "2 =", 2d0*(smh(1)*sqrt(0.001d0*q * (1 + af(1)*(0.001d0*q))))**2
    write(*,*) "3= ", qd**2
    write(*,*) "swk = ", swk(1,6,1)*gamma

    write(*,'(/,"Temperature = ",f5.1,"[K]",/,"Band gap energy =",f8.4,"[eV]  Offset =",f8.4,"[eV]",/, &
            &"Energy gap btween gamma and L valleys =",f8.4,"[eV]",/, &
            &"Ground satate energy =",f8.4,"[eV]",/,"Carrier dentisy =",e12.4,"[/m-3]",/, &
            &"Carrier number per particle =",f8.5)') &
            temp, bgap, offset, ec(2), shoot_en,cden0,epp

    do i = 1,inum
        spin = spin + pt(i,7)/ (h/2d0 * inum) * 100d0
        energy = energy + &
                 (-1d0+sqrt(1d0+2d0*af(ipt(i))*h**2*(pt(i,1)**2+pt(i,2)**2+pt(i,3)**2) / am(ipt(i)))) &
                 / (2d0 * af(ipt(i))) / q / inum  !キャリアの平均エネルギー[eV]
        momentum = momentum + h * pt(i,1)         !代表点のx方向の運動量の総和
    end do


    write(*,'(/,"Time =",f6.2,"[ps],   inum =",i8,/, &
            &"DP mechanism rate =",f9.4,"[%],  EY mechanism rate =",f9.4,"[%]",/, &
            &"BAP mechanism rate =",f9.4,"[%],  Spin =",f8.4,"[%]",/,&
            &"Momentum =",e13.5,"[kg·m/s],   Average energy =",f10.4,"[meV]")') &
           0d0, inum, 0d0,0d0,0d0, spin, momentum, energy*1d3


    write(10,'(e12.5,", ",e12.5,", ",e12.5,", ",e12.5,", ",e12.5,", ",e12.5,", ",e12.5,", ",i8)')&
           t*1.0e12, spin, &
           0d0,0d0,0d0, &
           momentum, energy, inum


    do jt = 1, jtl

        spin = 0d0
        momentum = 0d0
        energy = 0d0

        call emcd                                  !粒子のドリフトと散乱をΔtだけ進める関数
        call renew                                 !再結合した粒子(キャリア)を消去
        call charge                                !電荷分布を計算
        call poisson                               !ポアソン方程式ををといてポテンシャルを求める

        t = dt * real(jt)

        dpall  = dpall + dpef(jt)
        eyall  = eyall + eyef(jt)
        bapall = bapall + bapef(jt)

        do i = 1,inum
            spin = spin + pt(i,7) / (h/2.0d0 * inum) * 100d0
            energy = energy + &
                    (-1d0+sqrt(1d0+2d0*af(ipt(i))*h**2*(pt(i,1)**2+pt(i,2)**2+pt(i,3)**2) / am(ipt(i)))) &
                    / (2d0 * af(ipt(i))) / q / inum
            momentum = momentum + h * pt(i,1)   !代表点のx方向の運動量の総和
        end do

        if (mod(jt,100) == 0) then
            write(*,'(/,"Time =",f6.2,"[ps],   inum =",i8,/,&
                    &"DP mechanism rate =",f9.4,"[%],   EY mechanism rate =",f9.4,"[%]",/, &
                    &"BAP mechanism rate =",f9.4,"[%],   Spin =",f8.4,"[%]",/,&
                    &"Momentum =",e13.5,"[kg·m/s],   Average energy =",f10.4,"[meV]")') &
                    jt*dt*1d12, inum, dpall/(dpall+eyall+bapall)*100,eyall/(dpall+eyall+bapall)*100, &
                    bapall/(dpall+eyall+bapall)*100,spin, momentum, energy*1d3
        end if 


        if (dpall/=0d0 .and. eyall/=0d0 .and. bapall/=0d0) then
            write(10,'(e12.5,", ",e12.5,", ",e12.5,", ",e12.5,", ",e12.5,", ",e12.5,", ",e12.5,", ",i8)')&
                    t*1.0e12, spin, &
                    dpall/(dpall+eyall+bapall)*100,eyall/(dpall+eyall+bapall)*100,&
                    bapall/(dpall+eyall+bapall)*100, &
                    momentum, energy, inum
        else
            write(10,'(e12.5,", ",e12.5,", ",e12.5,", ",e12.5,", ",e12.5,", ",e12.5,", ",e12.5,", ",i8)')&
                    t*1.0e12, spin, &
                    0d0,0d0,&
                    0d0, &
                    momentum, energy, inum
        end if
              

        if (inum <= 0) then
            write(*,*) 'inum <= 0. inmu =' ,inum
            close(10)
            close(11)
            stop
        end if

        if (int(spin) > 100e0 .or. 0e0 > int(spin)) then !単精度使用はまるめ誤差対策のため
            write(*,*) 'Error : spin value is not proper!' ,spin, jt
            close(10)
            close(11)
            stop
        end if

    end do 

    do jt = 1, jtl 
        if (jt > range .or. jt <= jtl - range) then
            dpav  = sum(dpef(jt-range:jt+range)) / (2d0 * range + 1d0)
            eyav  = sum(eyef(jt-range:jt+range)) / (2d0 * range + 1d0)
            bapav = sum(bapef(jt-range:jt+range)) / (2d0 * range + 1d0)
        else if (jt <=5) then
            dpav  = sum(dpef(1:jt+range)) / (jt + range)
            eyav  = sum(eyef(1:jt+range)) / (jt + range)
            bapav = sum(bapef(1:jt+range)) / (jt + range)
        else
            dpav  = sum(dpef(jt-range:jtl)) / (jtl - jt + 1d0 + range)
            eyav  = sum(eyef(jt-range:jtl)) / (jtl - jt + 1d0 + range)
            bapav = sum(bapef(jt-range:jtl)) / (jtl - jt + 1d0 + range)
        end if 

        if (.not. (dpav==0d0 .and. eyav==0d0 .and. bapav==0d0)) then
            write(11,'(e12.5,", ",e12.5,", ",e12.5,", ",e12.5)') &
                    jt*dt*1d12, &
                    dpav/(dpav+eyav+bapav)*100,eyav/(dpav+eyav+bapav)*100,bapav/(dpav+eyav+bapav)*100
        else
            write(11,'(e12.5,", ",e12.5,", ",e12.5,", ",e12.5)') &
                    jt*dt*1d12, &
                    0d0,0d0,0d0
        end if

    end do


    close(10)
    close(11)

    write(*,'("DP mechanism rate = ",f9.4,"[%]",/,"EY mechanism rate = ",f9.4,"[%]",/,&
            "BAP mechanism rate = ",f9.4,"[%]")') &
            dpall/(dpall+eyall+bapall)*100,eyall/(dpall+eyall+bapall)*100,bapall/(dpall+eyall+bapall)*100

    stop
    end program main