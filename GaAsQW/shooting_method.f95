module shooting_method_mod
    use config_mod
    implicit none

    !This module gives band gaps and ground state energy of electron in a quantum well 

    integer ::  max_ix, sign
    real(double) :: eg0(2),eg(2), alpha(2), beta(2), m_semi(2), a, bowing, eg_a,&
                    shoot_x, shoot_dx, shoot_en, shoot_den, well_l, well_r, offset,&
                    max_psi, max_x, allow, allow_param,gstate_en
    real(double), allocatable :: f(:), shoot_psi(:)


    ! Shooting Method to calculate the ground state enrgey of a electron in a quantum well
    !
    !   __________________   ____________________
    !                     | | 
    !                     | |
    !    semiconductor 2  | |  semiconductor 2       offset = (band gap of 2 - band gap of 1)  / 2
    !                     | |
    !                     | |
    !                     |_|
    !                      ^
    !                semiconductor 1

    contains

    subroutine varshni  !This subroutine gives band gap energy using Varshni's law
        implicit none
    
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 参考文献 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !   I. Vurgaftman, J. R. Meyer, L. R. Ram-Mohan, J. Appl. Phys. 89 (11), 5815 (2001)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
        eg0(1) = 1.519d0                                     !0 K におけるガンマ谷のバンドギャプ [eV]
        alpha(1) = 5.41d-4                                   !varshni則の定数 [eV/K]
        beta(1) = 204d0                                      !varshni則の定数 [eV/K]
        eg(1) = eg0(1) - alpha(1) * temp**2 / (temp + beta(1)) !T K におけるバンドギャップ [eV]
        bgap = eg(1)
    
        eg0(2) = 3.099d0                                     !Al(x)Ga(1-x)Asのガンマ谷における計算
        alpha(2) = 8.85d-4
        beta(2) = 530d0
        eg(2) = eg0(2) - alpha(2) * temp**2 / (temp + beta(2))
    
        bowing = -0.127 + 1.31 * a                           !AlGaAsのbowing parameter
    
        eg_a = a*eg(2) + (1-a)*eg(1) - a*(1-a)*bowing        !Al(a)Ga(1-a)Asのバンドギャップ

        eg0_l = 1.815d0                                      !GaAsのL谷における計算
        alpha_l = 6.05d-4
        beta_l = 204
        eg_l = eg0_l - alpha_l * temp**2 / (temp + beta_l) 

        ec(1) = 0d0                                          !ガンマ谷の伝導体の底のエネルギー（eV）
        ec(2) = (eg_l - eg(1)) / 2d0                         !L谷の伝導体の底のエネルギー(eV)
    
        return
    end subroutine varshni



    subroutine shooting_method  !This subroutine gives Zero Point Eergy in the quantum well

        integer :: i, j


    open(unit=30, file='output_psi_graph')
    open(unit=31, file='output_ground_state_energy')

    a = 0.3                                       !Al(a)Ga(1-a)As
    max_x  = 62.2d-9                              !右端のx座標 [m]
    shoot_dx = 0.1d-9                             !刻み幅 [m]
    max_ix = int(max_x / shoot_dx + 0.5)            !節のの数(1~)
    well_l = 30d-9                                !井戸の左端のx座標 [m]
    iwell_l = int(well_l/shoot_dx + 0.5)
    well_r = well_l + wd                          !井戸の右端のx座標 [m]
    iwell_r = int(well_r/shoot_dx + 0.5)
    allow_param = 0.001d0                         !収束判定のパラメータ

    allocate(f(max_ix))                           !波動関数の傾き
    allocate(shoot_psi(max_ix))                   !波動関数
    
    sign = 1                                      !波動関数の正負を調べるのに使う
    f(1) = 0.10d0                                 !傾きの初期値
    shoot_psi(1) = 0d0                            !波動関数の初期値
    shoot_en = 1d-3                               !エネルギー[eV]
    shoot_den = 1d-3                              !エネルギーの刻みはば[eV]

    m_semi(1) = am0 * 0.067d0                     !GaAsのガンマ谷での有効質量
    m_semi(2) = am0 * (a*0.15 + (1-a)*0.067d0 )   !Al(a)Ga(1-a)Asのガンマ谷での有効質量----線形性を仮定

    call varshni
    offset = ( eg_a - eg(1) ) / 2d0         !バンドオフセット（井戸の壁の高さ）

    do i = 1, MAX_LOOP

        if(mod(i, 1000000) == 0) then
            write(*,'("loop =",i12,", den =",e12.5)') i, shoot_den
        end if

        max_psi = 0d0

        do j = 2, max_ix
            if (j< iwell_l .or. j > iwell_r) then
                f(j) = 2 * m_semi(2) * (offset - shoot_en)*q * shoot_psi(j - 1) * shoot_dx / h**2 + f(j - 1)
                shoot_psi(j) = f(j - 1) * shoot_dx + shoot_psi(j - 1)
            else if (j == iwell_l) then !境界条件
                f(j) = am(1)/am(2) &
                     * 2 * m_semi(2) * (offset - shoot_en)*q * shoot_psi(j - 1) * shoot_dx / h**2 + f(j - 1)
                shoot_psi(j) = f(j - 1) * shoot_dx + shoot_psi(j - 1)
            else if (j == iwell_r) then !境界条件
                f(j) = am(2)/am(1) &
                       *2 * m_semi(1) * (-shoot_en)*q * shoot_psi(j - 1) * shoot_dx / h**2 + f(j - 1)
                shoot_psi(j) = f(j - 1) * shoot_dx + shoot_psi(j - 1)
            else
                f(j) = 2 * m_semi(1) * (-shoot_en)*q * shoot_psi(j - 1) * shoot_dx / h**2 + f(j - 1)
                shoot_psi(j) = f(j - 1) * shoot_dx + shoot_psi(j - 1)
            end if
        end do

        do j = 2, max_ix
            if (abs(shoot_psi(j)) > max_psi) then
                max_psi = abs(shoot_psi(j))
            end if
        end do

        allow = max_psi * allow_param

        if (abs(shoot_psi(max_ix)) <= allow) then
            do j = 1, max_ix
                shoot_x = shoot_dx * j
                write(30,'(e12.5,", ",e12.5)') shoot_x, shoot_psi(j)
            end do
            write(31,'("ground satate energy =",e12.5,"[eV]")') shoot_en 
            close(30)
            close(31)
            gstate_en = shoot_en
            return
        end if

        if (sign * shoot_psi(max_ix) < 0) then
            sign = sign * (-1)
            shoot_en = shoot_en - shoot_den
            shoot_den = shoot_den * 0.1d0
        else
            shoot_en = shoot_en + shoot_den
        end if

        if (i == MAX_LOOP) then
            write(*,*) 'Error : the loop in shooting method didnt end!'
            close(30)
            close(31)
            close(10)
            close(11)
            stop
        end if

    end do

    close(30)
    close(31)

    end subroutine shooting_method

end module shooting_method_mod







    