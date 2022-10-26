module config_mod
  use initialize_random_mod
  implicit none

  !This module gives basic configuration

  real(double), parameter :: pi  = 3.14159d0   !円周率
  real(double), parameter :: q   = 1.60219d-19 !電気素量
  real(double), parameter :: h   = 1.05459d-34 !ディラック定数
  real(double), parameter :: kb  = 1.38066d-23 !ボルツマン定数
  real(double), parameter :: ep0 = 8.85419d-12 !真空の誘電率
  real(double), parameter :: am0 = 9.10953d-31 !電子の質量
  real(double), parameter :: mol = 6.022d23    !アボガドロ定数

  real(double):: dx,dy,dt,temp,pw,ec(2),am(2),efel,hz,lay_r,bgap,param_bap,&
             cden0,epp,xmax,ymax,puls_en,sim_en,ex_cden,eq_cden,hwo,hwij,lattice,&
             hwe,de,gamma,smh(2),hhm(2),af(2),hm(2),del,no,nij,ne,aco,ak,qd,&
             bimp,cl,da,dij,ef,opa,ope,poa,poe,qmax,qmin,rou,sef,sv,we,wij,wo,&
             z2,ei,eps,eqe,ep,epf,qq,tr,exr,eqr,wd,n1,n2,t1,tau,tdt,t_s,t_bap,&
             t_dp,t_r,alpha_dp,eta_dp,gamma_dp,t_p,eg0_l,alpha_l,beta_l,eg_l,t_ey,&
             t,spin,momentum,energy,iwell_l,iwell_r,wp(2),amh,ameh(2),fij,f0,la,lb,lc,&
             gae,ase,ionz,gaw,asw,gaasd,gaasn,dpall,eyall,bapall,dpav,eyav,bapav,gammai,gammaf
             
  real(double), allocatable :: swk(:,:,:),cden(:,:), fx(:,:), fy(:,:),phi(:,:),an(:,:), pt(:,:)
  
  real(double), allocatable :: dpef(:),eyef(:),bapef(:)

  integer :: inum,pmax,MAX_LOOP,mx,my,mx1,my1,jtl,ncon,i,ie,iemax,jt,range
  integer, allocatable :: ipt(:)

  contains

  subroutine config
    implicit none

    !---(Basic Parameters)---!

    temp      = 5d0                          !温度 [K]
    tr        = 1.64d-9                       !再結合時間 [s]
    t_p       = 0.846d-12                      !運動量緩和時間 [s] 
    param_bap = 4d20                        !BAP効果のパラメータ
    wd        = 2.2d-9                        !井戸の幅 [m]
    pw        = 60d-3                         !励起光強度（W）
    lay_r     = 1d-3                          !レーザーの半径（m）
    hz        = 800d6                         !周波数（Hz）
    n1        = 3.55d0                        !Al(0.3)Ga(0.7)Asの光子エネルギー1.7eVにおける屈折率
    n2        = 3.72d0                        !GaAsの光子エネルギー1.7eVにおける屈折率
    del       = 0.0132d0                      !Spin-orbit spilitting energy [eV]
    range     = 700                             !スピン緩和機構の寄与の平均をとる際の時間の範囲

    !---(Modeling Parameters)---!

    dx        = 1e-6                           !刻み幅 [m]
    dy        = dx                             !dx and dy should be the same for the poisson subruoutine
    mx        = 50                             !dxの数
    my        = 50
    dt        = 2e-15                          !刻み幅 [s]
    jtl       = 25000                           !dtの数
    de        = 0.001d0                        !エネルギーの刻み幅[eV]
    iemax     = 100                            !整数化したエネルギー
    xmax      = dx*(real(mx))                  !サンプル端のX座標
    ymax      = dy*(real(my))                  !サンプル端のy座標
    mx1       = mx+1                           !ｘ軸方向の格子点数
    my1       = my+1                           !ｚ軸方向の格子点数
    MAX_LOOP  = 1000000                       !ルーブの上限
    pmax      = 1000000                        !代表キャリア数の最大値
    ncon      = 30                             !１格子あたりの粒子数（計算負荷軽減のためのキャリアの代表点の数）

    !---(Material Parameters)---!

    am(1)    = 0.067d0*am0                     !GaAsのガンマ谷での有効質量[kg]
    am(2)    = 0.350d0*am0                     !GaAsのL谷での有効質量
    amh      = 0.51d0*am0                      !ホールの有効質量
    ameh(1)  = am(1) * amh / (am(1) + amh)     !電子とホールの換算質量
    ameh(2)  = am(2) * amh / (am(2) + amh)
    eps      = 12.90d0*ep0                     !半導体の誘電率(2.88)
    epf      = 10.92d0*ep0                     !光学的誘電率(2.83)
    ep       = 1d0/(1d0/epf-1d0/eps)           !P.48(2.89)
    lattice  = 5.65d-10                        !GaAsの格子定数

    !---(Phonon Parameters)---!
    
    rou     = 5360d0                           !P.38(2.53)
    sv      = 5240d0                           !結晶中の音速P.38(2.53)
    cl      = rou*sv*sv                        !決勝の弾性定数P.38
    z2      = 4d0                              !L谷の数
    da      = 7d0*q                            !変異ポテンシャルP.39(2.54)
    dij     = 1.0d11*q                         !バンド間散乱の変位ポテンシャルP.45(2.78)
    hwo     = 0.03536d0                        !ブリルアルゾーン中心付近での光学フォノンのエネルギー(eV)
    hwij    = 0.03d0                           !ブリルアルゾーン境界付近での光学フォノンのエネルギー(ev)
    hwe     = hwij                             !Lバンド内部での散乱に寄与する光学フォノンのエネルギー(eV)
    wo      = hwo*q/h                          !対応する角振動数
    wij     = hwij*q/h
    we      = hwe*q/h
    no      = 1.0d0/(exp(hwo*q / kb*temp)-1.0d0) !フォノンの数（ボース分布）
    nij     = 1.0d0/(exp(hwij*q / kb*temp)-1.0d0)
    ne      = 1.0d0/(exp(hwe*q / kb*temp)-1.0d0)

    !---(Coefficients)---!
    
    smh(1)  = sqrt(2d0 * am(1)) / h
    smh(2)  = sqrt(2d0 * am(2)) / h
    hhm(1)  = h**2 / (am(1) * 2d0)
    hhm(2)  = h**2 / (am(2) * 2d0)
    hm(1)   = h / am(1)
    hm(2)   = h / am(2)
    
    poe    = q**2 * wo / (8d0 * pi * ep) * (no+1d0)
    poa    = poe * no / (1d0+no)
    aco    = 2d0*pi * da**2 * kb*temp / (h * cl)       !P.61(2.134)
    ope    = pi * dij**2 / (rou * wij) * (nij + 1d0)   !P.46(2.80)　係数 emission
    opa    = ope * nij / (nij + 1d0)                   !absorption
    eqe    = pi * dij**2 / (we * rou) * (ne + 1d0)

    return
  end subroutine config

end module config_mod
