module simulatemaps
  !module SpinHarmonicsは http://cosmologist.info/polar/ から。
  !使い方は
  !     lambda = SH_lambdaslm(s,l,m,x) 
  !where sYlm = lambda*e^{i*m*phi}, s is the spin, and x is cos(theta). 
  !Sign conventions and documentation in astro-ph/0106536.
  !module pix_tools は Healpix 3.20 から。Makefile参照
  use SpinHarmonics
  use pix_tools
  implicit none

  !定数を定義
  integer, parameter :: dl = KIND(1d0) ! selected_real_kind(15,50)
  real(dl), parameter :: pi=3.14159265358979323846264338328_dl, &
           twopi=2*pi, fourpi=4*pi
  real(dl), parameter :: rad_to_arcmin = 360.0d0*60.0d0/twopi ![arcmin/rad]
  complex(dl), parameter :: img=(0._dl,1._dl) !純虚数

  !model parameter settings
  !ここでファイルのフォーマットとか、Healpixのパラメタなどを
  !設定する
  integer, parameter :: nalm=3           !1=T, 2=E, 3=B
  integer, parameter :: ncl=4            !1=TT,2=EE,3=BB,4=TE
  integer, parameter :: Nside=4          !Healpix resolution
  integer, parameter :: npix=12*Nside**2 !number of pixel
  integer, parameter :: nlmax0=2*Nside   !default lmax for the map
!!$  integer, parameter :: nlmax0=3*Nside-1 !testing (nyquist frequency)


  !smoothing scale
!!$  real(dl), parameter :: fwhm_arcmin = sqrt(fourpi/npix)*2.5d0*(rad_to_arcmin)  !2.5d0 times pixcel size
  real(dl), parameter :: fwhm_arcmin = 2220.0d0 !Nside4

  !読み込んだClを入れておく配列
  integer, parameter :: lmin=2          !Cl fileのlmin
  integer, parameter :: lmax_scal=1750  !KK2011 Cl fileのlmax
  integer, parameter :: lmax_tens=1500  !Cl fileのlmax
  real(dl), dimension(0:lmax_scal,1:ncl) :: ScalCls,TensCls,TotCls

  !lensed Clを用いるかどうかとKK2011時代のClを使うかどうか。
  logical, parameter :: read_lensed_Cls = .true.
  logical, parameter :: KK2011= .false.

  real(dl),allocatable,dimension(:,:,:) :: &
       QQ_signal_m_summed_E,QU_signal_m_summed_E,&
       UU_signal_m_summed_E,UQ_signal_m_summed_E,&
       QQ_signal_m_summed_B,QU_signal_m_summed_B,&
       UU_signal_m_summed_B,UQ_signal_m_summed_B

contains

  subroutine ReadCls

    !cambで計算したClファイルを読み込む
    !単位は (muK)^2
    real(dl), dimension(lmax_scal) :: ScalClTT,ScalClEE,ScalClTE,&
         ScalClphiphi,ScalClTphi,ScalClBB
    real(dl), dimension(lmax_tens) :: TensClTT,TensClEE,TensClBB,&
         TensClTE
    character(LEN=:), allocatable :: filename
    integer i,dummy
    real(dl) fact

    !初期化
    ScalCls(:,:) = 0._dl
    TensCls(:,:) = 0._dl
    TotCls(:,:) = 0._dl
    
    if(read_lensed_Cls) then
       !lensed scalar modeの読み込み
       if(KK2011) then
          filename = './Cls/lblike_lensedCls.dat'
       else
          filename = './Cls/planck15lens_lensedCls.dat'
       end if

       open(11,file=filename)
       do i = 2,lmax_scal
          read(11,*) dummy,ScalClTT(i),ScalClEE(i),ScalClBB(i),ScalClTE(i)
          fact=dble(i*(i+1))/twopi

          ScalCls(i,1) = ScalClTT(i)/fact
          ScalCls(i,2) = ScalClEE(i)/fact
          ScalCls(i,3) = ScalClBB(i)/fact
          ScalCls(i,4) = ScalClTE(i)/fact
       enddo
       close(11)
    else
       !unlensed scalar modeの読み込み
       if(KK2011) then
          filename = './Cls/lblike_scalCls.dat'
       else
          filename = './Cls/PlanckTTlowP_scalCls.dat'
       end if

       open(11,file=filename)
       do i = 2,lmax_scal
          if(KK2011) then
             read(11,*) dummy,ScalClTT(i),ScalClEE(i),ScalClTE(i)
          else
             read(11,*) dummy,ScalClTT(i),ScalClEE(i),ScalClTE(i)&
                  ,ScalClphiphi(i),ScalClTphi(i)
          end if

          !input sを掛けるのを止めた
          !!testing
          stop'do not come here'
          fact=dble(i*(i+1))/twopi
          ScalCls(i,1) = ScalClTT(i)/fact
          ScalCls(i,2) = ScalClEE(i)/fact
          ScalCls(i,3) = 0._dl
          ScalCls(i,4) = ScalClTE(i)/fact
       enddo
       close(11)
    end if


    !unlensed tensor modeの読み込み
    if(KK2011) then
       filename = './Cls/lblike_tensCls.dat'
    else
       filename = './Cls/planck15lens_tensCls.dat'
    end if

    open(11,file=filename)
    do i = 2,lmax_tens
       read(11,*) dummy,TensClTT(i),TensClEE(i),TensClBB(i)&
            ,TensClTE(i)

       !input rを掛けるのを止めた
       fact=dble(i*(i+1))/twopi

       TensCls(i,1) = TensClTT(i)/fact
       TensCls(i,2) = TensClEE(i)/fact
       TensCls(i,3) = TensClBB(i)/fact
       TensCls(i,4) = TensClTE(i)/fact
    enddo
    close(11)

  end subroutine ReadCls


  subroutine Calculate_signal_cov_matrix(covmat_signal,r_in,s_in)
    use SpinHarmonics
    use alm_tools
    real(dl), intent(in) :: r_in,s_in
    real(dl), dimension(0:2*npix-1,0:2*npix-1), intent(out) :: covmat_signal

    integer, parameter :: nlmax=nlmax0 !maximum multipole to calculate covmat
    real(dl), dimension(0:npix-1,0:npix-1) :: covQQ,covQU,covUQ,covUU
    !angular power spectra
    real(dl), dimension(0:nlmax) :: ClEE,ClBB
    !pixel & beam window functions
    real(dl), dimension(0:nlmax,1:nalm) :: beam,w_ell
    !indices for directions (i,j), ell
    integer i,j,ell

    !setting pixel & beam window functions
    call pixel_window(w_ell,Nside)  !get w_ell from Healpix
    call gaussbeam(fwhm_arcmin,nlmax,beam) !get beam from Healpix

    !$omp parallel do private(ell)
    do ell=2,nlmax
       !defalt
       ClEE(ell) = (s_in*ScalCls(ell,2)+r_in*TensCls(ell,2))&
            *w_ell(ell,2)**2*beam(ell,2)**2

       ClBB(ell) = (s_in*ScalCls(ell,3)+r_in*TensCls(ell,3))&
            *w_ell(ell,3)**2*beam(ell,3)**2
    end do

!!$    !testing output
!!$    open(60,file='signal_cl_scalar_lens.txt')
!!$    do ell=2,nlmax
!!$       write(60,*) ell,ClEE(ell),ClBB(ell)
!!$    end do
!!$    stop

    !初期化
    covQQ(:,:)=0._dl
    covQU(:,:)=0._dl
    covUQ(:,:)=0._dl
    covUU(:,:)=0._dl
       
    !$omp parallel do default(shared) 
    do i=0,npix-1
       do j=0,npix-1
          do ell=2,nlmax
             !noise covariance
             covQQ(i,j)=covQQ(i,j)&
                  +ClEE(ell)*QQ_signal_m_summed_E(ell,i,j)&
                  +ClBB(ell)*QQ_signal_m_summed_B(ell,i,j)

             covQU(i,j)=covQU(i,j)&
                  +ClEE(ell)*QU_signal_m_summed_E(ell,i,j)&
                  +ClBB(ell)*QU_signal_m_summed_B(ell,i,j)

             covUQ(i,j)=covUQ(i,j)&
                  +ClEE(ell)*UQ_signal_m_summed_E(ell,i,j)&
                  +ClBB(ell)*UQ_signal_m_summed_B(ell,i,j)

             covUU(i,j)=covUU(i,j)&
                  +ClEE(ell)*UU_signal_m_summed_E(ell,i,j)&
                  +ClBB(ell)*UU_signal_m_summed_B(ell,i,j)
          end do !ell loop
       end do !j loop
    end do!i loope

    !$omp parallel do default(shared) 
    do i=0,npix-1
       do j=0,npix-1
          covmat_signal(i,j)           =covQQ(i,j)
          covmat_signal(i,j+npix)      =covQU(i,j)
          covmat_signal(i+npix,j)      =covUQ(i,j)
          covmat_signal(i+npix,j+npix) =covUU(i,j)
       end do
    end do

  end subroutine Calculate_signal_cov_matrix


  subroutine Output_signal_cov_matrix(filename,covmat_signal)
    real(dl), dimension(0:2*npix-1,0:2*npix-1), intent(in) :: covmat_signal
    character(LEN=1024), intent(in) :: filename

    open(67,file=filename,form='binary',access='sequential')
    write(67) covmat_signal
    close(67)

  end subroutine Output_signal_cov_matrix

  subroutine Set_signal_cov_matrix(filename,covmat_signal)
    real(dl), dimension(0:2*npix-1,0:2*npix-1), intent(out) :: covmat_signal
    character(LEN=1024), intent(in) :: filename

    open(67,file=filename,form='binary',access='sequential')
    read(67) covmat_signal
    close(67)

  end subroutine Set_signal_cov_matrix

  subroutine SetXX_with_m_summed(filename)
    character(LEN=1024),intent(in) :: filename

    !計算したxx_signal_m_summedをファイルから取り戻す-----------------
    open(67,file=filename,form='binary',access='sequential')
    read(67) QQ_signal_m_summed_E,&
             QU_signal_m_summed_E,&
             UQ_signal_m_summed_E,&
             UU_signal_m_summed_E,&
             QQ_signal_m_summed_B,&
             QU_signal_m_summed_B,&
             UQ_signal_m_summed_B,&
             UU_signal_m_summed_B
    close(67)

  end subroutine SetXX_with_m_summed

  subroutine Allocate_arrays
    !QQ,QU,UU_with_m_summed
    allocate(QQ_signal_m_summed_E(0:nlmax0,0:npix-1,0:npix-1))
    allocate(QU_signal_m_summed_E(0:nlmax0,0:npix-1,0:npix-1))
    allocate(UU_signal_m_summed_E(0:nlmax0,0:npix-1,0:npix-1))
    allocate(UQ_signal_m_summed_E(0:nlmax0,0:npix-1,0:npix-1))
    allocate(QQ_signal_m_summed_B(0:nlmax0,0:npix-1,0:npix-1))
    allocate(QU_signal_m_summed_B(0:nlmax0,0:npix-1,0:npix-1))
    allocate(UU_signal_m_summed_B(0:nlmax0,0:npix-1,0:npix-1))
    allocate(UQ_signal_m_summed_B(0:nlmax0,0:npix-1,0:npix-1))
  end subroutine Allocate_arrays

  subroutine deAllocate_arrays
    !QQ,QU,UU_with_m_summed
    deallocate(QQ_signal_m_summed_E)
    deallocate(QU_signal_m_summed_E)
    deallocate(UU_signal_m_summed_E)
    deallocate(UQ_signal_m_summed_E)
    deallocate(QQ_signal_m_summed_B)
    deallocate(QU_signal_m_summed_B)
    deallocate(UU_signal_m_summed_B)
    deallocate(UQ_signal_m_summed_B)

  end subroutine deAllocate_arrays


end module simulatemaps


