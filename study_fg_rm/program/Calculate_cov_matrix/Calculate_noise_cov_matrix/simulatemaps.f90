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
  integer, parameter :: nalm=3          !1=T, 2=E, 3=B
  integer, parameter :: ncl=4           !1=TT,2=EE,3=BB,4=TE
  integer, parameter :: Nside=4         !Healpix resolution
  integer, parameter :: npix=12*Nside**2!number of pixel
  integer, parameter :: nlmax0=2*Nside  !default lmax for the map
  !testing 2.5d0 times pixcel size
  !real(dl), parameter :: fwhm_arcmin = sqrt(fourpi/npix)*2.5d0*(rad_to_arcmin)
  real(dl), parameter :: fwhm_arcmin = 2220.0d0 !for Nside4
  !real(dl), parameter :: fwhm_arcmin = 1110.0d0 !for Nside8

  logical smoothing
  namelist /params/ smoothing

  !QQ,QU,UU_with_m_summed
  real(dl), dimension(0:nlmax0,0:npix-1,0:npix-1) :: &
       QQ_signal_m_summed_E,QU_signal_m_summed_E,&
       UU_signal_m_summed_E,UQ_signal_m_summed_E,&
       QQ_signal_m_summed_B,QU_signal_m_summed_B,&
       UU_signal_m_summed_B,UQ_signal_m_summed_B

contains

  subroutine Calculate_noise_cov_matrix(covmat_noise,sig_noise)
    use SpinHarmonics
    use alm_tools
    real(dl), intent(in) :: sig_noise
    real(dl), dimension(0:2*npix-1,0:2*npix-1), intent(out) :: covmat_noise

    integer, parameter :: nlmax=nlmax0 !maximum multipole to calculate covmat
    real(dl), dimension(0:npix-1,0:npix-1) :: noiseQQ,noiseQU,noiseUQ,noiseUU
    real(dl), dimension(0:nlmax0) :: NlEE,NlBB
    !pixel & beam window functions
    real(dl), dimension(0:nlmax,1:nalm) :: beam,w_ell
    !indices for directions (i,j), ell
    integer i,j,ell

    real(dl) noise_per_pixel    !testing

    !setting pixel & beam window functions
    call pixel_window(w_ell,Nside)  !get w_ell from Healpix
    call gaussbeam(fwhm_arcmin,nlmax,beam) !get beam from Healpix
    call NoisePowerSpectrum(sig_noise,nlmax,NlEE)
    call NoisePowerSpectrum(sig_noise,nlmax,NlBB)

    !$omp parallel do private(ell)
    do ell=2,nlmax
       !noiseはmap baseでシミュレーションしているので
       !beam windowだけ掛ける
       NlEE(ell) = NlEE(ell)*beam(ell,2)**2
       NlBB(ell) = NlBB(ell)*beam(ell,3)**2
    end do


    !初期化
    noiseQQ(:,:)=0._dl
    noiseQU(:,:)=0._dl
    noiseUQ(:,:)=0._dl
    noiseUU(:,:)=0._dl
       
    !$omp parallel do default(shared) 
    do i=0,npix-1
       do j=0,npix-1
          do ell=2,nlmax
             !noise covariance
             noiseQQ(i,j)=noiseQQ(i,j)&
                  +NlEE(ell)*QQ_signal_m_summed_E(ell,i,j)&
                  +NlBB(ell)*QQ_signal_m_summed_B(ell,i,j)

             noiseQU(i,j)=noiseQU(i,j)&
                  +NlEE(ell)*QU_signal_m_summed_E(ell,i,j)&
                  +NlBB(ell)*QU_signal_m_summed_B(ell,i,j)

             noiseUQ(i,j)=noiseUQ(i,j)&
                  +NlEE(ell)*UQ_signal_m_summed_E(ell,i,j)&
                  +NlBB(ell)*UQ_signal_m_summed_B(ell,i,j)

             noiseUU(i,j)=noiseUU(i,j)&
                  +NlEE(ell)*UU_signal_m_summed_E(ell,i,j)&
                  +NlBB(ell)*UU_signal_m_summed_B(ell,i,j)
          end do !ell loop
       end do !j loop
    end do!i loope

    !$omp parallel do default(shared) 
    do i=0,npix-1
       do j=0,npix-1
          covmat_noise(i,j)           =noiseQQ(i,j)
          covmat_noise(i,j+npix)      =noiseQU(i,j)
          covmat_noise(i+npix,j)      =noiseUQ(i,j)
          covmat_noise(i+npix,j+npix) =noiseUU(i,j)
       end do
    end do

    if(smoothing .eq. .false.) then
       !testing (diagnal noise covariance)
       covmat_noise(:,:) = 0.0d0
       noise_per_pixel = NoisePerPix(sig_noise,npix)
       do i=0,2*npix-1
          covmat_noise(i,i) = noise_per_pixel**2
       end do
    end if

    !testing output
    write(*,*)"diag element after smoothing and noise_per_pixel**2 comparison"
    do i=0,1 !2*npix-1
       noise_per_pixel = NoisePerPix(sig_noise,npix)
       write(*,'(3E15.5)') covmat_noise(i,i),noise_per_pixel**2,&
                           noise_per_pixel/sqrt(covmat_noise(i,i))
    end do

  end subroutine Calculate_noise_cov_matrix

  subroutine NoisePowerSpectrum(sig_noise,nlmax,Nell)
    real(dl), intent(in) :: sig_noise
    integer, intent(in) :: nlmax
    real(dl), intent(out) :: Nell(0:nlmax)
    integer l

    Nell(:) = 0._dl !初期化
    do l=2,nlmax
       Nell(l) = ((pi/10800._dl)*sig_noise)**2
    end do
      
  end subroutine NoisePowerSpectrum


  function NoisePerPix(sig_noise,npixel) result(noise_per_pix)
    integer npixel
    real(dl) sig_noise
    real(dl) noise_per_pix
    real(dl) solid_angle
         
    solid_angle = (fourpi/npixel)
    noise_per_pix = ((pi/10800._dl)*sig_noise)/sqrt(solid_angle)
         
  end function NoisePerPix

  subroutine Output_noise_cov_matrix(filename,covmat_noise)
    real(dl), dimension(0:2*npix-1,0:2*npix-1), intent(in) :: covmat_noise
    character(LEN=1024), intent(in) :: filename
    open(67,file=filename,form='binary',access='sequential')
    write(67) covmat_noise
    close(67)

  end subroutine Output_noise_cov_matrix

  subroutine Set_noise_cov_matrix(filename,covmat_noise)
    real(dl), dimension(0:2*npix-1,0:2*npix-1), intent(out) :: covmat_noise
    character(LEN=1024), intent(in) :: filename
    open(67,file=filename,form='binary',access='sequential')
    read(67) covmat_noise
    close(67)

  end subroutine Set_noise_cov_matrix

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
    !-----------------------------------------------------------
    close(67)
  end subroutine SetXX_with_m_summed

  function mytrim(file)
    character(LEN=*), intent(in) :: file
    character(LEN=:), allocatable:: mytrim
    
    mytrim = trim(adjustl(file))
  end function mytrim
end module simulatemaps


