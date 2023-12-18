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

  !model parameter settings
  !ここでファイルのフォーマットとか、Healpixのパラメタなどを設定する
  integer, parameter :: ncl=4             !1=TT,2=EE,3=BB,4=TE
  integer, parameter :: Nside=4           !Healpix resolution
  integer, parameter :: npix=12*Nside**2  !number of pixel
  integer, parameter :: nlmax0=2*Nside   !default lmax 

  !QQ,QU,UU_with_m_summed
  real(dl), allocatable, dimension(:,:,:) :: &
       QQ_signal_m_summed_E,QU_signal_m_summed_E,&
       UU_signal_m_summed_E,UQ_signal_m_summed_E,&
       QQ_signal_m_summed_B,QU_signal_m_summed_B,&
       UU_signal_m_summed_B,UQ_signal_m_summed_B

contains

  subroutine CalculateXX_with_m_summed
    use SpinHarmonics
    use alm_tools
    complex(dl), parameter :: img=(0._dl,1._dl)!純虚数

    !functions defined in KK2011
    complex(dl) Wlm_i,Xlm_i,Wlm_j,Xlm_j
    complex(dl) Yp2lm_i,Ym2lm_i,Yp2lm_j,Ym2lm_j
    complex(dl) sumQQ_E,sumQQ_B,sumQU_E,sumQU_B,&
         sumUQ_E,sumUQ_B,sumUU_E,sumUU_B

      
    !arrays for spherical coordinates (phi,cos(theta))
    real(dl), dimension(0:npix-1) :: phi_arr,mu_arr

    !colatitude in radians and its cosine for each direction n_i and n_j
    real(dl) theta,mu_i,mu_j
    !longitude in radians
    real(dl) phi,phi_i,phi_j               
    !indices for directions (i,j), ell and m
    integer i,j,ell,m

    !初期化
    QQ_signal_m_summed_E(:,:,:)=0._dl
    QU_signal_m_summed_E(:,:,:)=0._dl
    UQ_signal_m_summed_E(:,:,:)=0._dl
    UU_signal_m_summed_E(:,:,:)=0._dl

    QQ_signal_m_summed_B(:,:,:)=0._dl
    QU_signal_m_summed_B(:,:,:)=0._dl
    UQ_signal_m_summed_B(:,:,:)=0._dl
    UU_signal_m_summed_B(:,:,:)=0._dl


    write(*,*)"setting theta & phi"
    !$omp parallel do private(i,theta,phi)
    do i=0,npix-1
       call pix2ang_nest(Nside,i,theta,phi)
       phi_arr(i)=phi
       mu_arr(i)=cos(theta)
    end do

    do i=0,npix-1
       phi_i=phi_arr(i)
       mu_i = mu_arr(i)

       !$omp parallel do default(shared) &
       !$omp & private(j,mu_j,phi_j,Yp2lm_i,Yp2lm_j,Ym2lm_i,Ym2lm_j,Wlm_i,Wlm_j,Xlm_i,Xlm_j) &
       !$omp & reduction(+:sumQQ_E,sumQU_E,sumUQ_E,sumUU_E,sumQQ_B,sumQU_B,sumUQ_B,sumUU_B)  
       do j=0,npix-1
          phi_j=phi_arr(j)
          mu_j = mu_arr(j)

          do ell=2,nlmax0
             !初期化
             sumQQ_E=0._dl*img
             sumQU_E=0._dl*img
             sumUQ_E=0._dl*img
             sumUU_E=0._dl*img
             sumQQ_B=0._dl*img
             sumQU_B=0._dl*img
             sumUQ_B=0._dl*img
             sumUU_B=0._dl*img

             do m=-ell,ell
                Yp2lm_i = SH_lambdaslm( 2,ell,m,mu_i)*exp(img*m*phi_i)
                Ym2lm_i = SH_lambdaslm(-2,ell,m,mu_i)*exp(img*m*phi_i)
                Yp2lm_j = SH_lambdaslm( 2,ell,m,mu_j)*exp(img*m*phi_j)
                Ym2lm_j = SH_lambdaslm(-2,ell,m,mu_j)*exp(img*m*phi_j)

                Wlm_i=(-1._dl)*(Yp2lm_i+Ym2lm_i)/2
                Xlm_i=(-img)  *(Yp2lm_i-Ym2lm_i)/2
                Wlm_j=(-1._dl)*(Yp2lm_j+Ym2lm_j)/2
                Xlm_j=(-img)  *(Yp2lm_j-Ym2lm_j)/2

                sumQQ_E=sumQQ_E +  Wlm_i*conjg(Wlm_j)
                sumQQ_B=sumQQ_B +  Xlm_i*conjg(Xlm_j)

                sumQU_E=sumQU_E +(-Wlm_i)*conjg(Xlm_j)
                sumQU_B=sumQU_B +  Xlm_i*conjg(Wlm_j)

                sumUQ_E=sumUQ_E +(-Xlm_i)*conjg(Wlm_j)
                sumUQ_B=sumUQ_B +  Wlm_i*conjg(Xlm_j)

                sumUU_E=sumUU_E +  Xlm_i*conjg(Xlm_j)
                sumUU_B=sumUU_B +  Wlm_i*conjg(Wlm_j)
             end do !m loop
             
             QQ_signal_m_summed_E(ell,i,j)=real(sumQQ_E)
             QU_signal_m_summed_E(ell,i,j)=real(sumQU_E)
             UQ_signal_m_summed_E(ell,i,j)=real(sumUQ_E)
             UU_signal_m_summed_E(ell,i,j)=real(sumUU_E)

             QQ_signal_m_summed_B(ell,i,j)=real(sumQQ_B)
             QU_signal_m_summed_B(ell,i,j)=real(sumQU_B)
             UQ_signal_m_summed_B(ell,i,j)=real(sumUQ_B)
             UU_signal_m_summed_B(ell,i,j)=real(sumUU_B)

          end do !ell loop
       end do !j loop
    end do!i loop

  end subroutine CalculateXX_with_m_summed

  subroutine Output_XX_with_m_summed(filename)
    character(LEN=1024), intent(in) :: filename
    !計算したxx_signal_m_summedをファイルに保存する ---------
    open(67,file=filename,form='binary',access='sequential')
    write(67) QQ_signal_m_summed_E,&
              QU_signal_m_summed_E,&
              UQ_signal_m_summed_E,&
              UU_signal_m_summed_E,&
              QQ_signal_m_summed_B,&
              QU_signal_m_summed_B,&
              UQ_signal_m_summed_B,&
              UU_signal_m_summed_B
    close(67)

  end subroutine Output_XX_with_m_summed

    
  subroutine SetXX_with_m_summed(filename)
    character(LEN=1024),intent(in):: filename
    !計算したxx_signal_m_summedをファイルから取り戻す--------
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


