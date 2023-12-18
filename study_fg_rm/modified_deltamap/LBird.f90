module LBird
  use Utils
  implicit none

  !global variables
  !model parameter settings
  integer Nside,npix  !Healpix resolution
!!$  integer nvar_fore   !new number of foreground parameters
  integer num_source  !=5 (mean(beta_s) + delta(beta_s) + mean(beta_d, T_d)+delta(beta_d)+ delta(T_d)
  logical noiseless   !if maps are generated withoug noise
  character(LEN=128) algorithm
  real(dl) sig_noise  !noise per sqrt(Omega_pix)=1 arcmin pixel (diag noise)
  real(dl) Td_Kelvin  !dust temperature
  namelist /Params/ Nside,noiseless,sig_noise,algorithm,Td_Kelvin!,nvar_fore

  !LiteBIRDの周波数 [GHz]
  integer numchan
  real(dl), allocatable,dimension(:) :: nu_array !new CMB+Foreground frequencies
  namelist /numCMB/ numchan !number of CMB channels
  namelist /bands/ nu_array 

  !signal & noise covariance matrices
  real(dl), allocatable, dimension(:,:) :: covmat_signal_scalar,covmat_signal_tensor
  real(dl), allocatable, dimension(:,:,:) :: covmat_noise_nu !new
  
  !mapをファイルから読むための配列
  integer, parameter :: nmap=3  !1=T, 2=Q, 3=U
  !maps, dimension(0:npix-1,1:nmap,1:numchan) 
  real(dl), allocatable, dimension(:,:,:) :: map_chan !new
  !diag noiseのmap
  real(dl), allocatable, dimension(:,:,:) :: diagnoise_map_tot_all !new

  !mask map (偏光のmask(1種類)のみ),dimension(0:npix-1,1:nmasks)
  integer, parameter :: nmasks=1 
  real(dl), allocatable, dimension(:,:) :: mask  

  !maskの情報を用いて配列を並び替えるためのvector
  integer,allocatable,dimension(:) :: reorder_vector_tot_all
  integer,allocatable,dimension(:) :: reorder_vector_for_Dmatrix
  
  !maskされたpixelの数とされてない数 npix_unmasked = npix-npix_masked
  integer npix_masked,npix_unmasked
  integer npix_tot_all !new (numchan*2*npix) main.f90でset
  integer npix_unmasked_tot_all,npix_masked_tot_all !new (=numchan*2*npix_unmaskedなど)subroutine WriteReOrderVectorでset
  integer npix_unmasked_source_all,npix_masked_source_all !new (=num_source*2*npix_unmaskedなど)subroutine WriteReOrderVectorでset
  real(dl) Td_GHz

contains

  subroutine Allocate_global_arrays

    allocate(reorder_vector_tot_all(0:npix_tot_all-1)) !new
    allocate(reorder_vector_for_Dmatrix(0:2*npix*num_source-1)) !new
    allocate(nu_array(1:numchan)) !new

    allocate(covmat_signal_scalar(0:2*npix-1,0:2*npix-1))
    allocate(covmat_signal_tensor(0:2*npix-1,0:2*npix-1))

    allocate(covmat_noise_nu(0:2*npix-1,0:2*npix-1,1:numchan))!new

    allocate(map_chan(0:npix-1,1:nmap,1:numchan)) !new
    allocate(diagnoise_map_tot_all(0:npix-1,1:nmap,1:numchan)) !new

    allocate(mask(0:npix-1,1:nmasks))

  end subroutine Allocate_global_arrays

  subroutine deAllocate_global_arrays

    deallocate(covmat_signal_scalar)
    deallocate(covmat_signal_tensor)

    deallocate(covmat_noise_nu) !new

    deallocate(mask)

    deallocate(reorder_vector_tot_all) !new
    deallocate(reorder_vector_for_Dmatrix) !new
    deallocate(nu_array) !new
    deallocate(diagnoise_map_tot_all) !new
    deallocate(map_chan) !new

  end subroutine deAllocate_global_arrays
  

  subroutine ReadMaskFits(filename)
    use fitstools
    logical anynull
    real(dl) nullval
    character(LEN=:), allocatable :: filename
    call read_bintab(filename,mask,npix,nmasks,nullval,anynull)
  end subroutine ReadMaskFits

  subroutine CountMaskedPixels
    integer i,j
    !maskされていないpixelの数を数える
    j=0
    do i=0,npix-1
       if(mask(i,1) == 1) j=j+1
    end do

    npix_unmasked = j                  !global変数へ
    npix_masked = npix - npix_unmasked !global変数へ

    npix_unmasked_tot_all = numchan*2*npix_unmasked !global変数へ(factor 2 includes Q and U)
    npix_masked_tot_all = numchan*2*npix_masked     !global変数へ

    npix_unmasked_source_all = num_source*2*npix_unmasked !global変数へ(factor 2 includes Q and U)
    npix_masked_source_all = num_source*2*npix_masked     !global変数へ
    
  end subroutine CountMaskedPixels

  subroutine WriteReOrderVector
    use Utils
    integer i,j,k
    integer,allocatable,dimension(:) :: reorder_vector_back_tot_all !new

    allocate(reorder_vector_back_tot_all(0:npix_tot_all-1)) !new

    !初期化
    reorder_vector_tot_all(:) = -1._dl       !global
    reorder_vector_back_tot_all(:) = -1._dl  !local

    !まず、maskされていないpixelを拾って並べる
    j=0
    do i=0,npix-1
          if(mask(i,1) == 1) then
             !new---
             do k = 1,numchan    !ふたつ目のmapの分も考慮する
                reorder_vector_tot_all(j+    2*(k-1)*npix_unmasked) = i+    2*(k-1)*npix !Q
                reorder_vector_tot_all(j+(2*(k-1)+1)*npix_unmasked) = i+(2*(k-1)+1)*npix !U
         
                !逆変換
                reorder_vector_back_tot_all(i+    2*(k-1)*npix)=j+    2*(k-1)*npix_unmasked !Q
                reorder_vector_back_tot_all(i+(2*(k-1)+1)*npix)=j+(2*(k-1)+1)*npix_unmasked !U
             end do
             !-----
             
             j=j+1
          end if
       end do

    !次に、maskされているpixelを拾って並べる
    j=0
    do i=0,npix-1
       if(mask(i,1) == 0) then
             !new----
             do k = 1,numchan    !ふたつ目のmapの分も考慮する
                reorder_vector_tot_all(npix_unmasked_tot_all+j+2*(k-1)*npix_masked)=i+2*(k-1)*npix
                reorder_vector_tot_all(npix_unmasked_tot_all+j+(2*(k-1)+1)*npix_masked)=i+(2*(k-1)+1)*npix

                !逆変換
                reorder_vector_back_tot_all(i+2*(k-1)*npix) = npix_unmasked_tot_all+j+2*(k-1)*npix_masked
                reorder_vector_back_tot_all(i+(2*(k-1)+1)*npix) = npix_unmasked_tot_all+j+(2*(k-1)+1)*npix_masked
             end do
             !------
             j=j+1
          end if
       end do

!!$      !testing output
!!$      do i=0,npix_tot_all-1
!!$         write(80,*)"reorder_vector(",i,")=",reorder_vector_tot_all(i),reorder_vector_back_tot_all(i)
!!$      end do
!!$      stop
      !--------------

       deallocate(reorder_vector_back_tot_all)

  end subroutine WriteReOrderVector


    subroutine WriteReOrderVector_for_Dmatrix
    use Utils
    integer i,j,k

    !初期化
    reorder_vector_for_Dmatrix(:) = -1._dl       !global

    !まず、maskされていないpixelを拾って並べる
    j=0
    do i=0,npix-1
          if(mask(i,1) == 1) then
             !new---
             do k = 1,num_source    !ふたつ目のmapの分も考慮する
                reorder_vector_for_Dmatrix(j+    2*(k-1)*npix_unmasked) = i+    2*(k-1)*npix !Q
                reorder_vector_for_Dmatrix(j+(2*(k-1)+1)*npix_unmasked) = i+(2*(k-1)+1)*npix !U
             end do
             !-----
             
             j=j+1
          end if
       end do

    !次に、maskされているpixelを拾って並べる
    j=0
    do i=0,npix-1
       if(mask(i,1) == 0) then
             !new----
             do k = 1,num_source    !ふたつ目のmapの分も考慮する
                reorder_vector_for_Dmatrix(npix_unmasked_source_all+j+2*(k-1)*npix_masked)=i+2*(k-1)*npix
                reorder_vector_for_Dmatrix(npix_unmasked_source_all+j+(2*(k-1)+1)*npix_masked)=i+(2*(k-1)+1)*npix
             end do
             !------
             j=j+1
          end if
       end do

     end subroutine WriteReOrderVector_For_Dmatrix


  subroutine CalculateCovMatrix(r_est,s_est,n1_est,covmat_tot,covmat_signal,covmat_noise)
    !alpha1,alpha2の二つを含むdelta map用
    use Utils
    implicit none

    real(dl), intent(in) :: r_est,s_est,n1_est
    real(dl), dimension(0:npix_tot_all-1,0:npix_tot_all-1), intent(out) :: covmat_tot,covmat_signal,covmat_noise

    !noise per pix for artificial noise diagnal matrix
    real(dl) noise_per_pixel

    !indices for directions (i,j), ell and m
    integer i,j
    integer index_i,index_j,index_ip,index_jp

    !初期化
    covmat_signal(:,:)=0._dl
    covmat_noise(:,:)=0._dl
    covmat_tot(:,:)=0._dl

    !CMB signal covariance matrix
    !fileから読み込んだ(r=1,s=1)のsignal covarianceに推定するr,sを掛けるだけ
    covmat_signal(0:2*npix-1,0:2*npix-1) = covmat_signal_scalar(0:2*npix-1,0:2*npix-1)*s_est &
                                         + covmat_signal_tensor(0:2*npix-1,0:2*npix-1)*r_est

    !signal covarianceをnumchan個のmapについて求める。
    do i=1,numchan
       index_i  = 2*(i-1)*npix
       index_ip = 2*i*npix - 1

       do j=1,numchan
          index_j  = 2*(j-1)*npix 
          index_jp = 2*j*npix - 1

          covmat_signal(index_i:index_ip,index_j:index_jp)=covmat_signal(0:2*npix-1,0:2*npix-1)

       end do
    end do

    !noiseは今回は対角ブロックのみ
    do i=1,numchan
       index_i  = 2*(i-1)*npix
       index_ip = 2*i*npix-1 
       covmat_noise(index_i:index_ip,index_i:index_ip)=covmat_noise_nu(0:2*npix-1,0:2*npix-1,i)
    end do
    covmat_noise=n1_est*covmat_noise

    if(noiseless) covmat_noise(:,:)=0._dl

    !singular matrixを避けるためのdiagnal noise matrix
    !noise covarianceに足しておく
    noise_per_pixel = NoisePerPix(sig_noise,npix)

    !$omp parallel do default(shared) 
    do i=0,npix_tot_all-1
       covmat_noise(i,i) = covmat_noise(i,i) + noise_per_pixel**2         
    end do

    covmat_tot = covmat_signal + covmat_noise

  end subroutine CalculateCovMatrix


  subroutine ReOrderCovMatrix(mat,mat_reordered,reorder_vector)
    use Utils
    real(dl), dimension(0:npix_tot_all-1,0:npix_tot_all-1), intent(in) :: mat !input
    real(dl), dimension(0:npix_tot_all-1,0:npix_tot_all-1), intent(out) :: mat_reordered !output
    integer, dimension(0:npix_tot_all-1), intent(in) :: reorder_vector
   
    integer i,j

    mat_reordered(:,:) = 0._dl

    !$omp parallel do default(shared) 
    do i=0,npix_tot_all-1
       do j=0,npix_tot_all-1
          mat_reordered(i,j) = mat(reorder_vector(i),reorder_vector(j))
       end do
    end do
   
  end subroutine ReOrderCovMatrix

    subroutine ReOrderDMatrix(mat,mat_reordered,reorder_vector)
    use Utils
    real(dl), dimension(0:npix_tot_all-1,0:2*npix*num_source-1), intent(in) :: mat !input
    real(dl), dimension(0:npix_tot_all-1,0:2*npix*num_source-1), intent(out) :: mat_reordered !output
    integer, dimension(0:npix_tot_all-1), intent(in) :: reorder_vector
   
    integer i,j

    mat_reordered(:,:) = 0._dl

    !$omp parallel do default(shared) 
    do i=0,npix_tot_all-1
       do j=0,2*npix*num_source-1
          mat_reordered(i,j) = mat(reorder_vector(i),reorder_vector_for_Dmatrix(j))
       end do
    end do
   
  end subroutine ReOrderDMatrix

  subroutine CalcDmatrix(beta_s,beta_d,gamma_d,nu,Dmatrix)
    !this subroutine sets tilde(D) matrix in Eq(59) 
    use Utils
    implicit none
    real(dl), intent(in) :: beta_s,beta_d,gamma_d
    real(dl), dimension(1:numchan), intent(in) :: nu
    real(dl), dimension(0:npix_tot_all-1,0:2*npix*num_source-1), intent(inout) :: Dmatrix

    real(dl), dimension(1:numchan) :: gnu
    real(dl), dimension(1:num_source,1:numchan) :: inmat
    real(dl), dimension(0:2*npix-1,0:2*npix-1) :: E
    integer i,j,k
    integer index_i,index_ip,index_j,index_jp

    real(dl), parameter :: nu_s = 60._dl
    real(dl), parameter :: nu_d = 280._dl
    
    !set identity matrix
    E(:,:) = 0._dl
    do i = 0, 2*npix-1
       E(i,i) = 1._dl
    end do
    
    do i = 1,numchan
       gnu(i) = GetThermoFactor(nu(i))
    end do

    if(algorithm .eq. "modified_delta_map") then
       do i = 1,numchan
!!$         inmat(1,i) = gnu(i)*nu(i)**beta_s
!!$         inmat(2,i) = gnu(i)*nu(i)**beta_s * log(nu(i))
!!$         inmat(3,i) = gnu(i)*nu(i)**beta_d / (exp(nu(i)/gamma_d)-1._dl)
!!$         inmat(4,i) = inmat(3,i) * (nu(i)/gamma_d) * exp(nu(i)/gamma_d)/(exp(nu(i)/gamma_d)-1._dl)
!!$         inmat(5,i) = inmat(3,i) * log(nu(i))

          inmat(1,i) = gnu(i)*(nu(i)/nu_s)**beta_s
          inmat(2,i) = gnu(i)*(nu(i)/nu_s)**beta_s * log(nu(i)/nu_s)
          inmat(3,i) = gnu(i)*(nu(i)/nu_d)**beta_d / (exp(nu(i)/gamma_d)-1._dl) * (exp(nu_d/gamma_d)-1._dl)
          inmat(4,i) = inmat(3,i) * ((nu(i)/gamma_d) * exp(nu(i)/gamma_d)/(exp(nu(i)/gamma_d)-1._dl)-(nu_d/gamma_d) * exp(nu_d/gamma_d)/(exp(nu_d/gamma_d)-1._dl))
          inmat(5,i) = inmat(3,i) * log(nu(i)/nu_d)
       end do
    end if

    if(algorithm .eq. "PowLawDust") then
       do i = 1,numchan
          inmat(1,i) = gnu(i)*(nu(i)/nu_d)**beta_d
       end do
    end if

    if(algorithm .eq. "PowLawSync") then
       do i = 1,numchan
          inmat(1,i) = gnu(i)*(nu(i)/nu_s)**beta_s
       end do
    end if

    if(algorithm .eq. "PowLawDustSync") then
       do i = 1,numchan
          inmat(2,i) = gnu(i)*(nu(i)/nu_s)**beta_s
          inmat(1,i) = gnu(i)*(nu(i)/nu_d)**beta_d
       end do
    end if

    if(algorithm .eq. "psm162_deltamap") then
       do i = 1,numchan
          inmat(1,i) = gnu(i)*(nu(i)/nu_s)**beta_s
          inmat(2,i) = gnu(i)*(nu(i)/nu_s)**beta_s* log(nu(i)/nu_s)
          inmat(3,i) = gnu(i)*(nu(i)/nu_d)**beta_d
          inmat(4,i) = gnu(i)*(nu(i)/nu_d)**beta_d* log(nu(i)/nu_d)
       end do
    end if

    
    
    do i=1,numchan
       index_i = 2*(i-1)*npix
       index_ip = 2*i*npix-1

       do j=1,num_source
          index_j = 2*(j-1)*npix
          index_jp = 2*j*npix-1

          Dmatrix(index_i:index_ip,index_j:index_jp) = inmat(j,i)*E(0:2*npix-1,0:2*npix-1)

       end do
    end do
    
  end subroutine CalcDmatrix
  
  subroutine CalcAlphas(beta_s,beta_d,gamma_d,nu,alphas)
    use Utils
    implicit none
    real(dl), intent(in) :: beta_s,beta_d,gamma_d
    real(dl), dimension(1:numchan), intent(in) :: nu
    real(dl), dimension(1:numchan), intent(inout) :: alphas
    real(dl), dimension(1:numchan) :: gnu
    real(dl), parameter :: nu_s = 100._dl
    real(dl), parameter :: nu_d = 400._dl
    integer i

    do i = 1,numchan
       gnu(i) = GetThermoFactor(nu(i))
    end do

    if(algorithm .eq. "modified_delta_map") then
       do i = 1,numchan
          alphas(1) = gnu(i)*(nu(i)/nu_s)**beta_s
          alphas(2) = gnu(i)*(nu(i)/nu_s)**beta_s * log(nu(i)/nu_s)
          alphas(3) = gnu(i)*(nu(i)/nu_d)**beta_d / (exp(nu(i)/gamma_d)-1._dl) * (exp(nu_d/gamma_d)-1._dl)
          alphas(4) = alphas(3) * ((nu(i)/gamma_d) * exp(nu(i)/gamma_d)/(exp(nu(i)/gamma_d)-1._dl)-(nu_d/gamma_d) * exp(nu_d/gamma_d)/(exp(nu_d/gamma_d)-1._dl))
          alphas(5) = alphas(3) * log(nu(i)/nu_d)
       end do
    end if

    if(algorithm .eq. "PowLawDust") then
       do i = 1,numchan
          alphas(i) = gnu(i)*(nu(i)/nu_d)**beta_d
       end do
    end if
    if(algorithm .eq. "PowLawSync") then
       do i = 1,numchan
          alphas(i) = gnu(i)*(nu(i)/nu_s)**beta_s
       end do
    end if
    if(algorithm .eq. "PowLawDustSync") then
       do i = 1,numchan
          alphas(i) = gnu(i)*(nu(i)/nu_d)**beta_d
       end do
    end if
  end subroutine CalcAlphas
  
    function GetThermoFactor(nu) result(thermo)
      use Utils
      real(dl) thermo
      real(dl), intent(in) :: nu
      
      real(dl) x
      ! Convert units from the antenna temperature to the thermodynamic
      ! temperature in muK. (1000.0 is to convert from mK to muK)
      x = (nu/56.780d0) !nu is in units of GHz
      thermo=(exp(x)-1.0d0)**2.0d0/(x**2.0d0*exp(x)) * 1000.0d0 

    end function GetThermoFactor

    function CalculateMaskedskyLikelihood(cov_s,cov_n,Dmat,map) result(lnLike)
      use Utils
      use fitstools
      use omp_lib !時間計測
      real(dl), dimension(0:npix_tot_all-1,0:npix_tot_all-1), intent(in) :: cov_s,cov_n
      real(dl), dimension(0:npix_tot_all-1,0:2*npix*num_source-1), intent(in) :: Dmat
      real(dl), dimension(0:npix-1,1:nmap,1:numchan),intent(in) :: map
      real(dl) lnLike,lnlike1,lnlike2,logdet,logdet2

      !とりあえず並び替えたもの
      real(dl), allocatable, dimension(:) :: map_vector_reordered
      real(dl), allocatable, dimension(:,:) :: cov_s_reordered,cov_n_reordered
      real(dl), allocatable, dimension(:,:) :: Dmat_reordered
      
      !並び替えて、maskされてない部分だけ抜き出したもの
      real(dl), allocatable, dimension(:) :: map_vector_unmasked,map_t
      real(dl), allocatable, dimension(:,:) :: cov_tot_unmasked
      real(dl), allocatable, dimension(:,:) :: Dmat_unmasked,Dmat_unmasked_T,Dmat_tmp
      real(dl), allocatable, dimension(:) :: map_source,map_source_t
      real(dl), allocatable, dimension(:,:) :: covmat_source

      
      !testing
!!$      real(dl), dimension(1:numchan), intent(in) :: alphas
      real(dl),allocatable,dimension(:,:) :: cov_tmp,cov_tmp2
      real(dl) tmp,tmp1,tmp2
      real(dl), allocatable, dimension(:) :: map_t_2,map_source_2 ,map_source_t_2
      real(dl) lnlike1_2,lnlike2_2
      integer, parameter :: nvar_fore = 0
      !----

      
      integer iK,iM,iN
      integer i,j,k


      allocate(map_vector_reordered(0:npix_tot_all-1))
      allocate(cov_s_reordered(0:npix_tot_all-1,0:npix_tot_all-1))
      allocate(cov_n_reordered(0:npix_tot_all-1,0:npix_tot_all-1))      
      allocate(Dmat_reordered(0:npix_tot_all-1,0:2*npix*num_source-1))
      allocate(map_vector_unmasked(0:npix_unmasked_tot_all-1))
      allocate(cov_tot_unmasked(0:npix_unmasked_tot_all-1,0:npix_unmasked_tot_all-1))

      allocate(Dmat_unmasked(0:npix_unmasked_tot_all-1,0:npix_unmasked_source_all-1))
      allocate(Dmat_unmasked_T(0:npix_unmasked_source_all-1,0:npix_unmasked_tot_all-1))
      allocate(Dmat_tmp(0:npix_unmasked_tot_all-1,0:npix_unmasked_source_all-1))
      allocate(map_t(0:npix_unmasked_tot_all-1))
      allocate(map_source(0:npix_unmasked_source_all-1))
      allocate(covmat_source(0:npix_unmasked_source_all-1,0:npix_unmasked_source_all-1))      
      allocate(map_source_t(0:npix_unmasked_source_all-1))
      
      allocate(map_t_2(0:npix_unmasked_tot_all-1)) !testing
      allocate(map_source_2(0:npix_unmasked_source_all-1)) !testing
      allocate(map_source_t_2(0:npix_unmasked_source_all-1)) !testing
      allocate(cov_tmp(0:npix_unmasked_tot_all-1,0:npix_unmasked_tot_all-1))!testing
      allocate(cov_tmp2(0:npix_unmasked_tot_all-1,0:npix_unmasked_tot_all-1))!testing      

      !maskされているのを後ろに持ってきて、cov_n_reorderedに格納
      call ReorderCovMatrix(cov_n,cov_n_reordered,reorder_vector_tot_all)
      call ReorderCovMatrix(cov_s,cov_s_reordered,reorder_vector_tot_all)
      call ReorderDMatrix(Dmat,Dmat_reordered,reorder_vector_tot_all)
      
      !CMB用のバンドのmapを一つにまとめ、maskされているのを後ろに持ってくる
      call SetMapVectorReordered(map,map_vector_reordered) 

      
      !maskされていない部分だけ抜き出す
      map_vector_unmasked(0:npix_unmasked_tot_all-1) = map_vector_reordered(0:npix_unmasked_tot_all-1)
      cov_tot_unmasked(0:npix_unmasked_tot_all-1,0:npix_unmasked_tot_all-1)=&
           + cov_s_reordered(0:npix_unmasked_tot_all-1,0:npix_unmasked_tot_all-1)&
           + cov_n_reordered(0:npix_unmasked_tot_all-1,0:npix_unmasked_tot_all-1)
      Dmat_unmasked(0:npix_unmasked_tot_all-1,0:npix_unmasked_source_all-1) = Dmat_reordered(0:npix_unmasked_tot_all-1,0:npix_unmasked_source_all-1)
      Dmat_unmasked_T = transpose(Dmat_unmasked)

      !DGEMMの使い方
      !alphas(iM,iN)=invmat(iM,iK)*alpha0(iK,iN)
      !call DGEMM('N','N',iM,iN,iK,1.0d0,invmat,iM,alpha0,iK,0.0d0,alphas,iM)
      cov_tmp2 = cov_tot_unmasked
      cov_tot_unmasked = InvCovMatrix(cov_tot_unmasked,logdet) !cov_tot_unmasked is now (S+N)^{-1}

      iM = npix_unmasked_tot_all
      iK = npix_unmasked_tot_all
      iN = 1
      call DGEMM('N','N',iM,iN,iK,1.0d0,cov_tot_unmasked,iM,map_vector_unmasked,iK,0.0d0,map_t,iM)

      lnLike1 = dot_product(map_vector_unmasked,map_t) !first line of Eq.(61)

      !map_t = (S+N)^{-1} m = vecter(npix_unmasked_tot_all)
      !D^T(npix_unmasked_source_all,npix_unmasked_tot_all)
      !map_source(npix_unmasked_source_all) = D^T(npix_unmasked_source_all,npix_unmasked_tot_all)*map_t(npix_unmasked_tot_all,1)
      iM = npix_unmasked_source_all
      iK = npix_unmasked_tot_all
      iN = 1
      call DGEMM('N','N',iM,iN,iK,1.0d0,Dmat_unmasked_T,iM,map_t,iK,0.0d0,map_source,iM) !map_source = D^T * (S+N)^{-1} * m

      !D^T (S+N)^{-1} D
      !D^T(2*npix*(nvar_fore+2),npix_unmasked_tot_all)*(S+N)^{-1}(npix_unmasked_tot_all,npix_unmasked_tot_all)*D(npix_unmasked_tot_all,2*npix*(nvar_fore+2))
      iM = npix_unmasked_tot_all
      iK = npix_unmasked_tot_all
      iN = npix_unmasked_source_all!2*npix*(nvar_fore+2)
      call DGEMM('N','N',iM,iN,iK,1.0d0,cov_tot_unmasked,iM,Dmat_unmasked,iK,0.0d0,Dmat_tmp,iM)

      
      iM = npix_unmasked_source_all!2*npix*(nvar_fore+2)
      iK = npix_unmasked_tot_all
      iN = npix_unmasked_source_all!2*npix*(nvar_fore+2)
      call DGEMM('N','N',iM,iN,iK,1.0d0,Dmat_unmasked_T,iM,Dmat_tmp,iK,0.0d0,covmat_source,iM) !covmat_source = D^T * (S+N)^{-1} * D

      
      !covmat_source = D^T (S+N)^{-1} D
      covmat_source = InvCovMatrix(covmat_source,logdet2) !covmat_source is now covmat_source^{-1}

      
      !map_source_t(2*npix*(nvar_fore+2) = covmat_source(2*npix*(nvar_fore+2),2*npix*(nvar_fore+2)*map_source(2*npix*(nvar_fore+2)
      iM = npix_unmasked_source_all!2*npix*(nvar_fore+2)
      iK = npix_unmasked_source_all!2*npix*(nvar_fore+2)
      iN = 1
      call DGEMM('N','N',iM,iN,iK,1.0d0,covmat_source,iM,map_source,iK,0.0d0,map_source_t,iM) !map_source_t = [D^T * (S+N)^{-1} * D]^{-1} D^T * (S+N)^{-1} * m

      lnlike2 = dot_product(map_source,map_source_t) !second line in Eq.(61)

      lnlike = lnlike1/2 + logdet/2 -lnlike2/2 &
              +logdet2/2 !<-newly added by Minami-san

      
      deallocate(map_vector_reordered)
      deallocate(cov_s_reordered)
      deallocate(cov_n_reordered)      
      deallocate(Dmat_reordered)
      deallocate(map_vector_unmasked)
      deallocate(cov_tot_unmasked)

      deallocate(Dmat_unmasked)
      deallocate(Dmat_unmasked_T)
      deallocate(Dmat_tmp)
      deallocate(map_t)
      deallocate(map_source)
      deallocate(covmat_source)      
      deallocate(map_source_t)
      
      deallocate(map_t_2) !testing
      deallocate(map_source_2) !testing
      deallocate(map_source_t_2) !testing
      deallocate(cov_tmp)!testing
      deallocate(cov_tmp2)!testing
      
    end function CalculateMaskedskyLikelihood




    function CalculateMaskedskyLikelihood3(cov_s,cov_n,Dmat,alphas,map) result(lnLike)
      use Utils
      use fitstools
      use omp_lib !時間計測
      real(dl), dimension(0:npix_tot_all-1,0:npix_tot_all-1), intent(in) :: cov_s,cov_n
      real(dl), dimension(0:npix_tot_all-1,0:2*npix*num_source-1), intent(in) :: Dmat
      real(dl), dimension(0:npix-1,1:nmap,1:numchan),intent(in) :: map
      real(dl) lnLike,lnlike1,lnlike2,logdet

      !とりあえず並び替えたもの
      real(dl), allocatable, dimension(:) :: map_vector_reordered
      real(dl), allocatable, dimension(:,:) :: cov_s_reordered,cov_n_reordered
      real(dl), allocatable, dimension(:,:) :: Dmat_reordered
      
      !並び替えて、maskされてない部分だけ抜き出したもの
      real(dl), allocatable, dimension(:) :: map_vector_unmasked,map_t
      real(dl), allocatable, dimension(:,:) :: cov_tot_unmasked
      real(dl), allocatable, dimension(:,:) :: Dmat_unmasked,Dmat_unmasked_T,Dmat_tmp
      real(dl), allocatable, dimension(:) :: map_source,map_source_t
      real(dl), allocatable, dimension(:,:) :: covmat_source

      
      !testing
            real(dl), dimension(1:numchan), intent(in) :: alphas
      real(dl),allocatable,dimension(:,:) :: cov_tmp,cov_tmp2
      real(dl) tmp,tmp1,tmp2
      real(dl), allocatable, dimension(:) :: map_t_2,map_source_2 ,map_source_t_2
      real(dl) lnlike1_2,lnlike2_2
      !----

      
      integer iK,iM,iN
      integer i,j,k


      allocate(map_vector_reordered(0:npix_tot_all-1))
      allocate(cov_s_reordered(0:npix_tot_all-1,0:npix_tot_all-1))
      allocate(cov_n_reordered(0:npix_tot_all-1,0:npix_tot_all-1))      
      allocate(Dmat_reordered(0:npix_tot_all-1,0:2*npix*num_source-1))
      allocate(map_vector_unmasked(0:npix_unmasked_tot_all-1))
      allocate(cov_tot_unmasked(0:npix_unmasked_tot_all-1,0:npix_unmasked_tot_all-1))
      allocate(Dmat_unmasked(0:npix_unmasked_tot_all-1,0:npix_unmasked_source_all-1))
      allocate(Dmat_unmasked_T(0:npix_unmasked_source_all-1,0:npix_unmasked_tot_all-1))
      allocate(Dmat_tmp(0:npix_unmasked_tot_all-1,0:npix_unmasked_source_all-1))
      allocate(map_t(0:npix_unmasked_tot_all-1))
      allocate(map_source(0:npix_unmasked_source_all-1))
      allocate(covmat_source(0:npix_unmasked_source_all-1,0:npix_unmasked_source_all-1))      
      allocate(map_source_t(0:npix_unmasked_source_all-1))
      
      allocate(map_t_2(0:npix_unmasked_tot_all-1)) !testing
      allocate(map_source_2(0:npix_unmasked_source_all-1)) !testing
      allocate(map_source_t_2(0:npix_unmasked_source_all-1)) !testing
      allocate(cov_tmp(0:npix_unmasked_tot_all-1,0:npix_unmasked_tot_all-1))!testing
      allocate(cov_tmp2(0:npix_unmasked_tot_all-1,0:npix_unmasked_tot_all-1))!testing      

      !maskされているのを後ろに持ってきて、cov_n_reorderedに格納
      call ReorderCovMatrix(cov_n,cov_n_reordered,reorder_vector_tot_all)
      call ReorderCovMatrix(cov_s,cov_s_reordered,reorder_vector_tot_all)
      call ReorderDMatrix(Dmat,Dmat_reordered,reorder_vector_tot_all)
      
      !CMB用のバンドのmapを一つにまとめ、maskされているのを後ろに持ってくる
      call SetMapVectorReordered(map,map_vector_reordered) 

      
      !maskされていない部分だけ抜き出す
      map_vector_unmasked(0:npix_unmasked_tot_all-1) = map_vector_reordered(0:npix_unmasked_tot_all-1)
      cov_tot_unmasked(0:npix_unmasked_tot_all-1,0:npix_unmasked_tot_all-1)=&
           + cov_s_reordered(0:npix_unmasked_tot_all-1,0:npix_unmasked_tot_all-1)&
           + cov_n_reordered(0:npix_unmasked_tot_all-1,0:npix_unmasked_tot_all-1)
      Dmat_unmasked(0:npix_unmasked_tot_all-1,0:npix_unmasked_source_all-1) = Dmat_reordered(0:npix_unmasked_tot_all-1,0:npix_unmasked_source_all-1)
      Dmat_unmasked_T = transpose(Dmat_unmasked)

!!$      write(*,*)"npix_unmasked",npix_unmasked
!!$      write(*,*)"npix_unmasked_tot_all",npix_unmasked_tot_all
!!$      write(*,*)"npix_unmasked_source_all",npix_unmasked_source_all
      
      !DGEMMの使い方
      !alphas(iM,iN)=invmat(iM,iK)*alpha0(iK,iN)
      !call DGEMM('N','N',iM,iN,iK,1.0d0,invmat,iM,alpha0,iK,0.0d0,alphas,iM)
      cov_tmp2 = cov_tot_unmasked
      cov_tot_unmasked = InvCovMatrix(cov_tot_unmasked,logdet) !cov_tot_unmasked is now (S+N)^{-1}

      iM = npix_unmasked_tot_all
      iK = npix_unmasked_tot_all
      iN = 1
      call DGEMM('N','N',iM,iN,iK,1.0d0,cov_tot_unmasked,iM,map_vector_unmasked,iK,0.0d0,map_t_2,iM)
      map_t = matmul(cov_tot_unmasked,map_vector_unmasked)

!!$      do i=0,npix_unmasked_tot_all-1
!!$         do j=0,npix_unmasked_tot_all-1
!!$            write(92,'(2I5,E25.15)') i,j,cov_tot_unmasked(i,j)
!!$         end do
!!$      end do
!!$      do i=0,npix_unmasked_tot_all-1
!!$         write(91,'(I5,3E25.15)') i,map_t(i),map_vector_unmasked(i)
!!$      end do
      
      lnLike1 = dot_product(map_vector_unmasked,map_t)
      lnLike1_2 = dot_product(map_vector_unmasked,map_t_2) !testing

!!$      write(*,*)"lnlike1 comparison:",lnLike1,lnLike1_2
!!$      stop
      
      !map_t = (S+N)^{-1} m = vecter(npix_unmasked_tot_all)
      !D^T(npix_unmasked_source_all,npix_unmasked_tot_all)
      !map_source(npix_unmasked_source_all) = D^T(npix_unmasked_source_all,npix_unmasked_tot_all)*map_t(npix_unmasked_tot_all,1)
!!$      iM = npix_unmasked_source_all
!!$      iK = npix_unmasked_tot_all
!!$      iN = 1
!!$      call DGEMM('N','N',iM,iN,iK,1.0d0,Dmat_unmasked_T,iM,map_t,iK,0.0d0,map_source,iM)
      map_source = matmul(Dmat_unmasked_T,map_t)
      map_source_2 =matmul(Dmat_unmasked_T,map_t_2) !testing

!!$      do i=0,npix_unmasked_source_all-1
!!$         write(92,'(I5,3E25.15)') i,map_source(i),map_source_2(i),map_source(i)-map_source_2(i)
!!$      end do
      
!!$      write(*,*)"map_source comparison:",dot_product(map_source,map_source),dot_product(map_source_2,map_source_2)


      !D^T (S+N)^{-1} D
      !D^T(2*npix*(nvar_fore+2),npix_unmasked_tot_all)*(S+N)^{-1}(npix_unmasked_tot_all,npix_unmasked_tot_all)*D(npix_unmasked_tot_all,2*npix*(nvar_fore+2))
!!$      iM = npix_unmasked_tot_all
!!$      iK = npix_unmasked_tot_all
!!$      iN = 2*npix*(nvar_fore+2)
!!$      call DGEMM('N','N',iM,iN,iK,1.0d0,cov_tot_unmasked,iM,Dmat_unmasked,iK,0.0d0,Dmat_tmp,iM)
      Dmat_tmp = matmul(cov_tot_unmasked,Dmat_unmasked)
      
!!$      iM = 2*npix*(nvar_fore+2)
!!$      iK = npix_unmasked_tot_all
!!$      iN = 2*npix*(nvar_fore+2)
!!$      call DGEMM('N','N',iM,iN,iK,1.0d0,Dmat_unmasked_T,iM,Dmat_tmp,iK,0.0d0,covmat_source,iM)
      covmat_source = matmul(Dmat_unmasked_T,Dmat_tmp)

!!$      !testing
!!$      do i=0,npix_unmasked_source_all-1
!!$         do j=0,npix_unmasked_source_all-1
!!$            write(70,"(2I5,E25.15)") i,j,covmat_source(i,j)
!!$         end do
!!$         write(70,*) " "
!!$      end do
!!$      stop
!!$
!!$      do i=0,npix_unmasked_source_all-1
!!$         do j=0,npix_unmasked_tot_all
!!$            write(71,"(2I5,2E15.5)") i,j,Dmat_unmasked_T(i,j),Dmat_unmasked(j,i)
!!$         end do
!!$         write(71,*) " "
!!$      end do

      covmat_source = InvCovMatrix(covmat_source) !covmat_source is now covmat_source^{-1}

!!$      do i=0,npix_unmasked_source_all-1
!!$         do j=0,npix_unmasked_source_all-1
!!$            write(71,"(2I5,E25.15)") i,j,covmat_source(i,j)
!!$         end do
!!$         write(71,*) " "
!!$      end do

      
      !map_source_t(2*npix*(nvar_fore+2) = covmat_source(2*npix*(nvar_fore+2),2*npix*(nvar_fore+2)*map_source(2*npix*(nvar_fore+2)

!!$      iM = 2*npix*(nvar_fore+2)
!!$      iK = 2*npix*(nvar_fore+2)
!!$      iN = 1
!!$      call DGEMM('N','N',iM,iN,iK,1.0d0,covmat_source,iM,map_source,iK,0.0d0,map_source_t,iM)

      map_source_t = matmul(covmat_source,map_source)
      lnlike2 = dot_product(map_source,map_source_t)

      map_source_t_2 = matmul(covmat_source,map_source_2)!testing
      lnlike2_2 = dot_product(map_source_2,map_source_t_2)!testing

!!$      write(*,*)"map_source_t comparison:",dot_product(map_source_t,map_source_t),dot_product(map_source_t_2,map_source_t_2)
      
!!$      do i=0,npix_unmasked_source_all-1
!!$         write(83,'(I5,3E25.15)') i,map_source_t(i),map_source_t_2(i),map_source_t(i)-map_source_t_2(i)
!!$         write(84,'(I5,3E25.15)') i,map_source(i),map_source_2(i),map_source(i)-map_source_2(i)
!!$      end do
      
!!$      write(*,*)"lnlike2 comparison:",lnlike2,lnlike2_2

!!$      tmp = 0._dl
!!$      do i=0,npix_unmasked_source_all-1
!!$         tmp = tmp + covmat_source(30,i)*map_source(i)
!!$         write(85,'(I5,3E25.15)') i,covmat_source(30,i),map_source(i),tmp
!!$      end do
!!$      write(*,*) "tmp, map_source_t(30)=",tmp,map_source_t(30)
!!$      stop
      
!!$      do i=0,npix_unmasked_source_all-1
!!$         write(72,*) i,map_source(i),map_source_t(i)
!!$      end do
!!$      do i=0,npix_unmasked_source_all-1
!!$         do j=0,npix_unmasked_source_all-1
!!$            write(71,*) i,j,covmat_source(i,j)
!!$         end do
!!$      end do
!!$
!!$      do i=0,npix_unmasked_source_all-1
!!$         write(70,*) i,map_source_t(i),map_source_t_2(i)
!!$      end do
!!$      stop
      
      !DGEMMの使い方
      !alphas(iM,iN)=invmat(iM,iK)*alpha0(iK,iN)
      !call DGEMM('N','N',iM,iN,iK,1.0d0,invmat,iM,alpha0,iK,0.0d0,alphas,iM)

!!$      write(*,*)"lnlike1,lnlike2,logdet:"
!!$      write(*,'(4E25.15)') lnlike1,lnlike2,logdet,lnlike1-lnlike2
      lnlike = lnlike1/2 +logdet/2-lnlike2/2
!!$      stop
      
      deallocate(map_vector_reordered)
      deallocate(cov_s_reordered)
      deallocate(cov_n_reordered)      
      deallocate(Dmat_reordered)
      deallocate(map_vector_unmasked)
      deallocate(cov_tot_unmasked)
      deallocate(Dmat_unmasked)
      deallocate(Dmat_unmasked_T)
      deallocate(Dmat_tmp)
      deallocate(map_t)
      deallocate(map_source)
      deallocate(covmat_source)      
      deallocate(map_source_t)
      
      deallocate(map_t_2) !testing
      deallocate(map_source_2) !testing
      deallocate(map_source_t_2) !testing
      deallocate(cov_tmp)!testing
      deallocate(cov_tmp2)!testing
      
    end function CalculateMaskedskyLikelihood3

    function CalculateMaskedskyLikelihood4(cov_s,cov_n,Dmat,alphas,map) result(lnLike)
      use Utils
      use fitstools
      use omp_lib !時間計測
      real(dl), dimension(0:npix_tot_all-1,0:npix_tot_all-1), intent(in) :: cov_s,cov_n
      real(dl), dimension(0:npix_tot_all-1,0:2*npix*num_source-1), intent(in) :: Dmat
      real(dl), dimension(0:npix-1,1:nmap,1:numchan),intent(in) :: map
      real(dl) lnLike,lnlike1,lnlike2,logdet
      !とりあえず並び替えたもの
      real(dl), allocatable, dimension(:) :: map_vector_reordered
      real(dl), allocatable, dimension(:,:) :: cov_s_reordered,cov_n_reordered
      real(dl), allocatable, dimension(:,:) :: Dmat_reordered
      
      !並び替えて、maskされてない部分だけ抜き出したもの
      real(dl), allocatable, dimension(:) :: map_vector_unmasked,map_t
      real(dl), allocatable, dimension(:,:) :: cov_tot_unmasked,cov_tot_unmasked2,cov_tot_unmasked3
      real(dl), allocatable, dimension(:,:) :: Dmat_unmasked,Dmat_unmasked_T,Dmat_tmp,Dmat_tmp_T
      real(dl), allocatable, dimension(:) :: map_source,map_source_t
      real(dl), allocatable, dimension(:,:) :: covmat_source

      real(dl), allocatable, dimension(:,:) :: DtransCov,DtransCov_T
      
      integer iK,iM,iN
      integer i,j


      !testing variables
      real(dl) g1D1,g2D2,chisq!,CMB_Q,CMB_U
      real(dl), dimension(0:npix_unmasked-1) :: CMB_Q,CMB_U
      real(dl), dimension(0:2*npix_unmasked-1) :: CMBmap
      real(dl), dimension(1:numchan), intent(in) :: alphas
      real(dl), dimension(0:2*npix_unmasked-1,0:2*npix_unmasked-1) :: cov_CMB
      !-----
      
      allocate(map_vector_reordered(0:npix_tot_all-1))
      allocate(cov_s_reordered(0:npix_tot_all-1,0:npix_tot_all-1))
      allocate(cov_n_reordered(0:npix_tot_all-1,0:npix_tot_all-1))      
      allocate(Dmat_reordered(0:npix_tot_all-1,0:2*npix*num_source-1))
      allocate(map_vector_unmasked(0:npix_unmasked_tot_all-1))
      allocate(cov_tot_unmasked(0:npix_unmasked_tot_all-1,0:npix_unmasked_tot_all-1))
      allocate(cov_tot_unmasked2(0:npix_unmasked_tot_all-1,0:npix_unmasked_tot_all-1))
      allocate(cov_tot_unmasked3(0:npix_unmasked_tot_all-1,0:npix_unmasked_tot_all-1))
      allocate(Dmat_unmasked(0:npix_unmasked_tot_all-1,0:npix_unmasked_source_all-1))
      allocate(Dmat_unmasked_T(0:npix_unmasked_source_all-1,0:npix_unmasked_tot_all-1))
      allocate(Dmat_tmp(0:npix_unmasked_tot_all-1,0:npix_unmasked_source_all-1))
      allocate(Dmat_tmp_T(0:npix_unmasked_source_all-1,0:npix_unmasked_tot_all-1))
      allocate(DtransCov(0:npix_unmasked_source_all-1,0:npix_unmasked_tot_all-1))
      allocate(DtransCov_T(0:npix_unmasked_tot_all-1,0:npix_unmasked_source_all-1))
      
      allocate(map_t(0:npix_unmasked_tot_all-1))
      allocate(map_source(0:npix_unmasked_source_all-1))
      allocate(covmat_source(0:npix_unmasked_source_all-1,0:npix_unmasked_source_all-1))      
      allocate(map_source_t(0:npix_unmasked_source_all-1))
      
      
      !maskされているのを後ろに持ってきて、cov_n_reorderedに格納
      call ReorderCovMatrix(cov_n,cov_n_reordered,reorder_vector_tot_all)
      call ReorderCovMatrix(cov_s,cov_s_reordered,reorder_vector_tot_all)
      call ReorderDMatrix(Dmat,Dmat_reordered,reorder_vector_tot_all)
      
      !CMB用のバンドのmapを一つにまとめ、maskされているのを後ろに持ってくる
      call SetMapVectorReordered(map,map_vector_reordered) 

      
      !maskされていない部分だけ抜き出す
      map_vector_unmasked(0:npix_unmasked_tot_all-1) = map_vector_reordered(0:npix_unmasked_tot_all-1)
      cov_tot_unmasked(0:npix_unmasked_tot_all-1,0:npix_unmasked_tot_all-1)=&
           + cov_s_reordered(0:npix_unmasked_tot_all-1,0:npix_unmasked_tot_all-1)&
           + cov_n_reordered(0:npix_unmasked_tot_all-1,0:npix_unmasked_tot_all-1)
      Dmat_unmasked(0:npix_unmasked_tot_all-1,0:npix_unmasked_source_all-1) = Dmat_reordered(0:npix_unmasked_tot_all-1,0:npix_unmasked_source_all-1)
      Dmat_unmasked_T = transpose(Dmat_unmasked)


      !testing
!!$      write(*,*)"npix_unmasked_tot_all=",npix_unmasked_tot_all
!!$      write(*,*)"npix_unmasked_source_all=",npix_unmasked_source_all
      
      !DGEMMの使い方
      !alphas(iM,iN)=invmat(iM,iK)*alpha0(iK,iN)
      !call DGEMM('N','N',iM,iN,iK,1.0d0,invmat,iM,alpha0,iK,0.0d0,alphas,iM)
      cov_tot_unmasked = InvCovMatrix(cov_tot_unmasked,logdet) !cov_tot_unmasked is now (S+N)^{-1}

      !calc D^T(S+N)^{-1}
      DtransCov = matmul(Dmat_unmasked_T,cov_tot_unmasked)
      !calc ((S+N)^{-1})^T D
      DtransCov_T = transpose(DtransCov)
      !calc D^T(S+N)^{-1} D
      covmat_source = matmul(DtransCov,Dmat_unmasked)
      covmat_source = InvCovMatrix(covmat_source) !cov_source is now [D^T(S+N)^{-1} D]^{-1}

      !calc DtransCov_T [D^T(S+N)^{-1} D]^{-1} D^T(S+N)^{-1}
      Dmat_tmp_T = matmul(covmat_source,DtransCov)
      cov_tot_unmasked2 = matmul(DtransCov_T, Dmat_tmp_T)

      cov_tot_unmasked3 = cov_tot_unmasked - cov_tot_unmasked2
!!$      do i=0,npix_unmasked_tot_all-1
!!$         do j=0,npix_unmasked_tot_all-1
!!$            write(50,'(2I5,3E25.15)') i,j,cov_tot_unmasked(i,j),cov_tot_unmasked2(i,j),cov_tot_unmasked3(i,j)
!!$         end do
!!$         write(50,*)""
!!$      end do
!!$      stop

      !testing
      !alpha_1 = g_{nu1} D_{nu1}
      !alpha_2 = g_{nu2} D_{nu2}

!!$      write(*,*) "cov_s(0,0)=" ,cov_s_reordered(0,0)
!!$      write(*,*) "cov_s(0,1)=" ,cov_s_reordered(0,1)
!!$      write(*,*) "cov_s(1,1)=" ,cov_s_reordered(1,1)
!!$      write(*,*) "cov_cmb(0,0)=" ,cov_cmb(0,0)
!!$      write(*,*) "cov_cmb(0,1)=" ,cov_cmb(0,1)
!!$      write(*,*) "cov_cmb(1,0)=" ,cov_cmb(1,0)
!!$      write(*,*) "cov_cmb(1,1)=" ,cov_cmb(1,1)
!!$      write(*,*) "npix_unmasked=",npix_unmasked
      
!!$      g1D1=alphas(1)
!!$      g2D2=alphas(2)
!!$      cov_cmb(0:2*npix_unmasked-1,0:2*npix_unmasked-1)=cov_s_reordered(0:2*npix_unmasked-1,0:2*npix_unmasked-1)&
!!$                      +(g2D2*g2D2*cov_n_reordered(0:2*npix_unmasked-1,0:2*npix_unmasked-1)+g1D1*g1D1*cov_n_reordered(2*npix_unmasked:4*npix_unmasked-1,2*npix_unmasked:4*npix_unmasked-1))/(g2D2-g1D1)/(g2D2-g1D1)
!!$
!!$
!!$      do j=0,2*npix_unmasked-1
!!$         write(101,*) cov_cmb(0,j)
!!$      end do
!!$
!!$
!!$      cov_cmb = InvCovMatrix(cov_cmb)
      
!!$      do i=0,npix_unmasked-1
!!$         CMB_Q(i) = g2D2*map(i,2,1)-g1D1*map(i,2,2)
!!$         CMB_U(i) = g2D2*map(i,3,1)-g1D1*map(i,3,2)
!!$         CMB_Q(i) = CMB_Q(i)/(g2D2-g1D1)
!!$         CMB_U(i) = CMB_U(i)/(g2D2-g1D1)
!!$         CMBmap(i) = CMB_Q(i)
!!$         CMBmap(i+npix_unmasked) = CMB_U(i)
!!$      end do

!!$      CMBmap(0:2*npix_unmasked-1) = g2D2*map_vector_unmasked(0:2*npix_unmasked-1)-g1D1*map_vector_unmasked(2*npix_unmasked:4*npix_unmasked-1)
!!$      CMBmap(0:2*npix_unmasked-1) = CMBmap(0:2*npix_unmasked-1)/(g2D2-g1D1)
      
!!$      write(*,*) "chisq for single pixel=",CMB_Q*cov_cmb(0,0)*CMB_Q + CMB_U*cov_cmb(1,1)*CMB_U + 2*CMB_Q*cov_cmb(0,1)*CMB_U

!!$      chisq=0._dl
!!$      do i=0,2*npix_unmasked-1
!!$         do j=0,2*npix_unmasked-1
!!$            chisq=chisq+CMBmap(i)*cov_cmb(i,j)*CMBmap(j)
!!$         end do
!!$      end do
!!$      write(*,*)"chisq=",chisq

      
!!$      write(*,*) "map(0,2,1)=",map(0,2,1) !Q
!!$      write(*,*) "map(0,3,2)=",map(0,2,2) !Q
!!$      write(*,*) "map(0,2,1)=",map(0,3,1) !U
!!$      write(*,*) "map(0,3,2)=",map(0,3,2) !U
!!$      write(*,*) "invcov_cmb(0,0)=" ,cov_cmb(0,0)
!!$      write(*,*) "invcov_cmb(0,1)=" ,cov_cmb(0,1)
!!$      write(*,*) "invcov_cmb(1,1)=" ,cov_cmb(1,1)
!!$      write(*,*) "CMB_Q, CMB_U=",CMB_Q,CMB_U

!!$
!!$      write(*,*) "map_vector_unmasked(0:3)=" , map_vector_unmasked(0:3)
      
      !testing
      map_t = matmul(cov_tot_unmasked,map_vector_unmasked)
      lnlike1 = dot_product(map_vector_unmasked,map_t)
      write(*,*)"breakdown lnlike1:",lnlike1
      map_t = matmul(cov_tot_unmasked2,map_vector_unmasked)
      lnlike2 = dot_product(map_vector_unmasked,map_t)
      write(*,*)"breakdown lnlike2:",lnlike2
      write(*,*)"expected total lnlike:",lnlike1-lnlike2
      !---
      map_t = matmul(cov_tot_unmasked3,map_vector_unmasked)
      lnLike = dot_product(map_vector_unmasked,map_t)
      write(*,*)"         total lnlike:",lnLike
      stop
      
      !DGEMMの使い方
      !alphas(iM,iN)=invmat(iM,iK)*alpha0(iK,iN)
      !call DGEMM('N','N',iM,iN,iK,1.0d0,invmat,iM,alpha0,iK,0.0d0,alphas,iM)
      lnlike = lnlike/2 +logdet/2
      lnlike = chisq/2 +logdet/2!testing (determinantの次元がdeltamapとgeneral deltamapで異なる)

      
      deallocate(map_vector_reordered)
      deallocate(cov_s_reordered)
      deallocate(cov_n_reordered)      
      deallocate(Dmat_reordered)
      deallocate(map_vector_unmasked)
      deallocate(cov_tot_unmasked)
      deallocate(cov_tot_unmasked2)
      deallocate(cov_tot_unmasked3)
      deallocate(Dmat_unmasked)
      deallocate(Dmat_unmasked_T)
      deallocate(Dmat_tmp)
      deallocate(Dmat_tmp_T)
      deallocate(DtransCov)
      deallocate(DtransCov_T)
      
      deallocate(map_t)
      deallocate(map_source)
      deallocate(covmat_source)      
      deallocate(map_source_t)

!!$      deallocate(map_vector_reordered)
!!$      deallocate(cov_s_reordered)
!!$      deallocate(cov_n_reordered)      
!!$      deallocate(Dmat_reordered)
!!$      deallocate(map_vector_unmasked)
!!$      deallocate(cov_tot_unmasked)
!!$      deallocate(Dmat_unmasked)
!!$      deallocate(Dmat_unmasked_T)
!!$      deallocate(Dmat_tmp)
!!$      deallocate(map_t)
!!$      deallocate(map_source)
!!$      deallocate(covmat_source)      
!!$      deallocate(map_source_t)
      

    end function CalculateMaskedskyLikelihood4


    function CalculateMaskedskyLikelihood5(cov_s,cov_n,Dmat,map) result(lnLike)
      use Utils
      use fitstools
      use omp_lib !時間計測
      real(dl), dimension(0:npix_tot_all-1,0:npix_tot_all-1), intent(in) :: cov_s,cov_n
      real(dl), dimension(0:npix_tot_all-1,0:2*npix*num_source-1), intent(in) :: Dmat
      real(dl), dimension(0:npix-1,1:nmap,1:numchan),intent(in) :: map
      real(dl) lnLike,lnlike1,lnlike2,logdet
      !とりあえず並び替えたもの
      real(dl), allocatable, dimension(:) :: map_vector_reordered
      real(dl), allocatable, dimension(:,:) :: cov_s_reordered,cov_n_reordered
      real(dl), allocatable, dimension(:,:) :: Dmat_reordered
      
      !並び替えて、maskされてない部分だけ抜き出したもの
      real(dl), allocatable, dimension(:) :: map_vector_unmasked,map_t
      real(dl), allocatable, dimension(:,:) :: cov_tot_unmasked,cov_tot_unmasked2,cov_tot_unmasked3
      real(dl), allocatable, dimension(:,:) :: Dmat_unmasked,Dmat_unmasked_T,Dmat_tmp,Dmat_tmp_T
      real(dl), allocatable, dimension(:) :: map_source,map_source_t
      real(dl), allocatable, dimension(:,:) :: covmat_source

      real(dl), allocatable, dimension(:,:) :: DtransCov,DtransCov_T
      
      integer iK,iM,iN
      integer i,j


      allocate(map_vector_reordered(0:npix_tot_all-1))
      allocate(cov_s_reordered(0:npix_tot_all-1,0:npix_tot_all-1))
      allocate(cov_n_reordered(0:npix_tot_all-1,0:npix_tot_all-1))      
      allocate(Dmat_reordered(0:npix_tot_all-1,0:2*npix*num_source-1))
      allocate(map_vector_unmasked(0:npix_unmasked_tot_all-1))
      allocate(cov_tot_unmasked(0:npix_unmasked_tot_all-1,0:npix_unmasked_tot_all-1))
      allocate(cov_tot_unmasked2(0:npix_unmasked_tot_all-1,0:npix_unmasked_tot_all-1))
      allocate(cov_tot_unmasked3(0:npix_unmasked_tot_all-1,0:npix_unmasked_tot_all-1))
      allocate(Dmat_unmasked(0:npix_unmasked_tot_all-1,0:npix_unmasked_source_all-1))
      allocate(Dmat_unmasked_T(0:npix_unmasked_source_all-1,0:npix_unmasked_tot_all-1))
      allocate(Dmat_tmp(0:npix_unmasked_tot_all-1,0:npix_unmasked_source_all-1))
      allocate(Dmat_tmp_T(0:npix_unmasked_source_all-1,0:npix_unmasked_tot_all-1))
      allocate(DtransCov(0:npix_unmasked_source_all-1,0:npix_unmasked_tot_all-1))
      allocate(DtransCov_T(0:npix_unmasked_tot_all-1,0:npix_unmasked_source_all-1))
      
      allocate(map_t(0:npix_unmasked_tot_all-1))
      allocate(map_source(0:npix_unmasked_source_all-1))
      allocate(covmat_source(0:npix_unmasked_source_all-1,0:npix_unmasked_source_all-1))      
      allocate(map_source_t(0:npix_unmasked_source_all-1))
      
      
      !maskされているのを後ろに持ってきて、cov_n_reorderedに格納
      call ReorderCovMatrix(cov_n,cov_n_reordered,reorder_vector_tot_all)
      call ReorderCovMatrix(cov_s,cov_s_reordered,reorder_vector_tot_all)
      call ReorderDMatrix(Dmat,Dmat_reordered,reorder_vector_tot_all)
      
      !CMB用のバンドのmapを一つにまとめ、maskされているのを後ろに持ってくる
      call SetMapVectorReordered(map,map_vector_reordered) 

      
      !maskされていない部分だけ抜き出す
      map_vector_unmasked(0:npix_unmasked_tot_all-1) = map_vector_reordered(0:npix_unmasked_tot_all-1)
      cov_tot_unmasked(0:npix_unmasked_tot_all-1,0:npix_unmasked_tot_all-1)=&
           + cov_s_reordered(0:npix_unmasked_tot_all-1,0:npix_unmasked_tot_all-1)&
           + cov_n_reordered(0:npix_unmasked_tot_all-1,0:npix_unmasked_tot_all-1)
      Dmat_unmasked(0:npix_unmasked_tot_all-1,0:npix_unmasked_source_all-1) = Dmat_reordered(0:npix_unmasked_tot_all-1,0:npix_unmasked_source_all-1)
      Dmat_unmasked_T = transpose(Dmat_unmasked)


      !DGEMMの使い方
      !alphas(iM,iN)=invmat(iM,iK)*alpha0(iK,iN)
      !call DGEMM('N','N',iM,iN,iK,1.0d0,invmat,iM,alpha0,iK,0.0d0,alphas,iM)
      cov_tot_unmasked = InvCovMatrix(cov_tot_unmasked,logdet) !cov_tot_unmasked is now (S+N)^{-1}

      !calc D^T(S+N)^{-1}
      DtransCov = matmul(Dmat_unmasked_T,cov_tot_unmasked)
      !calc ((S+N)^{-1})^T D
      DtransCov_T = transpose(DtransCov)
      !calc D^T(S+N)^{-1} D
      covmat_source = matmul(DtransCov,Dmat_unmasked)
      covmat_source = InvCovMatrix(covmat_source) !cov_source is now [D^T(S+N)^{-1} D]^{-1}

      !calc DtransCov_T [D^T(S+N)^{-1} D]^{-1} D^T(S+N)^{-1}
      Dmat_tmp_T = matmul(covmat_source,DtransCov)
      cov_tot_unmasked2 = matmul(DtransCov_T, Dmat_tmp_T)

      cov_tot_unmasked3 = cov_tot_unmasked - cov_tot_unmasked2
      
      !testing
!!$      map_t = matmul(cov_tot_unmasked,map_vector_unmasked)
!!$      lnlike1 = dot_product(map_vector_unmasked,map_t)
!!$      write(*,*)"breakdown lnlike1:",lnlike1
!!$      map_t = matmul(cov_tot_unmasked2,map_vector_unmasked)
!!$      lnlike2 = dot_product(map_vector_unmasked,map_t)
!!$      write(*,*)"breakdown lnlike2:",lnlike2
!!$      write(*,*)"expected total lnlike:",lnlike1-lnlike2
      !---
      map_t = matmul(cov_tot_unmasked3,map_vector_unmasked)
      lnLike = dot_product(map_vector_unmasked,map_t)
      write(*,*)"         total lnlike:",lnLike
      !stop
      
      !DGEMMの使い方
      !alphas(iM,iN)=invmat(iM,iK)*alpha0(iK,iN)
      !call DGEMM('N','N',iM,iN,iK,1.0d0,invmat,iM,alpha0,iK,0.0d0,alphas,iM)
      lnlike = lnlike/2 +logdet/2

      
      deallocate(map_vector_reordered)
      deallocate(cov_s_reordered)
      deallocate(cov_n_reordered)      
      deallocate(Dmat_reordered)
      deallocate(map_vector_unmasked)
      deallocate(cov_tot_unmasked)
      deallocate(cov_tot_unmasked2)
      deallocate(cov_tot_unmasked3)
      deallocate(Dmat_unmasked)
      deallocate(Dmat_unmasked_T)
      deallocate(Dmat_tmp)
      deallocate(Dmat_tmp_T)
      deallocate(DtransCov)
      deallocate(DtransCov_T)
      
      deallocate(map_t)
      deallocate(map_source)
      deallocate(covmat_source)      
      deallocate(map_source_t)


    end function CalculateMaskedskyLikelihood5

    
    
    subroutine SetMapVectorReordered(map_in,map_out)
      real(dl),dimension(0:npix-1,1:nmap,1:numchan),intent(in) :: map_in
      real(dl),dimension(0:npix_tot_all-1),intent(inout) :: map_out

      real(dl),allocatable,dimension(:) :: map_tmp
      integer i,j,index_im,index_i,index_ip

      !CMB用のバンドのmapを一つにまとめる
      allocate(map_tmp(0:npix_tot_all-1))
      do i=1,numchan
         j=2*(i-1)
         index_im = j*npix
         index_i  = (j+1)*npix
         index_ip = (j+2)*npix
         map_tmp(index_im:index_i-1) = map_in(0:npix-1,2,i)  !Q成分 (CMB(i))
         map_tmp(index_i:index_ip-1) = map_in(0:npix-1,3,i)  !U成分 (CMB(i))
      end do

      !maskされてないピクセルを配列の前方へ移動
      do i=0,npix_tot_all-1
         map_out(i) = map_tmp(reorder_vector_tot_all(i))
      end do
      deallocate(map_tmp)

    end subroutine SetMapVectorReordered


       function NoisePerPix(sig_noise,npixel) result(noise_per_pix)
         use Utils
         integer npixel
         real(dl) sig_noise
         real(dl) noise_per_pix
         real(dl) solid_angle
         
         solid_angle = (fourpi/npixel)
         noise_per_pix = ((pi/10800._dl)*sig_noise)/sqrt(solid_angle)
         
       end function NoisePerPix

       subroutine SetNoiseCovMatrices
         use Utils
         implicit none
         !nuCMB,nu1,nu2,nu3,nu4,nu5,nu6からファイル名を推測して代入する。
         character(LEN=1024) filename,bandstr
         integer i,j
         integer index_i,index_ip

         character(LEN=256) noise_covmat_prefix
         namelist /noisecov/ noise_covmat_prefix
 
         !パラメタファイルからfileのprefixを読みとる。
         open(unit=100,file='input.txt')
         read(100,noisecov)
         close(100)

         !new ---------------------
         do i=1,numchan
            write(bandstr,*) int(nu_array(i))            
            filename=trim(noise_covmat_prefix)//mytrim(bandstr)//'GHz.dat'
            call Set_noise_cov_matrix(filename,covmat_noise_nu(:,:,i))
         end do
         !-------------------------
       end subroutine SetNoiseCovMatrices

       subroutine Set_noise_cov_matrix(filename,covmat_noise)
         use Utils
         real(dl), dimension(0:2*npix-1,0:2*npix-1), intent(out) :: covmat_noise
         character(LEN=1024), intent(in) :: filename
         integer(8) nreclen    

!!$         INQUIRE(IOLENGTH=nreclen) covmat_noise
!!$         open(67,file=filename,access='direct',recl=nreclen)
         open(67,file=filename,form='binary',access='sequential')
         read(67) covmat_noise
         close(67)

       end subroutine Set_noise_cov_matrix

       subroutine SetSignalCovMatrices

         character(LEN=1024) cov_scal_file,cov_tens_file
         namelist /signalcov/ cov_scal_file,cov_tens_file
         open(unit=100,file='input.txt')
         read(100,signalcov)
         close(100)

         call Set_signal_cov_matrix(cov_scal_file,covmat_signal_scalar)
         call Set_signal_cov_matrix(cov_tens_file,covmat_signal_tensor)

       end subroutine SetSignalCovMatrices

       subroutine Set_signal_cov_matrix(filename,covmat)
         use Utils
         real(dl), dimension(0:2*npix-1,0:2*npix-1), intent(out) :: covmat
         character(LEN=1024), intent(in) :: filename
         integer(8) nreclen    

!!$         INQUIRE(IOLENGTH=nreclen) covmat
!!$         open(67,file=filename,access='direct',recl=nreclen)
         open(67,file=filename,form='binary',access='sequential')
         read(67) covmat
         close(67)

       end subroutine Set_signal_cov_matrix

       
subroutine myfunc(val, n, x, grad, need_gradient, f_data)
  use Utils
  integer i
  !nlopt関係の変数---------------------
  double precision val, x(n), grad(n)
  integer n, need_gradient
  integer f_data !unused (nloptの仕様)
  !------------------------------------

  real(dl) r_est,s_est,n_est                 !ratio parameters to be estimated
  real(dl) beta_s_est,beta_d_est,gamma_d_est !foreground params
  real(dl) chisq,chisq_prior

  !mapに掛ける係数alphaとcov matrix----------------
  real(dl), allocatable, dimension(:,:) :: covmat_tot,covmat_signal,covmat_noise
  real(dl), allocatable, dimension(:,:) :: Dmatrix  !new
  real(dl), allocatable, dimension(:,:,:) :: mapout !observed map + diag_noise

  real(dl),dimension(1:numchan) :: alphas !testing
  
  allocate(covmat_tot(0:npix_tot_all-1,0:npix_tot_all-1))
  allocate(covmat_signal(0:npix_tot_all-1,0:npix_tot_all-1))
  allocate(covmat_noise(0:npix_tot_all-1,0:npix_tot_all-1)) 
  allocate(Dmatrix(0:npix_tot_all-1,0:2*npix*num_source-1))
  allocate(mapout(0:npix-1,1:nmap,1:numchan))
  !------------------------------------------------

  r_est   = x(1)
  s_est   = x(2)
  n_est   = x(3)
  beta_s_est = x(4)
  beta_d_est = x(5)
  gamma_d_est = x(6)

  call CalculateCovMatrix(r_est,s_est,n_est,covmat_tot,covmat_signal,covmat_noise)

  !mapに後から人工的なdiag noiseを加える
  mapout = map_chan + diagnoise_map_tot_all

!!$  write(*,*)"map_chan(1,2,1)=",map_chan(1,2,1)
!!$  write(*,*)"x=",x(1:6)
!!$  do i = 0,npix-1
!!$     write(100,*) mapout(i,2,6)
!!$  end do

  call CalcDmatrix(beta_s_est,beta_d_est,gamma_d_est,nu_array,Dmatrix)
  !call CalcAlphas(beta_s_est,beta_d_est,gamma_d_est,nu_array,alphas) !testing
  
  !Masked Sky Likelihood
  !val = CalculateMaskedskyLikelihood5(covmat_signal,covmat_noise,Dmatrix,mapout)
  val = CalculateMaskedskyLikelihood(covmat_signal,covmat_noise,Dmatrix,mapout)
  !val = CalculateMaskedskyLikelihood4(covmat_signal,covmat_noise,Dmatrix,alphas,mapout) 
  !val = CalculateMaskedskyLikelihood(covmat_signal,covmat_noise,Dmatrix,mapout)
  
  !val = CalculateMaskedskyLikelihood2(covmat_signal,covmat_noise,mapout_CMB,reorder_vector,chisq) 
  !full sky likelihood
  !val = CalculateFullskyLikelihood(covmat_tot,mapout_CMB(:,:,1)) 


  !priorを課す場合
!!$  chisq_prior = CalculatePrior(n,x)
!!$  val = val+chisq_prior/2
  !---

  if (need_gradient.ne.0) then
     grad(1) = 0.0
     grad(2) = 0.5 / dsqrt(x(2))
     stop 'dont come here'
  endif
  !val = dsqrt(x(2))

  deallocate(covmat_tot)
  deallocate(covmat_signal)
  deallocate(covmat_noise) 
  deallocate(mapout)
  deallocate(Dmatrix)
  
end subroutine myfunc

subroutine myconstraint(val, n, x, grad, need_gradient, d)
  integer need_gradient
  integer n
  double precision val, x(n), grad(n), d(2), a, b
  a = d(1)
  b = d(2)
  if (need_gradient.ne.0) then
     grad(1) = 3. * a * (a*x(1) + b)**2
     grad(2) = -1.0
  endif
  val = (a*x(1) + b)**3 - x(2)
end subroutine myconstraint



end module LBird


