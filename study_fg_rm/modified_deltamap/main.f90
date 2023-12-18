module module
  use LBird
  implicit none
  !nlopt関係の宣言--------------------------------------------
  !external myfunc, myconstraint
  !double precision d1(4), d2(4)
  integer nvar
  namelist /nlo_nvar/ nvar
  !value, step, lower & upper bounds of parameters
  double precision, dimension(:), allocatable :: x,dx,lb,ub 
  double precision minf  !the value at the minimum
  integer*8 opt
  integer ires
  include 'nlopt.f'

contains

  subroutine Allocate_local_arrays
    use LBird
    allocate(x(nvar))
    allocate(dx(nvar))
    allocate(lb(nvar))
    allocate(ub(nvar))
  end subroutine Allocate_local_arrays

  subroutine deAllocate_local_arrays
    deallocate(x)
    deallocate(dx)
    deallocate(lb)
    deallocate(ub)
  end subroutine deAllocate_local_arrays

end module module

program main
  use LBird
  use rngmod, only: rand_init,planck_rng,rand_gauss
  use module
  implicit none

  !diag noiseを足すための変数---------------------------------
  type(planck_rng) :: rng_handleQ,rng_handleU
  real(dl) noise_per_pixel

  ! CMB関係の宣言---------------------------------------------
  integer i,j,mapnum
  character(LEN=1024) maskfilename,mapdir,resultdir
  character(LEN=1024) result_filename(1:2)
  namelist /io_filenames/ maskfilename,mapdir,resultdir
  character(LEN=:), allocatable :: filename


  ! misc----------------------------
  real(dl) start,finish !cpu時間を測る
  !-------------------------------


  !reading parameters from input file ------------------------------
  open(unit=100,file='input.txt')
  read(100,nlo_nvar) !get nvar
  read(100,Params)   !get Nside, sig_noise, Tdust etc.
  read(100,numCMB)   !channelの数 (numCMBchan & numForegroundchan)

  !-----------------------------------------
  !読み込んだcmb channelの数に従って配列を定義
  !-----------------------------------------
  call SetPixelNumber !setting npix,npix_tot_all,num_source
  Td_GHz = Td_Kelvin * Kelvin_to_GHz
  call Allocate_global_arrays
  read(100,bands)!周波数の値の情報を読み込む
  read(100,io_filenames)
  close(100)

  !set result_filenames
  call SetResultFileName(resultdir,nu_array,result_filename)

  !define local variables
  call Allocate_local_arrays

  !-----------------------------------------
  !使用する mask を読み込む      
  !-----------------------------------------
  filename = mytrim(maskfilename)
  call ReadMaskFits(filename)

  !-----------------------------------------
  !global変数であるnpix_(un)masked,
  !npix_(un)masked_totをsetする。
  !-----------------------------------------
  call CountMaskedPixels

  !-----------------------------------------
  !maskされたpixel情報からcov matrixを
  !並べ替えるためのvectorを準備する
  !-----------------------------------------
  call WriteReOrderVector
  call WriteReOrderVector_for_Dmatrix
  
  !-----------------------------------------
  !covarianceは ./Calculate_(noise/signal)_cov_matrix/
  !で計算し、r_est,s_estを掛けるだけにしておく。
  !fileの場所はsubroutineの中でinput.txtから読み込む
  !-----------------------------------------
  call SetNoiseCovMatrices 
  call SetSignalCovMatrices

  !-------------------------------------
  !mapに人工的に加えるdiag noiseの初期化
  !-------------------------------------
  call rand_init(rng_handleQ,1427)
  call rand_init(rng_handleU,2349)
  noise_per_pixel = NoisePerPix(sig_noise,npix)

  !open output files
  open(unit=60,file=result_filename(1))
  open(unit=61,file=result_filename(2))

  !解析するCMB mapについてのループ
  do mapnum=1,100
     write(*,*)"analysing map: mupnum=",mapnum

     !new
     do j=1,numchan
        call read_map_at_XXGHz(mapdir,nu_array(j),mapnum,map_chan(0:,1:,j))
     end do

     
     !mapoutに加えるための人工的なdiag noiseの用意
     !new----------------------------------------
     do j=1,numchan
        do i=0,npix-1
           diagnoise_map_tot_all(i,1,j) = 0._dl !temperature
           diagnoise_map_tot_all(i,2,j) = noise_per_pixel*rand_gauss(rng_handleQ)
           diagnoise_map_tot_all(i,3,j) = noise_per_pixel*rand_gauss(rng_handleU)

           !testing(CMB channelにだけdiag noiseをのせる)---------
!!$           if(j>1) then
!!$              diagnoise_map_tot_all(i,2,j) = 0._dl
!!$              diagnoise_map_tot_all(i,3,j) = 0._dl
!!$           end if
           !-----------------------------------------------------
        end do
     end do
     !-------------------------------------------

     !set algorithm and dimensionality -----------------
     ! Non-derivative algorithm are:
     ! LN_SBPLX,LN_NELDERMEAD(あまり良くない),LN_PRAXIS
     ! LN_NEWUOA, LN_BOBYQA, LN_COBYLA(時間かかる)
     ! LN_PRAXISはdimensionality=1だと失敗。
     !--------------------------------------------------     
     call nlo_create(opt, NLOPT_LN_BOBYQA, nvar)

     !精度の設定
     call nlo_set_xtol_rel(ires,opt,1.D-5) 

     !初期条件、上限、下限、最初のステップ幅
     call SetInitialValuesAndBounds(ires,opt,lb,ub,x,dx,nvar)
     call nlo_set_min_objective(ires, opt, myfunc, 0)

     !解析開始
     call cpu_time(start)
     call nlo_optimize(ires,opt,x,minf)
     if (ires.lt.0) stop 'nlopt failed!'
     call cpu_time(finish)
     !結果の記録
     !write(*,*)"time score:",finish-start
     call OutputResults(mapnum,nvar,x,minf)

     call nlo_destroy(opt)
  end do

  close(60)
  close(61)

  call deAllocate_global_arrays
  call deAllocate_local_arrays


end program main

subroutine OutputResults(mapnum,n,x,minf)
  use LBird
  use fitstools
  implicit none
  
  integer, intent(in) :: mapnum
  integer, intent(in) :: n
  real(dl), dimension(1:n), intent(in) :: x
  real(dl), intent(in) :: minf

  real(dl) r_est,s_est,n_est
  real(dl) val
  real(dl), allocatable, dimension(:,:) :: Dmatrix  !new
  real(dl), allocatable, dimension(:,:) :: covmat_tot,covmat_signal,covmat_noise
  real(dl), allocatable, dimension(:,:,:) :: mapout
  real(dl) beta_s_est,beta_d_est,gamma_d_est
  integer nout


  allocate(covmat_tot(0:npix_tot_all-1,0:npix_tot_all-1))
  allocate(covmat_signal(0:npix_tot_all-1,0:npix_tot_all-1))
  allocate(covmat_noise(0:npix_tot_all-1,0:npix_tot_all-1)) 
  allocate(mapout(0:npix-1,1:nmap,1:numchan))
  allocate(Dmatrix(0:npix_tot_all-1,0:2*npix*num_source-1))
  
  r_est   = x(1)
  s_est   = x(2)
  n_est   = x(3)
  beta_s_est  = x(4)
  beta_d_est  = x(5)
  gamma_d_est = x(6)

  call CalculateCovMatrix(r_est,s_est,n_est,covmat_tot,covmat_signal,covmat_noise)
  mapout = map_chan + diagnoise_map_tot_all

!!!! Dmatrixを計算してlikelihoodを計算
    call CalcDmatrix(beta_s_est,beta_d_est,gamma_d_est,nu_array,Dmatrix)
    !!!! calculateMaskedskyLikelihood3を作成すべし
  !val = CalculateMaskedskyLikelihood3(covmat_signal,covmat_noise,mapout,reorder_vector,chisq) 
    !val = CalculateMaskedskyLikelihood4(covmat_signal,covmat_noise,Dmatrix,alphas,mapout)
    !val = CalculateMaskedskyLikelihood(covmat_signal,covmat_noise,Dmatrix,mapout) 
  val = CalculateMaskedskyLikelihood5(covmat_signal,covmat_noise,Dmatrix,mapout)
  
  write(*,'("  -lnLike =",2f25.15)') val,minf
  write(*,'("        r =",f25.15)') r_est!x(1)
  write(*,'("   beta_s =",f25.15)') beta_s_est!x(2)
  write(*,'("   beta_d =",f25.15)') beta_d_est!x(3)
  write(*,'("  gamma_d =",f25.15)') gamma_d_est/Kelvin_to_GHz
  write(*,'("        s =",f15.5)') x(2)
  write(*,'("  c_noise =",f15.5)') x(3)


  nout = 2+n
  write(60,'(<nout>E25.15)') minf,r_est,s_est,n_est,beta_s_est,beta_d_est,gamma_d_est/Kelvin_to_GHz,val


  deallocate(covmat_tot)
  deallocate(covmat_signal)
  deallocate(covmat_noise) 
  deallocate(mapout)
  deallocate(Dmatrix)
end subroutine OutputResults

subroutine read_map_at_XXGHz(mapdir,XX,mapnum,map)
  !解析に用いるmapを読み込む。ディレクトリ構造は以下のようにしておく
  ! mapdir/XXGHz/mapnum.fits (for mapnum = 1,2,...,100)
  use LBird

  character(LEN=*), intent(in) :: mapdir
  real(dl), intent(in) :: XX !GHz
  integer, intent(in) :: mapnum
  character(LEN=:), allocatable :: file
  character(LEN=64) bandstr,mapnumber
  real(dl), dimension(0:npix-1,1:nmap), intent(out) :: map

  write(mapnumber,*) mapnum
  write(bandstr,*) int(XX)
  file = mytrim(mapdir)//mytrim(bandstr)//'GHz/'//mytrim(mapnumber)//'.fits'
  call read_map(file,map)
  
end subroutine read_map_at_XXGHz

subroutine read_map(mapfilename,mapout)
  use LBird
  use fitstools
  implicit none
  character(LEN=*), intent(in) :: mapfilename
  logical anynull
  real(dl) nullval
  real(dl), dimension(0:npix-1,1:nmap), intent(out) :: mapout

  call read_bintab(mapfilename,mapout,npix,nmap,nullval,anynull)         

end subroutine read_map

subroutine SetInitialValuesAndBounds(ires,opt,lb,ub,x,dx,nvar)
  use LBird
  implicit none
  integer i,nvar
  integer ires
  integer*8 opt
  double precision x(nvar),dx(nvar)  !value and step of parameters
  double precision lb(nvar),ub(nvar) !lower & uapper bounds

  !get default -(+)infty lower (upper) values first
  call nlo_get_lower_bounds(ires, opt, lb)
  call nlo_get_upper_bounds(ires, opt, ub)

  !initial guess values and initial step
  !for [r,s,noise]
  x(1:3) = [1.0d-1,1.0d0,1.0d0]
  x(4) = -3.2d0 !beta_s
  x(5) = 2.65d0 !beta_d
  x(6) = Td_GHz !gamma_d

  !for [r,s,noise]
  dx(1:3) = [1.0d-3,0.1d0,0.1d0] 
  dx(4) = 0.3d0     !beta_s
  dx(5) = 0.1d0    !beta_d
  dx(6) = Td_GHz/10.0 !gamma_d

  !lower bounds
  lb(1:3) = 0._dl ![r,s,noise]
  lb(4) = -4._dl  !beta_s
  lb(5) =  1._dl  !beta_d
  lb(6) =  5._dl*Kelvin_to_GHz  !gamma_d

  !upper bounds
  ub(4) = -2._dl  !beta_s
  ub(5) =  3._dl  !beta_d
  ub(6) =  40._dl*Kelvin_to_GHz  !gamma_d

  !scalar normalization(s)とnoise normalization(n)は固定
  lb(2) = 1._dl
  ub(2) = lb(2) !x(2)=sは動かない

  lb(3) = 1._dl
  ub(3) = lb(3) !x(3)=noise normalizationは動かない


  if(algorithm .eq. "PowLawDust") then
     lb(4)=-3._dl
     ub(4)=lb(4)

     lb(6)=Td_GHz !19.15*Kelvin_to_GHz
     ub(6)=lb(6)
  end if

    if(algorithm .eq. "PowLawSync") then
     lb(5)=1._dl
     ub(5)=lb(5)

     lb(6)=Td_GHz !19.15*Kelvin_to_GHz
     ub(6)=lb(6)
  end if

  if(algorithm .eq. "PowLawDustSync") then
     lb(6)=19.15*Kelvin_to_GHz
     ub(6)=lb(6)
  end if

  if(algorithm .eq. "psm162_deltamap") then
     lb(6)=19.15*Kelvin_to_GHz
     ub(6)=lb(6)
  end if

  call nlo_set_lower_bounds(ires, opt, lb)
  call nlo_set_upper_bounds(ires, opt, ub)
  call nlo_set_initial_step(ires, opt, dx)
end subroutine SetInitialValuesAndBounds


subroutine SetResultFileName(dirname,nu,filename)
  use LBird
    implicit none
      real(dl), dimension(1:numchan),intent(in) :: nu
      character(LEN=1024), intent(out) ::  filename(1:2)
      character(LEN=1024), intent(in) ::  dirname
      character(LEN=1024) char1,char2
      integer i

      char2 = ''
      do i = 1,numchan
         write(char1,*) int(nu(i))
         char2 = mytrim(char2)//'_'//mytrim(char1)
      end do

      filename(1) = mytrim(dirname)//mytrim(algorithm)//'_'//mytrim(char2)//'.dat'
      filename(2) = mytrim(dirname)//'alphavalues_'//mytrim(char2)//'.dat'

end subroutine SetResultFileName

subroutine SetPixelNumber
  use LBird
  implicit none
  npix=12*Nside**2
  npix_tot_all = numchan*2*npix

 !default deltamapではnvar_foreは
 !(beta_s,beta_d,T_d)でnvar_fore=3
 !sourceとしてはbar(beta_s)とd(beta_s)でsynchが2つ
 !bar(beta_d,T_d)とd(beta_d)&d(T_d)でnum_source=5
  if(algorithm .eq. "modified_delta_map")   num_source =5 
  if(algorithm .eq. "PowLawDust")   num_source =1 !spacially homogeneous spectral index
  if(algorithm .eq. "PowLawSync")   num_source =1 !spacially homogeneous spectral index
  if(algorithm .eq. "PowLawDustSync")   num_source =2 !spacially homogeneous spectral indices
  if(algorithm .eq. "psm162_deltamap")   num_source =4 !spacially varying spectral indices
end subroutine SetPixelNumber
