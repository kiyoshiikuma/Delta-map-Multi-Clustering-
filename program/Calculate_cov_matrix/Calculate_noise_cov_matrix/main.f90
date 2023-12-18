program main
  !cov matrixを作る際に必要なClファイルは subroutine ReadCls で指定
  use simulatemaps
  implicit none

  integer ioerr
  !output file名
  character(LEN=1024) output_Dirname,output_filename,bandstr
  namelist /OutputDirname/ output_Dirname
  !XX covarianceが入ったファイル
  character(LEN=1024) XX_filename
  namelist /XXfile/ XX_filename

  integer, parameter :: num_band = 10
  integer, dimension(1:num_band) :: band_array
  real(dl),dimension(1:num_band) :: noise_array
  namelist /bands/ band_array,noise_array

  !ノイズの大きさ(sig_noise)
  !testingのためにcovmatをscale_param倍するためのパラメタ
  real(dl) scale_param,sig_noise
  namelist /inputparam/ scale_param

  real(dl), dimension(0:2*npix-1,0:2*npix-1) :: covmat_noise

  integer i

  !-------------------------------------------------------
  !計算結果を入れるファイルの名前をinput.txtより得る
  !-------------------------------------------------------
  open(unit=100,file='input.txt',iostat=ioerr)
  read(100,params)
  read(100,bands)          !get band_array,noise_array
  read(100,inputparam)     !scale_param
  read(100,OutputDirname)
  read(100,XXfile)
  close(100)

  !------------------------------------------------------
  !covariance matrixを計算するために必要な
  !spin spherical harmonicsのmについての和まで
  !計算したファイルを読み込む。cov matrixはこの
  !結果にclを掛けてellで和を取る
  !------------------------------------------------------
  call SetXX_with_m_summed(XX_filename)

  do i=1,num_band
     write(bandstr,*) band_array(i)
     output_filename=mytrim(output_dirname)//mytrim(bandstr)//'GHz.dat'
     sig_noise = noise_array(i)

     call Calculate_noise_cov_matrix(covmat_noise,sig_noise)

     covmat_noise = scale_param*covmat_noise

     call Output_noise_cov_matrix(output_filename,covmat_noise)
  end do

end program main

