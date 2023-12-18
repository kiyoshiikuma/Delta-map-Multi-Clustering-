program main
  !cov matrixを作る際に必要なClファイルは subroutine ReadCls で指定
  use simulatemaps
  implicit none

  integer ioerr

  !output file名
  character(LEN=1024) output_filename_scalar,output_filename_tensor
  namelist /OutputFilename/ output_filename_scalar,output_filename_tensor

  !XX covarianceが入ったファイル
  character(LEN=1024) XX_filename
  namelist /XXfile/ XX_filename

  real(dl), dimension(0:2*npix-1,0:2*npix-1) :: covmat_signal
  real(dl) r_in,s_in

  !-------------------------------------------------------
  !計算結果を入れるファイルの名前をinput.txtより得る
  !-------------------------------------------------------
  open(unit=100,file='input.txt',iostat=ioerr)
  read(100,OutputFilename)
  read(100,XXfile)
  close(100)

  call Allocate_arrays

  !---------------------------------------------------------------
  !camb で計算した cl を読み込む 
  !file名はLBird.f90, subroutine ReadClsで定義している。
  !---------------------------------------------------------------
  call ReadCls

  !------------------------------------------------------
  !covariance matrixを計算するために必要な
  !spin spherical harmonicsのmについての和まで
  !計算したファイルを読み込む。cov matrixはこの
  !結果にclを掛けてellで和を取る
  !------------------------------------------------------
  call SetXX_with_m_summed(XX_filename)

  r_in=0._dl
  s_in=1._dl
  call Calculate_signal_cov_matrix(covmat_signal,r_in,s_in)
  call Output_signal_cov_matrix(output_filename_scalar,covmat_signal)

  r_in=1._dl
  s_in=0._dl
  call Calculate_signal_cov_matrix(covmat_signal,r_in,s_in)
  call Output_signal_cov_matrix(output_filename_tensor,covmat_signal)

  call deAllocate_arrays
end program main

