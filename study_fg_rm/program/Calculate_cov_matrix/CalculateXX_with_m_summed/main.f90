program main
  !cov matrixを作る際に必要なClファイルは subroutine ReadCls で指定
  use simulatemaps
  implicit none

  integer ioerr

  !input.txtでoutput fileの名前を設定する。
  character(LEN=1024) output_filename
  namelist /OutputFilename/ output_filename

  !-------------------------------------------------------
  !計算結果を入れるファイルの名前をinput.txtより得る
  !-------------------------------------------------------
  open(unit=100,file='input.txt',iostat=ioerr)
  read(100,OutputFilename)
  close(100)

  call Allocate_arrays

  !------------------------------------------------------
  !covariance matrixを計算するために必要な
  !spin spherical harmonicsのmについての和まで
  !計算してファイルに書き出す。cov matrixはこの
  !結果にclを掛けてellで和を取る
  !(CalculateCovMatrix)で得られる。
  !------------------------------------------------------
  call CalculateXX_with_m_summed
  call Output_XX_with_m_summed(output_filename)
  !call SetXX_with_m_summed(output_filename)

  write(*,*)"testing output:QQ_signal_m_summed_E(5,4,5)="
  write(*,'(E25.15)') QQ_signal_m_summed_E(5,4,5)
  write(*,*)"expected(nside4)=0.254692652345917"

  call deAllocate_arrays
end program main

