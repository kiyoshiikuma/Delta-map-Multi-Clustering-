module Utils
  implicit none

  integer, parameter :: dl = KIND(1d0) ! selected_real_kind(15,50)
  real(dl), parameter :: pi=3.14159265358979323846264338328_dl, &
           twopi=2*pi, fourpi=4*pi
  real(dl), parameter :: rad_to_arcmin = 360.0d0*60.0d0/twopi  

  ![http://physics.nist.gov/cgi-bin/cuu/Value?khz]
  real(dl), parameter :: Kelvin_to_GHz = 20.836612_dl !fixed 2016/10/07


contains
  function InvMatrix(A) result(Ainv)
    ! Returns the inverse of a matrix calculated by finding the LU
    ! decomposition.  Depends on LAPACK.
    real(dl), dimension(:,:), intent(in) :: A
    real(dl), dimension(size(A,1),size(A,2)) :: Ainv

    real(dl), dimension(size(A,1)) :: work  ! work array for LAPACK
    integer, dimension(size(A,1)) :: ipiv   ! pivot indices
    integer :: n, info
    ! External procedures defined in LAPACK
    external DGETRF
    external DGETRI

    ! Store A in Ainv to prevent it from being overwritten by LAPACK
    Ainv = A
    n = size(A,1)

    ! DGETRF computes an LU factorization of a general M-by-N matrix A
    ! using partial pivoting with row interchanges.
    call DGETRF(n, n, Ainv, n, ipiv, info)
    if (info /= 0) then
       stop 'Matrix is numerically singular!'
    end if

    ! DGETRI computes the inverse of a matrix using the LU factorization
    ! computed by DGETRF.
    call DGETRI(n, Ainv, n, ipiv, work, n, info)
    if (info /= 0) then
       stop 'Matrix inversion failed!'
    end if

  end function InvMatrix

  function InvCovMatrix(A,logdet) result(Ainv)
    ! Returns the inverse of a matrix calculated by finding the LU
    ! decomposition.  Depends on LAPACK.
    real(dl), dimension(:,:), intent(in) :: A
    real(dl), dimension(size(A,1),size(A,2)) :: Ainv
    real(dl), optional :: logdet
    integer :: i, n, info
    ! External procedures defined in LAPACK
    external DPOTRF
    external DPOTRI

    ! Store A in Ainv to prevent it from being overwritten by LAPACK
    Ainv = A
    n = size(A,1)

    !testing (対称行列かどうか)
    !call TestSymmetry(A)

    ! DGETRF computes an LU factorization of a general M-by-N matrix A
    ! using partial pivoting with row interchanges.
    call DPOTRF('U', n, Ainv, n, info)
    if (info /= 0) then
       stop 'Matrix is numerically singular!'
    end if

    ! log(|行列式|)も同時に求めておく
    if(present(logdet)) then
       logdet=0._dl
       do i=1,n
          if(Ainv(i,i) < 0) then
             write(*,*)"Ainv(i,i) < 0 !: Ainv(i,i)=",Ainv(i,i),"i=",i
             stop
          else
             logdet=logdet+log(Ainv(i,i))
          end if
       end do
       logdet=logdet*2
    end if
    !------------------------------------

    ! DGETRI computes the inverse of a matrix using the LU factorization
    ! computed by DPOTRF ---------------
    call DPOTRI('U',n, Ainv, n, info)
    ! 三角行列の下半分を埋める
    do i=1,n
       Ainv(i+1:n,i) = Ainv(i,i+1:n)
    end do
    ! ----------------------------------
    if (info /= 0) then
       stop 'Matrix inversion failed!'
    end if
  end function InvCovMatrix

  function GetLogdetFromSymmetricMatrix(A) result(logdet)
    real(dl), dimension(:,:), intent(in) :: A
    real(dl), dimension(size(A,1),size(A,2)) :: Ainv
    real(dl) logdet
    integer :: i, n, info

    ! External procedures defined in LAPACK
    external DPOTRF
    external DPOTRI

    ! Store A in Ainv to prevent it from being overwritten by LAPACK
    Ainv = A
    n = size(A,1)

    !testing (対称行列かどうか)
    !call TestSymmetry(A)
    
    ! DGETRF computes an LU factorization of a general M-by-N matrix A
    ! using partial pivoting with row interchanges.
    call DPOTRF('U', n, Ainv, n, info)
    if (info /= 0) then
       stop 'Matrix is numerically singular!'
    end if
    
    ! log(|行列式|)を求める
    logdet=0._dl
    do i=1,n
       if(Ainv(i,i) < 0) then
          write(*,*)"Ainv(i,i) < 0 !: Ainv(i,i)=",Ainv(i,i),"i=",i
          stop
       else
          logdet=logdet+log(Ainv(i,i))
       end if
    end do
    logdet=logdet*2
    !------------------------------------

  end function GetLogdetFromSymmetricMatrix

  subroutine TestSymmetry(A)
    real(dl), dimension(:,:), intent(in) :: A
    integer i,j,n
         
    n = size(A,1)
    do i=1,n
       do j=1,n
          if(A(i,j) .ne. A(j,i)) then
             write(*,*)"i,j=",i,j
             write(*,*)"A(i,j) .ne. A(j,i)",A(i,j),A(j,i),(A(i,j)-A(j,i))/(A(i,j)+A(j,i))
             stop 'matrix is not symmetric'
          end if
       end do
    end do
    write(*,*)"matrix is symmetric"
  end subroutine TestSymmetry

  subroutine MatrixSymmetrise(A)
    real(dl), dimension(:,:), intent(inout) :: A
    real(dl), dimension(size(A,1),size(A,2)) :: Atmp
    integer i,j,n
         
    n = size(A,1)
    Atmp = A
         
    !$omp parallel do default(shared) 
    do i=1,n
       do j=1,n
          A(i,j) = (Atmp(i,j)+Atmp(j,i))/2
       end do
    end do
    
  end subroutine MatrixSymmetrise

  function mytrim(file)
    character(LEN=*), intent(in) :: file
    character(LEN=:), allocatable:: mytrim
    
    mytrim = trim(adjustl(file))
  end function mytrim

end module Utils
