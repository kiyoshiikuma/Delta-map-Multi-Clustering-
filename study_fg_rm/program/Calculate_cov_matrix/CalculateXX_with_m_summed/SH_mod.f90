! Routines to compute the theta part of the spin harmonics and their equal m
! integrals from cos(theta) = x to 1. 
! The integral returned for l=l2 is therefore always *negative*. 
!
! Antony Lewis Oct 2002 (http://cosmologist.info)
! Oct 02: fixed bug in SH_initialize
! Mar08 added DH fix for overflow at high l
! ** this version has a thread safe versions of the GetIntegral routines **
! Integrals for separate m can be computed in parallel
! Ref: astro-ph/0106536

! Note that overlap integrals become very small for m >> l sin(theta) = l sqrt(1-x^2)
! and may give NaN results in extremes (e.g. x=0.9995, l=2500, m> 1000). Only generate
! them where they are significantly non-zero

!We use {}_s \lambda_{lm}, where {}_s Y_{lm} = {}_s \lambda_{lm} e^{im\phi}


module SpinHarmonics

  implicit none
  private

  integer  :: lmax
  integer, parameter :: dl = KIND(1d0) ! selected_real_kind(15,50)
  integer, parameter :: i8 = Selected_Int_Kind(8)
  integer  :: spin  ! only checked for spin=0,2
  integer  :: sp    ! abs(spin)
  real(dl), dimension (:,:), allocatable :: Cslm, Harmonic
  real(dl), dimension (:,:,:), allocatable :: Integral

  real(dl), parameter :: pi=3.14159265358979323846264338328_dl, &
      twopi=2*pi, fourpi=4*pi

  Type SH_mValues
   real(dl), dimension (:,:), pointer:: mIntegralP, mIntegralM
   integer :: mVal
  end  Type SH_mValues


!If spin=0, W, otherwise two spin two matrixes W_+ and W_-
  Type SH_Couplings  
   real(dl), dimension (:,:), pointer :: W, WP, WM
  end Type SH_Couplings

  integer :: lastlmax, lastspin
  real(dl) :: xVal, lastx

  ! The first four routines compute individual values, the remaineder
  ! are for computing arrays

  public SH_lambdaslm, SH_GetAmsss, SH_GetDiagInt, SH_CalcOffDiagInt, &
      SH_GetDiagIntArray,SH_initialize,  SH_FreeMem, SH_GetHarmonics, &
      SH_Lambda, SH_GetInt, SH_GetIntegrals, SH_GetmIntegrals, SH_mValues, &
      SH_GetCouplingMatrices, SH_FreeCouplingMatrices, SH_Couplings,SH_InitHarmonics,&
      SH_dlmn, SH_LambdaArray,SH_NullCouplingMatrices

contains

  function SH_dlmn(l,m,n,x)
    integer, intent(in) :: l,m,n
    real(dl), intent(in) :: x
    real(dl) SH_dlmn
    
    SH_dlmn = sqrt(fourpi/(2*l+1))*SH_lambdaslm(n,l,-m,x)
    if (mod(m,2)/=0) SH_dlmn = -SH_dlmn
  end  function SH_dlmn

  function SH_lambdaslm(ss,n,mm,x)
    !Get one-off lambda at any s (ss),l (n),m (mm), x=cos(theta)
    !This is not designed to return accurate answers for |x| very
    ! close to 1 or 0
    real(dl) SH_lambdaslm
    integer, intent(IN) :: ss,n,mm
    real(dl), intent(IN):: x
    real(dl) Pn,Pmm,pmmp1,fact
    integer i,l,m,s

    s=ss
    m=mm

    if (abs(m)>n .or. n<abs(s)) then
      SH_lambdaslm=0._dl
    else
      Pmm=1._dl
      if (abs(m)<abs(s)) then
        s=mm
        m=ss
        if (IAND(m+s,1) == 1) Pmm=-Pmm
      end if

      If (m<0) then
        m=-m
        s=-s
        if (IAND(m+s,1) == 1) Pmm=-Pmm
      end if

      fact=1._dl

      do i=m+abs(s)+1,2*m
        fact =fact*i/real(4*(2*m+1-i),dl)
      end do
      fact =sqrt(fact*(2*m+1)/fourpi)
      if (m/=-s) Pmm=Pmm*(1-x)**((m+s)/2._dl)
      if (m/=s) Pmm=Pmm*(1+x)**((m-s)/2._dl)
      Pmm=Pmm/real(2**abs(s),dl)*fact
      if (IAND(m,1) ==1) Pmm=-Pmm

      if (n==m) then
        SH_lambdaslm=Pmm 
      else
        pmmp1=(s+x*(m+1))/sqrt( ((m+1)**2-s**2)/real(2*m+3,dl) )*Pmm
        if (n==m+1) then
          SH_lambdaslm=pmmp1
        else
          do l=m+2,n
            Pn= ((l*(l-1)*x + s*m)*pmmp1 -l*Pmm*sqrt(((l-1)**2-m**2)*real((l-1)**2-s**2,dl)/(2*l-3)/(2*l-1))) &
                *sqrt((4*l*l-1)/real(l**2-m**2,dl)/(l**2-s**2))/(l-1)
            Pmm=pmmp1
            pmmp1=Pn
          end do
          SH_lambdaslm=Pn

        end if
      end if
    end if

  end function SH_lambdaslm

  

  subroutine SH_LambdaArray(ss,mm,x,arr)
  !Get array of lambdas from l_min up to lmax
    implicit none
    integer, intent(IN)  :: mm,ss
    real(dl), intent(IN) :: x
    real(dl) :: arr(0:lmax)
    real(dl) Pn,pmmp1, Pmm,fact
    integer i,l,m,s

    if (.not. allocated(Cslm)) stop 'Must initialize SpinHarmonics'

    s=ss
    m=mm
    Pmm=1._dl

    if (abs(m)<abs(s)) then
      s=mm
      m=ss
      if (IAND(m+s,1) ==1) Pmm=-Pmm
    end if

    If (m<0) then
      m=-m
      s=-s
      if (IAND(m+s,1) ==1) Pmm=-Pmm
    end if
    if (abs(s) /= sp) stop 'Array routines assume spin used to initialize'

    fact=1.0
    do i=m+abs(s)+1,2*m
      fact =fact*i/real(4*(2*m+1-i),dl)
    end do
    fact =sqrt(fact*(2*m+1)/fourpi)
    if (m/=s) Pmm=Pmm*(1+x)**((m-s)/2._dl)
    if (m/=-s) Pmm=Pmm*(1-x)**((m+s)/2._dl)
    Pmm=Pmm/(2**abs(s))*fact
    if (IAND(m,1) ==1) Pmm=-Pmm  

    arr(m)=Pmm
    if (m <= lmax) then 
      pmmp1=(s+x*(m+1))/sqrt( ((m+1)**2-s**2)/real(2*m+3,dl))*Pmm
      arr(m+1)=pmmp1

      do l=m+2,lmax
        Pn= ((x + s*m/real(l*(l-1),dl))*pmmp1 -Pmm/Cslm(l-1,m))*Cslm(l,m)
        Pmm=pmmp1
        pmmp1=Pn
        arr(l)=Pn
      end do
    end if

  end subroutine SH_LambdaArray

  function SH_CalcOffDiagInt(s,m,l,l2,x)
    real(dl)  SH_CalcOffDiagInt
    integer, intent(IN) :: s,m,l,l2
    real(dl), intent(IN):: x
    real(dl) tmp
    integer mm

    if (allocated(Harmonic) .and. s==lastspin .and. x==lastx .and. lastlmax+1>=max(l,l2) ) then
      mm=abs(m)
      if (sp==0) then
        tmp = x*SH_Lambda(l,m)*SH_Lambda(l2,m) 
      else
        tmp =(x-real(spin*m,dl)/(l*l2))*SH_Lambda(l,m)*SH_Lambda(l2,m) 
      end if
      if (l2 > mm .and. l2 > sp) tmp=tmp + (2*l2+1)/Cslm(l2,mm)*SH_Lambda(l,m)*SH_Lambda(l2-1,m)/(l-l2) 
      if (l> mm .and. l > sp) tmp = tmp- SH_Lambda(l-1,m)*SH_Lambda(l2,m)*(2*l+1)/Cslm(l,mm)/(l-l2)
      SH_CalcOffDiagInt=tmp*twopi/(l+l2+1)

    else
      SH_CalcOffDiagInt= twopi*sqrt((1-x)*(1+x))/(l2*(l2+1)-l*(l+1))* &
          ( sqrt(real(l2*(l2+1)-m*(m+1),dl))*SH_lambdaslm(s,l,m,x)*SH_lambdaslm(s,l2,m+1,x) - &
          sqrt(real(l*(l+1)-m*(m+1),dl))*SH_lambdaslm(s,l,m+1,x)*SH_lambdaslm(s,l2,m,x) )
    end if

  end  function SH_CalcOffDiagInt

  recursive function SH_GetAmsss(s,m,x) result(R)
    implicit none
    real(dl) R
    integer, intent(IN) :: s,m
    real(dl), intent(IN):: x

    if (s==m) then
      R = ((x-1)/2)**(2*s+1)
    else if (s==-m) then
      R = (((x+1)/2)**(2*s+1)-1._dl)
    else
      R = SH_GetAmsss(s-1,m,x) + x*twopi/(2*s+1) * &
          SH_lambdaslm(m,s,s,x)**2 + m/sqrt(real((2*s+1)*(s**2-m**2),dl)) * &
          SH_CalcOffDiagInt(m,s-1,s,s-1,x)
    end if

  end function SH_GetAmsss

  recursive function SH_GetDiagInt(m,l,x) result(R)
    !For use when we don't want to store all m values. |m|=l only.
    real(dl) R
    integer, intent(IN)  :: m, l
    real(dl), intent(IN) :: x
    integer mm
    real(dl) tmp

    mm=abs(m)

    if (l<sp .or. mm/=l) then
      stop 'SH_GetDiagInt with strange parameters'
    else if (l==sp .and. m==sp) then 
      R = ((x-1)/2)**(2*sp+1)
    else if (l==sp .and. m==-sp) then
      R =  (((x+1)/2)**(2*sp+1)-1._dl)   
    else if (l==m) then
      tmp = SH_GetDiagInt(m-1,l-1,x) + x*twopi/(2*m+1)* &
          SH_lambdaslm(spin,m,m,x)**2 
      if (sp /=0) tmp=tmp + &
          spin/sqrt((2*m+1)*real(m**2-spin**2,dl)) * &
          SH_CalcOffDiagInt(spin,m-1,m,m-1,x)
      R = tmp
    else if (l==-m) then
      tmp =  SH_GetDiagInt(-mm+1,l-1,x) + (x*twopi* &
          SH_lambdaslm(spin,mm,-mm,x)**2)/(2*mm+1) 
      if (spin /=0) tmp=tmp - &
          spin/sqrt((2*mm+1)*real(m**2-spin**2,dl)) * &
          SH_CalcOffDiagInt(spin,-mm+1,mm-1,mm,x)
      R = tmp
    else if (mm < sp) then
      R = SH_GetAmsss(sp,m,x)
    end if

  end function SH_GetDiagInt

  subroutine SH_GetDiagIntArray(m,arr,x)
    integer, intent(IN)   :: m
    real(dl), intent(OUT) :: arr(0:lmax)
    real(dl), intent(IN)  :: x
    integer mm,l
    real(dl)tmp

    if (.not. allocated(Cslm)) stop 'Must initialize SpinHarmonics'

    call SH_GetHarmonics(x)

    mm=abs(m)

    if (mm>=sp) then
      arr(mm) = SH_GetDiagInt(m,mm,x)
    else
      arr(sp) = SH_GetAmsss(sp,m,x)
    end if

    do l=max(mm+1,sp+1),lmax
      tmp = arr(l-1)  + & 
          Cslm(l,mm)/Cslm(l+1,mm)*SH_CalcOffDiagInt(spin,m,l+1,l-1,x) 
      if ((l-2 >= sp) .and. (l-2 >=mm)) tmp  = tmp- &
          Cslm(l,mm)/Cslm(l-1,mm)*SH_CalcOffDiagInt(spin,m,l,l-2,x)
      if ( sp > 0) tmp  = tmp &
          + 2*spin*m/real(l*l-1,dl)/l*Cslm(l,mm)*SH_CalcOffDiagInt(spin,m,l,l-1,x)
      arr(l)=tmp
    end do
  end subroutine SH_GetDiagIntArray


  subroutine SH_initialize(almax, aspin)
    integer, intent(IN) :: almax, aspin
    integer l,m
    real(dl) l2

    if (aspin /=0 .and. aspin /=2) print *,'Warning: not tested for spin /= 0, 2'
    spin=aspin
    sp=abs(aspin)
    lmax=almax
    call SH_FreeMem

    allocate(Cslm(0:lmax+1,0:lmax+1))

    do m=0, lmax+1
      do l=max(sp,m)+1,lmax+1
        l2=l*l
        Cslm(l,m)=sqrt(real(l2*(4*l2-1),dl)/((l2-m*m)*(l2-spin**2)))
      end do
    end do

  end subroutine SH_initialize


  subroutine SH_GetHarmonics(x)
    real(dl), intent(IN) :: x
    integer m,l

    if (allocated(Harmonic)) then
      if (x==lastx .and. lmax==lastlmax .and. spin==lastspin) return
      deallocate(Harmonic)
    end if
    allocate(Harmonic(sp:(lmax+1),-(lmax+1):(lmax+1))) 
    lastlmax=lmax
    lastspin=spin
    lastx=x

    do m=-lmax-1,-sp
      call GetLambdaArray(m,x)
    end do

    do m=sp,lmax+1
      call GetLambdaArray(m,x)
    end do

    do m=-sp+1,sp-1
      do l=sp,lmax+1  
        Harmonic(l,m) = SH_lambdaslm(spin,l,m,x)
      end do
    end do

  end  subroutine SH_GetHarmonics


  subroutine GetLambdaArray(mm,x)
    implicit none
    integer, intent(IN)  :: mm
    real(dl), intent(IN) :: x
    real(dl) Pn,pmmp1, Pmm,fact
    integer i,l,m,s
    integer :: factexp


    s=spin
    m=mm
    Pmm=1._dl

    if (abs(m)<abs(s)) then
      s=mm
      m=spin
      if (IAND(m+s,1) ==1) Pmm=-Pmm
    end if

    If (m<0) then
      m=-m
      s=-s
      if (IAND(m+s,1) ==1) Pmm=-Pmm
    end if

    fact=1._dl
    factexp=0
    do i=m+abs(s)+1,2*m
      fact =fact*i/real(4*(2*m+1-i),dl)
!Following overlow fix by Duncan Hanson Mar08
      factexp = factexp + exponent(fact)
      fact = set_exponent(fact, 0)
    end do
    fact = set_exponent(fact, factexp)

    fact =sqrt(fact*(2*m+1)/fourpi)
    if (m/=s) Pmm=Pmm*(1+x)**((m-s)/2._dl)
    if (m/=-s) Pmm=Pmm*(1-x)**((m+s)/2._dl)
    Pmm=Pmm/(2**abs(s))*fact
    if (IAND(m,1) ==1) Pmm=-Pmm  

    Harmonic(m,mm)=Pmm
    if (m <= lmax) then 
      pmmp1=(s+x*(m+1))/sqrt( ((m+1)**2-s**2)/real(2*m+3,dl))*Pmm
      Harmonic(m+1,mm)=pmmp1

      do l=m+2,lmax+1
        Pn= ((x + s*m/real(l*(l-1),dl))*pmmp1 -Pmm/Cslm(l-1,m))*Cslm(l,m)
        Pmm=pmmp1
        pmmp1=Pn
        Harmonic(l,mm)=Pn
      end do
    end if

  end subroutine GetLambdaArray

  function SH_Lambda(l,m)
    real(dl) SH_Lambda
    integer, intent(IN) :: l,m

    if (abs(m)>l .or. l < sp) then
      SH_Lambda = 0._dl
    else
      SH_Lambda = Harmonic(l,m)
    end if

  end function SH_Lambda


  function SH_GetInt(I,l,l2,m)
    implicit none
    Type(SH_mValues) I
    real(dl) SH_GetInt
    integer, intent(IN) :: m,l,l2

    if (abs(m) > min(l,l2) .or. (min(l,l2)<sp)) then

      SH_GetInt= 0._dl

    else                
      if (max(l,l2)> lmax) then
        if (l==l2) stop 'l=l2 too big'
        SH_GetInt = SH_CalcOffDiagInt(spin,m,l,l2,xVal)
      else
        if (allocated(Integral)) then
          if (l<=l2) then            
            SH_GetInt=Integral(m,l,l2)
          else
            SH_GetInt=Integral(m,l2,l)
          end if
        else
          if (abs(m)/=I%mVal) stop 'wrong m :'
          if (m<0) then
            SH_GetInt=I%mIntegralM(min(l2,l),max(l,l2))
          else
            SH_GetInt=I%mIntegralP(min(l2,l),max(l,l2))
          end if
        end if
      end if
    end if

  end function SH_GetInt


  subroutine GetOffDiags(Arr,m,x)
    integer, intent(IN)   :: m
    real(dl), intent(IN)  :: x
    real(dl), intent(OUT) :: Arr(sp:lmax,sp:lmax)
    integer mm,l,l2
    real(dl) tmp

    mm=abs(m)
    do l=max(sp,mm),lmax
      do l2=l,lmax
        if (l/=l2) then
          if (sp==0) then
            tmp = x*SH_Lambda(l,m)*SH_Lambda(l2,m) 
          else
            tmp =(x-real(spin*m,dl)/(l*l2))*SH_Lambda(l,m)*SH_Lambda(l2,m) 
          end if
          if (l2 > mm .and. l2 > sp) tmp=tmp + (2*l2+1)/Cslm(l2,mm)*SH_Lambda(l,m)*SH_Lambda(l2-1,m)/(l-l2) 
          if (l> mm .and. l > sp) tmp = tmp- SH_Lambda(l-1,m)*SH_Lambda(l2,m)*(2*l+1)/Cslm(l,mm)/(l-l2)
          Arr(l,l2) =tmp*twopi/(l+l2+1)
        end if
      end do
    end do

  end  subroutine GetOffDiags

  subroutine SH_FreeMem

    if (allocated(Cslm)) deallocate(Cslm)
    if (allocated(Harmonic)) deallocate(Harmonic)
    if (allocated(Integral)) deallocate(Integral)

  end  subroutine SH_FreeMem


  subroutine SH_FreemIntegrals(I)
    Type(SH_mValues) I
   
    if (associated(I%mIntegralP)) deallocate(I%mIntegralP)
    if (associated(I%mIntegralM)) deallocate(I%mIntegralM)
    nullify(I%mIntegralM,I%mIntegralP)

  end subroutine SH_FreemIntegrals

  subroutine SH_InitHarmonics(x)
    real(dl), intent(IN) :: x

    xVal = x
    call SH_GetHarmonics(x)
     
  end subroutine SH_InitHarmonics

  subroutine SH_GetmIntegrals(I,m,x)
    implicit none
    Type(SH_mValues) I
    integer, intent(IN)  :: m
    real(dl), intent(IN) :: x

    integer l,mm,minl

    if (x /= xVal) stop 'Must call SH_InitHarmonics before SH_GetmIntegrals with correct x'
    mm=abs(m)
    I%mVal=mm
    minl=max(sp,mm)

    if (associated(I%mIntegralP)) deallocate(I%mIntegralP)
    if (associated(I%mIntegralM)) deallocate(I%mIntegralM)
    Allocate(I%mIntegralP(sp:lmax,sp:lmax))
    Allocate(I%mIntegralM(sp:lmax,sp:lmax))

    call GetOffDiags(I%mIntegralP(:,:),mm,x)
    call GetOffDiags(I%mIntegralM(:,:),-mm,x)

    if (mm>=sp) then
      I%mIntegralP(mm,mm) = SH_GetDiagInt(mm,mm,x)
      I%mIntegralM(mm,mm) = SH_GetDiagInt(-mm,mm,x)
    else
      I%mIntegralP(sp,sp) = SH_GetAmsss(sp,mm,x)
      I%mIntegralM(sp,sp) = SH_GetAmsss(sp,-mm,x)
    end if


    do l=max(mm+1,sp+1),lmax

      I%mIntegralP(l,l) = I%mIntegralP(l-1,l-1) + & 
          Cslm(l,mm)/Cslm(l+1,mm)*SH_GetInt(I,l+1,l-1,mm) 
      if ((l-2 >= sp) .and. (l-2 >=mm))  I%mIntegralP(l,l) = I%mIntegralP(l,l)- &
          Cslm(l,mm)/Cslm(l-1,mm)*SH_GetInt(I,l,l-2,mm)
      if ( sp > 0) I%mIntegralP(l,l) = I%mIntegralP(l,l) + 2*spin*mm/real(l*l-1,dl)/l*Cslm(l,mm)*SH_GetInt(I,l,l-1,mm)

      I%mIntegralM(l,l) = I%mIntegralM(l-1,l-1) + & 
          Cslm(l,mm)/Cslm(l+1,mm)*SH_GetInt(I,l+1,l-1,-mm) 
      if ((l-2 >= sp) .and. (l-2 >=mm))  I%mIntegralM(l,l) = I%mIntegralM(l,l)- &
          Cslm(l,mm)/Cslm(l-1,mm)*SH_GetInt(I,l,l-2,-mm)
      if ( sp > 0) I%mIntegralM(l,l) = I%mIntegralM(l,l) - 2*spin*mm/real(l*l-1,dl)/l*Cslm(l,mm)*SH_GetInt(I,l,l-1,-mm)

    end do

  end subroutine SH_GetmIntegrals

  subroutine SH_GetIntegrals(x)
    implicit none
    real(dl), intent(IN) :: x
    integer l,m,s,mm
    Type(SH_mValues) I

    I%mVal = -1 !I is dummy here

    xVal = x

    call SH_GetHarmonics(x)


    s=sp
    if (allocated(Integral)) deallocate(Integral)
    Allocate(Integral(-lmax:lmax,sp:lmax,sp:lmax))

    do m=-lmax,lmax
      call GetOffDiags(Integral(m,:,:),m,x)
    end do

    Integral(s,s,s)  = ((x-1)/2)**(2*s+1)
    Integral(-s,s,s) = (((x+1)/2)**(2*s+1)-1._dl)


    do m=s+1, lmax
      Integral(m,m,m) = Integral(m-1,m-1,m-1) + x*twopi/(2*m+1)*SH_Lambda(m,m)**2 
      if (spin /=0) Integral(m,m,m)=Integral(m,m,m) + &
          spin/sqrt((2*m+1)*real(m*m-spin*spin,dl)) *SH_GetInt(I,m-1,m,m-1)

    end do
    do m=s+1, lmax
      Integral(-m,m,m) = Integral(-m+1,m-1,m-1) + (x*twopi*SH_Lambda(m,-m)**2)/(2*m+1) 
      if (spin /=0) Integral(-m,m,m)=Integral(-m,m,m)-&
          spin/sqrt((2*m+1)*real(m*m-spin*spin,dl)) *SH_GetInt(I,m-1,m,-m+1)

    end do

    do m=-s+1,s-1
      Integral(m,s,s) = SH_GetAmsss(s,m,x)
    end do


    do m=-lmax,lmax
      mm=abs(m)  
      do l=max(mm+1,s+1),lmax
        if (mm /= l) then
          Integral(m,l,l) = Integral(m,l-1,l-1) + & 
              Cslm(l,mm)/Cslm(l+1,mm)*SH_GetInt(I,l+1,l-1,m) 
          if ((l-2 >= s) .and. (l-2 >=mm))  Integral(m,l,l) = Integral(m,l,l)- &
              Cslm(l,mm)/Cslm(l-1,mm)*SH_GetInt(I,l,l-2,m)
          if ( s > 0) Integral(m,l,l) = Integral(m,l,l) + 2*spin*m/real(l*l-1,dl)/l*Cslm(l,mm)*SH_GetInt(I,l,l-1,m)

        end if
      end do

    end do



  end  subroutine SH_GetIntegrals


  subroutine SH_InitmValues(I)
   Type(SH_mValues) I

   nullify(I%mIntegralM,I%mIntegralP)
    
  end subroutine SH_InitmValues


  subroutine SH_GetCouplingMatrices(W,twocaps,m,x)
   logical, intent(in) :: twocaps
   integer, intent(in) :: m
   real(dl), intent(in) :: x
   integer l,l2,n, lmin
   Type(SH_mValues) I
   Type(SH_Couplings) W


  call SH_InitmValues(I)

  lmin=max(abs(m),sp)
  n=lmax-lmin+1
  call SH_GetmIntegrals(I,m,x)
  

 if (sp==0) then
   if (associated(W%W)) deallocate(W%W)
   allocate(W%W(n,n))
     do l=lmin,lmax
      do l2=lmin,lmax
      !DONT FORGET SIGN OF W%WP!
      if (twocaps) then
        if (IAND(l+l2,1) == 1) then
         W%W(l-lmin+1,l2-lmin+1)=0
        else
         W%W(l-lmin+1,l2-lmin+1)=-2*SH_GetInt(I,l,l2,m)
        end if
      else
        W%W(l-lmin+1,l2-lmin+1)=-SH_GetInt(I,l,l2,m)
      end if
      end do
     enddo

 else if (sp==2) then

     if (associated(W%WP)) deallocate(W%WP,W%WM)
     allocate(W%WP(n,n))
     allocate(W%WM(n,n))
 
     do l=lmin,lmax
      do l2=lmin,lmax
      !DONT FORGET SIGN OF W%WP!
      if (twocaps) then
        if (IAND(l+l2,1) == 1) then
         W%WP(l-lmin+1,l2-lmin+1)=0
         W%WM(l-lmin+1,l2-lmin+1)=-(SH_GetInt(I,l,l2,m)-SH_GetInt(I,l,l2,-m))
        else
         W%WP(l-lmin+1,l2-lmin+1)=-(SH_GetInt(I,l,l2,m)+SH_GetInt(I,l,l2,-m))
         W%WM(l-lmin+1,l2-lmin+1)=0
        end if
      else
        W%WP(l-lmin+1,l2-lmin+1)=-(SH_GetInt(I,l,l2,m)+SH_GetInt(I,l,l2,-m))/2
        W%WM(l-lmin+1,l2-lmin+1)=-(SH_GetInt(I,l,l2,m)-SH_GetInt(I,l,l2,-m))/2
      end if
      end do
     end do
     
 else
  stop 'SH_GetCouplingMatrices: Unsupported spin'
 end if

 call SH_FreemIntegrals(I)

  end subroutine SH_GetCouplingMatrices

  subroutine SH_FreeCouplingMatrices(W)
    Type(SH_Couplings) W

    if (associated(W%WP)) deallocate(W%WM,W%WP)
    nullify(W%WM,W%WP)
    if (associated(W%W)) deallocate(W%W)
    nullify(W%W)
    
  end subroutine SH_FreeCouplingMatrices

  subroutine SH_NullCouplingMatrices(W)
    Type(SH_Couplings) W

    nullify(W%WM,W%WP,W%W)
    
  end subroutine SH_NullCouplingMatrices



end  module SpinHarmonics
