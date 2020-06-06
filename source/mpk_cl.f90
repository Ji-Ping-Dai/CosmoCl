    !Module storing observed matter power spectrum datasets, their points and window functions
    !and routines for computing the likelihood

    !This code is based on that in cmbdata.f90
    !and on Sam Leach's incorporation of Max Tegmark's SDSS code
    !
    !Originally SLB Sept 2004
    !AL April 2006: added covariance matrix support (following 2df 2005)
    !LV_06 : incorporation of LRG DR4 from Tegmark et al . astroph/0608632
    !AL: modified LV SDSS to do Q and b^2 or b^2*Q marge internally as for 2df
    !BR09: added model LRG power spectrum.
    !AL Oct 20: switch to Ini_Read_xxx_File; fortran compatibility changes

    !JD 29/07/2013 removed LRG stuff
    !JD 03/08/2013 fixed compute_scaling_factor and associated functions
    !to work with w_a/=0

    !JD 09/13: Replaced compute_scaling_factor routines with routines that use CAMB's
    !          built in D_V function.

    !JD 02/14  CosmoTheory changes;  Added MPK_Common

    module MPK_Common
    use settings
    use CosmologyTypes
    use CosmoTheory
    use Calculator_Cosmology
    use Likelihood_Cosmology
    implicit none
    private

    Type TPKLikelihoodCommon
        real(mcp), pointer, dimension(:,:) :: mpk_win, mpk_invcov, mpk_W
        real(mcp), pointer, dimension(:) :: mpk_P,  mpk_sdev, mpk_k
    end type TPKLikelihoodCommon

    type, extends(TCosmoCalcLikelihood) :: TCosmologyPKLikelihood
        real(mcp) DV_fid   !Fiducial D_V
    contains
    procedure :: compute_scaling_factor
    end type TCosmologyPKLikelihood

    public TCosmologyPKLikelihood, TPKLikelihoodCommon
    contains

    function compute_scaling_factor(this,z,CMB)
    Class(TCosmologyPKLikelihood) this
    Class(CMBParams) CMB
    real(mcp), intent(in) :: z
    real(mcp) :: compute_scaling_factor

    !We use H_0*D_V because we dont care about scaling of h since
    !k is in units of h/Mpc
    compute_scaling_factor = this%DV_fid/(CMB%H0*this%Calculator%BAO_D_v(z))

    end function compute_scaling_factor

    end module MPK_Common

    module mpk
    use settings
    use CosmologyTypes
    use CosmoTheory
    use likelihood
    use MatrixUtils
    use MPK_Common
    use CAMB, only : grhoc,grhob,grhom,ScalarPower
    implicit none
    private

    type, extends(TCosmologyPKLikelihood) :: MPKLikelihood
        logical :: use_set
        integer :: num_mpk_points_use ! total number of points used (ie. max-min+1)
        integer :: min_mpk_points_use ! min ponit used
        integer :: max_mpk_points_use ! max ponit used
        integer :: num_mpk_kbands_use ! total number of kbands used (ie. max-min+1)
        logical :: use_scaling
        logical :: Q_marge, Q_flat
        logical :: CMBXLSS, LSSXLSS
        real(mcp) :: Q_mid, Q_sigma, Ag
        Type(TPKLikelihoodCommon) :: PKData
    contains
    procedure :: LogLike => MPK_Lnlike
    procedure :: ReadIni => MPK_ReadIni
    end type MPKLikelihood
    
    integer, parameter :: mpk_d = kind(1.d0)
    logical :: Use_mpk = .true.
    integer, parameter :: trace_cmb=0,trace_nvss=1
    double precision, parameter :: fpi = 3.1415926535897932384626433832795, fpi4=4*fpi
    double precision, parameter :: fG=6.6726d-11, kappa=2*fpi4*fG

    public use_mpk, MPKLikelihood, MPKLikelihood_Add

    contains

    subroutine MPKLikelihood_Add(LikeList, Ini)
    use IniObjects
    use settings
    implicit none
    class(TLikelihoodList) :: LikeList
    class(TSettingIni) :: ini
    type(MPKLikelihood), pointer :: this
    integer nummpksets, i

    use_mpk = (Ini%Read_Logical('use_mpk',.false.))

    if(.not. Use_mpk) return

    nummpksets = Ini%Read_Int('mpk_numdatasets',0)
    do i= 1, nummpksets
        allocate(this)
        this%LikelihoodType = 'MPK'
        this%needs_powerspectra = .true.
        this%num_z = 1
        call this%ReadDatasetFile(Ini%ReadFileName(numcat('mpk_dataset',i)) )
        call LikeList%Add(this)
    end do
    if (Feedback>1 .and. nummpksets>0) write(*,*) 'read MPK data sets'

    end subroutine MPKLikelihood_Add

    subroutine MPK_ReadIni(this,Ini)
    use MatrixUtils
    class(MPKLikelihood) this
    class(TSettingIni) :: Ini
    character(LEN=:), allocatable :: measurements_file, cov_file, windows_file

    integer i,iopb
    real(mcp) :: ell
    integer :: num_mpk_points_full ! actual number of bandpowers in the infile
    integer :: max_mpk_points_use ! in case you don't want the smallest scale modes (eg. sdss)
    integer :: min_mpk_points_use ! in case you don't want the largest scale modes

    Type(TTextFile) :: F

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!correlation!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    this%use_set =.true.
    if (Feedback > 0) write (*,*) 'reading: '//trim(this%name)

    min_mpk_points_use = Ini%Read_Int('min_mpk_points_use',1)
    max_mpk_points_use = Ini%Read_Int('max_mpk_points_use',1)
    num_mpk_points_full = max_mpk_points_use - min_mpk_points_use +1
    this%num_mpk_points_use = num_mpk_points_full
    this%min_mpk_points_use = min_mpk_points_use
    this%max_mpk_points_use = max_mpk_points_use

    this%CMBXLSS = Ini%Read_Logical('CMB_cross_LSS',.True.)
    this%LSSXLSS = Ini%Read_Logical('LSS_cross_LSS',.True.)

    allocate(this%PKData%mpk_P(num_mpk_points_full))
    
    measurements_file  = Ini%ReadFileName('measurements_file')
    if (Feedback > 0) write (*,*) 'reading: '//trim(measurements_file)
    call F%Open(measurements_file)
    this%PKData%mpk_P=0.    

    do i = 1, num_mpk_points_full
        read (F%unit,*, iostat=iopb) ell, this%PKData%mpk_P(i)
    end do

    cov_file  = Ini%ReadFileName('cov_file')
    if (cov_file /= '') then
        allocate(this%PKData%mpk_invcov(num_mpk_points_full,num_mpk_points_full))
        call File%ReadTextMatrix(cov_file,this%PKData%mpk_invcov,num_mpk_points_full,num_mpk_points_full)
        call Matrix_Inverse(this%PKData%mpk_invcov) !inverse if you input cov rather than incov
    end if
    
    windows_file  = Ini%ReadFileName('windows_file')
    if (windows_file /= '') then
        allocate(this%PKData%mpk_win(200,2))
        call File%ReadTextMatrix(windows_file,this%PKData%mpk_win,200,2)
    end if

    end subroutine MPK_ReadIni


    function MPK_LnLike(this,CMB,Theory,DataParams) 
    implicit none
    Class(CMBParams) CMB
    Class(MPKLikelihood) :: this
    Class(TCosmoTheoryPredictions), target :: Theory
    
    real(mcp) :: DataParams(:)
    real(mcp) MPK_LnLike,LnLike
    double precision, dimension(:), allocatable :: mpk_Pth, mpk_WPth
    double precision :: covth(this%num_mpk_points_use)
    double precision :: normV
  
    integer :: i,j
    
    integer nstep
    integer, parameter :: NMAX=50,NSTPMX=1001,nvar=2
    double precision xx(NSTPMX),y(NMAX,NSTPMX)
    double precision x1,x2,vstart(nvar)
    common /path/ xx,y

    double precision :: omev,omk,wq,wqpr
    common /para/ omev,omk,wq,wqpr

    allocate(mpk_Pth(this%num_mpk_points_use))
    allocate(mpk_WPth(this%num_mpk_points_use))

    wq = CMB%w
    wqpr = CMB%wa
    omk = CMB%omk
    omev = CMB%omv

    x1=1./(1.+1100.)
    x2=1./(1.)
    nstep=NSTPMX-1
    vstart(1)=1.d0
    vstart(2)=0.d0
    call rkdumb(vstart,x1,x2,nstep,GF_derivs)

    if (this%CMBXLSS) then
        call CMBXCl(this,CMB,Theory,mpk_Pth)
    else if (this%LSSXLSS) then
        call LSSXCl(this,CMB,Theory,mpk_Pth) 
    end if
    if (associated(this%PKData%mpk_invcov)) then
        mpk_WPth = mpk_Pth - this%PKData%mpk_P
        covth = matmul(this%PKData%mpk_invcov,mpk_WPth)
        normV = sum(mpk_WPth*covth)
    end if

    LnLike = normV/2

    MPK_LnLike=LnLike

    if (Feedback>1) write(*,*) 'mpk chi-sq:', LnLike*2

    end function MPK_LnLike


    function derivs_com_dis(z)
        implicit double precision (a-h,o-z)
        implicit integer (i-n)
        double precision :: omev,omk,wq,wqpr
        common /para/ omev,omk,wq,wqpr

        omem=1.0-omev-omk 
        derivs_com_dis=(omk*(1+z)**2+omem*(1+z)**3+omev)**(-0.5)
    !    derivs_com_dis=derivs_com_dis*2997.9            !Mpc/h
        return
    end function derivs_com_dis
    
    subroutine LSSXCl(this,CMB,Theory,Cl_th)
        implicit none
        Class(CMBParams) CMB 
        Class(MPKLikelihood) :: this
        Class(TCosmoTheoryPredictions), target :: Theory
        Type(TCosmoTheoryPK), pointer :: PK
        real(mcp) :: com_dis, dcom_dis, dVc
        real(mcp) :: a, GFatA, GFat0, fkh, fmpk, dz, Cl_th(this%num_mpk_points_use)
        integer :: i, j
        
        double precision :: omev,omk,wq,wqpr
        common /para/ omev,omk,wq,wqpr 

        integer, parameter :: NMAX=50,NSTPMX=1001,nvar=2
        double precision xx(NSTPMX),y(NMAX,NSTPMX)
        common /path/ xx,y

        PK=>Theory%MPK
        Cl_th=0
        
        dz=(this%PKData%mpk_win(200,1)-this%PKData%mpk_win(1,1))/199

        !open(file='pk.dat',unit=188)
        do i=1,this%num_mpk_points_use
            do j=1,200
                a=1/(1+this%PKData%mpk_win(j,1))
                call ENLGR(xx,y(1,:),NSTPMX,a,GFatA)
                GFat0=y(1,1001)
                GFatA=a*GFatA

                call qromb2(derivs_com_dis,0.d0,this%PKData%mpk_win(j,1),1d-5,1d-5,com_dis)
                com_dis=com_dis*2997.9 !omk=0
                dcom_dis=derivs_com_dis(this%PKData%mpk_win(j,1))*2997.9
                dVc=(com_dis)**2*dcom_dis
                fkh=(this%min_mpk_points_use-1+i+0.5)/com_dis
                fmpk=PK%PowerAt(fkh,0._mcp) !pk at z=0
                !for convenience, I set bias to be a constant, you should model bias more accurate see:
                !10.1016/j.physrep.2017.12.002 tabel.7
                !I also did not include fnl term. you can modify the code to calculate different models.
                Cl_th(i)=Cl_th(i)+fmpk*(this%PKData%mpk_win(j,2)*CMB%bias(1)*GFatA/GFat0)**2/dVc*dz
            end do
            !write(188,"(1I,1e18.9)") this%min_mpk_points_use-1+i,Cl_th(i)
        end do
        !close(188)
        return
    end subroutine

    subroutine CMBXCl(this,CMB,Theory,Cl_th)
        implicit none
        Class(CMBParams) CMB
        Class(MPKLikelihood) :: this
        Class(TCosmoTheoryPredictions), target :: Theory
        Type(TCosmoTheoryPK), pointer :: PK 
        real(mcp) :: com_dis, Tcmb=2.7255
        real(mcp) :: a, GFatA, dGFatA, GFat0, fkh, fmpk, dz, Cl_th(this%num_mpk_points_use)
        integer :: i, j

        double precision :: omev,omk,wq,wqpr 
        common /para/ omev,omk,wq,wqpr 
        integer, parameter :: NMAX=50,NSTPMX=1001,nvar=2
        double precision xx(NSTPMX),y(NMAX,NSTPMX)
        common /path/ xx,y 

        PK=>Theory%mpk
        Cl_th=0 
        dz=(this%PKData%mpk_win(200,1)-this%PKData%mpk_win(1,1))/199 
        do i=1,this%num_mpk_points_use
            do j=1,200
                a=1/(1+this%PKData%mpk_win(j,1))
                call ENLGR(xx,y(1,:),NSTPMX,a,GFatA)
                GFat0=y(1,1001)
                GFatA=a*GFatA
                call ENLGR(xx,y(2,:),NSTPMX,a,dGFatA)
                dGFatA = -a**2*dGFatA
                call qromb2(derivs_com_dis,0.d0,this%PKData%mpk_win(j,1),1d-5,1d-5,com_dis)
                com_dis=com_dis*2997.9 !omk=0 
                fkh=(this%min_mpk_points_use-1+i+0.5)/com_dis
                fmpk=PK%PowerAt(fkh,0._mcp)
                Cl_th(i)=Cl_th(i)+fmpk*this%PKData%mpk_win(j,2)*CMB%bias(1)*GFatA/GFat0*dGFatA/derivs_com_dis(this%PKData%mpk_win(j,1))*dz
            end do
            Cl_th(i)=Cl_th(i)*3*(1-omev-omk)*Tcmb/2997.9**3/(this%min_mpk_points_use-1+i+0.5)**2*CMB%h**3
        end do
        return
    end subroutine

    subroutine ENLGR2(x1a,x2a,ya,m,n,x1,x2,y)
    integer m,n,nmax,mmax
    real(mcp):: x1,x2,y,x1a(m),x2a(m),ya(m,n)
    parameter (nmax=100,mmax=100)
    integer j,k
    real(mcp) :: ymtmp(nmax),yntmp(mmax)
    do j=1,m
        do k=1,n
            yntmp(k)=ya(j,k)
        end do
        call ENLGR(x2a,yntmp,n,x2,ymtmp(j))
    end do
    call ENLGR(x1a,ymtmp,m,x1,y)
    end subroutine ENLGR2



	subroutine rk4(y,dydx,x,h,yout,derivs)
	implicit double precision (a-h,o-z)
	integer n,NumMAX
	parameter (NumMAX=50,n=2)
	double precision h,x,dydx(n),yout(n),y(n)
	external derivs
	integer i
	double precision h6,hh,xh,dym(NumMAX),dyt(NumMAX),yt(NumMAX)
	
	hh=h*0.5
	h6=h/6
	xh=x+hh

	do i=1,n
		yt(i)=y(i)+hh*dydx(i)
	end do
	call derivs(xh,yt,dyt)
	do i=1,n
		yt(i)=y(i)+hh*dyt(i)
	end do
	call derivs(xh,yt,dym)
	do i=1,n
		yt(i)=y(i)+h*dym(i)
		dym(i)=dyt(i)+dym(i)
	end do
	call derivs(x+h,yt,dyt)
	do i=1,n
		yout(i)=y(i)+h6*(dydx(i)+dyt(i)+2.0*dym(i))
	end do

	return
	end subroutine rk4

    subroutine rkdumb(vstart,x1,x2,nstep,derivs)
	implicit double precision (a-h,o-z)
	integer nstep  ,nvar,NMAX,NSTPMX
	parameter (NMAX=50,NSTPMX=1001,nvar=2)
	double precision x1,x2,vstart(nvar),xx(NSTPMX),y(NMAX,NSTPMX)
	external derivs
	common /path/ xx,y
	integer i,k
	double precision h,x,dv(NMAX),v(NMAX)

	do i=1,nvar
		v(i)=vstart(i)
		y(i,1)=v(i)
	end do
	xx(1)=x1
	x=x1
	h=(x2-x1)/nstep
	do k=1,nstep
		call derivs(x,v,dv)
		call rk4(v,dv,x,h,v,derivs)
		if (x+h.eq.x) pause 'stepsize not significant in rkdumb'
		x=x+h
		xx(k+1)=x
		do i=1,nvar
			y(i,k+1)=v(i)			
		end do
	end do

	return
	end subroutine rkdumb

	subroutine GF_derivs(x,y,dydx)
        use Precision
        use ModelParams

	implicit double precision (a-h,o-z)
	double precision x,y(*),dydx(*)
	double precision dtauda
        double precision wa,omx,omxa,wq,wqpr,z,omev,omk
	external dtauda
        common /para/ omev,omk,wq,wqpr
!	common /bg/ omx,wq,wqpr
        
!        omx=1-(grhoc+grhob)/grhom
	omx = omev
	z=1./x-1.d0

        wa=-1.+1./(1.-wq*(1.+z)**(-wqpr)/(1.+wq))
        omxa=((1.+wq)*(1.+z)**(wqpr)-wq)**(3./wqpr)*omx*(dtauda(1./(1.+z))/(1.+z)**2/dtauda(1.d0))**2
	dydx(1)=y(2)
        dydx(2)=-(3.5-1.5*omxa*wa)*y(2)/x-1.5*(1-omxa*wa-(1-omxa))*y(1)/x**2

	return
	end subroutine GF_derivs

	SUBROUTINE ENLGR(X,Y,N,T,Z)
	dimension X(N),Y(N)
	DOUBLE PRECISION X,Y,T,Z,S
	integer :: i,j,k,m,n

	Z=0.0
	IF (N.LE.0) RETURN
	IF (N.EQ.1) THEN
	  Z=Y(1)
	  RETURN
	END IF
	IF (N.EQ.2) THEN
	  Z=(Y(1)*(T-X(2))-Y(2)*(T-X(1)))/(X(1)-X(2))
	  RETURN
	END IF
	I=1
10	IF (X(I).LT.T) THEN
	  I=I+1
	  IF (I.LE.N) GOTO 10
	END IF
	K=I-4
	IF (K.LT.1) K=1
	M=I+3
	IF (M.GT.N) M=N
	DO 30 I=K,M
	  S=1.0
	  DO 20 J=K,M
	    IF (J.NE.I) THEN
	      S=S*(T-X(J))/(X(I)-X(J))
	    END IF
20	  CONTINUE
	  Z=Z+S*Y(I)
30	CONTINUE
	RETURN
	END SUBROUTINE ENLGR

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

           subroutine GF_ValsAta(a,gf)
            implicit none
            !Do interpolation for phi and phidot at a
        integer NMAX,NSTPMX
        parameter (NMAX=50,NSTPMX=1001)
        double precision xx(NSTPMX),y(NMAX,NSTPMX)
            double precision da,a, gf
            double precision a0,b0,ho2o6,a03,b03
            common /path/ xx,y
            integer ix
           
            da=1./(NSTPMX-1.)
            ix = int(a/da)+1
            a0 = (ix*da -a)/da
            b0 = (a-(ix-1)*da)/da
            ho2o6 = da**2/6.
            a03=(a0**3-a0)
            b03=(b0**3-b0)
            gf=a0*y(1,ix)+b0*y(1,ix+1)+(a03*y(2,ix)+b03*y(2,ix+1))*ho2o6
        end subroutine GF_ValsAta

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! num rec routines
!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE qromb2(func,a,b,epsabs,epsrel,ss)
! The numerical recipes routine, but modified so that is decides
! it's done when either the relative OR the absolute accuracy has been attained.
! The old version used relative errors only, so it always failed when
! when the integrand was near zero.
! epsabs = epsrel = 1e-6 are canonical choices.
  INTEGER JMAX,JMAXP,K,KM
  real(mpk_d) a,b,func,ss,epsabs,epsrel
  EXTERNAL func
  PARAMETER (JMAX=20, JMAXP=JMAX+1, K=5, KM=K-1)
                                !    USES polint,trapzd
  INTEGER j
  real(mpk_d) dss,h(JMAXP),s(JMAXP)
  h(1)=1.d0
  do j=1,JMAX
     call trapzd(func,a,b,s(j),j)
     if (j.ge.K) then
        call polint(h(j-KM),s(j-KM),K,0.d0,ss,dss)
        if (abs(dss).le.epsrel*abs(ss)) return
        if (abs(dss).le.epsabs) return
     endif
     s(j+1)=s(j)
     h(j+1)=0.25d0*h(j)
  ENDDO
  print *,'Too many steps in qromb'
      
  RETURN 
END SUBROUTINE qromb2
  
SUBROUTINE polint(xa,ya,n,x,y,dy) ! From Numerical Recipes
  INTEGER n,NMAX
  real(mpk_d) dy,x,y,xa(n),ya(n)
  PARAMETER (NMAX=10)
  INTEGER i,m,ns
  real(mpk_d) den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)
  ns=1
  dif=abs(x-xa(1))
  do  i=1,n
     dift=abs(x-xa(i))
     if (dift.lt.dif) then
        ns=i
        dif=dift
     endif
     c(i)=ya(i)
     d(i)=ya(i)
  enddo
  y=ya(ns)
  ns=ns-1
  do  m=1,n-1
     do  i=1,n-m
        ho=xa(i)-x
        hp=xa(i+m)-x
        w=c(i+1)-d(i)
        den=ho-hp
        if(den.eq.0.) then
           print*, 'failure in polint'
           stop
        endif
        den=w/den
        d(i)=hp*den
        c(i)=ho*den
     enddo
     if (2*ns.lt.n-m)then
        dy=c(ns+1)
     else
        dy=d(ns)
        ns=ns-1
     endif
     y=y+dy
  enddo
  return
END SUBROUTINE polint
      
SUBROUTINE trapzd(func,a,b,s,n) ! From Numerical Recipes
  INTEGER n
  real(mpk_d) a,b,s,func
  EXTERNAL func
  INTEGER it,j
  real(mpk_d) del,sum,tnm,x
  if (n.eq.1) then
     s=0.5*(b-a)*(func(a)+func(b))
  else
     it=2**(n-2)
     tnm=it
     del=(b-a)/tnm
     x=a+0.5*del
     sum=0.
     do  j=1,it
        sum=sum+func(x)
        x=x+del
     enddo
     s=0.5*(s+(b-a)*sum/tnm)
  endif
  return
END SUBROUTINE trapzd


    end module mpk
