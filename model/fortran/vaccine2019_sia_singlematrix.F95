! Title: vaccine mathematical model measles transmission
! Author: Ed Jones, London School of Hygiene & Tropical Medicine;
!      Adapted by Kevin van Zandvoort and Petra Klepac, London School of Hygiene & Tropical Medicine
! Date: 2019/11/25
! Description: Reads parameters from R-script. Age-stratified model is ran for given amount of years

!SETUP

module MeaslesFunctions
	implicit none

contains

	function vaceffByAge(a,b,t)

		integer, parameter :: dp = selected_real_kind(15, 307) !double precision
		real(dp) :: vaceffByAge
		real(dp) :: a           ! regression intercept
		real(dp) :: b           ! regression lope
		integer :: t            ! age in weeks
		real(dp) :: m           ! age in months

		!convert age to months
		if(t < 157) then
			!t is in weeks
			m = real(t, dp)/(real(52, dp)/real(12))
		else
			!t is in years
			m = 12*(3+t-156)
		endif

		! vacc_eff = a + b* age -> linear regression from Hughes et al paper
		vaceffByAge = a+b*m
		if( vaceffByAge .gt. real(0.98, dp) ) then
			vaceffByAge = real(0.98, dp)
		endif

		return

	end function vaceffByAge

	subroutine SIAcov(a, b, m1cov, scov, agefrom, ageto, pop, susceptibles, unvacI, unvacR, sia1, sia2)

		integer, parameter :: dp = selected_real_kind(15, 307) !double precision
		!real, dimension(2)  :: SIAcov     ! two outputs - coverage in zero reach, and coverage in everyone else
		real(dp), intent(in) :: a        ! regression intercept
		real(dp), intent(in) :: b        ! regression slope
		real(dp), intent(in) :: m1cov     ! input - mcv1 coverage
		real(dp), intent(in) :: scov     ! input - overall SIA coverage
		integer, intent(in)  :: agefrom
		integer, intent(in)  :: ageto
		real(dp)             :: tempreach
		real(dp)             :: doses
		real(dp)             :: zeroreach
		real(dp), dimension(:), intent(in)  :: pop
		real(dp), dimension(:), intent(in)  :: susceptibles
		real(dp), dimension(:), intent(in)  :: unvacI
		real(dp), dimension(:), intent(in)  :: unvacR
		real(dp), intent(out)  :: sia1
		real(dp), intent(out)  :: sia2
		real(dp)             :: zero_dose_children
		real(dp)             :: non_zero_dose_children

    ! --------------------------------------------------------------------------
		! find what the SIA reach should be from the simple linear regression
		! tempreach = a - b*m1cov
		! tempreach = a + b*m1cov
		! if this number is greater than 1 set it to 1
		! if(tempreach .gt. 1) then
		!	  tempreach = real(1, dp)
		! endif
		! --------------------------------------------------------------------------
		
		! --------------------------------------------------------------------------
		! find what the SIA reach should be from the weighted logistic regression
		! tempreach = exp(a + b * m1cov) / (1 + exp(a + b * m1cov))
		tempreach = exp(a + b * scov) / (1 + exp(a + b * scov))
		! --------------------------------------------------------------------------
		
		doses = scov*sum(pop(agefrom:ageto))
		
		! --------------------------------------------------------------------------
    zero_dose_children = sum(susceptibles(agefrom:ageto)) &
                        + sum(unvacI(agefrom:ageto)) &
                        + sum(unvacR(agefrom:ageto))
                        
    zeroreach = tempreach * zero_dose_children
		
		if( doses .le. zeroreach ) then
		  sia1 = doses / zero_dose_children
			sia2 = 0
		else
		  sia1 = tempreach
		  sia2 = (doses-zeroreach)/(sum(pop(agefrom:ageto)) - zero_dose_children)
		endif
		
		if( sia1 .gt. 1) then
		  sia1 = 1
		endif
		
		if( sia2 .gt. 1) then
		  sia2 = 1
		endif
    ! --------------------------------------------------------------------------   
    
    
    ! --------------------------------------------------------------------------
    ! zeroreach = tempreach*(sum(susceptibles(agefrom:ageto)) &
		!   				+ sum(unvacI(agefrom:ageto)) &
		!					+ sum(unvacR(agefrom:ageto)))
							
		! if( doses .le. zeroreach ) then
		!   sia1 = scov
		!	  sia2 = 0
		! else
		!   sia1 = tempreach
		! 	sia2 = (doses-zeroreach)/(sum(pop(agefrom:ageto))&
		!																 -(sum(susceptibles(agefrom:ageto)) &
		!																		+ sum(unvacI(agefrom:ageto)) &
		!																		+ sum(unvacR(agefrom:ageto))))
		! endif
		! --------------------------------------------------------------------------

    ! --------------------------------------------------------------------------
    ! zero_dose_children = sum(susceptibles(agefrom:ageto)) &
    !                    + sum(unvacI(agefrom:ageto)) &
    !                    + sum(unvacR(agefrom:ageto))
    
    ! non_zero_dose_children = sum(pop(agefrom:ageto)) - zero_dose_children
    
    ! if ((doses * tempreach) .le. zero_dose_children) then
    !   sia1 = (doses * tempreach) / zero_dose_children
    !   sia2 = (doses * (1 - tempreach)) / non_zero_dose_children
    ! else
    !   sia1 = 1
    !   sia2 = (doses - zero_dose_children) / non_zero_dose_children
    ! endif
    ! --------------------------------------------------------------------------

		return

	end subroutine SIAcov

end module MeaslesFunctions



program measles
  use MeaslesFunctions
	!options
  implicit none

	!declare and initialize variables
		!integer, parameter :: dp = selected_real_kind(6, 37) !single precision
		integer, parameter :: dp = selected_real_kind(15, 307) !double precision

		!parameters which cannot change
		!integer, parameter :: uu = 102
		!integer, parameter :: tmax = 101000

		integer, parameter :: amax = 254

		character(len=3) :: get_os
		character(len=1) :: sep
		character(len=4) :: get_uu
		character(len=5) :: get_tstep
		character(len=10) :: get_scenarioname
		character(len=20) :: get_foldername
		character(len=3) :: get_country_iso3
		character(len=4) :: get_psa_run
		character(len=1) :: get_debug_compartments
		character(len=1) :: get_debug_debug
		character(len=1) :: get_debug_age
		character(len=1) :: get_debug_timepoints
		character(len=1) :: get_debug_relative

		integer :: uu
		integer :: tmax
		integer :: psa_run

		integer :: debug_compartments
		integer :: debug_debug
		integer :: debug_age
		integer :: debug_timepoints
		integer :: debug_relative

		!other values used
		integer :: maternal_immunity !months of maternal immunity
		integer :: t
		integer :: tstep
		integer :: leng
		integer :: j
		integer :: a
		integer :: k
		integer :: is
		integer :: counter
		character(len=81) :: filename
		character(len=97) :: outfile
		character(len=81) :: files
		integer :: hh
		real(dp) :: pop(amax+1)
		real(dp) :: bet(amax,amax+1)
		real(dp), allocatable :: coverage(:)
		real(dp) :: cyc(amax+1)
		integer, allocatable :: t_sia(:)
		integer, allocatable :: t_run(:)
		integer :: h
		integer :: x
		real(dp), allocatable :: M(:,:)
		real(dp), allocatable :: S(:,:)
		real(dp), allocatable :: I(:,:)
		real(dp), allocatable :: R(:,:)
		real(dp), allocatable :: VS(:,:)
		real(dp), allocatable :: VI(:,:)
		real(dp), allocatable :: VR(:,:)
		real(dp), allocatable :: V2S(:,:)
		real(dp), allocatable :: V2I(:,:)
		real(dp), allocatable :: V2R(:,:)
		real(dp), allocatable :: V3S(:,:)
		real(dp), allocatable :: V3I(:,:)
		real(dp), allocatable :: V3R(:,:)
		real(dp), allocatable :: MCV1_eff(:,:)
		real(dp), allocatable :: newinfect(:,:)
		real(dp), allocatable :: newvac(:,:)
		real(dp), allocatable :: N(:,:)
		real(dp), allocatable :: infected(:,:)
		real(dp), allocatable :: infecteds(:)
		real(dp), allocatable :: infecteds1(:)
		real(dp), allocatable :: infectedlt5(:)
		real(dp), allocatable :: infectedgt5(:)
		real(dp), allocatable :: heights(:)
		real(dp) :: amp
		real(dp) :: tcycle!Seasonality in Beta
		real(dp) :: dI
		real(dp) :: dVI
		real(dp) :: dV2I
		real(dp) :: dV3I
		real(dp) :: dR
		real(dp) :: dVR
		real(dp) :: dV2R
		real(dp) :: dV3R
		real(dp) :: degree1
		real(dp) :: degree2
		real(dp) :: degree3
		real(dp) :: g
		real(dp) :: cov(amax)
		real(dp) :: cov2
		real(dp), allocatable :: sia(:)
		real(dp) :: take1
		real(dp) :: take2
		real(dp) :: take3
		real(dp) :: take13
		real(dp) :: take23
		integer, allocatable  :: a0(:)
		integer, allocatable  :: a1(:)
!		real :: dM_S
    real(dp) :: dS_VS
    real(dp) :: dS_VR
    real(dp) :: dI_VI
    real(dp) :: dR_VR
	  real(dp) :: dS_V2S
		real(dp) :: dS_V2R
		real(dp) :: dI_V2I
		real(dp) :: dR_V2R
		real(dp) :: dVS_V2S
		real(dp) :: dVS_V2R
		real(dp) :: dVI_V2I
		real(dp) :: dVR_V2R
		real(dp) :: dVS_V3S
		real(dp) :: dVS_V3R
		real(dp) :: dVI_V3I
		real(dp) :: dVR_V3R
		real(dp) :: dV2S_V3S
		real(dp) :: dV2S_V3R
		real(dp) :: dV2I_V3I
		real(dp) :: dV2R_V3R
		real(dp) :: tempp
		real(dp) :: r0target
		real(dp) :: doSIA
		real(dp) :: t_end
		real(dp) :: betta
		integer :: t_cov
		integer :: sia_n
		integer :: temp
		integer :: tstart
		integer :: tend
		integer :: t_scale
		integer :: t_increment
		integer :: astart
		integer :: aend
		real(dp) :: agerate_annual
		real(dp) :: agerate_week
		real(dp) :: totpop
!		real(dp) :: vaceffByAge
!    real(dp), dimension(2) :: SIAcov
!		real(dp), dimension(2) :: SIAoutput
		real(dp) :: sia1
		real(dp) :: sia2
		real(dp) :: take_inc
		real(dp) :: ve_a
		real(dp) :: ve_b
		real(dp) :: SIA_a
    real(dp) :: SIA_b
		real(dp) :: vac_numerator
		real(dp) :: vac_denominator
		real(dp), allocatable :: tempUnvacS(:)
		real(dp), allocatable :: tempUnvacI(:)
		real(dp), allocatable :: tempUnvacR(:)

	!start model
	print*, 'Running MeaslesAuto'

	!initialize variables
	call getarg(1, get_uu)				!number of years (string)
	read(get_uu,*) uu					!number of years (integer)

	call getarg(2, get_tstep)			!tstep (string)
	read(get_tstep,*) tstep				!tstep (integer)

	call getarg(3, get_scenarioname) !scenarioname (string, 10 characters)

	call getarg(4, get_foldername)		!foldername (string)

	call getarg(5, get_country_iso3)	!country (string)

	call getarg(6, get_psa_run)			!psa_run (string)
	read(get_psa_run,*) psa_run			!psa_run (integer)

	call getarg(7, get_os)				!OS, used for seperator
	if( get_os == "WIN" ) then
		sep = "\"
	else
		sep = "/"
	endif

	call getarg(8, get_debug_debug)			!1 or 0, to export cases per timestep
	read(get_debug_debug,*) debug_debug

	call getarg(9, get_debug_compartments)	!1 or 0, to export cases per timestep
	read(get_debug_compartments,*) debug_compartments

	call getarg(10, get_debug_age)		!1 or 0, to export cases per timestep
	read(get_debug_age,*) debug_age

	call getarg(11, get_debug_timepoints)	!0, 1 or 2, to export cases per timestep
	read(get_debug_timepoints,*) debug_timepoints

	call getarg(12, get_debug_relative)		!1 or 0, to export cases per timestep
	read(get_debug_relative,*) debug_relative

	tmax = uu*tstep
	
	! ----------------------------------------------------------------------------
	uu = uu + 1

	! allocate(coverage(uu))
	! allocate(t_sia(uu))
	! allocate(t_run(uu))
	! allocate(sia(uu))
	! allocate(a0(uu))
	! allocate(a1(uu))
	
	! memory allocation should be large enough to handle multiple SIAs over multiple years
	allocate(coverage(uu*2))
	allocate(t_sia(uu*2))
	allocate(t_run(uu*2))
	allocate(sia(uu*2))
	allocate(a0(uu*2))
	allocate(a1(uu*2))
	! ----------------------------------------------------------------------------
	
	allocate(M(amax+1,tmax+1))
	allocate(S(amax+1,tmax+1))
	allocate(I(amax+1,tmax+1))
	allocate(R(amax+1,tmax+1))
	allocate(VS(amax+1,tmax+1))
	allocate(VI(amax+1,tmax+1))
	allocate(VR(amax+1,tmax+1))
	allocate(V2S(amax+1,tmax+1))
	allocate(V2I(amax+1,tmax+1))
	allocate(V2R(amax+1,tmax+1))
	allocate(V3S(amax+1,tmax+1))
	allocate(V3I(amax+1,tmax+1))
	allocate(V3R(amax+1,tmax+1))

	!individuals for who the MCV1 vaccine was effective (VR who became immune following vaccination, not following natural immunity)
	allocate(MCV1_eff(amax+1,tmax+1))

	allocate(newinfect(amax+1,tmax+1))
	allocate(newvac(amax+1,tmax+1))
	allocate(N(amax+1,tmax+1))
	allocate(infected(amax+1,tmax+1))
	allocate(infecteds(tmax+1))
	allocate(infecteds1(tmax+1))
	allocate(infectedlt5(tmax+1))
	allocate(infectedgt5(tmax+1))
	allocate(heights(tmax+1))

	maternal_immunity	=	6 !months of maternal immunity
	sia_n				=	0
	sia					=	real(0, dp)
	t_sia				=	0
	t_run				=	0
	cov					=	real(0, dp)
	a0					=	0
	a1					=	0
	newinfect			=	real(0, dp)

	!vaccine efficacy parameters for function vaceffByAge
	!ve_a = real(0.64598,dp)
	!ve_b = real(0.01485,dp)

  ! ----------------------------------------------------------------------------
	! SIA reach regression parameters for function SIAcov 
  ! From Allison's paper: 
	!   From:
	! SIA_a = real(1.141, dp)
	! SIA_b = real(-0.679, dp)
	! SIA_a = real(1, dp)
	! SIA_b = real(-0.5, dp)
  ! SIA_a = real(1.0644, dp)
  ! SIA_b = real(-0.5983, dp)
  
  ! intercept and slope values from weighted logistic regression
  ! SIA_a = real(2.485921, dp)
  ! SIA_b = real(-2.074887, dp)
  
  ! intercept and slope values from weighted logistic regression
  SIA_a = real(-2.621733, dp)
  SIA_b = real(5.238249, dp)
  
	!
	!  To; (just for sensitivity analysis and error checking, reverting back to above)
	!  I think this will make sure all the doses go to zero dose kids until there are none left.
	!
	! SIA reach regression parameters for function SIAcov (amended 2020/05/31 to best case scenario)
  ! SIA_a = real(1, dp)
  ! SIA_b = real(0, dp)
	!
	! ----------------------------------------------------------------------------

	!initialize compartments
	M	= real(0, dp)
	S	= real(0, dp)
	MCV1_eff = real(0, dp)
	VS	= real(0, dp)
	V2S = real(0, dp)
	V3S = real(0, dp)
	I	= real(0, dp)
	VI	= real(0, dp)
	V2I = real(0, dp)
	V3I = real(0, dp)
	R	= real(0, dp)
	VR	= real(0, dp)
	V2R = real(0, dp)
	V3R = real(0, dp)

	!introduce 0.95 susceptible people in all agegroups at first timepoint
	S(:,1) = real(0.95, dp)
	!introduce 0.05 infectious people in all agegroups at first timepoint
	I(:,1) = real(0.05, dp)

	hh=2
	h=1

	!read parameters for first year from file
	counter = 1
	if(psa_run.gt.0) then
		write(files,'(A2,A1,A2,A1,A7,A1,A10,A1,A20,A1,A5,A1,A3,I3.3,A1,A3,A1,I3.3,A15)') '..',sep,'..',sep,'outcome',sep,get_scenarioname,sep,get_foldername,sep,'input',sep,'run',psa_run,sep,get_country_iso3,'_',counter,'measle_data.txt'
	else
		write(files,'(A2,A1,A2,A1,A7,A1,A10,A1,A20,A1,A5,A1,A3,A1,I3.3,A15)') '..',sep,'..',sep,'outcome',sep,get_scenarioname,sep,get_foldername,sep,'input',sep,get_country_iso3,'_',counter,'measle_data.txt'
	endif
	open(1999,file=files,status='unknown')
		!skip first two lines
		read(1999,*)
		read(1999,*)
		read(1999,*)g
		read(1999,*)tstep
		!vaccine efficacy
		read(1999,*)take1
		read(1999,*)take2
		read(1999,*)take3
		read(1999,*)degree1
		read(1999,*)degree2
		read(1999,*)degree3
		read(1999,*)leng
		do 999 j=1,leng
			read(1999,*)a0(j)
		999 continue
		read(1999,*)r0target
		read(1999,*)amp
		read(1999,*)doSIA
		read(1999,*)leng
		do 2012 j=1,leng
			read(1999,*)a1(j)
		2012 continue
		read(1999,*)t_end
		read(1999,*)t_cov
		read(1999,*)leng
		do 110 j=1,leng
			read(1999,*)t_sia(j)
		110 continue
		read (1999,*)leng
		do 170 j=1,leng
			read(1999,*)cov(j)
		170 continue
		!do 970 j=1,leng
			!read(1999,*)cov2(j)
		!970 continue
		read(1999,*)cov2
		read(1999,*)leng
		do 570 j=1,leng
			read(1999,*)sia(j)
		570 continue
		read(1999,*)leng
		do 210 j=1,leng
			read(1999,*)t_run(j)
		210 continue
		do 100 k=1,(amax)
			read(1999,*)pop(k)
		100 continue
		temp=sum(pop)
	close(1999)
	
	!vaccine efficacy parameters for function vaceffByAge
	ve_a = take1
	ve_b = take2
	
	!dont need to do this in the loop anymore
	!contact rates
	if(psa_run.gt.0) then
		write(filename,'(A2,A1,A2,A1,A7,A1,A10,A1,A20,A1,A5,A1,A3,I3.3,A1,A3,A1,I3.3,A8)') '..',sep,'..',sep,'outcome',sep,get_scenarioname,sep,get_foldername,sep,'input',sep,'run',psa_run,sep,get_country_iso3,'_',1,'beta.txt'
	else
		write(filename,'(A2,A1,A2,A1,A7,A1,A10,A1,A20,A1,A5,A1,A3,A1,I3.3,A8)') '..',sep,'..',sep,'outcome',sep,get_scenarioname,sep,get_foldername,sep,'input',sep,get_country_iso3,'_',1,'beta.txt'
	endif

	open(2000,file=filename,status='unknown')
		rewind(2000)
		do 101 k=1,amax
			read(2000,*,end=777)(bet(k,is),is=1,amax)
		101 continue
		777 continue
	close(2000)
	
	!substract t_cov from t_sia
	t_sia=t_sia-t_cov
	!set t_sia AFTER last timepoint where t_sia is negative (will never run, happens when SIA is disabled on LINUX)
	where(t_sia<0) t_sia=tmax+1

	agerate_annual	= real(1, dp)/real(tstep, dp)
	agerate_week 	= real(52, dp)/real(tstep, dp)

	!spin up the system (may need to check some countries that take longer to get to equilibrium)
	!loop through number of years in t_cov
	do 1979 t=1,(t_cov-1), 1
		h = t
		hh = t+1
		!calculate seasonality in contact-rate
		tcycle=sin(real(t, dp)*real(2, dp)*real(3.1415926535, dp)/(tstep)) !Seasonality in Beta
		tcycle=real(1, dp)+tcycle*amp
		!cyc= tcycle*(I(:,h)+VI(:,h)+V2I(:,h)+V3I(:,h))*pop(1:amax+1)
		cyc= tcycle*(I(:,h)+VI(:,h)+V2I(:,h)+V3I(:,h))
		totpop = sum(pop(1:amax+1)) ! total population size
		!loop through all age-groups
		do 1982 a=1,amax
			!calculate number of infected individuals
			betta= sum(bet(a,:)*cyc) ! removed - change: should be multiplied by the a-th column in bet matrix and divided by total population (frequency dependent mixing)
			dI		=	S(a,h)*betta
			dVI 	=	VS(a,h)*(real(1, dp) - degree1)*betta
			dV2I	=	V2S(a,h)*(real(1, dp) - degree2)*betta
			dV3I	=	V3S(a,h)*(real(1, dp) - degree3)*betta
			!calculate proportion of recovered individuals
		    !dM_S = (1/(real(maternal_immunity)/12*real(tstep)))*M(a,h)
			dR		=	g*I(a,h)
			dVR		=	g*VI(a,h)
			dV2R	=	g*V2I(a,h)
			dV3R	=	g*V3I(a,h)
			!calculate number of infected individuals by multiplying with population-size
			if (debug_relative == 0) then
				newinfect(a,h) = (dI+dVI+dV2I+dV3I)*pop(a)
			else
				newinfect(a,h) = (dI+dVI+dV2I+dV3I)
			endif
			!apply mutations to compartments
			!M(a,hh)=M(a,h)-dM_S
			!S(a,hh)=S(a,h)-dI +dM_S
			S(a,hh)		=	S(a,h)						&
							-dI

			I(a,hh)		=	I(a,h)						&
							+dI							&
							-dR

			R(a,hh)		=	R(a,h)						&
							+dR

			VS(a,hh)	=	VS(a,h)						&
							-dVI

			VI(a,hh)	=	VI(a,h)						&
							+dVI						&
							-dVR

			VR(a,hh)	=	VR(a,h)						&
							+dVR

			MCV1_eff(a,hh) = MCV1_eff(a,h)

			V2S(a,hh)	=	V2S(a,h)					&
							-dV2I

			V2I(a,hh)	=	V2I(a,h)					&
							+dV2I						&
							-dV2R

			V2R(a,hh)	=	V2R(a,h)					&
							+dV2R

			V3S(a,hh)	=	V3S(a,h)					&
							-dV3I

			V3I(a,hh)	=	V3I(a,h)					&
							+dV3I						&
							-dV3R

			V3R(a,hh)	=	V3R(a,h)					&
							+dV3R
		1982 continue
		!#Aging (except newborns) !! Maternal immunity is
		!!varied by changing the average maternal immunity e.g 0.002 = 6 months, 0.04 = 3 months ave) 12/6 * 1/tstep = 0.002
		!loop through age-groups, from oldest until those 3 years old
		do 1937 a=(amax+1),157,-1
			!move individuals with age a-1 to those with age a at timepoint t
			!M(a,hh)		=M(a,hh)	-M(a,hh)*1/real(tstep)	+M(a-1,hh)*1/real(tstep)
			!S(a,hh)		=S(a,hh)	-S(a,hh)*1/real(tstep)	+S(a-1,hh)*1/real(tstep)
			M(a,hh)		=	M(a,h)						&
							-agerate_annual*M(a,h)		&
							+agerate_annual*M(a-1,h)	&
							-real(0.002, dp)*M(a,h)

			S(a,hh)		=	S(a,hh)						&
							-agerate_annual*S(a,hh)		&
							+agerate_annual*S(a-1,hh)	&
							+real(0.002, dp)*M(a,h)

			I(a,hh)		=	I(a,hh)						&
							-agerate_annual*I(a,hh)		&
							+agerate_annual*I(a-1,hh)

			R(a,hh)		=	R(a,hh)						&
							-agerate_annual*R(a,hh)		&
							+agerate_annual*R(a-1,hh)

			VS(a,hh)	=	VS(a,hh)					&
							-agerate_annual*VS(a,hh)	&
							+agerate_annual*VS(a-1,hh)

			VI(a,hh)	=	VI(a,hh)					&
							-agerate_annual*VI(a,hh)	&
							+agerate_annual*VI(a-1,hh)

			VR(a,hh)	=	VR(a,hh)					&
							-agerate_annual*VR(a,hh)	&
							+agerate_annual*VR(a-1,hh)

			MCV1_eff(a, hh) = MCV1_eff(a,hh)				&
							-agerate_annual*MCV1_eff(a,hh)	&
							-agerate_annual*MCV1_eff(a-1,hh)

			V2S(a,hh)	=	V2S(a,hh)					&
							-agerate_annual*V2S(a,hh)	&
							+agerate_annual*V2S(a-1,hh)

			V2I(a,hh)	=	V2I(a,hh)					&
							-agerate_annual*V2I(a,hh)	&
							+agerate_annual*V2I(a-1,hh)

			V2R(a,hh)	=	V2R(a,hh)					&
							-agerate_annual*V2R(a,hh)	&
							+agerate_annual*V2R(a-1,hh)

			V3S(a,hh)	=	V3S(a,hh)					&
							-agerate_annual*V3S(a,hh)	&
							+agerate_annual*V3S(a-1,hh)

			V3I(a,hh)	=	V3I(a,hh)					&
							-agerate_annual*V3I(a,hh)	&
							+agerate_annual*V3I(a-1,hh)

			V3R(a,hh)	=	V3R(a,hh)					&
							-agerate_annual*V3R(a,hh)	&
							+agerate_annual*V3R(a-1,hh)
		1937 continue
		!loop through age-groups, from 156 weeks old to 2 weeks old
		do 1977 a=156,2,-1
			!should maternal immunity be 0.002 here?
			!this shifts those aged 54 weeks to 55 weeks, but rate 52
			!M(a,hh)		=M(a,hh)	-M(a,hh)*52/real(tstep)		+M(a-1,hh)*52/real(tstep)
			!S(a,hh)		=S(a,hh)	-S(a,hh) *52/real(tstep)	+S(a-1,hh)*52/real(tstep)
			M(a,hh)		=	M(a,h)						&
							-agerate_week*M(a,h)		&
							+agerate_week*M(a-1,h)		&
							-real(0.002, dp)*M(a,h)

			S(a,hh)		=	S(a,hh)						&
							-agerate_week*S(a,hh)		&
							+agerate_week*S(a-1,hh)		&
							+real(0.002, dp)*M(a,h)

			I(a,hh)		=	I(a,hh)						&
							-agerate_week*I(a,hh)		&
							+agerate_week*I(a-1,hh)

			R(a,hh)		=	R(a,hh)						&
							-agerate_week*R(a,hh)		&
							+agerate_week*R(a-1,hh)

			VS(a,hh)	=	VS(a,hh)					&
							-agerate_week*VS(a,hh)		&
							+agerate_week*VS(a-1,hh)

			VI(a,hh)	=	VI(a,hh)					&
							-agerate_week*VI(a,hh)		&
							+agerate_week*VI(a-1,hh)

			VR(a,hh)	=	VR(a,hh)					&
							-agerate_week*VR(a,hh)		&
							+agerate_week*VR(a-1,hh)

			MCV1_eff(a,hh) = MCV1_eff(a,hh)						&
							-agerate_week*MCV1_eff(a,hh)		&
							+agerate_week*MCV1_eff(a-1,hh)

			V2S(a,hh)	=	V2S(a,hh)					&
							-agerate_week*V2S(a,hh)		&
							+agerate_week*V2S(a-1,hh)

			V2I(a,hh)	=	V2I(a,hh)					&
							-agerate_week*V2I(a,hh)		&
							+agerate_week*V2I(a-1,hh)

			V2R(a,hh)	=	V2R(a,hh)					&
							-agerate_week*V2R(a,hh)		&
							+agerate_week*V2R(a-1,hh)

			V3S(a,hh)	=	V3S(a,hh)					&
							-agerate_week*V3S(a,hh)		&
							+agerate_week*V3S(a-1,hh)

			V3I(a,hh)	=	V3I(a,hh)					&
							-agerate_week*V3I(a,hh)		&
							+agerate_week*V3I(a-1,hh)

			V3R(a,hh)	=	V3R(a,hh)					&
							-agerate_week*V3R(a,hh)		&
							+agerate_week*V3R(a-1,hh)
		1977 continue
		!Births, newborn aging and routine vaccination (starts at t.cov)
		!NB birth rate per timestep = 1/real(tstep)
		!Maternal immunity in proportion to those infected and vaccinated (this is more complex in nature than in current model)
		M(1,hh) 	=	M(1,hh-1)																		&
							-agerate_week*M(1,hh-1)															&
							+agerate_week*sum(																&
								(																			&
									R(172:189,hh)+VR(172:189,hh)+V2R(172:189,hh)+V3R(172:189,hh)			&
								)*pop(172:189)																&
							)/sum(																			&
								(																			&
									R(172:189,hh)+VR(172:189,hh)+V2R(172:189,hh)+V3R(172:189,hh)			&
									+S(172:189,hh)+VS(172:189,hh)+V2S(172:189,hh)+V3S(172:189,hh)			&
								)*pop(172:189) 																&
							)

			S(1,hh) 	=	S(1,hh)																			&
							-agerate_week*S(1,hh)															&
							+agerate_week*sum(																&
								(																			&
									S(172:189,hh)+VS(172:189,hh)+V2S(172:189,hh)+V3S(172:189,hh)			&
								)*pop(172:189)																&
							)/sum(																			&
								(																			&
									R(172:189,hh)+VR(172:189,hh)+V2R(172:189,hh)+V3R(172:189,hh)			&
									+S(172:189,hh)+VS(172:189,hh)+V2S(172:189,hh)+V3S(172:189,hh)			&
								)*pop(172:189) 																&
							)


						I(1,hh)		=	I(1,hh)																				&
										-agerate_week*I(1,hh)

						R(1,hh)		=	R(1,hh)																				&
										-agerate_week*R(1,hh)

						VS(1,hh)	=	VS(1,hh)																			&
										-agerate_week*VS(1,hh)

						VI(1,hh)	=	VI(1,hh)																			&
										-agerate_week*VI(1,hh)

						VR(1,hh)	=	VR(1,hh)																			&
										-agerate_week*VR(1,hh)

						MCV1_eff(1, hh) = MCV1_eff(1,hh)																	&
										-agerate_week*MCV1_eff(1,hh)

						V2S(1,hh)	=	V2S(1,hh)																			&
										-agerate_week*V2S(1,hh)

						V2I(1,hh)	=	V2I(1,hh)																			&
										-agerate_week*V2I(1,hh)

						V2R(1,hh)	=	V2R(1,hh)																			&
										-agerate_week*V2R(1,hh)

						V3S(1,hh)	=	V3S(1,hh)																			&
										-agerate_week*V3S(1,hh)

						V3I(1,hh)	=	V3I(1,hh)																			&
										-agerate_week*V3I(1,hh)

						V3R(1,hh)	=	V3R(1,hh)																			&
										-agerate_week*V3R(1,hh)
		!calculate total population size for all ages
		do 1961 a=1,(amax)
			N(a,hh)=pop(a)
		1961 continue

		temp=hh

	1979 continue

	if( debug_debug == 1 .OR. debug_debug == 3 ) then
		if(psa_run.gt.0) then
			write(outfile,'(A2,A1,A2,A1,A7,A1,A10,A1,A20,A1,A6,A1,A3,I3.3,A1,A3,A15)') '..',sep,'..',sep,'outcome',sep,get_scenarioname,sep,get_foldername,sep,'output',sep,'run',psa_run,sep,get_country_iso3,'_debug_spin.txt'
		else
			write(outfile,'(A2,A1,A2,A1,A7,A1,A10,A1,A20,A1,A6,A1,A3,A15)') '..',sep,'..',sep,'outcome',sep,get_scenarioname,sep,get_foldername,sep,'output',sep,get_country_iso3,'_debug_spin.txt'
		endif
		open(9915,file=outfile,status='NEW')

		if( debug_timepoints == 0 ) then
			t_increment = tstep - 1
			t_scale = 1
		else
			t_increment = 1
			if(debug_timepoints == 1) then
				t_scale = 4
			else
				t_scale = 1
			endif
		endif
		do 1155 t=1, int(t_cov/t_scale), t_increment
			tstart = t
			if( debug_timepoints == 0 ) then
				tend = tstart+tstep
			else
				tend = tstart
			endif
			do 1144 a=1, amax, 1
				astart = a
				!move to next loop if yearly age-strata are used for first 3 years
				if ( debug_age == 0 .AND. a < 157 ) then
					if( a == 1 .OR. a == 53 .OR. a == 105 ) then
						aend = astart + 51
					else
						cycle
					endif
				elseif ( debug_age == 2) then
					if( a == 1 ) then
						aend = amax
					else
						cycle
					endif
				else
					aend = astart
				endif
				if( debug_compartments == 0 ) then
					write(9915,1070) tstart, tend, astart, aend, "cases", sum(newinfect(astart:aend,tstart:tend))
				else
					write(9915,1070) tstart, tend, astart, aend, "M", 	sum(M(astart:aend,tstart:tend))
					write(9915,1070) tstart, tend, astart, aend, "S", 	sum(S(astart:aend,tstart:tend))
					write(9915,1070) tstart, tend, astart, aend, "I", 	sum(I(astart:aend,tstart:tend))
					write(9915,1070) tstart, tend, astart, aend, "R", 	sum(R(astart:aend,tstart:tend))
					write(9915,1070) tstart, tend, astart, aend, "VS", 	sum(VS(astart:aend,tstart:tend))
					write(9915,1070) tstart, tend, astart, aend, "VI", 	sum(VI(astart:aend,tstart:tend))
					write(9915,1070) tstart, tend, astart, aend, "VR", 	sum(VR(astart:aend,tstart:tend))
	        write(9915,1070) tstart, tend, astart, aend, "MCV1_eff", 	sum(MCV1_eff(astart:aend,tstart:tend))
					write(9915,1070) tstart, tend, astart, aend, "V2S",	sum(V2S(astart:aend,tstart:tend))
					write(9915,1070) tstart, tend, astart, aend, "V2I",	sum(V2I(astart:aend,tstart:tend))
					write(9915,1070) tstart, tend, astart, aend, "V2R",	sum(V2R(astart:aend,tstart:tend))
					write(9915,1070) tstart, tend, astart, aend, "V3S",	sum(V3S(astart:aend,tstart:tend))
					write(9915,1070) tstart, tend, astart, aend, "V3I",	sum(V3I(astart:aend,tstart:tend))
					write(9915,1070) tstart, tend, astart, aend, "V3R",	sum(V3R(astart:aend,tstart:tend))
				endif
			1144 continue
		1155 continue
		close(9915)
	endif

	h=temp

	!initialize year to 0
	x=0
	!substract t_cov from t_run
	t_run=t_run-t_cov
	!set first t_run to h, which was memory-position that in arrays that is used in last-run in reaching equilibrium
	t_run(1)=1
	!set compartment_sizes at memory-position for t = 1 same as in last run in spin-up

	newinfect(:,1)	= newinfect(:,h-1)
	M(:,1)			= M(:,h)
	S(:,1)			= S(:,h)
	I(:,1)			= I(:,h)
	R(:,1)			= R(:,h)
	MCV1_eff(:,1)	= MCV1_eff(:,h)
	VS(:,1)			= VS(:,h)
	VI(:,1)			= VI(:,h)
	VR(:,1)			= VR(:,h)
	V2S(:,1)		= V2S(:,h)
	V2I(:,1)		= V2I(:,h)
	V2R(:,1)		= V2R(:,h)
	V3S(:,1)		= V3S(:,h)
	V3I(:,1)		= V3I(:,h)
	V3R(:,1)		= V3R(:,h)
	N(:,1)			= N(:,h)

	newinfect(:,2:tmax+1)	= real(0,dp)
	M(:,2:tmax+1)			= real(0,dp)
	S(:,2:tmax+1)			= real(0,dp)
	I(:,2:tmax+1)			= real(0,dp)
	R(:,2:tmax+1)			= real(0,dp)
  MCV1_eff(:,2:tmax+1)	= real(0,dp)
	VS(:,2:tmax+1)			= real(0,dp)
	VI(:,2:tmax+1)			= real(0,dp)
	VR(:,2:tmax+1)			= real(0,dp)
	V2S(:,2:tmax+1)			= real(0,dp)
	V2I(:,2:tmax+1)			= real(0,dp)
	V2R(:,2:tmax+1)			= real(0,dp)
	V3S(:,2:tmax+1)			= real(0,dp)
	V3I(:,2:tmax+1)			= real(0,dp)
	V3R(:,2:tmax+1)			= real(0,dp)
	N(:,2:tmax+1)			= real(0,dp)

	!set t_run to last run-of-model in case it is negative
	where(t_run<0) t_run=tmax
	counter = 0

	!used in goto, interpret as loop
		!go to next year
		1234 x=x+1
		!only read in new-data if year is not first year, otherwise use data thats read-in before
		if(x.gt.1)then
			!read data from current year
			counter = counter+1

			if(psa_run.gt.0) then
				write(files,'(A2,A1,A2,A1,A7,A1,A10,A1,A20,A1,A5,A1,A3,I3.3,A1,A3,A1,I3.3,A15)') '..',sep,'..',sep,'outcome',sep,get_scenarioname,sep,get_foldername,sep,'input',sep,'run',psa_run,sep,get_country_iso3,'_',counter,'measle_data.txt'
			else
				write(files,'(A2,A1,A2,A1,A7,A1,A10,A1,A20,A1,A5,A1,A3,A1,I3.3,A15)') '..',sep,'..',sep,'outcome',sep,get_scenarioname,sep,get_foldername,sep,'input',sep,get_country_iso3,'_',counter,'measle_data.txt'
			endif
			open(counter+1000,file=files,status='old')
				read(counter+1000,*)
				read(counter+1000,*)
				read(counter+1000,*)g
				read(counter+1000,*)tstep
				!take1 and 2 are never overwritten, so these will always be the same as in the first year!
				read(counter+1000,*)tempp	!take1
				read(counter+1000,*)tempp	!take2
				read(counter+1000,*)take3
				read(counter+1000,*)degree1
				read(counter+1000,*)degree2
				read(counter+1000,*)degree3
				read(counter+1000,*)leng
				do 888 j=1,leng
					read(counter+1000,*)a0(j)
				888 continue
				read(counter+1000,*)r0target
				read(counter+1000,*)amp
				read(counter+1000,*)doSIA
				read(counter+1000,*)leng
				do 7777 j=1,leng
					read(counter+1000,*)a1(j)
				7777 continue
				read(counter+1000,*)t_end
				read(counter+1000,*)t_cov
				read(counter+1000,*)leng
				do 1100 j=1,leng
					read(counter+1000,*)t_sia(j)
				1100 continue
				read (counter+1000,*)leng
				do 1170 j=1,leng
					read(counter+1000,*)cov(j)
				1170 continue
				!do 670 j=1,leng
				read(counter+1000,*)cov2!(j)
				!670 continue
				!read(counter+1000,*)cov2
				read(counter+1000,*)leng
				do 1135 j=1,leng
					read(counter+1000,*)sia(j)
				1135 continue
				read(counter+1000,*)leng
				do 1210 j=1,leng
					read(counter+1000,*)temp !!!t_run(j)
				1210 continue
				do 1000 k=1,(amax)
					read(counter+1000,*)pop(k)
				1000 continue
			close(counter+1000)
			!read contact-data
			!if(psa_run.gt.0) then
			!	write(filename,'(A2,A1,A2,A1,A7,A1,A10,A1,A20,A1,A5,A1,A3,I3.3,A1,A3,A1,I3.3,A8)') '..',sep,'..',sep,'outcome',sep,get_scenarioname,sep,get_foldername,sep,'input',sep,'run',psa_run,sep,get_country_iso3,'_',1,'beta.txt'
			!else
			!	write(filename,'(A2,A1,A2,A1,A7,A1,A10,A1,A20,A1,A5,A1,A3,A1,I3.3,A8)') '..',sep,'..',sep,'outcome',sep,get_scenarioname,sep,get_foldername,sep,'input',sep,get_country_iso3,'_',1,'beta.txt'
			!endif
			!open(counter,file=filename,status='old')
			!	do 1111 k=1,amax
			!		read(counter,*,end=787)(bet(k,is),is=1,amax)
			!	1111 continue
			!	787 continue
			!close(counter)
			!substract t_cov from t_sia
			t_sia=t_sia-t_cov
			!set t_sia AFTER last timepoint where t_sia is negative (will never run, happens when SIA is disabled on LINUX)
			where(t_sia<0) t_sia=tmax+1
		endif

		take13=(take3-take1)/(real(1, dp)-take1)
		take23=(take3-take2)/(real(1, dp)-take2)

		!set time to start of year x, loop over until start of next year
		do 1980 t=t_run(x),t_run(x+1)-1!tmax
			hh=t+1
			!take1=real(0.85, dp)
			!take2=real(0.85, dp)
			tcycle=sin(real(t, dp)*real(2, dp)*real(3.1415926535, dp)/(tstep)) !Seasonality in Beta
			tcycle=real(1, dp)+tcycle*amp
			!cyc= tcycle*(I(:,t)+VI(:,t)+V2I(:,t)+V3I(:,t))*pop(1:amax+1)
			cyc= tcycle*(I(:,t)+VI(:,t)+V2I(:,t)+V3I(:,t))
			totpop = sum(pop(1:amax+1))
			!loop through all age-groups
			do 1978 a=1,amax, 1
				!calculate number of infected individuals
				betta	=	sum(bet(a,:)*cyc)
				dI		=	S(a,t)*betta
				dVI		=	VS(a,t) *(real(1, dp)-degree1)*betta
				dV2I	=	V2S(a,t)*(real(1, dp)-degree2)*betta
				dV3I	=	V3S(a,t)*(real(1, dp)-degree3)*betta
				!calculate new removed
				!dM_S = (1/(real(maternal_immunity)/12*real(tstep)))*M(a,t)
				dR	=g*I(a,t)
				dVR =g*VI(a,t)
				dV2R=g*V2I(a,t)
				dV3R=g*V3I(a,t)
				!calculate number new infections
				if (debug_relative == 0) then
					newinfect(a,hh)	=	(dI+dVI+dV2I+dV3I)*pop(a)
				else
					newinfect(a,hh)	=	(dI+dVI+dV2I+dV3I)
				endif
				!update compartments with mutations
				!M(a,hh)		=M(a,t)		-dM_S
				!S(a,hh)		=S(a,t)		-dI		+dM_S
				S(a,hh)		=	S(a,t)			&
								-dI

				I(a,hh)		=	I(a,t)			&
								+dI				&
								-dR

				R(a,hh)		=	R(a,t)			&
								+dR

				VS(a,hh)	=	VS(a,t)			&
								-dVI

				VI(a,hh)	=	VI(a,t)			&
								+dVI			&
								-dVR

				VR(a,hh)	=	VR(a,t)			&
								+dVR

				V2S(a,hh)	=	V2S(a,t)		&
								-dV2I

	    	MCV1_eff(a,hh) = MCV1_eff(a,t)

				V2I(a,hh)	=	V2I(a,t)		&
								+dV2I			&
								-dV2R

				V2R(a,hh)	=	V2R(a,t)		&
								+dV2R

				V3S(a,hh)	=	V3S(a,t)		&
								-dV3I

				V3I(a,hh)	=	V3I(a,t)		&
								+dV3I			&
								-dV3R

				V3R(a,hh)	=	V3R(a,t)		&
								+dV3R
			1978 continue
			!loop through age-groups, oldest until 3 years, people ageing
			do 1927 a=(amax+1),157,-1
				!M(a,hh)	=M(a,hh) 	-M(a,hh)*1/real(tstep)		+M(a-1,hh)*1/real(tstep)
				!S(a,hh)	=S(a,hh) 	-S(a,hh)*1/real(tstep)		+S(a-1,hh)*1/real(tstep)
				M(a,hh)		=	M(a,hh-1)						&
								-agerate_annual*M(a,hh-1)		&
								+agerate_annual*M(a-1,hh-1)		&
								-real(0.002, dp)*M(a,hh-1)

				S(a,hh)		=	S(a,hh)							&
								-agerate_annual*S(a,hh)			&
								+agerate_annual*S(a-1,hh)		&
								+real(0.002, dp)*M(a,hh-1)

				I(a,hh)		=	I(a,hh)							&
								-agerate_annual*I(a,hh)			&
								+agerate_annual*I(a-1,hh)

				R(a,hh)		=	R(a,hh)							&
								-agerate_annual*R(a,hh)			&
								+agerate_annual*R(a-1,hh)

				VS(a,hh)	=	VS(a,hh)						&
								-agerate_annual*VS(a,hh)		&
								+agerate_annual*VS(a-1,hh)

				VI(a,hh)	=	VI(a,hh)						&
								-agerate_annual*VI(a,hh)		&
								+agerate_annual*VI(a-1,hh)

				VR(a,hh)	=	VR(a,hh)						&
								-agerate_annual*VR(a,hh)		&
								+agerate_annual*VR(a-1,hh)

				MCV1_eff(a,hh)	=	MCV1_eff(a,hh)				&
								-agerate_annual*MCV1_eff(a,hh)	&
								+agerate_annual*MCV1_eff(a-1,hh)

				V2S(a,hh)	=	V2S(a,hh)						&
								-agerate_annual*V2S(a,hh)		&
								+agerate_annual*V2S(a-1,hh)

				V2I(a,hh)	=	V2I(a,hh)						&
								-agerate_annual*V2I(a,hh)		&
								+agerate_annual*V2I(a-1,hh)

				V2R(a,hh)	=	V2R(a,hh)						&
								-agerate_annual*V2R(a,hh)		&
								+agerate_annual*V2R(a-1,hh)

				V3S(a,hh)	=	V3S(a,hh)						&
								-agerate_annual*V3S(a,hh)		&
								+agerate_annual*V3S(a-1,hh)

				V3I(a,hh)	=	V3I(a,hh)						&
								-agerate_annual*V3I(a,hh)		&
								+agerate_annual*V3I(a-1,hh)

				V3R(a,hh)	=	V3R(a,hh)						&
								-agerate_annual*V3R(a,hh)		&
								+agerate_annual*V3R(a-1,hh)
			1927 continue
			!old: loop through age-groups, 156 weeks until 105 weeks (3yo)
			!do 19167 a =(156),105,-1
			!update: loop through age-groups, 156 weeks until 73 weeks
			do 19167 a =(156),73,-1
				!M(a,hh)		=M(a,hh)	-M(a,hh)*52/real(tstep) 	+M(a-1,hh)*52/real(tstep)
				!S(a,hh)		=S(a,hh) 	-S(a,hh)*52/real(tstep) 	+(1-cov(a-1))*S(a-1,hh)*52/real(tstep)

				M(a,hh)		=	M(a,hh-1)													&
								-agerate_week*M(a,hh-1)										&
								+agerate_week*M(a-1,hh-1)									&
								-real(0.002, dp)*M(a,hh-1)

				S(a,hh)		=	S(a,hh)														&
								-agerate_week*S(a,hh)										&
								+agerate_week*(real(1,dp)-cov(a-1))*S(a-1,hh)				&
								+real(0.002, dp)*M(a,hh-1)

				I(a,hh)		=	I(a,hh)														&
								-agerate_week*I(a,hh)										&
								+agerate_week*(real(1,dp)-cov(a-1))*I(a-1,hh)

				R(a,hh)		=	R(a,hh)														&
								-agerate_week*R(a,hh)										&
								+agerate_week*(real(1,dp)-cov(a-1))*R(a-1,hh)

				VS(a,hh)	=	VS(a,hh)													&
								-agerate_week*VS(a,hh)										&
								+agerate_week*VS(a-1,hh)									&
								+agerate_week*(real(1,dp)-vaceffByAge(ve_a,ve_b,a-1))*cov(a-1)*S(a-1,hh)

				VI(a,hh)	=	VI(a,hh)													&
								-agerate_week*VI(a,hh)										&
								+agerate_week*VI(a-1,hh)									&
								+agerate_week*cov(a-1)*I(a-1,hh)

				VR(a,hh)	=	VR(a,hh)													&
								-agerate_week*VR(a,hh)										&
								+agerate_week*VR(a-1,hh)									&
								+agerate_week*vaceffByAge(ve_a,ve_b,a-1)*cov(a-1)*S(a-1,hh)	&
								+agerate_week*cov(a-1)*R(a-1,hh)

				MCV1_eff(a,hh)	=	MCV1_eff(a,hh)											&
								-agerate_week * MCV1_eff(a,hh)								&
								+agerate_week * MCV1_eff(a-1,hh)							&
								+agerate_week * vaceffByAge(ve_a,ve_b,a-1) * cov(a-1) * S(a-1,hh)

				V2S(a,hh)	=	V2S(a,hh)													&
								-agerate_week*V2S(a,hh)										&
								+agerate_week*V2S(a-1,hh)

				V2I(a,hh)	=	V2I(a,hh)													&
								-agerate_week*V2I(a,hh)										&
								+agerate_week*V2I(a-1,hh)

				V2R(a,hh)	=	V2R(a,hh)													&
								-agerate_week*V2R(a,hh)										&
								+agerate_week*V2R(a-1,hh)

				V3S(a,hh)	=	V3S(a,hh)													&
								-agerate_week*V3S(a,hh)										&
								+agerate_week*V3S(a-1,hh)


				V3I(a,hh)	=	V3I(a,hh)													&
								-agerate_week*V3I(a,hh)										&
								+agerate_week*V3I(a-1,hh)

				V3R(a,hh)	=	V3R(a,hh)													&
								-agerate_week*V3R(a,hh)										&
								+agerate_week*V3R(a-1,hh)
			19167 continue
			!old: those that are 104weeks (exactly 2yo)
			!a=104
			!update: those that are 72 weeks (15-18 months)
			a=72			
			!calculate new vaccinations
			!M(a,hh)		=M(a,hh)	-M(a,hh)*52/real(tstep) 	+M(a-1,hh)*52/real(tstep)
			!S(a,hh)		=S(a,hh)	-S(a,hh)*52/real(tstep)		+(1-cov2)*(1-cov(a-1))*S(a-1,hh)*52/real(tstep)
			M(a,hh)		=	M(a,hh-1)																	&
							-agerate_week*M(a,hh-1)														&
							+agerate_week*M(a-1,hh-1)													&
							-real(0.002, dp)*M(a,hh-1)

			S(a,hh)		=	S(a,hh)																		&
							-agerate_week*S(a,hh)														&
							+agerate_week*(real(1, dp)-cov2)*(real(1, dp)-cov(a-1))*S(a-1,hh)			&
							+real(0.002, dp)*M(a,hh-1)

			I(a,hh)		=	I(a,hh)																		&
							-agerate_week*I(a,hh)														&
							+agerate_week*(real(1, dp)-cov2)*(real(1, dp)-cov(a-1))*I(a-1,hh)

			R(a,hh)		=	R(a,hh)																		&
							-agerate_week*R(a,hh)														&
							+agerate_week*(real(1, dp)-cov2)*(real(1, dp)-cov(a-1))*R(a-1,hh)

			VS(a,hh)	=	VS(a,hh)																	&
							-agerate_week*VS(a,hh)														&
							+agerate_week*(real(1, dp)-cov2)*VS(a-1,hh)									&
							+agerate_week*(real(1, dp)-vaceffByAge(ve_a,ve_b,a-1))*(cov(a-1)+cov2-cov(a-1)*cov2)*S(a-1,hh)

			VI(a,hh)	=	VI(a,hh)																	&
							-agerate_week*VI(a,hh)														&
							+agerate_week*(real(1, dp)-cov2)*VI(a-1,hh)									&
							+agerate_week*(cov(a-1)+cov2-cov(a-1)*cov2)*I(a-1,hh)

			VR(a,hh)	=	VR(a,hh)																	&
							-agerate_week*VR(a,hh)														&
							+agerate_week*(real(1, dp)-cov2)*VR(a-1,hh)									&
							+agerate_week*vaceffByAge(ve_a,ve_b,a-1)*(cov(a-1)+cov2-cov(a-1)*cov2)*S(a-1,hh)					&
							+agerate_week*(cov(a-1)+cov2-cov(a-1)*cov2)*R(a-1,hh)

			MCV1_eff(a,hh)	=	MCV1_eff(a,hh)															&
							-agerate_week*MCV1_eff(a,hh)												&
							+agerate_week*(real(1, dp)-cov2)*MCV1_eff(a-1,hh)							&
							+agerate_week*vaceffByAge(ve_a,ve_b,a-1)*(cov(a-1)+cov2-cov(a-1)*cov2)*S(a-1,hh)

			vac_numerator = MCV1_eff(a,hh)
			vac_denominator = (VS(a,hh) + VI(a,hh) + VR(a,hh))
			if(vac_denominator.gt.real(0, dp)) then
				take_inc = (take3 - (vac_numerator/vac_denominator))/(real(1, dp) - (vac_numerator/vac_denominator))
			else
				take_inc = (take3 - real(0, dp))/(real(1, dp) - real(0, dp))
			endif

			V2S(a,hh)	=	V2S(a,hh)																	&
							-agerate_week*V2S(a,hh)														&
							+agerate_week*(real(1, dp)-cov2)*V2S(a-1,hh)								&
							!vaccine efficacy no longer increases after 2nd dose
							+agerate_week*(real(1, dp)-take_inc)*cov2*VS(a-1,hh)

			V2I(a,hh)	=	V2I(a,hh)																	&
							-agerate_week*V2I(a,hh)														&
							+agerate_week*(real(1, dp)-cov2)*V2I(a-1,hh)								&
							+agerate_week*cov2*VI(a-1,hh)

			V2R(a,hh)	=	V2R(a,hh)																	&
							-agerate_week*V2R(a,hh)														&
							+agerate_week*(real(1, dp)-cov2)*V2R(a-1,hh)								&
							+agerate_week*take_inc*cov2*VS(a-1,hh)										&
							+agerate_week*cov2*VR(a-1,hh)

			take_inc 	=	real(0, dp)

			V3S(a,hh)	=	V3S(a,hh)																	&
							-agerate_week*V3S(a,hh)														&
							+agerate_week*(real(1, dp)-take_inc)*cov2*V2S(a-1,hh)						&
							+agerate_week*V3S(a-1,hh)

			V3I(a,hh)	=	V3I(a,hh)																	&
							-agerate_week*V3I(a,hh)														&
							+agerate_week*cov2*V2I(a-1,hh)												&
							+agerate_week*V3I(a-1,hh)

			V3R(a,hh)	=	V3R(a,hh)																	&
							-agerate_week*V3R(a,hh)														&
							+agerate_week*take_inc*cov2*V2S(a-1,hh)										&
							+agerate_week*cov2*V2R(a-1,hh)												&
							+agerate_week*V3R(a-1,hh)
			!old: loop through age-groups 103-55 weeks old
			!do 1967 a=103,55,-1
			!loop through age-groups 71-55 weeks old
			do 1967 a=71,55,-1
				!M(a,hh)		=M(a,hh)	-M(a,hh)*52/real(tstep) 		+(1-cov(a-1))*M(a-1,hh)*52/real(tstep)
				!S(a,hh)		=S(a,hh) 	-S(a,hh)*52/real(tstep) 		+(1-cov(a-1))*S(a-1,hh)*52/real(tstep)
				M(a,hh)		=	M(a,hh-1)                                                        	&
								-agerate_week*M(a,hh-1)                                          	&
								+agerate_week*(real(1, dp)-cov(a-1))*M(a-1,hh-1)                 	&
								-real(0.002, dp)*M(a,hh-1)

				S(a,hh)		=	S(a,hh)                                                          	&
								-agerate_week*S(a,hh)                                            	&
								+agerate_week*(real(1, dp)-cov(a-1))*S(a-1,hh)                   	&
								+real(0.002, dp)*M(a,hh-1)

				I(a,hh)		=	I(a,hh)                                                          	&
								-agerate_week*I(a,hh)                                            	&
								+agerate_week*(real(1, dp)-cov(a-1))*I(a-1,hh)

				R(a,hh)		=	R(a,hh)                                                          	&
								-agerate_week*R(a,hh)                                            	&
								+agerate_week*(real(1, dp)-cov(a-1))*R(a-1,hh)

				VS(a,hh)	=	VS(a,hh)                                                         	&
								-agerate_week*VS(a,hh) 												&
								+agerate_week*VS(a-1,hh)                 							&
								+agerate_week*(real(1, dp)-vaceffByAge(ve_a,ve_b,a-1))*cov(a-1)*(M(a-1,hh-1)+S(a-1,hh))

				VI(a,hh)	=	VI(a,hh)                                            				&
								-agerate_week*VI(a,hh)												&
								+agerate_week*VI(a-1,hh)    										&
								+agerate_week*cov(a-1)*I(a-1,hh)

				VR(a,hh)	=	VR(a,hh)                                            				&
								-agerate_week*VR(a,hh)												&
								+agerate_week*VR(a-1,hh)    										&
								+agerate_week*cov(a-1)*R(a-1,hh)                    				&
								+agerate_week*vaceffByAge(ve_a,ve_b,a-1)*cov(a-1)*(M(a-1,hh-1)+S(a-1,hh))

				MCV1_eff(a,hh)	=	MCV1_eff(a,hh)                                            			&
								-agerate_week*MCV1_eff(a,hh)										&
								+agerate_week*MCV1_eff(a-1,hh)    									&
								+agerate_week*vaceffByAge(ve_a,ve_b,a-1)*cov(a-1)*(M(a-1,hh-1)+S(a-1,hh))

				V2S(a,hh)	=	V2S(a,hh)                                           				&
								-agerate_week*V2S(a,hh)												&
								+agerate_week*V2S(a-1,hh)

				V2I(a,hh)	=	V2I(a,hh)                                           				&
								-agerate_week*V2I(a,hh)												&
								+agerate_week*V2I(a-1,hh)

				V2R(a,hh)	=	V2R(a,hh)                                           				&
								-agerate_week*V2R(a,hh)												&
								+agerate_week*V2R(a-1,hh)

				V3S(a,hh)	=	V3S(a,hh)                                           				&
								-agerate_week*V3S(a,hh)												&
								+agerate_week*V3S(a-1,hh)

				V3I(a,hh)	=	V3I(a,hh)                                           				&
								-agerate_week*V3I(a,hh)												&
								+agerate_week*V3I(a-1,hh)

				V3R(a,hh)	=	V3R(a,hh)                                           				&
								-agerate_week*V3R(a,hh)												&
								+agerate_week*V3R(a-1,hh)
			1967 continue
			!loop through ages 54-2 weeks old (1 year and 2 weeks - 2 weeks)
			do 1957 a=54,2,-1
				!M(a,hh)		=M(a,hh) 	-M(a,hh)*52/real(tstep) 			+(1-cov(a-1))*M(a-1,hh)*52/real(tstep)
				!S(a,hh)		=S(a,hh) 	-S(a,hh)*52/real(tstep) 			+(1-cov(a-1))*S(a-1,hh)*52/real(tstep)

				M(a,hh)		=	M(a,hh-1)													&
								-agerate_week*M(a,hh-1)										&
								+agerate_week*(real(1, dp)-cov(a-1))*M(a-1,hh-1)			&
								-real(0.002, dp)*M(a,hh-1)

				S(a,hh)		=	S(a,hh)														&
								-agerate_week*S(a,hh)										&
								+agerate_week*(real(1, dp)-cov(a-1))*S(a-1,hh)				&
								+real(0.002, dp)*M(a,hh-1)

				I(a,hh)		=	I(a,hh)														&
								-agerate_week*I(a,hh)										&
								+agerate_week*(real(1, dp)-cov(a-1))*I(a-1,hh)

				R(a,hh)		=	R(a,hh)														&
								-agerate_week*R(a,hh)										&
								+agerate_week*(real(1, dp)-cov(a-1))*R(a-1,hh)

				VS(a,hh)	=	VS(a,hh)													&
								-agerate_week*VS(a,hh)										&
								+agerate_week*VS(a-1,hh)									&
								+agerate_week*(1-vaceffByAge(ve_a,ve_b,a-1))*cov(a-1)*(M(a-1,hh-1)+S(a-1,hh))

				VI(a,hh)	=	VI(a,hh)													&
								-agerate_week*VI(a,hh)										&
								+agerate_week*VI(a-1,hh)									&
								+agerate_week*cov(a-1)*I(a-1,hh)

				VR(a,hh)	=	VR(a,hh)													&
								-agerate_week*VR(a,hh)										&
								+agerate_week*VR(a-1,hh)									&
								+agerate_week*vaceffByAge(ve_a,ve_b,a-1)*cov(a-1)*(M(a-1,hh-1)+S(a-1,hh))		&
								+agerate_week*cov(a-1)*R(a-1,hh)

				MCV1_eff(a,hh)	=	MCV1_eff(a,hh)											&
								-agerate_week*MCV1_eff(a,hh)								&
								+agerate_week*MCV1_eff(a-1,hh)								&
								+agerate_week*vaceffByAge(ve_a,ve_b,a-1)*cov(a-1)*(M(a-1,hh-1)+S(a-1,hh))

				V2S(a,hh)	=	V2S(a,hh)													&
								-agerate_week*V2S(a,hh)										&
								+agerate_week*V2S(a-1,hh)

				V2I(a,hh)	=	V2I(a,hh)													&
								-agerate_week*V2I(a,hh)										&
								+agerate_week*V2I(a-1,hh)

				V2R(a,hh)	=	V2R(a,hh)													&
								-agerate_week*V2R(a,hh)										&
								+agerate_week*V2R(a-1,hh)

				V3S(a,hh)	=	V3S(a,hh)													&
								-agerate_week*V3S(a,hh)										&
								+agerate_week*V3S(a-1,hh)

				V3I(a,hh)	=	V3I(a,hh)													&
								-agerate_week*V3I(a,hh)										&
								+agerate_week*V3I(a-1,hh)

				V3R(a,hh)	=	V3R(a,hh)													&
								-agerate_week*V3R(a,hh)										&
								+agerate_week*V3R(a-1,hh)
			1957 continue
			!those that are 1 week old
			M(1,hh) 	=	M(1,hh-1)																		&
							-agerate_week*M(1,hh-1)															&
							+agerate_week*sum(																&
								(																			&
									R(172:189,hh)+VR(172:189,hh)+V2R(172:189,hh)+V3R(172:189,hh)			&
								)*pop(172:189)																&
							)/sum(																			&
								(																			&
									R(172:189,hh)+VR(172:189,hh)+V2R(172:189,hh)+V3R(172:189,hh)			&
									+S(172:189,hh)+VS(172:189,hh)+V2S(172:189,hh)+V3S(172:189,hh)			&
								)*pop(172:189) 																&
							)

			S(1,hh) 	=	S(1,hh)																			&
							-agerate_week*S(1,hh)															&
							+agerate_week*sum(																&
								(																			&
									S(172:189,hh)+VS(172:189,hh)+V2S(172:189,hh)+V3S(172:189,hh)			&
								)*pop(172:189)																&
							)/sum(																			&
								(																			&
									R(172:189,hh)+VR(172:189,hh)+V2R(172:189,hh)+V3R(172:189,hh)			&
									+S(172:189,hh)+VS(172:189,hh)+V2S(172:189,hh)+V3S(172:189,hh)			&
								)*pop(172:189) 																&
							)

			I(1,hh) 	=	I(1,hh)																			&
							-agerate_week*I(1,hh)

			R(1,hh) 	=	R(1,hh)																			&
							-agerate_week*R(1,hh)

			VS(1,hh) 	=	VS(1,hh)																		&
							-agerate_week*VS(1,hh)

			VI(1,hh) 	=	VI(1,hh)																		&
							-agerate_week*VI(1,hh)

			VR(1,hh) 	=	VR(1,hh)																		&
							-agerate_week*VR(1,hh)

			MCV1_eff(1,hh)	=	MCV1_eff(1,hh)																&
							-agerate_week*MCV1_eff(1,hh)

			V2S(1,hh) 	=	V2S(1,hh)																		&
							-agerate_week*V2S(1,hh)

			V2I(1,hh) 	=	V2I(1,hh)																		&
							-agerate_week*V2I(1,hh)

			V2R(1,hh) 	=	V2R(1,hh)																		&
							-agerate_week*V2R(1,hh)

			V3S(1,hh) 	=	V3S(1,hh)																		&
							-agerate_week*V3S(1,hh)

			V3I(1,hh) 	=	V3I(1,hh)																		&
							-agerate_week*V3I(1,hh)

			V3R(1,hh) 	=	V3R(1,hh)																		&
							-agerate_week*V3R(1,hh)

			!perform SIA if applicable at this time-point
			if(any(t.eq.t_sia))then
				!increase number of sia (iteration)
				sia_n=sia_n+1
				!loop through ages included in the SIA campaign
				do 1965 a=a0(sia_n),a1(sia_n)
					tempUnvacS = S(:,hh)
					tempUnvacI = I(:,hh)
					tempUnvacR = R(:,hh)
					! call function SIAcov(a, b, m1cov, scov, agefrom, ageto, pop, susceptibles)
					call SIAcov(SIA_a, SIA_b, sum(cov), sia(sia_n), a0(sia_n), a1(sia_n), pop, tempUnvacS, tempUnvacI, tempUnvacR, sia1, sia2)
          ! sia1 is sia coverage in unvaccinated (zero doses)
					!sia1 = SIAoutput(1)
					! sia2 is coverage in everyone else
					!sia2 = SIAoutput(2)
					!mutations between compartments
					if(a.gt.0) then
						!unvaccinated
						dS_VS		=	S(a,hh)*( real(1, dp) -vaceffByAge(ve_a,ve_b,a) )*sia1
						dS_VR		=	S(a,hh) * vaceffByAge(ve_a,ve_b,a) * sia1
						dI_VI		=	I(a,hh)*sia1
						dR_VR		=	R(a,hh)*sia1

						!vaccinated < age 1
						!estimate the incremental efficacy derived from MCV2 to achieve an overall efficacy of take3. Estimate current MCV1 efficacy from those for who the vaccine worked in all who received MCV1 only
						vac_numerator = MCV1_eff(a,hh)
						vac_denominator = (VS(a,hh) + VI(a,hh) + VR(a,hh))
						if(vac_denominator.gt.real(0, dp)) then
							take_inc = (take3 - (vac_numerator/vac_denominator))/(real(1, dp) - (vac_numerator/vac_denominator))
						else
							take_inc = (take3 - real(0, dp))/(real(1, dp) - real(0, dp))
						endif
						dVS_V2S		=	VS(a,hh)*(real(1, dp)-take_inc)*sia2
						dVS_V2R		=	VS(a,hh)*take_inc*sia2
						dVI_V2I		=	VI(a,hh)*sia2
						dVR_V2R		=	VR(a,hh)*sia2

						!vaccinated > age1
						take_inc		=	real(0, dp) !vaccine efficacy for 3rd dose is 0 (overall VE does no longer increase)
						dV2S_V3S	=	V2S(a,hh)*(real(1, dp)-take_inc)*sia2
						dV2S_V3R	=	V2S(a,hh)*take_inc*sia2
						dV2I_V3I	=	V2I(a,hh)*sia2
						dV2R_V3R	=	V2R(a,hh)*sia2

						!update compartments
						S(a,hh)		=	S(a,hh)				&
										-dS_VS				&
										-dS_VR

						I(a,hh)		=	I(a,hh)				&
										-dI_VI

						R(a,hh)		=	R(a,hh)				&
										-dR_VR

						VS(a,hh)	=	VS(a,hh)			&
										+dS_VS				&
										-dVS_V2S				&
										-dVS_V2R

						VI(a,hh)	=	VI(a,hh)			&
										+dI_VI				&
										-dVI_V2I

						VR(a,hh)	=	VR(a,hh)			&
										+dR_VR				&
										+dS_VR				&
										-dVR_V2R

						MCV1_eff(a,hh)	=	MCV1_eff(a,hh)	&
											+dS_VR			&
											-MCV1_eff(a,hh)*sia2

						V2S(a,hh)	=	V2S(a,hh)			&
										-dV2S_V3S			&
										-dV2S_V3R			&
										+dVS_V2S

						V2I(a,hh)	=	V2I(a,hh)			&
										-dV2I_V3I			&
										+dVI_V2I

						V2R(a,hh)	=	V2R(a,hh)			&
										-dV2R_V3R			&
										+dVS_V2R 			&
										+dVR_V2R

						V3S(a,hh)	=	V3S(a,hh)			&
										+dV2S_V3S

						V3I(a,hh)	=	V3I(a,hh)			&
										+dV2I_V3I

						V3R(a,hh)	=	V3R(a,hh)			&
										+dV2S_V3R			&
										+dV2R_V3R
					endif
				1965 continue
			endif
			do 196 a=1,(amax)
				N(a,hh)=pop(a)
			196 continue

		1980 continue
	!go to next year if tmax is not yet reached
	if(t.lt.tmax) then
		goto 1234
	endif

	if(psa_run.gt.0) then
		write(outfile,'(A2,A1,A2,A1,A7,A1,A10,A1,A20,A1,A6,A1,A3,I3.3,A1,A3,A32)') '..',sep,'..',sep,'outcome',sep,get_scenarioname,sep,get_foldername,sep,'output',sep,'run',psa_run,sep,get_country_iso3,'_age_stratified_cases_byyear.txt'
	else
		write(outfile,'(A2,A1,A2,A1,A7,A1,A10,A1,A20,A1,A6,A1,A3,A32)') '..',sep,'..',sep,'outcome',sep,get_scenarioname,sep,get_foldername,sep,'output',sep,get_country_iso3,'_age_stratified_cases_byyear.txt'
	endif
	open(1212,file=outfile,status='NEW')
	if(psa_run.gt.0) then
		write(outfile,'(A2,A1,A2,A1,A7,A1,A10,A1,A20,A1,A6,A1,A3,I3.3,A1,A3,A34)') '..',sep,'..',sep,'outcome',sep,get_scenarioname,sep,get_foldername,sep,'output',sep,'run',psa_run,sep,get_country_iso3,'_age_stratified_popsize_byyear.txt'
	else
		write(outfile,'(A2,A1,A2,A1,A7,A1,A10,A1,A20,A1,A6,A1,A3,A34)') '..',sep,'..',sep,'outcome',sep,get_scenarioname,sep,get_foldername,sep,'output',sep,get_country_iso3,'_age_stratified_popsize_byyear.txt'
	endif
	open(1213,file=outfile,status='NEW')
	!creates a small file
	!write data for each year
	do 1544 x=1,uu-1,1
		tstart = 1+(x-1)*tstep
		tend = x*tstep
		write(1212,1069)	sum(newinfect(1:52,tstart:tend)),				sum(newinfect(53:104,tstart:tend)),			&
							sum(newinfect(105:156,tstart:tend)),			sum(newinfect(157,tstart:tend)),			&
							sum(newinfect(158,tstart:tend)),				sum(newinfect(159,tstart:tend)),			&
							sum(newinfect(160,tstart:tend)),				sum(newinfect(161,tstart:tend)),			&
							sum(newinfect(162,tstart:tend)),				sum(newinfect(163,tstart:tend)),			&
							sum(newinfect(164,tstart:tend)),				sum(newinfect(165,tstart:tend)),			&
							sum(newinfect(166,tstart:tend)),				sum(newinfect(167,tstart:tend)),			&
							sum(newinfect(168,tstart:tend)),				sum(newinfect(169,tstart:tend)),			&
							sum(newinfect(170,tstart:tend)),				sum(newinfect(171,tstart:tend)),			&
							sum(newinfect(172,tstart:tend)),				sum(newinfect(173,tstart:tend)),			&
							sum(newinfect(174,tstart:tend)),				sum(newinfect(175,tstart:tend)),			&
							sum(newinfect(176,tstart:tend)),				sum(newinfect(177,tstart:tend)),			&
							sum(newinfect(178,tstart:tend)),				sum(newinfect(179,tstart:tend)),			&
							sum(newinfect(180,tstart:tend)),				sum(newinfect(181,tstart:tend)),			&
							sum(newinfect(182,tstart:tend)),				sum(newinfect(183,tstart:tend)),			&
							sum(newinfect(184,tstart:tend)),				sum(newinfect(185,tstart:tend)),			&
							sum(newinfect(186,tstart:tend)),				sum(newinfect(187,tstart:tend)),			&
							sum(newinfect(188,tstart:tend)),				sum(newinfect(189,tstart:tend)),			&
							sum(newinfect(190,tstart:tend)),				sum(newinfect(191,tstart:tend)),			&
							sum(newinfect(192,tstart:tend)),				sum(newinfect(193,tstart:tend)),			&
							sum(newinfect(194,tstart:tend)),				sum(newinfect(195,tstart:tend)),			&
							sum(newinfect(196,tstart:tend)),				sum(newinfect(197,tstart:tend)),			&
							sum(newinfect(198,tstart:tend)),				sum(newinfect(199,tstart:tend)),			&
							sum(newinfect(200,tstart:tend)),				sum(newinfect(201,tstart:tend)),			&
							sum(newinfect(202,tstart:tend)),				sum(newinfect(203,tstart:tend)),			&
							sum(newinfect(204,tstart:tend)),				sum(newinfect(205,tstart:tend)),			&
							sum(newinfect(206,tstart:tend)),				sum(newinfect(207,tstart:tend)),			&
							sum(newinfect(208,tstart:tend)),				sum(newinfect(209,tstart:tend)),			&
							sum(newinfect(210,tstart:tend)),				sum(newinfect(211,tstart:tend)),			&
							sum(newinfect(212,tstart:tend)),				sum(newinfect(213,tstart:tend)),			&
							sum(newinfect(214,tstart:tend)),				sum(newinfect(215,tstart:tend)),			&
							sum(newinfect(216,tstart:tend)),				sum(newinfect(217,tstart:tend)),			&
							sum(newinfect(218,tstart:tend)),				sum(newinfect(219,tstart:tend)),			&
							sum(newinfect(220,tstart:tend)),				sum(newinfect(221,tstart:tend)),			&
							sum(newinfect(222,tstart:tend)),				sum(newinfect(223,tstart:tend)),			&
							sum(newinfect(224,tstart:tend)),				sum(newinfect(225,tstart:tend)),			&
							sum(newinfect(226,tstart:tend)),				sum(newinfect(227,tstart:tend)),			&
							sum(newinfect(228,tstart:tend)),				sum(newinfect(229,tstart:tend)),			&
							sum(newinfect(230,tstart:tend)),				sum(newinfect(231,tstart:tend)),			&
							sum(newinfect(232,tstart:tend)),				sum(newinfect(233,tstart:tend)),			&
							sum(newinfect(234,tstart:tend)),				sum(newinfect(235,tstart:tend)),			&
							sum(newinfect(236,tstart:tend)),				sum(newinfect(237,tstart:tend)),			&
							sum(newinfect(238,tstart:tend)),				sum(newinfect(239,tstart:tend)),			&
							sum(newinfect(240,tstart:tend)),				sum(newinfect(241,tstart:tend)),			&
							sum(newinfect(242,tstart:tend)),				sum(newinfect(243,tstart:tend)),			&
							sum(newinfect(244,tstart:tend)),				sum(newinfect(245,tstart:tend)),			&
							sum(newinfect(246,tstart:tend)),				sum(newinfect(247,tstart:tend)),			&
							sum(newinfect(248,tstart:tend)),				sum(newinfect(249,tstart:tend)),			&
							sum(newinfect(250,tstart:tend)),				sum(newinfect(251,tstart:tend)),			&
							sum(newinfect(252,tstart:tend)),				sum(newinfect(253,tstart:tend)),			&
							sum(newinfect(254,tstart:tend))
		!calculate the average population-size in a year
		write(1213,1069)	sum(N(1:52,tstart:tend))/real(tstep, dp),		sum(N(53:104,tstart:tend))/real(tstep, dp),	&
							sum(N(105:156,tstart:tend))/real(tstep, dp),	sum(N(157,tstart:tend))/real(tstep, dp),	&
							sum(N(158,tstart:tend))/real(tstep, dp),		sum(N(159,tstart:tend))/real(tstep, dp),	&
							sum(N(160,tstart:tend))/real(tstep, dp),		sum(N(161,tstart:tend))/real(tstep, dp),	&
							sum(N(162,tstart:tend))/real(tstep, dp),		sum(N(163,tstart:tend))/real(tstep, dp),	&
							sum(N(164,tstart:tend))/real(tstep, dp),		sum(N(165,tstart:tend))/real(tstep, dp),	&
							sum(N(166,tstart:tend))/real(tstep, dp),		sum(N(167,tstart:tend))/real(tstep, dp),	&
							sum(N(168,tstart:tend))/real(tstep, dp),		sum(N(169,tstart:tend))/real(tstep, dp),	&
							sum(N(170,tstart:tend))/real(tstep, dp),		sum(N(171,tstart:tend))/real(tstep, dp),	&
							sum(N(172,tstart:tend))/real(tstep, dp),		sum(N(173,tstart:tend))/real(tstep, dp),	&
							sum(N(174,tstart:tend))/real(tstep, dp),		sum(N(175,tstart:tend))/real(tstep, dp),	&
							sum(N(176,tstart:tend))/real(tstep, dp),		sum(N(177,tstart:tend))/real(tstep, dp),	&
							sum(N(178,tstart:tend))/real(tstep, dp),		sum(N(179,tstart:tend))/real(tstep, dp),	&
							sum(N(180,tstart:tend))/real(tstep, dp),		sum(N(181,tstart:tend))/real(tstep, dp),	&
							sum(N(182,tstart:tend))/real(tstep, dp),		sum(N(183,tstart:tend))/real(tstep, dp),	&
							sum(N(184,tstart:tend))/real(tstep, dp),		sum(N(185,tstart:tend))/real(tstep, dp),	&
							sum(N(186,tstart:tend))/real(tstep, dp),		sum(N(187,tstart:tend))/real(tstep, dp),	&
							sum(N(188,tstart:tend))/real(tstep, dp),		sum(N(189,tstart:tend))/real(tstep, dp),	&
							sum(N(190,tstart:tend))/real(tstep, dp),		sum(N(191,tstart:tend))/real(tstep, dp),	&
							sum(N(192,tstart:tend))/real(tstep, dp),		sum(N(193,tstart:tend))/real(tstep, dp),	&
							sum(N(194,tstart:tend))/real(tstep, dp),		sum(N(195,tstart:tend))/real(tstep, dp),	&
							sum(N(196,tstart:tend))/real(tstep, dp),		sum(N(197,tstart:tend))/real(tstep, dp),	&
							sum(N(198,tstart:tend))/real(tstep, dp),		sum(N(199,tstart:tend))/real(tstep, dp),	&
							sum(N(200,tstart:tend))/real(tstep, dp),		sum(N(201,tstart:tend))/real(tstep, dp),	&
							sum(N(202,tstart:tend))/real(tstep, dp),		sum(N(203,tstart:tend))/real(tstep, dp),	&
							sum(N(204,tstart:tend))/real(tstep, dp),		sum(N(205,tstart:tend))/real(tstep, dp),	&
							sum(N(206,tstart:tend))/real(tstep, dp),		sum(N(207,tstart:tend))/real(tstep, dp),	&
							sum(N(208,tstart:tend))/real(tstep, dp),		sum(N(209,tstart:tend))/real(tstep, dp),	&
							sum(N(210,tstart:tend))/real(tstep, dp),		sum(N(211,tstart:tend))/real(tstep, dp),	&
							sum(N(212,tstart:tend))/real(tstep, dp),		sum(N(213,tstart:tend))/real(tstep, dp),	&
							sum(N(214,tstart:tend))/real(tstep, dp),		sum(N(215,tstart:tend))/real(tstep, dp),	&
							sum(N(216,tstart:tend))/real(tstep, dp),		sum(N(217,tstart:tend))/real(tstep, dp),	&
							sum(N(218,tstart:tend))/real(tstep, dp),		sum(N(219,tstart:tend))/real(tstep, dp),	&
							sum(N(220,tstart:tend))/real(tstep, dp),		sum(N(221,tstart:tend))/real(tstep, dp),	&
							sum(N(222,tstart:tend))/real(tstep, dp),		sum(N(223,tstart:tend))/real(tstep, dp),	&
							sum(N(224,tstart:tend))/real(tstep, dp),		sum(N(225,tstart:tend))/real(tstep, dp),	&
							sum(N(226,tstart:tend))/real(tstep, dp),		sum(N(227,tstart:tend))/real(tstep, dp),	&
							sum(N(228,tstart:tend))/real(tstep, dp),		sum(N(229,tstart:tend))/real(tstep, dp),	&
							sum(N(230,tstart:tend))/real(tstep, dp),		sum(N(231,tstart:tend))/real(tstep, dp),	&
							sum(N(232,tstart:tend))/real(tstep, dp),		sum(N(233,tstart:tend))/real(tstep, dp),	&
							sum(N(234,tstart:tend))/real(tstep, dp),		sum(N(235,tstart:tend))/real(tstep, dp),	&
							sum(N(236,tstart:tend))/real(tstep, dp),		sum(N(237,tstart:tend))/real(tstep, dp),	&
							sum(N(238,tstart:tend))/real(tstep, dp),		sum(N(239,tstart:tend))/real(tstep, dp),	&
							sum(N(240,tstart:tend))/real(tstep, dp),		sum(N(241,tstart:tend))/real(tstep, dp),	&
							sum(N(242,tstart:tend))/real(tstep, dp),		sum(N(243,tstart:tend))/real(tstep, dp),	&
							sum(N(244,tstart:tend))/real(tstep, dp),		sum(N(245,tstart:tend))/real(tstep, dp),	&
							sum(N(246,tstart:tend))/real(tstep, dp),		sum(N(247,tstart:tend))/real(tstep, dp),	&
							sum(N(248,tstart:tend))/real(tstep, dp),		sum(N(249,tstart:tend))/real(tstep, dp),	&
							sum(N(250,tstart:tend))/real(tstep, dp),		sum(N(251,tstart:tend))/real(tstep, dp),	&
							sum(N(252,tstart:tend))/real(tstep, dp),		sum(N(253,tstart:tend))/real(tstep, dp),	&
							sum(N(254,tstart:tend))/real(tstep, dp)
	1544 continue

	if( debug_debug > 1 ) then
		if(psa_run.gt.0) then
			write(outfile,'(A2,A1,A2,A1,A7,A1,A10,A1,A20,A1,A6,A1,A3,I3.3,A1,A3,A15)') '..',sep,'..',sep,'outcome',sep,get_scenarioname,sep,get_foldername,sep,'output',sep,'run',psa_run,sep,get_country_iso3,'_debug_model.txt'
		else
			write(outfile,'(A2,A1,A2,A1,A7,A1,A10,A1,A20,A1,A6,A1,A3,A16)') '..',sep,'..',sep,'outcome',sep,get_scenarioname,sep,get_foldername,sep,'output',sep,get_country_iso3,'_debug_model.txt'
		endif
		open(9916,file=outfile,status='NEW')

		if( debug_timepoints == 0 ) then
			t_increment = tstep - 1
			t_scale = 1
		else
			t_increment = 1
			if(debug_timepoints == 1) then
				t_scale = 4
			else
				t_scale = 1
			endif
		endif
		do 9950 t=1, int(t_cov/t_scale), t_increment
			tstart = t
			if( debug_timepoints == 0 ) then
				tend = tstart+tstep
			else
				tend = tstart
			endif
			do 9951 a=1, amax, 1
				astart = a
				!move to next loop if yearly age-strata are used for first 3 years
				if ( debug_age == 0 .AND. a < 157 ) then
					if( a == 1 .OR. a == 53 .OR. a == 105 ) then
						aend = astart + 51
					else
						cycle
					endif
				elseif ( debug_age == 2) then
					if( a == 1 ) then
						aend = amax
					else
						cycle
					endif
				else
					aend = astart
				endif
				if( debug_compartments == 0 ) then
					write(9916,1070) tstart, tend, astart, aend, "cases", sum(newinfect(astart:aend,tstart:tend))
				else
					write(9916,1070) tstart, tend, astart, aend, "M", 	sum(M(astart:aend,tstart:tend))
					write(9916,1070) tstart, tend, astart, aend, "S", 	sum(S(astart:aend,tstart:tend))
					write(9916,1070) tstart, tend, astart, aend, "I", 	sum(I(astart:aend,tstart:tend))
					write(9916,1070) tstart, tend, astart, aend, "R", 	sum(R(astart:aend,tstart:tend))
					write(9916,1070) tstart, tend, astart, aend, "VS", 	sum(VS(astart:aend,tstart:tend))
					write(9916,1070) tstart, tend, astart, aend, "VI", 	sum(VI(astart:aend,tstart:tend))
					write(9916,1070) tstart, tend, astart, aend, "VR", 	sum(VR(astart:aend,tstart:tend))
					write(9916,1070) tstart, tend, astart, aend, "MCV1_eff", sum(MCV1_eff(astart:aend,tstart:tend))
					write(9916,1070) tstart, tend, astart, aend, "V2S",	sum(V2S(astart:aend,tstart:tend))
					write(9916,1070) tstart, tend, astart, aend, "V2I",	sum(V2I(astart:aend,tstart:tend))
					write(9916,1070) tstart, tend, astart, aend, "V2R",	sum(V2R(astart:aend,tstart:tend))
					write(9916,1070) tstart, tend, astart, aend, "V3S",	sum(V3S(astart:aend,tstart:tend))
					write(9916,1070) tstart, tend, astart, aend, "V3I",	sum(V3I(astart:aend,tstart:tend))
					write(9916,1070) tstart, tend, astart, aend, "V3R",	sum(V3R(astart:aend,tstart:tend))
				endif
			9951 continue
		9950 continue
		close(9916)
	endif

	1069 format(101F16.5)
	1070 format(4I10,A10,F15.5)

	return


end program measles



! subroutine SIAcov(a, b, m1cov, scov, agefrom, ageto, pop, susceptibles, unvacI, unvacR, sia1, si)
!
! 	integer, parameter :: dp = selected_real_kind(15, 307) !double precision
! 	!real, dimension(2)  :: SIAcov     ! two outputs - coverage in zero reach, and coverage in everyone else
! 	real(dp)            :: a        ! regression intercept
! 	real(dp)            :: b        ! regression slope
! 	real(dp)            :: mcv1     ! input - mcv1 coverage
! 	real(dp)            :: scov     ! input - overall SIA coverage
! 	integer             :: agefrom
! 	integer             :: ageto
! 	real(dp)            :: tempreach
! 	real(dp) 	          :: doses
! 	real(dp)            :: zeroreach
! 	real(dp), dimension(:) :: pop
! 	real(dp), dimension(:) :: susceptibles
! 	real(dp), dimension(:) :: unvacI
! 	real(dp), dimension(:) :: unvacR
!
! 	! find what the SIA reach should be from the regression
! 	tempreach = a - b*m1cov
! 	! if this number is greater than 1 set it to 0.98
! 	if(tempreach .gt. 1) then
! 		tempreach = real(0.98, dp)
! 	endif
!
! 	doses = scov*sum(pop(agefrom:ageto))
! 	zeroreach = tempreach*(sum(susceptibles(agefrom:ageto)) &
! 							+ sum(unvacI(agefrom:ageto)) &
! 							+ sum(unvacR(agefrom:ageto)))
!
! 	if( doses .le. zeroreach ) then
! 		SIAcov(1) = scov
! 		SIAcov(2) = 0
! 	else
! 		SIAcov(1) = tempreach
! 		SIAcov(2) = (doses-zeroreach)/(sum(pop(agefrom:ageto))&
! 																	 -(sum(susceptibles(agefrom:ageto)) &
! 																			+ sum(unvacI(agefrom:ageto)) &
! 																			+ sum(unvacR(agefrom:ageto))))
! 	endif
!
! 	return
!
! end function SIAcov
