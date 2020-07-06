{ Fill in the following sections (removing comment marks ! if necessary),
  and delete those that are unused.}
TITLE 'Design Project Owen Bruce'     { the problem identification }
COORDINATES cartesian2  { coordinate system, 1D,2D,3D, etc }
VARIABLES        { system variables }
  u (threshold=1e-3)
  v (threshold=1e-3)
SELECT         { method controls }
	ngrid = 1
DEFINITIONS    { parameter definitions }
	MoonRad = 1737.4e3
    
	Mfuel = 729465
	Mr = 104535
    Mp = 11e3
	
	Astar = 1 !Rocket Parameters
    Pt = 1*100e3
    Tt = 3315.6+273.15
    gamma = 1.4
    Rg = 8.314
    p0 = 3*(10e-15)*100e3
    Ae = 45
    G = 6.67e-11
    
    Mm = 7.346e22
    
    Me = 5.768
    
    
	mdot = Astar*pt/sqrt(Tt)*sqrt(gamma/Rg)*((gamma+1)/2)^(-(gamma+1)/(2*(gamma-1)))   !Thrust Equations
    Te = Tt*(1+(gamma-1)/2*Me^2)^(-1)
    pe = pt*(1+(gamma-1)/2*Me^2)^(-gamma/(gamma-1))
    ve = Me*sqrt(gamma*Rg*Te)
    
    Tfuel = Mfuel/mdot
    
    M = if t<tfuel then Mfuel + Mr + Mp - tfuel*mdot else Mr+Mp
    
	Ft = if t<tfuel then mdot*Ve+(pe-p0)*Ae else (pe-p0)*Ae   !Force Equations
    Fg = -G*M*mm/u^2
    Fnet = Ft + Fg
    a = Fnet/M
INITIAL VALUES
u = MoonRad
EQUATIONS        { PDE's, one for each variable }
	u: dt(u) = v
    v: dt(v) = a
! CONSTRAINTS    { Integral constraints }
BOUNDARIES       { The domain definition }
  REGION 1       { For each material region }
    START(0,0)   { Walk the domain boundary }
    LINE TO (1,0) TO (1,1) TO (0,1) TO CLOSE
TIME 0 TO 1000 halt(u-MoonRad>3e3)    { if time dependent }
MONITORS         { show progress }
PLOTS            { save result displays }
for time = endtime
  history(u) at (0,0)
  history(Fg) at (0,0)
  history(Ft) at (0,0)
Summary
	report(if t<tfuel then tintegral(mdot) else mdot*tfuel) as 'Total mass of fuel consumed'
END
