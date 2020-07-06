TITLE 'Design Project Owen Bruce'     { the problem identification }
COORDINATES cartesian2  { coordinate system, 1D,2D,3D, etc }
VARIABLES        { system variables }
  u (threshold=1e-3)
  v (threshold=1e-3)
SELECT         { method controls }
	ngrid = 1
DEFINITIONS    { parameter definitions }
	MoonRad = 1737.4e3   !Constants
    Mm = 7.346e22	!Moon Mass
    Mr = 104535   !Rocket Mass
    Mp = 11e3	!Payload Mass
    gamma = 1.4
    Rg = 8.314
    p0 = 3*(10e-15)*100e3	!External Pressure
    G = 6.67e-11
 
 
	Astar = 0.8177873014855989 !Rocket Parameters
    Ae = 9.817787301485597
    Pt = 9088936.507427996
    Tt = 4088.9365074279945
    Mfuel =230792.39386
    Me = 4.12776289190494
 
 
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
TIME 0 TO 3600 halt(u-MoonRad>3e3 or u<(MoonRad-1e3))    { if time dependent }
MONITORS         { show progress }
PLOTS            { save result displays }
for time = endtime
  history(u) at (0,0)
  history(Fg) at (0,0)
  history(Ft) at (0,0)
  history(if val(u,0,0)<MoonRad+2.9e3 then 0 else if t<tfuel then tintegral(mdot) else mdot*tfuel)	Export Format '#t#b#1' file = 'FuelConsumed.txt'
Summary
	report(if t<tfuel then tintegral(mdot) else mdot*tfuel) as 'Total mass of fuel consumed'
END

