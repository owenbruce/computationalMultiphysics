TITLE 'New Problem'     { the problem identification }
COORDINATES cartesian2  { coordinate system, 1D,2D,3D, etc }
VARIABLES        { system variables }
  Voltage              { choose your own names }
  rho_free
SELECT         { method controls }
quadratic
errlim = 1e-12
DEFINITIONS    { parameter definitions }
hb = 1e-3
ht = 1.5e-3
hm = 3e-3
Lx = 12e-3
Lt = 1e-3

Vmax = 3

rhoe = 80e-9
sigmae = 1/rhoe

epr = 1
ep = 8.85e-12*epr

Efield = -grad(Voltage)
Dfield = ep*Efield
J = sigmae*Efield

Qtheory = ep*(Lx-Lt)/hm*Vmax
! INITIAL VALUES
EQUATIONS        { PDE's, one for each variable }
  voltage: div(J)=0 { one possibility }
  rho_free: div(Dfield)=rho_free
! CONSTRAINTS    { Integral constraints }
BOUNDARIES       { The domain definition }
  REGION 'Bot Contact'       { For each material region }
    START(0,0)   { Walk the domain boundary }
    LINE TO (Lx,0)
    line to (Lx,(hb+hm+ht)/2)
    line to (Lx-Lt, (hb + hm + ht)/2)
    line to (Lx-Lt, hb)
    line to (0, hb)
    value(Voltage) = 0
    line to close
   
REGION 'Top Contact'
	Start(0,hb+hm+ht)
    value(Voltage) = Vmax
    line to (0,hb + hm)
    load(Voltage) = 0
    line to (Lx-lt, hb + hm)
    line to (Lx-lt, (hb + hm + ht)/2)
    line to (Lx, (hb + hm + ht)/2)
    line to (Lx, hb + hm + ht)
    line to close
Region 'Top Dielectric'
	sigmae = 1e-12
    epr = 1000
	Start(0,hb+hm)
    line to (Lx-Lt, hb+hm)
    line to (Lx-Lt, (hb + hm + ht)/2)
    line to (0, (hb + hm + ht)/2)
    line to close
    
Region 'Bot Dielectric'
	sigmae = 1e-12
    epr = 1000
	Start(0,hb)
    line to (Lx-Lt, hb)
    line to (Lx-Lt, (hb + hm + ht)/2)
    line to (0, (hb + hm + ht)/2)
    line to close
    
! TIME 0 TO 1    { if time dependent }
MONITORS         { show progress }
PLOTS            { save result displays }
  CONTOUR(Voltage) painted
  contour(rho_free) painted
  report val(Voltage,Lx-Lt/2,(hb + hm + ht)/2)
  vector(J)
  vector(Dfield)
  summary
  report (integral(rho_free,'Top Contact') + integral(rho_free,'Top Dielectric'))
  report (integral(rho_free,'Bot Contact') + integral(rho_free,'Bot Dielectric'))
  report Qtheory
END
