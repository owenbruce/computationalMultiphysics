{ Fill in the following sections (removing comment marks ! if necessary),
  and delete those that are unused.}
TITLE 'New Problem'     { the problem identification }
COORDINATES cartesian2  { coordinate system, 1D,2D,3D, etc }
VARIABLES        { system variables }
  Temp(threshold=1e-6)
! SELECT         { method controls }
DEFINITIONS    { parameter definitions }
k
cp
cp2 = cp*0.2
rho 
Tempi
qdot=k*grad(Temp)
AreaIron = area_integral(1,'Iron')
AreaBacon = area_integral(1,'Bacon')
AvgIronTemp = area_integral(Temp,'Iron')/AreaIron
AvgBaconTemp = area_integral(Temp,'Bacon')/AreaBacon
AvgBaconTemp2 = time_integral(Line_Integral(normal(qdot),'contact'))/area_integral(rho*cp,'bacon') + 4
qdotv = 0
INITIAL VALUES
Temp = Tempi
EQUATIONS        { PDE's, one for each variable }
  div(k*grad(Temp)) + qdotv = rho * cp*dt(Temp)
! CONSTRAINTS    { Integral constraints }
BOUNDARIES       { The domain definition }
  REGION 'Iron'       { For each material region }
  k = 80.4
  cp = 450
  rho = 7874
  Tempi = 4
  qdotv = if(AvgIronTemp<180) then 7.15E6 else 0
    START(0,0)
    LINE TO (0.1,0)
    line TO (0.1,0.008)
    line TO (0,0.008)
    line TO CLOSE
    start 'baconline' (0.01,0.008) line to (0.09,0.008)
  REGION 'Bacon'
  k = 0.4
  cp = 4200
  rho = 1000
  Tempi = 4
  	START 'Contact' (0.01,0.008)
    Line to (0.09,0.008)
    LINE TO (0.09,0.028)
    LINE TO (0.01,0.028)
    LINE TO CLOSE
    
TIME 0 TO 10*60  { if time dependent }
MONITORS         { show progress }
PLOTS            { save result displays }
  for time = 10*60
  history(AvgIronTemp) at (0,0)
  history(AvgBaconTemp,AvgBaconTemp2) at (0.05, 0.02)
  contour(Temp) painted
  surface(Temp)
  
SUMMARY
  report val(AvgBaconTemp,0,0)
  report val(AvgBaconTemp2,0,0)
  !report val(line_integral(div(qdot),'Contact'),0,0)
END