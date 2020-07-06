{ Fill in the following sections (removing comment marks ! if necessary),
  and delete those that are unused.}
TITLE 'New Problem'     { the problem identification }
COORDINATES cartesian2  { coordinate system, 1D,2D,3D, etc }
VARIABLES        { system variables }
  Me              { choose your own names }
! SELECT         { method controls }
DEFINITIONS    { parameter definitions }
Ae = 45
Astar = 1
gamma = 1.4
INITIAL VALUES
	Me = 1
EQUATIONS        { PDE's, one for each variable }
  Ae/Astar = ((gamma+1)/2)^(-(gamma+1)/(2*(gamma-1)))*(1+(gamma-1)/2*Me^2)^(-(gamma+1)/(2*(gamma-1)))/Me
! CONSTRAINTS    { Integral constraints }
BOUNDARIES       { The domain definition }
  REGION 1       { For each material region }
    START(0,0)   { Walk the domain boundary }
    LINE TO (1,0) TO (1,1) TO (0,1) TO CLOSE
! TIME 0 TO 1    { if time dependent }
MONITORS         { show progress }
PLOTS            { save result displays }
Summary
report val(Me,0,0)
END
