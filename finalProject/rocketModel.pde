{ Fill in the following sections (removing comment marks ! if necessary),
  and delete those that are unused.}
TITLE 'New Problem'     { the problem identification }
COORDINATES cartesian3  { coordinate system, 1D,2D,3D, etc }
VARIABLES        { system variables }
  u              { choose your own names }
! SELECT         { method controls }
! DEFINITIONS    { parameter definitions }
! INITIAL VALUES
EQUATIONS        { PDE's, one for each variable }
  div(grad(u))=0 { one possibility }
! CONSTRAINTS    { Integral constraints }
EXTRUSION
	surface "bottom" z = 0
    	layer "body"
    surface "narrow" z = 2
    	layer "cone"
    surface "top" z = 3
BOUNDARIES       { The domain definition }
  REGION 1       { For each material region }
  layer "body"
    START(1,0)   { Walk the domain boundary }
    arc(center=0,0) angle 360
	REGION 2
    layer "cone"
    START (2,0)
    arc(center=0,0) angle 360
! TIME 0 TO 1    { if time dependent }
MONITORS         { show progress }
PLOTS            { save result displays }
	grid(x,y,z)
END
