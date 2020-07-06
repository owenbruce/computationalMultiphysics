TITLE 'New Problem'     { the problem identification }
COORDINATES cartesian1  { coordinate system, 1D,2D,3D, etc }
VARIABLES        { system variables }
  psi              { choose your own names }
SELECT         { method controls }
modes = 3
DEFINITIONS    { parameter definitions }
l = 5e-9
hbar = 1.05457e-34
mp = 1.6726219E-27
mn = 1.6749274E-27
m = 2*(mn+mp)

qe = 1.602e-19
E = lambda/qe
Eth = mode^2*pi^2*hbar^2/(2*m*L^2)/qe
U

psinorm=psi/sqrt(integral(psi^2))
! INITIAL VALUES
EQUATIONS        { PDE's, one for each variable }
  (-hbar^2/(2*m))*del2(psi) + U*psi = lambda*psi
! CONSTRAINTS    { Integral constraints }
BOUNDARIES       { The domain definition }
region 'outside'
	U = 8*qe!2E-5*qe !Example Height
    start(-2*L) line to (2*L)
region 'well'
	U = 0
    start(-L/2) line to (L/2)
! TIME 0 TO 1    { if time dependent }
MONITORS         { show progress }
PLOTS            { save result displays }
  !CONTOUR(u)
  elevation(psinorm^2,U*globalmax(psinorm^2)/globalmax(U),E*qe*globalmax(psinorm^2)/globalmax(U))
  !elevation(U,E*qe)
  !surface(u)
  SUMMARY
  report E
  report Eth
END
