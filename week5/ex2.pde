TITLE 'New Problem'     { the problem identification }
COORDINATES cartesian1  { coordinate system, 1D,2D,3D, etc }
VARIABLES        { system variables }
  psi              { choose your own names }
SELECT         { method controls }
modes = 25
DEFINITIONS    { parameter definitions }
l = 1e-9
hbar = 1.05457e-34
m = 9.1e-31

qe = 1.602e-19
E = lambda/qe
Eth = mode^2*pi^2*hbar^2/(2*m*L^2)
U

psinorm=psi/sqrt(integral(psi^2))

probout = 1-integral(psinorm^2,'well')-integral(psinorm^2,'well2')
probmid = bintegral(psinorm^2,'middle')
probwell1 = integral(psinorm^2,'well')
probwell2 = integral(psinorm^2,'well2')
!INITIAL VALUES
EQUATIONS        { PDE's, one for each variable }
  (-hbar^2/(2*m))*del2(psi) + U*psi = lambda*psi
! CONSTRAINTS    { Integral constraints }
BOUNDARIES       { The domain definition }
region 'outside'
U = 5*qe
start(-L*20) line to (L*20)
  REGION 'well'       { For each material region }
  U= (12*l+x)^2*(9*l+x)^2*1E15
    start(-L*15) line to (-L*6)
    REGION 'well2'
    U = (12*l-x)^2*(9*l-x)^2*1E15
    start(L*6) line to (L*15)
    
    Feature
    Start 'middle' (-L*6) line to (L*6)
! TIME 0 TO 1    { if time dependent }
MONITORS         { show progress }
PLOTS            { save result displays }
  !CONTOUR(u)
  elevation(psinorm^2,U*globalmax(psinorm^2)/globalmax(U),E*qe*globalmax(psinorm^2)/globalmax(U))
  !elevation(U,E*qe)
  !surface(u)
  SUMMARY
	report min(probwell1,probwell2) as 'min:'
END
