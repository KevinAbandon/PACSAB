




#default values
file9 = "nativain.pdb"
file10 = "res"
file11 = "distancia.dat"
file12 = "energia.dat"
file15 = "res"
file19 = "output.pdb"
file20 = "snapcg.pdb"
file21 = "snapca.pdb"
file7 = "topologia.dat"
file16 = "atomtypes.dat"
file17 = "potentials.dat"
pi= atan(1)*4
a = 1*(10^-10)
tmin = 1 * (10^-30)
dijmin = 1 *(10^-4)
kkk=2381
temp = 300.0
nbloc = 10000
iwr = 0
iprint = 1
tact = 2 * (10^-14)
tene = 1*(10^-12)
tsnap=1 *(10^-11)
rcutgo = 8.0
sigma  = 0.05
sigmago = 0.1
idab = 0
irig = 0
fterm = 4.0
ebond = 1000.0
dstep = 1 *(10^-4)
tpush = 5 * (10^-4)
isolv = 1
iterm = 1
factm = 1.0
icm = 0
rbox = 0.0
rshake = 50.0
rpot = 50.0
dcut = 10.0
icons = 1

#interacciones
eps = 16.5
fvdw = 8.0
fsolv = 15.0
xlamb = 3.5
facthc = 0.8
factr = 0.9

#Puntos de hidrogeno
isec = 0
ehbc = 4.0
ehb = 3
rohmin  = 1.75
roha = 2.15
rohb = 2.34
rohmax = 2.50
rnomin = 2.75
rnoa = 3.1
rnob = 3.2
rnomax = 3.50
rchmin = 2.90
rcha = 3.27
rchb = 3.4
rchmax = 3.75
rsolv = 3.5
asolv = 10
bsolv = 0.5
dwat = 6.0

####Lectura sobre stdinput los parametros INPUT
####Lectura sobre stinput los parametros DISTANCIES

tene = tene > tsnap ? tsnap : tene
tact = tact > tene ? tene: tact
if rbox < 1*(10^-10)
	icm = 1
	rbox = 300.0
end

rbox2 = 0.5 *rbox
rshake2 = rshake*rshake
rpot2 = rpot*rpot

####call random_seed()

#Lectura de coordenadas

