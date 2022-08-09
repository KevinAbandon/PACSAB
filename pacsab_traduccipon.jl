#####Useful thins

#Funcion de lectura de archivos
function read(unit, file, end):
end

#Arreglo de atomos
atoms = []

#Funcion de escritura
write(unit)
end

######



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
File9 = read(unit=9, file=file9)
kk = 0
im = 0
for n in File9:
	c1, j, c2, c3, c4, k, x, y, z = read(9, 1000, end=51)
	atom[n] = c2
	res[n] = c3
	cad[n] = c4
	ind1[n] = k
	if  k < ind1[n-1]
		kk = kk + ind1[n-1]
	end
	if cad[n] != cad[n-1]
		im = im + 1
	end
	ind2[n] = k + kk
	imol[n] = im 
	k1 = ind2[n]
	r[n, 1] = x
	r[n, 2] = y
	r[n, 3] = z
	if atom[n] == "N"
		in[k1] = n
	elseif atom[n] == "H"
		ih[k1] = n
	elseif atom[n] == "CA"
		ica[k1] = n
	elseif atom[n] == "C"
		ico[k1] = n
	end
end
natom = n -1
nres = ind2[natom]
if natom > natmax
	write(6)
end

##### Asigna un tipo a cada atomo

File16 = open(uint=7, file=file16)

for i = 1:natom
	if atom[i] == "N" || atom[i] == "H" || atom[i] == "C" || atom[i] == "O" || atom[i] == "OXT"
		nat[i] = 1
	else
		for i2 in File16
			c2, c3, c4 = read(7, 71)
			if atom[i] == c2 && res[i] == c3
				nat[i] = nat[i]+1
				j = nat[i]
				atp[i, j] = c4
			end
		end
	end
end

for i = 1:natom
	if atom[i] == "N"
		atp[i, 1] = "nh"
	elseif atom[i] == "H"
		atp[i, 1] = "h"
	elseif atom[i] == "C"
		atp[i, 1] == "co"
	else
		atp[i, 1] = "oc"
	end
end

####Carga los parametros de tipo de atomo

File17 = open(unit=7, file=file17)

for i = 1:natom
	for j= 1:nat(i)
		c1,xq,xfree,xvol,xevdw,xrvdw,xrhc,xmassa = read(7)
		if atp[i,j] == c1
			qa[i, j] = xq
			gfreea[i, j] = xfree
			va[i, h] = xvol
			evdwa[i, j] = xevdw
			rvdwa[i, j] = xrvdw
			rhca[i, j] = 0.8 * xrvdw
			xma[i, j] = xmassa
		end
	end
end

####Propiedades de las esferas
xmassa = 0.0
for i = 1:natom
	xm[i] = 0.0
	qq[i] = 0.0
	vol[i] = 0.0
	gfree[i] = 0.0
	evdw[i] = 0.0
	sumrhc = 0.0
	sumrvdw = 0.0
	if nat[i] == 1
		xm[i] = xma[i, 1]
		qq[i] = qa[i, 1]
		vol[i] = va[i, 1]
		gfree[i] = gfreea[i, 1]
		evdw[i] = evdwa[i, 1]
		rvdw[i] = rvdwa[i, 1]
		rhc[i] = 0.8*rvdw[i]
	else
		for j= 1:nat[i]
			xm[i] = xm[i] + xma[i, j]
			qq[i] = qq[i] + qa[i, j]
			vol[i] = vol[i] + va[i,j]
			gfree[i] = gfree[i] + gfreea[i, j]
			evdw[i] = evdw[i] + evdwa[i, j]
			sumrhc = sumrhc + rhca[i, j]^3
			sumrvdw = sumrvdw + rvdwa[i, j]^3
		end
		rvdw[i] = factr * sumrvdw^0.3333
		rhc[i] = 0.8 * rvdw[i]
	end
	xmassa = xmassa + xm[i]
	v[i, 1] = 0.0
	v[i, 2] = 0.0
	v[i, 3] = 0.0
end

#### Lee la matriz de topologia

for i=1:natom-1
	for j = i+1:natom
		icov[i, j] = 0
		nstep[i, j] = 0
		inter[i, j] = 0
		istruct[i, j] = 0
	end
end

File7 = open(unit=7, file=file7)
for k = 1:File7
	i, j, rij = read(7)
	icov(i,j) = 1
	ibound(k, 1) = i
	ibound(k, 2) = rij
end

npair = k

dmd = open(unit=8, file="dmd.out")
write(8, INPUT)

c1 = "ATOM"

#### Reconoce estructura secundaria y establece puntos de hidrogeno

nbh = 0

for i = 1:natom
	ihb[i] = 0
end

for i = 1:nres-4
	ii = io[i]
	for j = i+4:nres
		if res[ica[j]] != "PRO"
			jj = ih[j]
			if ihb[ii] == 0 &&  ihb[jj] == 0
				n1 = io[i]
				n2 = ih[k]
				rij1 = dbox[n2, n1, 1]
				rij2 = dbox[n2, n1, 2]
				rij3 = dbox[n2, n1, 3]
				roh = sqrt(rij1*rij1+rij2*rij2+rij3*rij3)
				if roh < rohmax && roh > rohmin
					n1 = io[i]
					n2 = in[j]
					rij1 = dbox[n2, n1, 1]
					rij2 = dbox[n2, n1, 2]
					rij3 = dbox[n2, n1, 3]
					rno = sqrt(rij1*rij1+rij2*rij2+rij3*rij3)
					if rno < rnomax && rno > rnomin
						n1 = io[i]
						n2 = in[j]
						rij1 = dbox[n2, n1, 1]
						rij2 = dbox[n2, n1, 2]
						rij3 = dbox[n2, n1, 3]
						rch = sqrt(rij1*rij1+rij2*rij2+rij3*rij3)
						if rch < rchmax && rch > rchmin
							nhb = nhb +1 
							n1 = io[i]
							n2 = ih[j]
							ihb[n1] = 1
							ihb[n2] = 1
							write(6, "HBOND", atom(ii), res(ii), ind2(ii), atom(jj), res(jj), ind2(jj))
							write(8, "HBOND", atom(ii), res(ii), ind2(ii), atom(jj), res(jj), ind2(jj))
						end
					end
				end
			end
		end
	end
end

potencial(natom)

####Lista de solapamientos posibles

for i=1:natom-1
	ishk[i] = 0
	ipot[i] = 0
	for j=i+1:natom
		if icov[i,j] == 0
			rij1 = dbox[i, j, 1]
			rij2 = dbox[i, j, 2]
			rij3 = dbox[i, j, 3]
			rmod2 = rij1*rij1+rij2*rij2+rij3*rij3
			if rmod2 < rshake2
				ishk[i] = ishk[i]+1
				k = ishk[i]
				nshk[i, k] = j
			end
			if istruct[i, j] != 1
				if rmod2 < rpot2
					ipot[i] = ipot[i] +1
					k = ipot[i]
					npot[i, k] = j
				end
			end
		end
	end
end
if isolv != 0
	enchufa(natom, dcut)
end

#### Asigna la región donde se encuentra la interacción entre dos partículas

for i=1:natom-1
	for j=i+1:natom
		if icov[i, j] == 0 
			rij1 = dbox[j, i, 1]
			rij2 = dbox[j, i, 2]
			rij3 = dbox[j, i, 3]
			rmod2 = rij1*rij1+rij2*rij2+rij3+rij3
			rij= sqrt(rmod2)
			k=1
			while rij > rstep[i,j, k] && k < nstep[i, j]
				k = k+1
			end
			ireg[i, j] = k
		end
	end
end

#### Suma la energia potencial de la confromacion inicial

epot0 = 0.0
epotmol0 = 0.0
epothb0 = 0.0
epothmol0 = 0.0
for i=1:natom-1
	for j=i+1_natom
		rij1 = dbox[i, j, 1]
		rij2 = dbox[i, j, 2]
		rij3 = dbox[i, j, 3]
		rmod2=rij1*rij1+rij2*rij2+rij3*rij3
		dist = sqrt(rmod2)
		if inter[i, j] == 1
			k = nstep[i, j]
			while dist < rstep[i, j, k] && k > 0
				epot0 = epot0 - estep[i, j, k]
				if imol[i] != imol[j]
					epotmol0 = epotmol0 - estep[i, j, k]
				end
				if istruct[i, j] == 1
					epothb0 = epothb0 - estep[i, j, k]
				if imol[i] != imol[j] 
					epothmol0 = epothmol0 - estep[i, j, k]
				end
			end
		end
	end
end

#### Asigna velocidades aleatorias

for j=1:3
	vcm[j] = 0.0
	for i=1:natom
		random_number(fi)
		v[i, j] = fi
		vcm[j] = vcm[j] + xm[i]*v[i, j]
	end
	vcm[j] = vcm[j]/xmassa
end

#### Ajusta la energia cinetica a la temperatura requerida

ekin = 0.0
for j=1:3
	for i=1:natom
		v[i, j] = v[i, j] - vcm[j]
		ekin = ekin + 0.5 * xm[i] *(v[i, j]* a)^2
	end
end
sto = 1.5*8.314*natom*temp/ekin
ekin0 = 0
for j = 1:3
	for i = 1: natom
		v[i, j] = v[i, j] * sqrt(sto)
		ekin0 = ekin0 + 0.5 * xm[i] * (v[i, j] * a)^2
	end
end
ekin0 = ekin0/facte
etot0 = epot0 + ekin0

#### Busca el CM

for j = 1:3
	rcm[j] = 0.0
	for i=1:natom
		rcm[j] = rcm[j] + xm[i]*r[i, j]
	end
	rcm[j] = rcm[j] / xmassa
end

ibloc = 0

#### Escribe las coordenadas en el SRCM

#### Lista de solapamientos posibles
for i=1:natom-1 
	ishk[i] = 0
	ipot[i] = 0
	for j=i+1:natom
		if icov[i, j] == 0
			rij1 = dbox[i, j, 1]
			rij2 = dbox[i, j, 2]
			rij3 = dbox[i, j, 3]
			rmod2=rij1*rij1+rij2*rij2+rij3*rij3
			if rmod2 < rshake2
				ishk[i] = ishk[i] +1 
				k = ishk[i]
				nshk[i, k] = j
			end
			if istruct[i, j] != 1
				if rmod2 < rpot2
					ipot[i] = ipot[i] +1
					k = ipot[i]
					npot[i, k] = j
				end
			end
		end
	end
end

mem1 = 0
mem2 = 0

iev = 0
ierr = 0
ierr2 = 0

temps0 = temps

for i=1:natom
	for j=1:3
		r[i, j] = rant[i, j]
	end
end

tacene = 0.0

while tacene < tene
	tacact = 0.0
	icint = 0

	dmdshake(natom, npair, mem1, mem2)

	for i = 1:natom-1
		for j = i+1:natom
			if istruct[i, j] == 0
				inter[i, j] = 0
			end
		end
	end
end
if isolv != 0
	enchufa(natom, dcut)
end





