# -*- coding: utf-8 -*

from math import *
from numpy import *

class Solver_DE():
	'''Solver of a simple dark energy model (Carlo Schimdt)
	lamnda-cdm.'''

	def __init__(self, nb_pas_temps, nb_divisions_sphere, rayon_sphere, densite_background_initiale, souspression):
		'''Initialization of variables as soon as the instance
		is created. First we copy the arguments as attributes,
		then we calculate new attributes based on our model.'''
		# Number of time steps
		self.nb_pas_temps = nb_pas_temps
		# Number of inner-shells in our model
		self.nb_divisions_sphere = nb_divisions_sphere
		# Radius of the sphere representing the void
		self.rayon_sphere = rayon_sphere/1 #Mparsecs
		# Initial background density of the universe
		self.densite_background_initiale = densite_background_initiale/1
		# Density difference in the void compared to the background
		self.souspression = souspression
		# Hubble today in (Gigayears)-1
		self.H0 = 1/13.8
		# Time step in Gigayears
		self.pas_temps = 1/(self.H0*self.nb_pas_temps)
		# Percentage of matter in the universe today
		self.Omega_m0 = 0.3
		# Percentage of dark energy in the universe today
		self.Omega_lambda0 = 0.7
		# Cold dark matter ration
		self.Omega_m_0 = 0.3
		# ?
		self.q = 0.01
		# Gravitational constant
		self.G = 1
		# Dark energy state equation
		self.w_de = -1
		# Dark energy density (not accurate)
		self.rho_de = 1
		# Initial radiuses and masses of our void model
		(self.rayons_initiaux, self.masses) = self.calcul_rayons_masses_initiaux()
		# Scale factor function
		self.a = self.calcul_a()
		# Logarithm of the scale factor function
		self.lna = self.calcul_ln_a()
		# Hubble function
		self.H = self.calcul_Hubble()
		# Initial expansion rates of the shells of our void model
		self.vitesses_initiales = self.calcul_vitesses_initiales()
		# Hubble function lna-derivate
		self.dH = self.calcul_derivate_Hubble_lna()
		# Evolution of the ratio of dark-ernegy in the universe
		self.Omega_de_a = self.calcul_Omega_de_a()
		# Evolution of the expansion rate and the radius of the shells
		(self.rayons_evolution, self.vitesses_evolution) = self.calcul_vitesses_rayons_evolution()


	def calcul_a(self):
		'''Calculates the scale factor, a function of
		time. It is represented by a list.
		a(t) = (Omo/Olo)**1/3 * (sinh(t/tl))**2/3
		Here O is Omega and l is lambda'''
		a = []
		t_lambda = 2/(3*self.H0*sqrt(self.Omega_lambda0))
		for i in range(self.nb_pas_temps):
			a = a + [((self.Omega_m0/self.Omega_lambda0)**(1/3))*((sinh(self.pas_temps*(i+1)/t_lambda))**(2/3))]
		return(a)

	def calcul_ln_a(self):
		'''Calculates ln(a), it's necessary for the 
		resolution of the ODE later on.'''
		lna = []
		for i in range(self.nb_pas_temps):
			lna = lna + [log(self.a[i])]
		return(lna)

	def calcul_Hubble(self):
		'''Calculates the Hubble function (time
		dependant). H(t) = d(ln(a))/dt
		I will use finite differences. It's going to
		be ugly.'''
		H = []
		#First value of H will be calculated using only
		#lna[0] and lna[1]. It is unprecise
		#H[0] = (lna[1]-lna[0])/(t[1]-t[0])
		H = H + [(self.lna[1]-self.lna[0])/(self.pas_temps)]
		#For the other values of H that are neither the first
		#nor the last one, we are foinf to be more precise
		#H[i] = (lna[i+1]-lna[i-1])/(t[i+1]-t[i-1])
		for i in range(1,len(self.lna)-1):
			H = H + [(self.lna[i+1]-self.lna[i-1])/(2*self.pas_temps)]
		#For the last value of H, we use whatever is available
		#H[-1] = (lna[-1]-lna[-2])/(t[-1]-t[-2])
		H = H + [(self.lna[-1]-self.lna[-2])/self.pas_temps]
		#Check if H is the same length as lna
		if len(H)==len(self.lna):
			return(H)
		else:
			return("AAAAAAAAAAAAAA")

	def calcul_derivate_Hubble_lna(self):
		'''Calculates the Hubble function's lna-derivate
		using the method of finite differences, it's the same
		as in the calcul_Hubble() function'''
		dH = []
		dH = dH + [(self.H[1]-self.H[0])/(self.lna[1]-self.lna[0])]
		for i in range(1,len(self.H)-1):
			dH = dH + [(self.H[i+1]-self.H[i-1])/(self.lna[i+1]-self.lna[i-1])]
		dH = dH + [(self.H[-1]-self.H[-2])/(self.lna[-1]-self.lna[-2])]
		if len(dH)==len(self.H):
			return(dH)
		else:
			return("AAAAAAAAAAAAAA")

	def calcul_Omega_de_a(self):
		'''Calculates Omegade(a)'''
		Omega_de_a = []
		for i in range(len(self.H)):
			Omega_de_a = Omega_de_a + [(self.H0*self.Omega_lambda0)/self.H[i]]
		return(Omega_de_a)

	def calcul_rayons_masses_initiaux(self):
		'''Calculates the radius and masses of the shells
		representing the void. In order to simulate a 
		smooth transition between the background and the
		void, we modell twice as much shells in total as
		there are inside the void.'''
		rayons_initiaux = []
		masses = []
		for i in range(self.nb_divisions_sphere):
			rayons_initiaux = rayons_initiaux +[(self.rayon_sphere/self.nb_divisions_sphere)*(i+1)]
			masses = masses + [(4/3)*pi*(rayons_initiaux[i]**3)*self.densite_background_initiale*(1+self.souspression)]
		for i in range (self.nb_divisions_sphere):
			rayons_initiaux = rayons_initiaux + [(self.rayon_sphere/self.nb_divisions_sphere)*(self.nb_divisions_sphere+i+1)]
			masses = masses + [(4/3)*pi*(rayons_initiaux[self.nb_divisions_sphere+i]**3)*self.densite_background_initiale*(1+((self.rayon_sphere/rayons_initiaux[i+self.nb_divisions_sphere])**3)*self.souspression)]
		return(rayons_initiaux, masses)

	def calcul_vitesses_initiales(self):
		'''Calculate approximately the initial expansion
		speeds of each shell. I don't know what I'm doing
		so I'll just guess at this point. And I need this
		in order to do the rest.'''
		vitesses_initiales = []
		for i in range(self.nb_divisions_sphere*2):
			vitesses_initiales = vitesses_initiales + [self.rayons_initiaux[i]*self.H[0]]
		return(vitesses_initiales)

	def function_shear_vorticity_DE(self, y, v, t):
		'''Function that gives us the evolution of the
		shell radiusses. Here v = dy/dlna and we write
		dv/dlna'''
		#Convert into uint64 bits because of the size of
		#the numbers. They're too big (numpy function)
		exp_accel = uint64(-(self.dH[t]/self.H[t])*v - 0.5*self.Omega_de_a[t]*(1+3*self.w_de)*y + (4/3)*((self.Omega_m_0/(2*(self.a[t]**3)))*(y**3)-self.q*((self.Omega_m_0**2)/(4*(self.a[t]**6)))*(y**6)))
		return(exp_accel)

	def calcul_vitesses_rayons_evolution(self):
		'''Second-order non-linéar differential equation
		y'' + (H'/H)y' = -0.5*Omegade(a)*(1+3*wde)*y + 
		(4/3)*((Omegam0/2a**3)*y**3-q(Omegam0**2/4a**6)*y**6)
		I proceed as following : v = y'
		and I calculate first v using initial expansion
		rates and initial radiuses and only after can I
		calculate y.
		For the expansion rates and radiuses, I use RK4.'''
		rayons_evolution = [0]*len(self.rayons_initiaux)
		vitesses_evolution = [0]*len(self.vitesses_initiales)
		for i in range(len(rayons_evolution)):
			rayons_evolution[i] = [self.rayons_initiaux[i]]
			vitesses_evolution[i] = [self.vitesses_initiales[i]]
			for t in range(self.nb_pas_temps):
				h = self.lna[1]-self.lna[0] if t==0 else self.lna[t]-self.lna[t-1]
				k1 = h*self.function_shear_vorticity_DE(rayons_evolution[i][t], vitesses_evolution[i][t], t)
				k2 = h*self.function_shear_vorticity_DE(rayons_evolution[i][t], vitesses_evolution[i][t] + 0.5*k1, t)
				k3 = h*self.function_shear_vorticity_DE(rayons_evolution[i][t], vitesses_evolution[i][t] + 0.5*k2, t)
				k4 = h*self.function_shear_vorticity_DE(rayons_evolution[i][t], vitesses_evolution[i][t] + k3, t)
				vitesses_evolution[i] = vitesses_evolution[i] + [vitesses_evolution[i][t] + (1/6)*(k1+2*k2+2*k3+k4)]
				k1 = h*vitesses_evolution[i][-1]
				k2 = h*(vitesses_evolution[i][-1]+0.5*k1)
				k3 = h*(vitesses_evolution[i][-1]+0.5*k2)
				k4 = h*(vitesses_evolution[i][-1]+k3)
				rayons_evolution[i] = rayons_evolution[i] + [rayons_evolution[i][t]+(1/6)*(k1+2*k2+2*k3+k4)]
		return(rayons_evolution, vitesses_evolution)

	def get_rayons_adimension_t(self, pas):
		'''Donne le tableau des rayons (adimensionnés avec le rayon initial de la sphère)
		pour un certain pas de temps choisi'''
		#initialisation
		resultat = []
		#remplissage de tableau
		for i in range(len(self.rayons_evolution)):
			#concaténation
			resultat = resultat + [self.rayons_evolution[i][pas]/self.rayon_sphere]
		return(resultat)

	def get_vitesse_sur_Hubble_t(self, pas):
		'''Donne le tableau du rapport de la vitesse sur la valeur de la fonction d'Hubble
		en un certain temps choisi, pris sur tous les rayons des coquilles'''
		#initialisation
		resultat = []
		#remplissage de tableau
		for i in range(len(self.vitesses_evolution)):
			resultat = resultat + [self.vitesses_evolution[i][pas]/self.H[pas]]
		return(resultat)