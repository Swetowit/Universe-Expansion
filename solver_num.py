# -*- coding: utf-8 -*

from math import *

class Solver_Num():
	'''Solveur numérique pour le calcul de l'évolution des souspressions
	en faisant l'hypothèse d'un modèle sphérique.'''

	def __repr__(self):
		return("Ceci est un solveur numérique.")

	def __init__(self, rayon_sphere, nb_divisions_sphere, densite_background_initiale, souspression, nb_pas_temps):
		'''Calcul de toutes les grandeurs à l'instanciation de l'objet'''
		# Variables stockées
		self.rayon_sphere = rayon_sphere/10
		self.nb_divisions_sphere = nb_divisions_sphere
		self.densite_background_initiale = densite_background_initiale/(10**(-27))
		self.souspression = souspression
		self.nb_pas_temps = nb_pas_temps
		self.G = 1
		#Données calculées
		self.pas_t = 10**16/self.nb_pas_temps
		(self.rayons_initiaux, self.masses) = self.calcul_rayons_masses_initiaux()
		self.a = self.calcul_facteur_echelle()
		self.densite_background_evolution = self.calcul_densite()
		self.H = self.calcul_Hubble()
		self.energie =self.calcul_energie()
		self.rayons_evolution = self.calcul_evolution_rayons()
		self.vitesses_evolution = self.calcul_evolution_vitesse()

	def calcul_rayons_masses_initiaux(self):
		'''Calcule les rayons et les masses initiales des shells
		Le modèle est construit de telle sorte à ce qu'il y ait
		autant de shell à l'extérieur qu'à l'intérieur de la sphère
		afin qu'il y ait une transition continue entre les deux'''
		rayons_initiaux = []
		masses = []
		for i in range(self.nb_divisions_sphere):
			rayons_initiaux = rayons_initiaux +[(self.rayon_sphere/self.nb_divisions_sphere)*(i+1)]
			masses = masses + [(4/3)*pi*(rayons_initiaux[i]**3)*self.densite_background_initiale*(1+self.souspression)]
		for i in range (self.nb_divisions_sphere):
			rayons_initiaux = rayons_initiaux + [(self.rayon_sphere/self.nb_divisions_sphere)*(self.nb_divisions_sphere+i+1)]
			masses = masses + [(4/3)*pi*(rayons_initiaux[self.nb_divisions_sphere+i]**3)*self.densite_background_initiale*(1+((self.rayon_sphere/rayons_initiaux[i+self.nb_divisions_sphere])**3)*self.souspression)]
		return(rayons_initiaux, masses)

	def calcul_facteur_echelle(self):
		'''Donne un tableau du facteur d'échelle en fonction du temps'''
		a =[]
		for i in range(self.nb_pas_temps):
			a = a +[(self.pas_t*(i+1))**(2/3)]
		return(a)

	def calcul_densite(self):
		'''Donne le tableau de l'évolution de la densité du background
		en fonction du temps'''
		rho = []
		for i in range(self.nb_pas_temps):
			rho = rho + [self.densite_background_initiale/(self.a[i]**3)]
		return(rho)

	def calcul_Hubble(self):
		'''Donne la fonction d'Hubble en fonction du temps
		Ici H est adimensioné par H0, on calcule donc en fait H/H0'''
		H =[]
		for i in range(self.nb_pas_temps):
			H= H +[sqrt(8*pi*self.G*self.densite_background_evolution[i])/(3*self.densite_background_initiale)]
		return(H)

	def calcul_energie(self):
		'''Donne l'énergie contenue dans chaque shell'''
		energie = []
		for i in range(len(self.rayons_initiaux)):
			energie = energie + [0.5*(self.H[0]*self.rayons_initiaux[i])**2]
		return(energie)

	def equation_newton(self, masse, rayon, energie):
		'''Donnée de l'équation de newton qui sera utilisée pour
		trouver l'évolution des rayons.
		Elle fonctionne pour des nombres et non pas pour des tableaux'''
		return(sqrt(2*self.G*masse/rayon+2*energie))

	def calcul_evolution_rayons(self):
		'''Calcule l'évolution des rayons des shells
		grâce à la méthode RK4 en fonction du temps'''
		rayons_evolution = [0]*len(self.rayons_initiaux)
		# Calcul pour chaque shell (donc chaque rayon initiale)
		for i in range(len(rayons_evolution)):
			rayons_evolution[i] = [self.rayons_initiaux[i]]
			# Une fois la shell sélectionnée, on calcule l'évolution
			# de son rayon dans le temps
			for t in range(self.nb_pas_temps):
				k1 = self.equation_newton(self.masses[i], rayons_evolution[i][t], self.energie[i])
				k2 = self.equation_newton(self.masses[i], rayons_evolution[i][t] + (self.pas_t*k1)/2, self.energie[i])
				k3 = self.equation_newton(self.masses[i], rayons_evolution[i][t] + (self.pas_t*k2)/2, self.energie[i])
				k4 = self.equation_newton(self.masses[i], rayons_evolution[i][t] + self.pas_t*k3, self.energie[i])
				rayons_evolution[i] = rayons_evolution[i] + [rayons_evolution[i][t] + (self.pas_t/6)*(k1+2*k2+2*k3+k4)]
		return(rayons_evolution)

	def calcul_evolution_vitesse(self):
		'''Calcule l'évolution de la vitesse de l'expansion du rayon
		des shells en fonction du temps'''
		vitesse = [0]*len(self.rayons_evolution)
		for i in range(len(self.rayons_evolution)):
			vitesse[i] = []
			for t in range(len(self.rayons_evolution[0])):
				vitesse[i] = vitesse[i] + [self.equation_newton(self.masses[i], self.rayons_evolution[i][t], self.energie[i])]
		return(vitesse)

	def get_rayons_adimension_t(self, pas):
		'''Permet de récupérer le rayon des shells divisé par le 
		rayon de la sphère pour un certain instant t ou pas précisé'''
		res = []
		for i in range(len(self.rayons_initiaux)):
			res = res + [self.rayons_evolution[i][pas]/self.rayon_sphere]
		return(res)

	def get_vitesse_sur_Hubble_t(self, pas):
		'''Permet de récupérer la vitesse d'expansion des shells divisé
		par la valeur de la fonction d'Hublle a un certain pas/instant
		est adimensionné. H est déjà adimensionné.'''
		res = []
		for i in range(len(self.vitesses_evolution)):
			res = res + [self.vitesses_evolution[i][pas]/self.H[pas]]
		return(res)