# -*- coding: utf-8 -*

from math import *
from scipy.optimize import fsolve

class Solver_Anal():
	

	def __init__(self, rayon_sphere, nb_divisions_sphere, densite_background_initiale, souspression, nb_pas_temps):
		'''Constructeur'''
		#stockage des données récupérées à l'instanciation de l'objet
		self.rayon_sphere = rayon_sphere/10
		self.nb_divisions_sphere = nb_divisions_sphere
		self.densite_background_initiale = densite_background_initiale/(10**(-27))
		self.souspression = souspression
		self.nb_pas_temps = nb_pas_temps
		self.G = 1
		#nouvelles valeurs calculées par l'objet
		self.pas_t = 10**16/self.nb_pas_temps
		(self.rayons_initiaux, self.masses) = self.calcul_rayons_masses_initiaux()
		self.a = self.calcul_facteur_echelle()
		self.densite_background_evolution = self.calcul_densite()
		self.H = self.calcul_Hubble()
		self.energie =self.calcul_energie()

		self.tableau_temps = self.time()
		self.theta = self.calcul_theta()
		self.masse_sphere = (4/3)*pi*((self.rayon_sphere)**3)*self.densite_background_initiale*(1+self.souspression)
		self.A = self.calcul_A()
		self.B = self.calcul_B()
		self.rayons_evolution = self.calcul_rayons_evolution()
		self.vitesses_evolution = self.calcul_vitesse()

	'''Des fonctions qu'on pourrait faire hériter pour les solveurs car ça me soule de faire
	copier-coller '''

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

	def time(self):
		'''Fonction qui donne un tableau de temps différents,
		en relation avec les temps donnés dans le constructeur'''
		temps = []
		for i in range(1,self.nb_pas_temps+1):
			temps =  temps + [self.pas_t*i]
		return(temps)

	def calcul_A(self):
		'''Fonction qui calcule le A dans la solution analytique du problème
		A est différent pour chaque rayobn'''
		constante_H = self.H[0]
		A = []
		for i in range(len(self.rayons_initiaux)):
			A = A + [(self.G*self.masses[i])/((constante_H*self.rayons_initiaux[i])**2)]
		return(A)

	def calcul_B(self):
		'''Fonction qui calcule le B dans la solution analytique du problème
		B calculé avec le A du rayon/masse de la sphère de sous-pression'''
		n = self.nb_divisions_sphere
		return(sqrt((self.A[n]**3)/(self.G*self.masses[n])))

	def calcul_theta(self):
		'''fonction à calculer pour connaître theta : t - B(sinh(theta) - theta) = 0
		Résolution avec fsolve de scipy.optimize, je ne sais toujours pas comment ça
		fonctionne.'''
		theta = []
		for i in range(len(self.tableau_temps)):
			#creation de la fonction anonyme t - B(sinh(theta) - theta) pour chaque t
			#Problème avec le B qui est tout petit et qui fait tout sauter !
			fonction_resoudre = lambda x : self.tableau_temps[i]/(10**16) - (sinh(x)-x)
			#on trouve le 0 de la fonction précédente avec un paramètre 1 (ça marche lol)
			#puis conversion en un float pour le stocker dans un tableau sinon ça renvoie
			#array([trucmuche]) et ça c'est inutilisable
			theta = theta + [float(fsolve(fonction_resoudre,1))]
		return(theta)

	def calcul_rayons_evolution(self):
		'''Calcul du rayon analytique en fonction de theta
		R = A * (cosh(theta) - 1)'''
		rayons_evolution = [0]*len(self.rayons_initiaux)
		for i in range(len(rayons_evolution)):
			rayons_evolution[i] = []
			for j in range(self.nb_pas_temps):
				rayons_evolution[i] = rayons_evolution[i] + [self.A[i]*(cosh(self.theta[j])-1)]
		return(rayons_evolution)

	def equation_newton(self, masse, rayon, energie):
		'''Fonction appelée par le solveur
		y'=f(y) et correspond à la dérivée du rayon par rapport au temps'''
		res = sqrt(2*self.G*masse/rayon + 2*energie) #J'appelle G comme self.G
		return(res)

	def calcul_vitesse(self):
		'''Nous renvoie un tableau de l'évolution des rayons des coquilles,
		l'évolution des vitesses pour chaque rayon de coquille
		et s'il existe un shell-collapse, si oui à quel pas survient-il'''
		#Tableau de l'évolution des vitesses des coquilles
		vitesse_evolution = [0]*len(self.rayons_evolution)
		for i in range(len(vitesse_evolution)):
			vitesse_evolution[i] = [0]*len(self.rayons_evolution[0]) #création d'une colonne
			for t in range(len(self.rayons_evolution[0])):
				vitesse_evolution[i][t] = self.equation_newton(self.masses[i], self.rayons_evolution[i][t], self.energie[i]) #remplissage de chaque valeur
		return(vitesse_evolution)

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