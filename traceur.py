# -*- coding: utf-8 -*

import matplotlib.pyplot as plt
import random

class Tracer() :
	'''Va permettre de tracer et de visualiser les résultats du solver
	Présentation des résultats de façon adimensionnée'''

	def affichage_vitesse_Hubble(self, solveur, pas):
		'''Fonction qui affiche l'évolution du rapport de la vitesse sur 
		la constante d'Hubble adimensionné en fonction de l'évolution des
		rayons pour un pas de temps donné'''
		#tableau des abscisses
		abscisses = solveur.get_rayons_adimension_t(pas)
		#tableau des ordonnées
		ordonnées = solveur.get_vitesse_sur_Hubble_t(pas)
		#affichage du graphe
		plt.plot(abscisses, ordonnées)
		plt.xlabel('rayon')
		plt.ylabel('vitesse/Hubble')
		plt.title('Solveur Numérique')
		plt.savefig('vH_solv1.pdf')
		plt.show()
		return()

	def affichage_vitesse_Hubble_comparaison(self, solveur_anal, solveur_num, pas):
		'''Fonction qui affiche l'évolution du rapport de la vitesse sur 
		la constante d'Hubble adimensionné en fonction de l'évolution des
		rayons pour un pas de temps donné'''
		#tableau des abscisses
		abscisses = solveur_anal.get_rayons_adimension_t(pas)
		#tableau des ordonnées
		ordonnées1 = solveur_anal.get_vitesse_sur_Hubble_t(pas)
		ordonnées2 = solveur_num.get_vitesse_sur_Hubble_t(pas)
		#correction de l'échelle : je choisis un rayon au hasard dans
		#la liste, je regarde les rapports entre les solutions analytiques
		#et numériques pour ce rayon. Puis j'applique ce rapport à la
		#deuxième liste d'ordonnées (tous les rayons). De cette façon,
		#je remets mes solutions à la même échelle.
		rayon_choisi = random.randint(0,len(ordonnées1)-1)
		self.rapport = ordonnées1[rayon_choisi]/ordonnées2[rayon_choisi]
		nouvelle_ordonnée2 = []
		for i in range(len(ordonnées2)):
			nouvelle_ordonnée2 = nouvelle_ordonnée2 + [ordonnées2[i]*self.rapport]
		#affichage du graphe
		plt.plot(abscisses, ordonnées1, label='Solveur Analytique')
		plt.plot(abscisses, nouvelle_ordonnée2, label='Solveur Numérique')
		plt.legend(loc='upper left')
		plt.xlabel('rayon')
		plt.ylabel('vitesse/Hubble')
		plt.title('Comaparaison entre les solutions analytiques et numériques')
		plt.savefig('vH_solv.pdf')
		plt.show()
		return()