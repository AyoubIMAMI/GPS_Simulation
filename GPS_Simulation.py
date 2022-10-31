# -*- coding: utf-8 -*-
"""
Created on Thu Oct 15 15:42:21 2020

@author: Ayoub IMAMI
@author: Adame BORKI
"""
from operator import itemgetter
from math import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # Fonction pour la 3D
from mpl_toolkits import mplot3d
import numpy as np


RAAN0=32.847
a1=["A1",RAAN0+240,268.126]
a2=["A2",RAAN0+240,161.786]
a3=["A3",RAAN0+240,11.676]
a4=["A4",RAAN0+240,41.806]
b1=["B1",RAAN0+300,80.956]
b2=["B2",RAAN0+300,173.336]
b3=["B3",RAAN0+300,309.976]
b4=["B4",RAAN0+300,204.376]
c1=["C1",RAAN0,111.876]
c2=["C2",RAAN0,11.796]
c3=["C3",RAAN0,339.666]
c4=["C4",RAAN0,241.556]
d1=["D1",RAAN0+60,135.226]
d2=["D2",RAAN0+60,265.446]
d3=["D3",RAAN0+60,35.156]
d4=["D4",RAAN0+60,167.356]
e1=["E1",RAAN0+120,197.046]
e2=["E2",RAAN0+120,302.596]
e3=["E3",RAAN0+120,66.066]
e4=["E4",RAAN0+120,333.686]
f1=["F1",RAAN0+180,238.886]
f2=["F2",RAAN0+180,345.226]
f3=["F3",RAAN0+180,105.206]
f4=["F4",RAAN0+180,135.346]

satellites=[a1,a2,a3,a4,b1,b2,b3,b4,c1,c2,c3,c4,d1,d2,d3,d4,e1,e2,e3,e4,f1,f2,f3,f4]
Rterre=6400
Rsat=26600
c=300000


def angleradian(x):
    return ((x*pi)/180)
    
# Détermination des coordonnées cartésiennes du point M récepteur
def pointM(longitude,latitude,altitude):
    a=angleradian(longitude)
    b=angleradian(latitude)
    h=Rterre+altitude
    M=[h*np.cos(a)*np.cos(b),h*np.sin(a)*np.cos(b),h*np.sin(b)]
    return M

# Détermination de la position des satellites en coordonnées cartésiennes

def positionsatellite(t): 
    position=[]
    for i in range(len(satellites)) :
        omega = (2*pi)/(12*60*60)
        phi = angleradian(satellites[i][2])+(omega*t)%(2*pi)
        raan = angleradian(satellites[i][1])
        r = angleradian(55)
        # x,y,z déterminées à partir des figures planes de rotation
        x = Rsat*np.cos(raan)*np.cos(phi) - Rsat*np.sin(phi)*np.cos(r)*np.sin(raan)
        y = Rsat*np.sin(raan)*np.cos(phi) - Rsat*np.sin(phi)*np.cos(r)*np.cos(raan)
        z = Rsat*np.sin(phi)*np.sin(r)
        position.append([satellites[i][0],x,y,z])
    return position
    
    
# Figure avec la Terre,les satellites et les orbites 
# Les satellites de la même couleur sont sur la même orbite et qui est aussi de la même couleur
    
def FigureGlobal(t):  
    positioncarté=positionsatellite(t)
    couleur=['b','g','r','y','k','m']
    fig = plt.figure()
    ax = fig.gca(projection='3d')  # Affichage en 3D  
    phi, theta = np.mgrid[0.0:pi:100j, 0.0:2.0*pi:100j] 
    # Tracé de la sphère qui représente la terre
    a = Rterre*np.sin(phi)*np.cos(theta)
    b = Rterre*np.sin(phi)*np.sin(theta)
    c = Rterre*np.cos(phi)
    ax.plot_surface(a, b, c,  rstride=10, cstride=1, color='c', alpha=0.7, linewidth=0)
    for i in range (0,24,4):
        # Récupération des coordonnées des satellites
        for j in range(4) :
                x=positioncarté[i+j][1]
                y=positioncarté[i+j][2]
                z=positioncarté[i+j][3]
                # Tracé des satellites 
                if j==0:
                    ax.scatter(x, y, z,s=200,c=couleur[i//4],marker='o',label='Orbite '+str(satellites[i][0][0]))
                else :
                    ax.scatter(x, y, z,s=200,c=couleur[i//4],marker='o')
        # Tracé des orbites
        # On coupe l'orbite en 720 points chacun séparé de 60 secondes
        for omega in range(0,720) :
            phi = angleradian(satellites[i][2])+(omega*60*t)%(2*pi)
            raan = angleradian(satellites[i][1])
            r = angleradian(55)
            x = Rsat*np.cos(raan)*np.cos(phi) - Rsat*np.sin(phi)*np.cos(r)*np.sin(raan)
            y = Rsat*np.sin(raan)*np.cos(phi) - Rsat*np.sin(phi)*np.cos(r)*np.cos(raan)
            z = Rsat*np.sin(phi)*np.sin(r)
            # Les orbites sont en fait des suites de 720 points qui les composent
            ax.scatter(x, y, z,s=5,edgecolor = 'none',c=couleur[i//4],marker='.')
    plt.title("Représentation des orbites")
    plt.legend(loc=2)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    plt.tight_layout()
    plt.show()
    
# Détermination des satellites visibles par le récepteur
    
# Cette fonction affiche la liste des noms des satellites visibles 
# Return une liste du type ['nom du satellite',distance au point M]    
    
def satellitesvisibles(longitude,latitude,altitude,t):
    distancestriées=[]
    ldistance=[]
    lvisible=[]
    L=positionsatellite(t)
    M=pointM(longitude,latitude,altitude)
    # On calcule la distance récepteur-satellite avec un produit scalaire
    for s in range (len(L)) :
        # On forme le vecteur MS entre le récepteur et le satellite
        vecteurMS=[L[s][1]-M[0],L[s][2]-M[1],L[s][3]-M[2]]
        # Calcul du produit scalaire
        ps = M[0]*vecteurMS[0]+M[1]*vecteurMS[1]+M[2]*vecteurMS[2]
        # Si le produit scalaire est positif, alors le satellite est visible par le récepteur
        if ps>0 :
            lvisible.append(L[s][0])   
            # On calcule la norme du vecteur MS pour obtenir la distance entre le récepteur et le satellite
            distanceMS=sqrt((L[s][1]-M[0])**2 + (L[s][2]-M[1])**2 + (L[s][3]-M[2])**2)
            ldistance.append([L[s][0],distanceMS])
    # On trie la liste des distances par rapport à la deuxième variable (les diatances)
    distancestriées=sorted(ldistance, key=itemgetter(1))
    print("Les satellites visibles sont "+str(lvisible)+".")
    return (distancestriées)
  
# Figure avec la Terre, les satellites visibles 
# Affiachant en couleur les quatre plus proches, les autres en noir et le récepteur en blanc  
  
def FigureVisibles(longitude,latitude,altitude,t):
    satellitesvisu=satellitesvisibles(longitude,latitude,altitude,t)
    positioncarté=positionsatellite(t)
    M=pointM(longitude,latitude,altitude)
    couleur = ['r', 'g', 'c','y','k','k','k','k','k']
    ax = Axes3D(plt.figure())   # Affichage en 3D
    phi, theta = np.mgrid[0.0:pi:100j, 0.0:2.0*pi:100j]
    # Tracé de la sphère qui représente la Terre
    a = Rterre*np.sin(phi)*np.cos(theta)
    b = Rterre*np.sin(phi)*np.sin(theta)
    c = Rterre*np.cos(phi)
    ax.plot_surface(a, b, c,  rstride=10, cstride=1, color='c', alpha=0.7, linewidth=0)
    # Tracé du point récepteur
    ax.scatter(M[0],M[1],M[2],s=200,c='w',marker='o')
    for i in range (len(satellitesvisu)):
        for k in positioncarté :
            if satellitesvisu[i][0]==k[0]:
               x=k[1]
               y=k[2]
               z=k[3] 
               # Tracé des satellites visibles
               ax.scatter(x, y, z,c=couleur[i],s=200,marker='o',label=(satellitesvisu[i][0],satellitesvisu[i][1]))
               # Tracé des droites entre le récepteur et les différents satellites
               ax.plot3D([x,M[0]],[y,M[1]],[z,M[2]],c='grey')
   
    plt.title("Satellites visibles et les 4 plus proches")
    plt.legend(loc=2)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    plt.show()

# Détermination des durées de propagation de chaque signal des quatre satellites les plus proches

# Affiche une liste des noms des quatres satellites les plus proches
# Return une liste du type ['nom du satellite',durée de propagation]
      
def duréepropa(longitude, latitude,altitude,t):
    listedurée=[]
    listeplusproches=[]
    distancestriéesbis=satellitesvisibles(longitude, latitude,altitude,t)
    # On récupère les distances des quatre plus proches dans la liste des distances triées
    for i in range (4):
       durée=distancestriéesbis[i][1]/(c*10**(-3))
       listedurée.append([distancestriéesbis[i][0],durée]) 
       listeplusproches.append(distancestriéesbis[i][0])
    print("Les quatres satellites les plus proches sont "+str(listeplusproches)+".")
    return listedurée   
 
  
# Permet de convertir un nombre en bianire et de sélectionner en argument le nombre de bits souhaités dans la trame
  
def dec2bin(d,nb=0):
    if d==0:
        b="0"
    else:
        b=""
        while d!=0:
            b="01"[d&1]+b
            d=d>>1
    return b.zfill(nb)
    
    
# Détermination de l'instant de départ commun  
# Affiche l'instant d'émission en secondes
# Return l'instant d'émisiion en binaire sur une trame de 20 bits
    
    
def instant_départ(longitude,latitude,altitude,t):
    instant_récep = t
    satvisi = satellitesvisibles(longitude,latitude,altitude,t)
    # Récupértion de la distance du 4e satellite le plus proche
    instant_trajet = satvisi[3][1]/(c*10**(-3))
    instant_émis = int(instant_récep - instant_trajet)
    print(instant_émis)
    instant_émis_binaire = dec2bin(instant_émis,20)
    return instant_émis_binaire














