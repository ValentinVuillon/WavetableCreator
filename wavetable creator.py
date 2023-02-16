import numpy as np
import matplotlib.pyplot as plt
from scipy.io.wavfile import write

#a faire:
#-ajouter labels x, y, z sur plot matplotlib
#-creer un operateur qui prend une fonction t fait une serie de ces fonctions, avec rotation, en utilisant l'operateur rotation.
#-creer un operateur qui fait un filtrage->convolution avec une gaussienne.
#creer un autre code pour ça: créeer un code qui prend une wavetable de serum, la transforme en fonction 2D puis fait un filtrage avec un kernel. Ensuite en changeant le kernel je changerai la nature du filtre.

#ce code permet de creer des wavetables pour serum, ableton wavetable. 

#les axes sont
#      1 y
#       /
#      /
#     /
# -1 /
#  -1-----------1x
#
#et z est vers le haut

#---------------------------------------------------

#fonctions

#fonctions simples

#cosinus de pulsation w
def fs1(x, y,w):   
    return np.cos(w*x)

def fs2(x, y, w):  
    return np.cos(w*x)*np.sin(w*y)

def fs3(x, y, w):  
    return np.cos(w*(y+1)*2*x)

#gaussienne 2D d'ecart type s
def fs4(x,y,s):    
    return np.exp(-(s*x)**2-(s*y)**2)

 #sinus cardinal 2D de pulsation w
def fs5(x,y,w):  
    return np.sin(w*x)*np.sin(w*y)/(w**2*x*y)




#fonctions complexes (fonctions définies à partir d'autres fonctions)
def fc1(x,y):
    return signe(fs1(x,y,10))

def fc2(x,y,w):
    return signe(fs3(x, y, w))
    

def fc3(x,y):
    return fs4(x,y,2)*fc2(x,y,4)

def fc4(x,y,w,s):
    return fs4(x,y,s)*fc2(x,y, w)

def fc5(x,y):
    return (fc4(x,y,100, 4)+fc4(x-0.5,y-0.5, 80, 2)+fc4(x+0.5,y+0.5, 20, 4)+fc4(x+0.5,y-0.5, 60, 4)+fc4(x-0.5,y+0.5, 30, 4))/5
   
def fc6(x,y):
    return fc5(x,y)*fs4(x,y,2)

def fc7(x,y):
    return demultiplication_translation_somme_2args(5, 0.7, x, y, fs4, 10)

def fc8(x,y):
    return demultiplication_translation_somme_4args(5, 0.7, x, y, fc4, 10, 4)

def fc9(x,y):
    return demultiplication_translation_somme_2args(5, 0.7, x, y, fc5)

def fc10(x,y):
    return demultiplication_translation_somme_3args(5, 0.7,x, y, fs5, 20)

def fc11(x,y):
    x_warp_3args(x, y, fs5,20)

def fc12(x,y):
    return rotation_3args(np.pi/4,x, y, fs1,10)

def fc13(x,y):
    return rotation_2args(np.pi/4+10**(-2),x, y, fc11)

def fc14(x,y):
    return rotation_2args(np.pi/4+10**(-2),x, y, fc9)

    

#operateurs (agissent sur fonctions 2D)

#operateurs simples

#prend une fonction 2D et retourne son signe
def signe(arr):  
    arr1=arr
    for i in range(waveform_number):
        for j in range(sample_number_x):
            if arr1[i][j]>0:
                arr1[i][j]=1
            if arr1[i][j]<0:
                arr1[i][j]=-1
    return arr1


#prend une fonction 2D et la démultiplie. N: nombre de copies. d: définit distance entre les copies.
def demultiplication_translation_somme_2args(N,d, x,y, f):     
    arr1=np.zeros([waveform_number, sample_number_x])
    for a in np.linspace(-d,d,N):
        arr1+=f(x+a,y+a)
    if N%2==0:
        for a in np.linspace(-d,d,N):
            arr1+=f(x-a,y+a)
    else:
        for a in np.linspace(-d,d,N-1):
            arr1+=f(x-a,y+a)
      
    return arr1/4

def demultiplication_translation_somme_3args(N,d,x,y, f, extra_parameter1):
    arr1=np.zeros([waveform_number, sample_number_x])
    for a in np.linspace(-d,d,N):
        arr1+=f(x+a,y+a,extra_parameter1)
    if N%2==0:
        for a in np.linspace(-d,d,N):
            arr1+=f(x-a,y+a,extra_parameter1)
    else:
        for a in np.linspace(-d,d,N-1):
            arr1+=f(x-a,y+a,extra_parameter1)
      
    return arr1/4

def demultiplication_translation_somme_4args(N,d,x,y, f, extra_parameter1, extra_parameter2):
    arr1=np.zeros([waveform_number, sample_number_x])
    for a in np.linspace(-d,d,N):
        arr1+=f(x+a,y+a,extra_parameter1,extra_parameter2)
    if N%2==0:
        for a in np.linspace(-d,d,N):
            arr1+=f(x-a,y+a,extra_parameter1,extra_parameter2)
    else:
        for a in np.linspace(-d,d,N-1):
            arr1+=f(x-a,y+a,extra_parameter1,extra_parameter2)
      
    return arr1/4

#une transformation qui agit sur l'ensemble de valeurs x, y c'est comme les warps de serum (on pourrait faire les bend warp comme ça)
def x_warp_3args(x,y,f,extra_parameter1):  
    return f(-x**2,y,extra_parameter1)

#operateur qui fait une rotation suivant axe -z d'un angle theta en radian de la fonction donnée en argument
def rotation_2args(theta, x,y,f):  
    x_prime=np.cos(theta)*x-np.sin(theta)*y
    y_prime=np.sin(theta)*x+np.cos(theta)*y
    
    return f(x_prime, y_prime)

def rotation_3args(theta, x,y,f,extra_parameter1):  
    x_prime=np.cos(theta)*x-np.sin(theta)*y
    y_prime=np.sin(theta)*x+np.cos(theta)*y
    
    return f(x_prime, y_prime,extra_parameter1)

#operateurs complexes (operateurs définis à partir d'autres operateurs)


#---------------------------------------------------

prototypage=0

if prototypage==True:
    sample_number_x=60
    waveform_number=60 
else:
    sample_number_x=2048  #nb de samples dans une wavetable. Pour serum mettre 2048. Pour wavetable mettre 1024.
    waveform_number=256 


x = np.linspace(-1, 1, sample_number_x)
y = np.linspace(-1, 1, waveform_number)

X, Y = np.meshgrid(x, y)

Z = fc4(X,Y,10,4)  #changer ici la fonction, c'est elle qui formera la wavetable


#je coupe tout valeur au delà de -1 et 1, c'est ce que fait serum, ableton wavetable
for i in range(waveform_number):
    for j in range(sample_number_x):
        if Z[i][j]>1:
            Z[i][j]=1
        if Z[i][j]<-1:
            Z[i][j]=-1
            
            
#visualisation de la wavetable que serum créera
if prototypage==True:
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.plot_surface(X, Y, Z, rstride=1, cstride=1,cmap='jet', edgecolor = 'none')

print(Z.shape)


#je cree le fichier .wav
wav=np.ones(sample_number_x*waveform_number, dtype='float32') #dtype='float32' permet d'avoir des fichiers lisibles par ableton, ableton wavetable mais ça ne sert pas le cas de serum,
#ne pas mettre la commande crée des fichiers 64 bit que serum arrive à lire

for i in range(waveform_number):
    for j in range(sample_number_x):
        wav[j+i*sample_number_x]=Z[i][j]

rate = 44100 #je peux mettre n'importe quoi, c'est sans importance
if prototypage==False:
    write('wavetable.wav', rate, wav)

#consignes d'importation:
#serum: importer avec decoupe tous les 2048 samples
#ableton wavetable: importer puis mettre option "raw"
