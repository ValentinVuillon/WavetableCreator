import numpy as np
import matplotlib.pyplot as plt
from scipy.io.wavfile import write

#---------------------------------------------------

#functions

#simple functions

#cosine of pulsation w
def fs1(x, y,w):   
    return np.cos(w*x)

def fs2(x, y, w):  
    return np.cos(w*x)*np.sin(w*y)

def fs3(x, y, w):  
    return np.cos(w*(y+1)*2*x)

#2D Bell curve with parameter s 
def fs4(x,y,s):    
    return np.exp(-(s*x)**2-(s*y)**2)

 #2D cardinal sine with pulsation w
def fs5(x,y,w):  
    return np.sin(w*x)*np.sin(w*y)/(w**2*x*y)




#complex functions (functions defined with other functions)
def fc1(x,y):
    return sign(fs1(x,y,10))

def fc2(x,y,w):
    return sign(fs3(x, y, w))
    

def fc3(x,y):
    return fs4(x,y,2)*fc2(x,y,4)

def fc4(x,y,w,s):
    return fs4(x,y,s)*fc2(x,y, w)

def fc5(x,y):
    return (fc4(x,y,100, 4)+fc4(x-0.5,y-0.5, 80, 2)+fc4(x+0.5,y+0.5, 20, 4)+fc4(x+0.5,y-0.5, 60, 4)+fc4(x-0.5,y+0.5, 30, 4))/5
   
def fc6(x,y):
    return fc5(x,y)*fs4(x,y,2)

def fc7(x,y):
    return demultiplication_translation_sum_2args(5, 0.7, x, y, fs4, 10)

def fc8(x,y):
    return demultiplication_translation_sum_4args(5, 0.7, x, y, fc4, 10, 4)

def fc9(x,y):
    return demultiplication_translation_sum_2args(5, 0.7, x, y, fc5)

def fc10(x,y):
    return demultiplication_translation_sum_3args(5, 0.7,x, y, fs5, 20)

def fc11(x,y):
    x_warp_3args(x, y, fs5,20)

def fc12(x,y):
    return rotation_3args(np.pi/4,x, y, fs1,10)

def fc13(x,y):
    return rotation_2args(np.pi/4+10**(-2),x, y, fc11)

def fc14(x,y):
    return rotation_2args(np.pi/4+10**(-2),x, y, fc9)

    

#operators (they act on 2D functions)

#simple operators

#takes a function and returns its sign
def sign(arr):  
    arr1=arr
    for i in range(waveform_number):
        for j in range(sample_number):
            if arr1[i][j]>0:
                arr1[i][j]=1
            if arr1[i][j]<0:
                arr1[i][j]=-1
    return arr1


#takes a function and demultiplies it. N: number of copies. d: defines the distance between copies
def demultiplication_translation_sum_2args(N,d, x,y, f):     
    arr1=np.zeros([waveform_number, sample_number])
    for a in np.linspace(-d,d,N):
        arr1+=f(x+a,y+a)
    if N%2==0:
        for a in np.linspace(-d,d,N):
            arr1+=f(x-a,y+a)
    else:
        for a in np.linspace(-d,d,N-1):
            arr1+=f(x-a,y+a)
      
    return arr1/4

def demultiplication_translation_sum_3args(N,d,x,y, f, extra_parameter1):
    arr1=np.zeros([waveform_number, sample_number])
    for a in np.linspace(-d,d,N):
        arr1+=f(x+a,y+a,extra_parameter1)
    if N%2==0:
        for a in np.linspace(-d,d,N):
            arr1+=f(x-a,y+a,extra_parameter1)
    else:
        for a in np.linspace(-d,d,N-1):
            arr1+=f(x-a,y+a,extra_parameter1)
      
    return arr1/4

def demultiplication_translation_sum_4args(N,d,x,y, f, extra_parameter1, extra_parameter2):
    arr1=np.zeros([waveform_number, sample_number])
    for a in np.linspace(-d,d,N):
        arr1+=f(x+a,y+a,extra_parameter1,extra_parameter2)
    if N%2==0:
        for a in np.linspace(-d,d,N):
            arr1+=f(x-a,y+a,extra_parameter1,extra_parameter2)
    else:
        for a in np.linspace(-d,d,N-1):
            arr1+=f(x-a,y+a,extra_parameter1,extra_parameter2)
      
    return arr1/4

#an operator which acts on x 
def x_warp_3args(x,y,f,extra_parameter1):  
    return f(-x**2,y,extra_parameter1)

#an operator which make the function rotate around -z with an angle theta [rad]
def rotation_2args(theta, x,y,f):  
    x_prime=np.cos(theta)*x-np.sin(theta)*y
    y_prime=np.sin(theta)*x+np.cos(theta)*y
    
    return f(x_prime, y_prime)

def rotation_3args(theta, x,y,f,extra_parameter1):  
    x_prime=np.cos(theta)*x-np.sin(theta)*y
    y_prime=np.sin(theta)*x+np.cos(theta)*y
    
    return f(x_prime, y_prime,extra_parameter1)

#complex operators (operators defined from other operators)

#to be created

#---------------------------------------------------

prototyping=0

if prototyping==True:
    sample_number=60
    waveform_number=60 
else:
    sample_number=2048  #number of samples per wavetables.For serum 2048. For ableton wavetable 1024.
    waveform_number=256 


x = np.linspace(-1, 1, sample_number)
y = np.linspace(-1, 1, waveform_number)

X, Y = np.meshgrid(x, y)

Z = fc4(X,Y,10,4)  #change the function name here. This function will become the wavetable


#I cut any values below -1 and above 1 it's what serum and ableton wavetable doe
for i in range(waveform_number):
    for j in range(sample_number):
        if Z[i][j]>1:
            Z[i][j]=1
        if Z[i][j]<-1:
            Z[i][j]=-1
            
            
#visualizing the wavetable
if prototyping==True:
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.plot_surface(X, Y, Z, rstride=1, cstride=1,cmap='jet', edgecolor = 'none')

print(Z.shape)


#.wav file is created
wav=np.ones(sample_number*waveform_number, dtype='float32') #dtype='float32' is import for ableton to read the files.
#but serum is ok with reading 64 bit files so it's not necessary to use the command for serum

for i in range(waveform_number):
    for j in range(sample_number):
        wav[j+i*sample_number]=Z[i][j]

rate = 44100 # any number can be used
if prototyping==False:
    write('wavetable.wav', rate, wav)
