#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm_notebook


# In[2]:


# Aceleración
grav = 9.81

# Masa del dipolo
m = 385+250

# es la permeabilidad máxima del MetGlass
mu = 1e6

# Permeabilidad del espacio libre
mu_0 = 4*np.pi*1e-7

# Radio del anillo de corriente
a = 2

# Resistividad
R = 9e-5

# Densidad del aire
ro= 1.225

# Coeficiente de arrastre
C= 1.05

# Área de contacto del juego
A= (3**2)*np.pi-(2**2)*np.pi


# In[3]:


k = (9*(mu*mu_0)**2*a**4)/(4*R)


# In[4]:


# Aceleración causada por la resistencia del aire

def Ar(y):
    return ((1/2)*(C*ro*A*y**2))/m


# In[5]:


def p(z):
    return z**2/np.sqrt((z**2+a**2)**5)
def f(y):
    return y

# Aceleración sin resistencia
def g(x, y):
    return -grav-(k/m)*p(x)*y

# Aceleración con resistencia
def gr(x, y):
    return -grav-(k/m)*p(x)*y+Ar(y)


# In[6]:


# Tiempo inicial
t_0 = 0
# Timpo final
t_f = 10

# Altura Inicial z_0
x_0 = 67
# Velocidad Inicial z'_0
y_0 = 0


# In[7]:


# Paso horizontal
h = 0.001

# Número de pasos
# Calcula el número de pasos entre t_0 y t_f cuando tienes un paso horizontal de h
N = int ((t_f-t_0)/h)

# Arreglo con valores del tiempo:
# Usa np.linspace() para definir un arreglo que vaya de t_0 a t_f en N pasos
t = np.linspace(t_0, t_f, N)


# In[8]:


all_xs = [] 
all_ys = []

# Usar los valores iniciales para el primer paso
x_n = x_0
y_n = y_0

# Iterar para cada paso de la lista de tiempo t
# Usamos tqdm para generar una barra de progreso
for t_n in tqdm_notebook(t):
    
    kx1= f(y_n)
    ky1= g(x_n, y_n)
    
    kx2= f(y_n+(1/2)*h*ky1)
    ky2= g(x_n+(1/2)*h*kx1, y_n+(1/2)*h*ky1)
    
    kx3= f(y_n+(1/2)*h*ky2)
    ky3= g(x_n+(1/2)*h*kx2, y_n+(1/2)*h*ky2)
    
    kx4= f(y_n+h*ky3)
    ky4= g(x_n+h*kx3, y_n+h*ky3)
    
    x_nplus1= x_n+h/6*(kx1+2*kx2+2*kx3+kx4)
    y_nplus1= y_n+h/6*(ky1+2*ky2+2*ky3+ky4)
    
    all_xs.append(x_nplus1)
    all_ys.append(y_nplus1)
    
    x_n = x_nplus1
    y_n = y_nplus1


# In[9]:


all_xrs = [] 
all_yrs = []

x_n = x_0
y_n = y_0

for t_n in tqdm_notebook(t):
    
    kx1= f(y_n)
    ky1= gr(x_n, y_n)
    
    kx2= f(y_n+(1/2)*h*ky1)
    ky2= gr(x_n+(1/2)*h*kx1, y_n+(1/2)*h*ky1)
    
    kx3= f(y_n+(1/2)*h*ky2)
    ky3= gr(x_n+(1/2)*h*kx2, y_n+(1/2)*h*ky2)
    
    kx4= f(y_n+h*ky3)
    ky4= gr(x_n+h*kx3, y_n+h*ky3)
    
    x_nplus1= x_n+h/6*(kx1+2*kx2+2*kx3+kx4)
    y_nplus1= y_n+h/6*(ky1+2*ky2+2*ky3+ky4)
    
    all_xrs.append(x_nplus1)
    all_yrs.append(y_nplus1)
    
    x_n = x_nplus1
    y_n = y_nplus1


# In[10]:


all_xs = np.array(all_xs)
all_ys = np.array(all_ys)

fz = -k*p(all_xs)*all_ys
az = +fz/m - grav


# In[11]:


all_xrs = np.array(all_xrs)
all_yrs = np.array(all_yrs)

frz = -k*p(all_xrs)*all_yrs + Ar(all_yrs)*m
arz = +frz/m - grav


# In[12]:


plt.figure(figsize=(18,12))

plt.subplot(4,2,1)
plt.plot(t,all_xs)
plt.title("Altura de Dipolo con Frenado Magnético")
plt.xlabel("Tiempo [s]")
plt.ylabel("Altura [m]")

plt.subplot(4,2,2)
plt.plot(t,all_ys)
plt.title("Velocidad de Dipolo con Frenado Magnético")
plt.xlabel("Tiempo [s]")
plt.ylabel("Velocidad [m/s]")

plt.subplot(4,2,3)
plt.plot(t,az)
plt.title("Aceleración del dipolo")
plt.xlabel("Tiempo [s]")
plt.ylabel("Aceleracion [m/s^2]")

plt.subplot(4,2,4)
plt.plot(t,fz)
plt.title("Fuerza del campo magnético inducido sobre el dipolo")
plt.xlabel("Tiempo [s]")
plt.ylabel("Fuerza [N]")

plt.subplot(4,2,5)
plt.plot(t,all_xrs)
plt.title("Altura de Dipolo con Frenado Magnético y resistencia del aire")
plt.xlabel("Tiempo [s]")
plt.ylabel("Altura [m]")

plt.subplot(4,2,6)
plt.plot(t,all_yrs)
plt.title("Velocidad de Dipolo con Frenado Magnético y resistencia del aire")
plt.xlabel("Tiempo [s]")
plt.ylabel("Velocidad [m/s]")

plt.subplot(4,2,7)
plt.plot(t,arz)
plt.title("Aceleración del dipolo con resistencia del aire")
plt.xlabel("Tiempo [s]")
plt.ylabel("Aceleracion [m/s^2]")

plt.subplot(4,2,8)
plt.plot(t,frz)
plt.title("Fuerza de la resistencia del aire y campo magnético inducido sobre el dipolo")
plt.xlabel("Tiempo [s]")
plt.ylabel("Fuerza [N]")
plt.tight_layout()


# In[ ]:




