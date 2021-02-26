import numpy as np
import matplotlib.pyplot as plt
from numpy.core.numeric import full
import cmath
b_old =[12.4727649,0.005325123]
a_old=[0.190193385,0.197129232]


def get_data(element):
    data = np.genfromtxt("peak_data/{}.txt".format(element))
    
    channel = data[:,0]
    channel_error = data[:,1]
    lit_energy = data[:,2]
  
    return channel,channel_error,lit_energy
def cal_data(channel,channel_error,a,b,c,lit_energy):
    energy =[]
    for ch in channel:
        d = (b[0]**2) - (4*a[0]*(c[0]-ch))
        #coeff = [b[0],a[0],c[0]-ch]
        
        #energy = (channel - a[0])/b[0]
        energy.append((-b[0]+cmath.sqrt(d))/(2*a[0]))
    energy_error = np.sqrt(np.square(np.sqrt(np.square(channel_error)+np.square(a[1]))/(channel-a[0]))+np.square(b[1]/b[0]))*energy
   
    residual = lit_energy-energy
    
    return residual,energy_error,channel,energy
elements = ["Ag","Ba","Cu","Mo","Rb","Tb"]
atomic_number = np.asarray([47,56,29,42,37,65])
element_data =[]

x_allele=np.array([])
y_allele=np.array([])
for e in elements:
    
    data = get_data(e)
    element_data.append(data)
    
    
    
    x_allele =np.concatenate((x_allele,data[2]),axis=None)
    y_allele =np.concatenate((y_allele,data[0]),axis=None)
    

    

coeff,cov = np.polyfit(x_allele,y_allele,2,cov=True,full=False)
error = np.sqrt(np.diag(cov))
b=[coeff[0],error[0]]
a=[coeff[1],error[1]]
c=[coeff[2],error[2]]

# b = [np.mean(b_s),np.std(b_s)/np.sqrt(len(b_s))]   
# a = [np.mean(a_s),np.std(a_s)/np.sqrt(len(a_s))] 
#print("calibration constants:"+str(b)+"---"+str(a))
plt.title("Energy v/s Channel")
plt.suptitle("linear fit : "+"y=%.6fx+(%.6f)"%(coeff[0],coeff[1]))

plt.ylabel("Energy(keV)")
plt.xlabel("Channel No.")
plt.legend(["Ag","Ba","Cu","Mo","Rb","Tb"])
 # fig.suptitle('Horizontally stacked subplots')
for count,e in enumerate(element_data):
    plt.plot(e[2],e[0],"o",label = elements[count])

# calc the trendline
p = np.poly1d(coeff)
print(p)
linear_fit = p(x_allele)

plt.plot(x_allele,linear_fit,"r-",label="Linear Fit")
# the line equation:
print("Calibration Const.Eq")
print("y=%.6fx+(%.6f)"%(coeff[0],coeff[1]))
plt.legend(loc="lower right")
# x =int(np.around(p(np.sqrt((get_data("X")[0][0]-a[0])/b[0]))))
# print("Element X has atomic number :"+ str(x))
plt.show()
plt.close()

res_data = []
for data in element_data:
    res_data.append(cal_data(data[0],data[1],a,b,c,data[2]))
    """Plotting Residual
    """

plt.title("Residual Plot : Energy Difference against Channels")
plt.axhline(y=0)
plt.ylabel("Residual")
plt.xlabel("Channel")
plt.legend(["Ag","Ba","Cu","Mo","Rb","Tb"])
 # fig.suptitle('Horizontally stacked subplots')

for count,e in enumerate(res_data):
    plt.errorbar(e[2],e[0],e[1],0,"o",label = elements[count])
plt.legend(loc="upper right")
plt.show()
plt.close()






"""Moseley law
K-alpha
"""
k_a_energy =np.asarray([
    (res_data[0][3][0]),
    (res_data[1][3][0]),
    (res_data[2][3][0]),
    (res_data[3][3][0]),
    (res_data[4][3][0]),
    (res_data[5][3][2])
    ])
k_a_energy_sqrt = np.sqrt(k_a_energy)

k_a_energy_sqrt_error =0.5*k_a_energy_sqrt*(np.asarray([
(res_data[0][1][0]),
(res_data[1][1][0]),
(res_data[2][1][0]),
(res_data[3][1][0]),
(res_data[4][1][0]),
(res_data[5][1][2])])
 /np.asarray([(res_data[0][3][0]),(res_data[1][3][0]),(res_data[2][3][0]),
              (res_data[3][3][0]),(res_data[4][3][0]),
              (res_data[5][3][2])]))

i=0
for x,y,z in sorted(zip(k_a_energy_sqrt,k_a_energy_sqrt_error,atomic_number)):
    k_a_energy_sqrt[i]=x
    k_a_energy_sqrt_error[i]=y
    atomic_number[i]=z
    i+=1


plt.title("Moselet Law Plot:Atomic number against Square Root of Energy(k-alpha)")

plt.ylabel("Sqrt Energy(K-alpha - keV)")
plt.xlabel("Atomic Number")

 # fig.suptitle('Horizontally stacked subplots')
plt.errorbar(atomic_number,k_a_energy_sqrt,k_a_energy_sqrt_error,None,"o")
# calc the trendline
z,errorka = np.polyfit(atomic_number,k_a_energy_sqrt, 1,cov=True,full=False)
p = np.poly1d(z)
mosley_fit = p(atomic_number)

plt.plot(atomic_number,mosley_fit,"r--")
# the line equation:
print("Line Equation for K-alpha:x is atomic number and y is sqrt Energy")
print("y=%.6fx+(%.6f)"%(z[0],z[1]))
print(np.sqrt(np.diag(errorka)))
# x =int(np.around(p(np.sqrt((get_data("X")[0][0]-a[0])/b[0]))))
# print("Element X has atomic number :"+ str(x))
plt.show()
plt.close()

#Mosley Residual

m_residual = mosley_fit - k_a_energy_sqrt

plt.title("Moselet Law Residual Plot: Residual of Sqrt Energy(k-alpha) v/s Atomic Number ")
plt.axhline(y=0)
plt.ylabel("Residual")
plt.xlabel("Atomic Number")
plt.errorbar(atomic_number,m_residual,k_a_energy_sqrt_error,None,"ro")
plt.show()



"""Moseley law
K-beta
"""
k_b_energy=np.asarray([
    (res_data[0][3][1]),
    (res_data[1][3][1]),
    (res_data[2][3][1]),
    (res_data[3][3][1]),
    (res_data[4][3][1]),
    (res_data[5][3][1])])
print(k_b_energy)
k_b_energy_sqrt = np.sqrt(k_b_energy)

k_b_energy_sqrt_error =0.5*k_b_energy_sqrt*(np.asarray([(res_data[0][1][1]),(res_data[1][1][1]),(res_data[2][1][1]),
              (res_data[3][1][1]),
              (res_data[4][1][1]),(res_data[5][1][1])])/np.asarray([(res_data[0][3][1]),(res_data[1][3][1]),(res_data[2][3][1]),
              (res_data[3][3][1]),
              (res_data[4][3][1]),(res_data[5][3][1])]))

i=0
for x,y in sorted(zip(k_b_energy_sqrt,k_b_energy_sqrt_error)):
    k_b_energy_sqrt[i]=x
    k_b_energy_sqrt_error[i]=y
    i+=1


plt.title("Moselet Law Plot:Square Root of Energy(k-beta) against Atomic Number")

plt.ylabel("Sqrt Energy(K-beta)")
plt.xlabel("Atomic Number")

 # fig.suptitle('Horizontally stacked subplots')
plt.errorbar(atomic_number,k_b_energy_sqrt,k_b_energy_sqrt_error,None,"o")
# calc the trendline
z,errorkb = np.polyfit(atomic_number,k_b_energy_sqrt, 1,cov=True,full=False)
p = np.poly1d(z)
mosley_fit = p(atomic_number)

plt.plot(atomic_number,mosley_fit,"r--")
# the line equation:
print("Line Equation for K-beta:x is atomic number and y is sqrt Energy")
print("y=%.6fx+(%.6f)"%(z[0],z[1]))
print(np.sqrt(np.diag(errorkb)))
# x =int(np.around(p(np.sqrt((get_data("X")[0][0]-a[0])/b[0]))))
# print("Element X has atomic number :"+ str(x))
plt.show()
plt.close()

#Mosley Residual

m_residual = mosley_fit - k_b_energy_sqrt

plt.title("Moselet Law Residual Plot: Residual of Sqrt Energy(k-beta) v/s Atomic Number ")
plt.axhline(y=0)
plt.ylabel("Residual")
plt.xlabel("Atomic Number")
plt.errorbar(atomic_number,m_residual,k_b_energy_sqrt_error,None,"ro")
plt.show()



"""Moseley law
L
"""
atomic_number = [47,56,65]
k_l_energy=np.asarray([(res_data[0][3][2]),(res_data[1][3][3]),
              (res_data[5][3][6])])
k_l_energy_sqrt = np.sqrt(k_l_energy)

k_l_energy_sqrt_error =0.5*k_l_energy_sqrt*(np.asarray([(res_data[0][1][2]),(res_data[1][1][3]),(res_data[5][1][6])])/k_l_energy)

# i=0
# for x,y in sorted(zip(k_l_energy_sqrt,k_l_energy_sqrt_error)):
#     k_l_energy_sqrt[i]=x
#     k_l_energy_sqrt_error[i]=y
#     i+=1


plt.title("Moselet Law Plot:Atomic number against Square Root of Energy(k   -l)")

plt.ylabel("Sqrt Energy(K-l)")
plt.xlabel("Atomic Number")

 # fig.suptitle('Horizontally stacked subplots')
plt.errorbar(atomic_number,k_l_energy_sqrt,k_l_energy_sqrt_error,None,"o")
# calc the trendline
z,errorla = np.polyfit(atomic_number,k_l_energy_sqrt, 1,cov=True,full=False)
p = np.poly1d(z)
mosley_fit = p(atomic_number)

plt.plot(atomic_number,mosley_fit,"r--")
# the line equation:
print("Line Equation for K-l:x is atomic number and y is sqrt Energy")
print("y=%.6fx+(%.6f)"%(z[0],z[1]))
print(np.sqrt(np.diag(errorla)))
# x =int(np.around(p(np.sqrt((get_data("X")[0][0]-a[0])/b[0]))))
# print("Element X has atomic number :"+ str(x))
plt.show()
plt.close()

#Mosley Residual

m_residual = mosley_fit - k_l_energy_sqrt

plt.title("Moselet Law Residual Plot: Sqrt Energy(k-l Atomic Number ")
plt.axhline(y=0)
plt.ylabel("Residual")
plt.xlabel("Atomic Number")
plt.errorbar(atomic_number,m_residual,k_l_energy_sqrt_error,None,"ro")
plt.show()