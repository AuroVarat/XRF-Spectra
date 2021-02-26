import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks
from scipy.signal._peak_finding import peak_prominences
import pylab as plb
from scipy.optimize import curve_fit
from scipy import asarray as ar,exp
b = 1/12.4727649
a = 0.190193385

def get_data(file_handle,element_no):
    count =[]
    channel=[]
    for line in file_handle:
        linedata = line.split()
        channel.append(float(linedata[0]))
        count.append(float(linedata[element_no+1]))
    count = np.asarray(count)
    return count,channel

"""
    Estimating Ka and Kb peaks
"""
file = open("data2.txt",'r')
file.readline() #remove the title column
filedata = file.readlines()

elements_token = ["Ag","Ba","Cu","Mo","X","Tb"]
element_no_input = str(input("Specify element for Specific energy estimate?"))
# elements =[]
# for element in elements:
#     element = np.zeros(len(filedata))
#     elements.append(element)
# Counts and Channels for all elements
counts = []
channels =[]
for number,element in enumerate(elements_token):
    count,channel = get_data(filedata,number)
    counts.append(count)
    channels.append(channel)

if 0 <= elements_token.index(element_no_input) <= 5:
    
    element_no = elements_token.index(element_no_input)
    #Estimate Peak energies for provided element
    peaks = find_peaks(counts[element_no],height=0)
    prominence = peaks[1]["peak_heights"]
    print(len( peaks[1]["peak_heights"]))
    sorlist = np.sort(prominence)
    first_max = np.where(prominence == sorlist[-1])[0][0]
    second_max = int(np.where(prominence == sorlist[-2])[0][0])

    if element_no == 5:
        second_max = int(np.where(prominence == sorlist[-4])[0][0])
    sortcount = sorted(count)
    ka_channel = channel[int(peaks[0][first_max])]
    kb_channel = channel[int(peaks[0][second_max])]
    print("Energies of Element:"+elements_token[element_no])
    print("ka_energy "+"= "+str((ka_channel-a)*b))
    print("kb_energy "+"= "+str((kb_channel-a)*b))
"""
    Gaussian Fit
"""


# plt.plot(channels[element_no],counts[element_no],".")
# for xc in peaks[0]:
#     plt.axvline(x=xc,color ='r',linewidth=1,linestyle =':')
# plt.show()
# plt.close()
# start = int(input())
# end = int(input())
# channel_selected = channels[element_no]
# x = channel_selected[start:end]
# y= counts[element_no][start:end]
# n = len(x)                          #the number of data
# mean = sum(x*y)/n                   #note this correction
# sigma = sum(y*(x-mean)**2)/n        #note this correction

# def gaus(x,a,x0,sigma):
#     return a*exp(-(x-x0)**2/(2*sigma**2))

# popt,pcov = curve_fit(gaus,x,y,p0=[1,mean,sigma])

# plt.plot(x,y,'b+:',label='data')
# plt.plot(x,gaus(x,*popt),'ro:',label='fit')
# plt.legend()
# plt.title('Fig. 3 - Fit for Time Constant')
# plt.xlabel('Time (s)')
# plt.ylabel('Voltage (V)')
# plt.show()



"""
    Comparsion of Auxillary Peaks
"""
plt.close()
show_compare = input("Show Comparsion?")
if show_compare == "Y":
        # fig, axes = plt.subplots(len(elements))
        plt.title("Count v/s Channel Number")

        plt.ylabel("Count")
        plt.xlabel("Channel")
        plt.yticks(np.arange(0))
        plt.legend(elements_token)
        # fig.suptitle('Horizontally stacked subplots')
        for number,element in enumerate(elements_token):
                
                plt.plot(channels[number], counts[number]+number*1500,label =element)
        plt.legend(loc="upper right")
        plt.show()

