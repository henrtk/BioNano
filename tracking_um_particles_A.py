import xml.etree.ElementTree as ET
from matplotlib import pyplot as plt
from math import log, pi
import numpy as np
import os
import seaborn as sns
from seaborn.palettes import color_palette
# Globals
t = 100 / 7.5 #s
pixel_size = 2.5 * 10**(-7) #m/pixel
T = 273 + 20 #Kelvin
K_b = 1.38 * 10 **(-23) #J/K
nu = 1.0016 * 10**(-3) #N/m^2 *s
plt.style.use("seaborn")
# Globals

def parse_file(file):
    tree = ET.parse(file)
    root = tree.getroot()
    return tree, root

def remove_noise(file): #removes particles tracked in 10 or less frames
    tree, root = parse_file(file)
    for particle in root.findall('particle'):
        if int(particle.attrib['nSpots'] ) <= 10:
            root.remove(particle)
    tree.write(file)

def move(i,para):
    return float(i[-1].attrib[para]) -float(i[0].attrib[para])

def make_D_list(particles): #calculates D from displacement
    D_coeff_list = []
    for i in particles:
        xmove = move(i,"x") 
        ymove = move(i,"y") 
        tmove = move(i,"t") 
        displace2 = (xmove**2+ymove**2)*pixel_size**2
        D = 10**12*displace2/(4*tmove*t)
        if displace2 > 16*pixel_size**2:
            D_coeff_list.append(D)
    return D_coeff_list

def newMakeDList(particles):
    dList = []
    for i in particles:
        displaceSq, N = meanDisplacement(i)
        D = displaceSq*7.5/(4)
        if N > 10 and D>10**(-13):
            dList.append(D) # displace/(4t)
    return dList
 
def make_R_list(D_list): #calculates R from Stokes-Einstain eq and calculated Ds
    radii = []
    for elem in D_list:
        try:
            radii.append(K_b*T/(elem*6*pi*nu)) #radii in um
        except:
            ZeroDivisionError  
    return radii

def avg_D_const(D_list):
    return sum(D_list)/len(D_list)

def avg_R(R_list): #lagde vist to funksjoner for å regne gjennomsnitt:)
    return sum(R_list)/len(R_list)

def SD(datalist):
    calc_list = []
    avg = avg_D_const(datalist)

    for data in datalist:
        calc_list.append((data-avg)**2)
    return np.sqrt(sum(calc_list)/(len(datalist)-1))

def dataDict(filenames : tuple, dlistfunc) -> dict:
    # create loopable structure
    fileDict = dict({})
    for file in filenames:
        tree, root = parse_file(file)
        print(type(root))
        dList = np.asarray(dlistfunc(root))
        rList = np.asarray(make_R_list(dList))
        #filenaming should be more general....
        fileDict[file.name[23:25]] = dict({"Dlist" :  dList*10**(12),
                                      "Rlist" :  rList*10**(6),
                                      "avgD"  :  avg_D_const(dList),
                                      "avgR"  :  avg_R(rList),
                                      "particles"  :  root})
    return fileDict     

def meanDisplacement(particle):
    N = len(particle)
    tot = 0
    for i in range(N-1):
        tot += ((float(particle[i+1].attrib["y"]) - float(particle[i].attrib["y"]))*pixel_size)**2
        tot += ((float(particle[i+1].attrib["x"]) - float(particle[i].attrib["x"]))*pixel_size)**2
    return tot/N , N

def meanDrift(particle, param):
    N = len(particle)
    drift = (float(particle[-1].attrib[param]) - float(particle[0].attrib[param]))/N   
    return drift, N
     
def calculateMeanDrift(particles):
    data = []
    for particle in particles:
        x, N = meanDrift(particle,"x")
        y, _ = meanDrift(particle,"y")
        if N>10:
            data.append(np.asarray([x, y]))
    return np.asarray(data)

def plotDrift(particles):
    drift = calculateMeanDrift(particles)
    plt.scatter(drift[:,0],-drift[:,1], s =1 ) 
    plt.scatter(drift[:,0].mean(),drift[:,1].mean(), color = "r", label = "Mean drift")
    plt.legend()
    plt.show()
    return

def plotPaths(particles):
    plt.xlim(0,1360)
    plt.ylim(-1024,0)
    for particle in particles:
        data = np.zeros((len(particle),2))
        for i, detection in enumerate(particle):
            data[i] = np.asarray([detection.attrib["x"],detection.attrib["y"]])
        plt.plot(data[:,0],-data[:,1])
    plt.show()

def plotColours():
    x = np.linspace(0,2*np.pi*2,200)
    for i, el in enumerate(sns.color_palette("muted")):
        plt.plot(x,np.sin(x+i*np.pi/6), color = el,label = str(i) )
    plt.legend()
    plt.show()

def plotAttribs(fileDict, param = "D",cols=1):
    fig, axarr = plt.subplots(len(fileDict)//cols,cols, sharex=True)
    axarr = axarr.flatten()
    colors = sns.color_palette("muted")
    fig.subplots_adjust(right = 0.6)
    for index, sample in enumerate(fileDict.keys()):
        params = fileDict[sample][param+"list"]
        print(f"Average {param} in {sample} is {params.mean()} +- {SD(params)}")
        axarr[index].hist(fileDict[sample][param+"list"], bins=40, label = f" {sample}", color = colors[0])
        axarr[index].axvline(params.mean(), color='k', linestyle='dashed', linewidth=1, label = "$\langle D_{} \\rangle $")
        #axarr[index].set_ylabel("Particle counts")
        handles, labels = axarr[index].get_legend_handles_labels()
        axarr[index].legend(fontsize=14)
        labels, handles = zip(*sorted(zip(labels, handles), key=lambda t: t[0][1]))
        axarr[index].legend(handles, labels)
        #if param == "R":
            #axarr[index].set_xlabel("Radius " + " $\mathrm{\mu m}$")
        #else: 
            #pass#axarr[index].set_xlabel("Diff. coeff. " + " $[\mathrm{\mu m^2/s}]$")
    plt.figtext(0.3,0.01,"Particle radius $[\mathrm{\mu m}]$", wrap = True, fontsize = 14,horizontalalignment = "center")
    plt.figtext(0.07,0.38,"Particle counts",rotation = "vertical", wrap = True, fontsize = 14)
    plt.show()
    return
def plotall(filedict,param,title):
    tot = np.asarray([0.93])
    for i in filedict.keys():
        tot = np.concatenate((tot,filedict[i][param+"list"]))
    fig, ax = plt.subplots(1)
    ax.hist(tot,bins=40, color = color_palette("muted")[0],label = title)
    ax.axvline(tot.mean(), color='k', linestyle='dashed', linewidth=1, label = "$\langle r_{tot} \\rangle $")
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(fontsize=14)
    labels, handles = zip(*sorted(zip(labels, handles), key=lambda t: t[0][1]))
    ax.legend(handles, labels)
    plt.show()

