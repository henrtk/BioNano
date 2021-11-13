import xml.etree.ElementTree as ET
from matplotlib import pyplot as plt
from math import log, pi
import numpy as np
import os
# Globals
t = 100 / 7.5 #s
pixel_size = 2.5 * 10**(-7) #m/pixel
T = 273 + 20 #Kelvin
K_b = 1.38 * 10 **(-23) #J/K
nu = 1.0016 * 10**(-3) #Pa*s
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


treeA1, rootA1 = parse_file('DF40x_7_Tracks_A1.xml')

#treeA2, rootA2 = parse_file('DF40x_7_Tracks_A2 - Copy.xml')

#treeA3, rootA3 = parse_file('DF40x_7_Tracks_A3 - Copy.xml')

def move(i,para):
    return float(i[-1].attrib[para]) -float(i[0].attrib[para])

def make_D_list(particles): #calculates D from displacement
    D_coeff_list = []
    for i in particles:
        xmove = move(i,"x") 
        ymove = move(i,"y") 
        tmove = move(i,"t") 
        displace2 = (xmove**2+ymove**2)*pixel_size**2
        D = displace2/(4*tmove*t)
        if displace2 > 16*pixel_size**2:
            D_coeff_list.append(D)
    return D_coeff_list

def newMakeDList(particles):
    dList = []
    for i in particles:
        displaceSq, N = meanDisplacement(i)
        if N > 10:
            dList.append(displaceSq*7.5/(4)) # frames
    return dList
 
def make_R_list(D_list): #calculates R from Stokes-Einstain eq and calculated Ds
    radii = []
    for elem in D_list:
        try:
            radii.append(K_b*T*10**6/(elem*6*pi*nu)) #radii in um
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

def dataDict(filenames : tuple, dfunc) -> dict:
    # create loopable structure
    fileDict = dict({})
    for file in filenames:
        tree, root = parse_file(file)
        dList = np.asarray(dfunc(root))
        rList = np.asarray(make_R_list(dList))
        fileDict[file[-6:-4]] = dict({"Dlist" :  dList,
                                      "Rlist" :  rList,
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
    for j, particle in enumerate(particles):
        data = np.zeros((len(particle),2))
        for i, detection in enumerate(particle):
            data[i] = np.asarray([detection.attrib["x"],detection.attrib["y"]])
        plt.plot(data[:,0],-data[:,1])
    plt.show()

def plotAttribs(fileDict, param = "D",cols=1):
    fig, axarr = plt.subplots(len(fileDict)//cols,cols, sharex=True)
    axarr = axarr.flatten()
    for index, sample in enumerate(fileDict.keys()):
        params = fileDict[sample][param+"list"]
        print(f"Average {param} in {sample} is {params.mean()} +- {SD(params)}")
        axarr[index].hist(fileDict[sample][param+"list"], bins=40, color = 'c',edgecolor='k', label = param+f"$_"+"{" + f"{sample}" +"}$")
        axarr[index].axvline(params.mean(), color='k', linestyle='dashed', linewidth=1)
        axarr[index].set_ylabel("Number of particles")
        axarr[index].legend()
        if param == "R":
            axarr[index].set_xlabel(param+ " $\mu m$")
        else: 
            axarr[index].set_xlabel(param + " $m^2/s$")
    plt.show()
    return

os.listdir("A kopi")

files = dataDict(('ResultofDF40x_7_Tracks.xml','DF40x_7_Tracks_A1.xml','ResultofDF40x_5-5_A2.xml'),newMakeDList)
#plotAttribs(files, param="D")
plotDrift(files["A2"]["particles"])
plotDrift(files["A1"]["particles"])
plotPaths(files["A2"]["particles"])

"""
treeA1, rootA1 = parse_file('ResultofDF40x_7_Tracks.xml')
treeA2, rootA2 = parse_file("ResultofDF40x_A1_higherT copy.xml")
treeA3, rootA3 = parse_file('DF40x_7_Tracks_A1.xml')
plotPaths(treeA2)
fig, ax = plt.subplots(3,1, sharex= False, constrained_layout = False )
figD , axD = plt.subplots(3,1 , sharex= True)

D_coeff_list_A1 = make_D_list(rootA1)
R_list_A1 = make_R_list(D_coeff_list_A1)
avg_D_A1 = avg_D_const(D_coeff_list_A1)
avg_R_A1 = avg_R(R_list_A1)
print(f'Average D in A1 is {avg_D_A1} +- {SD(D_coeff_list_A1)}')
print(f'Average radius in A1 is {avg_R_A1} +- {SD(R_list_A1)}')

D_coeff_list_A2 = newMakeDList(rootA1)
R_list_A2 = make_R_list(D_coeff_list_A2)
avg_D_A2 = avg_D_const(D_coeff_list_A2)
avg_R_A2 = avg_R(R_list_A2)
print(f'Average D in A2 is {avg_D_A2} +- {SD(D_coeff_list_A2)}')
print(f'Average radius in A2 is {avg_R_A2} +- {SD(R_list_A2)}')

D_coeff_list_A3 = make_D_list(rootA3)
R_list_A3 = make_R_list(D_coeff_list_A3)
avg_D_A3 = avg_D_const(D_coeff_list_A3)
avg_R_A3 = avg_R(R_list_A3)
print(f'Average D in A3 is {avg_D_A3} +- {SD(D_coeff_list_A3)}')
print(f'Average radius in A3 is {avg_R_A3} +- {SD(R_list_A3)}') 

totlistD = D_coeff_list_A1 + D_coeff_list_A2 + D_coeff_list_A3
totlistR = R_list_A1 + R_list_A2 + R_list_A3
'''
print(f"Total avg D in A: {avg_D_const(totlistD)} +- {SD(totlistD)} ")
print(f"total avg R in A: {avg_R(totlistR)} +- {SD(totlistR)}")

print(np.percentile(totlistR, [25, 50, 75], interpolation= 'midpoint'))
print(np.percentile(totlistD, [25, 50, 75], interpolation= 'midpoint'))

print(min(totlistD), max(totlistD))
print(min(totlistR), max(totlistR))


RA1_sort = R_list_A1.copy()
RA2_sort = R_list_A2.copy()
RA3_sort = R_list_A3.copy()

DA1_sort = D_coeff_list_A1.copy()
DA2_sort = D_coeff_list_A2.copy()
DA3_sort = D_coeff_list_A3.copy()
#husker ikke helt hvorfor jeg dreiv med den kopieringen og sorteringa
DA1_sort.sort()
DA2_sort.sort()
DA3_sort.sort()
RA1_sort.sort()
RA2_sort.sort()
RA3_sort.sort()

_, bins_r = np.histogram(RA1_sort, bins = 40)
logbins_r = np.logspace(np.log10(bins_r[0]), np.log10(bins_r[-1]), len(bins_r))
_, bins_D = np.histogram(DA1_sort, bins = 40)
logbins_D = np.logspace(np.log10(bins_D[0]), np.log10(bins_D[-1]), len(bins_D))

ax[0].hist(RA1_sort, bins = logbins_r)
ax[0].set_xscale('log')
ax[1].hist(RA2_sort, bins = logbins_r)
ax[1].set_xscale('log')
ax[2].hist(RA3_sort, bins = logbins_r)
ax[2].set_xscale('log')

axD[0].hist(DA1_sort, bins = logbins_D)
axD[0].set_xscale('log')
axD[1].hist(DA2_sort, bins = logbins_D)
axD[1].set_xscale('log')
axD[2].hist(DA3_sort, bins = logbins_D)
axD[2].set_xscale('log')

ax[0].annotate('A1', (0.5, 1.05), xycoords = 'axes fraction')
ax[2].annotate('A3', (0.5, 1.05), xycoords = 'axes fraction')
fig.text(0.05, 0.5, 'number of tracked particles', va='center', rotation='vertical')
fig.text(0.5, 0.03, 'particle radius in um', ha='center')

axD[0].annotate('A1', (0.5, 1.05), xycoords = 'axes fraction')
axD[1].annotate('A2', (0.5, 1.05), xycoords = 'axes fraction')
axD[2].annotate('A3', (0.5, 1.05), xycoords = 'axes fraction')
figD.text(0.05, 0.5, 'number of tracked particles', va='center', rotation='vertical')
figD.text(0.5, 0.03, 'Diffusion coefficient [m^2/s]', ha='center')
'''
for ax in ax.flat:
    ax.label_outer()

for axD in axD.flat:
    ax.label_outer()   
plt.show()
#fig.savefig('histR_A.jpg')
#figD.savefig('histD_A.jpg')


"""
