from scipy import signal
from scipy.fftpack import fft
from pylab import *
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab

mu,sigma = 100,15
data = []
x = []
y = []
normalized_data = []
rrtotal = rrmean = rrnum = rr_rrmean = rr_rrmean_total =  0.0
sdnn_base = 0
sdnn_index_base = 0.0
sdnn_index = 0.0
rmssd_base = 0
sdann_index = 0.0
sdann_loop = 0
sdann_base = 0.0
sdann_tot = 0.0
n = 0
num = 0
nn50 = 0.0
pnn50 = 0.0
nn20 = 0.0
pnn20 =0.0
vlf_num = 0.0
lf_num = 0.0
hf_num = 0.0
sample_rate =250                                          #sample rate
time_step = 1/sample_rate                                 #time step 
ulfpower_freq_rate = 0.003/(sample_rate/2)
vlfpower_freq_rate = 0.04/(sample_rate/2) #0.04/sample_rate
lfpower_freq_rate = 0.15/(sample_rate/2) #0.15/sample_rate
hfpower_freq_rate = 0.4/(sample_rate/2) #0.4/sample_rate
ulf_power = 0.0
vlf_power = 0.0
lf_power = 0.0
hf_power = 0.0
tot_power = 0.0      
rr_max = 0.0
rr_min = 0.0


f = open('chf03.txt',"r")

for l in f:
	rr = int(l)
	rrtotal += rr	
	data.append(float(l))

num = len(data)                                            #sample size
rrmean = rrtotal/num
#print data
#print len(data)
#print rrtotal



############################# time domain ####################################
for i in range( 0, len(data)-2,1) :                        #time domain
	
	sdnn_base += np.square(data[i] - rrmean)           #sdnn
	rmssd_base += np.square(data[i+1] - data[i])       #rmssd
	if np.abs(data[i] - data[i+1]) > 50 :              #nn50
		nn50 += 1
	if np.abs(data[i] - data[i+1]) > 20 :              #nn20
		nn20 += 1
nn20 -= nn50
sdnn_base += np.square(data[len(data)-1] - rrmean)
pnn50 = nn50/(num-1) * 100                                 #pnn50
pnn20 = nn20/(num-1) * 100                                 #pnn20

sdann_loop = num / (5*60*sample_rate+1)                    #sdann index
for i in range(0,sdann_loop,1) :
	n += 1
	for j in range( i *(5*60*sample_rate+1) , (i+1)*(5*60*sample_rate+1) ,1) :
		sdann_base += data[j]
		sdnn_index_base += np.square(data[j] - rrmean)
	sdann_tot += sdann_base / (5*60*sample_rate)
	sdnn_index += np.sqrt(sdnn_index_base/(5*60*sample_rate-1))


sdann_index = sdann_tot / n
sdnn_index = sdnn_index / n
######################### normalize ##########################################
rr_max = rr_min = data[0]
for i in range(0,len(data),1) :
	
	if data[i] > rr_max :
		rr_max = data[i]
	if data[i] < rr_min :
		rr_min = data[i]

print "rrmax= ",rr_max,"rrmin= ",rr_min
for i in range(0, len(data), 1) :
	normalized_data.append((data[i] - rr_min)/(rr_max - rr_min))
#	print data[i]

######################### FREQ ################################################
fft = np.fft.fft(normalized_data)                                         #freq tomain
amplitude = np.abs(fft) #**2 #/1000000                          #ampitude  ms^2 => s^2
#amplitude = 10*np.log10(np.abs(fft))
#print amplitude
#print "amplitude[14] = ",amplitude[14]
freq = np.fft.fftfreq(fft.size, 0.004)                            #frequency
#print freq, len(freq)
#print freq[0:fft.size/2]

ulf_num = int(num/2 * ulfpower_freq_rate)
vlf_num = int(num/2 * vlfpower_freq_rate)
lf_num = int(num/2 * lfpower_freq_rate)
hf_num = int(num/2 * hfpower_freq_rate)
#print vlf_num,lf_num,hf_num

if ulf_num < 1 :
	ulf_num = 2
for i in range(1,ulf_num,1) :	                             #ulf power
	ulf_power += amplitude[i] * freq[i]

for i in range(1, vlf_num-1, 1) :                            #vlf power
	vlf_power += amplitude[i]  * freq[i]
	#print amplitude[i]
for i in range(vlf_num-1, lf_num, 1) :                       #lf power	
	#print amplitude[i] 
	lf_power += amplitude[i]  * freq[i] 

for i in range(lf_num+1, hf_num-1, 1) :                        #hf power	
	hf_power += amplitude[i]   * freq[i]
	#print amplitude[i]
for i in range(0, hf_num-1, 1) :                             #total power
	tot_power += amplitude[i]   * freq[i]

#################### nonlinear #################################################

print "sample number: ",num
print "//////////////////////////time domain/////////////////////////"
print "mean RR= ",rrmean,"ms"
print "SDNN=",np.sqrt(sdnn_base/(len(data)-1)),"ms"
print "SDNN index= ",sdnn_index,"ms"
print "SDANN= ",rrmean,"ms"
print "SDANN index= ",sdann_index,"ms"
print "RMSSD= ",np.sqrt(rmssd_base/(len(data)-1)),"ms"
print "NN50= ",nn50
print "pNN50= ",pnn50,"%"
print "NN20= ",nn20
print "pNN20= ",pnn20,"%"
print "TINN= "
print "Trangular index= "
print "//////////////////////////freq domain/////////////////////////"
print "ulf power= ",ulf_power
print "vlf power = ",vlf_power
print "lf power= ",lf_power
print "hf power= ",hf_power
print "tot_power = ",tot_power
print "lf/hf = ",lf_power/hf_power
print " lf % = ",(vlf_power+lf_power)/tot_power *100,"%"
print " hf % = ",hf_power/tot_power *100,"%"
print "nlf=",lf_power / (tot_power - vlf_power)
print "nhf=",hf_power / (tot_power - vlf_power)
print "////////////////////////nonlear/////////////////////////////////"
print "SD1= ",np.sqrt(0.5*(rmssd_base/(len(data)-1)))
print "SD2= ",np.sqrt(2*(sdnn_base/(len(data)-1)) - 0.5*(rmssd_base/(len(data)-1)))
#for i in range(0,len(data)/2+1,1) :
#	print freq[i]
#print vlf_num,lf_num,hf_num
#for lf in amplitude([0,1000])
#	lf += lf
#print lf


# db = 10 * np.log10(power)

#fftx = np.linspace(0.0 , 125.0 , 40000, endpoint=True, retstep=True)

############################ plot fft ###################################
#plt.plot(freq[0:fft.size/2],amplitude[0:fft.size/2],'r')
#plt.plot(freq[0:vlf_num],amplitude[0:vlf_num],'r')
#plt.plot(freq[vlf_num:lf_num],amplitude[vlf_num:lf_num],'m')
#plt.plot(freq[lf_num:hf_num],amplitude[lf_num:hf_num],'y')
#plt.plot(freq[hf_num:fft.size/2],amplitude[hf_num:fft.size/2])
#plt.xlim([0,0.5])
#plt.ylim([0,300])
#xlabel('Freqency(Hz)')
#ylabel('Amplitude ms^2')
#plt.show()


############################plot a histogram ##########################
#a, bins, patches = plt.hist(data, 70, normed=0, facecolor='green',alpha=0.75)
#y = mlab.normpdf( bins, mu, sigma)
#plt.axis([0,1100,0,1000])
#plt.xlabel('RR time(ms)')
#plt.ylabel('number of intervals')
#plt.grid(True)
#plt.show()

########################### nonlinear #####################################3
for i in range(0,len(data),1) :
	
		if (i%2 != 0) :                #odd  1.3.5.7.9....
			y.append(int(data[i]))
		else :                         #even 0.2.4.6.8....
			x.append(int(data[i]))

#print len(x)
#print len(y)
plt.plot(x,y,"c*")
plt.show()








#plt.plot(amplitude[0:39999])
#plt.plot(fftx,amplitude[0:39999])
#plt.xlim([0,1000])

#plt.xlim([0,0.5])
#plt.ylim([-0.5, 0.5])

#signal = np.array([int(x,10) for x in f.read().split()], dtype=int)
#sample numbers n
#n = signal.size

#print signal
#print "n=",n

#the frenquency via fft is a complex number
#fourier = np.fft.fft(signal)

#print "fft = ",fourier
#
#amplitude = np.abs(fourier)
#power = np.square(np.abs(fourier))


#print power

#print type(power[0])

#plot(power[0:n/2])

#print "amplitude = ",amplitude
#print "power = ",power
#print "------------------------------------------------"
#timestep = 0.4
#frenquency = np.fft.fftfreq(n, d=timestep) #sample spacing
#frenquency = np.log(frenquency)
#print "freqency = ",frenquency
#print type(frenquency[0])



#fftx = np.linspace(0.0,0.5,n)

#subplot(313)
#plot(fftx,power)

#plot(frenquency[0:len(frenquency)/2])

#xlabel('freqency(Hz)')
#ylabel('power ms2')
#title('frequency domain graph')

#plt.show()

