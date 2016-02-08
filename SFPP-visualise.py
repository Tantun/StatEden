import random
import pylab
import math


file_name = 'results_50kx100k-uniform'    # File to get the data from

def CDF(seq):    # Computing the tail distribution (1 - cummulative distribution function) of the sorted sample
    n = len(seq)
    cdf = []
    i = 0
    j = 0
    while i <= seq[-1]:
        value = n - 0.0 - j
        cdf.append(value/n)
        while j < n and seq[j] <= i:
            j += 1
        i += 1
    return cdf[:]
    




Heights, LeftDs, RightDs, Widths, Disps = [], [], [], [], []
data_file = open(file_name, 'r')   # Get and extract the data from file_name
for line in data_file:
    line_parsed = line.split()
    Heights.append(int(line_parsed[0]))
    LeftDs.append(int(line_parsed[1]))
    RightDs.append(int(line_parsed[2]))
    Widths.append(int(line_parsed[3]))
    Disps.append(max((-1)*LeftDs[-1], RightDs[-1]))
data_file.close()


# Computing tail distributions for various statistics, processed for inferring the exponent gamma in P(X>t) ~ t^(-gamma); We compute both log P(X>t) and log P(X>t)/log t

Heights.sort()    
Heights_cdf = CDF(Heights)
Heights_cdf_log = []    
for num in Heights_cdf:
    Heights_cdf_log.append(math.log(num))
log_xh = list(math.log(i+1) for i in range(len(Heights_cdf)))

Heights_cdf_log_exp = []
for i in xrange(1,len(log_xh)):
    Heights_cdf_log_exp.append(Heights_cdf_log[i]/log_xh[i])
    
        
Widths.sort()
Widths_cdf = CDF(Widths)[1:]
Widths_cdf_log = []
for num in Widths_cdf:
    Widths_cdf_log.append(math.log(num))
log_xw = list(math.log(i+1) for i in range(len(Widths_cdf)))

Widths_cdf_log_exp = []
for i in xrange(1,len(log_xw)):
    Widths_cdf_log_exp.append(Widths_cdf_log[i]/log_xw[i])


RightDs.sort()
RightDs_cdf = CDF(RightDs)
RightDs_cdf_log = []
for num in RightDs_cdf:
    RightDs_cdf_log.append(math.log(num))
log_xrd = list(math.log(i+1) for i in range(len(RightDs_cdf)))

RightDs_cdf_log_exp = []
for i in xrange(1,len(log_xrd)):
    RightDs_cdf_log_exp.append(RightDs_cdf_log[i]/log_xrd[i])


        
Disps.sort()
Disps_cdf = CDF(Disps)
Disps_cdf_log = []
for num in Disps_cdf:
    Disps_cdf_log.append(math.log(num))
log_xd = list(math.log(i+1) for i in range(len(Disps_cdf)))

Disps_cdf_log_exp = []
for i in xrange(1,len(log_xd)):
    Disps_cdf_log_exp.append(Disps_cdf_log[i]/log_xd[i])


    
# Plotting tail distributions for various statistics and the conjectured exponent values P(X>t) ~ t^(-1) when X = Width and X = Displacement and P(X>t) ~ t^(-2/3) when X = Height
    
    
Fig_H = pylab.figure('Heights')
pylab.plot(log_xh, Heights_cdf_log, label = 'Heights')
pylab.plot([0,log_xh[-1]], [0,(-2)*log_xh[-1]/3.0], label = "-2/3")
pylab.legend()
Fig_H.show()

        
Fig_H = pylab.figure('Heights exponents')
pylab.plot(log_xh[1:], Heights_cdf_log_exp, label = 'Height exponents')
pylab.plot([0,log_xh[-1]], [-2.0/3.0,-2.0/3.0], label = "-2/3")
pylab.legend()
Fig_H.show()


Fig_W = pylab.figure('Widths')
pylab.plot(log_xw, Widths_cdf_log, label = 'Widths')
pylab.plot(log_xd, Disps_cdf_log, label = 'Displacement')
pylab.plot([0,max(log_xw[-1],log_xd[-1])], [0,(-1)*max(log_xw[-1], log_xd[-1])], label = "-1")
pylab.legend()
Fig_W.show()

Fig_W = pylab.figure('Widths exponents')
pylab.plot(log_xw[1:], Widths_cdf_log_exp, label = 'Width exponents')
pylab.plot(log_xd[1:], Disps_cdf_log_exp, label = 'Displacement exponents')
pylab.plot([0,max(log_xw[-1],log_xd[-1])], [-1,-1], label = "-1")
pylab.legend()
Fig_W.show()

    
    

