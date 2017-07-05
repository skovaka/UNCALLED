import numpy as np                               

readfile = open("events_40525_ch264_read510.txt")
readfile.readline()                              
read_means = list()                              
read_stdvs = list()                              
for line in readfile:
    i, mean, stdv = map(float, line.strip().split())
    read_means.append(mean)
    read_stdvs.append(stdv)

modelfile = open("r9_250bps.nucleotide.6mer.template.model")
model_means = list()
model_stdvs = list()
for line in modelfile:
    if line[0] == "#" or line[:4] == "kmer":
        continue
    mean, stdv = line.split()[3:5]
    model_means.append(float(mean))
    model_stdvs.append(float(stdv))

scale = np.std(read_stdvs) / np.std(model_means)
shift = np.mean(read_stdvs) - np.mean(model_means)*scale

print "Events:", np.mean(read_stdvs), np.std(read_stdvs), np.var(read_stdvs)
print "Model:", np.mean(model_means), np.std(model_means), np.var(model_means) 
print "Scale/Shift:",scale, shift

scaled_means = list(map(lambda m: m*scale, model_means))
print "Adjusted:",np.mean(scaled_means), np.std(scaled_means), np.var(scaled_means) 

