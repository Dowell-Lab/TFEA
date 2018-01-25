import numpy as np
import time

a = [i for i in range(30000)]

numpytime = time.time()
simES = list()
for i in range(1000):
	Eval = 0.0
	ES = list()
	np.random.shuffle(a)
	for distance in a:
		if distance < 0.5:
			Eval += distance
			ES.append(Eval)
		else:
			Eval += -distance
			ES.append(Eval)
	simES.append(max(ES,key=abs))
print simES[:10]
print "Numpy: ", time.time()-numpytime

othertime = time.time()
simES = list()
for i in range(1000):
        Eval = 0.0
        maximum = 0.0
	minimum = 0.0
        np.random.shuffle(a)
        for distance in a:
                if distance < 0.5:
                        Eval += distance
                        if Eval > maximum:
				maximum = Eval
                else:
                        Eval += -distance
                        if Eval < minimum:
				minimum = Eval
        simES.append(max(maximum,minimum,key=abs))
print simES[:10]
print "Other time: ", time.time()-othertime
