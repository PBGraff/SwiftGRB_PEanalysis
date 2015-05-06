import sys
import pickle
import numpy as np
from sklearn import tree,ensemble

fp1 = open('AdaBoost_best.pkl','r')
AB = pickle.load(fp1)
fp1.close()

data = np.loadtxt(sys.argv[1])

ypred1 = AB.predict_proba(data)

np.savetxt(sys.argv[2], ypred1[:,1], fmt='%0.9lf')
