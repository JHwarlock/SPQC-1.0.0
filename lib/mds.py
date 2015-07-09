"""
	  Programed by Rong ZhenqQin (rongzhenqqin@genomics.org.cn  or  zju3351689@gmail.com)
"""
import sys
import numpy as np
import time
import struct
def resampling1(wvector,numselect):
	"""
	This is a method of resampling with replacement based on the  probability created by weight matrix.
	Here, we used a Cumulative probability method to resampling.

	Input:
		wvector:	the weight vector of samples (is a nx1 list vector)
		numselect:	the number of samples resampling

	Output:
		Index:	the index of resampling samples (a list)
	"""

	wvector = map(float,wvector)
	pro = wvector/np.sum(wvector)
	procum = np.cumsum(pro)
	randseq = np.random.rand(numselect,1)
	Index = []
	for i in xrange(numselect):
		temp = np.random.rand()<procum
		temp = list(temp)
		Index.append(temp.index(True))
	return Index


def resampling2(wvector,numselect):
	"""
	This is a method of resampling without replacement based on the probability created by weight matrix.
	It is similar to resampling with replacement
	It is also could be used to sort randomly (where the wvector is ones vector, and the len(wvector) == numselect)
	
	Input:
		wvector:	the weight vector of samples (is a 1 x n list vector)
		numselect:	the number of samples resampling

	Output:
		Index:	the index of resampling samples (a list)
	"""

	wvector = map(float,wvector)
	pro = wvector/np.sum(wvector)
	procum = np.cumsum(pro)
	randseq = np.random.rand(numselect,1)
	Index = []

	def weightas0(pro,ind,procum):
		"""
		This is a sub-function to be easy to resampling.
		"""

		prosum = procum[-1]-pro[ind]
		pro[ind] = 0;
#		print pro
#		print prosum
		pro = pro/prosum;
		procum = np.cumsum(pro)
		return pro,procum

	for i in xrange(numselect-1):
		temp = np.random.rand()<procum
		temp = list(temp)
		ind = temp.index(True)
		Index.append(ind)
		pro,procum = weightas0(pro,ind,procum)
	i = numselect-1;temp = np.random.rand()<procum;temp = list(temp);ind = temp.index(True);Index.append(ind);
	#print len(Index)
	#print len(set(Index))
	return Index



class mdsoutput:
	def __init__(self):
		self.w = 0
		self.v = 0
		self.p = 0
def check_data(X_SNPs):
	if X_SNPs.dtype == np.float64 or X_SNPs.dtype == np.float32:
		xtype = X_SNPs.dtype
		pass
	else:
		sys.stderr.write("""The format of X_SNPs matrix should be numpy.float32 or numpy.float64, please check it.
		If the memory is sufficient, we suggest you use the numpy.float64\n""")
		exit(1)

import mutilstats
def mds_ps(X_raw,nvs_output=10):
	"""
	Here, a mutil demensional scale method was used to calculate the population structure.
	I think this must be improve, we will check the reliability of method for population structure analysis, soon.
	
	Input:
		nvs_output: is number of reduced demensions of raw data
		X_SNPs: is the same as it in plsgwas
	Output:
		w : a list (len(list) = nvs_output) of eigenvalue
		v : a matrix of eigenvector corresponds to the eigenvalue
	"""
	X_SNPs = X_raw.copy()
	if X_SNPs.dtype == np.float64 or X_SNPs.dtype == np.float32:
		xtype = X_SNPs.dtype
		pass
	else:
		sys.stderr.write("""The format of X_SNPs matrix should be numpy.float32 or numpy.float64, please check it.
		If the memory is sufficient, we suggest you use the numpy.float64\n""")
		exit(1)
	nx,px = X_SNPs.shape
	X_SNPs = np.asmatrix(X_SNPs)
	if nvs_output>nx:
		sys.stderr.write('too many nvs_output, it must be smaller than number of samples, we have changed auto\n')
	nvs_output = min(nx,nvs_output)
	mutilstats.centring(X_SNPs)
	mutilstats.normalize(X_SNPs)
	#print X_SNPs
	dist = np.asmatrix(np.zeros((nx,nx)))
	for i in xrange(nx):
		temp = X_SNPs - X_SNPs[i,:]
		temp = np.power(temp,2)
		dist[:,i] = np.power(np.sum(temp,axis=1),0.5)
	I = np.asmatrix(np.eye(nx))
	I_n = np.asmatrix(np.ones((nx,nx)))
	dist = -1*(I-(1.0/nx)*I_n)*dist*(I-(1.0/nx)*I_n)/2
	del I_n
	del I
	w,v=np.linalg.eig(dist)
	del dist
	idx = np.argsort(w)[::-1]
	w = w[idx]
	v = v[:,idx]
	precent = np.cumsum(w)/np.sum(w) * 100
	mds_output = mdsoutput()
	mds_output.p = precent[0:nvs_output]
	mds_output.w = w[0:nvs_output]
	mds_output.v = v[:,0:nvs_output]
	"""
	w=list(w)
	wtemp=w[:]
	wtemp.sort()
	last=-1
	vector_ind = []
	return_v = np.asmatrix(np.zeros((nx,nvs_output)))
	while nvs_output:
		vector_ind.append(w.index(wtemp[last]))
		last -= 1
		nvs_output -= 1
	return_w = []
	while vector_ind:
		ind = vector_ind.pop(0)
		return_w.append(w[ind])
		return_v[:,nvs_output] = v[:,ind]
		nvs_output += 1
	"""
	return mds_output


def randomsort(len_seq):
	wvector = np.ones(len_seq)
	ind = resampling2(wvector,len_seq)
	return ind

def dotvector(X,Xt):
	return Xt*X.transpose()


def rbfnorm2(X,Xt):
	n0 = X.shape[0]
	nx = Xt.shape[0]
	RBFnorm2 = np.zeros((nx,n0));
	for i in xrange(nx):
		for j in xrange(n0):
			RBFnorm2[i,j] = np.linalg.norm(Xt[i,:]-X[j,:])**2
	return RBFnorm2


def crossvalidate(Xtrain,Ytrain):
	"""
	Data sorted randomly for crossvalidation
	"""
	X = Xtrain.copy()
	Y = Ytrain.copy()
	len_seq=len(Y)
	ind = randomsort(len_seq)
	X = X[ind,:]
	Y = Y[ind,:]
	return X,Y

def centring(X,axis=0,Xmean = None):
	ret = None
	if Xmean == None:
		Xmean = np.mean(X,axis=axis)
	X -= Xmean
	return Xmean

def normalize(X,axis=0,Xstd=None):
	ret = None
	if Xstd == None:
		Xstd =np.std(X,ddof=1,axis=axis)
	X /= Xstd
	return Xstd

