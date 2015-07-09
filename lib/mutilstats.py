#!/usr/bin/env python
# -*- coding: UTF-8 -*-
import compress
import os
import numpy as np
import scipy as sp
from scipy import stats
import sys
import time
import ctypes
import itertools
import statplot
class SampleInfo(object):
	def __init__(self):
		#SN Files   samplename  classid classname
		self.samplenum = 0
		
		self.classlabels = []
		self.uniqclasslabel = []
		
		self.samplecolors = []
		self.classcolors = []
		self.samplelines = []
		self.classlines = []
		self.samplemarkers = []  ## 'o','^','v','+'
		self.classmarkers = []  ##  'o','o','^','^'
		
		self.uniqcolor = []
		
		self.classnums =[] 
		self.uniqclassnum = []
		
		self.sns = []#第一列
		self.samplenames=[]#第三列
		self.traits = []#第四列
		self.uniqline = []
		self.uniqmarker = []
		
		self.files = []
		self.hidx = {}
	def parse_sampleinfo(self,sampleinfo,phenotype='category'):
		fh = file(sampleinfo,"r")
		classnumOrtraits = []
		for line in fh:
			if line.startswith("#") or line.startswith("\n") or line.startswith(" ") or line.startswith("\t"):continue
			nameid,filedetail,samplename,trait,other = line.rstrip("\n").split("\t",4)
			self.sns.append(nameid)
			files = filedetail.rstrip(",").split(",")
			for filename in files:
				self.files.append(nameid+"///"+filename)
			classname = other.split("\t")[0]
			self.samplenum += 1
			self.classlabels.append(classname)
			if classname in self.uniqclasslabel:pass
			else:self.uniqclasslabel.append(classname)
			classnumOrtraits.append(trait)
			self.samplenames.append(samplename)
		fh.close()
		## use sample number to get iteration colors, markers and lines
		self.uniqcolor,self.uniqline,self.uniqmarker = statplot.styles(len(self.uniqclasslabel))
		self.samplecolors,self.samplelines,self.samplemarkers = statplot.styles(self.samplenum)
		#self.uniqclasslabel = list(set(self.classlabels))
		self.uniqclassnum = range(len(self.uniqclasslabel))
		h={}
		tmpclasslabels = np.asarray(self.classlabels)
		for i in xrange(len(self.uniqclasslabel)):
			tmplabel = self.uniqclasslabel[i]
			self.hidx[tmplabel] = tmpclasslabels == tmplabel
			h[tmplabel] = i
		for classlabel in self.classlabels:
			self.classnums.append(h[classlabel])
			self.classcolors = [self.uniqcolor[i] for i in self.classnums]
			self.classmarkers = [self.uniqmarker[i] for i in self.classnums]
			self.classlines  = [self.uniqline[i] for i in self.classnums]
		#print self.classcolors
		#print self.classmarkers
		#print self.classlines
		
		if phenotype == "category":
			self.traits = np.transpose(np.asmatrix(np.float64(classnumOrtraits[:])))
		elif phenotype == "quantificat":
			self.traits = np.transpose(np.asmatrix(np.float64(classnumOrtraits[:])))
		else:return 1
		#print self.traits
		return 0

class MatrixAnno(object):
	##np.asarray(a[:,0].T)[0].tolist()
	def __init__(self):
		self.p = 0
		self.n = 0
		self.data  = None
		self.anno = []
		self.anno1 = []
		self.anno2 = []
	def parse_matrix_anno(self,fmatrixanno,cutoff=-10000000.0,precent=0.5,addtolog=1,log2tr=0):
		fh = compress.gz_file(fmatrixanno,"r")
		t0 = time.time()	
		sys.stderr.write('[INFO] Start to Build data ...\n')
		for line in fh:
			if line.startswith("#") or line.startswith("\n") or line.startswith(" ") or line.startswith("\t"):
				continue
			else:
				arr = line.rstrip("\n").split("\t")
				self.n = len(arr[2:])
				break
		fh.seek(0)
		t0 = time.time()
		num = int(self.n * precent)
		for line in fh:
			if line.startswith("#") or line.startswith("\n") or line.startswith(" ") or line.startswith("\t"):continue
			else:
				arr = line.rstrip("\n").split("\t")
				assert self.n == len(arr[2:])
				try:
					tmpdata = np.float64(arr[2:])
				except:
					sys.stderr.write("[ERROR] n is not same as exprsnums\n")
					print arr
					exit(1)
				if np.std(tmpdata,ddof=1) <=0:continue## filter the no var data
				if np.sum(tmpdata > cutoff) <= num: 
					#sys.stderr.write("[INFO] data filtered: %s\n"%(arr[0]+"\t"+arr[1]))
					continue
				if log2tr:
					tmpdata = np.log2(tmpdata+addtolog)
				self.p += 1
				if self.data == None:
					self.data = tmpdata
				else:
					self.data = np.concatenate((self.data,tmpdata))
				self.anno.append(arr[0] + "\t" + arr[1])
				self.anno1.append(arr[0])
				self.anno2.append(arr[1])
		self.data = np.asmatrix(np.transpose(self.data.reshape(self.p,self.n)))
		fh.close()
		assert len(self.anno) == self.p
		sys.stderr.write('[INFO] Data Built done! cost %.2fs\n'%(time.time()-t0))
		return 0

class FactorFrame(object):
	def __init__(self):
		self.fnm = []##factor name
		self.snm = []##must same as sample infos
		self.lvs = 0## number of variables
		self.var = []##
		self.levels = []
	def parse_factor(self,factorfile):
		f = compress.gz_file(factorfile,"r")
		for line in f:
			if line.startswith("##"):continue
			if line.startswith("#"):
				self.fnm = line.rstrip("\n").split("\t")[1:]
				self.lvs = len(self.fnm)
				self.levels = [0,]*self.lvs
				continue
			arr = line.rstrip("\n").split("\t")
			self.snm.append(arr[0])
			self.var.append(map(str,arr[1:]))
		f.close()
		self.var = np.asarray(self.var)
		for i in xrange(self.lvs):
			self.levels[i] = len(set(self.var[:,i].tolist()))
		print self.levels
		self.var = np.float64(self.var)
		return 0

def datacheck(X):
	ret = 1
	if isinstance(X,np.matrix) or isinstance(X,np.array):pass
	else:
		sys.stderr.write('[ERROR] centring only support matrix or array data\n')
		return ret
	if X.dtype == np.float64 or X.dtype == np.float32:pass
	else:
		sys.stderr.write('[ERROR] centring only support float64 or float32 data\n')
		return ret
	return 0

def centring(X,axis=0):
	ret = 1
	if datacheck(X):
		return ret
	Xmean = np.mean(X,axis=axis)
	X -= Xmean
	return 0

def normalize(X,axis=0):
	ret = 1
	if datacheck(X):
		return ret
	Xstd =np.std(X,ddof=1,axis=axis)
	X /= Xstd
	return 0

def quantile(Xnp):
	ranks=[]
	n, p = Xnp.shape
	Xnormalize = np.asmatrix(np.zeros((n,p)))
	for i in xrange(n):
		ranks.append(np.int32(stats.rankdata(Xnp[i,:],"min")).tolist())
	Xnptmp = Xnp.copy()
	Xnptmp.sort()
	ranks = np.asarray(ranks)
	Xnptmpmean = np.asarray(np.mean(Xnptmp,axis = 0))
	for i in xrange(n):
		Xnormalize[i] = Xnptmpmean[0][ranks[i]-1]
	return Xnormalize

def comb_replace(datalist,num):
	return list(itertools.combinations_with_replacement([datalist], num))

