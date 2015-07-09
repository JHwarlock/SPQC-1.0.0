# -*- coding: UTF-8 -*-
import sys
import numpy as np
import scipy as sp
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from  utils import dirDetectCreate 
import matplotlib as mpl
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm
from matplotlib import font_manager as fm
import itertools

#print cm._cmapnames
#['Spectral', 'copper', 'RdYlGn', 'Set2', 'summer', 'spring', 'gist_ncar', 'terrain', 'OrRd', 'RdBu', 'autumn', 'gist_earth', 'Set1', 'PuBu', 'Set3', 'brg', 'gnuplot2', 'gist_rainbow', 'pink', 'binary', 'winter', 'jet', 'BuPu', 'Dark2', 'prism', 'Oranges', 'gist_yarg', 'BuGn', 'hot', 'PiYG', 'YlOrBr', 'PRGn', 'Reds', 'spectral', 'bwr', 'RdPu', 'cubehelix', 'Greens', 'rainbow', 'Accent', 'gist_heat', 'YlGnBu', 'RdYlBu', 'Paired', 'flag', 'hsv', 'BrBG', 'seismic', 'Blues', 'Purples', 'cool', 'Pastel2', 'gray', 'coolwarm', 'Pastel1', 'gist_stern', 'gnuplot', 'GnBu', 'YlGn', 'Greys', 'RdGy', 'ocean', 'YlOrRd', 'PuOr', 'PuRd', 'gist_gray', 'CMRmap', 'PuBuGn', 'afmhot', 'bone']

# for projection='3d'
from mpl_toolkits.mplot3d import Axes3D

###动态设置字体
#mpl.rcParams['font.family'] = 'sans-serif'
#mpl.rcParams['font.sans-serif'] = ['Arial']
# font size for figures
mpl.rcParams.update({'font.size': 10})
mpl.rcParams['legend.fontsize'] = 9
# Arial font
#mpl.rc('font',**{'family':'sans-serif','sans-serif':['Arial']})

###
from matplotlib.patches import Polygon
# to get kmeans and scipy.cluster.hierarchy
from scipy.cluster.vq import *
from scipy.cluster.hierarchy import *

###
from matplotlib.colors import LogNorm

##kmeans归一化处理 from scipy.cluster.vq import whiten
from scipy.cluster.vq import whiten
mpl.style.use('ggplot')

def styles(num,cm=0):
	#	plot(x, y, color='green', linestyle='dashed', marker='o',
	#			         markerfacecolor='blue', markersize=12).
	# 因此，产生3个元素列表，for  colorstyle， linestyle， makerstyle  ##
	#1247 'b'         blue
	#1248 'g'         green
	#1249 'r'         red
	#1250 'c'         cyan
	#1251 'm'         magenta
	#1252 'y'         yellow
	#1253 'k'         black
	#1271     ================    ===============================
	#1272     character           description
	#1273     ================    ===============================
	#1274     ``'-'``             solid line style
	#1275     ``'--'``            dashed line style
	#1276     ``'-.'``            dash-dot line style
	#1277     ``':'``             dotted line style
	#1278     
	#1279     
	#1280     ``'.'``             point marker
	#1281     ``','``             pixel marker
	#1282     ``'o'``             circle marker
	#1283     ``'v'``             triangle_down marker
	#1284     ``'^'``             triangle_up marker
	#1285     ``'<'``             triangle_left marker
	#1286     ``'>'``             triangle_right marker
	#1287     ``'1'``             tri_down marker
	#1288     ``'2'``             tri_up marker
	#1289     ``'3'``             tri_left marker
	#1290     ``'4'``             tri_right marker
	#1291     ``'s'``             square marker
	#1292     ``'p'``             pentagon marker
	#1293     ``'*'``             star marker
	#1294     ``'h'``             hexagon1 marker
	#1295     ``'H'``             hexagon2 marker
	#1296     ``'+'``             plus marker
	#1297     ``'x'``             x marker
	#1298     ``'D'``             diamond marker
	#1299     ``'d'``             thin_diamond marker
	#1300     ``'|'``             vline marker
	#1301     ``'_'``             hline marker
	#1302 
	#1303 marker: [ ``7`` | ``4`` | ``5`` | ``6`` | ``'o'`` | ``'D'`` | ``'h'`` | ``'H'`` | ``'_'`` | ``''`` | ``'None'`` | ``' '`` | ``None`` | ``'8'`` | ``'p'``      | ``','`` | ``'+'`` | ``'.'`` | ``'s'`` | ``'*'`` | ``'d'`` | ``3`` | ``0`` | ``1`` | ``2`` | ``'1'`` | ``'3'`` | ``'4'`` | ``'2'`` | ``'v'`` | ``'<'``       | ``'>'`` | ``'^'`` | ``'|'`` | ``'x'`` | ``'$...$'`` | *tuple* | *Nx2 array* ]
	
	#color_raw = itertools.cycle(itertools.chain(['b','r','g','c','m','y','k']))
	color_raw = itertools.cycle(itertools.chain(plt.rcParams['axes.color_cycle']))
	lines_raw = itertools.cycle(itertools.chain(['-','--','-.',':']))
	marker_raw = itertools.cycle(itertools.chain(['o','^','s','*','+','D','v','1','2','x','3','4','s','p','>','<','h','d','H']))
	ret_color = []
	ret_lines = []
	ret_marker = []
	for i in xrange(num):
		ret_color.append(color_raw.next())
		ret_lines.append(lines_raw.next())
		ret_marker.append(marker_raw.next())
	return ret_color,ret_lines,ret_marker
def color_grad(num,colorgrad=cm.Set3):
	color_class = cm.Set3(np.linspace(0, 1, num))

def test_iter(num):
	fig = plt.figure(dpi=300)
	x = 1
	y = 1
	ax = fig.add_subplot(111)
	ret_color,ret_lines,ret_marker = styles(num)
	for i in xrange(num):
		ax.plot([x,x+1,x+2,x+3,x+4],[y,y,y,y,y],color=ret_color[i],linestyle=ret_lines[i],marker=ret_marker[i],markeredgecolor=ret_color[i],markersize=12,alpha=0.8)
		y += 1
	plt.savefig("test_style.png",format='png',dpi=300)	
	plt.clf()
	plt.close()
	return 0

def plot_enrich(resultmark,resultothers,fig_prefix,xlabel,ylabel):
	fig = plt.figure(figsize=(8,6),dpi=300)
	num = len(resultmark) + 1
	ret_color,ret_lines,ret_marker = styles(num)
	ax = fig.add_subplot(111)
	for i in xrange(num-1):
		#ax.plot(resultmark[i][1],resultmark[i][2],ret_color[i]+ret_marker[i],label=resultmark[i][0],markeredgecolor=ret_color[i],markersize=8,alpha=0.7)
		ax.plot(resultmark[i][1],resultmark[i][2],color=ret_color[i],linestyle='',marker=ret_marker[i],label=resultmark[i][0],markeredgecolor=ret_color[i],markersize=10,alpha=0.7)
	xarr = []
	yarr = []
	for ret in resultothers:
		xarr.append(ret[0])
		yarr.append(ret[1])
	ax.plot(xarr,yarr,'ko',label="others",markeredgecolor='k',markersize=3,alpha=0.5)
	art = []
	lgd = ax.legend(bbox_to_anchor=(1.02, 1),loc=0,borderaxespad=0,numpoints=1,fontsize=6)
	art.append(lgd)
	ax.set_xlabel(xlabel)
	ax.set_ylabel(ylabel)
	plt.savefig(fig_prefix+".png",format='png',additional_artists=art,bbox_inches="tight",dpi=300)
	plt.savefig(fig_prefix+".svg",format='svg',additional_artists=art,bbox_inches="tight",dpi=300)
	plt.clf()
	plt.close()
	return 0

def kdensity(var_arr,num = 200,fun='pdf',cdfstart=-np.inf):
	"""
	plot theory distribution
	y = P.normpdf( bins, mu, sigma)
	l = P.plot(bins, y, 'k--', linewidth=1.5)
	"""
	if fun not in ['cdf','pdf']:
		sys.stderr.write("kdensity Fun should be 'cdf' or 'pdf'")
		exit(1)
	kden = stats.gaussian_kde(var_arr)
	#kden.covariance_factor = lambda : .25	
	#kden._compute_covariance()
	min_a = np.min(var_arr)
	max_a = np.max(var_arr)
	xnew = np.linspace(min_a, max_a, num)
	if fun == 'cdf':
		ynew = np.zeros(num)
		ynew[0] = kden.integrate_box_1d(cdfstart,xnew[0])
		for i in xrange(1,num):
			ynew[i] = kden.integrate_box_1d(cdfstart,xnew[i])
	else: ynew = kden(xnew)
	return xnew,ynew
def hcluster(Xnp,samplenames,fig_prefix):
	linkage_matrix = linkage(Xnp,'ward','euclidean')
	fig = plt.figure(dpi=300)
	ax = fig.add_subplot(111)
	#dendrogram(linkage_matrix,labels=samplenames,leaf_label_rotation=45) ## new version of scipy
	dendrogram(linkage_matrix,labels=samplenames,orientation='right')
	ax.grid(visible=False)
	plt.savefig(fig_prefix+"_hcluster.png",format='png',dpi=300)
	plt.savefig(fig_prefix+"_hcluster.svg",format='svg',dpi=300)
	plt.clf()
	plt.close()
	return 0

def plot_hmc_curve(X,Y,colors,classlabels,figname_prefix="out",scale=0):
	#调和曲线生成Harmonic curve
	#X = n x p   Y is list, colors is list
	n,p = X.shape
	if n == len(Y) and len(Y) == len(colors):pass
	else: return 1
	
	if scale ==1:
		X = whiten(X)
	step = 100
	t = np.linspace(-np.pi, np.pi, num=step)
	f = np.zeros((n,step))
	for i in xrange(n):
		f[i,:] = X[i,0]/np.sqrt(2)
		for j in xrange(1,p):
			if j%2 == 1:
				f[i,:] += X[i,j]*np.sin(int((j+1)/2)*t)
			else:
				f[i,:] += X[i,j]*np.cos(int((j+1)/2)*t)
	fig = plt.figure(dpi=300)
	ax = fig.add_subplot(111)
	uniq_colors = list(set(colors))
	idx = [colors.index(color) for color in uniq_colors]
	labels = [classlabels[i] for i in idx]
	for i in idx:
		ax.plot(t,f[i,:],colors[i])
	ax.legend(labels,loc=0)
	for i in xrange(n):
		ax.plot(t,f[i,:],colors[i])
	ax.set_xlabel("$t(-\pi,\ \pi)$",fontstyle='italic')
	ax.set_ylabel("$f(t)$",fontstyle='italic')
	plt.savefig(figname_prefix+".png",format='png',dpi=300)
	plt.savefig(figname_prefix+".svg",format='svg',dpi=300)
	plt.clf()
	plt.close()
	return 0

def plot_linear_regress(X,Y,xlabel,ylabel,classnum,h_uniq_colors,h_uniq_classlabels,figname_prefix="out"):
	##h_uniq_classlabels = {0:'class1',1:'class2'} , 0 and 1 must be the classnum
	##h_uniq_colors = {0:'r^',1:'b.'}
	plt.style.use('grayscale')
	if X.size != Y.size != len(classnum):
		sys.stderr("Error: X, Y should be same dimensions")
		return 1
	slope,intercept,rvalue,pvalue,stderr = stats.linregress(X,Y)
	tmpX = np.linspace(np.min(X),np.max(X),num=50)
	tmpY = tmpX*slope+intercept
	uniq_classnum = list(set(classnum))
	np_classnum = np.array(classnum)
	fig = plt.figure()
	ax = fig.add_subplot(111)
	ax.plot(tmpX,tmpY,'k--')
	ax.grid(True,color='k',alpha=0.5,ls=':')
	for i in uniq_classnum:
		try:
			color = h_uniq_colors[i]
			label = h_uniq_classlabels[i]
		except:
			plt.clf()
			plt.close()
			sys.stderr("Error: key error")
			return 1
		idx = np.where(np_classnum == i)
		ax.plot(X[idx],Y[idx],color,label=label,alpha=0.6)
	ax.legend(loc=0,numpoints=1)
	ax.set_xlabel(xlabel)
	ax.set_ylabel(ylabel)
	ax.set_title('slope:%.3g,intercept:%.3g,r:%.3g,p:%.3g,stderr:%.3g'%(slope,intercept,rvalue,pvalue,stderr))
	plt.savefig(figname_prefix+".png",format='png',dpi=300)
	plt.savefig(figname_prefix+".svg",format='svg',dpi=300)
	plt.clf()
	plt.close()
	return 0

def plot_boxplot(Xnp,fig_prefix,xlabel,ylabel,xticks_labels,outshow=1,colors=None,ylim=1):
	fig = plt.figure(figsize=(10,8),dpi=300)
	ax1 = fig.add_subplot(111)
	if outshow == 1:
		bp = ax1.boxplot(Xnp.T)
		plt.setp(bp['boxes'], color='white')
		plt.setp(bp['whiskers'], color='black')
		plt.setp(bp['fliers'], color='red', marker='+')
	else:
		bp = ax1.boxplot(Xnp.T,0,'')
	if colors:
		n,p = Xnp.shape
		for i in xrange(n):
			box = bp['boxes'][i]
			boxX = box.get_xdata().tolist()
			boxY = box.get_ydata().tolist()
			boxCoords = zip(boxX,boxY)
			boxPolygon = Polygon(boxCoords, facecolor=colors[i])
			ax1.add_patch(boxPolygon)
	ax1.set_xticklabels(xticks_labels,rotation=45)
	ax1.set_xlabel(xlabel)
	ax1.set_ylabel(ylabel)
	if ylim:
		ax1.set_ylim(-10,10)
	plt.savefig(fig_prefix+".png",format='png',dpi=300)
	plt.savefig(fig_prefix+".svg",format='svg',dpi=300)
	plt.clf()
	plt.close()
	return 0

def plot_Xscore(Xnp,classnums,uniqclassnum,uniqcolor,uniqmarker,uniqclasslabel,fig_prefix,xlabel,ylabel,zlabel=None,dim=2):
	plt.style.use('grayscale')
	leng = len(uniqclassnum)
	Xnp = np.asarray(Xnp)
	fig = plt.figure(figsize=(10,8),dpi=300)
	if dim == 3:
		ax1 = fig.add_subplot(111,projection ='3d')
	elif dim==2:
		ax1 = fig.add_subplot(111)
	else:
		sys.stderr.write("[ERROR] Dim '%d' plot failed\n"%dim)
		return 1
	for i in xrange(leng):
		tmpclassidx = np.array(classnums) == uniqclassnum[i]
		tmplabel = uniqclasslabel[i]
		tmpcolor = uniqcolor[i%(len(uniqcolor))]
		tmpmarker = uniqmarker[i%(len(uniqmarker))]
		if dim == 2:
			ax1.plot(Xnp[tmpclassidx,0],Xnp[tmpclassidx,1],ls='',markerfacecolor=tmpcolor,marker=tmpmarker,label=tmplabel,markeredgecolor = tmpcolor,alpha=0.7,markersize=10)
		else:ax1.plot(Xnp[tmpclassidx,0],Xnp[tmpclassidx,1],Xnp[tmpclassidx,2],ls='',markerfacecolor=tmpcolor,marker=tmpmarker,label=tmplabel,markeredgecolor = tmpcolor,alpha=0.7,markersize=10)
	ax1.legend(loc=0,numpoints=1)
	ax1.grid()
	ax1.set_xlabel(xlabel)
	ax1.set_ylabel(ylabel)
	if dim == 3 and zlabel !=None:
		ax1.set_zlabel(zlabel)
	plt.savefig(fig_prefix+".png",format='png',dpi=300)
	plt.savefig(fig_prefix+".svg",format='svg',dpi=300)
	plt.clf()
	plt.close()
	return 0
def plot_XYscore(Xnp,Y,classnums,uniqclassnum,uniqcolor,uniqmarker,uniqclasslabel,fig_prefix,xlabel,ylabel,zlabel=None,dim=2):
	Xnp[:,dim-1] = Y[:,0]
	return plot_Xscore(Xnp,classnums,uniqclassnum,uniqcolor,uniqmarker,uniqclasslabel,fig_prefix,xlabel,ylabel,zlabel,dim)

def plot_markxy(X1,Y1,X2,Y2,xlabel,ylabel,fig_prefix):
	fig = plt.figure(dpi=300)
	ax = fig.add_subplot(111)
	ax.plot(X1,Y1,'b+')
	ax.plot(X2,Y2,'ro')
	ax.set_xlabel(xlabel)
	ax.set_ylabel(ylabel)
	plt.savefig(fig_prefix+".png",format='png',dpi=300)
	plt.savefig(fig_prefix+".svg",format='svg',dpi=300)
	plt.clf()
	plt.close()
	return 0

def plotline(Xvector,Ys,fig_prefix,xlabel,ylabel,colors,legends=None,title=None,xlimmax = None,ylimmax = None):
	n,p = Ys.shape
	fig = plt.figure(dpi=300)
	ax = fig.add_subplot(111)
	if legends:
		leng = len(legends)
	else:
		leng = 0
	for i in xrange(n):
		if i < leng:
			tmplabel = legends[i]
			ax.plot(Xvector,Ys[i,:],colors[i],label=tmplabel)
		else:
			ax.plot(Xvector,Ys[i,:],colors[i])
	if legends != None:
		ax.legend(loc=0)
	#ax.grid()
	ax.set_xlabel(xlabel)
	ax.set_ylabel(ylabel)
	if title != None:
		ax.set_title(title)
	if ylimmax:
		ax.set_ylim(0,ylimmax)
	if xlimmax:
		ax.set_xlim(0,p)
	plt.savefig(fig_prefix+".png",format='png',dpi=300)
	plt.savefig(fig_prefix+".svg",format='svg',dpi=300)
	plt.clf()
	plt.close()
	return 0
def barh_dict_class(hdata,fig_prefix,xlabel,ylabel,title = "",width=0.4,legends=[],colors=[],fmt="%.2f",ylog=0,rotation=0,plot_txt = 1):
	data = []
	yticklabels = []
	classnames = []
	classnumbers = [0] * len(hdata.keys())
	if not colors:
		color_class = cm.Set3(np.linspace(0, 1, len(hdata.keys())))
	else:
		color_class = colors
	idx = 0
	plot_idx = []
	plot_start = 0
	for classname in sorted(hdata.keys()):
		classnames.append(classname)
		for key in hdata[classname]:
			if hdata[classname][key] <=0:continue
			yticklabels.append(key)
			classnumbers[idx] += 1
			data.append(hdata[classname][key])
		plot_idx.append([plot_start,len(data)])
		plot_start += len(data)-plot_start
		idx += 1
	if len(data) > 16:
		fig = plt.figure(figsize=(5,15),dpi=300)
		fontsize_off = 2
	else:
		fig = plt.figure(figsize=(5,7),dpi=300)
	ax = fig.add_subplot(111)
	linewidth = 0
	alpha=0.8
	ylocations = np.asarray(range(len(data)))+width*2
	rects = []
	for i in xrange(len(plot_idx)):
		s,e = plot_idx[i]
		rect = ax.barh(ylocations[s:e],np.asarray(data[s:e]),width,color=color_class[i],linewidth=linewidth,alpha=alpha,align='center')
		rects.append(rect)
	ax.set_yticks(ylocations)
	ax.set_ylabel(ylabel)
	ax.set_xlabel(xlabel)
	ylabelsL = ax.set_yticklabels(yticklabels)
	ax.set_ylim(0,ylocations[-1]+width*2)
	tickL = ax.yaxis.get_ticklabels()
	for t in tickL:
		t.set_fontsize(t.get_fontsize() - 2)
	ax.xaxis.grid()
	ax.legend(classnames,loc=0,fontsize=8)
	#print fig.get_size_inches()
	fig.set_size_inches(10,12)
	fig.tight_layout()
	plt.savefig(fig_prefix+".png",format='png',dpi=300)
	plt.savefig(fig_prefix+".svg",format='svg',dpi=300)
	plt.clf();plt.close();
	return 0
def bar_dict_class(hdata,fig_prefix,xlabel,ylabel,title = "",width=0.4,legends=[],colors=[],fmt="%.2f",ylog=0,rotation=0,plot_txt = 1):
	data = []
	xticklabels = []
	classnames = []
	classnumbers = [0] * len(hdata.keys())
	if not colors:
		color_class = cm.Set3(np.linspace(0, 1, len(hdata.keys())))
	else:
		color_class = colors
	idx = 0
	plot_idx = []
	plot_start = 0
	for classname in sorted(hdata.keys()):
		classnames.append(classname)
		for key in hdata[classname]:
			if hdata[classname][key] <=0:continue
			xticklabels.append(key)
			classnumbers[idx] += 1
			data.append(hdata[classname][key])
		plot_idx.append([plot_start,len(data)])
		plot_start += len(data)-plot_start
		idx += 1
	fontsize_off = 0
	if len(data) > 16:
		fig = plt.figure(figsize=(10,5),dpi=300)
		fontsize_off = 2
	else:
		fig = plt.figure(figsize=(7,5),dpi=300)
	ax = fig.add_subplot(111)
	if ylog:
		ax.set_yscale("log",nonposy='clip')
	linewidth = 0
	alpha=0.8
	xlocations = np.asarray(range(len(data)))+width*2
	#rects = ax.bar(xlocations,np.asarray(data),width,color=plot_colors,linewidth=linewidth,alpha=alpha,align='center')
	rects = []
	for i in xrange(len(plot_idx)):
		s,e = plot_idx[i]
		rect = ax.bar(xlocations[s:e],np.asarray(data[s:e]),width,color=color_class[i],linewidth=linewidth,alpha=alpha,align='center')
		rects.append(rect)
	max_height = 0
	if plot_txt:
		for rk in rects:
			for rect in rk:
				height = rect.get_height()
				if height < 0.1:continue
				ax.text(rect.get_x()+rect.get_width()/2., 1.01*height, fmt%float(height),ha='center', va='bottom',fontsize=(8-fontsize_off))
	ax.set_xticks(xlocations)
	ax.set_ylabel(ylabel)
	ax.set_xlabel(xlabel)
	if rotation == 0 or rotation == 90:hafmt="center"
	else:hafmt="right"
	xlabelsL = ax.set_xticklabels(xticklabels,ha=hafmt,rotation=rotation)
	#print xlocations
	ax.set_xlim(0,xlocations[-1]+width*2)
	tickL = ax.xaxis.get_ticklabels()
	for t in tickL:
		t.set_fontsize(t.get_fontsize() - 2)
	ax.yaxis.grid()
	if ylog:
		ax.set_ylim(0.99,np.max(data)*2)
	else:
		ax.set_ylim(0,np.max(data)*1.35)
	ax.legend(classnames,loc=0,fontsize=8-fontsize_off)
	#ax.xaxis.set_major_locator(plt.NullLocator())
	plt.tick_params(axis='x',          # changes apply to the x-axis
				    which='both',      # both major and minor ticks are affected
					bottom='off',      # ticks along the bottom edge are off
					top='off',         # ticks along the top edge are off
					labelbottom='on') # labels along the bottom edge are off
	fig.tight_layout()
	plt.savefig(fig_prefix+".png",format='png',dpi=300)
	plt.savefig(fig_prefix+".svg",format='svg',dpi=300)
	plt.clf();plt.close();
	return 0

def bar_dict(hdata,fig_prefix,xlabel,ylabel,title = "",width=0.4,legends=[],colors=[],fmt="%.2f",ylog=0,hlist=None,rotation=0):
	data = []
	xticklabels = []
	if hlist == None:
		for key in hdata:
			if hdata[key] <=0:
				continue
			xticklabels.append(key)
			data.append(hdata[key])
	else:
		for key in hlist:
			if hdata[key] <=0:
				continue
			xticklabels.append(key)
			data.append(hdata[key])
	fig = plt.figure(figsize=(7,5),dpi=300)
	ax = fig.add_subplot(111)
	if ylog:
		ax.set_yscale("log",nonposy='clip')
	linewidth = 0
	alpha=0.5
	if not colors:
		colors = cm.Set3(np.linspace(0, 1, len(data)))
	xlocations = np.asarray(range(len(data)))+width*2
	rects = ax.bar(xlocations,np.asarray(data),width,color=colors,linewidth=linewidth,alpha=alpha,align='center')
	for rect in rects:
		height = rect.get_height()
		if height < 0.1:continue
		ax.text(rect.get_x()+rect.get_width()/2., 1.01*height, fmt%float(height),ha='center', va='bottom',fontsize=8)
	ax.set_xticks(xlocations)
	ax.set_ylabel(ylabel)
	ax.set_xlabel(xlabel)
	if rotation == 0 or rotation == 90:
		hafmt='center'
	else:
		hafmt = 'right'
	xlabelsL = ax.set_xticklabels(xticklabels,ha=hafmt,rotation=rotation)
	#if rotation:
	#	for label in xlabelsL:
	#		label.set_rotation(rotation)
	ax.set_title(title)
	ax.set_xlim(0,xlocations[-1]+width*2)
	tickL = ax.xaxis.get_ticklabels()
	for t in tickL:
		t.set_fontsize(t.get_fontsize() - 2)
	ax.yaxis.grid()
	#ax.set_adjustable("datalim")
	if ylog:
		ax.set_ylim(0.99,np.max(data)*2)
	else:
		ax.set_ylim(0,np.max(data)*1.5)
	#ax.set_ylim(ymin=0)
	fig.tight_layout()
	plt.savefig(fig_prefix+".png",format='png',dpi=300)
	plt.savefig(fig_prefix+".svg",format='svg',dpi=300)
	plt.clf();plt.close();
	return 0

def stackv_bar_plot(data,xticks_labels,fig_prefix,xlabel,ylabel,title="",width=0.3,legends=[],colors=[],scale=0):
	n,p = data.shape
	ind = np.arange(p)
	tmp = np.float64(data.copy())
	#tmp = np.cumsum(data,0)
	#print tmp - data
	if scale:
		tmp = tmp/np.sum(tmp,0)*100
		print tmp
	#tmp = np.cumsum(tmp,0)
	if not colors:
		colors = cm.Set3(np.linspace(0, 1, n))
	fig = plt.figure(figsize=(14,8),dpi=300)
	ax = fig.add_subplot(111)
	linewidth = 0
	alpha=0.8
	for i in xrange(n):
		if i:
			cumtmp = cumtmp + np.asarray(tmp[i-1,:])[0]
			rects = ax.bar(ind,np.asarray(tmp[i,:])[0],width,color=colors[i],linewidth=linewidth,alpha=alpha,bottom=cumtmp,align='center',label=legends[i])
		else:
			cumtmp = 0
			print ind,np.asarray(tmp[i,:])[0]
			rects = ax.bar(ind,np.asarray(tmp[i,:])[0],width,color=colors[i],linewidth=linewidth,alpha=alpha,align='center',label=legends[i])
	ax.legend(loc=0)
	ax.set_ylabel(ylabel)
	ax.set_xlabel(xlabel)
	ax.set_xticks(ind)
	ax.set_xticklabels(xticks_labels)
	ax.set_title(title)
	if scale:
		ax.set_ylim(0,100)
	plt.savefig(fig_prefix+".png",format='png',dpi=300)
	plt.savefig(fig_prefix+".svg",format='svg',dpi=300)
	plt.clf()
	plt.close()
	return 0

def err_bar_group(data,error,group_label,xticklabel,xlabel,ylabel,colors,fig_prefix,title=None,width=0.3):
	num_groups,p = data.shape
	assert num_groups == len(group_label)
	fig = plt.figure(dpi=300)
	ax = fig.add_subplot(111)
	xlocations = np.arange(p)
	for i in xrange(num_groups):
		ax.bar(xlocations+width*i, data[i,:],yerr=error[i,:], width=width,linewidth=0,color=colors[i],ecolor=colors[i],capsize=5,alpha=0.6,label=group_label[i])
	ax.legend(group_label,loc=0)
	ax.set_xticks(xlocations+width/2*num_groups)
	ax.set_xticklabels(xticklabel)
	ax.set_ylabel(ylabel)
	ax.set_xlabel(xlabel)
	if title <> None:ax.set_title(title)
	fig.tight_layout()
	plt.savefig(fig_prefix+".png",format='png',dpi=300)
	plt.savefig(fig_prefix+".svg",format='svg',dpi=300)
	plt.clf()
	plt.close()
	return 0
	
def err_bar(data,error,xlabel,ylabel,fig_prefix,title=None,mark_sig=None,mark_range=[[0,1],],width=0.3):
	num = len(data)
	assert num == len(error) == len(xlabel)
	#colors = cm.Set3(np.linspace(0, 1, len(xlabel)))
	colors = ["black","gray"]
	fig = plt.figure(dpi=300)
	ax = fig.add_subplot(111)
	xlocations = np.asarray(range(len(data)))+width
	ax.bar(xlocations, data, yerr=error, width=width,linewidth=0,ecolor='black',capsize=5,color=colors,alpha=0.5)
	ax.set_xticks(xlocations+width/2)
	ax.set_xticklabels(xlabel)
	ax.set_ylabel(ylabel)
	ax.set_xlim(0, xlocations[-1]+width*2)
	if title <> None:ax.set_title(title)
	if mark_sig <> None:
		xlocations = xlocations+width/2
		ybin = np.max(np.asarray(data)+np.asarray(error))
		step = ybin/20
		offset = ybin/40
		assert len(mark_sig) == len(mark_range)
		for i in xrange(len(mark_range)):
			mark_r = mark_range[i]
			sig_string = mark_sig[i]
			xbin = np.asarray(mark_r)
			ybin += step
			ax.plot([xlocations[mark_r[0]],xlocations[mark_r[1]]],[ybin,ybin],'k-')
			ax.text((xlocations[mark_r[0]]+xlocations[mark_r[1]])/2,ybin+offset,sig_string)
		ax.set_ylim(0,ybin+step*2.5)
	fig.tight_layout()
	plt.savefig(fig_prefix+".png",format='png',dpi=300)
	plt.savefig(fig_prefix+".svg",format='svg',dpi=300)
	plt.clf()
	plt.close()
	return 0

def trendsiglabel(Xvec,Yvec,meansdata,totmean,color,xticklabel,fig_prefix="trend"):
	num = len(Xvec)
	ngenes_sig,p = meansdata.shape
	ngenes_tot,p = totmean.shape
	assert num == len(Yvec) == len(xticklabel) == p
	
	fig = plt.figure(dpi=300)
	ax = fig.add_subplot(111)
	#ax.plot(Xvec,Yvec,color+'^-',markeredgecolor='None',markersize = 12)
	for i in xrange(ngenes_tot):
		#print i
		ax.plot(Xvec,totmean[i,:],'g-',lw=0.5,alpha=0.2)
	for i in xrange(ngenes_sig):
		ax.plot(Xvec,meansdata[i,:],'b-',lw=0.8,alpha=0.6)
	ax.plot(Xvec,Yvec,color+'^-',markeredgecolor=color,markersize = 5)
	ax.set_xticks(np.arange(num))
	xlabelsL = ax.set_xticklabels(xticklabel)
	#clean y
	#ax.get_yaxis().set_ticks([])
	#min_a = np.min(Yvec)
	#max_a = np.max(Yvec)
	#ax.set_ylim(min_a-1,max_a+1)
	fig.tight_layout()
	plt.savefig(fig_prefix+".png",format='png',dpi=300)
	plt.savefig(fig_prefix+".svg",format='svg',dpi=300)
	plt.clf()
	plt.close()
	return 0

def twofactor_diff_plot(Xmeanarr,Xstdarr,xticklabel,fig_prefix="Sigplot",title=None,xlabel=None,ylabel="Expression",width=0.3,labels=None):
	num = Xmeanarr.shape[-1]
	fmts = ['o-','^--','x-.','s--','v-.','+-.']
	ecolors = ['r','b','g','c','m','y','k']
	assert num == Xstdarr.shape[-1] == len(xticklabel)
	fig = plt.figure(dpi=300)
	ax = fig.add_subplot(111)
	xlocations = np.asarray(range(num))+width
	n,p = Xmeanarr.shape
	for i in xrange(n):
		ax.errorbar(xlocations, Xmeanarr[i,:], yerr=Xstdarr[i,:],fmt=fmts[i],ecolor='r',markeredgecolor='None')
	if labels:
		ax.legend(labels,loc=0,numpoints=1)
	ax.set_xticks(xlocations)
	ax.set_xticklabels(xticklabel)
	ax.set_ylabel(ylabel)
	if xlabel: ax.set_xlabel(xlabel)
	ax.set_xlim(0, xlocations[-1]+width*2)
	if title <> None:ax.set_title(title)
	fig.tight_layout()
	plt.savefig(fig_prefix+".png",format='png',dpi=300)
	plt.savefig(fig_prefix+".svg",format='svg',dpi=300)
	plt.clf()
	plt.close()
	return 0

def onefactor_diff_plot(Xmeanarr,Xstdarr,xticklabel,fig_prefix="Sigplot",title=None,xlabel=None,ylabel="Expression",width=0.3):
	num = len(Xmeanarr)
	assert num == len(Xstdarr) == len(xticklabel)
	fig = plt.figure(dpi=300)
	ax = fig.add_subplot(111)
	xlocations = np.asarray(range(len(Xmeanarr)))+width
	ax.errorbar(xlocations, Xmeanarr, yerr=Xstdarr,fmt='o-',ecolor='r')
	ax.set_xticks(xlocations)
	ax.set_xticklabels(xticklabel)
	ax.set_ylabel(ylabel)
	if xlabel: ax.set_xlabel(xlabel)
	ax.set_xlim(0, xlocations[-1]+width*2)
	if title <> None:ax.set_title(title)
	fig.tight_layout()
	plt.savefig(fig_prefix+".png",format='png',dpi=300)
	plt.savefig(fig_prefix+".svg",format='svg',dpi=300)
	plt.clf()
	plt.close()
	return 0

def bar_plot(data,xticks_labels,fig_prefix,xlabel,ylabel,title="",width=0.3,rotation=0,fmt='%.2f',ylog=0,colors=None):
	ind = np.arange(len(data))
	fig = plt.figure()
	ax = fig.add_subplot(111)
	if ylog:
		ax.set_yscale("log",nonposy='clip')
	linewidth = 0
	alpha=0.5
	if not colors:
		colors = 'k'
	rects = ax.bar(ind,data,width,color=colors,linewidth=linewidth,alpha=alpha,align='center')
	ax.set_ylabel(ylabel)
	ax.set_xlabel(xlabel)
	ax.set_xticks(ind)
	ax.yaxis.grid()
	#ax.set_xticks(ind+width/2)
	if rotation == 0 or rotation == 90:hafmt='center'
	else:hafmt = 'right'
	xlabelsL = ax.set_xticklabels(xticks_labels,ha=hafmt,rotation=rotation)
	#rotate labels 90 degrees
	if rotation:
		for label in xlabelsL:
			label.set_rotation(rotation)
	ax.set_title(title)
	for rect in rects:
		height = rect.get_height()
		if height < 0.1:continue
		ax.text(rect.get_x()+rect.get_width()/2., 1.01*height, fmt%float(height),ha='center', va='bottom',fontsize=8)
	fig.tight_layout()
	plt.savefig(fig_prefix+".png",format='png',dpi=300)
	plt.savefig(fig_prefix+".svg",format='svg',dpi=300)
	plt.clf()
	plt.close()
	return 0

def MA_vaco_plot(avelogFC,logFC,totavelogFC,totlogFC,fig_prefix,xlabel,ylabel,title="MAplot"):
	fig = plt.figure(figsize=(5,4),dpi=300)
	ax = fig.add_subplot(111)
	ax.plot(avelogFC,logFC,'ro',markersize = 5,alpha=0.5,markeredgecolor='r')
	ax.plot(totavelogFC,totlogFC,'bo',markersize = 3,alpha=0.5,markeredgecolor='b')
	ax.set_xlabel(xlabel)
	ax.set_ylabel(ylabel)
	ax.set_title(title)
	fig.tight_layout()
	plt.savefig(fig_prefix+".png",format='png',dpi=300)
	plt.savefig(fig_prefix+".svg",format='svg',dpi=300)
	plt.clf()
	plt.close()
	return 0

def vaco_plot(X,Y,Xcut,Ycut,fig_prefix,xlabel,ylabel,title=None):
	# X is rho or fc
	Xcutx = [np.min(X),np.max(X)]
	Ycuts = [Ycut,Ycut]
	idx1 = (Y > Ycut) & (np.abs(X) > Xcut)
	idx2 = ~idx1
	fig = plt.figure(figsize=(5,4),dpi=300)
	ax = fig.add_subplot(111)
	ax.plot(X[idx1],Y[idx1],'ro',markersize = 5,alpha=0.5,markeredgecolor='None')
	ax.plot(X[idx2],Y[idx2],'bo',markersize = 5,alpha=0.5,markeredgecolor='None')
	ax.plot(Xcutx,Ycuts,'r--')
	ax.set_xlabel(xlabel)
	ax.set_ylabel(ylabel)
	#ax.set_xlim(-6,6)
	if title != None:
		ax.set_title(title)
	plt.savefig(fig_prefix+".png",format='png',dpi=300)
	plt.savefig(fig_prefix+".svg",format='svg',dpi=300)
	plt.clf()
	plt.close()
	return 0

def baohedu_plot(genes,reads,samples,fig_prefix,xlabel="number of reads",ylabel="number of detected genes",title=None):
	n1,p1 = genes.shape
	n2,p2 = reads.shape
	assert n1==n2 and p1==p2
	"saturability"
	#types = ['ro-','b^--','gs-.','kv:','c^-.','m*--','yp:']
	ret_color,ret_lines,ret_marker = styles(n1)
	fig = plt.figure(figsize=(8,6),dpi=300)
	ax = fig.add_subplot(111)
	for i in xrange(n1):
		x = reads[i,:]
		y = genes[i,:]
		ax.plot(x,y,color=ret_color[i],linestyle=ret_lines[i],marker=ret_marker[i],markeredgecolor=ret_color[i],markersize = 4,alpha=0.7,label=samples[i])
	ax.set_xlabel(xlabel)
	ax.set_ylabel(ylabel)
	if title != None: ax.set_title(title)
	ax.legend(loc=0,numpoints=1)
	plt.savefig(fig_prefix+".png",format='png',dpi=300)
	plt.savefig(fig_prefix+".svg",format='svg',dpi=300)
	plt.tight_layout()
	plt.clf()
	plt.close()
	return 0

def plotyy(Xvector,Y1np,Y2np,fig_prefix,xlabel,ylabel1,ylabel2,title=None):
	Y1np = np.asarray(Y1np)
	Y2np = np.asarray(Y2np)
	fig = plt.figure(figsize=(10,8),dpi=300)
	ax1 = fig.add_subplot(111)
	try:
		n1,p1 = Y1np.shape
	except ValueError:
		n1 = 1
	try:
		n2,p2 = Y2np.shape
	except ValueError:
		n2 = 1
	for i in xrange(n1):
		if n1 == 1:
			ax1.plot(Xvector,Y1np, 'b-')
			break
		if i == 0:
			ax1.plot(Xvector,Y1np[i,:], 'b-')
		else:
			ax1.plot(Xvector,Y1np[i,:], 'b--')
	ax1.set_xlabel(xlabel)
	ax1.set_ylabel(ylabel1, color='b')
	if title: ax1.set_title(title)
	for tl in ax1.get_yticklabels():
		tl.set_color('b')
	ax2 = ax1.twinx()
	for i in xrange(n2):
		if n2 == 1:
			ax2.plot(Xvector,Y2np, 'r-')
			break
		if i == 0:
			ax2.plot(Xvector,Y2np[i,:], 'r-')
		else:
			ax2.plot(Xvector,Y2np[i,:], 'r-.')
	ax2.set_ylabel(ylabel2, color='r')
	for tl in ax2.get_yticklabels():
		tl.set_color('r')
	ax1.grid()
	plt.savefig(fig_prefix+".png",format='png',dpi=300)
	plt.savefig(fig_prefix+".svg",format='svg',dpi=300)
	plt.clf()
	plt.close()
	return 0

def clean_axis(ax):
	"""Remove ticks, tick labels, and frame from axis"""
	ax.get_xaxis().set_ticks([])
	ax.get_yaxis().set_ticks([])
	for spx in ax.spines.values():
		spx.set_visible(False)
def density_plt(Xarr,colors,legendlabel,figname_prefix="density",xlabel=None,ylabel=None,fun="pdf"):
	"""not at the same scale  """
	fig = plt.figure(dpi=300)
	ax = fig.add_subplot(111)
	n = len(Xarr)
	assert len(colors) == len(legendlabel)
	for i in xrange(n):
		dat = Xarr[i]
		xp,yp = kdensity(np.asarray(dat),num = 200,fun=fun)
		ax.plot(xp,yp,colors[i],label=legendlabel[i],markeredgecolor='None')
	ax.legend(loc=0,numpoints=1)
	if xlabel: ax.set_xlabel(xlabel)
	if ylabel: ax.set_ylabel(ylabel)
	plt.savefig(figname_prefix+".png",format='png',dpi=300)
	plt.savefig(figname_prefix+".svg",format='svg',dpi=300)
	plt.clf()
	plt.close()
	return 0
def exprs_density(Xnp,colors,classlabels,figname_prefix="out",xlabel=None,ylabel=None,fun="cdf"):
	n,p = Xnp.shape
	fig = plt.figure(dpi=300)
	ax = fig.add_subplot(111)
	uniq_colors = list(set(colors))
	idx = [colors.index(color) for color in uniq_colors]
	labels = [classlabels[i] for i in idx]
	for i in idx:
		if fun == "cdf":
			xp,yp = kdensity(np.asarray(Xnp[i,:]),fun="cdf")
		elif fun == "pdf":
			xp,yp = kdensity(np.asarray(Xnp[i,:]),fun="pdf")	
		ax.plot(xp,yp,colors[i])
	ax.legend(labels,loc=0)
	for i in xrange(n):
		if fun == "cdf":
			xp,yp = kdensity(np.asarray(Xnp[i,:]),fun="cdf")
		elif fun == "pdf":
			xp,yp = kdensity(np.asarray(Xnp[i,:]),fun="pdf")
		ax.plot(xp,yp,colors[i])
	ax.set_xlabel(xlabel)
	ax.set_ylabel(ylabel)
	ax.set_xlim(0,10)
	fig.tight_layout()
	plt.savefig(figname_prefix+".png",format='png',dpi=300)
	plt.savefig(figname_prefix+".svg",format='svg',dpi=300)
	plt.clf()
	plt.close()
	return 0

def hist_groups(data,labels,xlabel,fig_prefix,bins=25,alpha=0.7,normed=True,colors=None):
	n = len(data)
	assert n == len(labels)
	if not colors:
		if n >7:
			colors = cm.Set1(np.linspace(0, 1,n))
		else:
			colors = plt.rcParams['axes.color_cycle'][0:n]
	if normed:ylabel = "Density"
	else:ylabel = "Frequency"
	fig = plt.figure(dpi=300)
	ax = fig.add_subplot(111)
	for i in xrange(n):
		xp,yp = kdensity(data[i],fun="pdf")
		ax.hist(data[i],histtype="stepfilled",bins=bins, alpha=alpha,normed=normed,label=labels[i],color=colors[n-i-1])
		ax.plot(xp,yp,color=colors[n-i-1],linestyle='--',lw=1.0)
	ax.set_xlabel(xlabel)
	ax.set_ylabel(ylabel)
	ax.legend(loc=0)
	plt.savefig(fig_prefix+".png",format='png',dpi=300)
	plt.savefig(fig_prefix+".svg",format='svg',dpi=300)
	plt.clf()
	plt.close()
	return 0

def exprs_RLE(Xnp,mean="median",fig_prefix=None,samplenames=None,colors=None):
	###!!!!用median 保持robust
	#在同一组实验中，即使是相互比较的对照组与实验组之间，大部分基因的表达量还是应该保持一致的，何况平行实验之间。当我们使用相对对数表达（Relative Log Expression(RLE)）的的箱线图来控制不同组之间的实验质量时，我们会期待箱线图应该在垂直中央相类的位置（通常非常接近0）。如果有一个芯片的表现和其它的平行组都很不同，那说明它可能出现了质量问题。
	n,p = Xnp.shape
	if mean == "median":
		Xmean =np.median(Xnp,axis=0)	
	elif mean == "mean":
		Xmean =np.mean(Xnp,axis=0)
	plot_boxplot(Xnp-Xmean,fig_prefix,"","Relative Log Expression",samplenames,colors=colors,ylim=0)
	
	return 0

def exprs_NUSE():
	#1. hist 
	#2. julei 
	#3. RLE
	#array corr
	#
	#相对标准差（Normalized Unscaled Standard Errors(NUSE)）
	#是一种比RLE更为敏感 的质量检测手段。如果你在RLE图当中对某组芯片的质量表示怀疑，那当你使用NUSE图时，这种怀疑就很容易被确定下来。NUSE的计算其实也很简单，它是某芯片基因标准差相对于全组标准差的比值。我们期待全组芯片都是质量可靠的话，那么，它们的标准差会十分接近，于是它们的NUSE值就会都在1左右。然而，如果有实验芯片质量有问题的话，它就会严重的偏离1，进而影响其它芯片的NUSE值偏向相反的方向。当然，还有一种非常极端的情况，那就是大部分芯片都有质量问题，但是它们的标准差却比较接近，反而会显得没有质量问题的芯片的NUSE值会明显偏离1，所以我们必须结合RLE及NUSE两个图来得出正确的结论
	return 0

def exprs_corrarray(Xnp,samplenames,figname_prefix):
	corr_coef = np.abs(np.corrcoef(Xnp))
	n,p = corr_coef.shape
	fig = plt.figure(figsize=(8,6),dpi=300)
	ax = fig.add_subplot(111)
	clean_axis(ax)
	image_instance = ax.imshow(corr_coef,interpolation='nearest',aspect='auto',cmap=cm.autumn,alpha=0.8,origin='lower')
	ax.set_yticks(np.arange(p))
	ax.set_yticklabels(samplenames)
	ax.set_xticks(np.arange(n))
	xlabelsL = ax.set_xticklabels(samplenames)
	for label in xlabelsL:
		label.set_rotation(90)
	cb = fig.colorbar(image_instance,ax=ax)
	cb.set_label("correlation")
	cb.outline.set_linewidth(0)
	ax.grid(visible=False)
	fig.tight_layout()
	plt.savefig(figname_prefix+".png",format='png',dpi=300)
	plt.savefig(figname_prefix+".svg",format='svg',dpi=300)
	plt.clf()
	plt.close()
	return corr_coef

def pie_plot(sizes,labels,fig_prefix="pie_plot",autopct='%1.1f%%',colors=None,explode=None,shadow=False, startangle=90,radius=1):
	fig = plt.figure(figsize=(6,6),dpi=300)
	ax5 = fig.add_subplot(111)
	if not colors:
		colors = cm.Set3(np.linspace(0, 1, len(labels)))
	#patches, texts, autotexts = ax5.pie(sizes,explode,labels=labels, colors=colors,autopct=autopct, shadow=shadow, startangle=startangle,radius=radius)
	patches, texts = ax5.pie(sizes,explode,colors=colors, shadow=shadow, startangle=startangle,radius=radius)
	tmplabels = []
	total = sum(sizes)
	for i in xrange(len(labels)):
		lable = labels[i]
		size = float(sizes[i])/total*100
		#print lable+"("+ autopct+")"
		tmplabels.append((lable+"("+ autopct+")")%size)
	ax5.legend(patches,tmplabels,loc='best')
	for w in patches:
		w.set_linewidth(0.2)
		w.set_edgecolor('white')
	##plt.legend(patches, labels, loc="best")
	#proptease = fm.FontProperties()
	#proptease.set_size('xx-small')
	#plt.setp(autotexts, fontproperties=proptease)
	#plt.setp(texts, fontproperties=proptease)
	plt.axis('equal')
	plt.tight_layout()
	plt.savefig(fig_prefix+".png",format='png',dpi=300)
	plt.savefig(fig_prefix+".svg",format='svg',dpi=300)
	plt.clf()
	plt.close()
	return 0

def cluster_heatmap(Xnp,samplenames,annonames,fig_prefix="test_cluster_heatmap",colornorm = True,nosample=False,nogene=False):
	n,p = Xnp.shape
	assert n == len(samplenames) and p == len(annonames)
	# make norm
	if colornorm:
		vmin = np.floor(np.min(Xnp))
		vmax = np.ceil(np.max(Xnp))
		vmax = max([vmax,abs(vmin)]) # choose larger of vmin and vmax
		#vmin = vmax * -1
		my_norm = mpl.colors.Normalize(vmin, vmax)
	else:my_norm = None
	# heatmap with row names
	if len(samplenames)/3 <=7:
		rightx = 7
	else:
		rightx = len(samplenames)/3
	if len(annonames)/3 <=8:
		leftx = 8
	else:
		leftx = int(len(annonames)/4)
	leftx = min(20,leftx)
	sys.stderr.write("[INFO] plot size is %dX%d\n"%(leftx,rightx))
	fig = plt.figure(figsize=(rightx,leftx))
	heatmapGS = gridspec.GridSpec(2,2,wspace=0.0,hspace=0.0,width_ratios=[0.25,1],height_ratios=[0.15,1])
	### col dendrogram ### col is sample cluster
	if not nosample:
		col_pairwise_dists = sp.spatial.distance.squareform(sp.spatial.distance.pdist(Xnp))
		col_clusters = linkage(col_pairwise_dists,method='ward')
		col_denAX = fig.add_subplot(heatmapGS[0,1])
		col_denD = dendrogram(col_clusters)
		col_denAX.set_axis_off()
	###  fcluster(col_clusters,0.7*max(col_clusters[:,2]),'distance') 
	###  to return the index of each sample for each cluster
	
	### row dendrogram ### row is anno cluster
	if not nogene:
		row_pairwise_dists = sp.spatial.distance.squareform(sp.spatial.distance.pdist(Xnp.T))
		row_clusters = linkage(row_pairwise_dists,method='ward')
		row_denAX = fig.add_subplot(heatmapGS[1,0])
		row_denD = dendrogram(row_clusters,orientation='right')
		row_denAX.set_axis_off()
	### heatmap ####
	heatmapAX = fig.add_subplot(heatmapGS[1,1])
	if nogene:
		Xtmp = Xnp.T.copy()
	else:
		Xtmp = Xnp.T[row_denD['leaves'],:]
	if nosample:
		pass
	else:
		Xtmp = Xtmp[:,col_denD['leaves']]
	axi = heatmapAX.imshow(Xtmp,interpolation='nearest',aspect='auto',origin='lower',norm=my_norm,cmap=cm.coolwarm)
	clean_axis(heatmapAX)
	## row labels ##
	if not nogene:
		t_annonames = [annonames[i] for i in row_denD['leaves']]
	else:
		t_annonames = annonames
	heatmapAX.set_yticks(np.arange(p))
	heatmapAX.yaxis.set_ticks_position('right')
	heatmapAX.set_yticklabels(t_annonames)
	## col labels ##
	if not nosample:
		t_samplenames = [samplenames[i] for i in col_denD['leaves']]
	else:
		t_samplenames = samplenames
	heatmapAX.set_xticks(np.arange(n))
	xlabelsL = heatmapAX.set_xticklabels(t_samplenames)
	#rotate labels 90 degrees
	for label in xlabelsL:
		label.set_rotation(90)
	#remove the tick lines
	for l in heatmapAX.get_xticklines() + heatmapAX.get_yticklines():
		l.set_markersize(0)
	heatmapAX.grid(False)
	### scale colorbar ###
	#scale_cbGSSS = gridspec.GridSpecFromSubplotSpec(1,2,subplot_spec=heatmapGS[0,0],wspace=0.0,hspace=0.0)
	#scale_cbAX = fig.add_subplot(scale_cbGSSS[0,1])
	scale_cbGSSS = gridspec.GridSpecFromSubplotSpec(1,1,subplot_spec=heatmapGS[0,0],wspace=0.0,hspace=0.0)
	scale_cbAX = fig.add_subplot(scale_cbGSSS[0,0])
	scale_cbAX.set_axis_off()
	#cb = fig.colorbar(axi,ax=scale_cbAX,orientation='horizontal')
	cb = fig.colorbar(axi,ax=scale_cbAX)
	cb.set_label('Measurements')
	cb.ax.yaxis.set_ticks_position('left')
	cb.ax.yaxis.set_label_position('left')
	cb.outline.set_linewidth(0)
	tickL = cb.ax.yaxis.get_ticklabels()
	for t in tickL:
		t.set_fontsize(t.get_fontsize() - 3)
	fig.tight_layout()
	plt.savefig(fig_prefix+".png",format='png',dpi=300)
	plt.savefig(fig_prefix+".svg",format='svg',dpi=300)
	plt.clf()
	plt.close()
	return 0

