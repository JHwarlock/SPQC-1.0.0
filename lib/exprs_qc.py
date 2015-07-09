# -*- coding: UTF-8 -*-
import sys
import numpy as np
import scipy as sp
from scipy import stats
import mutilstats
import statplot
import mhtml
import mds
class exprsqc_opt():
	def __init__(self):
		self.phenotype = ""
		self.log2tr = 0
		self.noise = -1
		self.sample_info = ""
		self.matrixdata = ""
		#self.outdir = ""
		self.addbg = 0
		self.normalize = 0
	def destroy(self):
		del self.phenotype,self.log2tr,self.noise,self.sample_info,self.matrixdata,self.addbg
		return 0

def fmtout(a):
	 return "%.3f"%a

def qcrun(opt):
	# init 
	fsample = opt.sample_info
	fmatrix = opt.matrixdata

	sampleinfo = mutilstats.SampleInfo()
	ret = sampleinfo.parse_sampleinfo(fsample)
	if ret <> 0:
		sys.stderr.write("[ERROR] Sample information failed to parse, please check the sampleinfo file\n")
		return 1
	## process
	data = mutilstats.MatrixAnno()
	if opt.log2tr == 1:
		ret = data.parse_matrix_anno(fmatrix,addtolog=opt.addbg,log2tr=1,cutoff=opt.noise)
	else:
		assert opt.log2tr == 0
		ret = data.parse_matrix_anno(fmatrix,cutoff=opt.noise)
	if opt.normalize:
		ret = mutilstats.normalize(data.data)
	if ret <> 0:
		sys.stderr.write("[ERROR] Data parse failed, please check the matrix file\n")
		return 1
	
	#1 CDF
	statplot.exprs_density(data.data,sampleinfo.classcolors,sampleinfo.classlabels,"exprs_CDF","expression level","cumulative distribution","cdf")
	#2 PDF
	statplot.exprs_density(data.data,sampleinfo.classcolors,sampleinfo.classlabels,"exprs_pdf","expression level","probability density distribution","pdf")
	#3 boxplot
	statplot.plot_boxplot(data.data,"exprs_boxplot","","expression level",sampleinfo.samplenames,colors=sampleinfo.classcolors,ylim=0)
	#4 RLE
	statplot.exprs_RLE(data.data,"mean","RLE_plot",sampleinfo.samplenames,colors=sampleinfo.classcolors)

	#5 corr 
	corr_matrix = statplot.exprs_corrarray(data.data,sampleinfo.samplenames,"corrarray")
	## out corr_matrix
	foutcorr = file("corr_matrix.xls","w")
	foutcorr.write("## correlation coefficient matrix for samples' expression level\n")
	foutcorr.write("\t".join(["#correlation"]+sampleinfo.samplenames)+"\n")
	for i in xrange(len(sampleinfo.samplenames)):
		foutcorr.write("\t".join([sampleinfo.samplenames[i],]+map(fmtout,corr_matrix[i,:].tolist()))+"\n")
	foutcorr.close()

	#6 cluster
	statplot.hcluster(data.data,sampleinfo.samplenames,"hcluster")

	#7 MDS

	#statdir = "Exprs"
	html_main = mhtml.simple_main(title="样本表达质量控制结果",css="../CSS")
	html_main.add_head("样本表达质量控制结果")
	html_main.add_enter()
	#html_main.add_back1()
	#html_main.add_enter()
	html_main.add_head("1. 样本信息列表",2)
	html_main.add_line()
	html_main.add_enter()
	tmptable,tmpnote = mhtml.xls2table("%s"%fsample)
	html_main.add_content(tmptable)
	html_main.add_precontent(tmpnote)
	html_main.add_enter()
	html_main.add_head("2. 样本表达质量控制结果,数据可靠性分析",2)
	html_main.add_line()
	html_main.add_enter()
	html_main.add_head("a. 样本表达水平概率密度分布",3)
	html_main.add_enter()
	if opt.log2tr == 1:
		strlog = "采用log<sub>2</sub>变换，并"
	else:
		strlog = ""
	html_main.add_content("""对各样本表达水平，%s计算概率密度。查看各样本及各组间表达水平的分布情况。概率密度估计采用Kernel density estimation， implementation in python with scipy（http://www.scipy.org/）"""%strlog)
	html_main.add_content("""<img src="./exprs_pdf.png" width="50%" /><a href="./exprs_pdf.svg">SVG矢量图版本</a>""")
	html_main.add_enter()
	html_main.add_head("b. 样本表达水平累积概率密度分布",3)
	html_main.add_enter()
	html_main.add_content("""对各样本表达水平，%s计算累积概率密度。查看各样本及各组间表达水平的分布情况。"""%strlog)
	html_main.add_content("""<img src="./exprs_CDF.png" width="50%" /><a href="./exprs_CDF.svg">SVG矢量图版本</a>""")
	html_main.add_enter()
	html_main.add_head("c. 样本表达水平箱式图",3)
	html_main.add_content("""对各样本表达水平，%s绘制箱式图。查看各样本及各组间表达水平的分布情况。"""%strlog)
	html_main.add_content("""<img src="./exprs_boxplot.png" width="50%" /><a href="./exprs_boxplot.svg">SVG矢量图版本</a>""")
	html_main.add_enter()

	html_main.add_head("d. 样本间表达相关性分析",3)
	html_main.add_enter()
	html_main.add_content("""计算样本两两间的fisher 相关系数， 将相关系数矩阵按实验分组形式，绘制成热图。样本处理组间，大部分表达具有相关性，主要是因为维持生命基本活动的大部分基因均不差异表达，只有少部分为差异表达（当处理条件十分剧烈时，实验组和处理组间可能并不满足此假设，但组内样本应满足此假设）。因此，各样本间，表达水平相关性应较高。若图中存在特异性的样本，或实验条件不统一且未做校正时，该特殊样本与其他样本的表达相关性会非常低。""")
	html_main.add_content("""<img src="./corrarray.png" width="50%" /><a href="./corrarray.svg">SVG矢量图版本</a>""")
	tmptable,tmpnote = mhtml.xls2table("corr_matrix.xls")
	html_main.add_content(tmptable)
	html_main.add_enter()
	html_main.add_head("e. 样本相对表达水平比较",3)
	html_main.add_enter()
	html_main.add_content("""在同一组实验中，即使是相互比较的对照组与实验组之间，大部分基因的表达量还是应该保持一致的。当使用相对对数表达水平（Relative Log Expression(RLE)）的箱线图来控制不同组之间的实验质量时，箱线图应该在垂直中央相类的位置（通常接近0）。如果有一个样本的表现和其它的平行组都很不同，那说明它可能出现了质量问题。""")
	html_main.add_content("""<img src="./RLE_plot.png" height="450" width="550" /><a href="./RLE_plot.svg">SVG矢量图版本</a>""")
	html_main.add_enter()
	html_main.add_head("f. 样本聚类结果",3)
	html_main.add_enter()
	html_main.add_content("""基于表达水平数据的样本聚类，计算样本间欧式距离，采用离差平方和法（wald法）进行层次聚类，验证聚类结果是否同实验设计基本一致。若聚类结果明显不一致，则样本间存在着明显的其他未知的因素，而不仅仅是实验处理效应。""")
	html_main.add_content("""<img src="./hcluster_hcluster.png" width="50%" /><a href="./hcluster_hcluster.svg">SVG矢量图版本</a>""")
	html_main.add_enter()

	html_main.add_head("g. 样本多维尺度分析",3)
	html_main.add_content("多维尺度分析(Multi Dimensional Scaling, MDS)是一种将多维空间的研究对象简化到低维空间进行定位、分析和归类，同时又保留对象间原始关系的数据分析方法。此处我们采用样本间欧式距离反映样本间的差异，选择前3个本征值最大的维度，绘制样本在前三个维度上的分布，若实验处理因素为表达差异的主要因素，则一般而言，样本组内差异应小于组间差异。")

	sinfo = sampleinfo
	snnum = len(sinfo.sns)

	mdsout = mds.mds_ps(data.data,10)
	xlabel = "Number of dimensions"
	ylabel = "Variation percentage"
	statplot.plotline(np.arange(0,len(mdsout.p)+1),np.asarray([[0,]+mdsout.p.tolist(),]),"Variation_percentage",xlabel,ylabel,['r^-'],xlimmax=10+1,ylimmax=102)

	if snnum == 2:
		ret = statplot.plot_Xscore(mdsout.v,sinfo.classnums,sinfo.uniqclassnum,sinfo.uniqcolor,sinfo.uniqmarker,sinfo.uniqclasslabel,"MDS_samples_distribution","1st dimension","2nd dimension","3rd dimension",dim=2)
		html_main.add_content("前n个维度累积解释变异的百分比图(见下图)。其中，前两个维度，累积解释变异的百分比: %.2f%%, %.2f%%"%(mdsout.p[0],mdsout.p[1]))
	elif snnum >= 3:
		ret = statplot.plot_Xscore(mdsout.v,sinfo.classnums,sinfo.uniqclassnum,sinfo.uniqcolor,sinfo.uniqmarker,sinfo.uniqclasslabel,"MDS_samples_distribution","1st dimension","2nd dimension","3rd dimension",dim=3)
		html_main.add_content("前n个维度累积解释变异的百分比图(见下图)。其中，前三个维度，累积解释变异的百分比: %.2f%%, %.2f%%, %.2f%%"%(mdsout.p[0],mdsout.p[1],mdsout.p[2]))
	html_main.add_content("""<img src="./Variation_percentage.png" width="50%" /><a href="./Variation_percentage.svg">SVG矢量图版本</a>""")
	html_main.add_enter()
	html_main.add_content("""样本在前三个维度中的空间分布图""")
	html_main.add_content("""<img src="./MDS_samples_distribution.png" width="50%" /><a href="./MDS_samples_distribution.svg">SVG矢量图版本</a>""")
	html_main.add_enter()


	f = file("exprs_samples_qc.html","w")
	f.write(str(html_main))
	f.close()
	return 0
