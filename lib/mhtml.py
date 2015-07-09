# -*- coding: UTF-8 -*-
import sys


class table_html(object):
	def __init__(self,listhead,table_title=""):
		#self.head = """<table frame="hsides" border="1" cellpadding="2" align="left" cellspacing="0">\n"""+"""<caption align="top">%s</caption>"""%table_title
		self.head = """<table frame="hsides" border="1" cellpadding="2" cellspacing="0">\n"""+"""<caption align="top">%s</caption>"""%table_title
		self.end = """</table>\n"""
		self.h_head = """<tr>"""
		self.h_end = """</tr>"""
		self.v_head = """<td>"""
		self.v_end = """</td>"""
		self.listhead = listhead[:]
		self.length = len(self.listhead)
		self.content = []
		strcontent = "<tr>"
		for i in xrange(self.length):
			strcontent += "<th>"+str(listhead[i])+"</th>"
		strcontent += "</tr>"
		self.content.append(strcontent)
	def add_row(self,content):
		assert len(content) == self.length
		strcontent = "<tr>"
		for i in xrange(self.length):
			strcontent += "<td><pre>"+str(content[i])+"</pre></td>"
		strcontent += "</tr>"
		self.content.append(strcontent)
	def add_noprerow(self,content):
		assert len(content) == self.length
		strcontent = "<tr>"
		for i in xrange(self.length):
			strcontent += "<td>"+str(content[i])+"</td>"
		strcontent += "</tr>"
		self.content.append(strcontent)
	def __str__(self):
		return "\n".join([self.head,"\n".join(self.content),self.end])

class simple_main(object):
	def __init__(self,title="Result report",css="CSS"):
		#start = """<html><body><h1>KEGG Pathway Enrichment Analysis</h1><hr />"""
		self.start = """<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN"  "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html>
<head>
<meta http-equiv="Content-Type" content="text/html; charset=utf-8">
<title>%s</title>\n
<link href="%s/report.css" rel="stylesheet" type="text/css">
</head>
<body>"""%(title,css) ###also can add 
#<link href="Files/Pages/css/report.css" rel="stylesheet" type="text/css">
#<script language="javascript" src="Files/Pages/js/jquery.js"></script>
#<script language="javascript" src="Files/Pages/js/report.js"></script>
		self.end = """</body></html>"""
		self.content = ""
	def add_head(self,content,headlevel=1):
		self.content += "<h%d>%s</h%d>\n"%(headlevel,content,headlevel)
	def add_content(self,content):#content is str
		self.content += "<p>"+content+"</p>\n"
	def add_precontent(self,content):
		self.content += "<p><pre>"+content+"</pre></p>\n"
	def add_line(self):
		self.content += "<hr />\n"
	def add_back1(self,content="返回上一页"):
		self.content += """<input type="button" name="Submit" onclick="javascript:history.back(-1);" value="返回上一页">\n"""
		#self.content += """<a href="<a href="javascript :history.back(-1)">返回上一页</a>\n"""
	def add_enter(self):
		self.content += "<br />\n"
	def __str__(self):
		if self.content == "":
			sys.stderr.write("[ERROR] html page has no content!\n")
			return None
		return "\n".join([self.start,self.content,self.end])

def get_sampleinfo(sampleinfo):
	f = file(sampleinfo,"r")
	header = ["#","Samplename","Filename","SN","Category","Category name"]
	for line in f:
		if line.startswith("##"):continue
		if line.startswith("#"):
			header = ["#",] + line.rstrip("\n")[1:].split("\t")
			break
	f.seek(0)
	#print header
	table = table_html(header,"Table: sample information list")
	idx = 1
	samplename = []
	files = []
	for line in f:
		if line.startswith("##"):continue
		if line.startswith("#"):
			header = ["#",] + line.rstrip("\n")[1:].split("\t")
			continue
		arr = line.rstrip("\n").split()
		sn = arr[0]
		files.append(arr[1].split(","))
		#print [str(idx),"""<a href="#%s">%s</a>"""%(sn,sn)] + arr[1:]
		table.add_row([str(idx),"""<a href="#%s">%s</a>"""%(sn,sn)] + arr[1:])
		idx += 1
		samplename.append(sn)
	return str(table),samplename,files

def xls2table(fxls,title=None,header=None):
	f = file(fxls,"r")
	if header == None:
		header = "";
	if not title:
		title = "";
	for line in f:
		if line.startswith("##"):
			title = "Table: "+line[2:].rstrip()
			continue
		elif line.startswith("#"):
			header = line.rstrip("\n").split("\t")
			continue
		else:break
	f.seek(0)
	table = table_html(header,title)
	note = ""
	while 1:
		line = f.readline()
		if not line: break
		if line.startswith("#"):continue
		if line.startswith("Note:"):
			note = line+f.read()
			break
		table.add_row(line.rstrip("\n").split("\t"))
	f.close()
	return str(table),note
if __name__  == "__main__":
	print get_sampleinfo(sys.argv[1])
