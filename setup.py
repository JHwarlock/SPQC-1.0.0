import sys
from distutils.core import setup

long_description = \
"""sample quality control for exprs data.
"""

def get_setup_args():
	#SPQC_VERSION = None
	py_modules = ["exprsQC"]
	setup_args = dict(
			name = "exprsQC",
			version = "1.0.0",
			description = long_description,
			author = "Rong ZhengQin",
			author_email = "rongzhengqin@honortech.cn",
			platforms = ["Linux","Mac OS-X","UNIX"],
			#url = "",
			# Description of the modules and packages in the distribution
			package_dir = {"exprsQC": "lib"},
			packages = ["exprsQC"],
			#scripts=[],
			license = "BSD",
			data_files = [],
			ext_modules = [],
			#entry_points = {
			#	    'console_scripts' : [
			#			'SPQC = '],
			#	}
			scripts=['bin/SPQC'],
			)
	return setup_args

if __name__ == '__main__':
	setup(**get_setup_args())

