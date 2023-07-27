from setuptools import setup, find_packages

setup(name='ContessUtils',
	version = "0.1.0",
	packages = find_packages(include = [ "ContessUtils", "ContessUtils.*" ]),

	author = "Shayan Hoshyari, Chenxi Liu",
	description = "",

	url="",

	project_urls=dict(),

	license = "PRIVATE",

	classifiers=[
	],

	python_requires=">=3.7",


	entry_points = {
		'console_scripts': [
			# 'banach=Banach.Scripts.main:main_cmd',
		]
	},

)
