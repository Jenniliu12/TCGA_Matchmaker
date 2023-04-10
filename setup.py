from setuptools import setup

setup(name = "TCGA_Matchmaker",
	version = "0.1.0",
	packages = ["TCGA_code"], 
    entry_points = {
		'console_scripts': [
		'test_run = demo_pkg.test:main',
        'another_module_run = demo_pkg.another_module:main'
		]
	}
	)