###################################################################################################
###---------------------------------------------------------------------------------------------###
### This file is part of bibench (Biclustering Benchmarking)			                        ###
### Copyright (c) 2011,                                                                         ###
### By:    Kemal Eren,                                                                          ###
##         Mehmet Deveci,                                                                       ###
###        Onur Kucuktunc,                                                                      ###
###        Umit V. Catalyurek                                                                   ###
###                                                                                             ###
###---------------------------------------------------------------------------------------------###
### For license info, please see the README.txt and LICENSE.txt files in the main directory.    ###
###---------------------------------------------------------------------------------------------###
###################################################################################################
# -*- coding: utf-8 -*-
from distutils.core import setup

setup(
    name='BiBench',
    version='0.1',
    packages=['bibench',
              'bibench.algorithms',
              'bibench.datasets',
              'bibench.validation',
              'bibench.test'],

    scripts=['bin/generate_dataset.py',
             'bin/run_algorithm.py',
             'bin/util.py'],

    install_requires=['numpy',
                      'rpy2',
                      'decorator'],

    extras_require = {
        "test" : ["nose"],
        "doc" : ["sphinx"]
        },

    #metadata
    author=u'Kemal Eren, Mehmet Deveci, Umit Catalyurek',
    author_email='ekemal@bmi.osu.edu, mdeveci@bmi.osu.edu, umit@bmi.osu.edu',
    description='Biclustering framework',
    long_description=open('README').read(),
    license='GPL license, see LICENSE',
    keywords = 'biclustering',
    url='',
    download_url = '',
)
