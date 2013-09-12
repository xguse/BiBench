.. include:: global.txt

============
Installation
============

Instructions for installing BiBench and its dependencies. BiBench is
written for Python 2. It has been tested with both Python 2.6 and 2.7.

----------------------
Installation Procedure
----------------------

++++++++++++++++++++++++++++++++
Automatic Installation with pip
++++++++++++++++++++++++++++++++

`pip <http://www.pip-installer.org/>`_ and `virtualenv
<http://www.virtualenv.org>`_ are currently the best way
to install and manage Python packages. This method is also recommended
because pip should automatically install BiBench's
dependencies. (Note: One of BiBench's dependencies is rpy2, which
requires that `R <http://www.r-project.org/>`_ be installed. R must
have been built as a library. i.e., configured with
``--enable-R-shlib``, and libR.so must be available in
LD_LIBRARY_PATH. For more information, see the `rpy FAQ
<http://rpy.sourceforge.net/rpy_faq.html>`_).

If you want pip to install BiBench without installing dependencies,
use ``pip install --no-deps`` in the following steps.

If it is not installed, install pip. If necessary, you can use pip to install
virtualenv and virtualenvwrapper::

    pip install virtualenv virtualenvwrapper

To use virtualenvwrapper you must source virtualenvwrapper.sh, so that
its commands are available in Bash. One option is to add the following
to your .bashrc::

    source `which virtualenvwrapper.sh`

To create a new virtual environment called bibench_env and install
BiBench and its dependencies to it:

.. parsed-literal::

    mkvirtualenv --no-site-packages bibench_env
    pip install |downloadurl|

Pip can also install from a local package:

.. parsed-literal::

    wget |downloadurl|
    pip install BiBench-0.1.tar.gz

BiBench is now installed and ready to go.

The BiBench package is now installed into its own virtual environment,
``bibench_env``, which is an isolated Python environment. So before
using BiBench you must switch to that environment::

    workon bibench_env

When you are done working in that environment, simply run the
deactivate command::

    deactivate

    
+++++++++++++++++++
Manual installation
+++++++++++++++++++

It is possible to manually install BiBench. In this case, all of
BiBench's Python dependencies must be installed seperately.

BiBench is packaged using `distutils
<http://docs.python.org/library/distutils.html>`_.  To install BiBench
manually, simply download it, unpack it, and run the distutils setup
script (depending on your Python setup, you may need to be root for
this step):

.. parsed-literal::
    wget |downloadurl|
    tar xzf BiBench-0.1.tar.gz
    cd BiBench-0.1
    python setup.py install


+++++++++++++++++++++++++++
Building the documentation
+++++++++++++++++++++++++++

This documentation may be compiled into a number of formats, using
sphinx. To generate html::

    make html

For a list of all possible targets, run::

    make help

------------
Dependencies
------------

The BiBench package is now installed, but most of its functionality
will not be available until extra dependencies are available.

In particular, BiBench does not provide implementations for biclustering algorithms; they must be installed seperately.

+++++++++++++++++++
Python Dependencies
+++++++++++++++++++

* `NumPy <http://numpy.scipy.org/>`_
* `rpy2 <http://rpy.sourceforge.net/rpy2.html>`_
* `decorator <http://pypi.python.org/pypi/decorator>`_
* `nose <http://code.google.com/p/python-nose/>`_ (optional: for running unit tests)
* `sphinx <http://sphinx.pocoo.org/>`_ (optional: to make the documentation)

If you install BiBench using pip, you should not need to install these
packages manually; pip should automatically handle dependencies.

If you chose not to install BiBench's dependencies before, you can install them with pip::

    pip install numpy rpy2 decorator

To install optional dependencies simply run::

    pip install nose sphinx


++++++++++++++
R Dependencies
++++++++++++++

Much of BiBench's functionality depends on `R
<http://www.r-project.org/>`_. Some algorithms that BiBench supports
are available as packages for R; BiBench also relies on R for some
visualization methods and other functionality. To communicate with R,
the Python package rpy2 must be installed.

rpy2 requires that R be compiled as a shared library using the
configure option ``--enable-R-shlib``. Also libR.so must be available
in LD_LIBRARY_PATH. For more information, see the `rpy FAQ
<http://rpy.sourceforge.net/rpy_faq.html>`_.

Assuming that R is installed, entering the following commands **in R**
should install all of BiBench's R dependencies::

    install.packages(c('biclust', 'isa2', 'MASS'))
    source("http://bioconductor.org/biocLite.R")
    biocLite()
    biocLite(c('fabia',
               'GEOquery',
	       'GEOmetadb',
	       'GOstats',
	       'GO.db',
	       'multtest',
	       'pcaMethods'))

The the rest of this section provides links to those dependencies.

The following may be installed with the R command ``install.packages()``:

* `biclust <http://cran.r-project.org/package=biclust>`_
* `isa2 <http://cran.r-project.org/web/packages/isa2/index.html>`_
* `MASS <http://cran.r-project.org/web/packages/MASS/index.html>`_

The rest of the packages require `Bioconductor
<http://www.bioconductor.org/>`_. Here are the `installation
directions <http://www.bioconductor.org/install/>`_ for Bioconductor
and its packages.

* `fabia <http://www.bioconductor.org/packages/release/bioc/html/fabia.html>`_

For downloading GDS data:

* `GEOquery <http://www.bioconductor.org/packages/1.8/bioc/html/GEOquery.html>`_
* `GEOmetadb <http://www.bioconductor.org/packages/2.2/bioc/html/GEOmetadb.html>`_

For `Gene Ontology <http://www.geneontology.org/>`_ enrichment analysis:

* `GOstats <http://www.bioconductor.org/packages/release/bioc/html/GOstats.html>`_
* `GO.db <http://www.bioconductor.org/packages/release/data/annotation/html/GO.db.html>`_
* `multtest <http://www.bioconductor.org/packages/release/bioc/html/multtest.html>`_
* The correct `AnnotationData <http://www.bioconductor.org/packages/2.6/data/annotation/>`_ package for your organism.

For missing data imputation:

* `pcaMethods <http://www.bioconductor.org/packages/1.9/bioc/html/pcaMethods.html>`_


++++++++++++++++++++++++++++
Other Algorithm Dependencies
++++++++++++++++++++++++++++

These algorithms are not available for R or Python. They must be
manually built and installed, and their binaries must be available on
the PATH, in order for the appropriate module in
``bibench.algorithms`` to work.

* `BBC <http://www.people.fas.harvard.edu/~junliu/BBC/>`_
* `COALESCE <http://huttenhower.sph.harvard.edu/sleipnir/>`_ (as part of the Sleipnir software)
* `CPB <http://bmi.osu.edu/hpc/software/cpb>`_
* `OPSM <http://bmi.osu.edu/hpc/software/OPSM.tar.gz>`_ (modified for command-line use from `BicAT <http://www.tik.ethz.ch/sop/bicat/>`_)
* `QUBIC <http://csbl.bmb.uga.edu/~maqin/bicluster/>`_
