from setuptools import setup, find_packages
from distutils.extension import Extension
try:
    import numpy
except ImportError:
    print "-------------------------------------------------"
    print "Please install numpy"
    print "-------------------------------------------------"

setup(
        name='spanky',
        version='0.1.2',
        #py_modules=['modules/spanky_parse_utils,modules/spanky_utils'],
        packages = ['spanky'],
        package_data={'spanky': ['data/*.txt']},
        include_package_data=True,
   		scripts = ['bin/astacomp','bin/junccomp','bin/merge_jtabs_all','bin/sim_transcripts','bin/spankyjunc','bin/spankysplice','bin/splicecomp','bin/quickjunc','bin/simrnaseq_make_models'],
        author='David Sturgill',
        description="A splicing analysis toolkit",
        author_email='davidsturgill@niddk.nih.gov',
        url='none',
        install_requires=[
        	"fisher >= 0.1.4",
         	"scipy >= 0.1.10",
        	"scikits.statsmodels >= 0.3.1",
        	"biopython >= 1.50",
        	"pyfasta >= 0.4.4",
        	"argparse >= 1.2.1",
         	"pysam >= 0.5"]
    )
