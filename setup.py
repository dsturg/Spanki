from setuptools import setup, find_packages
from distutils.extension import Extension
try:
    import numpy
except ImportError:
    print "-------------------------------------------------"
    print "Please install numpy"
    print "-------------------------------------------------"

setup(
        name='spanki',
        version='0.5.0',
        #py_modules=['modules/spanki_parse_utils,modules/spanki_utils'],
        packages = ['spanki'],
        package_data={'spanki': ['data/*.txt']},
        include_package_data=True,
   		scripts = ['bin/junccomp','bin/merge_jtabs','bin/spankisim_transcripts','bin/spankijunc','bin/spankisplice','bin/splicecomp','bin/quickjunc','bin/spankisim_models','bin/make_curated_jtab','bin/annotate_junctions'],
        author='David Sturgill',
        description="A splicing analysis toolkit",
        author_email='dave.sturgill@gmail.com',
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
