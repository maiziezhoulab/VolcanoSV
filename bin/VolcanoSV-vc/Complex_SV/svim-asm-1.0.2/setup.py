import os
from setuptools import setup, find_packages

here = os.path.abspath(os.path.dirname(__file__))

# Get the long description from the relevant file
with open(os.path.join(here, 'README.rst')) as f:
      long_description = f.read()

def get_version(string):
      """ Parse the version number variable __version__ from a script. """
      import re
      version_re = r"^__version__ = ['\"]([^'\"]*)['\"]"
      version_str = re.search(version_re, string, re.M).group(1)
      return version_str


setup(name='svim-asm',
      version=get_version(open('src/svim_asm/svim-asm').read()),
      description='A structural variant caller for genome-genome alignments.',
      long_description=long_description,
      url='https://github.com/eldariont/svim-asm',
      author='David Heller',
      author_email='heller_d@molgen.mpg.de',
      license='GPLv3',
      classifiers=[
      'Development Status :: 5 - Production/Stable',
      'Environment :: Console',
      'Intended Audience :: Science/Research',
      'Topic :: Scientific/Engineering :: Bio-Informatics',
      'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
      'Programming Language :: Python :: 3.6'
      ],
      keywords='svim-asm SV assembly structural variation caller',
      packages = find_packages("src"),
      package_dir = {"": "src"},
      data_files = [("", ["LICENSE"])],
      zip_safe=False,
      install_requires=['pysam', 'numpy', 'scipy', 'matplotlib', 'edlib'],
      scripts=['src/svim_asm/svim-asm'])
