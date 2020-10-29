import glob
from setuptools import setup

setup(name='biota',
      packages = ['biota'],
      version='0.2',
      description='Tools to generate aboveground biomass forest change maps from ALOS PALSAR/PALSAR-2 mosaic data.',
      url='https://bitbucket.org/sambowers/biota',
      data_files=[('./cfg/',glob.glob('./cfg/*'))],
      author='Samuel Bowers',
      author_email='sam.bowers@ed.ac.uk',
      license='GNU General Public License',
      zip_safe=False)
