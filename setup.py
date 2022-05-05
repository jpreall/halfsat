import pathlib
from setuptools import setup, find_packages

from halfsat import __version__

HERE = pathlib.Path(__file__).parent

VERSION = __version__
PACKAGE_NAME = 'halfsat'
AUTHOR = 'Jonathan Preall and Brian He'
AUTHOR_EMAIL = 'jpreall@cshl.edu; bhe@cshl.edu'
URL = 'https://github.com/jpreall/halfsat'

LICENSE = 'BSD 3-Clause License'
DESCRIPTION = 'Scrapes information from web_summary.html. Creates half-sat curves and tables for informational reference.'
LONG_DESCRIPTION = (HERE / "README.md").read_text()
LONG_DESC_TYPE = "text/markdown"

INSTALL_REQUIRES = [
      'bs4',
      'matplotlib',
      'argparse',
      'scanpy',
      'h5py',
      'numpy',
      'pandas',
      'scipy'
]

setup(name=PACKAGE_NAME,
      version=VERSION,
      description=DESCRIPTION,
      long_description=LONG_DESCRIPTION,
      long_description_content_type=LONG_DESC_TYPE,
      author=AUTHOR,
      license=LICENSE,
      author_email=AUTHOR_EMAIL,
      url=URL,
      install_requires=INSTALL_REQUIRES,
      packages=find_packages()
      )
