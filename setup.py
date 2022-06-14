import pathlib
from setuptools import setup, find_packages

#from halfsat import __version__

HERE = pathlib.Path(__file__).parent

#VERSION = __version__
VERSION = '1.0.0'
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

docs_extras = [
    'Sphinx >= 3.0.0',  # Force RTD to use >= 3.0.0
    'docutils',
    'pylons-sphinx-themes >= 1.0.8',  # Ethical Ads
    'pylons_sphinx_latesturl',
    'repoze.sphinx.autointerface',
    'sphinxcontrib-autoprogram',
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
      packages=find_packages(),
      extras_require={'docs': docs_extras}
      )
