from setuptools import setup

setup(name='tfea',
      version='3.0',
      description='Transcription Factor Enrichment Analysis',
      url='https://github.com/jdrubin91/TFEA.git',
      author='Jonathan Rubin',
      author_email='jonathan.rubin@colorado.edu',
      license='CU Boulder Dowell Lab',
      packages=['TFEA'],
      scripts=['bin/TFEA', 'bin/TFEA-annotate'],
      install_requires=[
          'htseq',
          'configparser',
          'argparse',
          'matplotlib',
          'scipy',
          'numpy',
          'pybedtools',
          'psutil'
      ],
      zip_safe=False)
