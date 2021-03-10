from setuptools import setup


__version__ = 'unknown'
exec(open('treeomics/version.py').read())

setup(
      name='treeomics',                                  # package name
      packages=setuptools.find_packages(),
      #packages=['treeomics', 'tests'],  # 'treeomics.phylogeny', 'treeomics.plots', 'treeomics.utils'],
      version=__version__,
      description='Treeomics infers metastatic seeding patterns from DNA sequencing data.',
      install_requires=['numpy', 'scipy', 'pandas', 'matplotlib', 'seaborn', 'networkx<2.0',
                        'wkhtmltopdf', 'pdfkit', 'ete3'],
      setup_requires=['pytest-runner', 'flake8'],
      tests_require=['pytest', 'pytest-cov'],
      extras_require={'plotting': ['matplotlib', 'seaborn', 'jupyter', 'etetoolkit']},
      url='https://github.com/reiterlab/treeomics',
      author='Johannes Reiter',
      author_email='johannes.reiter@stanford.edu',
      license='GNUv3',
      classifiers=[
        'Programming Language :: Python :: 3.6',
      ],
      test_suite='tests',
      entry_points={
        'console_scripts': [
            'treeomics = treeomics.__main__:main'
        ]
      }
)
