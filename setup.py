from setuptools import setup


__version__ = 'unknown'
exec(open('treeomics/version.py').read())

setup(
      name='treeomics',                                  # package name
      # packages=find_packages(),
      packages=['treeomics', 'tests'],  # 'treeomics.phylogeny', 'treeomics.plots', 'treeomics.utils'],
      version=__version__,
      description='Treeomics infers metastatic seeding patterns from DNA sequencing data.',
      install_requires=['numpy', 'scipy', 'pandas', 'matplotlib==1.5', 'seaborn==10', 'networkx<2.0',
                        'wkhtmltopdf', 'pdfkit'],
      setup_requires=['pytest-runner', 'flake8'],
      tests_require=['pytest', 'pytest-cov'],
      extras_require={'plotting': ['matplotlib', 'seaborn', 'jupyter']},
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
