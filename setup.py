from setuptools import setup


__version__ = 'unknown'
exec(open('treeomics/version.py').read())

setup(
      name='treeomics',                                  # package name
      # packages=find_packages(),
      packages=['treeomics', 'tests'],  # 'treeomics.phylogeny', 'treeomics.plots', 'treeomics.utils'],
      version=__version__,
      description='Treeomics infers metastatic seeding patterns from DNA sequencing data.',
      install_requires=['numpy', 'scipy', 'pandas', 'matplotlib<2.0', 'seaborn', 'networkx<2.0',
                        'wkhtmltopdf', 'pdfkit'],
      url='https://github.com/johannesreiter/treeomics',
      author='Johannes Reiter',
      author_email='johannes.reiter@stanford.edu',
      license='GNUv3',
      classifiers=[
        'Programming Language :: Python :: 3.4',
      ],
      test_suite='tests',
      entry_points={
        'console_scripts': [
            'treeomics = treeomics.__main__:main'
        ]
      }
)
