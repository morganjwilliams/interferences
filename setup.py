from setuptools import setup, find_packages
import versioneer

tests_require = ['pytest',
                 'pytest-runner',
                 'pytest-cov',
                 'coverage',
                 'coveralls']

dev_require = ['versioneer',
               'nbstripout',
               'nbdime',
               'twine',
               'sphinx_rtd_theme',
               'sphinx-autodoc-annotation',
               #'mock',
               ]

with open('README.md', 'r') as src:
    LONG_DESCRIPTION = src.read()

setup(name='interferences',
      description="Prediction of elemental interferences in ICP-MS.",
      long_description=LONG_DESCRIPTION,
      long_description_content_type='text/markdown',
      version=versioneer.get_version(),
      url='https://github.com/morganjwilliams/interferences',
      author='Morgan Williams',
      author_email='morgan.williams@csiro.au',
      classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Natural Language :: English',
        'Programming Language :: Python :: 3 :: Only',
        'Topic :: Scientific/Engineering',
        'Topic :: Software Development :: Libraries :: Python Modules',
        ],
      packages=find_packages(exclude=['test*']),
      install_requires=['pathlib',
                        'numpy',
                        'scipy',
                        'pandas',
                        'matplotlib',
                        'periodictable',
                        ],
      extras_require={'dev': dev_require + tests_require,
      },

      tests_require=tests_require,
      test_suite="test",
      package_data={'interferences': ['data/*']},
      include_package_data=True,
      license='CSIRO Modifed MIT/BSD',
      cmdclass=versioneer.get_cmdclass()
)
