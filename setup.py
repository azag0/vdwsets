from setuptools import setup


setup(
    name='vdwsets',
    version='0.1',
    description='van der Waals datasets',
    author='Jan Hermann',
    author_email='dev@janhermann.cz',
    url='https://github.com/azag0/vdwsets',
    packages=['vdwsets'],
    package_data={'vdwsets': ['data/*/energies.csv', 'data/*/geoms/*.xyz']},
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: Mozilla Public License 2.0 (MPL 2.0)',
        'Natural Language :: English',
        'Operating System :: Unix',
        'Programming Language :: Python :: 3.6',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Topic :: Scientific/Engineering :: Physics',
    ],
    license='Mozilla Public License 2.0',
    install_requires=['pandas'],
)
