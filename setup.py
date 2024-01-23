from setuptools import setup, find_packages

setup(
    name='Accessive',
    version='0.1.0',
    author='William Max Alexander',
    author_email='max@alexander.bio',
    description='A bioinformatics accession mapping tool',
    url='https://github.com/maxalex/Accessive',
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        'pandas',
        'psycopg2-binary'
    ],
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
    ],
    python_requires='>=3.7',
)
