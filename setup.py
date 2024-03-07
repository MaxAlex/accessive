from setuptools import setup, find_packages

setup(
    name="Accessive",
    version="0.1.0",
    author="William Max Alexander",
    author_email="accessive@alexander.bio",
    description="A Python library for converting between different bioinformatic accession types",
    long_description=open('README.md').read(),
    long_description_content_type="text/markdown",
    url="https://github.com/maxalex/accessive",
    packages=find_packages(),
    install_requires=[
        'requests',
        'pandas>=2.1.0'
    ],
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: MIT License',
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
    ],
    python_requires='>=3.6',
    package_data={"accessive": ["db_build/*"]},
    include_package_data=True,
)
