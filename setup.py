from setuptools import setup, find_packages

setup(
    name='scrna_atlas',
    version='0.1',
    packages=find_packages(),
    install_requires=[
        'numpy',
        'pandas',
        's3fs',
        'anndata',
        'scanpy',
    ],
    # author='Your Name',
    # author_email='your@email.com',
    # description='A simple package for single-cell RNA sequencing data exploration',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
)
