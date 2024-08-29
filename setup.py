from setuptools import setup, find_packages

setup(
    name="scrna_atlas",
    version="0.1",
    packages=find_packages(),
    install_requires=[
        "numpy",
        "pandas",
        "s3fs",
        "anndata",
        "scanpy",
    ],
)
