from setuptools import setup, find_packages

setup(
    name="NM_interaction",
    version="1.0.0",
    packages=find_packages(),
    include_package_data=True,  # Include files specified in MANIFEST.in
    install_requires=[
        "pandas",
        "matplotlib",
        "numpy",
    ],
)