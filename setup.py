import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="dyntripy", # Replace with your own username
    version="3.2.3",
    author="Naidan YUN",
    author_email="yunnaidan@gmail.com",
    description="A package for detecting dynamic triggering automatically",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/yunndlalala/dynamic_earthquake_triggering",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires=[
        "obspy>=1.1.0",
        "matplotlib>=3.1.2",
        "pandas>=0.24.2",
        "scipy>=1.1.0",
        "numpy>=1.16.3"
    ],
    python_requires='>=3.6.2',
)
