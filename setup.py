import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="lumpyrem", 
    version="0.1.1",
    author="Rui Hugman",
    author_email="rthugman@gmail.com",
    description="This package reads, writes and runs LUMPREM2 and associated utility programs",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/rhugman/lumpyrem",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)