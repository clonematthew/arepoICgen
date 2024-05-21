import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="arepoICgen",
    author="Matt Cusack",
    author_email="cusackmt@cardiff.ac.uk",
    description="Package for generating initial conditions for AREPO",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/clonematthew/arepoICgen",
    package_dir={'':'src'},
    packages=setuptools.find_packages(where='src')
)