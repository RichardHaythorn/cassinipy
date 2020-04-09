import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="cassinipy-richardhaythorn", # Replace with your own username
    version="0.1.0",
    author="Richard Haythornthwaite",
    author_email="richardhaythornthwaite@hotmail.co.uk",
    description="Code for analysing Cassini in-situ data",
    url="https://github.com/RichardHaythorn/cassinipy",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
    ],
    python_requires='>=3.6',
)