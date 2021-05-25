from setuptools import setup

extra_requirements = {
    "tests": ["pytest", "coverage", "pytest-cov"],
    "docs": ["sphinx", "sphinx-autodoc-typehints"],
}

setup(
    name="framol",
    author="Sarai D. Folkestad",
    description="Description of package",
    install_requires=["numpy", "scipy"],
    extras_require=extra_requirements,
)
