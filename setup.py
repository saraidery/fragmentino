from setuptools import setup

extra_requirements = {
    "tests": ["pytest", "coverage", "pytest-cov"],
    "docs": [
        "sphinx",
        "sphinx-autodoc-typehints",
        "furo",
        "sphinx-copybutton",
        "nbsphinx",
        "pandoc",
    ],
}

setup(
    name="fragmentino",
    author="Sarai D. Folkestad",
    description="Description of package",
    install_requires=["numpy", "scipy", "plotly"],
    extras_require=extra_requirements,
)
