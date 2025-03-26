def get_bibtex():
    bibtex_entries = {
        "trimesh": """@article{trimesh,
    author  = {Dawson-Haggerty, M. and others},
    title   = {trimesh: Python library for loading and using triangular meshes},
    year    = {2019},
    note    = {Available at \\url{https://github.com/mikedh/trimesh}}
}""",
        "scikit-image": """@article{scikit-image,
    author  = {Van der Walt, S. and Sch\\"{o}nberger, J. L. and Nunez-Iglesias, J. and Boulogne, F. and Warner, J. D. and others},
    title   = {scikit-image: image processing in Python},
    journal = {PeerJ},
    volume  = {2},
    pages   = {e453},
    year    = {2014},
    doi     = {10.7717/peerj.453}
}""",
        "scipy": """@article{scipy,
    author  = {Virtanen, P. and Gommers, R. and Oliphant, T. E. and Haberland, M. and Reddy, T. and others},
    title   = {{SciPy} 1.0: fundamental algorithms for scientific computing in Python},
    journal = {Nature Methods},
    volume  = {17},
    pages   = {261--272},
    year    = {2020},
    doi     = {10.1038/s41592-019-0686-2}
}"""
    }
    return "\n\n".join(bibtex_entries.values())
