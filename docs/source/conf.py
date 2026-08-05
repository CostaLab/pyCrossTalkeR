
# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information


project = 'pyCrossTalkeR'
copyright = '2025, James Nagai'
author = 'James Nagai'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration


version = '0.1.8.8'

extensions = [
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.intersphinx',
    "sphinx.ext.napoleon",
    'sphinxcontrib.bibtex',
    'sphinx_design',
    "myst_parser",
    "nbsphinx",
    "sphinx_rtd_size",
]

autosummary_generate = True
bibtex_bibfiles = ['references.bib']
bibtex_default_style = 'alpha'

intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'sphinx': ('https://www.sphinx-doc.org/en/master/', None),
}
intersphinx_disabled_domains = ['std']

templates_path = ['_templates']

#html_theme = 'sphinx_rtd_theme'
html_theme = 'sphinx_book_theme'


html_logo = "_static/logo.png"

html_theme_options = {
    "repository_url": "https://github.com/CostaLab/pyCrossTalkeR",
    "use_repository_button": True,
    "use_download_button": True,
    "use_fullscreen_button": True,
    "collapse_navbar": True,
    "show_toc_level": 4,
    "show_navbar_depth": 2,

    "logo": {
        "text": "<b>pyCrossTalkeR</b>",
    }
}


sphinx_rtd_size_width = "85%"

nbsphinx_prolog = """
.. raw:: html

    <style>
        div.nbinput.container div.prompt,
        div.nboutput.container div.prompt,
        span.prompt {
            display: none !important;
            min-width: 0 !important;
            padding: 0 !important;
        }
    </style>
"""


# -- Options for EPUB output
epub_show_urls = 'footnote'    

html_js_files = [
    "require.min.js",
    "custom.js",
]

html_static_path = ["_static"]

html_css_files = [
    "custom.css",
]