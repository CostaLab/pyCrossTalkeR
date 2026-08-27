
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
    'nbsphinx', 
    'sphinx_design',
    "myst_parser", 
    'sphinx_book_theme' ,
    "sphinx_rtd_size",
]
myst_enable_extensions = []
autosummary_generate = True # Auto-generates individual API doc pages from summary tables

# -- BibTeX citation settings
bibtex_bibfiles = ['references.bib']
bibtex_default_style = 'alpha'

intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'sphinx': ('https://www.sphinx-doc.org/en/master/', None),
}
intersphinx_disabled_domains = ['std']
templates_path = ['_templates']

sphinx_rtd_size_width = "85%"


# -- It hides the input/output prompt numbers (In [1]: / Out [1]:).
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

# -- Options for HTML output -------------------------------------------------

html_theme = 'sphinx_book_theme'
html_logo = "_static/logo.png"


# Specific parameters passed directly to the 'sphinx_book_theme'.
html_theme_options = {
    # URL of the GitHub repository
    "repository_url": "https://github.com/CostaLab/pyCrossTalkeR",
    # GitHub button in the top navigation bar
    "use_repository_button": True,
    # Show a download button (e.g., PDF or Markdown/RST)
    "use_download_button": True,
    # Full-screen reading mode
    "use_fullscreen_button": True,
    # Collapse inactive subsections in the left sidebar to keep it tidy
    "collapse_navbar": True,
    # Set the maximum depth of heading levels (H1 to H4) shown in the right-hand page TOC
    "show_toc_level": 4,
    # Control how many navigation levels deep are automatically expanded in the left sidebar
    "show_navbar_depth": 2,
    "logo": {
        "text": "<b>pyCrossTalkeR</b>",
    }
}

# -- Static files configuration
html_js_files = [
    "require.min.js",
    "custom.js",
]

html_static_path = ["_static"]

# Link custom CSS (logo sizing and styling, located in '_static')
html_css_files = [
    "custom.css",
]