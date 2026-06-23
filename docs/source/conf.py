import os
import sys
sys.path.insert(0, os.path.abspath('../..'))

from mnished._version import __version__

project = 'MNiShed'
copyright = '2019-2026, MNiMORPH'
author = 'Andrew D. Wickert'
release = __version__

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.intersphinx',
    'sphinx.ext.viewcode',
    'sphinx.ext.mathjax',
    'sphinx.ext.graphviz',
    'sphinx_rtd_theme',
]

# Render Graphviz diagrams (e.g. the ET-pathway figure) as scalable SVG.
graphviz_output_format = 'svg'

# Keep the sidebar/TOC to section headings only — do not add every Python
# object (class/method) to the nav. This stops the API reference from nesting
# the ``Buckets`` class under the ``Buckets`` section heading. (Sphinx >= 5.2.)
toc_object_entries = False

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
html_theme_options = {
    'logo_only': False,
    'prev_next_buttons_location': 'bottom',
    'style_external_links': False,
    'collapse_navigation': True,
    'sticky_navigation': True,
    'navigation_depth': 4,
}

# NumPy-style docstrings
napoleon_numpy_docstring = True
napoleon_google_docstring = False
napoleon_include_init_with_doc = True

autodoc_member_order = 'bysource'

intersphinx_mapping = {
    'python': ('https://docs.python.org/3', None),
    'numpy': ('https://numpy.org/doc/stable/', None),
    'pandas': ('https://pandas.pydata.org/docs/', None),
}
