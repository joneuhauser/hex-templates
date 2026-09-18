"""Shared static navigation for the mesh atlas and its page generators."""
from html import escape

PAGES=[('index.html','Overview'),('pyramid.html','Pyramid'),('geode.html','Geode'),
       ('geode-sides.html','Side walls'),('periodic.html','Periodic rotation'),('methodology.html','Methodology')]


def header(current):
    links=''.join(f'<a href="{path}"'+(' aria-current="page"' if current==path else '')+f'>{escape(label)}</a>' for path,label in PAGES)
    return ('<header class="site-header"><a class="brand" href="index.html">'
            '<img src="assets/figures/mark.svg" width="25" height="28" alt="">Hex template library</a>'
            '<nav aria-label="Main navigation">'+links+'</nav></header>')


def search_summary(rows):
    """Render labeled search settings; descriptions contain trusted page HTML."""
    return ('<section class="search-summary" aria-labelledby="search-setup">'
            '<h2 id="search-setup">Search setup and references</h2><dl>'
            + ''.join(f'<dt>{escape(label)}</dt><dd>{description}</dd>'
                      for label, description in rows)
            + '</dl></section>')
