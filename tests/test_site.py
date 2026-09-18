"""Static Pages contract: links, catalog paths, and archive contents agree."""
from pathlib import Path
from html.parser import HTMLParser
from urllib.parse import urlsplit, unquote
import hashlib, json, zipfile
import xml.etree.ElementTree as ET
ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / 'docs'
class Page(HTMLParser):
    def __init__(self, text):
        super().__init__(); self.links=[]; self.ids=set(); self.feed(text)
    def handle_starttag(self, tag, attrs):
        a=dict(attrs)
        if 'id' in a:
            assert a['id'] not in self.ids, a['id']
            self.ids.add(a['id'])
        for key in ('href','src'):
            if key in a:self.links.append(a[key])
pages={p:Page(p.read_text()) for p in DOCS.glob('*.html')}
for name, catalog_path in [('pyramid', 'catalog.json'), ('geode', 'geode/catalog.json'),
                           ('periodic', 'periodic/current-catalog.json'),
                           ('refinement', 'refinement/catalog.json')]:
    # Catalog cards and periodic class headings are inserted by catalog.js.
    page = pages[DOCS/f'{name}.html']
    for entry in json.loads((DOCS/'assets'/catalog_path).read_text())['templates']:
        page.ids.add(entry['id'])
        if entry.get('topology_class'):
            page.ids.add(entry['topology_class'])
for path,page in pages.items():
    for link in page.links:
        u=urlsplit(link)
        if u.scheme or u.netloc:continue
        target=(path.parent/unquote(u.path)).resolve() if u.path else path
        assert target.is_relative_to(DOCS.resolve()), (path,link)
        assert target.is_file(), (path,link)
        if u.fragment and target in pages:assert u.fragment in pages[target].ids,(path,link)
for path in (DOCS/'assets/figures').glob('*.svg'):
    ET.parse(path)
catalog=json.loads((DOCS/'assets/catalog.json').read_text())
with zipfile.ZipFile(DOCS/'assets/templates.zip') as archive:
    assert archive.testzip() is None
    for r in catalog['templates']:
        for key in ('mesh','certificate'):
            p=DOCS/r[key]
            assert p.is_file() and archive.read(r[key])==p.read_bytes()
    assert archive.read('catalog.json')==(DOCS/'assets/catalog.json').read_bytes()
    assert archive.read('quality.csv')==(DOCS/'assets/quality.csv').read_bytes()
for p in DOCS.rglob('*'):
    if p.suffix in ('.html','.js','.json','.csv'):
        assert '/tmp/center-test' not in p.read_text(), p
geode=json.loads((DOCS/'assets/geode/catalog.json').read_text())
assert len(geode['templates'])==29
assert {r['id'] for r in geode['templates'] if r.get('side','') and r['side'].startswith('Q7')}=={
    'G-Q7-06-H42-02', 'G-Q7-06-H44-01', 'G-Q7-06-H44-02', 'G-Q7-06-flipped-H44-02'}
assert {r['id'] for r in geode['templates'] if r.get('side')=='Q6-01'}=={
    'G-Q6-01-H42-05', 'G-Q6-01-H44-09'}
assert [r['hexes'] for r in geode['templates']]==sorted(r['hexes'] for r in geode['templates'])
with zipfile.ZipFile(DOCS/'assets/geode/geodes.zip') as archive:
    assert archive.testzip() is None
    for r in geode['templates']:
        for key in ('mesh','certificate','connectivity','vtu'):
            p=DOCS/r[key]
            assert p.is_file() and archive.read(r[key])==p.read_bytes()
        cert=json.loads((DOCS/r['certificate']).read_text())
        assert cert['accepted'] and cert['all_jacobians_positive'] and cert['translation_matching']
        connectivity=json.loads((DOCS/r['connectivity']).read_text())
        assert connectivity['index_base']==0 and len(connectivity['cells'])==r['hexes']
        assert connectivity['mesh_sha256']==cert['mesh_sha256']==r['sha256']
        assert hashlib.sha256((DOCS/r['mesh']).read_bytes()).hexdigest()==r['sha256']
        assert len(r['metrics']['cell_min_sj'])==len(r['symmetry']['cell_orbits'])==r['hexes']
        for row in connectivity['cells']:
            assert len(set(row))==8 and min(row)>=0 and max(row)<r['vertices']
    assert archive.read('catalog.json')==(DOCS/'assets/geode/catalog.json').read_bytes()
print(f"SITE_OK relative links, fragments, {len(catalog['templates'])} pyramid and {len(geode['templates'])} Geode meshes, archives, and portable data")
