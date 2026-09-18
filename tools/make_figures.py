"""Deterministic vector figures for the atlas and its mathematical methodology."""
from pathlib import Path
import json, math, sys
from html import escape
import numpy as np
from mesh_tools import read_mesh, FACES
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'docs/assets/figures';OUT.mkdir(parents=True,exist_ok=True)
INK='#233c48';BLUE='#467a91';ORANGE='#b65a32';MUTED='#6d7f83';LINE='#cfdbd1'

def svg(name,width,height,body,title):
 (OUT/name).write_text(f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" role="img" aria-label="{escape(title, quote=True)}"><title>{escape(title)}</title><style>text{{font-family:system-ui,sans-serif;fill:{INK};font-size:13px}}.small{{font-size:11px;fill:{MUTED}}}.label{{font-size:11px;letter-spacing:1px;fill:{ORANGE}}}</style>'+body+'</svg>')
def text(x,y,t,cls='',anchor='middle'):return f'<text x="{x}" y="{y}" class="{cls}" text-anchor="{anchor}">{escape(str(t))}</text>'
def line(x,y,a,b,color=LINE,width=1,dash=''):return f'<line x1="{x}" y1="{y}" x2="{a}" y2="{b}" stroke="{color}" stroke-width="{width}"'+(f' stroke-dasharray="{dash}"' if dash else '')+'/>'
def rect(x,y,w,h,fill='#eef3ec',stroke=LINE):return f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="3" fill="{fill}" stroke="{stroke}"/>'
def circle(x,y,r=7,fill=BLUE):return f'<circle cx="{x}" cy="{y}" r="{r}" fill="{fill}"/>'
svg('mark.svg',40,44,'<g stroke="#245a82" stroke-width="2" stroke-linejoin="round"><path d="M20 3 37 13 20 23 3 13Z" fill="#eaf2f6"/><path d="M3 13 20 23 20 41 3 31Z" fill="#cfdfe8"/><path d="M20 23 37 13 37 31 20 41Z" fill="#93b5c8"/></g>','Hexahedral cell mark')
p,_,q=read_mesh(ROOT/'solver/input/pyramid.mesh')
camera=np.array([4.,-6.,3.5]);direction=camera/np.linalg.norm(camera);right=np.cross(direction,[0,0,1]);right/=np.linalg.norm(right);up=np.cross(right,direction)
def project(v):return np.array([282+v@right*132,254-v@up*132])
xy=np.array([project(v) for v in p]);body=''
for k in range(-3,4):body+=line(100+k*34,305,330+k*34,175,'#dbe4dc',.6)
for k in range(-3,4):body+=line(110+k*34,175,360+k*34,305,'#dbe4dc',.6)
for face in sorted(q,key=lambda f:float(np.mean(p[f],axis=0)@direction)):
 pts=' '.join(f'{x:.2f},{y:.2f}' for x,y in xy[face]);normal=np.cross(p[face[1]]-p[face[0]],p[face[3]]-p[face[0]])
 fill='#d8e6db' if normal@direction>0 else '#edf2e9'
 body+=f'<polygon points="{pts}" fill="{fill}" fill-opacity=".93" stroke="#7c9994" stroke-width="1.1"/>'
for v in sorted(set(q.flat)):
 x,y=xy[v];body+=circle(round(x,2),round(y,2),2.4,ORANGE)
body+=line(330,113,429,78,MUTED,.8)+text(435,73,'3 quads per side','small','start')
body+=rect(42,60,52,52,'#edf2e9','#77958c')+line(68,60,68,112,'#77958c')+line(42,86,94,86,'#77958c')+text(68,133,'4 base quads','small')
body+=text(292,326,'Fixed boundary · free interior','small')
svg('pyramid.svg',580,350,body,'Quadrangulated pyramid and four-cell base inset')

body=''
labels=[('01','Enumerate','Cell connectivities'),('02','Embed & improve','Fixed boundary and mirrors'),('03','Certify','Topology + whole-cell det J'),('04','Publish','Validated meshes')]
for i,(n,title,sub) in enumerate(labels):
 x=10+i*190;body+=rect(x,35,172,115)+text(x+16,57,n,'label','start')+text(x+86,85,title)+text(x+86,109,sub,'small')
 if i<3:body+=text(x+181,95,'→')
body+=text(385,184,'Only certified meshes enter the collection.','small')
svg('pipeline.svg',780,205,body,'Search, geometry optimization, exact certification, and publication pipeline')

body=''
for panel in range(3):
 x=panel*250;cx=x+125;cy=120
 body+=rect(x+10,18,230,218,'#fafcf7')+line(cx,47,cx,186,ORANGE,1.4,'5 4')
 if panel==0:body+=rect(cx-57,82,114,68,'#c4d9dd',BLUE)
 elif panel==1:
  for a in [-1,1]:body+=rect(cx+(13 if a==1 else -78),82,65,68,'#c4d9dd',BLUE)
 else:
  body+=line(x+44,cy,x+207,cy,ORANGE,1.4,'5 4')
  for a in [-1,1]:
   for b in [-1,1]:body+=rect(cx+(10 if a==1 else -66),cy+(10 if b==1 else -55),56,45,'#c4d9dd',BLUE)
 body+=text(cx,211,['One fixed cell','Two exchanged cells','Four exchanged cells'][panel])+text(cx,259,['|orbit| = 1','|orbit| = 2','|orbit| = 4'][panel],'small')
svg('orbits.svg',750,278,body,'Top views of fixed and exchanged cell orbits; dashed lines are mirror planes')

body=''
for offset,invalid in [(20,False),(400,True)]:
 pts=[(offset+70,60),(offset+210,60),(offset+210,175),(offset+70,175)]
 for i in range(4):body+=line(*pts[i],*pts[(i+1)%4],BLUE,2)
 if invalid:body+=line(*pts[0],*pts[2],ORANGE,3)+text(offset+170,105,'same color','small')
 for i,(x,y) in enumerate(pts):body+=circle(x,y,15,BLUE if i%2==0 else '#d59465')+f'<text x="{x}" y="{y+4}" text-anchor="middle" style="fill:white">{i%2}</text>'
 body+=text(offset+140,222,'Consistent quad cycle' if not invalid else 'An odd cycle cannot be repaired')
svg('parity.svg',760,248,body,'Opposite-color edge constraints reject an odd cycle as soon as it forms')

body='';center=(198,143);pts=[(center[0]+94*math.cos(-math.pi/2+i*2*math.pi/5),center[1]+85*math.sin(-math.pi/2+i*2*math.pi/5)) for i in range(5)]
for i in range(5):body+=line(*pts[i],*pts[(i+1)%5],ORANGE,2)
for i,(x,y) in enumerate(pts):body+=circle(x,y,17,BLUE)+f'<text x="{x}" y="{y+4}" text-anchor="middle" style="fill:white">q{i+1}</text>'
body+=text(198,28,'Incompatibility graph')+text(198,264,'Packing finds 2; a full coloring needs 3.','small')+text(380,145,'→')
for i,(label,faces) in enumerate([('Cell A','q1, q3'),('Cell B','q2, q4'),('Cell C','q5')]):
 x=421+i*103;body+=rect(x,91,90,104,'#e2ece7')+text(x+45,127,label)+text(x+45,160,faces,'small')
body+=text(572,237,'Each group must also fit a cube.','small')
svg('cover.svg',750,288,body,'Illustrative five-cycle incompatibility graph and a three-cell cover')

body='';pts=np.array([[85,185],[190,185],[240,135],[135,135],[85,80],[190,80],[240,30],[135,30]])
for a,b in [[0,1],[1,2],[2,3],[3,0],[4,5],[5,6],[6,7],[7,4],[0,4],[1,5],[2,6],[3,7]]:
 body+=line(*pts[a],*pts[b],BLUE if max(a,b)<6 else '#c89678',1.7,'' if max(a,b)<6 else '5 4')
for i,(x,y) in enumerate(pts):body+=circle(x,y,6,BLUE if i<6 else '#fffaf0')+(f'<circle cx="{x}" cy="{y}" r="6" fill="none" stroke="{ORANGE}"/>' if i>=6 else '')+text(x+13,y-10,str(i),'small')
body+=text(160,232,'6 assigned corners','small')+text(324,119,'→')
body+=rect(382,43,326,143)+text(545,74,'An optimistic successor')+text(545,103,'Remove every face that might be consumed.','small')+text(545,129,'Keep only unavoidable new faces.','small')+text(545,155,'Reflect both the faces and the constraints.','small')
body+=text(545,231,'Reject only when even this relaxation cannot fit.','small')
svg('partial.svg',750,262,body,'Partial cube assignment and the conservative successor used by the cover bound')

# Standard plotting tools produce the two scientific figures as standalone SVGs.
import matplotlib
matplotlib.use('Agg');import matplotlib.pyplot as plt
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'axes.edgecolor':LINE,'text.color':INK,'axes.labelcolor':INK,'xtick.color':MUTED,'ytick.color':MUTED,'svg.fonttype':'none','figure.facecolor':'#fcfdf9','axes.facecolor':'#fcfdf9'})
fig,ax=plt.subplots(figsize=(8.2,3.25));t=np.linspace(0,1,400);ax.plot(t,(t-.5)**2+.02,color=BLUE,lw=2.3,label='Positive polynomial');ax.plot([0,.5,1],[.27,-.23,.27],':',color=ORANGE,label='Original Bernstein coefficients');ax.scatter([0,.5,1],[.27,-.23,.27],color=ORANGE,s=28)
ax.plot([0,.25,.5,.75,1],[.27,.02,.02,.02,.27],'--',color='#6f9977',label='Coefficients after subdivision');ax.axhline(0,color=LINE,lw=1);ax.axvline(.5,color=LINE,lw=1);ax.set_xlabel('Reference coordinate t');ax.set_ylabel('Value');ax.legend(frameon=False,fontsize=8,loc='lower right');fig.tight_layout();fig.savefig(OUT/'bernstein.svg');plt.close(fig)
catalog=[r for r in json.loads((ROOT/'docs/assets/catalog.json').read_text())['templates'] if r['hexes']<=48];fig,ax=plt.subplots(figsize=(8.2,3.7))
for count in sorted({r['hexes'] for r in catalog}):
 rows=[r for r in catalog if r['hexes']==count];x=np.zeros(len(rows));y=[r['metrics']['min_scaled_jacobian'] for r in rows];ax.scatter(np.array(x)+count,y,s=23,color=BLUE,alpha=.7);best=int(np.argmax(y));ax.scatter([count+x[best]],[y[best]],s=60,facecolors='none',edgecolors=ORANGE,linewidths=1.6)
ax.set_xticks([36,40,42,44,46,48]);ax.set_xlabel('Hexahedra');ax.set_ylabel('Sampled minimum scaled Jacobian');ax.set_ylim(0,.405);ax.grid(axis='y',alpha=.22);fig.tight_layout();fig.savefig(OUT/'quality.svg');plt.close(fig)
for path in OUT.glob('*.svg'):
 path.write_text('\n'.join(line.rstrip() for line in path.read_text().splitlines())+'\n')
print('Generated 8 SVG figures and the atlas mark')

from make_search_figure import generate as generate_search_figure
generate_search_figure()
