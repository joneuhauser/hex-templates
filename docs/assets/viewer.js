import * as THREE from 'three';
import {OrbitControls} from 'three/addons/controls/OrbitControls.js';

export const FACES=[[0,3,2,1],[4,5,6,7],[0,1,5,4],[1,2,6,5],[2,3,7,6],[3,0,4,7]];
export const EDGES=[[0,1],[1,2],[2,3],[3,0],[4,5],[5,6],[6,7],[7,4],[0,4],[1,5],[2,6],[3,7]];
const assets=new Map();
export function parseMesh(text){
  const tokens=text.replace(/#[^\n]*/g,'').trim().split(/\s+/);let i=0,points=[],cells=[],quads=[];
  const integer=()=>{const n=Number(tokens[i++]);if(!Number.isInteger(n)||n<0)throw new Error('Invalid mesh integer');return n;};
  while(i<tokens.length){const key=tokens[i++];if(key==='End')break;
    if(key==='MeshVersionFormatted'){integer();continue;}if(key==='Dimension'){if(integer()!==3)throw new Error('Expected a 3D mesh');continue;}
    const width={Vertices:3,Hexahedra:8,Quadrilaterals:4}[key];if(!width)throw new Error('Unsupported mesh section: '+key);
    const count=integer();const rows=[];
    for(let n=0;n<count;n++){const row=[];for(let j=0;j<width;j++){const value=Number(tokens[i++]);if(!Number.isFinite(value))throw new Error('Invalid mesh coordinate');row.push(key==='Vertices'?value:value-1);}integer();rows.push(row);}
    if(key==='Vertices')points=rows;else if(key==='Hexahedra')cells=rows;else quads=rows;
  }
  if(!points.length||!cells.length)throw new Error('Empty mesh');
  for(const c of [...cells,...quads])for(const v of c)if(!Number.isInteger(v)||v<0||v>=points.length)throw new Error('Mesh index out of range');
  if(cells.some(c=>new Set(c).size!==8))throw new Error('Degenerate cell indexing');
  if(!quads.length){const faces=new Map();for(const cell of cells)for(const f of FACES){const q=f.map(i=>cell[i]),key=[...q].sort((a,b)=>a-b).join(',');const entry=faces.get(key);if(entry)entry.count++;else faces.set(key,{q,count:1});}quads=[...faces.values()].filter(f=>f.count===1).map(f=>f.q);}
  return {points,cells,quads};
}
function loadMesh(url){if(!assets.has(url))assets.set(url,fetch(url).then(r=>{if(!r.ok)throw new Error('Mesh HTTP '+r.status);return r.text();}).then(parseMesh).catch(e=>{assets.delete(url);throw e;}));return assets.get(url);}
const palette=[new THREE.Color('#cc623d'),new THREE.Color('#e3b581'),new THREE.Color('#9dbabc'),new THREE.Color('#3e728b')];
function qualityColor(value){const t=THREE.MathUtils.clamp(value/.4,0,1)*3;const a=Math.min(2,Math.floor(t));return palette[a].clone().lerp(palette[a+1],t-a);}
function facePoint(p,u,v){return p[0].clone().multiplyScalar((1-u)*(1-v)).addScaledVector(p[1],u*(1-v)).addScaledVector(p[2],u*v).addScaledVector(p[3],(1-u)*v);}
function faceNormal(p,u,v){const a=p[1].clone().sub(p[0]).multiplyScalar(1-v).addScaledVector(p[2].clone().sub(p[3]),v);const b=p[3].clone().sub(p[0]).multiplyScalar(1-u).addScaledVector(p[2].clone().sub(p[1]),u);return a.cross(b).normalize();}

export class MeshViewer{
 constructor(card,template){this.card=card;this.template=template;this.disposed=false;this.repetition=1;this.repeatGeneration=0;this.listeners=[];this.dirty=true;}
 listen(element,type,fn){element.addEventListener(type,fn);this.listeners.push(()=>element.removeEventListener(type,fn));}
 async initialize(){
  this.data=await loadMesh(this.template.mesh);if(this.disposed)return;
  if(this.template.display_matrix){const m=this.template.display_matrix;this.data={...this.data,points:this.data.points.map(p=>m.map(row=>row.reduce((sum,a,i)=>sum+a*p[i],0)))};}
  const canvas=this.card.querySelector('canvas');this.renderer=new THREE.WebGLRenderer({canvas,antialias:true,alpha:true,powerPreference:'low-power'});
  this.renderer.setPixelRatio(Math.min(window.devicePixelRatio||1,2));this.renderer.setClearColor(0x000000,0);this.renderer.outputColorSpace=THREE.SRGBColorSpace;this.renderer.toneMapping=THREE.ACESFilmicToneMapping;this.renderer.toneMappingExposure=1.15;
  this.scene=new THREE.Scene();this.scene.add(new THREE.HemisphereLight(0xffffff,0x738b78,2.1));
  const key=new THREE.DirectionalLight(0xfff6e9,3.1);key.position.set(-3,-5,7);this.scene.add(key);const fill=new THREE.DirectionalLight(0xe0efff,1.3);fill.position.set(4,2,3);this.scene.add(fill);
  this.camera=new THREE.PerspectiveCamera(33,1,.05,80);this.camera.up.set(0,0,1);
  this.controls=new OrbitControls(this.camera,canvas);this.controls.enableDamping=false;this.controls.minDistance=2.3;this.controls.maxDistance=15;this.controls.enablePan=true;this.controls.target.set(0,0,.59);this.reset();
  this.controls.addEventListener('change',()=>{this.dirty=true;});
  this.material=new THREE.MeshStandardMaterial({vertexColors:true,roughness:.62,metalness:.02,side:THREE.DoubleSide});this.edgeMaterial=new THREE.LineBasicMaterial({color:0x294450,transparent:true,opacity:.6});
  this.mesh=new THREE.Mesh(new THREE.BufferGeometry(),this.material);this.edges=new THREE.LineSegments(new THREE.BufferGeometry(),this.edgeMaterial);this.scene.add(this.mesh,this.edges);this.makeBoundary();this.makePlanes();this.rebuild();
  this.listen(this.card.querySelector('.shrink'),'input',()=>this.rebuild());this.listen(this.card.querySelector('.color-mode'),'change',()=>this.rebuild());this.listen(this.card.querySelector('.reset-view'),'click',()=>this.reset());
  this.listen(this.card.querySelector('.show-boundary'),'change',e=>{this.boundary.visible=e.target.checked;this.dirty=true;});this.listen(this.card.querySelector('.show-planes'),'change',e=>{this.planes.visible=e.target.checked;this.dirty=true;});
  const repeat=this.card.querySelector('.repeat-mode');
  if(repeat)this.listen(repeat,'change',()=>this.setRepetition(Number(repeat.value)));
  this.listen(canvas,'keydown',e=>{if(e.key.toLowerCase()==='r'){this.reset();e.preventDefault();}});
  this.listen(canvas,'webglcontextlost',()=>{if(!this.disposed)this.card.querySelector('.viewer-status').textContent='3D context lost. Reload this page to restore the view.';});
  this.resize=new ResizeObserver(()=>{if(this.disposed)return;const box=this.card.querySelector('.viewer').getBoundingClientRect();this.renderer.setSize(box.width,box.height,false);this.camera.aspect=box.width/box.height;this.camera.updateProjectionMatrix();this.dirty=true;});this.resize.observe(this.card.querySelector('.viewer'));
  this.card.querySelector('.viewer-status').textContent='';this.card.querySelector('.viewer').dataset.ready='true';
  if(repeat&&Number(repeat.value)!==1)await this.setRepetition(Number(repeat.value));
  if(this.disposed)return;
  this.renderer.setAnimationLoop(()=>{if(this.disposed||!this.dirty)return;this.renderer.render(this.scene,this.camera);this.dirty=false;});
 }
 reset(){if(!this.camera)return;this.camera.position.set(...(this.repetition===2?this.template.repeated_camera_position:this.template.camera_position||[3.8,-5.3,3.35]));this.controls.target.set(...(this.repetition===2?this.template.repeated_camera_target:this.template.camera_target||[0,0,.59]));this.controls.update();this.dirty=true;}
 async setRepetition(count){
  const generation=++this.repeatGeneration,status=this.card.querySelector('.viewer-status');
  status.textContent='Loading periodic completion…';
  try{
   const data=await loadMesh(count===2?this.template.repeated_mesh:this.template.mesh);
   if(this.disposed||generation!==this.repeatGeneration)return;
   this.data=data;this.repetition=count;
   this.scene.remove(this.boundary);this.boundary.geometry.dispose();this.boundary.material.dispose();
   this.makeBoundary();this.rebuild();this.reset();status.textContent='';
  }catch(error){
   if(this.disposed||generation!==this.repeatGeneration)return;
   this.card.querySelector('.repeat-mode').value=String(this.repetition);
   status.textContent='Could not load the completion. Try switching views again.';
  }
 }
 rebuild(){
  if(this.disposed||!this.mesh)return;const shrink=1-Number(this.card.querySelector('.shrink').value)/100,mode=this.card.querySelector('.color-mode').value;
  const positions=[],normals=[],colors=[],lines=[];const divisions=3;
  this.data.cells.forEach((cell,id)=>{
   const vertices=cell.map(v=>new THREE.Vector3(...this.data.points[v]));const center=vertices.reduce((a,b)=>a.add(b),new THREE.Vector3()).multiplyScalar(1/8);vertices.forEach(v=>v.sub(center).multiplyScalar(shrink).add(center));
   const color=mode==='uniform'?new THREE.Color('#4b7f96'):mode==='orbits'?new THREE.Color().setHSL((this.template.symmetry.cell_orbits[id%this.template.hexes]*.61803398875+.53)%1,.34,.52):qualityColor(this.template.metrics.cell_min_sj[id%this.template.hexes]);
   const push=(p,u,v)=>{const point=facePoint(p,u,v),normal=faceNormal(p,u,v);positions.push(point.x,point.y,point.z);normals.push(normal.x,normal.y,normal.z);colors.push(color.r,color.g,color.b);};
   for(const f of FACES){const p=f.map(i=>vertices[i]);for(let u=0;u<divisions;u++)for(let v=0;v<divisions;v++){
    const a=u/divisions,b=v/divisions,c=(u+1)/divisions,d=(v+1)/divisions;push(p,a,b);push(p,c,b);push(p,c,d);push(p,a,b);push(p,c,d);push(p,a,d);
   }}
   for(const [a,b] of EDGES)lines.push(...vertices[a].toArray(),...vertices[b].toArray());
  });
  const geometry=new THREE.BufferGeometry();geometry.setAttribute('position',new THREE.Float32BufferAttribute(positions,3));geometry.setAttribute('normal',new THREE.Float32BufferAttribute(normals,3));geometry.setAttribute('color',new THREE.Float32BufferAttribute(colors,3));
  const edgeGeometry=new THREE.BufferGeometry();edgeGeometry.setAttribute('position',new THREE.Float32BufferAttribute(lines,3));this.mesh.geometry.dispose();this.edges.geometry.dispose();this.mesh.geometry=geometry;this.edges.geometry=edgeGeometry;this.dirty=true;
 }
 makeBoundary(){const seen=new Set(),lines=[];for(const q of this.data.quads)for(let i=0;i<4;i++){const a=q[i],b=q[(i+1)%4],key=[a,b].sort((a,b)=>a-b).join(',');if(seen.has(key))continue;seen.add(key);lines.push(...this.data.points[a],...this.data.points[b]);}const g=new THREE.BufferGeometry();g.setAttribute('position',new THREE.Float32BufferAttribute(lines,3));this.boundary=new THREE.LineSegments(g,new THREE.LineBasicMaterial({color:0x536c64,transparent:true,opacity:.24}));this.boundary.visible=this.card.querySelector('.show-boundary').checked;this.scene.add(this.boundary);}
 makePlanes(){this.planes=new THREE.Group();for(const action of this.template.symmetry.actions.filter(a=>a.kind==='mirror')){
   const d=action.id==='mx'?[0,1]:action.id==='my'?[1,0]:action.id==='md'?[1,1]:[1,-1];const length=Math.hypot(...d),a=d.map(x=>x/length*(this.template.plane_radius||1.25)),z=this.template.plane_range||[-.03,1.64];const points=[[-a[0],-a[1],z[0]],[a[0],a[1],z[0]],[a[0],a[1],z[1]],[-a[0],-a[1],z[1]]];
   const g=new THREE.BufferGeometry();g.setAttribute('position',new THREE.Float32BufferAttribute([0,1,2,0,2,3].flatMap(i=>points[i]),3));g.computeVertexNormals();const m=new THREE.MeshBasicMaterial({color:0xb3633e,side:THREE.DoubleSide,transparent:true,opacity:.12,depthWrite:false});this.planes.add(new THREE.Mesh(g,m));
  }this.planes.visible=this.card.querySelector('.show-planes').checked;this.scene.add(this.planes);}
 dispose(){this.disposed=true;this.listeners.forEach(remove=>remove());this.resize?.disconnect();this.controls?.dispose();if(this.renderer){this.renderer.setAnimationLoop(null);this.scene.traverse(o=>{o.geometry?.dispose();if(o.material)for(const m of Array.isArray(o.material)?o.material:[o.material])m.dispose();});this.renderer.dispose();this.renderer.forceContextLoss();const canvas=this.card.querySelector('canvas');canvas.replaceWith(canvas.cloneNode(false));}this.card.querySelector('.viewer').removeAttribute('data-ready');}
}
