/* Catalog rendering stays independent of WebGL so data remain usable if the viewer cannot load. */
(() => {
  const METRICS = [
    ['min_scaled_jacobian','Minimum sampled scaled Jacobian ↑',6],['nine_point_min_sj','Corner + center minimum SJ ↑',6],
    ['cell_min_sj_p05','5th percentile of cell minima ↑',6],['cell_min_sj_median','Median of cell minima ↑',6],
    ['max_condition','Worst Jacobian condition number ↓',3],['min_det','Minimum sampled determinant ↑','exp'],
    ['volume_min','Smallest cell volume','exp'],['volume_max','Largest cell volume',6],
    ['volume_ratio','Volume max / min',3],['volume_cv','Volume coefficient of variation',4],
    ['worst_cell_edge_ratio','Worst cell edge-length ratio ↓',3]
  ];
  let catalog, selected='all', sort='quality', viewerModule, observer;
  const gallery=document.getElementById('gallery'), status=document.getElementById('collection-status');
  const esc=s=>String(s).replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
  const number=(n,p)=>p==='exp'?n.toExponential(3):n.toFixed(p);
  const loaded=new Map();
  const requestViewer=()=>viewerModule||(viewerModule=import('./viewer.js').catch(error=>{viewerModule=null;throw error;}));
  function card(t){
    const s=t.symmetry, displayId=t.topology_class||t.id;
    const extraDownloads=['connectivity','vtu','periodic_vtu'].filter(key=>t[key]).map(key=>`<a href="${esc(t[key])}" download>${key==='connectivity'?(t.periodic?'Periodic JSON':'Connectivity (.json)'):key==='periodic_vtu'?'2×2 VTU':'VTU'}</a>`).join('');
    const kicker=t.periodic?`${t.hexes} cells per tile · ${t.vertices} periodic vertices`:`${t.vertices} vertices · ${t.interior_vertices} interior`;
    const description=t.periodic?`<strong>${t.cap_faces[0]} → ${t.cap_faces[1]} cap quads · height ${number(t.height,3)}</strong><br>${esc(t.origin)}`:`<span class="sym-label">Symmetry: ${s.group} (order ${s.order})</span><strong>${s.mirrors.length?s.mirrors.map(esc).join(' & ')+(s.mirrors.length===1?' mirror':' mirrors'):'No measured mirrors'}</strong>${s.rotations.length?' · '+s.rotations.map(esc).join(', ')+' rotation':''}<br>${s.orbit_count} cell orbits · ${esc(t.origin)}`;
    const repetition=t.periodic?`<label>View <select class="repeat-mode" aria-label="Periodic completion for ${esc(t.id)}"><option value="1">One tile</option><option value="2">2×2 tiles</option></select></label>`:'';
    const variants=t.variants?`<details class="geometry-variants"><summary>${t.variants.length} ${t.variants.length===1?'embedding':'embeddings'} · representative ${esc(t.id)}</summary><p>Same periodic connectivity; different stored coordinates or cell placements.</p><ul>${t.variants.map(v=>`<li><a href="${esc(v.connectivity)}" download>${esc(v.id)}</a><br>min SJ ${number(v.min_sj,4)} · condition ${number(v.condition,2)} · height ${number(v.height,3)} · <a href="${esc(v.certificate)}">certificate</a></li>`).join('')}</ul></details>`:'';
    const height=t.height_ratio===undefined?'':`<br><strong>Height / width ${number(t.height_ratio,3)}</strong> · ${esc({fixed:'fixed reference',sj:'best sampled SJ',condition:'lowest condition among trials'}[t.geometry_mode])}`;
    const metrics=METRICS.map(([key,label,precision])=>`<tr><th scope="row">${label}</th><td>${number(t.metrics[key],precision)}</td></tr>`).join('');
    const e=document.createElement('article');e.className='mesh-card';e.id=t.id;e.dataset.meshId=t.id;
    e.innerHTML=`<header class="card-header"><div><p class="card-kicker">${kicker}</p><h4${t.topology_class?' id="'+esc(t.topology_class)+'"':''}>${esc(displayId)}</h4></div>${t.reference_url?'<a class="publication-link" href="'+esc(t.reference_url)+'">Publication</a>':''}</header>
      <div class="viewer" aria-label="Interactive shrunk-element view of ${esc(displayId)}"><canvas aria-label="${esc(displayId)}: drag to rotate, scroll to zoom" tabindex="0"></canvas><div class="viewer-status" role="status">Loading mesh…</div><button class="reset-view" title="Reset view" aria-label="Reset view for ${esc(displayId)}">↺</button></div>
      <div class="viewer-controls"><label class="shrink-control" for="shrink-${t.id}">Shrink <input id="shrink-${t.id}" class="shrink" type="range" min="0" max="55" value="18" step="1"><output for="shrink-${t.id}">18%</output></label><label><span class="hidden">Color by</span><select class="color-mode" aria-label="Color ${esc(displayId)} by"><option value="quality">Cell quality</option><option value="orbits">${t.periodic?'Cell identity':'Cell orbits'}</option><option value="uniform">Uniform</option></select></label><label class="toggle" ${t.periodic?'hidden':''}><input class="show-planes" type="checkbox">Mirrors</label><label class="toggle"><input class="show-boundary" type="checkbox" checked>${t.periodic?'Caps':'Boundary'}</label>${repetition}</div>
      <div class="card-body"><p class="symmetry-line">${description}${height}</p>${variants}<table class="metrics"><caption>Quality measurements (21³ samples per cell)</caption><tbody>${metrics}</tbody></table><div class="legend"><span>Cell minimum SJ</span><i aria-hidden="true"></i><span>0 → 0.4+</span></div><div class="card-downloads"><a href="${esc(t.mesh)}" download>Mesh (.mesh)</a><a href="${esc(t.certificate)}" download>Validation report</a></div></div>`;
    e.querySelector('.card-downloads').insertAdjacentHTML('beforeend',extraDownloads);
    const range=e.querySelector('.shrink');range.addEventListener('input',()=>{e.querySelector('output').value=range.value+'%';});
    e.querySelector('.color-mode').addEventListener('change',event=>{e.querySelector('.legend').hidden=event.target.value!=='quality';});
    e.template=t;e.visible=false;e.view=null;e.generation=0;return e;
  }
  function stop(e){e.visible=false;e.generation++;if(e.view){e.view.dispose();e.view=null;}loaded.delete(e.id);}
  async function start(e){
    e.visible=true;const generation=++e.generation;
    try{
      const {MeshViewer}=await requestViewer();if(!e.visible||e.generation!==generation)return;
      // Browsers impose a small WebGL-context limit. Only nearby cards own contexts.
      if(loaded.size>=6){const first=loaded.values().next().value;stop(first);}
      const view=new MeshViewer(e,e.template);e.view=view;loaded.set(e.id,e);await view.initialize();
    }catch(error){if(e.generation!==generation)return;e.view?.dispose();e.view=null;loaded.delete(e.id);e.querySelector('.viewer-status').textContent='The 3D viewer could not load. The mesh, metrics, and certificate are available below. Reload to retry.';console.warn('Mesh viewer:',error);}
  }
  function render(){
    if(observer)observer.disconnect();gallery.querySelectorAll('.mesh-card').forEach(stop);gallery.replaceChildren();
    const list=catalog.templates.filter(t=>selected==='all'||t.hexes===Number(selected));
    const counts=[...new Set(list.map(t=>t.hexes))].sort((a,b)=>a-b);
    for(const count of counts){
      const templates=list.filter(t=>t.hexes===count);
      templates.sort((a,b)=>sort==='id'?String(a.topology_class||a.id).localeCompare(String(b.topology_class||b.id),undefined,{numeric:true}):sort==='condition'?a.metrics.max_condition-b.metrics.max_condition:b.metrics.min_scaled_jacobian-a.metrics.min_scaled_jacobian);
      const section=document.createElement('section');section.className='group';section.id='cells-'+count;section.setAttribute('aria-labelledby','heading-'+count);
      section.innerHTML=`<div class="group-heading"><h3 id="heading-${count}">${count} <span>elements</span></h3><span>${templates.length} ${catalog.geometry_count?(templates.length===1?'connectivity':'connectivities'):(templates.length===1?'mesh':'meshes')}${catalog.count_notes?.[count]?' · '+esc(catalog.count_notes[count]):catalog.partial_enumeration_counts.includes(count)?' · diagonal enumeration incomplete':!document.body.dataset.catalog&&count>48?' · published template':''}</span><div class="group-navigation"><button data-direction="-1" aria-label="Previous ${count}-cell mesh">←</button><button data-direction="1" aria-label="Next ${count}-cell mesh">→</button></div></div><div class="mesh-grid" tabindex="0" role="region" aria-label="${count}-cell meshes; scroll horizontally"></div>`;
      const grid=section.querySelector('.mesh-grid');for(const t of templates)grid.append(card(t));
      section.querySelectorAll('[data-direction]').forEach(button=>button.addEventListener('click',()=>grid.scrollBy({left:Number(button.dataset.direction)*(grid.querySelector('.mesh-card').getBoundingClientRect().width+20),behavior:matchMedia('(prefers-reduced-motion: reduce)').matches?'auto':'smooth'})));
      if(templates.length===1)section.querySelector('.group-navigation').hidden=true;
      gallery.append(section);
    }
    status.textContent=catalog.geometry_count?`${list.length} distinct periodic ${list.length===1?'connectivity':'connectivities'} · ${catalog.geometry_count} validated embeddings in the full collection`:`${list.length} ${list.length===1?'mesh':'meshes'}`;
    observer=new IntersectionObserver(entries=>{for(const entry of entries){const e=entry.target;if(entry.isIntersecting&&!e.visible)start(e);else if(!entry.isIntersecting&&e.visible)stop(e);}}, {rootMargin:'140px 0px'});
    gallery.querySelectorAll('.mesh-card').forEach(e=>observer.observe(e));
  }
  document.querySelectorAll('[data-count]').forEach(button=>button.addEventListener('click',()=>{
    selected=button.dataset.count;document.querySelectorAll('[data-count]').forEach(b=>{const active=b===button;b.classList.toggle('active',active);b.setAttribute('aria-pressed',String(active));});
    const url=new URL(location.href);if(selected==='all')url.searchParams.delete('cells');else url.searchParams.set('cells',selected);history.replaceState(null,'',url);if(catalog)render();
  }));
  document.getElementById('sort').addEventListener('change',e=>{sort=e.target.value;if(catalog)render();});
  let catalogPath=document.body.dataset.catalog||'assets/catalog.json';
  const geometry=document.getElementById('geometry');
  if(geometry){
    const paths={fixed:catalogPath,sj:document.body.dataset.catalogSj,condition:document.body.dataset.catalogCondition};
    const requested=new URL(location.href).searchParams.get('geometry');
    geometry.value=Object.hasOwn(paths,requested)?requested:'sj';
    catalogPath=paths[geometry.value];
    geometry.addEventListener('change',()=>{const url=new URL(location.href);url.searchParams.set('geometry',geometry.value);location.assign(url);});
  }
  fetch(catalogPath).then(r=>{if(!r.ok)throw new Error(r.status);return r.json();}).then(data=>{
    catalog=data;const count=new URL(location.href).searchParams.get('cells');const button=[...document.querySelectorAll('[data-count]')].find(b=>b.dataset.count===count);if(button)button.click();else render();
    if(location.hash)requestAnimationFrame(()=>document.getElementById(location.hash.slice(1))?.scrollIntoView());
  }).catch(error=>{status.textContent='The catalog could not load. Serve this directory over HTTP, or download the complete mesh archive above.';console.error(error);});
})();
