#include "lab_page.hpp"

namespace climate {

const char* climate_lab_js() {
    return R"LABJS(
/* ---------------- COMPLETE FORMULA LAB ---------------- */
const lab={groups:[],count:0,selected:new Set(),filter:''};
function labAll(){const out=[];lab.groups.forEach(g=>g.formulas.forEach(f=>out.push(f)));return out;}
function labLookup(id){for(const g of lab.groups){for(const f of g.formulas)if(f.id===id)return f;}return null;}
async function loadLab(){
  try{
    const d=await api('/api/lab');
    lab.groups=d.groups;lab.count=d.count;
    $('#lab-count').textContent=d.count;
    const all=labAll();
    ['physics::equilibrium_temperature_change','physics::evi','trend::sens_slope','spatial::getis_ord_gi']
      .forEach(id=>{if(all.some(f=>f.id===id))lab.selected.add(id);});
    renderLab();
  }catch(e){$('#lab-list').innerHTML=`<div class="err">${esc(e.message)}</div>`;}
}
function labVisible(){
  const q=lab.filter.trim().toLowerCase();
  return lab.groups.map(g=>{
    const formulas=q?g.formulas.filter(f=>`${f.name} ${f.module} ${f.returns} ${f.signature}`.toLowerCase().includes(q)):g.formulas;
    return {id:g.id,label:g.label,formulas};
  }).filter(g=>g.formulas.length);
}
function renderLab(){
  const groups=labVisible();
  $('#lab-list').innerHTML=groups.map(g=>{
    const all=g.formulas.every(f=>lab.selected.has(f.id));
    return `<div class="lab-group">
      <label class="lab-grouphead"><input type="checkbox" data-lab-group="${esc(g.id)}" ${all?'checked':''}><span>${esc(g.id)} · ${esc(g.label)}</span><span class="hint">${g.formulas.length}</span></label>
      ${g.formulas.map(f=>`<label class="lab-row"><input type="checkbox" data-lab-pick="${esc(f.id)}" ${lab.selected.has(f.id)?'checked':''}><span class="lab-name">${esc(f.name)}</span><span class="lab-meta">${esc(f.module)} · ${esc(f.returns)}</span></label>`).join('')}
    </div>`;
  }).join('')||'<div class="hint">No formula matches that filter.</div>';
  renderLabInputs();
  $('#lab-selected').textContent=`${lab.selected.size} selected of ${lab.count}`;
}
function renderLabInputs(){
  const ids=[...lab.selected];
  if(!ids.length){$('#lab-inputs').innerHTML='<div class="hint">Nothing selected yet.</div>';return;}
  const scalar=['number','complex','seed','enum'];
  $('#lab-inputs').innerHTML=ids.map(id=>{
    const f=labLookup(id);if(!f)return '';
    const nums=f.params.filter(p=>scalar.includes(p.kind));
    const shared=f.params.some(p=>!scalar.includes(p.kind));
    return `<div class="ph-card"><div class="ph-title">${esc(f.name)}<span class="hint" style="margin-left:auto">${esc(f.module)}</span></div>
      <div class="formula">${esc(f.signature)}</div>
      ${nums.length?`<div class="ph-grid">${nums.map(p=>`<div class="field"><label>${esc(p.name)}</label><input data-lab-arg="${esc(id)}|${esc(p.name)}" type="number" step="any" value="${p.default}"></div>`).join('')}</div>`:''}
      ${shared?`<div class="hint">Also reads the shared series, matrix and text inputs.</div>`:''}
    </div>`;
  }).join('');
}
async function runLab(){
  if(!lab.selected.size){toast('Select at least one formula.');return;}
  const argsById={};
  $$('[data-lab-arg]').forEach(input=>{
    const [id,name]=input.dataset.labArg.split('|');
    if(!argsById[id])argsById[id]={};
    const value=Number(input.value);
    if(Number.isFinite(value))argsById[id][name]=value;
  });
  const body={file:$('#lab-file').value,text:$('#lab-text').value,
    matrix_columns:Number($('#lab-cols').value)||3,
    calls:[...lab.selected].map(id=>({id,args:argsById[id]||{}}))};
  const button=$('#lab-run');button.disabled=true;button.textContent='Running…';
  try{
    const d=await api('/api/lab/run',body);
    $('#lab-summary').textContent=`${d.succeeded}/${d.requested} returned a value · ${d.series} · ${d.samples} samples`;
    $('#lab-output').innerHTML=d.results.map(r=>`<div class="lab-hit ${r.ok?'':'bad'}"><div class="ph-title"><b>${esc(r.name||r.id)}</b><span class="hint" style="margin-left:auto">${esc(r.module||'')}</span></div>${r.ok?`<pre>${esc(r.result)}</pre>`:`<pre class="err">${esc(r.error)}</pre>`}</div>`).join('');
  }catch(e){toast('Run failed: '+e.message);}
  finally{button.disabled=false;button.textContent='Run selected →';}
}
$('#lab-search')&&$('#lab-search').addEventListener('input',e=>{lab.filter=e.target.value;renderLab();});
$('#lab-list')&&$('#lab-list').addEventListener('change',e=>{
  const pick=e.target.closest('[data-lab-pick]');
  if(pick){if(pick.checked)lab.selected.add(pick.dataset.labPick);else lab.selected.delete(pick.dataset.labPick);renderLab();return;}
  const group=e.target.closest('[data-lab-group]');
  if(group){const target=labVisible().find(g=>g.id===group.dataset.labGroup);
    if(target)target.formulas.forEach(f=>{if(group.checked)lab.selected.add(f.id);else lab.selected.delete(f.id);});
    renderLab();}
});
$('#lab-select-visible')&&$('#lab-select-visible').addEventListener('click',()=>{labVisible().forEach(g=>g.formulas.forEach(f=>lab.selected.add(f.id)));renderLab();});
$('#lab-clear')&&$('#lab-clear').addEventListener('click',()=>{lab.selected.clear();renderLab();});
$('#lab-run')&&$('#lab-run').addEventListener('click',runLab);
)LABJS";
}

}  // namespace climate
