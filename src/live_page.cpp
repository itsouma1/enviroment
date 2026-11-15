#include "live_page.hpp"

namespace climate {

const char* climate_live_js() {
    return R"LIVEJS(
/* ---------------- LIVE NASA DATA ---------------- */
const LIVE={sources:[]};
async function loadLive(){
  try{
    const d=await api('/api/live/status');
    LIVE.sources=d.sources;
    $('#live-creds').textContent=d.earthaccess_credentials
      ? 'Earthdata credentials detected — Earthdata Library pulls are available'
      : 'No Earthdata credentials found — open endpoints work, Earthdata pulls will say so';
    renderLive();
  }catch(e){
    $('#live-list').innerHTML=`<div class="err">${esc(e.message)}</div>`;
  }
}
function renderLive(){
  $('#live-list').innerHTML=LIVE.sources.map(s=>{
    const c=s.cached;
    return `<div class="panel livecard">
      <div class="panelhead"><b>${esc(s.label)}</b><span class="hint">${esc(s.provider)}</span></div>
      <div class="panelbody">
        <div class="kv"><span class="k">endpoint</span><span class="v" style="word-break:break-all">${esc(s.url)}</span></div>
        <div class="kv"><span class="k">cadence</span><span class="v">${esc(s.cadence)}</span></div>
        <div class="kv"><span class="k">access</span><span class="v">${s.needs_auth?'Earthdata Login':'open'}</span></div>
        <div class="kv"><span class="k">cached</span><span class="v">${c?`${c.samples} samples · ${esc(c.retrieved_at)}`:'nothing yet'}</span></div>
        <button class="btn primary" data-live-pull="${esc(s.key)}">Pull from ${esc(s.key)} →</button>
      </div>
    </div>`;
  }).join('');
}
function liveBody(source){
  const body={source};
  if(source==='power'){
    body.lat=Number($('#live-lat').value);
    body.lon=Number($('#live-lon').value);
    body.start=$('#live-start').value;
    body.end=$('#live-end').value;
    body.parameter=$('#live-param').value;
  }
  if(source==='seaice') body.region=$('#live-region').value;
  if(source==='earthaccess'){
    body.short_name=$('#live-short').value;
    body.start=$('#live-start').value;
    body.end=$('#live-end').value;
    const user=$('#ea-user'),pass=$('#ea-pass'),token=$('#ea-token');
    if(user&&user.value.trim())body.earthdata_username=user.value.trim();
    if(pass&&pass.value)body.earthdata_password=pass.value;
    if(token&&token.value.trim())body.earthdata_token=token.value.trim();
  }
  return body;
}
async function pullLive(source){
  const out=$('#live-output');
  const started=Date.now();
  out.innerHTML=`<div class="panel"><div class="panelhead"><b>Contacting ${esc(source)}</b><span class="hint" id="live-tick">0s</span></div><div class="panelbody"><div class="hint">A real network request is in flight. Open endpoints answer in a few seconds; Earthdata downloads a granule and can take up to a minute.</div></div></div>`;
  const tick=setInterval(()=>{const el=document.getElementById('live-tick');if(el)el.textContent=Math.round((Date.now()-started)/1000)+'s';},1000);
  const controller=new AbortController();
  const timer=setTimeout(()=>controller.abort(),180000);
  try{
    const response=await fetch('/api/live/pull',{method:'POST',headers:{'Content-Type':'application/json'},
      body:JSON.stringify(liveBody(source)),signal:controller.signal});
    const raw=await response.text();
    let d;
    try{d=JSON.parse(raw);}catch(err){throw new Error(raw.slice(0,300)||('the server answered HTTP '+response.status));}
    if(!d.ok){
      out.innerHTML=`<div class="panel"><div class="panelhead"><b>Pull failed</b><span class="hint">${esc(source)}</span></div><div class="panelbody"><div class="err">${esc(d.error)}</div>
        ${d.hint?`<div class="note">${esc(d.hint)}</div>`:''}</div></div>`;
      return;
    }
    const s=d.series;
    out.innerHTML=`<div class="panel"><div class="panelhead"><b>${esc(s.variable)} pulled</b><span class="hint">${esc(s.retrieved_at)}</span></div><div class="panelbody">
      <div class="kv"><span class="k">source url</span><span class="v" style="word-break:break-all">${esc(s.url)}</span></div>
      <div class="kv"><span class="k">samples</span><span class="v">${s.samples} · ${esc(s.first_time)} → ${esc(s.last_time)}</span></div>
      <div class="kv"><span class="k">units</span><span class="v">${esc(s.units)}</span></div>
      <div class="kv"><span class="k">file</span><span class="v">${esc(s.file)}</span></div>
      <div class="kv"><span class="k">elapsed</span><span class="v">${((Date.now()-started)/1000).toFixed(1)} s</span></div>
      <div id="live-chart" style="margin-top:8px">loading series…</div>
    </div></div>`;
    const text=await (await fetch('/api/data/file?name='+encodeURIComponent(s.file.split('/').pop()))).text();
    const values=text.split(/\r?\n/).slice(1).map(line=>{const parts=line.split(',');const v=Number(parts[1]);return parts.length>1&&isFinite(v)?v:null;}).filter(v=>v!==null);
    $('#live-chart').innerHTML=chart([{name:s.variable,values,color:'#ffab54'}],{height:220});
    loadLive();
    toast(`${s.samples} real samples from ${source}`);
  }catch(e){
    const message=e.name==='AbortError'
      ? 'The request timed out after 180 s. The endpoint may be unreachable from this machine; try again or pick an open source.'
      : e.message;
    out.innerHTML=`<div class="panel"><div class="panelhead"><b>Pull failed</b><span class="hint">${esc(source)}</span></div><div class="panelbody"><div class="err">${esc(message)}</div></div></div>`;
  }finally{
    clearInterval(tick);clearTimeout(timer);
  }
}
async function verifyEarthdata(){
  const box=$('#ea-status');
  if(box)box.innerHTML='<div class="hint">Asking Earthdata Login to authenticate…</div>';
  try{
    const d=await api('/api/earthdata/check',{
      earthdata_username:($('#ea-user')||{}).value||'',
      earthdata_password:($('#ea-pass')||{}).value||'',
      earthdata_token:($('#ea-token')||{}).value||''});
    if(!box)return;
    box.innerHTML=d.authenticated
      ? `<div class="statusline"><span>LOGIN <b>ACCEPTED</b></span><span>METHOD <b>${esc(d.method||'—')}</b></span><span>USER <b>${esc(d.username||'—')}</b></span></div>
         <div class="note">Earthdata pulls can run now: the credentials stay in this tab and are sent once per pull.</div>`
      : `<div class="err">${esc(d.error||'login failed')}</div>${d.hint?`<div class="note">${esc(d.hint)}</div>`:''}`;
  }catch(e){
    if(box)box.innerHTML=`<div class="err">${esc(e.message)}</div>`;
  }
}
$('#ea-check')&&$('#ea-check').addEventListener('click',verifyEarthdata);
$('#live-reload')&&$('#live-reload').addEventListener('click',loadLive);
$('#live-list')&&$('#live-list').addEventListener('click',e=>{
  const b=e.target.closest('[data-live-pull]');
  if(b)pullLive(b.dataset.livePull);
});
)LIVEJS";
}

}  // namespace climate
