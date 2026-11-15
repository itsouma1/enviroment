#include "visual_page.hpp"

namespace climate {

const char* climate_visual_js() {
    return R"VISJS(
/* ---------------- VISUALIZE EVERYTHING ---------------- */
const SERIES_COLORS=['#f6e6c6','#ffab54','#45d6c9','#ff8f86','#b4e35f','#c6a8ff'];
const STAGE_COLORS=['#45d6c9','#86bcff','#ffc861','#ff8f86','#b4e35f','#c6a8ff','#f6e6c6'];
function labTitle(file){return file.replace(/\.csv$/,'').replace(/_monthly$/,'').replace(/_/g,' ');}
async function loadVisual(){
  try{
    const d=await api('/api/visualize');
    $('#vis-count').textContent=`${d.count} series`;
    $('#vis-grid').innerHTML=d.datasets.map((x,i)=>`
      <div class="panel viscard">
        <div class="panelhead"><b>${esc(labTitle(x.name))}</b><span class="hint">${esc(x.variable)} · ${x.samples} pts</span></div>
        <div class="panelbody">
          ${chart([{name:x.variable,values:x.values,color:SERIES_COLORS[i%SERIES_COLORS.length]}],{height:168})}
          <div class="kv"><span class="k">mean</span><span class="v">${fmt(x.mean,4)}</span></div>
          <div class="kv"><span class="k">range</span><span class="v">${fmt(x.minimum,3)} … ${fmt(x.maximum,3)}</span></div>
          <div class="kv"><span class="k">Sen slope</span><span class="v">${fmt(x.slope,5)}</span></div>
          <div class="kv"><span class="k">end-to-end</span><span class="v" style="color:${x.percent>=0?'#ffab54':'#45d6c9'}">${x.percent>=0?'+':''}${fmt(x.percent,2)}%</span></div>
        </div>
      </div>`).join('');
  }catch(e){toast('Could not load the series views: '+e.message);}
  try{
    const g=await api('/api/grid-trend');
    const slopes=g.slopes;
    const lo=Math.min(...slopes),hi=Math.max(...slopes);
    const span=Math.max(1e-9,Math.max(Math.abs(lo),Math.abs(hi)));
    const cells=slopes.map((s,i)=>{
      const hot=g.hotspots[i];
      const t=Math.min(1,Math.abs(s)/span);
      const hue=s>=0?28:176;
      const light=18+30*t;
      const sat=s>=0?88:62;
      return `<div class="gridcell" title="lat ${fmt(g.latitudes[Math.floor(i/g.columns)],1)} · lon ${fmt(g.longitudes[i%g.columns],1)} · Sen ${fmt(s,5)}${hot?' · Gi* hotspot':''}" style="background:hsl(${hue} ${sat}% ${light}%);${hot?'outline:1px solid #f6e6c6;outline-offset:-1px':''}"></div>`;
    }).join('');
    $('#vis-map').innerHTML=`<div class="gridmap">${cells}</div>
      <div class="tagrow" style="margin-top:7px"><span class="tag">COOLING / TEAL</span><span class="tag amber">WARMING / ORANGE</span><span class="tag">${g.hotspot_count} Gi* HOTSPOTS</span><span class="tag">${g.rows} × ${g.columns} CELLS</span></div>`;
  }catch(e){toast('Could not load the spatial field: '+e.message);}
}
$('#vis-reload')&&$('#vis-reload').addEventListener('click',loadVisual);
)VISJS";
}

}  // namespace climate
