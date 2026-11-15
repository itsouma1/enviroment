#include "ui_page.hpp"

#include <string>

#include "lab_page.hpp"
#include "live_page.hpp"
#include "visual_page.hpp"

namespace climate {

const char* climate_ui_html() {
    static const std::string page = [] {
        std::string html = R"CLIMATE_UI(<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Climate Workbench</title>
<style>
:root{
  --bg:#0b0806; --bg2:#140e09; --panel:#1b130d; --panel2:#241a12;
  --line:#4a3524; --line2:#6b4a30; --text:#f4e8d2; --muted:#cdb795;
  /* Hue plan: warm pole 28 deg (orange) against its complement 174 deg (teal);
     a triad of 28 / 140 / 258 keeps series separable; 347 deg rose is the
     accent complement of 167 deg. Neutrals stay warm so the shell reads as one
     material and every hue above carries the data. */
  --teal:#ff9f43; --sky:#e8c48a; --amber:#ffc35e; --rose:#ff8f86;
  --green:#7fe0a0; --violet:#b49cff; --bezel:#c8b294;
  --series-obs:#f7e8c8; --series-trend:#ff9e4a; --series-season:#3fd0c0; --series-rem:#ff7d92;
  --series-alt:#7fe0a0; --series-violet:#b49cff;
  --mono:"Courier New",Courier,"Liberation Mono",Menlo,monospace;
  --sans:"Courier New",Courier,"Liberation Mono",Menlo,monospace;
}
*{box-sizing:border-box;border-radius:0!important}
html,body{margin:0;padding:0}
html{overflow-y:scroll}
body{min-height:100%}
.wrap{overflow-x:clip}
body{background:radial-gradient(1200px 700px at 50% -10%,#241708,#0d0906 70%);color:var(--text);font-family:var(--mono);font-size:13px;line-height:1.5;
  padding:18px;text-shadow:0 0 8px rgba(255,168,80,.10);-webkit-font-smoothing:antialiased}
body:before{content:"";position:fixed;inset:0;pointer-events:none;z-index:9998;
  box-shadow:inset 0 0 150px rgba(0,0,0,.72),inset 0 0 40px rgba(255,160,70,.05)}
body:after{content:"";position:fixed;inset:0;pointer-events:none;z-index:9999;opacity:.09;
  background-image:repeating-linear-gradient(to bottom,rgba(0,0,0,.55) 0 1px,transparent 1px 3px)}
@media(prefers-reduced-motion:no-preference){
  .app{animation:phosphor 7s ease-in-out infinite}
  @keyframes phosphor{0%,100%{filter:brightness(1)}50%{filter:brightness(1.015)}}}
button,input,select,textarea{font:inherit;color:inherit}
button{cursor:pointer}
a{color:var(--amber)}
::selection{background:#6b4a30;color:#fff3dc}
.app{display:grid;grid-template-columns:228px minmax(0,1fr);min-height:calc(100vh - 36px);
  border:10px solid var(--bezel);outline:2px solid #54402e;
  box-shadow:inset 0 0 90px rgba(0,0,0,.78),0 0 0 1px #000,0 18px 60px rgba(0,0,0,.6);background:var(--bg)}
.sidebar{position:sticky;top:0;height:calc(100vh - 56px);max-height:calc(100vh - 36px);overflow-y:auto;align-self:start;display:flex;flex-direction:column;
  background:#1a120c;border-right:2px solid var(--line);padding:16px 12px 10px;gap:3px;z-index:5}
.brand{display:flex;align-items:center;gap:9px;padding:0 5px 16px;border-bottom:1px dashed var(--line2);margin-bottom:8px}
.mark{width:34px;height:34px;border:1px solid var(--bezel);background:#2b1e14;position:relative;flex:none;
  display:grid;place-items:center;color:var(--teal);font-size:15px;font-weight:700;
  box-shadow:inset 0 0 12px rgba(216,166,87,.16)}
.mark:after{content:"▚";position:absolute;font-size:20px;color:#d8c3a0}
.brand b{display:block;font-size:12px;letter-spacing:.18em;text-transform:uppercase;color:#f0dfc2}
.brand span{display:block;font-size:9.5px;color:var(--muted);letter-spacing:.12em;margin-top:1px}
.navlabel{font-size:9px;letter-spacing:.2em;color:#d3bd98;text-transform:uppercase;padding:14px 6px 5px}
.nav{display:flex;align-items:center;gap:8px;width:100%;border:1px solid transparent;
  background:transparent;color:#e9d6b6;padding:7px 8px;text-align:left;font-size:11.5px;
  letter-spacing:.04em;transition:background .1s,color .1s}
.nav:hover{background:#2a1e14;color:#f3e3c6;border-color:#3c2b1d}
.nav.active{background:#33241a;border-color:var(--line2);color:#ffe9c2}
.nav .ico{width:22px;height:22px;display:grid;place-items:center;flex:none;
  font-size:9px;font-weight:700;color:#e6c79a;background:#2b1d13;
  clip-path:polygon(50% 0,100% 25%,100% 75%,50% 100%,0 75%,0 25%)}
.nav.active .ico{color:#120d09;background:var(--teal);border-color:var(--teal)}
.nav[data-nav="co2"] .ico{color:#ffab54;border-color:#7a4a25}
.nav[data-nav="retrieval"] .ico{color:#c6a8ff;border-color:#5a4670}
.nav[data-nav="physics"] .ico{color:#ff8f86;border-color:#6e4242}
.nav[data-nav="info"] .ico{color:#b4e35f;border-color:#59632f}
.nav[data-nav="uncertainty"] .ico{color:#e8c46a;border-color:#6e5a2a}
.nav[data-nav="catalog"] .ico{color:#8fd0a8;border-color:#3f6b50}
.nav[data-nav="formulas"] .ico{color:#d5b3ff;border-color:#57436f}
.sidefoot{margin-top:auto;padding:12px 5px 2px;border-top:1px dashed var(--line);font-size:9.5px;color:#d3bd98;line-height:1.7}
.sidefoot .dot{display:inline-block;width:6px;height:6px;background:var(--green);
  box-shadow:0 0 6px var(--green);margin-right:5px}
.main{min-width:0;padding:0 22px 44px;background:#100b07}
.topbar{display:flex;align-items:center;justify-content:space-between;gap:16px;height:54px;
  border-bottom:1px solid var(--line2);margin-bottom:20px;position:sticky;top:0;z-index:4;
  background:#100b07}
.crumb{display:flex;align-items:center;gap:8px;color:var(--muted);font-size:10px;letter-spacing:.12em;text-transform:uppercase}
.crumb b{color:var(--text);font-weight:700}
.pill{display:inline-flex;align-items:center;gap:6px;border:1px solid var(--line2);background:#1b130d;
  color:#d6bd97;padding:4px 9px;font-size:9.5px;letter-spacing:.12em;text-transform:uppercase}
.pill i{width:6px;height:6px;background:var(--green);box-shadow:0 0 6px var(--green);display:inline-block;
  animation:blink 2.4s steps(1) infinite}
@keyframes blink{50%{opacity:.25}}
.page{display:none}
.page.active{display:block}
.hero{border:1px solid var(--line2);padding:26px 28px;position:relative;overflow:hidden;
  background:#1a120c;display:grid;grid-template-columns:1.25fr .75fr;gap:24px;align-items:center;
  box-shadow:inset 0 0 40px rgba(0,0,0,.5)}
.eyebrow{font-size:9.5px;letter-spacing:.22em;text-transform:uppercase;color:var(--teal);font-weight:700}
h1{font-size:clamp(24px,3vw,38px);line-height:1.1;letter-spacing:-.02em;margin:8px 0 10px;max-width:760px;font-weight:700;
  text-transform:uppercase}
h1 em{font-style:normal;color:var(--teal)}
.lede{color:#e4d0ae;max-width:650px;font-size:12.5px;margin:0 0 18px}
.actions{display:flex;gap:8px;flex-wrap:wrap}
.btn{border:1px solid var(--line2);background:#241a12;color:var(--text);padding:8px 13px;
  font-size:10.5px;font-weight:700;letter-spacing:.09em;text-transform:uppercase;transition:.1s;
  display:inline-flex;align-items:center;gap:7px}
.btn:hover{border-color:var(--bezel);background:#332419;color:#fff0d4}
.btn.primary{background:var(--teal);border-color:var(--teal);color:#160f0a}
.btn.primary:hover{background:#ffcf8a;border-color:#ffcf8a}
.btn.ghost{background:transparent}
.btn:disabled{opacity:.45;cursor:not-allowed}
.stats{display:grid;grid-template-columns:repeat(4,minmax(0,1fr));gap:1px;margin:14px 0 24px;
  border:1px solid var(--line);background:var(--line)}
.stat{background:#1a120c;padding:12px 14px}
.stat .k{font-size:9px;letter-spacing:.14em;text-transform:uppercase;color:#cdb794}
.stat .v{font-size:23px;font-weight:700;margin-top:3px;font-variant-numeric:tabular-nums;color:var(--teal)}
.stat .s{font-size:10px;color:#c2ac8b;margin-top:1px}
.sectionhead{display:flex;align-items:flex-end;justify-content:space-between;gap:16px;margin:22px 0 10px}
.sectionhead h2{font-size:13px;letter-spacing:.12em;text-transform:uppercase;margin:0;font-weight:700;color:#f0dfc2}
.sectionhead h2:before{content:"■ ";color:var(--teal)}
.sectionhead p{margin:3px 0 0;color:var(--muted);font-size:10.5px}
.grid{display:grid;gap:10px}
.g2{grid-template-columns:repeat(2,minmax(0,1fr))}
.g3{grid-template-columns:repeat(3,minmax(0,1fr))}
.g4{grid-template-columns:repeat(4,minmax(0,1fr))}
.card{border:1px solid var(--line);background:#1a120c;padding:14px}
.card h3{margin:0 0 3px;font-size:11.5px;font-weight:700;letter-spacing:.06em;text-transform:uppercase;color:#f0dfc2}
.card .desc{color:var(--muted);font-size:10.5px;margin:0 0 10px}
.module{cursor:pointer;transition:.1s;position:relative;overflow:hidden}
.module:hover{border-color:var(--line2);background:#241a12}
.module .tag{width:27px;height:27px;display:grid;place-items:center;font-size:9px;font-weight:700;
  color:var(--teal);background:#241a12;border:1px solid var(--line2);margin-bottom:9px}
.module .arrow{position:absolute;right:12px;top:12px;color:#bda57e;font-size:14px}
.tool{display:grid;grid-template-columns:350px minmax(0,1fr);gap:10px;align-items:start}
.panel{border:1px solid var(--line);background:#18110c;overflow:hidden}
.panelhead{padding:9px 12px;border-bottom:1px solid var(--line2);display:flex;align-items:center;justify-content:space-between;gap:10px;background:#22180f}
.panelhead b{font-size:10.5px;letter-spacing:.1em;text-transform:uppercase;color:#f0dfc2}
.panelhead .hint{font-size:9px;color:#c9b493;letter-spacing:.06em;text-transform:uppercase}
.panelbody{padding:12px}
.field{margin-bottom:9px}
.field:last-child{margin-bottom:0}
.field label,.label{display:block;font-size:9px;letter-spacing:.12em;text-transform:uppercase;color:#d8c2a1;margin-bottom:4px;font-weight:700}
input[type=text],input[type=number],select,textarea{width:100%;background:#0e0a07;border:1px solid #4a3524;
  padding:7px 8px;font-size:11.5px;outline:none;color:#ead8ba}
input[type="password"],input[type="text"],input[type="number"],select,textarea{background:#120c08;color:var(--text)}
input[type="password"]{-webkit-text-security:disc}
input::placeholder,textarea::placeholder{color:#a08a6c}
input:focus,select:focus,textarea:focus{border-color:var(--teal);box-shadow:0 0 0 1px rgba(216,166,87,.25)}
textarea{font-size:10.5px;line-height:1.5;resize:vertical;min-height:96px}
select{appearance:none;background-image:url("data:image/svg+xml,%3Csvg xmlns='http://www.w3.org/2000/svg' width='8' height='6'%3E%3Cpath d='M0 0h8L4 6z' fill='%23a38562'/%3E%3C/svg%3E");
  background-position:calc(100% - 10px) center;background-repeat:no-repeat;padding-right:24px}
.row{display:grid;grid-template-columns:1fr 1fr;gap:8px}
.row3{display:grid;grid-template-columns:repeat(3,1fr);gap:8px}
.check{display:flex;align-items:center;gap:7px;font-size:10.5px;color:#e4d0ae;margin:6px 0;cursor:pointer}
.check input{accent-color:var(--teal);width:13px;height:13px}
.run{width:100%;justify-content:center;margin-top:5px}
.result{min-height:180px}
.placeholder{min-height:160px;display:grid;place-items:center;text-align:center;color:#c2ac8b;padding:24px}
.placeholder .big{font-size:26px;color:#5a4330;margin-bottom:6px}
.placeholder p{max-width:340px;font-size:10.5px;margin:0}
.cards{display:grid;grid-template-columns:repeat(4,minmax(0,1fr));gap:1px;margin-bottom:10px;border:1px solid var(--line);background:var(--line)}
.metric{background:#1a120c;padding:8px 9px}
.metric .k{font-size:8.5px;text-transform:uppercase;letter-spacing:.12em;color:#cdb794}
.metric .v{font-size:17px;font-weight:700;margin-top:2px;font-variant-numeric:tabular-nums;color:#eedcbe}
.metric .v.good{color:var(--green)} .metric .v.warn{color:var(--amber)} .metric .v.bad{color:var(--rose)}
.metric .s{font-size:9.5px;color:#c2ac8b;margin-top:1px}
.verdict{border:1px solid var(--line2);background:#22180f;
  padding:8px 10px;font-size:11px;margin-bottom:10px;color:#f2e2c4}
.verdict b{color:var(--teal)}
.chartwrap{border:1px solid var(--line);background:#0d0906;padding:5px 4px 2px;margin-bottom:8px}
.chart{width:100%;height:auto;display:block;overflow:visible}
.chart text{fill:#c9b493;font-family:var(--mono);font-size:9px}
.chart .gridline{stroke:#3a2a1d;stroke-width:1;stroke-dasharray:2 5}
.chart .axis{stroke:#4a3524;stroke-width:1}
.chart .line{fill:none;stroke-width:1.15;stroke-linejoin:round;stroke-linecap:round;stroke-opacity:.92}
.chart .area{stroke:none;fill:rgba(255,154,77,.06)}
.chart .legend{font-family:var(--mono);font-size:9px;fill:#dcc6a4}
.ph-picks{display:grid;grid-template-columns:1fr 1fr;gap:4px;border:1px solid #3c2b1d;background:#0d0906;padding:5px;max-height:190px;overflow:auto}
.ph-picks label{display:flex;gap:5px;align-items:flex-start;font-size:9px;color:#cbb18d;padding:3px;cursor:pointer}
.ph-picks label:hover{background:#1d140d}
.ph-card{border:1px solid #3c2b1d;background:#120d08;padding:7px;margin-top:7px}
.ph-card .ph-title{font-size:9px;color:#ffe2b0;letter-spacing:.06em;text-transform:uppercase;margin-bottom:5px}
.ph-card .ph-grid{display:grid;grid-template-columns:1fr 1fr;gap:5px}
.lab-group{border:1px solid #33241a;background:#0d0906;margin-bottom:6px}
.lab-grouphead{display:flex;gap:6px;align-items:center;padding:4px 6px;background:#171009;border-bottom:1px solid #33241a;font-size:9px;letter-spacing:.08em;text-transform:uppercase;color:#ffe2b0;cursor:pointer}
.lab-grouphead .hint{margin-left:auto;color:#d8c2a1;letter-spacing:.04em}
.lab-row{display:flex;gap:6px;align-items:baseline;padding:2px 6px;font-size:9px;color:#cbb18d;cursor:pointer;border-bottom:1px solid #1a1209}
.lab-row:hover{background:#1d140d}
.lab-row .lab-name{color:#f7e8c8;min-width:170px}
.lab-row .lab-meta{color:#d3bd98;margin-left:auto;text-align:right}
.filebtn{display:inline-block;margin-top:7px;border:1px solid var(--line2);background:#22180f;color:#f0dfc0;padding:5px 9px;font-size:10px;letter-spacing:.08em;text-transform:uppercase;cursor:pointer}
.filebtn:hover{background:#2b1e14;color:#fff3dc}
.vhidden{position:absolute;width:1px;height:1px;padding:0;margin:-1px;overflow:hidden;clip:rect(0 0 0 0);white-space:nowrap;border:0}
.mission{border:1px solid var(--line2);background:#1a120c;padding:16px 18px;margin-bottom:10px;
  box-shadow:inset 0 0 40px rgba(0,0,0,.5);position:relative}
.mission .mtag{font-size:9px;letter-spacing:.2em;text-transform:uppercase;color:var(--teal);font-weight:700}
.mission p{margin:6px 0 10px;color:#e8d6b6;font-size:12px;max-width:1000px}
.mission .qs{display:grid;grid-template-columns:repeat(4,minmax(0,1fr));gap:8px}
.mq{border:1px solid var(--line);background:#120c08;padding:8px 10px}
.mq b{display:block;font-size:10.5px;letter-spacing:.14em;color:var(--amber)}
.mq span{font-size:9.5px;color:#cbb493}
.crow{display:grid;grid-template-columns:1.5fr 1fr 1fr 1fr auto;gap:8px;align-items:end}
.cxgrid{display:grid;grid-template-columns:1fr 1fr;gap:10px;margin-top:10px}
.cxcard{border:1px solid var(--line);background:#120c08;padding:10px}
.cxhead{display:flex;justify-content:space-between;gap:8px;align-items:center;margin-bottom:6px}
.cxhead b{font-size:11px;letter-spacing:.06em;color:#f4e8d2}
.cxsig{font-size:8.5px;letter-spacing:.12em;padding:2px 6px;border:1px solid}
.cxsig.yes{color:#0b0806;background:var(--amber);border-color:var(--amber);font-weight:700}
.cxsig.no{color:#d8c2a1;border-color:var(--line2)}
.cxarrow{font-size:12px;font-weight:700;color:#ffd79a;margin-bottom:6px}
.cxarrow span{color:#cbb493;font-weight:400;font-size:9.5px;letter-spacing:.06em;text-transform:uppercase}
.cxstats{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:6px;margin-bottom:8px}
.cxstats span{display:block;font-size:8.5px;letter-spacing:.1em;text-transform:uppercase;color:#c2ac8b}
.cxstats b{font-size:11.5px;color:#f7e8c8}
.cxverdict{border:1px solid var(--line2);background:#22180f;padding:9px 12px;font-size:11px;color:#f2dfbe;margin-top:10px}
.cxverdict.opposed{border-color:var(--amber);color:#ffe0a6}
summary.fhead{cursor:pointer;list-style:none}
summary.fhead::-webkit-details-marker{display:none}
.bridge-group{margin-bottom:10px}
.bridge-q{font-size:9px;letter-spacing:.16em;color:var(--amber);font-weight:700;margin-bottom:4px}
.bridge-group .etable{width:100%}
.coach{position:fixed;z-index:10000;max-width:320px;background:#2a1c12;border:1px solid var(--amber);
  box-shadow:0 14px 40px rgba(0,0,0,.7);padding:12px 13px 10px;color:#f7e8c8}
.coach[hidden]{display:none}
.coach-arrow{position:absolute;width:12px;height:12px;background:#2a1c12;border-left:1px solid var(--amber);
  border-bottom:1px solid var(--amber);transform:rotate(45deg)}
.coach-step{font-size:8.5px;letter-spacing:.18em;text-transform:uppercase;color:#ffc37f;margin-bottom:5px}
.coach-text{font-size:11px;line-height:1.55;color:#efdcc0}
.coach-actions{display:flex;justify-content:space-between;gap:8px;margin-top:10px}
.coach-actions .btn{padding:5px 9px;font-size:9px}
.coach-focus{outline:2px solid var(--amber)!important;outline-offset:2px;box-shadow:0 0 0 4px rgba(255,195,94,.18)!important}
.qbar{display:grid;grid-template-columns:repeat(4,minmax(0,1fr));gap:8px;margin-bottom:10px}
.qcell{border:1px solid var(--line2);background:#1a120c;padding:8px 10px}
.qcell span{display:block;font-size:8.5px;letter-spacing:.14em;color:var(--teal);font-weight:700}
.qcell b{display:block;font-size:11px;color:#f4e8d2;margin-top:3px;line-height:1.35;font-weight:400}
.qcell.set b{color:#ffe0a6}
.varcards{display:grid;grid-template-columns:repeat(auto-fill,minmax(320px,1fr));gap:10px;padding:10px}
.varcard{border:1px solid var(--line);background:#140e09;padding:11px 12px;display:flex;flex-direction:column;gap:7px}
.varcard .vtop{display:flex;justify-content:space-between;gap:8px;align-items:baseline}
.varcard .vid{font-size:11.5px;font-weight:700;color:#f7e8c8;letter-spacing:.04em}
.varcard .vname{font-size:9.5px;color:#c2ac8b}
.varcard .vchange{font-size:11px;color:#ffc37f;line-height:1.4}
.varcard .vmeta{font-size:9px;color:#bda98c;letter-spacing:.04em;text-transform:uppercase}
.varcard .vregions{font-size:9px;color:#a9917a}
.varcard .acts{display:flex;gap:6px;margin-top:2px}
.fetch-extra{margin-top:14px}
.fetch-extra .sectionhead{margin-top:10px}
@media(max-width:1100px){.qbar{grid-template-columns:1fr 1fr}}
@media(max-width:1100px){.crow{grid-template-columns:1fr 1fr}.cxgrid{grid-template-columns:1fr}
  .mission .qs{grid-template-columns:repeat(2,1fr)}}
.lab-results{display:flex;flex-direction:column;gap:5px}
.lab-results .lab-hit{border:1px solid #33241a;background:#120d08;padding:6px}
.lab-results .lab-hit.bad{border-color:#6e4242}
.lab-results .lab-hit .ph-title{display:flex;gap:8px;align-items:baseline}
.lab-results pre{margin:3px 0 0;font:11px/1.45 ui-monospace,Menlo,Consolas,monospace;color:#f3dfbd;white-space:pre-wrap;word-break:break-word}
.lab-results .err{color:#ff9b7d}
.vis-grid{display:grid;grid-template-columns:1fr 1fr;gap:10px;margin-top:10px}
.live-grid{display:grid;grid-template-columns:1fr 1fr;gap:10px;margin-top:10px}
.livecard .panelbody{padding:8px}
.visbar{display:grid;grid-template-columns:250px 1fr 54px;align-items:center;gap:9px;font-size:10px;color:#e9d6b6;margin:4px 0}
.visbar .bar{display:block;height:10px;background:#0d0906;border:1px solid #3c2b1d}
.visbar .bar i{display:block;height:100%}
.visbar b{text-align:right;color:#f7e8c8}
.viscard .chart{height:auto}
#vis-map .gridmap{grid-template-columns:repeat(24,1fr)}
@media(max-width:1100px){.vis-grid{grid-template-columns:1fr}.visbar{grid-template-columns:140px 1fr 44px}}
.fgroup{font-size:9px;letter-spacing:.1em;text-transform:uppercase;color:#ffc07a;border-bottom:1px solid #4a3524;margin:10px 0 4px;padding:3px 2px}
.etable{width:100%;border-collapse:collapse;font-family:var(--mono);font-size:9px}
.etable th{text-align:left;color:#cbb390;font-weight:400;letter-spacing:.08em;text-transform:uppercase;border-bottom:1px solid #3c2b1d;padding:4px 5px}
.etable td{padding:5px;border-bottom:1px dashed #2f2218;color:#cbb18d}
.etable tr.primary td{color:#ffe2b0;background:#24180e}
.etable .pick{color:#ffc07a;font-weight:700}
.etable .fallback{color:#cbb390}
.minicharts{display:grid;grid-template-columns:1fr 1fr;gap:8px}
.minicharts .chartwrap{margin:0}
.tagrow{display:flex;gap:4px;flex-wrap:wrap;margin-top:6px}
.tag{display:inline-flex;align-items:center;border:1px solid var(--line2);background:#22180f;color:#e9d6b6;
  padding:2px 6px;font-size:9px;letter-spacing:.04em;text-transform:uppercase}
.tag.teal{border-color:#7b5b32;color:#ffcf8a;background:#2a1d10}
.tag.amber{border-color:#6e4a2c;color:#ffc861;background:#2a1d10}
.warnlist{margin:8px 0 0;padding:0;list-style:none}
.warnlist li{font-size:10.5px;color:#d8b77e;background:#241a12;border-left:2px solid var(--amber);
  padding:5px 8px;margin-top:4px}
details{margin-top:8px;border-top:1px dashed #3c2b1d;padding-top:7px}
summary{cursor:pointer;color:#d8c2a1;font-size:10px;letter-spacing:.06em}
pre{font-family:var(--mono);font-size:10px;line-height:1.55;background:#0a0705;border:1px solid #3c2b1d;
  padding:9px;overflow:auto;max-height:300px;color:#f0dfc0;white-space:pre-wrap;overflow-wrap:anywhere}
.tabs{display:flex;gap:2px;border-bottom:1px solid var(--line2);padding:0 2px;margin-bottom:12px;flex-wrap:wrap}
.tab{background:transparent;border:1px solid transparent;border-bottom:0;color:#d8c2a1;padding:7px 10px;
  font-size:10px;font-weight:700;letter-spacing:.06em;text-transform:uppercase}
.tab:hover{color:var(--text);background:#22180f}
.tab.active{color:var(--teal);border-color:var(--line2);background:#22180f}
.tabpane{display:none}.tabpane.active{display:block}
.split{display:grid;grid-template-columns:1fr 1fr;gap:12px}
.bigvalue{font-size:31px;letter-spacing:-.03em;color:var(--teal);line-height:1;
  font-variant-numeric:tabular-nums;margin:2px 0 4px}
.unit{font-size:10px;color:var(--muted);text-transform:uppercase;letter-spacing:.08em}
.kv{display:grid;grid-template-columns:1fr auto;gap:5px 12px;font-size:10.5px;margin-top:9px}
.kv .k{color:#cdb794}.kv .v{color:#f2e2c4;text-align:right}
.coeffs{display:grid;grid-template-columns:repeat(4,1fr);gap:4px;margin-top:8px}
.coeff{border:1px solid var(--line);padding:5px 6px;background:#140e09}
.coeff b{display:block;font-size:8.5px;color:#c9b493;font-weight:400}
.coeff span{font-size:10px;color:#f2e2c4}
.table{width:100%;border-collapse:collapse;font-size:10.5px}
.table th{text-align:left;font-size:8.5px;letter-spacing:.12em;text-transform:uppercase;color:#cdb794;
  border-bottom:1px solid var(--line2);padding:6px;font-weight:700}
.table td{border-bottom:1px solid #2e2116;padding:6px;color:#e9d6b6;vertical-align:top}
.table tr:hover td{background:#22180f}
.table td.id{color:var(--teal);white-space:nowrap}
.table td.mut{color:#d8c2a1}
.searchrow{display:flex;gap:8px;margin-bottom:10px;flex-wrap:wrap}
.searchrow input{flex:1;min-width:160px}.searchrow select{width:auto;min-width:140px}
.empty{padding:28px;text-align:center;color:#c2ac8b;font-size:10.5px}
.spinner{display:inline-block;width:10px;height:10px;border:2px solid rgba(0,0,0,.3);
  border-top-color:#160f0a;animation:spin .7s linear infinite;vertical-align:-1px;margin-right:6px}
@keyframes spin{to{transform:rotate(360deg)}}
.toast{position:fixed;right:22px;bottom:22px;z-index:10000;max-width:380px;border:1px solid #8a4a32;
  background:#2d150f;color:#ffd8c8;padding:9px 12px;font-size:10.5px;box-shadow:0 0 24px rgba(0,0,0,.7);display:none}
.toast.show{display:block}
.formula{font-size:9.5px;color:#dcc6a4;background:#0d0906;border:1px dashed #3c2b1d;padding:6px 8px;margin-top:8px;line-height:1.6}
.note{font-size:9.5px;color:#c2ac8b;margin-top:7px;line-height:1.6}
.hero .scope{border:1px solid var(--line2);background:#0a0705;padding:7px;min-height:230px;display:flex;flex-direction:column}
.scopehead{display:flex;justify-content:space-between;font-size:8.5px;letter-spacing:.12em;color:#c9b493;
  border-bottom:1px dashed #3c2b1d;padding-bottom:5px;text-transform:uppercase}
.scope svg{width:100%;height:auto;flex:1}
.scopefoot{font-size:8.5px;color:#bda57e;border-top:1px dashed #3c2b1d;padding-top:4px;letter-spacing:.08em}
.workflow{display:grid;grid-template-columns:272px minmax(0,1fr) 290px;gap:10px;align-items:start}
.wleft,.wright{display:flex;flex-direction:column;gap:10px;position:sticky;top:56px}
.wcenter{display:flex;flex-direction:column;gap:8px;min-width:0}
.step{border:1px solid var(--line);background:#18110c}
.step-head{display:flex;align-items:center;gap:9px;padding:8px 10px;border-bottom:1px solid var(--line2);background:#22180f}
.step-no{width:26px;height:26px;display:grid;place-items:center;font-size:10px;
  font-weight:700;color:#ffcf8a;flex:none;background:#2b1d13;
  clip-path:polygon(30% 0,70% 0,100% 30%,100% 70%,70% 100%,30% 100%,0 70%,0 30%)}
.step-title{font-size:10.5px;letter-spacing:.1em;text-transform:uppercase;font-weight:700;color:#f0dfc2}
.step-status{margin-left:auto;font-size:8.5px;letter-spacing:.09em;text-transform:uppercase;color:#c9b493}
.step-body{padding:11px}
.step:not(.open) .step-body{display:none}
.step.locked{opacity:.55}
.step.done .step-no{background:var(--teal);color:#160f0a}
.step.done .step-status{color:var(--green)}
.facts{display:grid;grid-template-columns:repeat(3,1fr);gap:1px;background:var(--line);border:1px solid var(--line)}
.fact{background:#140e09;padding:7px 8px}
.fact .k{font-size:8px;letter-spacing:.1em;text-transform:uppercase;color:#c9b493}
.fact .v{font-size:11px;color:#e0c9a6;margin-top:2px}
.fact .s{font-size:9px;color:#c2ac8b}
.bar{height:6px;background:#2a1e14;border:1px solid #3c2b1d;margin:3px 0 7px}
.bar i{display:block;height:100%;background:var(--teal)}
.bar.warn i{background:var(--rose)}
.log{max-height:250px;overflow:auto;font-size:9.5px;line-height:1.55}
.logline{border-bottom:1px dotted #3c2b1d;padding:3px 0;color:#dcc6a4}
.logline b{color:var(--teal);font-weight:400}
.logline.warn{color:#ffc861}
.logline.err{color:#ff9b7d}
.artlist{margin:0;padding:0;list-style:none;font-size:9.5px;color:#e9d6b6}
.artlist li{border-left:2px solid var(--rose);padding:3px 7px;margin-bottom:4px;background:#22130e}
.varpick{max-height:300px;overflow:auto;border:1px solid var(--line);background:#0e0a07}
.varitem{display:flex;gap:7px;align-items:flex-start;padding:6px 7px;border-bottom:1px solid #2a1e14;cursor:pointer}
.varitem:hover{background:#22180f}
.varitem.active{background:#2b1e14;border-left:2px solid var(--teal)}
.varitem .code{font-size:8.5px;color:var(--teal);border:1px solid var(--line2);padding:1px 3px;flex:none;margin-top:1px}
.varitem .nm{display:block;font-size:10px;color:#f2e2c4}
.varitem .meta{display:block;font-size:8.5px;color:#c2ac8b;text-transform:uppercase;letter-spacing:.04em;margin-top:1px}
.sphere{font-size:8.5px;letter-spacing:.14em;text-transform:uppercase;color:#c9b493;background:#22180f;
  border-top:1px solid var(--line);border-bottom:1px solid var(--line);padding:4px 7px;margin-top:2px}
.cmdbar{display:flex;gap:6px}
.cmdbar input{flex:1}
.console{font-size:10px;line-height:1.55;background:#0a0705;border:1px solid #3c2b1d;padding:9px;
  min-height:260px;max-height:480px;overflow:auto;white-space:pre-wrap;color:#f0dfc0}
.console .ok{color:var(--green)} .console .bad{color:#ff9b7d} .console .cmd{color:var(--teal)}
.quick{display:flex;gap:4px;flex-wrap:wrap;margin-top:6px}
.quick button{font-size:9px;padding:3px 6px;background:#22180f;border:1px solid var(--line2);color:#e9d6b6;text-transform:uppercase}
.quick button:hover{color:#ffcf8a;border-color:var(--teal)}
.flist{max-height:none;overflow:visible}
.fitem{border-bottom:1px solid #2a1e14}
.fhead{display:flex;align-items:center;gap:7px;padding:5px 7px;cursor:pointer;font-size:10px;color:#f2e2c4}
.fhead:hover{background:#22180f}
.fhead .cnt{margin-left:auto;font-size:8.5px;color:#c9b493}
.fbody{display:none;padding:4px 7px 7px;background:#0e0a07}
.fitem.open .fbody{display:block}
.fsig{font-size:9px;color:#dcc6a4;padding:2px 0;border-bottom:1px dotted #2a1e14;word-break:break-word}
.gridmap{display:grid;gap:1px;background:#2a1e14;border:1px solid #4a3524;padding:1px}
.gridcell{height:13px;min-width:2px;background:#2a1e14}
.dsgrid{display:grid;grid-template-columns:repeat(2,1fr);gap:8px}
.dscard{border:1px solid var(--line);background:#1a120c;padding:10px}
.dscard h3{font-size:10.5px;margin:0 0 3px;color:#f0dfc2;text-transform:uppercase;letter-spacing:.05em}
.dscard .id{font-size:8.5px;color:var(--teal);word-break:break-all}
.dscard p{font-size:9.5px;color:#d8c2a1;margin:5px 0}
.dscard .meta{font-size:8.5px;color:#c9b493;text-transform:uppercase;letter-spacing:.05em}
.telemetry{font-size:8.5px;color:#c2ac8b;letter-spacing:.08em;text-transform:uppercase}
.statusline{display:flex;gap:12px;flex-wrap:wrap;border:1px solid var(--line);background:#0d0906;
  padding:5px 8px;font-size:8.5px;color:#c9b493;letter-spacing:.08em;text-transform:uppercase;margin-bottom:8px}
.statusline b{color:var(--green);font-weight:400}
.mapbox{border:1px solid var(--line2);background:#0a0705;padding:8px;min-height:190px;position:relative}
.mapbox .grid{position:absolute;inset:8px;background-image:url("data:image/svg+xml,%3Csvg xmlns='http://www.w3.org/2000/svg' width='24' height='24'%3E%3Cpath d='M24 0H0v24' fill='none' stroke='%232a1e14'/%3E%3C/svg%3E")}
::-webkit-scrollbar{width:10px;height:10px}
::-webkit-scrollbar-track{background:#0a0705;border-left:1px solid #2a1e14}
::-webkit-scrollbar-thumb{background:#4a3524;border:1px solid #6b4a30}
::-webkit-scrollbar-thumb:hover{background:#6b4a30}
@media(max-width:1450px){.workflow{grid-template-columns:240px minmax(0,1fr) 260px}}
@media(max-width:1200px){.workflow{grid-template-columns:250px minmax(0,1fr)}.wright{grid-column:1/-1;position:static;display:grid;grid-template-columns:1fr 1fr 1fr}.hero{grid-template-columns:1fr}.stats{grid-template-columns:repeat(2,1fr)}.g4{grid-template-columns:repeat(2,1fr)}}
@media(max-width:1100px){.workflow{grid-template-columns:1fr}.wleft,.wright{position:static}.wleft{display:grid;grid-template-columns:1fr 1fr;align-items:start}}
@media(max-width:900px){.app{grid-template-columns:1fr;border-width:6px}.sidebar{position:static;height:auto;
  flex-direction:row;flex-wrap:wrap;padding:8px;border-right:0;border-bottom:2px solid var(--line)}
  .brand{width:100%;padding-bottom:8px;margin-bottom:2px}.navlabel,.sidefoot{display:none}.nav{width:auto;padding:5px 6px;font-size:10px}
  .nav .ico{display:none}.main{padding:0 10px 30px}.workflow{grid-template-columns:1fr}.wleft,.wright{position:static}
  .wright{grid-template-columns:1fr}.tool,.split,.g2,.g3,.g4,.minicharts,.dsgrid{grid-template-columns:1fr}
  .cards{grid-template-columns:1fr 1fr}.stats{grid-template-columns:1fr 1fr}.hero{padding:18px}
  body{padding:6px}.topbar{height:44px}}

</style>
</head>
<body>
<div class="app">
<aside class="sidebar">
  <div class="brand"><div class="mark"></div><div><b>Trend Detective</b><span>NASA climate workbench</span></div></div>
  <div class="navlabel">The challenge</div>
  <button class="nav" data-nav="live"><span class="ico">1</span>Fetch data</button>
  <button class="nav" data-nav="catalog"><span class="ico">2</span>Select variables</button>
  <button class="nav" data-nav="visual"><span class="ico">3</span>Visualize over time</button>
  <button class="nav active" data-nav="workflow"><span class="ico">4</span>What / where / how much</button>
  <button class="nav" data-nav="significance"><span class="ico">5</span>Significance</button>
  <div class="navlabel">Deep dives</div>
  <button class="nav" data-nav="timeseries"><span class="ico">TS</span>Time series lab</button>
  <button class="nav" data-nav="co2"><span class="ico">CO</span>CO₂ curve fit</button>
  <button class="nav" data-nav="retrieval"><span class="ico">RS</span>Retrievals</button>
  <button class="nav" data-nav="physics"><span class="ico">PH</span>Formula calculator</button>
  <button class="nav" data-nav="info"><span class="ico">MI</span>Information theory</button>
  <button class="nav" data-nav="uncertainty"><span class="ico">σ</span>Uncertainty</button>
  <div class="navlabel">Engine</div>
  <button class="nav" data-nav="console"><span class="ico">&gt;_</span>Command console</button>
  <button class="nav" data-nav="formulas"><span class="ico">FX</span>Calculation modes</button>
  <button class="nav" data-nav="alllab"><span class="ico">Σ</span>All calculations</button>
  <button class="nav" data-nav="overview"><span class="ico">SY</span>System</button>
  <div class="sidefoot"><span class="dot"></span>LIVE NASA ENDPOINTS<br><span id="foot-catalog">—</span> catalog variables<br><span id="foot-formulas">—</span> public functions</div>
</aside>
<main class="main">
  <header class="topbar">
    <div class="crumb"><span>TREND DETECTIVE</span><span>/</span><b id="crumb">What / where / how much</b></div>
    <div class="pill"><i></i>change engine ready</div>
  </header>

  <!-- TREND DETECTIVE WORKFLOW -->
  <section class="page active" id="page-workflow">
    <div class="sectionhead" style="margin-top:0"><div><h2>What, where, how much, significant</h2><p>Change is always occurring in Earth's interconnected system, and a measured variable can trend one way in one region and the opposite way in another. Select a mission or model variable, see its time evolution, then read what is changing, where, how much, and whether it is significant. Nothing is tested until you confirm the data.</p></div><div class="telemetry">NASA / EARTH OBSERVATION / LIVE</div></div>
    <div class="qbar" id="wf-qbar">
      <div class="qcell"><span>WHAT IS CHANGING</span><b id="wf-q-what">run an analysis</b></div>
      <div class="qcell"><span>WHERE IS IT CHANGING</span><b id="wf-q-where">run an analysis</b></div>
      <div class="qcell"><span>HOW MUCH</span><b id="wf-q-howmuch">run an analysis</b></div>
      <div class="qcell"><span>IS IT SIGNIFICANT</span><b id="wf-q-sig">run an analysis</b></div>
    </div>

    <div class="panel" style="margin-bottom:10px" id="wf-bridge">
      <div class="panelhead"><b>Formulas behind each answer</b><span class="hint">many formulas, one question each</span></div>
      <div class="panelbody" id="wf-formula-bridge">
        <div class="empty">Run an analysis and every formula that contributes to WHAT, WHERE, HOW MUCH and SIGNIFICANT appears here with its own number.</div>
      </div>
    </div>

    <div class="panel" style="margin-bottom:10px" id="wf-plan-panel">
      <div class="panelhead"><b>Formulas selected for this variable</b><span class="hint" id="wf-plan-hint">chosen by the engine from the registry</span></div>
      <div class="panelbody" id="wf-plan">
        <div class="empty">The engine reads the variable, adds its retrieval formulas and the seasonality / trend / significance / power chain that applies to it, then runs them as one batch.</div>
      </div>
    </div>

    <div class="panel" style="margin-bottom:10px">
      <div class="panelhead"><b>Regional contrast · one variable, two places</b><span class="hint">live NASA POWER</span></div>
      <div class="panelbody">
        <div class="crow">
          <div class="field"><label>Variable (POWER parameter)</label><select id="cx-param">
            <option value="T2M">T2M · air temperature at 2 m</option>
            <option value="T2M_MAX">T2M_MAX · daily maximum temperature</option>
            <option value="T2M_MIN">T2M_MIN · daily minimum temperature</option>
            <option value="TS">TS · skin temperature</option>
            <option value="PRECTOTCORR">PRECTOTCORR · precipitation</option>
            <option value="WS2M">WS2M · wind speed at 2 m</option>
            <option value="RH2M">RH2M · relative humidity at 2 m</option>
            <option value="ALLSKY_SFC_SW_DWN">ALLSKY_SFC_SW_DWN · downward shortwave</option>
          </select></div>
          <div class="field"><label>Region A · lat / lon</label><div class="row"><input id="cx-lat1" type="number" value="72" step="0.5"><input id="cx-lon1" type="number" value="-40" step="0.5"></div></div>
          <div class="field"><label>Region B · lat / lon</label><div class="row"><input id="cx-lat2" type="number" value="-70" step="0.5"><input id="cx-lon2" type="number" value="0" step="0.5"></div></div>
          <div class="field"><label>Years</label><div class="row"><input id="cx-start" type="number" value="1990"><input id="cx-end" type="number" value="2024"></div></div>
          <button class="btn primary run" id="cx-run">Compare regions →</button>
        </div>
        <div class="note">Each region is a real NASA POWER monthly request (no credentials needed). Sen's slope and the seasonal Mann–Kendall test give each region its own answer, so an opposition between regions becomes measurable instead of hidden inside one global average.</div>
        <div id="cx-out"><div class="empty">Run the comparison to see whether the same variable rises in one region and falls in another.</div></div>
      </div>
    </div>

    <div class="workflow">
      <aside class="wleft">
        <div class="panel"><div class="panelhead"><b>Variable selection</b><span class="hint">Step 0</span></div><div class="panelbody">
          <div class="field"><label>Search catalog</label><input id="wf-search" type="text" placeholder="sst, ozone, soil moisture…"></div>
          <div class="field"><label>Earth sphere</label><select id="wf-sphere"><option value="">All spheres</option><option>ATMOSPHERE</option><option>OCEAN</option><option>CRYOSPHERE</option><option>LAND / BIOSPHERE</option></select></div>
          <div class="varpick" id="wf-vars"><div class="empty">Loading catalog…</div></div>
          <div class="field" style="margin-top:9px"><label>Region</label><select id="wf-region"></select><div class="note">For scalar demo series the region is a report tag; it does not subset the CSV. Gridded acquisitions carry their own latitude/longitude mask.</div></div>
          <div class="row"><div class="field"><label>Start year</label><input id="wf-start" type="number" value="2000" step="1"></div><div class="field"><label>End year</label><input id="wf-end" type="number" value="2024" step="1"></div></div>
          <button class="btn primary run" id="wf-load">Load data →</button>
          <input type="file" id="wf-upload" accept=".csv,text/csv" style="display:none">
          <div class="note">Local demo CSVs are synthetic fixtures. NASA product routes are real; live acquisition is planned on the NASA Data page.</div>
        </div></div>
        <div class="panel"><div class="panelhead"><b>Artifact advisory</b><span class="hint">Disclose</span></div><div class="panelbody"><ul class="artlist" id="wf-artifacts"><li>Select a variable to load its retrieval caveats.</li></ul></div></div>
      </aside>

      <div class="wcenter">
        <div class="step open" id="wf-step1">
          <div class="step-head"><div class="step-no">01</div><div class="step-title">Data acquisition</div><div class="step-status" id="wf-step1-status">AWAITING LOAD</div></div>
          <div class="step-body">
            <div id="wf-data"><div class="empty">Choose a variable and press LOAD DATA. Coverage, gaps and quality flags appear here before any statistic is computed.</div></div>
            <div id="wf-quality"></div>
            <label class="check"><input id="wf-proceed" type="checkbox" checked> Proceed with available data; gaps will be disclosed in the report</label>
            <button class="btn primary run" id="wf-confirm" disabled>Data confirmed → remove seasonality</button>
          </div>
        </div>

        <div class="step locked" id="wf-step2">
          <div class="step-head"><div class="step-no">02</div><div class="step-title">Seasonality removal</div><div class="step-status" id="wf-step2-status">LOCKED</div></div>
          <div class="step-body">
            <div class="statusline"><span>METHOD <b id="wf-method">STL</b></span><span>PERIOD <b id="wf-period">12</b></span><span>MODEL <b>Y = T + S + R</b></span></div>
            <div id="wf-decomp"><div class="empty">The decomposition appears after confirming the data.</div></div>
          </div>
        </div>

        <div class="step locked" id="wf-step3">
          <div class="step-head"><div class="step-no">03</div><div class="step-title">Trend + significance</div><div class="step-status" id="wf-step3-status">LOCKED</div></div>
          <div class="step-body">
            <div id="wf-trend"><div class="empty">Mann–Kendall and Sen's slope appear here.</div></div>
            <details><summary>Optional advanced tests (not required)</summary><div class="note">Acceleration, BEAST changepoints and field-significance corrections are implemented in the engine, but are deliberately not part of the core verdict. Run them from the Command Console.</div></details>
          </div>
        </div>

        <div class="step locked" id="wf-step4">
          <div class="step-head"><div class="step-no">04</div><div class="step-title">Spatial analysis</div><div class="step-status" id="wf-step4-status">LOCKED</div></div>
          <div class="step-body">
            <div id="wf-spatial"><div class="empty">A gridded field is required for map-wide inference.</div></div>
          </div>
        </div>

        <div class="step locked" id="wf-step5">
          <div class="step-head"><div class="step-no">05</div><div class="step-title">Final report</div><div class="step-status" id="wf-step5-status">LOCKED</div></div>
          <div class="step-body"><div id="wf-report"><div class="empty">Confirm the trend to assemble the report.</div></div></div>
        </div>
      </div>

      <aside class="wright">
        <div class="panel"><div class="panelhead"><b>Analysis log</b><span class="hint">Live</span></div><div class="panelbody"><div class="log" id="wf-log"><div class="logline">[BOOT] Workbench ready. Awaiting data.</div></div></div></div>
        <div class="panel"><div class="panelhead"><b>Active formula</b><span class="hint">Audit</span></div><div class="panelbody" id="wf-active"><div class="formula">Select a variable to see the retrieval and trend formulas.</div></div></div>
        <div class="panel"><div class="panelhead"><b>Inputs / outputs</b><span class="hint">Trace</span></div><div class="panelbody"><div class="kv" id="wf-io"><span class="k">STATUS</span><span class="v">IDLE</span></div></div></div>
      </aside>
    </div>
  </section>


  <!-- SIGNIFICANCE -->
  <section class="page" id="page-significance">
    <div class="sectionhead" style="margin-top:0"><div><h2>Is it significant?</h2><p>Step 5. A change is only reported when it survives the tests that could explain it away: serial correlation, the seasonal cycle, estimator disagreement, and the finite length of the record. Every test below is labelled with the question it answers, and its statistic and p-value are computed by the C++ engine from the actual series.</p></div><div class="telemetry">ALPHA / POWER / AUTOCORRELATION</div></div>
    <div class="tool">
      <div class="panel"><div class="panelhead"><b>Test setup</b><span class="hint">one variable at a time</span></div><div class="panelbody">
        <div class="field"><label>Variable</label><select id="sig-var"></select></div>
        <div class="field"><label>Significance level (alpha)</label><input id="sig-alpha" type="number" value="0.05" step="0.01" min="0.001" max="0.2"></div>
        <button class="btn primary run" id="sig-run">Run the significance suite →</button>
        <div class="note">The suite runs Mann–Kendall on the deseasonalised series, the seasonal Mann–Kendall test, the Sen slope with its rank-based confidence interval, the estimator ensemble, and the minimum detectable trend for a record of this length. Nothing is hard-coded: each number comes from the series you selected.</div>
      </div></div>
      <div id="sig-out"><div class="panel"><div class="placeholder"><div><div class="big">α</div><p>Choose a variable and run the suite to see which changes survive.</p></div></div></div></div>
    </div>
  </section>
  <!-- FORMULA REGISTRY -->
  <section class="page" id="page-formulas">
    <div class="sectionhead" style="margin-top:0"><div><h2>Calculation modes</h2><p>The same calculations the analysis tabs use, grouped by the question they answer: what is changing, how much, where, whether it is significant, and why. Every mode here is executable — pick one in <b>All calculations</b> to run it with your own inputs, or let the analysis tabs call it on a real series.</p></div><div class="telemetry" id="formula-count">—</div></div>
    <div class="searchrow"><input id="formula-search" type="text" placeholder="Search calculations by name or symbol…"><button class="btn" id="formula-refresh">Rescan</button></div>
    <div class="tagrow" id="formula-cats"></div>
    <div class="panel"><div class="flist" id="formula-list"><div class="empty">Loading registry…</div></div></div>
    <div class="note">The registry is generated from <span class="telemetry">src/*.hpp</span> at runtime. Use the Command Console to execute any command-line implementation; the Laboratories expose the retrieval and analysis paths interactively.</div>
  </section>

  <!-- COMMAND CONSOLE -->
  <section class="page" id="page-console">
    <div class="sectionhead" style="margin-top:0"><div><h2>Command console</h2><p>Every CLI formula and diagnostic, executed by the same binary.</p></div><div class="telemetry">STDIN → ENGINE</div></div>
    <div class="tool">
      <div class="panel"><div class="panelhead"><b>Transmit</b><span class="hint">climate &lt;command&gt;</span></div><div class="panelbody">
        <div class="field"><label>Command</label><div class="cmdbar"><input id="cmd-input" type="text" value="catalog --format json"><button class="btn primary" id="cmd-run">Run</button></div></div>
        <div class="field"><label>Batch · one command per line</label><textarea id="cmd-batch" rows="4" spellcheck="false">physics toa --forcing 3.7 --feedback -1.25 --delta-t 1.2
physics forcing-co2 --ppm 420 --reference-ppm 280
physics wind-stress --u10 10 --v10 0</textarea></div>
        <button class="btn" id="cmd-batch-run">Run batch →</button>
        <div class="quick">
          <button data-cmd="catalog --format json">catalog</button>
          <button data-cmd="physics mie --size-parameter 5">mie</button>
          <button data-cmd="analyze --file data/sst_monthly.csv --variable sst --units degC">trend</button>
          <button data-cmd="earth two-layer --years 150 --forcing 3.7">two-layer</button>
          <button data-cmd="ocean teos-pure --t 10 --p 1e7">teos-10</button>
          <button data-cmd="data list">nasa</button>
        </div>
        <div class="note">Commands run in the project directory. The <span class="telemetry">ui</span> command is blocked to prevent recursive servers. Output is capped at 200 kB.</div>
      </div></div>
      <div class="panel"><div class="panelhead"><b>Output</b><span class="hint">stdout + stderr</span></div><div class="panelbody"><div class="console" id="cmd-output">$ ready</div></div></div>
    </div>
    <div class="sectionhead"><div><h2>Command reference</h2><p>Complete help text, searchable.</p></div></div>
    <div class="searchrow"><input id="help-search" type="text" placeholder="Filter commands…"></div>
    <div class="panel"><div class="panelbody"><pre id="help-text">Loading help…</pre></div></div>
  </section>

  <!-- OVERVIEW -->
  <section class="page" id="page-overview">
    <div class="hero">
      <div>
        <div class="eyebrow">NASA climate methods / live satellite and model data</div>
        <h1>Every equation, <em>open to inspection.</em></h1>
        <p class="lede">A working bench for satellite retrievals, trend detection and climate diagnostics.
          No black boxes, no numerical libraries — the formulas from the notes, compiled and traceable.</p>
        <div class="actions">
          <button class="btn primary" data-nav="workflow">Open trend detective →</button>
          <button class="btn" data-nav="formulas">Browse all formulas</button>
        </div>
      </div>
      <div class="scope"><div class="scopehead"><span>SIGNAL MONITOR</span><span>DEMO SERIES</span></div><div id="hero-chart" style="flex:1;display:grid;place-items:center;color:#5a4330;font-size:9px">AWAITING DATASET</div><div class="scopefoot">CH 01 / MONTHLY / DEG C / NO SMOOTHING</div></div>
    </div>
    <div class="stats" id="overview-stats">
      <div class="stat"><div class="k">Catalog variables</div><div class="v">—</div><div class="s">satellite &amp; climate products</div></div>
      <div class="stat"><div class="k">Process groups</div><div class="v">—</div><div class="s">ocean, ice, land, atmosphere</div></div>
      <div class="stat"><div class="k">Demo datasets</div><div class="v">—</div><div class="s">ready to load locally</div></div>
      <div class="stat"><div class="k">Engine groups</div><div class="v">—</div><div class="s">retrieval, statistics, physics</div></div>
    </div>
    <div class="grid g4" id="module-grid"></div>

    <div class="sectionhead"><div><h2>Ask in plain language</h2><p>The deterministic intent compiler maps a question to catalog variables.</p></div></div>
    <div class="grid g2">
      <div class="card">
        <div class="field"><label for="query-text">Research question</label>
          <input id="query-text" type="text" value="is the arctic changing?" placeholder="e.g. how has sahel rainfall changed since 1980?"></div>
        <button class="btn primary run" id="query-run">Compile question →</button>
      </div>
      <div class="card" id="query-result"><div class="placeholder" style="height:auto;min-height:140px">
        <div><div class="big">⌘</div><p>Try “is the arctic changing?” or “how has sahel rainfall changed since 1980?”</p></div></div></div>
    </div>

    <div class="sectionhead"><div><h2>Designed for auditability</h2><p>What makes this different from a plotting notebook.</p></div></div>
    <div class="grid g3">
      <div class="card"><h3>Standard library only</h3><p class="desc" style="margin:0">No NetCDF, JSON, CSV or numerical library is linked. The parsers, SVD/least squares, FFT and PDE solvers are in this repository.</p></div>
      <div class="card"><h3>Stated formulas, not stubs</h3><p class="desc" style="margin:0">Every implemented method is exercised by the self-checking test binary, including closure checks and closed-form identities.</p></div>
      <div class="card"><h3>Corrections are documented</h3><p class="desc" style="margin:0">Dimensionally wrong notes formulas are replaced by their standard forms, with the change recorded in the README rather than hidden.</p></div>
    </div>
  </section>

  <!-- TIME SERIES -->
  <section class="page" id="page-timeseries">
    <div class="sectionhead" style="margin-top:0"><div><h2>Time series lab</h2><p>STL decomposition · Mann–Kendall · Sen’s slope with Hamed–Rao correction</p></div></div>
    <div class="tool">
      <div class="panel">
        <div class="panelhead"><b>Input series</b><span class="hint">CSV</span></div>
        <div class="panelbody">
          <div class="field"><label>Demo dataset</label><select id="ts-dataset"></select></div>
          <div class="field"><label>Or paste / upload CSV</label>
            <textarea id="ts-csv" rows="7" spellcheck="false"></textarea>
            <label class="filebtn" for="ts-file">Choose CSV file</label>
            <input type="file" id="ts-file" accept=".csv,text/csv" class="vhidden">
          </div>
          <div class="row">
            <div class="field"><label>Time column</label><select id="ts-time"></select></div>
            <div class="field"><label>Value column</label><select id="ts-value"></select></div>
          </div>
          <div class="row3">
            <div class="field"><label>Period</label><input id="ts-period" type="number" value="12" min="2" step="1"></div>
            <div class="field"><label>Alpha</label><input id="ts-alpha" type="number" value="0.05" min="0.001" max="0.5" step="0.01"></div>
            <div class="field"><label>Units</label><input id="ts-units" type="text" value=""></div>
          </div>
          <label class="check"><input id="ts-hamed" type="checkbox" checked> Hamed–Rao autocorrelation correction</label>
          <label class="check"><input id="ts-prewhiten" type="checkbox"> Trend-free prewhitening</label>
          <button class="btn primary run" id="ts-run">Analyze series →</button>
          <div class="note">Default demo series are synthetic and exist to exercise the pipeline. They are not NASA products.</div>
        </div>
      </div>
      <div id="ts-output"><div class="panel"><div class="placeholder"><div><div class="big">∿</div><p>Choose a dataset, then run the analysis. STL separates trend and seasonal components before the trend test.</p></div></div></div></div>
    </div>
  </section>

  <!-- CO2 -->
  <section class="page" id="page-co2">
    <div class="sectionhead" style="margin-top:0"><div><h2>CO₂ curve fit</h2><p>Polynomial + annual harmonics with iterative outlier rejection (Thoning et al.)</p></div></div>
    <div class="tool">
      <div class="panel">
        <div class="panelhead"><b>Input series</b><span class="hint">monthly</span></div>
        <div class="panelbody">
          <div class="field"><label>Demo dataset</label><select id="co2-dataset"></select></div>
          <div class="field"><label>Or paste / upload CSV</label>
            <textarea id="co2-csv" rows="7" spellcheck="false"></textarea>
            <label class="filebtn" for="co2-file">Choose CSV file</label>
            <input type="file" id="co2-file" accept=".csv,text/csv" class="vhidden">
          </div>
          <div class="row">
            <div class="field"><label>Time column</label><select id="co2-time"></select></div>
            <div class="field"><label>Value column</label><select id="co2-value"></select></div>
          </div>
          <div class="row3">
            <div class="field"><label>Polynomial</label><input id="co2-poly" type="number" value="2" min="1" max="5"></div>
            <div class="field"><label>Harmonics</label><input id="co2-harmonics" type="number" value="4" min="1" max="8"></div>
            <div class="field"><label>Units</label><input id="co2-units" type="text" value="ppm"></div>
          </div>
          <button class="btn primary run" id="co2-run">Fit curve →</button>
          <div class="note">The growth rate is the derivative of the polynomial trend, so the seasonal cycle does not contaminate it.</div>
        </div>
      </div>
      <div id="co2-output"><div class="panel"><div class="placeholder"><div><div class="big">CO₂</div><p>Fit a polynomial plus annual harmonics and read the growth rate at the ends of the record.</p></div></div></div></div>
    </div>
  </section>

  <!-- RETRIEVAL -->
  <section class="page" id="page-retrieval">
    <div class="sectionhead" style="margin-top:0"><div><h2>Retrieval lab</h2><p>Single-footprint and single-pixel calculations, using the same code as the CLI</p></div></div>
    <div class="tabs" id="retrieval-tabs">
      <button class="tab active" data-tab="sst">SST · NLSST</button>
      <button class="tab" data-tab="ndvi">Vegetation index</button>
      <button class="tab" data-tab="ice">Sea ice · NASA Team</button>
      <button class="tab" data-tab="precip">Precipitation · GPI</button>
    </div>
    <div class="tabpane active" id="tab-sst">
      <div class="split">
        <div class="panel"><div class="panelhead"><b>Brightness temperatures</b><span class="hint">Kelvin</span></div><div class="panelbody">
          <div class="row"><div class="field"><label>BT11</label><input id="sst-bt11" type="number" value="300.2" step="0.01"></div>
          <div class="field"><label>BT12</label><input id="sst-bt12" type="number" value="299.3" step="0.01"></div></div>
          <div class="row"><div class="field"><label>First guess</label><input id="sst-guess" type="number" value="300.0" step="0.01"></div>
          <div class="field"><label>Zenith angle °</label><input id="sst-zenith" type="number" value="35" step="0.1"></div></div>
          <div class="row"><div class="field"><label>Mirror side</label><input id="sst-mirror" type="number" value="0" min="0" max="1" step="1"></div>
          <div class="field"><label>Theta star °</label><input id="sst-theta" type="number" value="0" step="0.1"></div></div>
          <button class="btn primary run" id="sst-run">Retrieve SST →</button>
          <div class="formula">SST = a₀ + a₁BT₁₁ + a₂(BT₁₁−BT₁₂)Tₛ + a₃(secθ−1)(BT₁₁−BT₁₂) + a₄mirror + a₅θ* + a₆θ²</div>
        </div></div>
        <div id="sst-output"><div class="panel"><div class="placeholder"><div><div class="big">°C</div><p>The model is calibrated from the bundled matchup file at request time.</p></div></div></div></div>
      </div>
    </div>
    <div class="tabpane" id="tab-ndvi">
      <div class="split">
        <div class="panel"><div class="panelhead"><b>Reflectances</b><span class="hint">0–1</span></div><div class="panelbody">
          <div class="field"><label>Index</label><select id="vi-index"><option value="ndvi">NDVI</option><option value="evi">EVI</option><option value="savi">SAVI</option></select></div>
          <div class="row"><div class="field"><label>Red</label><input id="vi-red" type="number" value="0.12" step="0.001"></div>
          <div class="field"><label>NIR</label><input id="vi-nir" type="number" value="0.48" step="0.001"></div></div>
          <div class="row" id="vi-extra"><div class="field"><label>Blue (EVI)</label><input id="vi-blue" type="number" value="0.08" step="0.001"></div>
          <div class="field"><label>Soil factor L</label><input id="vi-soil" type="number" value="0.5" step="0.05"></div></div>
          <button class="btn primary run" id="vi-run">Compute index →</button>
          <div class="formula" id="vi-formula">NDVI = (ρNIR − ρRED) / (ρNIR + ρRED)</div>
        </div></div>
        <div id="vi-output"><div class="panel"><div class="placeholder"><div><div class="big">NDVI</div><p>An undefined ratio (both reflectances zero) is reported as non-finite rather than silently zero.</p></div></div></div></div>
      </div>
    </div>
    <div class="tabpane" id="tab-ice">
      <div class="split">
        <div class="panel"><div class="panelhead"><b>Passive-microwave footprint</b><span class="hint">Kelvin</span></div><div class="panelbody">
          <div class="row3"><div class="field"><label>19V</label><input id="ice-19v" type="number" value="240" step="0.1"></div>
          <div class="field"><label>19H</label><input id="ice-19h" type="number" value="228" step="0.1"></div>
          <div class="field"><label>37V</label><input id="ice-37v" type="number" value="242" step="0.1"></div></div>
          <button class="btn primary run" id="ice-run">Invert tie points →</button>
          <div class="formula">PR = (TB₁₉ⱽ−TB₁₉ᴴ)/(TB₁₉ⱽ+TB₁₉ᴴ) · GR = (TB₃₇ⱽ−TB₁₉ⱽ)/(TB₃₇ⱽ+TB₁₉ⱽ)</div>
          <div class="note">Uses the three-constituent physical mixture plus the calibrated polarization/gradient retrieval.</div>
        </div></div>
        <div id="ice-output"><div class="panel"><div class="placeholder"><div><div class="big">❄</div><p>The 15% extent convention is applied to the retrieved total concentration.</p></div></div></div></div>
      </div>
    </div>
    <div class="tabpane" id="tab-precip">
      <div class="split">
        <div class="panel"><div class="panelhead"><b>GPI inputs</b><span class="hint">IR threshold 235 K</span></div><div class="panelbody">
          <div class="row"><div class="field"><label>Cold-cloud fraction</label><input id="pr-cold" type="number" value="0.35" step="0.01" min="0" max="1"></div>
          <div class="field"><label>Precipitable water cm</label><input id="pr-wv" type="number" value="4.2" step="0.1"></div></div>
          <div class="row"><div class="field"><label>a₁</label><input id="pr-a1" type="number" value="0.1" step="0.01"></div>
          <div class="field"><label>a₂</label><input id="pr-a2" type="number" value="3.0" step="0.1"></div></div>
          <button class="btn primary run" id="pr-run">Compute GPI →</button>
          <div class="formula">P = 3 mm/h × f(cold) × exp[a₁(w − a₂)]</div>
        </div></div>
        <div id="pr-output"><div class="panel"><div class="placeholder"><div><div class="big">☂</div><p>Returns both the raw 3 mm/h cold-cloud estimate and the water-vapour adjusted value.</p></div></div></div></div>
      </div>
    </div>
  </section>

  <!-- INFO -->
  <section class="page" id="page-info">
    <div class="sectionhead" style="margin-top:0"><div><h2>Information theory</h2><p>Histogram MI, Miller–Madow correction and the KSG k-NN estimator side by side</p></div></div>
    <div class="tool">
      <div class="panel"><div class="panelhead"><b>Samples</b><span class="hint">comma-separated</span></div><div class="panelbody">
        <div class="field"><label>X</label><textarea id="mi-x" rows="4">1,2,3,4,5,6,7,8,9,10,11,12</textarea></div>
        <div class="field"><label>Y</label><textarea id="mi-y" rows="4">2,4,5,8,9,12,13,16,17,20,21,24</textarea></div>
        <div class="row"><div class="field"><label>Histogram bins</label><input id="mi-bins" type="number" value="6" min="2" max="64"></div>
        <div class="field"><label>KSG neighbours</label><input id="mi-k" type="number" value="3" min="1" max="20"></div></div>
        <button class="btn primary run" id="mi-run">Estimate MI →</button>
        <div class="formula">I(X;Y) = Σ p(x,y) log₂[p(x,y)/(p(x)p(y))]</div>
        <div class="note">The Miller–Madow term is (K−1)(L−1)/(2n ln 2). KSG uses the max-norm neighbourhood.</div>
      </div></div>
      <div id="mi-output"><div class="panel"><div class="placeholder"><div><div class="big">I</div><p>Mutual information is reported in bits. Use monotonically related series to see the estimators diverge on small samples.</p></div></div></div></div>
    </div>
  </section>

  <!-- UNCERTAINTY -->
  <section class="page" id="page-uncertainty">
    <div class="sectionhead" style="margin-top:0"><div><h2>SST uncertainty</h2><p>HadSST4-style uncorrelated, correlated and ensemble components</p></div></div>
    <div class="tool">
      <div class="panel"><div class="panelhead"><b>Error model</b><span class="hint">weights × errors</span></div><div class="panelbody">
        <div class="field"><label>Grid-cell weights</label><textarea id="un-w" rows="3">0.25,0.25,0.25,0.25</textarea></div>
        <div class="field"><label>Uncorrelated 1σ errors</label><textarea id="un-s" rows="3">0.1,0.2,0.15,0.25</textarea></div>
        <div class="row"><div class="field"><label>Correlation ρ</label><input id="un-rho" type="number" value="0.5" step="0.05" min="0" max="0.99"></div>
        <div class="field"><label>Ensemble values (optional)</label><input id="un-ens" type="text" value="0.1,0.2,0.15,0.25"></div></div>
        <button class="btn primary run" id="un-run">Propagate errors →</button>
        <div class="formula">σ²ᵤ = Σwᵢ²σᵢ² · σ²ₛ = σ²[ρ(Σwᵢ)² + (1−ρ)Σwᵢ²]</div>
        <div class="note">The correlated component uses the constant-covariance limit. Effective sample size is n·Σw²/[ρ(Σw)²+(1−ρ)Σw²].</div>
      </div></div>
      <div id="un-output"><div class="panel"><div class="placeholder"><div><div class="big">σ</div><p>Components are combined in quadrature. The effective sample size tells you how much independent information the weighted mean carries.</p></div></div></div></div>
    </div>
  </section>

  <!-- LIVE NASA DATA -->
  <section class="page" id="page-live">
    <div class="sectionhead" style="margin-top:0"><div><h2>Fetch data</h2><p>Step 1. Every series here is pulled from a live NASA endpoint at the moment you press the button: GISTEMP, NASA POWER, the NSIDC Sea Ice Index, and any Earthdata dataset you can log in for. Each pull records its source URL and retrieval time in <code>data/live_*.json</code>, and the pulled series then behaves like any other variable in the workbench. The NASA product planner below returns a dry-run manifest and the exact fetch command before anything is downloaded.</p></div><div class="telemetry" id="live-creds">checking…</div></div>
    <div class="panel"><div class="panelhead"><b>Pull parameters</b><span class="hint">applies to POWER and Earthdata</span></div><div class="panelbody">
      <div class="row3">
        <div class="field"><label>Latitude</label><input id="live-lat" type="number" step="0.5" value="72"></div>
        <div class="field"><label>Longitude</label><input id="live-lon" type="number" step="0.5" value="-40"></div>
        <div class="field"><label>POWER parameter</label><input id="live-param" type="text" value="T2M"></div>
      </div>
      <div class="row3">
        <div class="field"><label>Start</label><input id="live-start" type="text" value="2000-01-01"></div>
        <div class="field"><label>End</label><input id="live-end" type="text" value="2024-12-31"></div>
        <div class="field"><label>Sea-ice region</label><select id="live-region"><option value="N">North</option><option value="S">South</option></select></div>
      </div>
      <div class="row"><div class="field"><label>Earthdata short name</label><input id="live-short" type="text" value="MODISA_L3m_SST"></div><button class="btn" id="live-reload">Refresh status</button></div>
    </div></div>
    <div class="panel"><div class="panelhead"><b>Earthdata Login</b><span class="hint">optional · used once · never stored</span></div><div class="panelbody">
      <div class="row3">
        <div class="field"><label>Earthdata username</label><input id="ea-user" type="text" autocomplete="off" placeholder="your Earthdata login"></div>
        <div class="field"><label>Password</label><input id="ea-pass" type="password" autocomplete="off" placeholder="password"></div>
        <div class="field"><label>Or bearer token</label><input id="ea-token" type="password" autocomplete="off" placeholder="token instead of a password"></div>
      </div>
      <div class="row" style="margin-top:8px"><button class="btn" id="ea-check">Verify Earthdata login</button><span class="telemetry">nothing is downloaded by this check</span></div>
      <div id="ea-status"></div>
      <div class="note">Used only for the Earthdata pull you press: the browser keeps them in this tab's memory, the server passes them to the fetch helper as environment variables for that one child process, and nothing is written to disk, logs or the page. Leave blank to fall back to credentials already exported in your shell or a <span class="telemetry">~/.netrc</span> file. Open endpoints (GISTEMP, POWER, NSIDC sea ice) never need a login. A bearer token is checked directly here, with no Python and no extra package; a password check needs the optional Python helper, which the portable build does not ship. Earthdata accounts with multi-factor authentication cannot use a password at all: paste a bearer token.</div>
    </div></div>
    <div id="live-list" class="live-grid"><div class="hint">Loading sources…</div></div>
    <div id="live-output"></div>
  <div class="fetch-extra">
    <div class="sectionhead" style="margin-top:0"><div><h2>NASA Earthdata uplink</h2><p>Product registry, live CMR metadata, and API-backed acquisition routes for every variable class.</p></div><div class="telemetry">EARTHDATA / CMR / DAAC</div></div>
    <div class="statusline"><span>CMR API <b>LIVE / NO DOWNLOAD</b></span><span>DATA ACCESS <b>EARTHDATA LOGIN</b></span></div>
    <div class="dsgrid" id="nasa-ds"><div class="empty">Loading NASA product registry…</div></div>
    <div class="sectionhead"><div><h2>Fetch planner</h2><p>Dry run: validates the request without downloading or requiring credentials.</p></div></div>
    <div class="tool">
      <div class="panel"><div class="panelhead"><b>Request</b><span class="hint">CMR</span></div><div class="panelbody">
        <div class="field"><label>Dataset key</label><select id="nasa-key"></select></div>
        <div class="row"><div class="field"><label>Start</label><input id="nasa-start" type="text" value="2024-01-01"></div><div class="field"><label>End</label><input id="nasa-end" type="text" value="2024-12-31"></div></div>
        <div class="field"><label>Maximum granules</label><input id="nasa-max" type="number" value="4" min="1" max="100"></div>
        <div class="row"><button class="btn primary run" id="nasa-plan">Generate fetch plan →</button><button class="btn" id="nasa-search">Query CMR API →</button></div>
        <div class="note">CMR search is a live API call and downloads no data. Live arrays need a free NASA Earthdata Login, <span class="telemetry">earthaccess</span>, <span class="telemetry">xarray</span> and <span class="telemetry">netCDF4</span>. No credentials are stored by this application.</div>
      </div></div>
      <div class="panel"><div class="panelhead"><b>Plan / command output</b><span class="hint">JSON</span></div><div class="panelbody"><pre id="nasa-output">Generate a plan to see the exact fetch command and dry-run manifest.</pre></div></div>
    </div>
    <div class="sectionhead"><div><h2>Bundled NASA-pattern fixtures</h2><p>Small local files used by the readers and tests. Not freshly downloaded products.</p></div></div>
    <div class="panel"><div class="panelbody"><div id="nasa-files" class="tagrow"></div></div></div>
  </div>
  </section>


    <div class="sectionhead"><div><h2>How the retrieved variable changes over time</h2><p>A retrieval above gives one footprint. This runs the change analysis on the full record of the same quantity, so the number you retrieved and the change you report are the same variable.</p></div><div class="telemetry">SERIES → TREND</div></div>
    <div class="tool">
      <div class="panel"><div class="panelhead"><b>Record</b><span class="hint">real values from data/</span></div><div class="panelbody">
        <div class="field"><label>Variable</label><select id="ret-series">
          <option value="sst_monthly.csv|sst|degC">Sea surface temperature · NLSST record</option>
          <option value="ndvi_monthly.csv|value|unitless">Vegetation index · NDVI record</option>
          <option value="sea_ice_extent_monthly.csv|value|million km2">Sea ice extent · NASA Team record</option>
          <option value="precipitation_monthly.csv|value|mm/day">Precipitation · GPI record</option>
          <option value="live_gistemp.csv|value|degC anomaly">GISTEMP global anomaly · live pull</option>
          <option value="live_seaice_n.csv|value|million km2">Arctic sea ice · live pull</option>
        </select></div>
        <div class="field"><label>Alpha</label><input id="ret-alpha" type="number" value="0.05" step="0.01" min="0.001" max="0.2"></div>
        <button class="btn primary run" id="ret-run">Plot the record →</button>
        <div class="note">The chart shows the raw series, the STL trend and the Sen slope. The verdict beside it is the Mann–Kendall test on the deseasonalised series, so the same number is defended in the Significance tab.</div>
      </div></div>
      <div id="ret-out"><div class="panel"><div class="placeholder"><div><div class="big">∿</div><p>Pick a record to see the retrieval and its change together.</p></div></div></div></div>
    </div>
  </section>

  <!-- PHYSICS -->
  <section class="page" id="page-physics">
    <div class="sectionhead" style="margin-top:0"><div><h2>Formula laboratory</h2><p>Every public physics formula, selectable in any combination. Energy balance, radiative forcing, Mie scattering, wind-driven transport, geostrophy, vegetation indices.</p></div></div>
    <div class="tool">
      <div class="panel"><div class="panelhead"><b>Calculator</b><span class="hint">SI unless marked</span></div><div class="panelbody">
        <div class="field"><label>Formulas · select one or many</label><div class="ph-picks" id="ph-picks"></div></div>
        <div class="row"><button class="btn" id="ph-all">Select all</button><button class="btn" id="ph-clear">Clear</button></div>
        <div id="ph-fields"></div>
        <button class="btn primary run" id="ph-run">Calculate selected →</button>
        <div class="formula" id="ph-formula-desc">Each checked formula gets its own inputs and its own result. Batch runs share no hidden state.</div>
      </div></div>
      <div id="ph-output"><div class="panel"><div class="placeholder"><div><div class="big">ƒ</div><p>Pick a formula. The fields and the equation update with the selection.</p></div></div></div></div>
    </div>
  </section>

  <!-- COMPLETE FORMULA LAB -->
  <section class="page" id="page-alllab">
    <div class="sectionhead" style="margin-top:0"><div><h2>Complete formula lab</h2><p>Every public function in the engine — <span id="lab-count">—</span> of them — generated directly from the headers and executable from here. Select any number, set their inputs, and run them in one batch.</p></div><div class="telemetry" id="lab-selected">0 selected</div></div>
    <div class="tool">
      <div class="panel"><div class="panelhead"><b>Formulas</b><span class="hint">grouped in workflow order</span></div><div class="panelbody">
        <div class="field"><label>Series source for series-taking formulas</label><select id="lab-file"></select></div>
        <div class="row"><div class="field"><label>Matrix columns</label><input id="lab-cols" type="number" value="3" min="1" step="1"></div><div class="field"><label>Text argument</label><input id="lab-text" type="text" value="keyword"></div></div>
        <div class="field"><label>Filter</label><input id="lab-search" type="text" placeholder="name, module, return type…"></div>
        <div class="row"><button class="btn" id="lab-select-visible">Select visible</button><button class="btn" id="lab-clear">Clear</button></div>
        <div id="lab-list" style="margin-top:7px"><div class="hint">Loading formula bindings…</div></div>
      </div></div>
      <div>
        <div class="panel"><div class="panelhead"><b>Inputs</b><span class="hint">one card per selected formula</span></div><div class="panelbody" id="lab-inputs"><div class="hint">Nothing selected yet.</div></div></div>
        <button class="btn primary run" id="lab-run">Run selected →</button>
        <div class="panel" style="margin-top:10px"><div class="panelhead"><b>Results</b><span class="hint" id="lab-summary">—</span></div><div class="panelbody lab-results" id="lab-output"><div class="hint">Each formula returns its own value; failures are reported per formula and do not stop the batch.</div></div></div>
      </div>
    </div>
  </section>

  <!-- VISUALIZE EVERYTHING -->
  <section class="page" id="page-visual">
    <div class="sectionhead" style="margin-top:0"><div><h2>Visualize change over time</h2><p>Every series the workbench can reach, drawn from its real values: the raw record, the fitted trend and the sign of the change. Regional opposition shows up as two curves with different slopes, never as one averaged line.</p></div><div class="row"><span class="telemetry" id="vis-count">—</span><button class="btn" id="vis-reload">Reload</button></div></div>
    <div class="panel" style="margin-top:10px"><div class="panelhead"><b>Spatial trend field</b><span class="hint">synthetic 12 × 24 demonstration of per-cell Sen slope and Gi* hotspots — not a downloaded product</span></div><div class="panelbody" id="vis-map"><div class="hint">Loading field…</div></div></div>
    <div class="vis-grid" id="vis-grid"><div class="hint">Loading series…</div></div>
  </section>

  <!-- CATALOG -->
  <section class="page" id="page-catalog">
    <div class="sectionhead" style="margin-top:0"><div><h2>Select variables</h2><p>Step 2. A measured variable is only worth analysing if its record can resolve the change you are looking for, so every card states <b>what kind of change it measures</b>. Pick one and it is carried straight into the analysis and significance tabs.</p></div><div class="telemetry">WHAT KIND OF CHANGE</div></div>
    <div class="searchrow">
      <input id="cat-search" type="text" placeholder="Search name, retrieval, description…">
      <select id="cat-process"><option value="">All processes</option></select>
      <select id="cat-region"><option value="">All regions</option></select>
    </div>
    <div class="panel"><div id="cat-list"><div class="placeholder"><div><p>Loading catalog…</p></div></div></div></div>
  </section>
</main>
</div>
<div class="toast" id="toast"></div>
<div class="coach" id="coach" hidden>
  <div class="coach-arrow" id="coach-arrow"></div>
  <div class="coach-body">
    <div class="coach-step" id="coach-step">1 / 5</div>
    <div class="coach-text" id="coach-text"></div>
    <div class="coach-actions"><button class="btn" id="coach-skip">Skip the tour</button><button class="btn primary" id="coach-next">Next</button></div>
  </div>
</div>
<script>
const $ = s => document.querySelector(s);
const $$ = s => Array.from(document.querySelectorAll(s));
const state = {overview:null, catalog:[], physics:null, workflow:null, nasa:null, formulas:null, help:""};

function esc(s){return String(s ?? '').replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));}
function fmt(v,d=3){if(v===null||v===undefined||!isFinite(v))return '—';
  const a=Math.abs(v); if(a!==0&&(a<1e-3||a>=1e6))return v.toExponential(2);
  return Number(v).toFixed(d).replace(/(\.\d*?)0+$/,'$1').replace(/\.$/,'');}
function pval(v){if(!isFinite(v))return '—'; return v<0.0001?'< 0.0001':v.toFixed(4);}
function toast(msg){const t=$('#toast');t.textContent=msg;t.classList.add('show');clearTimeout(t._h);t._h=setTimeout(()=>t.classList.remove('show'),6000);}
async function api(path,body){
  const r=await fetch(path,{method:body?'POST':'GET',headers:body?{'Content-Type':'application/json'}:{},body:body?JSON.stringify(body):undefined});
  const text=await r.text(); let data;
  try{data=JSON.parse(text);}catch(e){throw new Error(text.slice(0,300)||'Invalid server response');}
  if(!r.ok){
    const detail = typeof data.error==='string' ? data.error : JSON.stringify(data.error||data);
    throw new Error(detail||('HTTP '+r.status));
  }
  return data;
}
function nav(name){
  $$('.page').forEach(p=>p.classList.toggle('active',p.id==='page-'+name));
  $$('.nav').forEach(b=>b.classList.toggle('active',b.dataset.nav===name));
  const label={workflow:'What / where / how much',significance:'Significance',formulas:'Calculation modes',console:'Command console',overview:'System',timeseries:'Time series lab',co2:'CO₂ curve fit',retrieval:'Retrievals',info:'Information theory',uncertainty:'Uncertainty',physics:'Formula calculator',alllab:'All calculations',visual:'Visualize over time',live:'Fetch data',catalog:'Select variables'}[name]||name;
  $('#crumb').textContent=label; window.scrollTo({top:0,behavior:'smooth'});
}
document.addEventListener('click',e=>{const b=e.target.closest('[data-nav]');if(b){nav(b.dataset.nav);return;}
  const g=e.target.closest('[data-goto]');if(g){nav(g.dataset.goto);window.scrollTo({top:0,behavior:'smooth'});}});
$$('#retrieval-tabs .tab').forEach(t=>t.addEventListener('click',()=>{
  $$('#retrieval-tabs .tab').forEach(x=>x.classList.remove('active'));t.classList.add('active');
  $$('.tabpane').forEach(p=>p.classList.remove('active'));$('#tab-'+t.dataset.tab).classList.add('active');
}));

function metric(k,v,s='',tone=''){return `<div class="metric"><div class="k">${esc(k)}</div><div class="v ${tone}">${v}</div>${s?`<div class="s">${esc(s)}</div>`:''}</div>`;}
function chart(series,opts={}){
  const W=880,H=opts.height||290,L=54,R=14,T=26,B=26;
  let vals=[];series.forEach(s=>s.values.forEach(v=>{if(v!==null&&isFinite(v))vals.push(v);}));
  if(!vals.length)return '<div class="empty">No finite values to plot.</div>';
  let min=Math.min(...vals),max=Math.max(...vals);
  if(min===max){min-=1;max+=1;}
  const pad=(max-min)*0.08;min-=pad;max+=pad;
  const n=Math.max(...series.map(s=>s.values.length));
  const x=i=>L+(n<=1?0:i*(W-L-R)/(n-1));
  const y=v=>T+(max-v)*(H-T-B)/(max-min);
  let g='';
  for(let i=0;i<=4;i++){const yy=T+i*(H-T-B)/4;const val=max-(max-min)*i/4;
    g+=`<line class="gridline" x1="${L}" y1="${yy}" x2="${W-R}" y2="${yy}"/><text x="${L-7}" y="${yy+3}" text-anchor="end">${fmt(val,2)}</text>`;}
  g+=`<line class="axis" x1="${L}" y1="${H-B}" x2="${W-R}" y2="${H-B}"/>`;
  const xl=opts.xlabels||[];
  if(xl.length>=2){for(const p of [0,Math.floor((xl.length-1)/2),xl.length-1]){
    g+=`<text x="${x(p)}" y="${H-8}" text-anchor="${p===0?'start':(p===xl.length-1?'end':'middle')}">${esc(xl[p])}</text>`;}}
  series.forEach(s=>{
    let points=[],segments=[];
    s.values.forEach((v,i)=>{if(v===null||!isFinite(v)){if(points.length)segments.push(points);points=[];return;}points.push({x:x(i),y:y(v)});});
    if(points.length)segments.push(points);
    let d='';
    segments.forEach(seg=>{
      if(seg.length===1){d+=`M${seg[0].x.toFixed(1)} ${seg[0].y.toFixed(1)} `;return;}
      d+=`M${seg[0].x.toFixed(1)} ${seg[0].y.toFixed(1)} `;
      const deltas=[],slopes=[];
      for(let i=0;i<seg.length-1;i++)deltas.push((seg[i+1].y-seg[i].y)/(seg[i+1].x-seg[i].x));
      for(let i=0;i<seg.length;i++){
        if(i===0)slopes.push(deltas[0]);
        else if(i===seg.length-1)slopes.push(deltas[deltas.length-1]);
        else if(deltas[i-1]*deltas[i]<=0)slopes.push(0);
        else slopes.push((deltas[i-1]+deltas[i])/2);
      }
      for(let i=0;i<seg.length-1;i++){
        const dx=(seg[i+1].x-seg[i].x)/3;
        const m1=Math.min(Math.abs(slopes[i]),3*Math.abs(deltas[i]))*Math.sign(slopes[i]);
        const m2=Math.min(Math.abs(slopes[i+1]),3*Math.abs(deltas[i]))*Math.sign(slopes[i+1]);
        d+=`C${(seg[i].x+dx).toFixed(1)} ${(seg[i].y+dx*m1).toFixed(1)} ${(seg[i+1].x-dx).toFixed(1)} ${(seg[i+1].y-dx*m2).toFixed(1)} ${seg[i+1].x.toFixed(1)} ${seg[i+1].y.toFixed(1)} `;
      }
    });
    g+=`<path class="line" d="${d}" stroke="${s.color||'#f7e8c8'}"/>`;
  });
  let legend='';let lx=L;
  series.forEach(s=>{legend+=`<line x1="${lx}" y1="8" x2="${lx+10}" y2="8" stroke="${s.color||'#f7e8c8'}" stroke-width="2"/><text class="legend" x="${lx+14}" y="11" fill="${s.color||'#f7e8c8'}">${esc(s.name)}</text>`;lx+=s.name.length*6+34;});
  return `<svg class="chart" viewBox="0 0 ${W} ${H}" role="img"><g>${g}</g><g>${legend}</g></svg>`;
}
function chartCard(title,series,opts={}){return `<div class="chartwrap"><div style="font-size:9.5px;color:#d8c2a1;padding:2px 8px 0;text-transform:uppercase;letter-spacing:.08em">${esc(title)}</div>${chart(series,opts)}</div>`;}

function loadOverview(){
  api('/api/overview').then(d=>{
    state.overview=d;
    $('#overview-stats').innerHTML=`
      <div class="stat"><div class="k">Catalog variables</div><div class="v">${d.catalog_count}</div><div class="s">satellite &amp; climate products</div></div>
      <div class="stat"><div class="k">NASA product routes</div><div class="v" id="stat-nasa">—</div><div class="s">CMR / DAAC registry</div></div>
      <div class="stat"><div class="k">Public C++ functions</div><div class="v" id="stat-functions">—</div><div class="s">indexed from headers</div></div>
      <div class="stat"><div class="k">Demo datasets</div><div class="v">${d.files.length}</div><div class="s">local synthetic fixtures</div></div>`;
    $('#foot-catalog').textContent=d.catalog_count;
    $('#module-grid').innerHTML=d.modules.map(m=>`<div class="card module" data-nav="${m.target}">
      <div class="tag">${esc(m.code)}</div><h3>${esc(m.name)}</h3><p class="desc" style="margin:0">${esc(m.desc)}</p><div class="arrow">↗</div></div>`).join('');
    const ts=$('#ts-dataset'),co=$('#co2-dataset');
    ts.innerHTML=d.files.map(f=>`<option value="${esc(f)}">data/${esc(f)}</option>`).join('');
    co.innerHTML=d.files.map(f=>`<option value="${esc(f)}">data/${esc(f)}</option>`).join('');
    const lf=$('#lab-file');
    if(lf){lf.innerHTML=d.files.map(f=>`<option value="${esc(f)}">data/${esc(f)}</option>`).join('');
      if(d.files.includes('sst_monthly.csv'))lf.value='sst_monthly.csv';}
    const sst=d.files.includes('sst_monthly.csv')?'sst_monthly.csv':d.files[0];
    const c=d.files.includes('co2_monthly.csv')?'co2_monthly.csv':d.files[0];
    if(sst){ts.value=sst;loadCsv('ts');}
    if(c){co.value=c;loadCsv('co2');}
    fillSelects();
    loadCatalog();
  }).catch(e=>toast('Could not load overview: '+e.message));
}
function fillSelects(){
  if(!state.overview)return;
  $('#cat-process').innerHTML='<option value="">All processes</option>'+state.overview.processes.map(p=>`<option>${esc(p)}</option>`).join('');
  $('#cat-region').innerHTML='<option value="">All regions</option>'+state.overview.regions.map(p=>`<option>${esc(p)}</option>`).join('');
}
function columns(text,timeSel,valSel,prefTime,prefVal){
  const line=text.split(/\r?\n/).find(l=>l.trim());
  if(!line)return;
  const hs=line.split(',').map(s=>s.trim().replace(/^"|"$/g,'')).filter(Boolean);
  const opts=hs.map(h=>`<option value="${esc(h)}">${esc(h)}</option>`).join('');
  timeSel.innerHTML=opts;valSel.innerHTML=opts;
  const ti=hs.findIndex(h=>h.toLowerCase()===prefTime); if(ti>=0)timeSel.selectedIndex=ti;
  const vi=hs.findIndex(h=>h.toLowerCase()===prefVal); if(vi>=0)valSel.selectedIndex=vi;
  else if(hs.length>1)valSel.selectedIndex=1;
}
async function loadCsv(prefix){
  const name=$('#'+prefix+'-dataset').value;
  try{
    const text=await (await fetch('/api/data/file?name='+encodeURIComponent(name))).text();
    $('#'+prefix+'-csv').value=text;
    columns(text,$('#'+prefix+'-time'),$('#'+prefix+'-value'),
      prefix==='ts'?'time':'time',prefix==='ts'?'sst':'co2');
    if(prefix==='ts')drawHeroChart(text);
  }catch(e){toast(e.message);}
}
$('#ts-dataset').addEventListener('change',()=>loadCsv('ts'));
$('#co2-dataset').addEventListener('change',()=>loadCsv('co2'));
$('#ts-csv').addEventListener('input',()=>columns($('#ts-csv').value,$('#ts-time'),$('#ts-value'),'time',''));
$('#co2-csv').addEventListener('input',()=>columns($('#co2-csv').value,$('#co2-time'),$('#co2-value'),'time','co2'));
function wireFile(input,textarea,prefix){
  $(input).addEventListener('change',e=>{const f=e.target.files[0];if(!f)return;const r=new FileReader();
    r.onload=()=>{$(textarea).value=r.result;columns(r.result,$('#'+prefix+'-time'),$('#'+prefix+'-value'),'time','');};r.readAsText(f);});
}
wireFile('#ts-file','#ts-csv','ts');wireFile('#co2-file','#co2-csv','co2');

function mkBlock(mk){
  if(!mk)return '<div class="empty">No test result.</div>';
  return `${metric('Mann–Kendall τ',fmt(mk.tau,4),mk.direction,mk.significant?'good':'')}
    ${metric('p-value',pval(mk.p_value),mk.significant?'significant':'not significant',mk.significant?'good':'warn')}
    ${metric('S statistic',fmt(mk.s,0),'variance '+fmt(mk.variance,1))}
    ${metric('Correction',fmt(mk.correction_factor,3)+'×',mk.variance_corrected?'Hamed–Rao':'tie-corrected')}`;
}
function renderAnalysis(d){
  const mk=d.trend.mann_kendall,sen=d.trend.sen;
  const verdict=mk.significant
    ? `<b>Significant ${esc(mk.direction)}</b> at α = ${d.alpha}. Sen slope ${fmt(sen.slope,4)} ${esc(d.units)} per time step (${fmt(d.trend.sen_per_decade,3)} per decade).`
    : `<b>No significant monotonic trend</b> at α = ${d.alpha}. The estimated Sen slope is ${fmt(sen.slope,4)} ${esc(d.units)} per time step.`;
  const ci=sen.computed?`95% CI ${fmt(sen.lower,3)} … ${fmt(sen.upper,3)}`:'CI unavailable';
  const labels=d.series.time.map(t=>fmt(t,1));
  $('#ts-output').innerHTML=`<div class="panel"><div class="panelhead"><b>${esc(d.variable)}</b><span class="hint">${d.samples} samples · ${d.decomposed?'STL':'climatology'}</span></div><div class="panelbody">
    <div class="verdict">${verdict}</div>
    <div class="cards">
      ${metric('Sen slope / decade',fmt(d.trend.sen_per_decade,3),ci,'good')}
      ${metric('MK p-value',pval(mk.p_value),mk.direction,mk.significant?'good':'warn')}
      ${metric('Seasonal strength',d.decomposition.seasonal_strength?fmt(d.decomposition.seasonal_strength*100,1)+'%':'—','variance explained')}
      ${metric('Trend strength',d.decomposition.trend_strength?fmt(d.decomposition.trend_strength*100,1)+'%':'—','of variance')}
    </div>
    ${chartCard('Observed series and STL trend',[{name:'observed',values:d.series.observed,color:'#f7e8c8'},{name:'trend',values:d.components.trend,color:'#ffab54'}],{xlabels:labels,height:300})}
    <div class="minicharts">
      ${chartCard('Seasonal component',[{name:'seasonal',values:d.components.seasonal,color:'#45d6c9'}],{xlabels:labels,height:210})}
      ${chartCard('Remainder',[{name:'remainder',values:d.components.remainder,color:'#ff8f86'}],{xlabels:labels,height:210})}
    </div>
    <div class="cards" style="margin-top:12px">${mkBlock(mk)}</div>
    <div class="kv">
      <span class="k">Record</span><span class="v">${fmt(d.start_time,3)} → ${fmt(d.end_time,3)}</span>
      <span class="k">Spacing</span><span class="v">${fmt(d.spacing,4)}</span>
      <span class="k">Mean ± SD</span><span class="v">${fmt(d.summary.mean,4)} ± ${fmt(d.summary.sd,4)} ${esc(d.units)}</span>
      <span class="k">Range</span><span class="v">${fmt(d.summary.minimum,4)} … ${fmt(d.summary.maximum,4)}</span>
      <span class="k">Missing filled</span><span class="v">${d.missing_filled}</span>
    </div>
    ${d.warnings.length?`<ul class="warnlist">${d.warnings.map(w=>`<li>${esc(w)}</li>`).join('')}</ul>`:''}
    <details><summary>Raw JSON</summary><pre>${esc(JSON.stringify(d,null,2))}</pre></details>
  </div></div>`;
}
$('#ts-run').addEventListener('click',async e=>{
  const b=e.target;b.disabled=true;b.innerHTML='<span class="spinner"></span>Analyzing…';
  try{
    const d=await api('/api/analyze',{csv:$('#ts-csv').value,time:$('#ts-time').value,value:$('#ts-value').value,
      period:+$('#ts-period').value,alpha:+$('#ts-alpha').value,units:$('#ts-units').value,
      hamed_rao:$('#ts-hamed').checked,prewhiten:$('#ts-prewhiten').checked});
    renderAnalysis(d);
  }catch(err){toast(err.message);$('#ts-output').innerHTML=`<div class="panel"><div class="empty">${esc(err.message)}</div></div>`;}
  b.disabled=false;b.textContent='Analyze series →';
});

function renderCo2(d){
  const g=d.growth;
  $('#co2-output').innerHTML=`<div class="panel"><div class="panelhead"><b>${esc(d.variable)}</b><span class="hint">${d.samples} samples</span></div><div class="panelbody">
    <div class="cards">
      ${metric('Growth at end',fmt(g.end_per_year,3),esc(d.units)+' / year','good')}
      ${metric('Mean growth',fmt(g.mean_per_year,3),esc(d.units)+' / year')}
      ${metric('R²',fmt(d.model.r_squared,5),'residual SD '+fmt(d.model.residual_sd,4))}
      ${metric('Seasonal amplitude',fmt(d.seasonal_amplitude,3),esc(d.units))}
    </div>
    ${chartCard('Observed and fitted curve',[{name:'observed',values:d.series.observed,color:'#f7e8c8'},{name:'fitted',values:d.series.fitted,color:'#ffab54'}],{xlabels:d.series.time.map(t=>fmt(t,1)),height:310})}
    <div class="kv">
      <span class="k">Model</span><span class="v">degree ${d.model.polynomial_order} + ${d.model.harmonic_count} harmonics</span>
      <span class="k">Growth at start</span><span class="v">${fmt(g.start_per_year,4)} ${esc(d.units)}/yr</span>
      <span class="k">Last five years</span><span class="v">${fmt(g.last_five_years_per_year,4)} ${esc(d.units)}/yr</span>
      <span class="k">Carbon rate</span><span class="v">${fmt(g.end_pgc_per_year,3)} PgC/yr</span>
      <span class="k">Rejected points</span><span class="v">${d.model.rejected_points}</span>
    </div>
    <div class="tagrow">${d.harmonics.map(h=>`<span class="tag">n=${h.number} · ${fmt(h.amplitude,3)} ${esc(d.units)} · ${fmt(h.phase_degrees,1)}°</span>`).join('')}</div>
    ${d.warnings.length?`<ul class="warnlist">${d.warnings.map(w=>`<li>${esc(w)}</li>`).join('')}</ul>`:''}
    <details><summary>Raw JSON</summary><pre>${esc(JSON.stringify(d,null,2))}</pre></details>
  </div></div>`;
}
$('#co2-run').addEventListener('click',async e=>{
  const b=e.target;b.disabled=true;b.innerHTML='<span class="spinner"></span>Fitting…';
  try{
    const d=await api('/api/co2',{csv:$('#co2-csv').value,time:$('#co2-time').value,value:$('#co2-value').value,
      poly:+$('#co2-poly').value,harmonics:+$('#co2-harmonics').value,units:$('#co2-units').value});
    renderCo2(d);
  }catch(err){toast(err.message);$('#co2-output').innerHTML=`<div class="panel"><div class="empty">${esc(err.message)}</div></div>`;}
  b.disabled=false;b.textContent='Fit curve →';
});

$('#sst-run').addEventListener('click',async()=>{
  try{
    const d=await api('/api/retrieve/sst',{bt11:+$('#sst-bt11').value,bt12:+$('#sst-bt12').value,
      first_guess:+$('#sst-guess').value,zenith:+$('#sst-zenith').value,mirror:+$('#sst-mirror').value,
      theta_star:+$('#sst-theta').value});
    const names=['a0','a1','a2','a3','a4','a5','a6'];
    $('#sst-output').innerHTML=`<div class="panel"><div class="panelhead"><b>Retrieved SST</b><span class="hint">NLSST</span></div><div class="panelbody">
      <div class="bigvalue">${fmt(d.result_c,3)} °C</div><div class="unit">${fmt(d.result_k,4)} K</div>
      <div class="cards" style="grid-template-columns:1fr 1fr;margin-top:14px">
        ${metric('Calibration RMSE',fmt(d.model.rmse,3),'K')}
        ${metric('Calibration bias',fmt(d.model.bias,3),'K')}
      </div>
      <div class="coeffs">${d.model.coefficients.map((c,i)=>`<div class="coeff"><b>${names[i]}</b><span>${fmt(c,5)}</span></div>`).join('')}</div>
      <div class="note">${d.model.samples} matchups · calibrated ${d.model.calibrated?'yes':'no'}</div>
    </div></div>`;
  }catch(err){toast(err.message);}
});
$('#vi-index').addEventListener('change',()=>{
  const v=$('#vi-index').value;
  $('#vi-formula').textContent={ndvi:'NDVI = (ρNIR − ρRED) / (ρNIR + ρRED)',
    evi:'EVI = 2.5 (ρNIR − ρRED) / (ρNIR + 6ρRED − 7.5ρBLUE + 1)',
    savi:'SAVI = (ρNIR − ρRED)(1 + L) / (ρNIR + ρRED + L)'}[v];
});
$('#vi-run').addEventListener('click',async()=>{
  try{
    const d=await api('/api/retrieve/index',{index:$('#vi-index').value,red:+$('#vi-red').value,
      nir:+$('#vi-nir').value,blue:+$('#vi-blue').value,soil:+$('#vi-soil').value});
    $('#vi-output').innerHTML=`<div class="panel"><div class="panelhead"><b>${esc(d.index.toUpperCase())}</b><span class="hint">unitless</span></div><div class="panelbody">
      <div class="bigvalue">${fmt(d.value,4)}</div><div class="unit">${d.defined?'defined':'undefined (zero denominator)'}</div>
      <div class="kv"><span class="k">Red</span><span class="v">${fmt(d.red,4)}</span><span class="k">NIR</span><span class="v">${fmt(d.nir,4)}</span>
      <span class="k">Blue</span><span class="v">${fmt(d.blue,4)}</span><span class="k">Soil factor</span><span class="v">${fmt(d.soil,3)}</span></div>
    </div></div>`;
  }catch(err){toast(err.message);}
});
$('#ice-run').addEventListener('click',async()=>{
  try{
    const d=await api('/api/retrieve/seaice',{tb19v:+$('#ice-19v').value,tb19h:+$('#ice-19h').value,tb37v:+$('#ice-37v').value});
    const card=(t,v)=>`<div class="metric"><div class="k">${t}</div><div class="v">${(v*100).toFixed(1)}%</div></div>`;
    $('#ice-output').innerHTML=`<div class="panel"><div class="panelhead"><b>NASA Team concentrations</b><span class="hint">mixture inversion</span></div><div class="panelbody">
      <div class="cards">${card('Open water',d.mixture.open_water)}${card('First-year',d.mixture.first_year)}${card('Multiyear',d.mixture.multiyear)}${card('Total ice',d.mixture.total)}</div>
      <div class="kv"><span class="k">PR(19)</span><span class="v">${fmt(d.pr,5)}</span><span class="k">GR(37V/19V)</span><span class="v">${fmt(d.gr,5)}</span>
      <span class="k">Empirical total</span><span class="v">${(d.empirical.total*100).toFixed(1)}%</span>
      <span class="k">Extent indicator</span><span class="v">${d.empirical.extent_indicator?'ice ≥ 15%':'below 15%'}</span>
      <span class="k">Model RMSE</span><span class="v">${fmt(d.model.rmse_cf,4)} / ${fmt(d.model.rmse_cm,4)}</span></div>
    </div></div>`;
  }catch(err){toast(err.message);}
});
$('#pr-run').addEventListener('click',async()=>{
  try{
    const d=await api('/api/retrieve/precip',{cold_fraction:+$('#pr-cold').value,water_vapour:+$('#pr-wv').value,
      a1:+$('#pr-a1').value,a2:+$('#pr-a2').value});
    $('#pr-output').innerHTML=`<div class="panel"><div class="panelhead"><b>GPI estimate</b><span class="hint">mm/day</span></div><div class="panelbody">
      <div class="bigvalue">${fmt(d.adjusted,3)}</div><div class="unit">mm/day after water-vapour adjustment</div>
      <div class="cards" style="grid-template-columns:1fr 1fr;margin-top:14px">
        ${metric('Raw GPI',fmt(d.gpi,3),'mm/day')}${metric('Adjustment',fmt(d.adjusted-d.gpi,3),'mm/day')}
      </div>
      <div class="kv"><span class="k">Cold fraction</span><span class="v">${fmt(d.cold_fraction,3)}</span>
      <span class="k">Water vapour</span><span class="v">${fmt(d.water_vapour,3)} cm</span>
      <span class="k">exp[a₁(w−a₂)]</span><span class="v">${fmt(d.factor,4)}</span></div>
    </div></div>`;
  }catch(err){toast(err.message);}
});

function parseList(s){return String(s).split(/[\s,;]+/).filter(Boolean).map(Number).filter(v=>isFinite(v));}
$('#mi-run').addEventListener('click',async()=>{
  try{
    const d=await api('/api/info',{x:parseList($('#mi-x').value),y:parseList($('#mi-y').value),
      bins:+$('#mi-bins').value,neighbours:+$('#mi-k').value});
    const max=Math.max(d.plugin,d.adjusted,d.ksg,0.001);
    const bar=(n,v,c)=>`<div style="margin:10px 0"><div style="display:flex;justify-content:space-between;font-size:11px;color:#9db1c8"><span>${n}</span><b style="font-family:var(--mono)">${fmt(v,4)} bits</b></div>
      <div style="height:5px;background:#122036;border-radius:3px;margin-top:5px"><div style="height:5px;width:${Math.max(2,v/max*100)}%;background:${c};border-radius:3px"></div></div></div>`;
    $('#mi-output').innerHTML=`<div class="panel"><div class="panelhead"><b>Mutual information</b><span class="hint">${d.samples} paired samples</span></div><div class="panelbody">
      <div class="bigvalue">${fmt(d.ksg,4)} <span style="font-size:13px;color:#9fc4ef">bits · KSG</span></div>
      ${bar('Histogram (plugin)',d.plugin,'#7f96b2')}
      ${bar('Miller–Madow corrected',d.adjusted,'#38bdf8')}
      ${bar('KSG k-NN',d.ksg,'#5eead4')}
      <div class="kv"><span class="k">Samples</span><span class="v">${d.samples}</span>
      <span class="k">Histogram</span><span class="v">${d.bins} × ${d.bins}</span>
      <span class="k">Neighbours</span><span class="v">k = ${d.neighbours}</span>
      <span class="k">MM correction</span><span class="v">${fmt(d.miller_madow,4)} bits</span></div>
      <div class="note">Negative MI estimates are sampling artefacts; the histogram estimator can return them.</div>
    </div></div>`;
  }catch(err){toast(err.message);}
});
$('#un-run').addEventListener('click',async()=>{
  try{
    const d=await api('/api/sstunc',{weights:parseList($('#un-w').value),sigmas:parseList($('#un-s').value),
      correlation:+$('#un-rho').value,ensemble:parseList($('#un-ens').value)});
    $('#un-output').innerHTML=`<div class="panel"><div class="panelhead"><b>Combined uncertainty</b><span class="hint">1σ</span></div><div class="panelbody">
      <div class="cards">
        ${metric('Total',fmt(d.total_sd,4),'σ units','good')}
        ${metric('Uncorrelated',fmt(d.uncorrelated_sd,4),'quadrature')}
        ${metric('Correlated',fmt(d.correlated_sd,4),'constant covariance')}
        ${metric('Effective n',fmt(d.effective_sample_size,2),'independent samples')}
      </div>
      ${chartCard('Variance budget',[{name:'uncorrelated',values:[d.uncorrelated_variance],color:'#38bdf8'},{name:'correlated',values:[d.correlated_variance],color:'#fbbf24'}],{height:180})}
      <div class="kv"><span class="k">Total variance</span><span class="v">${fmt(d.total_variance,6)}</span>
      <span class="k">Correlation ρ</span><span class="v">${fmt(d.correlation,3)}</span>
      <span class="k">Ensemble component</span><span class="v">${fmt(d.ensemble_sd,4)}</span>
      <span class="k">Weights sum</span><span class="v">${fmt(d.weight_sum,4)}</span>
      <span class="k">dσ²/dρ</span><span class="v">${fmt(d.correlation_sensitivity,6)}</span></div>
      <div class="note">The total includes the ensemble component in quadrature when at least two values are supplied.</div>
    </div></div>`;
  }catch(err){toast(err.message);}
});

const physics={
  'toa':{desc:'N = F + λ ΔT',fields:[['F','Effective forcing (W/m²)',3.7],['lambda','Feedback λ (W/m²/K)',-1.2],['delta_t','ΔT (K)',1.2]]},
  'forcing-co2':{desc:'ΔF = 5.35 ln(C / C₀)',fields:[['ppm','CO₂ concentration (ppm)',420],['reference_ppm','Reference (ppm)',280]]},
  'ecs':{desc:'ECS = −F₂ₓ / λ',fields:[['lambda','Feedback λ (W/m²/K)',-1.25],['forcing_2x','F₂ₓ (W/m²)',3.7]]},
  'mie':{desc:'Qext, Qsca and g from the Mie series',fields:[['size','Size parameter x',5],['real','Real refractive index',1.33],['imag','Imaginary index k',0.01]]},
  'wind':{desc:'τ = ρₐ C_d |U| U',fields:[['u','u₁₀ (m/s)',10],['v','v₁₀ (m/s)',0],['rho','Air density (kg/m³)',1.225],['cd','Drag coefficient',0.0013]]},
  'ekman':{desc:'M = (τy, −τx) / (ρ f)',fields:[['tau_x','τx (N/m²)',0.1],['tau_y','τy (N/m²)',0],['latitude','Latitude (°)',30],['rho','Water density (kg/m³)',1025]]},
  'aod':{desc:'τ = N π r² Qext L',fields:[['n','Number concentration (m⁻³)',1e8],['radius','Radius (µm)',0.1],['qext','Qext',2],['path','Path length (m)',1000]]},
  'equilibrium':{desc:'ΔT = −F / λ',fields:[['F','Forcing (W/m²)',3.7],['lambda','Feedback λ (W/m²/K)',-1.25]]},
  'coriolis':{desc:'f = 2 Ω sin φ',fields:[['latitude','Latitude (°)',30]]},
  'ekman-pumping':{desc:'w = (∂τy/∂x − ∂τx/∂y) / (ρ f)',fields:[['dtau_y_dx','∂τy/∂x (N/m²/m)',1e-8],['dtau_x_dy','∂τx/∂y (N/m²/m)',0],['rho','Water density (kg/m³)',1025],['latitude','Latitude (°)',30]]},
  'geostrophic':{desc:'u = −∂p/∂y/(ρf), v = ∂p/∂x/(ρf)',fields:[['rho','Density (kg/m³)',1025],['dp_dx','∂p/∂x (Pa/m)',1e-5],['dp_dy','∂p/∂y (Pa/m)',0],['latitude','Latitude (°)',30]]},
  'thermal-wind':{desc:'∂u/∂z = g ∂ρ/∂y/(fρ)',fields:[['rho','Density (kg/m³)',1025],['drho_dx','∂ρ/∂x (kg/m⁴)',1e-5],['drho_dy','∂ρ/∂y (kg/m⁴)',0],['latitude','Latitude (°)',30]]},
  'savi':{desc:'SAVI = (NIR−RED)/(NIR+RED+L)(1+L)',fields:[['red','RED reflectance',0.08],['nir','NIR reflectance',0.45],['soil','Soil factor L',0.5]]},
  'evi':{desc:'EVI = G(NIR−RED)/(NIR+C1 RED−C2 BLUE+L)',fields:[['red','RED reflectance',0.08],['nir','NIR reflectance',0.45],['blue','BLUE reflectance',0.04],['gain','Gain G',2.5],['c1','C1',6],['c2','C2',7.5],['canopy','Canopy L',1]]}
};
state.physicsSelected=new Set(['toa']);
function renderPhysicsPicks(){
  $('#ph-picks').innerHTML=Object.entries(physics).map(([k,v])=>`<label><input type="checkbox" data-ph-pick="${k}" ${state.physicsSelected.has(k)?'checked':''}><span>${esc(k.replace(/-/g,' '))}</span></label>`).join('');
}
function renderPhysics(){
  renderPhysicsPicks();
  const chosen=[...state.physicsSelected];
  $('#ph-fields').innerHTML=chosen.map(k=>{
    const def=physics[k];
    return `<div class="ph-card" data-ph-card="${k}"><div class="ph-title">${esc(k.replace(/-/g,' '))}</div><div class="formula">${esc(def.desc)}</div><div class="ph-grid">${def.fields.map(f=>`<div class="field"><label>${esc(f[1])}</label><input data-ph="${k}:${f[0]}" type="number" value="${f[2]}" step="any"></div>`).join('')}</div></div>`;
  }).join('')||'<div class="empty">Select at least one formula.</div>';
  $('#ph-formula-desc').textContent=chosen.length?`${chosen.length} formula${chosen.length===1?'':'s'} selected. Each card has independent inputs; all results are returned together.`:'Select one or more formulas.';
}
$('#ph-picks')&&$('#ph-picks').addEventListener('change',e=>{
  const box=e.target.closest('[data-ph-pick]');if(!box)return;
  if(box.checked)state.physicsSelected.add(box.dataset.phPick);else state.physicsSelected.delete(box.dataset.phPick);
  renderPhysics();
});
$('#ph-all')&&$('#ph-all').addEventListener('click',()=>{state.physicsSelected=new Set(Object.keys(physics));renderPhysics();});
$('#ph-clear')&&$('#ph-clear').addEventListener('click',()=>{state.physicsSelected=new Set();renderPhysics();});
$('#ph-run')&&$('#ph-run').addEventListener('click',async()=>{
  const cards=$$('[data-ph-card]');
  if(!cards.length){toast('Select at least one formula');return;}
  try{
    const results=await Promise.all(cards.map(async card=>{
      const body={formula:card.dataset.phCard};
      card.querySelectorAll('[data-ph]').forEach(i=>{const [,name]=i.dataset.ph.split(':');body[name]=+i.value;});
      return api('/api/physics',body);
    }));
    $('#ph-output').innerHTML=`<div class="panel"><div class="panelhead"><b>Batch results</b><span class="hint">${results.length} formulas</span></div><div class="panelbody">
      <table class="etable"><thead><tr><th>Formula</th><th>Result</th><th>Equation</th></tr></thead><tbody>${results.map(d=>`<tr><td>${esc(d.title)}</td><td>${fmt(d.primary.value,5)} ${esc(d.primary.unit)}</td><td>${esc(d.formula)}</td></tr>`).join('')}</tbody></table>
      ${results.map(d=>`<div class="ph-card"><div class="ph-title">${esc(d.title)}</div><div class="kv">${d.details.map(x=>`<span class="k">${esc(x.label)}</span><span class="v">${fmt(x.value,5)} ${esc(x.unit||'')}</span>`).join('')}</div></div>`).join('')}
    </div></div>`;
  }catch(err){toast(err.message);}
});

function renderCatalog(){
  const q=$('#cat-search').value.toLowerCase(),p=$('#cat-process').value,r=$('#cat-region').value;
  const list=state.catalog.filter(e=>(!p||e.processes.includes(p))&&(!r||e.regions.includes(r))&&
    (!q||(e.id+' '+e.name+' '+e.description+' '+e.retrieval).toLowerCase().includes(q)));
  if(!list.length){$('#cat-list').innerHTML='<div class="empty">No catalog variables match those filters.</div>';return;}
  $('#cat-list').innerHTML=`<div class="varcards">${list.map(e=>`
    <div class="varcard" data-pick="${esc(e.id)}">
      <div class="vtop"><span class="vid">${esc(e.id)}</span><span class="vname">${esc(e.name)} · ${esc(e.units||'—')}</span></div>
      <div class="vchange">${esc(e.change||'Trend in the measured variable')}</div>
      <div class="vmeta">${e.processes.map(x=>esc(x)).join(' · ')}</div>
      <div class="vregions">regions · ${e.regions.map(esc).join(', ')}</div>
      <div class="acts"><button class="btn primary" data-analyze="${esc(e.id)}">Analyse this variable →</button><button class="btn" data-plot="${esc(e.id)}">Plot over time</button></div>
    </div>`).join('')}</div>`;
}
function loadCatalog(){
  api('/api/catalog').then(d=>{state.catalog=d;renderCatalog();renderWfVars();sigOptions();})
    .catch(e=>{$('#cat-list').innerHTML=`<div class="empty">${esc(e.message)}</div>`;});
}
$('#cat-list').addEventListener('click',e=>{
  const a=e.target.closest('[data-analyze]');
  if(a){analyseVariable(a.dataset.analyze);return;}
  const pl=e.target.closest('[data-plot]');
  if(pl){nav('visual');toast('Series for '+pl.dataset.plot+' are in the grid below.');return;}
  const card=e.target.closest('[data-pick]');
  if(card){analyseVariable(card.dataset.pick);}
});
$('#cat-search').addEventListener('input',renderCatalog);
$('#cat-process').addEventListener('change',renderCatalog);
$('#cat-region').addEventListener('change',renderCatalog);

$('#query-run').addEventListener('click',async()=>{
  const text=$('#query-text').value.trim();if(!text)return;
  try{
    const d=await api('/api/query',{text});
    if(!d.valid){$('#query-result').innerHTML=`<div class="card"><div class="empty">Could not resolve: ${esc(d.error)}</div></div>`;return;}
    const it=d.intent||{};
    $('#query-result').innerHTML=`<div class="card"><h3>Compiled selection</h3>
      <div class="tagrow">${it.process?`<span class="tag teal">process: ${esc(it.process)}</span>`:''}
      ${it.region?`<span class="tag">region: ${esc(it.region)}</span>`:''}
      ${it.start_year?`<span class="tag">from ${it.start_year}</span>`:''}
      ${it.end_year?`<span class="tag">to ${it.end_year}</span>`:''}</div>
      <div class="kv" style="margin-top:12px">${d.variables.map(v=>`<span class="k">${esc(v.id)}</span><span class="v">${esc(v.name)}</span>`).join('')}</div>
      <div class="note">${d.reasons.map(esc).join(' ')}</div></div>`;
  }catch(err){toast(err.message);}
});
$('#query-text').addEventListener('keydown',e=>{if(e.key==='Enter')$('#query-run').click();});

/* ---------------- TREND DETECTIVE ---------------- */
const WFMETA={
  sst:{file:'sst_monthly.csv',value:'sst',sphere:'OCEAN',product:'MODIS Aqua L3 SST',pid:'MODISA_L3m_SST',provider:'PO.DAAC',res:'4 km / 8-day composite',grid:false,
    artifacts:['Skin vs. bulk SST: satellite radiometer measures the skin; not corrected in this demo.','Diurnal thermocline: nighttime product reduces, but does not remove, diurnal heating.','Sensor drift: MODIS Aqua is intercalibrated in the CDR, not in the synthetic demo.'],
    retrieval:'NLSST: SST = a0 + a1 BT11 + a2 (BT11-BT12) Ts + a3 (sec z-1)(BT11-BT12) + ...',analysis:'STL decomposition + Mann-Kendall + Sen slope'},
  air_temperature:{file:'live_gistemp.csv',value:'value',sphere:'ATMOSPHERE',product:'AIRS / MERRA-2 surface air temperature',pid:'AIRS L3 / M2T1NXSLV',provider:'GES DISC',res:'1 degree monthly',grid:false,
    artifacts:['This file is the real NASA GISTEMP v4 global land-ocean anomaly, not an AIRS or MERRA-2 download.','It is an anomaly against 1951-1980, not an absolute temperature.','Station coverage is spatially uneven, especially in the early decades.'],
    retrieval:'Gridded temperature retrieval / reanalysis variable',analysis:'STL decomposition + Mann-Kendall + Sen slope'},
  co2:{file:'live_co2_mlo.csv',value:'value',sphere:'ATMOSPHERE',product:'OCO-2 Lite XCO2',pid:'OCO2_L2_Lite_FP',provider:'GES DISC',res:'point retrievals, monthly composite',grid:false,
    artifacts:['This file is the real NOAA GML Mauna Loa in-situ monthly mean, not an XCO2 column product.','Mauna Loa is one high-altitude station; it is not a global mean.','Missing months are dropped, so the time axis is uneven.'],
    retrieval:'OCO-2 XCO2 retrieval',analysis:'Polynomial + harmonics (CO2 lab); STL + MK in this workflow'},
  ndvi:{file:'ndvi_monthly.csv',value:'value',sphere:'LAND / BIOSPHERE',product:'MODIS/Terra vegetation indices',pid:'MOD13A3',provider:'LP DAAC',res:'1 km monthly',grid:false,
    artifacts:['NDVI saturates in dense canopy.','Cloud and aerosol contamination persist in composites.','Land-cover change can mimic a trend.'],
    retrieval:'NDVI = (NIR - RED)/(NIR + RED)',analysis:'STL decomposition + Mann-Kendall + Sen slope'},
  precipitation:{file:'live_power_prectotcorr_65_25.csv',value:'value',sphere:'ATMOSPHERE',product:'GPCP monthly precipitation',pid:'GPCPMON',provider:'GES DISC',res:'0.5 degree monthly',grid:false,
    artifacts:['This file is the real NASA POWER monthly precipitation at 65N 25E, not a demo series.','It is one model grid cell, not a regional or global mean.','POWER reanalysis precipitation carries model and gauge bias.'],
    retrieval:'GPI cold-cloud index + GPCP inverse-variance merge',analysis:'STL decomposition + Mann-Kendall + Sen slope'},
  sea_ice_extent:{file:'live_seaice_n.csv',value:'value',sphere:'CRYOSPHERE',product:'NSIDC sea ice concentration CDR',pid:'NSIDC-0051',provider:'NSIDC DAAC',res:'25 km / monthly extent',grid:false,
    artifacts:['This file is the real NSIDC Sea Ice Index northern monthly extent, not a demo series.','Extent uses a 15% concentration threshold; it is not area.','Passive-microwave retrievals are weather-dependent.'],
    retrieval:'NASA Team concentration + 15% extent rule',analysis:'STL decomposition + Mann-Kendall + Sen slope'},
  sea_ice_concentration:{file:null,value:'value',sphere:'CRYOSPHERE',product:'NSIDC sea ice concentration CDR',pid:'NSIDC-0051',provider:'NSIDC DAAC',res:'25 km monthly',grid:true,
    artifacts:['Concentration is a fractional coverage, not thickness.','Melt ponds and snow cover affect the brightness temperatures.','The 15% extent record is a different variable (km2, thresholded) and is never substituted for concentration; this variable needs its own series.'],
    retrieval:'NASA Team / Bootstrap concentration',analysis:'Area-weighted anomalies + Mann-Kendall + Sen slope'},
  ch4:{file:'live_ch4_gl.csv',value:'value',sphere:'ATMOSPHERE',product:'OCO-2 / TROPOMI column methane',pid:'OCO2_L2_Lite_FP',provider:'GES DISC',res:'point retrievals, monthly composite',grid:false,
    artifacts:['This file is the real NOAA GML globally averaged marine-surface monthly mean, not a column product.','The seasonal cycle is large, so trend work needs deseasonalised data.','A few months are missing where the network was too sparse.'],
    retrieval:'Column-averaged CH4 retrieval',analysis:'Polynomial + harmonics (Thoning et al.); growth rate from the fitted trend'},
  ozone:{file:'ozone_monthly.csv',value:'value',sphere:'ATMOSPHERE',product:'OMI total column ozone',pid:'OMDOAO3e',provider:'GES DISC',res:'1 degree monthly',grid:false,
    artifacts:['Total column ozone is not surface ozone.','Cloud and snow/ice surfaces bias the retrieval.','Demo series is synthetic; it is not an OMI download.'],
    retrieval:'OMI DOAS total column ozone retrieval',analysis:'STL decomposition + Mann-Kendall + Sen slope'},
  aerosol_optical_depth:{file:'aerosol_monthly.csv',value:'value',sphere:'ATMOSPHERE',product:'OMI aerosol optical depth',pid:'OMAERO',provider:'GES DISC',res:'1 degree monthly',grid:false,
    artifacts:['Aerosol index and optical depth are not interchangeable.','Cloud contamination and surface albedo remain.','Demo series is synthetic; it is not an OMI download.'],
    retrieval:'OMI near-UV aerosol retrieval',analysis:'STL decomposition + Mann-Kendall + Sen slope'},
  soil_moisture:{file:'soil_moisture_monthly.csv',value:'value',sphere:'LAND / BIOSPHERE',product:'SMAP surface soil moisture',pid:'SPL3SMP_E',provider:'NSIDC DAAC',res:'9 km daily composite',grid:false,
    artifacts:['Surface soil moisture is not root-zone water.','Freeze/thaw and dense vegetation degrade the retrieval.','Demo series is synthetic; it is not an SMAP download.'],
    retrieval:'SMAP L-band brightness-temperature inversion',analysis:'STL decomposition + Mann-Kendall + Sen slope'},
  ocean_color:{file:'ocean_color_monthly.csv',value:'value',sphere:'OCEAN',product:'MODIS Aqua chlorophyll-a',pid:'MODISA_L3m_CHL',provider:'OB.DAAC',res:'4 km monthly',grid:false,
    artifacts:['Chlorophyll-a is a proxy for biomass, not carbon.','Clouds and coloured dissolved matter bias the blue/green ratio.','Demo series is synthetic; it is not a MODIS download.'],
    retrieval:'MODIS ocean-colour chlorophyll-a algorithm',analysis:'STL decomposition + Mann-Kendall + Sen slope'},
  ocean_heat_content:{sphere:'OCEAN',product:'Ocean heat content (Armor/EN4)',pid:'OHC',provider:'PO.DAAC',res:'gridded monthly',grid:true,
    artifacts:['OHC depends on the mapping method and climatology.','Deep-ocean coverage is sparse.','No local demo series is bundled.'],
    retrieval:'depth-integrated temperature anomaly',analysis:'Trend and acceleration diagnostics'},
  toa_energy_imbalance:{sphere:'ATMOSPHERE',product:'CERES TOA energy imbalance',pid:'CERES_EBAF',provider:'ASDC',res:'1 degree monthly',grid:true,
    artifacts:['TOA imbalance is a small residual of large upward and downward fluxes.','Instrument calibration and diurnal sampling matter.','No local demo series is bundled.'],
    retrieval:'net shortwave + net longwave at TOA',analysis:'Trend and acceleration diagnostics'},
  snow_cover:{sphere:'CRYOSPHERE',product:'MODIS snow cover',pid:'MOD10CM',provider:'NSIDC DAAC',res:'0.05 degree monthly',grid:true,
    artifacts:['Snow under cloud is not observed.','Extent is sensitive to the snow-cover threshold.','No local demo series is bundled.'],
    retrieval:'MOD10CM snow-cover fraction',analysis:'STL decomposition + Mann-Kendall + Sen slope'},
  sea_level:{sphere:'OCEAN',product:'Satellite altimetry sea level',pid:'SEA_LEVEL_CDR',provider:'PO.DAAC',res:'global mean monthly',grid:false,
    artifacts:['Global mean sea level includes steric and mass components.','Altimetry needs glacial isostatic adjustment correction.','No local demo series is bundled.'],
    retrieval:'altimetry range + geophysical corrections',analysis:'Trend and acceleration diagnostics'}
};
function sphereOf(e){
  const p=((e&&e.processes)||[]).join(' ').toLowerCase();
  if(/sea_ice|cryo|snow/.test(p))return 'CRYOSPHERE';
  if(/ocean|heat_content|sea_level|energy_balance/.test(p))return 'OCEAN';
  if(/vegetation|drought|land|biosphere/.test(p))return 'LAND / BIOSPHERE';
  return 'ATMOSPHERE';
}
function wfMeta(id){
  return WFMETA[id]||{sphere:'ATMOSPHERE',product:'Catalog variable',pid:id,provider:'—',res:'CSV series',grid:false,
    artifacts:['No retrieval-specific artifact note is registered for this variable.'],retrieval:'See catalog entry',analysis:'STL decomposition + Mann-Kendall + Sen slope'};
}
function wfLog(msg,kind=''){
  const el=$('#wf-log'); if(!el)return;
  const t=new Date().toLocaleTimeString('en-GB',{hour12:false});
  el.insertAdjacentHTML('beforeend',`<div class="logline ${kind}"><b>${t}</b> ${esc(msg)}</div>`);
  el.scrollTop=el.scrollHeight;
}
function wfStep(n,stateName){
  const el=$('#wf-step'+n); if(!el)return;
  el.classList.toggle('open',stateName!=='locked');
  el.classList.toggle('locked',stateName==='locked');
  el.classList.toggle('done',stateName==='done');
  const status=$('#wf-step'+n+'-status');
  if(status)status.textContent=stateName==='done'?'COMPLETE':stateName==='locked'?'LOCKED':'ACTIVE';
}
function selectWfVar(id){
  if(!state.catalog.length)return;
  const entry=state.catalog.find(e=>e.id===id)||state.catalog[0];
  const meta=wfMeta(entry.id);
  const changed=!state.workflow||!state.workflow.csv||state.workflow.entry.id!==entry.id;
  state.workflow={entry,meta};
  if(changed){
    $('#wf-confirm').disabled=true;
    $('#wf-data').innerHTML='<div class="empty">Data gate reset for the new variable. Press LOAD DATA.</div>';
    $('#wf-quality').innerHTML='';
    ['2','3','4','5'].forEach(n=>wfStep(n,'locked'));
    $('#wf-step1-status').textContent='AWAITING LOAD';
  }
  $$('#wf-vars .varitem').forEach(x=>x.classList.toggle('active',x.dataset.id===entry.id));
  const region=$('#wf-region');
  const regions=(entry.regions&&entry.regions.length)?entry.regions:['global'];
  region.innerHTML=regions.map(r=>`<option>${esc(r)}</option>`).join('');
  $('#wf-artifacts').innerHTML=meta.artifacts.map(a=>`<li>${esc(a)}</li>`).join('');
  $('#wf-active').innerHTML=`<div style="font-size:10px;color:#f0dfc2;text-transform:uppercase;letter-spacing:.06em">${esc(entry.id)} / ${esc(meta.product)}</div>
    <div class="formula">${esc(meta.retrieval)}<br><br>${esc(meta.analysis)}</div>`;
  $('#wf-io').innerHTML=`<span class="k">VARIABLE</span><span class="v">${esc(entry.id)}</span><span class="k">UNITS</span><span class="v">${esc(entry.units||'—')}</span><span class="k">REGION</span><span class="v">${esc(regions[0])}</span><span class="k">STATUS</span><span class="v">SELECTED</span>`;
  $('#wf-load').disabled=false;
}
function renderWfVars(){
  const box=$('#wf-vars'); if(!box)return;
  if(!state.catalog.length){box.innerHTML='<div class="empty">Catalog unavailable.</div>';return;}
  const q=($('#wf-search').value||'').toLowerCase();
  const sphere=$('#wf-sphere').value;
  const matches=state.catalog.filter(e=>{
    const hay=(e.id+' '+e.name+' '+e.description+' '+e.processes.join(' ')).toLowerCase();
    return (!q||hay.includes(q))&&(!sphere||sphereOf(e)===sphere);
  });
  const groups={};
  matches.forEach(e=>{const s=sphereOf(e);(groups[s]=groups[s]||[]).push(e);});
  let html='';
  ['ATMOSPHERE','OCEAN','CRYOSPHERE','LAND / BIOSPHERE'].forEach(s=>{
    if(!groups[s])return;
    html+=`<div class="sphere">${s}</div>`;
    html+=groups[s].map(e=>`<div class="varitem ${state.workflow&&state.workflow.entry.id===e.id?'active':''}" data-id="${esc(e.id)}">
      <span class="code">${esc(e.id.slice(0,3).toUpperCase())}</span><span><span class="nm">${esc(e.name)}</span><span class="meta">${esc(e.units||'—')} / ${esc(e.processes.slice(0,2).join(' · '))}</span></span></div>`).join('');
  });
  box.innerHTML=html||'<div class="empty">No variables match this filter.</div>';
  if(!state.workflow&&matches.length)selectWfVar(matches[0].id);
}
$('#wf-search')&&$('#wf-search').addEventListener('input',renderWfVars);
$('#wf-sphere')&&$('#wf-sphere').addEventListener('change',renderWfVars);
$('#wf-vars')&&$('#wf-vars').addEventListener('click',e=>{const item=e.target.closest('.varitem');if(item)analyseVariable(item.dataset.id);});
async function analyseVariable(id){
  if(!id)return;
  if(!state.catalog.length){try{state.catalog=await api('/api/catalog');}catch(e){toast(e.message);return;}}
  nav('workflow');
  selectWfVar(id);
  try{
    await loadWorkflowData();
    const confirm=$('#wf-confirm');
    if(confirm&&!confirm.disabled)confirm.click();
  }catch(e){toast(e.message);}
}
function parseCsv(text){
  const lines=text.split(/\r?\n/).filter(l=>l.trim());
  if(lines.length<2)throw new Error('CSV needs a header and at least one row');
  const header=lines[0].split(',').map(s=>s.trim().replace(/^"|"$/g,''));
  const rows=lines.slice(1).map(l=>l.split(',').map(s=>s.trim().replace(/^"|"$/g,'')));
  return {header,rows};
}
function isNum(s){return s!==''&&isFinite(Number(s));}
function gapRuns(values){
  const runs=[];let start=-1;
  values.forEach((v,i)=>{if(!isFinite(v)){if(start<0)start=i;}else if(start>=0){runs.push([start,i-1]);start=-1;}});
  if(start>=0)runs.push([start,values.length-1]);
  return runs;
}
function applyWorkflowCsv(text,source){
  const table=parseCsv(text);
  const entry=state.workflow.entry,meta=state.workflow.meta;
  const ti=Math.max(0,table.header.findIndex(h=>h.toLowerCase()==='time'));
  const lower=table.header.map(h=>h.toLowerCase());
  const declared=meta.value?String(meta.value):'';
  let vi=declared?lower.indexOf(declared.toLowerCase()):-1;
  if(vi<0)vi=lower.indexOf('value');
  if(declared&&vi<0)
    throw new Error('this file has no '+declared+" or 'value' column, so it is not a "+entry.id+' record');
  if(vi<0)vi=table.header.length>1?1:0;
  const start=Number($('#wf-start').value)||0,end=Number($('#wf-end').value)||9999;
  const rows=table.rows.filter(r=>{const t=Number(r[ti]);return isFinite(t)&&t>=start-0.0001&&t<=end+1;});
  if(rows.length<3)throw new Error('fewer than three observations in the selected period');
  const times=rows.map(r=>Number(r[ti])),values=rows.map(r=>Number(r[vi]));
  const filtered=[table.header.join(',')].concat(rows.map(r=>r.join(','))).join('\n');
  const missing=values.filter(v=>!isFinite(v)).length;
  const years=(end-start+1);
  const coverage=Math.min(100,rows.length/Math.max(1,years*12)*100);
  const runs=gapRuns(values);
  state.workflow={entry,meta,source,csv:filtered,times,values,missing,coverage,runs,ti,vi,start,end,
    timeName:table.header[ti]||'time',valueName:table.header[vi]||'value',
    min:Math.min(...values.filter(v=>isFinite(v))),max:Math.max(...values.filter(v=>isFinite(v)))};
  const finite=rows.length-missing;
  const gapText=runs.length?runs.slice(0,4).map(([a,b])=>`${fmt(times[a],2)} → ${fmt(times[b],2)} (${b-a+1} steps)`).join(' · '):'none detected';
  $('#wf-data').innerHTML=`<div class="statusline"><span>ACQUISITION <b>SUCCESS</b></span><span>SOURCE <b>${esc(source)}</b></span><span>RECORDS <b>${rows.length}</b></span></div>
    <div class="facts">
      <div class="fact"><div class="k">Source</div><div class="v">${esc(source)}</div><div class="s">${esc(meta.provider||'local')}</div></div>
      <div class="fact"><div class="k">NASA product route</div><div class="v">${esc(meta.product)}</div><div class="s">${esc(meta.pid)}</div></div>
      <div class="fact"><div class="k">Temporal range</div><div class="v">${fmt(times[0],2)} → ${fmt(times[times.length-1],2)}</div><div class="s">${fmt(rows.length,0)} observations</div></div>
      <div class="fact"><div class="k">Resolution</div><div class="v">${esc(meta.res)}</div><div class="s">native product spec</div></div>
      <div class="fact"><div class="k">Missing values</div><div class="v">${missing}</div><div class="s">${finite} finite samples</div></div>
      <div class="fact"><div class="k">Spatial support</div><div class="v">${meta.grid?'GRID FIELD':'SCALAR SERIES'}</div><div class="s">${meta.grid?'map inference available':'no map inference'}</div></div>
    </div>
    <div class="note">Temporal coverage</div><div class="bar ${coverage<90?'warn':''}"><i style="width:${coverage.toFixed(1)}%"></i></div>
    <div class="kv"><span class="k">COVERAGE</span><span class="v">${coverage.toFixed(1)}% of requested months</span><span class="k">GAPS</span><span class="v">${esc(gapText)}</span><span class="k">RANGE</span><span class="v">${fmt(state.workflow.min,3)} … ${fmt(state.workflow.max,3)} ${esc(entry.units||'')}</span></div>`;
  $('#wf-quality').innerHTML=`<div class="note" style="margin-top:8px">Quality flags encoded in the local CSV: <span class="telemetry">NONE</span>. Source-product flags are not present; do not infer cloud-mask or retrieval quality from this fixture.</div>`;
  $('#wf-confirm').disabled=false;
  $('#wf-step1-status').textContent='LOADED / CONFIRM';
  $('#wf-io').innerHTML=`<span class="k">VARIABLE</span><span class="v">${esc(entry.id)}</span><span class="k">SAMPLES</span><span class="v">${rows.length}</span><span class="k">MISSING</span><span class="v">${missing}</span><span class="k">STATUS</span><span class="v">AWAITING CONFIRMATION</span>`;
  wfLog(`Data loaded: ${entry.id}, ${rows.length} records, ${coverage.toFixed(1)}% coverage.`);
  wfLog('No statistic has been computed yet. Confirm the data gate to proceed.');
  wfStep(1,'active');
  ['2','3','4','5'].forEach(n=>wfStep(n,'locked'));
}
async function loadWorkflowData(){
  const entry=state.workflow&&state.workflow.entry; if(!entry)return toast('Select a variable first.');
  const meta=wfMeta(entry.id);
  const file=meta.file;
  if(!file){
    wfLog(`${entry.id} has no bundled series; no other variable's record is used in its place.`,'warn');
    toast(`${entry.id} needs its own CSV (columns: time, ${meta.value||'value'}).`);
    $('#wf-data').innerHTML='<div class="empty">This variable has no bundled series, and no other variable\'s record is used in its place. Load its own CSV.</div>';
    $('#wf-upload').click(); return;
  }
  if(!state.overview||!state.overview.files.includes(file)){
    wfLog(`The bundled file ${file} is missing; requesting upload.`,'warn');
    toast('No bundled demo for this variable. Choose a CSV file.');
    $('#wf-upload').click(); return;
  }
  try{
    const text=await (await fetch('/api/data/file?name='+encodeURIComponent(file))).text();
    applyWorkflowCsv(text,'LOCAL DEMO / '+file);
  }catch(e){toast(e.message);wfLog(e.message,'err');}
}
$('#wf-load')&&$('#wf-load').addEventListener('click',loadWorkflowData);
$('#wf-upload')&&$('#wf-upload').addEventListener('change',e=>{
  const f=e.target.files[0];if(!f)return;const r=new FileReader();
  r.onload=()=>{try{applyWorkflowCsv(r.result,'UPLOADED CSV / '+f.name);}catch(err){toast(err.message);}};r.readAsText(f);
});
function wfChart(d){
  const labels=d.series.time.map(t=>fmt(t,1));
  const title=d.method==='harmonic'?'Observed and fitted polynomial trend':'Observed and STL trend';
  return chartCard(title,[{name:'observed',values:d.series.observed,color:'#f7e8c8'},{name:'trend',values:d.components.trend,color:'#ffab54'}],{xlabels:labels,height:250})+
    `<div class="minicharts">${chartCard('Seasonal',[{name:'seasonal',values:d.components.seasonal,color:'#45d6c9'}],{xlabels:labels,height:210})}${chartCard('Remainder',[{name:'remainder',values:d.components.remainder,color:'#ff8f86'}],{xlabels:labels,height:210})}</div>`;
}
function wfTrend(d){
  const mk=d.trend.mann_kendall,sen=d.trend.sen;
  const verdict=mk.significant?`<b>SIGNIFICANT ${esc(mk.direction).toUpperCase()}</b> at alpha = ${d.alpha}.`:`<b>NO SIGNIFICANT MONOTONIC TREND</b> at alpha = ${d.alpha}.`;
  const ci=sen.computed?`${fmt(sen.lower,4)} … ${fmt(sen.upper,4)}`:'not available';
  const ens=d.trend.ensemble||{members:[],agreement:0,consensus:''};
  const rows=(ens.members||[]).filter(m=>m.computed).map(m=>`<tr class="${m.selected?'primary':''}"><td>${esc(m.name)}${m.selected?' <span class="pick">PRIMARY</span>':' <span class="fallback">FALLBACK</span>'}</td><td>${fmt(m.slope,4)}</td><td>${pval(m.p_value)}</td><td>${fmt(m.rmse,3)}</td><td>${fmt(m.score,3)}</td></tr>`).join('');
  return `<div class="verdict">${verdict}</div>
    <div class="cards">${metric('Sen slope',fmt(sen.slope,4),esc(d.units)+' / step')}${metric('Per decade',fmt(d.trend.sen_per_decade,3),esc(d.units)+' / decade')}${metric('MK p-value',pval(mk.p_value),mk.direction,mk.significant?'good':'warn')}${metric('95% CI',ci,'Sen slope','good')}</div>
    <div class="panel" style="margin-top:8px"><div class="panelhead"><b>Estimator ensemble</b><span class="hint">accuracy-weighted</span></div><div class="panelbody">
      <div class="statusline"><span>PRIMARY <b>${esc((ens.members||[]).find(m=>m.selected)?.name||'—')}</b></span><span>AGREEMENT <b>${Math.round(ens.agreement*100)}%</b></span><span>CONSENSUS <b>${esc(ens.consensus)}</b></span></div>
      <table class="etable"><thead><tr><th>Estimator</th><th>Slope</th><th>p</th><th>RMSE</th><th>Score</th></tr></thead><tbody>${rows}</tbody></table>
      <div class="note">Every estimator sees the same data. The primary is the highest accuracy-weighted score (residual RMSE, significance, sign agreement); the rest remain as auditable fallbacks. This is why a single method is not allowed to hide disagreement.</div>
    </div></div>
    <div class="cards">${metric('MK S',fmt(mk.s,0),'statistic')}${metric('Variance',fmt(mk.variance,0),'tie-corrected')}${metric('Z',fmt(mk.z,3),'normal score')}${metric('Correction',fmt(mk.correction_factor,3)+'x','Hamed-Rao factor')}</div>
    <div class="kv"><span class="k">SERIES TESTED</span><span class="v">deseasonalised</span><span class="k">SAMPLE / CORR</span><span class="v">${d.samples} / ${fmt(mk.correction_factor,3)}x</span><span class="k">RAW MK p</span><span class="v">${pval(d.raw_trend.mann_kendall.p_value)}</span><span class="k">SEASONAL MK p</span><span class="v">${d.seasonal_mann_kendall?pval(d.seasonal_mann_kendall.p_value):'—'}</span></div>`;
}
function wfSpatial(d){
  const meta=state.workflow.meta;
  return `${meta.grid?'':'<div class="verdict">The selected acquisition is a scalar series, so its own map cannot be drawn. The grid below is a separate synthetic demonstration of the spatial path — it is not the selected data.</div>'}
    <div class="statusline"><span>GRID <b>24 × 12 DEMO</b></span><span>TREND <b>SEN / CELL</b></span><span>HOTSPOT <b>GETIS-ORD Gi*</b></span></div>
    <div id="wf-grid"><div class="empty">Computing per-cell slopes…</div></div>
    <div class="note">Synthetic demonstration field, not a NASA download and not the selected series. It exercises the exact per-cell Sen + Getis-Ord code that a gridded acquisition would use. The Walker threshold is implemented in the spatial module; field-significance is optional, not part of the core verdict.</div>`;
}
async function loadGridDemo(){
  const box=$('#wf-grid');if(!box)return;
  try{
    const d=await api('/api/grid-trend');
    const maxAbs=Math.max(...d.slopes.map(v=>Math.abs(v)),1e-9);
    const cells=d.slopes.map((s,i)=>{
      const lat=d.latitudes[Math.floor(i/d.columns)],lon=d.longitudes[i%d.columns];
      const t=Math.min(1,Math.abs(s)/maxAbs);
      const hue=s>=0?30:175;
      const light=18+38*t;
      const sat=s>=0?78:55;
      const hot=d.hotspots[i];
      return `<div class="gridcell" title="lat ${fmt(lat,1)} · lon ${fmt(lon,1)} · Sen ${fmt(s,5)} / month${hot?' · Gi* hotspot':''}" style="background:hsl(${hue} ${sat}% ${light}%);${hot?'outline:1px solid #f3dfb5;outline-offset:-1px':''}"></div>`;
    }).join('');
    box.innerHTML=`<div class="gridmap" style="grid-template-columns:repeat(${d.columns},1fr)">${cells}</div>
      <div class="tagrow" style="margin-top:6px"><span class="tag">COOLING / TEAL</span><span class="tag amber">WARMING / ORANGE</span><span class="tag">${d.hotspot_count} Gi* HOTSPOTS</span></div>`;
  }catch(e){box.innerHTML=`<div class="empty">${esc(e.message)}</div>`;}
}
function wfQuestionBar(d){
  const e=state.workflow.entry,meta=state.workflow.meta,mk=d.trend.mann_kendall,sen=d.trend.sen;
  const region=$('#wf-region')?$('#wf-region').value:'global';
  const unit=esc(e.units||'');
  const set=(id,text)=>{const el=$(id);if(!el)return;el.textContent=text;el.parentElement.classList.add('set');};
  set('#wf-q-what',`${e.id} · ${mk.direction}`);
  set('#wf-q-where',`${meta.grid?'grid field · '+region:'series region tag: '+region}`);
  set('#wf-q-howmuch',`${fmt(sen.slope,4)} ${unit}/step · ${fmt(d.trend.sen_per_decade,2)} ${unit}/decade`);
  set('#wf-q-sig',mk.significant?`yes · p = ${pval(mk.p_value)}`:`no · p = ${pval(mk.p_value)}`);
}
function wfFormulaBridge(d){
  const box=$('#wf-formula-bridge');if(!box)return;
  const mk=d.trend.mann_kendall,sen=d.trend.sen,smk=d.seasonal_mann_kendall;
  const ens=d.trend.ensemble||{members:[]};
  const unit=esc(d.units||'');
  const rows=[];
  const row=(q,name,how,value)=>rows.push({q,name,how,value});
  row('WHAT','Mann–Kendall S','sign of the rank-sum S over all pairs',fmt(mk.s,0)+' ('+esc(mk.direction)+')');
  row('WHAT','Kendall τ','S normalised by n(n−1)/2',fmt(mk.tau,3));
  (ens.members||[]).filter(m=>/tau|kendall/i.test(m.name)).forEach(m=>row('WHAT',m.name,m.formula||'rank concordance',fmt(m.slope,4)+' '+unit));
  row('WHERE','Sign partition','two regional Sen slopes from /api/contrast','open Regional contrast');
  if(d.grid)row('WHERE','Getis-Ord Gi*','local hotspot statistic per cell','grid field');
  row('HOW MUCH','Sen slope','median of all pairwise slopes','0.0137'.replace('0.0137',fmt(sen.slope,4))+' '+unit+'/step');
  (ens.members||[]).filter(m=>/ols|least/i.test(m.name)).forEach(m=>row('HOW MUCH',m.name,m.formula||'least squares',fmt(m.slope,4)+' '+unit));
  row('HOW MUCH','Per decade','Sen slope × steps per decade',fmt(d.trend.sen_per_decade,3)+' '+unit);
  row('HOW MUCH','Relative change','total change ÷ |mean| × 100','see report');
  row('SIGNIFICANT','Mann–Kendall p','two-sided normal approximation of S',pval(mk.p_value));
  row('SIGNIFICANT','Hamed–Rao factor','variance inflated for serial correlation',fmt(mk.correction_factor,3)+'x');
  if(smk)row('SIGNIFICANT','Seasonal Mann–Kendall p','months compared across years',pval(smk.p_value));
  row('SIGNIFICANT','Ensemble agreement','fraction of estimators sharing the sign',Math.round((ens.agreement||0)*100)+'%');
  row('SIGNIFICANT','Minimum detectable trend','2.802 σ √(12/n) / span',(mdtPerYear(d)===null?'—':fmt(mdtPerYear(d),4)+' '+unit+'/yr'));
  const order=['WHAT','WHERE','HOW MUCH','SIGNIFICANT'];
  box.innerHTML=order.map(q=>`<div class="bridge-group"><div class="bridge-q">${q}</div><table class="etable"><thead><tr><th>Formula</th><th>How it calculates the variable</th><th>Value on this record</th></tr></thead><tbody>
    ${rows.filter(r=>r.q===q).map(r=>`<tr><td>${esc(r.name)}</td><td class="mut">${esc(r.how)}</td><td>${r.value}</td></tr>`).join('')}
  </tbody></table></div>`).join('')+`<div class="note">Each question is answered by several formulas, not one: the estimators agree or they do not, and the disagreement is printed rather than averaged. The same set is reachable as individual calculations in Calculation modes and runnable with your own inputs in All calculations.</div>`;
}
async function loadFormulaPlan(entryId){
  const box=$('#wf-plan');if(!box)return;
  box.innerHTML='<div class="empty">Selecting every formula that applies to this variable…</div>';
  try{
    const d=await api('/api/plan/run',{variable:entryId,csv:state.workflow.csv,time:state.workflow.timeName,value:state.workflow.valueName,period:12});
    const rows=(d.results||[]).map(r=>`<tr><td class="mut">${esc(r.stage)}</td><td>${esc(r.id)}</td><td class="mut">${esc(r.role)}</td><td>${r.skipped?`<span class="mut">method step · needs sensor inputs</span>`:(r.ok?`<b>${esc(String(r.result).slice(0,140))}</b>`:`<span class="warn">${esc(r.error||'rejected')}</span>`)}</td></tr>`).join('');
    const hint=$('#wf-plan-hint');
    if(hint)hint.textContent=`${d.selected} of ${d.registry_size} registry formulas · ${d.executed} computed · ${d.skipped} method steps need sensor inputs · ${d.rejected} rejected`;
    box.innerHTML=`<div class="statusline"><span>METHOD <b>${esc(String(d.method).toUpperCase())}</b></span><span>SELECTED <b>${d.selected}</b></span><span>EXECUTED <b>${d.executed}</b></span><span>SAMPLES <b>${d.samples}</b></span></div>
      <table class="etable"><thead><tr><th>Stage</th><th>Formula</th><th>Why it belongs to this variable</th><th>Value on this record</th></tr></thead><tbody>${rows}</tbody></table>
      <div class="note">${(d.notes||[]).map(n=>esc(n)).join(' · ')}</div>`;
  }catch(e){box.innerHTML=`<div class="empty">${esc(e.message)}</div>`;}
}
function renderWfReport(d){
  wfQuestionBar(d);
  wfFormulaBridge(d);
  const e=state.workflow.entry,meta=state.workflow.meta,mk=d.trend.mann_kendall,sen=d.trend.sen;
  $('#wf-report').innerHTML=`<div class="statusline"><span>REPORT <b>GENERATED</b></span><span>VARIABLE <b>${esc(e.id)}</b></span><span>REGION <b>${esc($('#wf-region').value)}</b></span></div>
    <div class="verdict">${mk.significant?'The series shows a statistically significant '+esc(mk.direction)+' trend.':'No statistically significant monotonic trend was detected.'}</div>
    <div class="facts"><div class="fact"><div class="k">What</div><div class="v">${esc(e.name)}</div><div class="s">${esc(e.units||'')}</div></div><div class="fact"><div class="k">Where</div><div class="v">${esc($('#wf-region').value)}</div><div class="s">catalog region tag</div></div><div class="fact"><div class="k">How much</div><div class="v">${fmt(d.trend.sen_per_decade,3)} / decade</div><div class="s">Sen slope</div></div><div class="fact"><div class="k">Significance</div><div class="v">p ${pval(mk.p_value)}</div><div class="s">Mann-Kendall</div></div><div class="fact"><div class="k">Coverage</div><div class="v">${state.workflow.coverage.toFixed(1)}%</div><div class="s">requested months</div></div><div class="fact"><div class="k">Gaps</div><div class="v">${state.workflow.runs.length}</div><div class="s">runs of missing data</div></div></div>
    <div class="note">Formulas applied · ${d.method==='harmonic'?'polyharmonic':'STL'} route chosen for ${esc(e.id)}</div><div class="formula">${esc(meta.retrieval)}<br>${d.method==='harmonic'?'Polynomial + harmonics (Thoning et al. 1989) with iterative outlier rejection; growth rate from the fitted trend':'STL: Y = T + S + R (Cleveland et al.)'}<br>Mann-Kendall S with tie correction + Hamed-Rao variance<br>Sen slope + rank-based 95% CI${d.method==='harmonic'?`<br>Growth at end of record: ${fmt(d.growth.end_per_year,3)} ${esc(e.units||'')}/yr`:''}</div>
    <div class="note">Artifacts to disclose</div><ul class="warnlist">${meta.artifacts.map(a=>`<li>${esc(a)}</li>`).join('')}</ul>
    <div class="actions"><button class="btn" data-wfact="copy">Copy JSON</button><button class="btn" data-wfact="csv">Download CSV</button><button class="btn" data-wfact="print">Print report</button></div>`;
}
async function runWorkflow(){
  if(!state.workflow||!state.workflow.csv)return toast('Load and confirm data first.');
  const b=$('#wf-confirm');b.disabled=true;b.innerHTML='<span class="spinner"></span>Computing…';
  try{
    wfStep(1,'done');wfStep(2,'active');
    $('#wf-step1-status').textContent='CONFIRMED';
    const e=state.workflow.entry;
    wfLog('Data gate confirmed. The engine selects the method from the variable.');
    const d=await api('/api/analyze',{csv:state.workflow.csv,variable:e.id,time:state.workflow.timeName,value:state.workflow.valueName,period:12,alpha:0.05,units:e.units||'',hamed_rao:true});
    state.analysis=d;
    $('#wf-method').textContent=d.method==='harmonic'?'POLYHARMONIC':(d.decomposed?'STL':'CLIMATOLOGY');
    wfLog(d.method==='harmonic'?'Polyharmonic fit selected for this greenhouse-gas record.':'STL decomposition selected for this record.');
    $('#wf-period').textContent=d.period;
    $('#wf-decomp').innerHTML=wfChart(d);
    wfStep(2,'done');wfStep(3,'active');
    $('#wf-trend').innerHTML=wfTrend(d);
    wfStep(3,'done');wfStep(4,'active');
    $('#wf-spatial').innerHTML=wfSpatial(d);
    loadGridDemo();
    wfStep(4,'done');
    wfStep(5,'active');
    renderWfReport(d);
    await loadFormulaPlan(e.id);
    wfStep(5,'done');
    $('#wf-io').innerHTML=`<span class="k">SEN</span><span class="v">${fmt(d.trend.sen.slope,5)} / step</span><span class="k">MK p</span><span class="v">${pval(d.trend.mann_kendall.p_value)}</span><span class="k">Z</span><span class="v">${fmt(d.trend.mann_kendall.z,3)}</span><span class="k">N</span><span class="v">${d.samples}</span>`;
    wfLog(`STL ${d.decomposed?'complete':'fallback'} (period ${d.period}); ${d.samples} samples.`);
    wfLog(`Mann-Kendall p=${pval(d.trend.mann_kendall.p_value)}; Sen=${fmt(d.trend.sen.slope,4)} ${e.units||''}/step.`);
    if(!d.decomposed)wfLog('STL skipped: fewer than two seasonal cycles. Climatology fallback used.','warn');
    if(d.warnings.length)d.warnings.forEach(w=>wfLog(w,'warn'));
    wfLog('Report assembled. Scalar map withheld; synthetic grid demo shown separately.','warn');
  }catch(err){toast(err.message);wfLog(err.message,'err');b.disabled=false;b.textContent='Retry data confirmation';}
  b.disabled=false;b.innerHTML='Data confirmed → remove seasonality';
}
$('#wf-confirm')&&$('#wf-confirm').addEventListener('click',runWorkflow);
$('#wf-report')&&$('#wf-report').addEventListener('click',async e=>{
  const act=e.target.closest('[data-wfact]');if(!act||!state.analysis)return;
  const d=state.analysis;
  if(act.dataset.wfact==='print'){window.print();return;}
  if(act.dataset.wfact==='copy'){await navigator.clipboard.writeText(JSON.stringify(d,null,2));toast('Analysis JSON copied.');return;}
  const rows=['time,observed,trend,seasonal,remainder'];
  d.series.time.forEach((t,i)=>rows.push([t,d.series.observed[i],d.components.trend[i],d.components.seasonal[i],d.components.remainder[i]].join(',')));
  const blob=new Blob([rows.join('\n')],{type:'text/csv'});const a=document.createElement('a');a.href=URL.createObjectURL(blob);a.download='trend-detective.csv';a.click();URL.revokeObjectURL(a.href);
});
function drawHeroChart(text){
  try{
    const t=parseCsv(text);const vi=t.header.length>1?1:0;
    const values=t.rows.map(r=>Number(r[vi])).filter(v=>isFinite(v)).slice(-160);
    $('#hero-chart').innerHTML=chart([{name:'sst',values,color:'#ffe08a'}],{height:210,xlabels:['start','','now']});
  }catch(e){}
}

/* ---------------- NASA DATA ---------------- */
async function loadNasa(){
  try{
    const d=await api('/api/nasa/datasets');state.nasa=d;
    if($('#stat-nasa'))$('#stat-nasa').textContent=d.datasets.length;
    $('#nasa-ds').innerHTML=d.datasets.map(x=>`<div class="dscard"><div class="id">${esc(x.key)}</div><h3>${esc(x.short_name)}</h3><div class="meta">${esc(x.provider)} / ${esc(x.mode)}</div><p>${esc(x.description)}</p><div class="meta">VAR ${esc(x.variable)} ${x.units?'['+esc(x.units)+']':''} / V ${esc(x.version||'latest')}</div><div class="actions" style="margin-top:6px"><button class="btn" data-plan="${esc(x.key)}">Plan fetch</button></div></div>`).join('');
    $('#nasa-key').innerHTML=d.datasets.map(x=>`<option value="${esc(x.key)}">${esc(x.key)} / ${esc(x.short_name)}</option>`).join('');
    $('#nasa-files').innerHTML=(d.files||[]).map(f=>`<span class="tag">${esc(f)}</span>`).join('')||'<span class="tag">none bundled</span>';
  }catch(e){$('#nasa-ds').innerHTML=`<div class="empty">${esc(e.message)}</div>`;}
}
$('#nasa-ds')&&$('#nasa-ds').addEventListener('click',e=>{const b=e.target.closest('[data-plan]');if(b){$('#nasa-key').value=b.dataset.plan;planNasa();}});
async function planNasa(){
  const out=$('#nasa-output');out.textContent='Planning…';
  try{
    const d=await api('/api/nasa/plan',{key:$('#nasa-key').value,start:$('#nasa-start').value,end:$('#nasa-end').value,max:+$('#nasa-max').value});
    out.textContent=`$ ${d.command}\n\n[exit ${d.exit_code}]\n${d.output}`;
    out.scrollIntoView({block:'nearest'});
  }catch(e){out.textContent='ERROR: '+e.message;out.scrollIntoView({block:'nearest'});}
}
$('#nasa-plan')&&$('#nasa-plan').addEventListener('click',planNasa);
async function searchNasa(){
  const out=$('#nasa-output');out.textContent='Querying CMR API…';
  try{
    const d=await api('/api/nasa/search',{key:$('#nasa-key').value,start:$('#nasa-start').value,end:$('#nasa-end').value,max:+$('#nasa-max').value});
    out.textContent=d.output;
  }catch(e){out.textContent='CMR query failed: '+e.message;}
}
$('#nasa-search')&&$('#nasa-search').addEventListener('click',searchNasa);

/* ---------------- FORMULA REGISTRY ---------------- */
async function loadFormulas(){
  try{
    const d=await api('/api/formulas');state.formulas=d;
    if($('#foot-formulas'))$('#foot-formulas').textContent=d.count;
    if($('#formula-count'))$('#formula-count').textContent=d.count+' PUBLIC FUNCTIONS';
    if($('#stat-functions'))$('#stat-functions').textContent=d.count;
    renderFormulas();
  }catch(e){$('#formula-list').innerHTML=`<div class="empty">${esc(e.message)}</div>`;}
}
const FORMULA_GROUPS=[
  ['01 · What is changing — observation, retrieval and Earth-system variables',['retrieval.hpp','aerosol.hpp','ozone.hpp','precip.hpp','cryo.hpp','ocean.hpp','ocean2.hpp','land.hpp','atmos.hpp','earth.hpp']],
  ['02 · How much is it changing — trend estimators, change detection and cycles',['series.hpp','stl.hpp','trend.hpp','multitrend.hpp','change.hpp','detect.hpp','harmonic.hpp','fft.hpp','ssa.hpp','loess.hpp','regime.hpp','hmm.hpp']],
  ['03 · Where is it changing — space, region, coherence and attribution',['spatial.hpp','geostat.hpp','regional.hpp','attrib.hpp','attrib2.hpp','causal.hpp','causal2.hpp','clustering.hpp','coherence.hpp','flow.hpp','landuse.hpp']],
  ['04 · Is it significant — uncertainty, power, validation and ensembles',['uncertainty.hpp','sstunc.hpp','validation.hpp','metrology.hpp','homogen.hpp','intercal.hpp','ensemble.hpp','quantile.hpp','relevance.hpp','power.hpp']],
  ['05 · Why it changes — physics, dynamics and transport',['physics.hpp','pde.hpp','assim.hpp','tempac.hpp','magnitude.hpp','maturity.hpp','copula.hpp','restrand.hpp']],
  ['06 · Numerics behind the answers — data, statistics and linear algebra',['json.hpp','csv.hpp','netcdf.hpp','readers.hpp','stats.hpp','linalg.hpp','ml.hpp','kernel.hpp']]
];
let formulaGroup='';
function formulaChips(){
  const box=$('#formula-cats');if(!box)return;
  const chip=(value,label)=>`<button class="tag ${formulaGroup===value?'amber':''}" data-fcat="${value}">${esc(label)}</button>`;
  box.innerHTML=[chip('','All questions')].concat(
    FORMULA_GROUPS.map(([label],i)=>chip(String(i),label.split(' · ')[0]+' · '+label.split(' — ')[1]))
  ).join('')+chip('leftover','Supporting helpers');
}
function renderFormulas(){
  if(!state.formulas)return;
  formulaChips();
  const q=($('#formula-search').value||'').toLowerCase();
  const byFile={};
  state.formulas.modules.forEach(m=>{byFile[m.file]=(m.functions||[]).filter(f=>!q||f.toLowerCase().includes(q));});
  const used=new Set();
  let html='';
  const block=(label,files,index)=>{
    const names=[];
    files.forEach(f=>{if(byFile[f]&&byFile[f].length){names.push(...byFile[f]);used.add(f);}});
    if(!names.length)return '';
    const open=(q||formulaGroup===String(index)||(!q&&!formulaGroup&&index===0))?' open':'';
    return `<details class="fitem"${open}><summary class="fhead"><span>${esc(label)}</span><span class="cnt">${names.length} calculations</span></summary><div class="fbody">${names.map(f=>`<div class="fsig">${esc(f)}</div>`).join('')}</div></details>`;
  };
  FORMULA_GROUPS.forEach(([label,files],i)=>{if(formulaGroup&&formulaGroup!==String(i)&&formulaGroup!=='leftover')return;html+=block(label,files,i);});
  if(!formulaGroup||formulaGroup==='leftover'){
    const leftovers=[];
    state.formulas.modules.forEach(m=>{if(!used.has(m.file)&&byFile[m.file]&&byFile[m.file].length)leftovers.push(...byFile[m.file]);});
    if(leftovers.length)html+=`<details class="fitem" open><summary class="fhead"><span>07 · Supporting helpers</span><span class="cnt">${leftovers.length} calculations</span></summary><div class="fbody">${leftovers.map(f=>`<div class="fsig">${esc(f)}</div>`).join('')}</div></details>`;
  }
  $('#formula-list').innerHTML=html||'<div class="empty">No calculation matches that filter.</div>';
}
$('#formula-cats')&&$('#formula-cats').addEventListener('click',e=>{
  const chip=e.target.closest('[data-fcat]');
  if(!chip)return;
  formulaGroup=chip.dataset.fcat;
  renderFormulas();
  window.scrollTo({top:0,behavior:'smooth'});
});
$('#formula-search')&&$('#formula-search').addEventListener('input',renderFormulas);
$('#formula-refresh')&&$('#formula-refresh').addEventListener('click',loadFormulas);

/* ---------------- COMMAND CONSOLE ---------------- */
async function loadHelp(){
  try{const d=await api('/api/help');state.help=d.text;renderHelp();}catch(e){$('#help-text').textContent=e.message;}
}
function renderHelp(){
  const q=($('#help-search').value||'').toLowerCase();
  const lines=state.help.split('\n');
  $('#help-text').textContent=q?lines.filter(l=>l.toLowerCase().includes(q)).join('\n'):state.help;
}
$('#help-search')&&$('#help-search').addEventListener('input',renderHelp);
function consoleWrite(html){const el=$('#cmd-output');el.innerHTML=html;el.scrollTop=el.scrollHeight;}
async function runConsole(){
  const cmd=$('#cmd-input').value.trim();if(!cmd)return;
  consoleWrite(`<span class="cmd">$ climate ${esc(cmd)}</span>\n`);
  try{
    const d=await api('/api/run',{command:cmd});
    consoleWrite(`<span class="cmd">$ climate ${esc(cmd)}</span>\n${esc(d.output)}\n<span class="${d.exit_code?'bad':'ok'}">[exit ${d.exit_code}]</span>`);
  }catch(e){consoleWrite(`<span class="cmd">$ climate ${esc(cmd)}</span>\n<span class="bad">${esc(e.message)}</span>`);}
}
async function runBatch(){
  const lines=$('#cmd-batch').value.split(/\r?\n/).map(s=>s.trim()).filter(Boolean);
  if(!lines.length){toast('Enter at least one command');return;}
  let out='';
  for(const line of lines){
    out+=`<span class="cmd">$ climate ${esc(line)}</span>\n`;
    try{
      const d=await api('/api/run',{command:line});
      out+=`${esc(d.output)}\n<span class="${d.exit_code?'bad':'ok'}">[exit ${d.exit_code}]</span>\n\n`;
    }catch(e){out+=`<span class="bad">${esc(e.message)}</span>\n\n`;}
  }
  consoleWrite(out);
}
$('#cmd-run')&&$('#cmd-run').addEventListener('click',runConsole);
$('#cmd-batch-run')&&$('#cmd-batch-run').addEventListener('click',runBatch);
$('#cmd-input')&&$('#cmd-input').addEventListener('keydown',e=>{if(e.key==='Enter')runConsole();});
$$('.quick button').forEach(b=>b.addEventListener('click',()=>{$('#cmd-input').value=b.dataset.cmd;runConsole();}));

/* ---------------- REGIONAL CONTRAST ---------------- */
function cxSpark(values,color){
  if(!values||values.length<2)return '';
  const w=260,h=54,pad=3;
  const min=Math.min(...values),max=Math.max(...values),span=(max-min)||1;
  const pts=values.map((v,i)=>[pad+(w-2*pad)*i/(values.length-1||1),h-pad-(h-2*pad)*(v-min)/span]);
  return `<svg viewBox="0 0 ${w} ${h}" preserveAspectRatio="none" style="width:100%;height:54px;display:block;background:#120c08;border:1px solid var(--line)"><polyline fill="none" stroke="${color}" stroke-width="1.4" points="${pts.map(p=>p[0].toFixed(1)+','+p[1].toFixed(1)).join(' ')}"/></svg>`;
}
async function runContrast(){
  const out=$('#cx-out');if(!out)return;
  const t0=Date.now();
  out.innerHTML='<div class="empty">Calling the NASA POWER API for both regions… <span id="cx-tick">0s</span></div>';
  const tick=setInterval(()=>{const el=document.getElementById('cx-tick');if(el)el.textContent=Math.round((Date.now()-t0)/1000)+'s';},500);
  const body={parameter:$('#cx-param').value,start:+$('#cx-start').value,end:+$('#cx-end').value,
    points:[{lat:+$('#cx-lat1').value,lon:+$('#cx-lon1').value,label:'Region A · lat '+$('#cx-lat1').value+', lon '+$('#cx-lon1').value},
            {lat:+$('#cx-lat2').value,lon:+$('#cx-lon2').value,label:'Region B · lat '+$('#cx-lat2').value+', lon '+$('#cx-lon2').value}]};
  try{
    const d=await api('/api/contrast',body);
    const cards=d.regions.map(r=>{
      const rise=r.slope_per_year>0;
      const curve=(d.curves&&d.curves[r.label]||[]).map(p=>p[1]);
      const sig=r.significant?'<span class="cxsig yes">SIGNIFICANT</span>':'<span class="cxsig no">NOT SIGNIFICANT</span>';
      return `<div class="cxcard ${rise?'rise':'fall'}">
        <div class="cxhead"><b>${esc(r.label)}</b>${sig}</div>
        <div class="cxarrow">${rise?'▲ RISING':'▼ FALLING'} <span>${esc(r.direction)}</span></div>
        <div class="cxstats">
          <div><span>Sen slope</span><b>${fmt(r.slope_per_year,4)} ${esc(d.units)}/yr</b></div>
          <div><span>Change ${fmt(r.start,0)}–${fmt(r.end,0)}</span><b>${r.change>=0?'+':''}${fmt(r.change,2)} ${esc(d.units)}</b></div>
          <div><span>Relative</span><b>${r.percent>=0?'+':''}${fmt(r.percent,1)}%</b></div>
          <div><span>Seasonal MK p</span><b>${pval(r.p_value)}</b></div>
          <div><span>Samples</span><b>${r.samples}</b></div>
          <div><span>Mean</span><b>${fmt(r.mean,2)} ${esc(d.units)}</b></div>
        </div>
        ${cxSpark(curve,rise?'#ffab54':'#ff8f86')}
        <div class="note">${esc(r.file)}</div>
      </div>`;
    }).join('');
    out.innerHTML=`<div class="cxverdict ${d.opposed?'opposed':''}">${esc(d.verdict)}</div>
      <div class="cxgrid">${cards}</div>
      <div class="note">Both regions were pulled live from the NASA POWER API for ${esc(d.parameter)}. Each verdict comes from its own seasonal Mann–Kendall test at alpha = ${d.alpha}; "significant" means a monotonic change this large is unlikely when there is no trend.</div>`;
    clearInterval(tick);
  }catch(e){
    clearInterval(tick);
    const message = typeof e==='string' ? e : (e&&e.message)||JSON.stringify(e);
    out.innerHTML=`<div class="panel"><div class="panelhead"><b>Contrast could not run</b><span class="hint">live NASA POWER</span></div><div class="panelbody">
      <div class="err">${esc(message)}</div>
      <div class="note">This panel pulls two live NASA POWER series, so it needs this machine to reach <span class="telemetry">power.larc.nasa.gov</span>. POWER needs no login — only Earthdata datasets do. If the network is blocked, open Fetch data to test the connection, then come back and retry.</div>
      <div class="row"><button class="btn primary" id="cx-retry">Retry the comparison</button><button class="btn" data-goto="live">Open Fetch data</button><button class="btn" data-goto="significance">Open Significance</button></div>
    </div></div>`;
    const retry=document.getElementById('cx-retry');
    if(retry)retry.addEventListener('click',runContrast);
  }
}
$('#cx-run')&&$('#cx-run').addEventListener('click',runContrast);

/* ---------------- SIGNIFICANCE ---------------- */
function sigOptions(){
  const box=$('#sig-var');if(!box)return;
  const list=state.catalog.length?state.catalog:[];
  box.innerHTML=list.map(e=>`<option value="${esc(e.id)}">${esc(e.id)} · ${esc(e.name)}</option>`).join('');
}
function mdtPerYear(d){
  const conv=d.spacing>0&&d.spacing<0.5?12:1;               // steps per year
  const sd=d.summary&&isFinite(d.summary.sd)?d.summary.sd:NaN;
  const n=d.samples;
  const spanYears=(d.end_time-d.start_time)/conv;
  if(!isFinite(sd)||!n||spanYears<=0)return null;
  return 2.802*sd*Math.sqrt(12/n)/spanYears;                // 80% power, alpha 0.05
}
async function sigRun(){
  const out=$('#sig-out');if(!out)return;
  const id=$('#sig-var').value;if(!id)return toast('Pick a variable first.');
  const meta=wfMeta(id),file=meta.file||(id+'_monthly.csv');
  out.innerHTML='<div class="panel"><div class="panelbody"><div class="hint">Loading the series and running the tests…</div></div></div>';
  try{
    const text=await (await fetch('/api/data/file?name='+encodeURIComponent(file))).text();
    const alpha=Math.min(0.2,Math.max(0.001,Number($('#sig-alpha').value)||0.05));
    const entry=state.catalog.find(e=>e.id===id)||{};
    const d=await api('/api/analyze',{csv:text,variable:id,alpha,period:12,units:entry.units||''});
    const mk=d.trend.mann_kendall,raw=d.raw_trend.mann_kendall,sen=d.trend.sen,smk=d.seasonal_mann_kendall;
    const ens=d.trend.ensemble||{members:[],agreement:0,consensus:''};
    const mdt=mdtPerYear(d);
    const rows=[
      ['Is there a monotonic trend at all?','Mann–Kendall on the deseasonalised series',`S = ${fmt(mk.s,0)}, z = ${fmt(mk.z,3)}`,pval(mk.p_value),mk.significant?'yes':'no'],
      ['Does the raw record already show it?','Mann–Kendall without removing the trend first',`tau = ${fmt(raw.tau,3)}`,pval(raw.p_value),raw.significant?'yes':'no'],
      ['Is the seasonal cycle hiding it?','Seasonal Mann–Kendall, month by month',smk?`z = ${fmt(smk.z,3)}`:'not available',smk?pval(smk.p_value):'—',smk?(smk.significant?'yes':'no'):'—'],
      ['Is autocorrelation inflating the test?','Hamed–Rao variance correction',`${fmt(mk.correction_factor,3)}x variance`,pval(mk.p_value),mk.correction_factor>1.2?'penalised':'no inflation'],
      ['How large is the change, with uncertainty?','Sen slope and its rank-based 95% interval',`${fmt(sen.slope,4)} [${fmt(sen.lower,4)}, ${fmt(sen.upper,4)}]`,'—',sen.lower>0?'firmly rising':sen.upper<0?'firmly falling':'straddles zero'],
      ['Do the estimators agree on the sign?','Accuracy-weighted ensemble',`${Math.round((ens.agreement||0)*100)}% agreement`,pval(mk.p_value),ens.consensus||'—'],
      ['Could a trend this size be detected at all?','Minimum detectable trend, 80% power',mdt===null?'—':`${fmt(mdt,4)} ${esc(d.units||'')}/yr`,'—',mdt!==null&&Math.abs(sen.slope*12)>mdt?'record resolves it':'record is short'],
    ];
    out.innerHTML=`<div class="panel"><div class="panelhead"><b>${esc(id)} · significance suite</b><span class="hint">alpha = ${fmt(alpha,3)} · ${d.samples} samples</span></div><div class="panelbody">
      <div class="verdict">${mk.significant?`<b>SIGNIFICANT ${esc(mk.direction).toUpperCase()}</b>`:'<b>NO SIGNIFICANT MONOTONIC TREND</b>'} at alpha = ${fmt(alpha,3)}. ${mk.significant?'The change survives autocorrelation, the seasonal cycle and estimator disagreement.':'The change cannot be separated from noise at this alpha, so it is not reported as a trend.'}</div>
      <table class="etable"><thead><tr><th>Question</th><th>Test</th><th>Statistic</th><th>p</th><th>Verdict</th></tr></thead><tbody>
        ${rows.map(r=>`<tr><td>${esc(r[0])}</td><td class="mut">${esc(r[1])}</td><td>${r[2]}</td><td>${r[3]}</td><td class="${/yes|firmly|resolves|penalised|no inflation/.test(r[4])?'good':'warn'}">${esc(r[4])}</td></tr>`).join('')}
      </tbody></table>
      <div class="kv"><span class="k">RANGE</span><span class="v">${fmt(d.start_time,2)} → ${fmt(d.end_time,2)}</span><span class="k">MISSING FILLED</span><span class="v">${d.missing_filled}</span><span class="k">ALPHA</span><span class="v">${fmt(d.alpha,3)}</span></div>
      <div class="note">Every row is computed by the C++ engine from the same series. The last row is the smallest trend this record could have detected, so a "no" above can be read as "no detectable trend" instead of "no change".</div>
    </div></div>`;
  }catch(e){
    out.innerHTML=`<div class="panel"><div class="panelbody"><div class="err">${esc(e.message)}</div></div></div>`;
  }
}
$('#sig-run')&&$('#sig-run').addEventListener('click',sigRun);

/*__LAB_JS__*/
/*__VIS_JS__*/
/*__LIVE_JS__*/

/* ---------------- RETRIEVAL RECORD ---------------- */
async function retRun(){
  const out=$('#ret-out');if(!out)return;
  const [file,column,units]=$('#ret-series').value.split('|');
  out.innerHTML='<div class="panel"><div class="panelbody"><div class="hint">Loading the record…</div></div></div>';
  try{
    const text=await (await fetch('/api/data/file?name='+encodeURIComponent(file))).text();
    const alpha=Math.min(0.2,Math.max(0.001,Number($('#ret-alpha').value)||0.05));
    const d=await api('/api/analyze',{csv:text,value:column,alpha,period:12,units});
    const labels=d.series.time.map(t=>fmt(t,1));
    out.innerHTML=`<div class="panel"><div class="panelhead"><b>${esc(file)}</b><span class="hint">${d.samples} samples · ${esc(units)}</span></div><div class="panelbody">
      ${chartCard('Observed and STL trend',[{name:'observed',values:d.series.observed,color:'#f7e8c8'},{name:'trend',values:d.components.trend,color:'#ffab54'}],{xlabels:labels,height:250})}
      <div class="cards">${metric('Sen slope',fmt(d.trend.sen.slope,4),esc(units)+' / step')}${metric('Per decade',fmt(d.trend.sen_per_decade,3),esc(units)+' / decade')}${metric('Mann–Kendall p',pval(d.trend.mann_kendall.p_value),d.trend.mann_kendall.direction,d.trend.mann_kendall.significant?'good':'warn')}${metric('Change over record',fmt(d.trend.sen.slope*(d.end_time-d.start_time),2),esc(units))}</div>
      <div class="verdict">${d.trend.mann_kendall.significant?`<b>SIGNIFICANT ${esc(d.trend.mann_kendall.direction.toUpperCase())}</b>`:'<b>NO SIGNIFICANT MONOTONIC TREND</b>'} at alpha = ${fmt(d.alpha,3)}.</div>
    </div></div>`;
  }catch(e){
    const message=typeof e==='string'?e:(e&&e.message)||JSON.stringify(e);
    out.innerHTML=`<div class="panel"><div class="panelbody"><div class="err">${esc(message)}</div></div></div>`;
  }
}
$('#ret-run')&&$('#ret-run').addEventListener('click',retRun);

/* ---------------- COACH MARKS ---------------- */
const COACH_KEY='trend-detective-tour-v1';
const COACH_STEPS=[
  {sel:'[data-nav="live"]',text:'Start here. Every series is pulled live from a NASA endpoint the moment you press the button — GISTEMP, POWER, NSIDC, or your own Earthdata dataset.'},
  {sel:'[data-nav="catalog"]',text:'Each variable card states what kind of change it measures, so you can choose the variable before choosing a method.'},
  {sel:'[data-nav="workflow"]',text:'Here the same record answers the four questions: what is changing, where, how much, and whether it is significant.'},
  {sel:'#cx-run',text:'Compare one variable in two places. Two live POWER series, two Sen slopes, two Mann–Kendall verdicts — an opposition becomes visible instead of averaged away.'},
  {sel:'[data-nav="significance"]',text:'Only changes that survive autocorrelation, the seasonal cycle and estimator disagreement are reported as significant.'}
];
let coachIndex=0,coachTarget=null;
function coachShow(index){
  const box=$('#coach');if(!box)return;
  if(index>=COACH_STEPS.length){coachEnd();return;}
  coachIndex=index;
  const step=COACH_STEPS[index];
  if(step.sel==='#cx-run'&&$('.page.active').id!=='page-workflow')nav('workflow');
  const target=document.querySelector(step.sel);
  if(!target){coachShow(index+1);return;}
  if(coachTarget)coachTarget.classList.remove('coach-focus');
  coachTarget=target;target.classList.add('coach-focus');
  const rect=target.getBoundingClientRect();
  box.hidden=false;
  $('#coach-step').textContent=(index+1)+' / '+COACH_STEPS.length;
  $('#coach-text').textContent=step.text;
  const boxW=320;
  let left=rect.right+16, arrowSide='left';
  if(left+boxW>window.innerWidth-16){left=Math.max(16,rect.left-boxW-16);arrowSide='right';}
  let top=Math.min(Math.max(12,rect.top-8),window.innerHeight-190);
  box.style.left=left+'px';box.style.top=top+'px';
  const arrow=$('#coach-arrow');
  arrow.style.top=Math.max(10,rect.top-top+rect.height/2-6)+'px';
  arrow.style.left=arrowSide==='left'?'-6px':'auto';
  arrow.style.right=arrowSide==='right'?'-6px':'auto';
  arrow.style.transform='rotate('+(arrowSide==='left'?'45deg':'-135deg')+')';
}
function coachEnd(){
  const box=$('#coach');if(box)box.hidden=true;
  if(coachTarget)coachTarget.classList.remove('coach-focus');
  try{localStorage.setItem(COACH_KEY,'done');}catch(e){}
}
function coachStart(){
  let seen=null;try{seen=localStorage.getItem(COACH_KEY);}catch(e){}
  if(seen)return;
  setTimeout(()=>coachShow(0),900);
}
$('#coach-next')&&$('#coach-next').addEventListener('click',()=>coachShow(coachIndex+1));
$('#coach-skip')&&$('#coach-skip').addEventListener('click',coachEnd);
window.addEventListener('resize',()=>{if(!$('#coach').hidden)coachShow(coachIndex);});
coachStart();

/* ---------------- BOOT ---------------- */
renderPhysics();
loadOverview();
loadNasa();
loadFormulas();
loadHelp();
loadLab();
loadVisual();
loadLive();
</script>
</body>
</html>)CLIMATE_UI";
        const std::string marker = "/*__LAB_JS__*/";
        const std::size_t at = html.find(marker);
        if (at != std::string::npos) html.replace(at, marker.size(), climate_lab_js());
        const std::string vis_marker = "/*__VIS_JS__*/";
        const std::size_t vis_at = html.find(vis_marker);
        if (vis_at != std::string::npos)
            html.replace(vis_at, vis_marker.size(), climate_visual_js());
        const std::string live_marker = "/*__LIVE_JS__*/";
        const std::size_t live_at = html.find(live_marker);
        if (live_at != std::string::npos)
            html.replace(live_at, live_marker.size(), climate_live_js());
        return html;
    }();
    return page.c_str();
}

}  // namespace climate
