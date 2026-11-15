#include "ui_page.hpp"

namespace climate {

const char* climate_ui_html() {
    return R"CLIMATE_UI(<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Climate Workbench</title>
<style>
:root{
  --bg:#0b0806; --bg2:#140e09; --panel:#1b130d; --panel2:#241a12;
  --line:#4a3524; --line2:#6b4a30; --text:#e8d6b8; --muted:#a68b6c;
  --teal:#d8a657; --sky:#c98b4b; --amber:#e0b15f; --rose:#c96a4a;
  --green:#a9b87a; --violet:#b58a6a; --bezel:#b7a184;
  --series-obs:#f0d9b5; --series-trend:#ff9a4d; --series-season:#5fc2c2; --series-rem:#d98a8a;
  --series-alt:#b8c97a; --series-violet:#c3a6e0;
  --mono:"Courier New",Courier,"Liberation Mono",Menlo,monospace;
  --sans:"Courier New",Courier,"Liberation Mono",Menlo,monospace;
}
*{box-sizing:border-box;border-radius:0!important}
html,body{margin:0;padding:0}
body{background:#080605;color:var(--text);font-family:var(--mono);font-size:13px;line-height:1.5;
  padding:18px;text-shadow:0 0 7px rgba(216,166,87,.10);-webkit-font-smoothing:antialiased}
body:before{content:"";position:fixed;inset:0;pointer-events:none;z-index:9998;
  box-shadow:inset 0 0 190px rgba(0,0,0,.86)}
body:after{content:"";position:fixed;inset:0;pointer-events:none;z-index:9999;opacity:.16;
  background-image:url("data:image/svg+xml,%3Csvg xmlns='http://www.w3.org/2000/svg' width='4' height='4'%3E%3Crect width='4' height='1' fill='%23000'/%3E%3C/svg%3E")}
button,input,select,textarea{font:inherit;color:inherit}
button{cursor:pointer}
a{color:var(--teal)}
::selection{background:#6b4a30;color:#fff3dc}
.app{display:grid;grid-template-columns:228px minmax(0,1fr);min-height:calc(100vh - 36px);
  border:10px solid var(--bezel);outline:2px solid #54402e;
  box-shadow:inset 0 0 90px rgba(0,0,0,.78),0 0 0 1px #000;background:var(--bg)}
.sidebar{position:sticky;top:0;height:calc(100vh - 56px);display:flex;flex-direction:column;
  background:#1a120c;border-right:2px solid var(--line);padding:16px 12px 10px;gap:3px;z-index:5}
.brand{display:flex;align-items:center;gap:9px;padding:0 5px 16px;border-bottom:1px dashed var(--line2);margin-bottom:8px}
.mark{width:34px;height:34px;border:1px solid var(--bezel);background:#2b1e14;position:relative;flex:none;
  display:grid;place-items:center;color:var(--teal);font-size:15px;font-weight:700;
  box-shadow:inset 0 0 12px rgba(216,166,87,.16)}
.mark:after{content:"▚";position:absolute;font-size:20px;color:#d8c3a0}
.brand b{display:block;font-size:12px;letter-spacing:.18em;text-transform:uppercase;color:#f0dfc2}
.brand span{display:block;font-size:9.5px;color:var(--muted);letter-spacing:.12em;margin-top:1px}
.navlabel{font-size:9px;letter-spacing:.2em;color:#8a6c4f;text-transform:uppercase;padding:14px 6px 5px}
.nav{display:flex;align-items:center;gap:8px;width:100%;border:1px solid transparent;
  background:transparent;color:#c7ab87;padding:7px 8px;text-align:left;font-size:11.5px;
  letter-spacing:.04em;transition:background .1s,color .1s}
.nav:hover{background:#2a1e14;color:#f3e3c6;border-color:#3c2b1d}
.nav.active{background:#2b1e14;border-color:var(--line2);border-left:3px solid var(--teal);color:#ffe9c2}
.nav .ico{width:22px;height:22px;display:grid;place-items:center;flex:none;
  font-size:9px;font-weight:700;color:#c9ad86;background:#241a12;border:1px solid #4a3524}
.nav.active .ico{color:#120d09;background:var(--teal);border-color:var(--teal)}
.sidefoot{margin-top:auto;padding:12px 5px 2px;border-top:1px dashed var(--line);font-size:9.5px;color:#8a6c4f;line-height:1.7}
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
.hero:before{content:"";position:absolute;left:0;top:0;bottom:0;width:3px;background:var(--teal)}
.eyebrow{font-size:9.5px;letter-spacing:.22em;text-transform:uppercase;color:var(--teal);font-weight:700}
h1{font-size:clamp(24px,3vw,38px);line-height:1.1;letter-spacing:-.02em;margin:8px 0 10px;max-width:760px;font-weight:700;
  text-transform:uppercase}
h1 em{font-style:normal;color:var(--teal)}
.lede{color:#bfa582;max-width:650px;font-size:12.5px;margin:0 0 18px}
.actions{display:flex;gap:8px;flex-wrap:wrap}
.btn{border:1px solid var(--line2);background:#241a12;color:var(--text);padding:8px 13px;
  font-size:10.5px;font-weight:700;letter-spacing:.09em;text-transform:uppercase;transition:.1s;
  display:inline-flex;align-items:center;gap:7px}
.btn:hover{border-color:var(--bezel);background:#332419;color:#fff0d4}
.btn.primary{background:var(--teal);border-color:var(--teal);color:#160f0a}
.btn.primary:hover{background:#f0c078;border-color:#f0c078}
.btn.ghost{background:transparent}
.btn:disabled{opacity:.45;cursor:not-allowed}
.stats{display:grid;grid-template-columns:repeat(4,minmax(0,1fr));gap:1px;margin:14px 0 24px;
  border:1px solid var(--line);background:var(--line)}
.stat{background:#1a120c;padding:12px 14px}
.stat .k{font-size:9px;letter-spacing:.14em;text-transform:uppercase;color:#9c7c59}
.stat .v{font-size:23px;font-weight:700;margin-top:3px;font-variant-numeric:tabular-nums;color:var(--teal)}
.stat .s{font-size:10px;color:#8e7153;margin-top:1px}
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
.module .arrow{position:absolute;right:12px;top:12px;color:#7a5c40;font-size:14px}
.tool{display:grid;grid-template-columns:350px minmax(0,1fr);gap:10px;align-items:start}
.panel{border:1px solid var(--line);background:#18110c;overflow:hidden}
.panelhead{padding:9px 12px;border-bottom:1px solid var(--line2);display:flex;align-items:center;justify-content:space-between;gap:10px;background:#22180f}
.panelhead b{font-size:10.5px;letter-spacing:.1em;text-transform:uppercase;color:#f0dfc2}
.panelhead .hint{font-size:9px;color:#947554;letter-spacing:.06em;text-transform:uppercase}
.panelbody{padding:12px}
.field{margin-bottom:9px}
.field:last-child{margin-bottom:0}
.field label,.label{display:block;font-size:9px;letter-spacing:.12em;text-transform:uppercase;color:#a38562;margin-bottom:4px;font-weight:700}
input[type=text],input[type=number],select,textarea{width:100%;background:#0e0a07;border:1px solid #4a3524;
  padding:7px 8px;font-size:11.5px;outline:none;color:#ead8ba}
input:focus,select:focus,textarea:focus{border-color:var(--teal);box-shadow:0 0 0 1px rgba(216,166,87,.25)}
textarea{font-size:10.5px;line-height:1.5;resize:vertical;min-height:96px}
select{appearance:none;background-image:url("data:image/svg+xml,%3Csvg xmlns='http://www.w3.org/2000/svg' width='8' height='6'%3E%3Cpath d='M0 0h8L4 6z' fill='%23a38562'/%3E%3C/svg%3E");
  background-position:calc(100% - 10px) center;background-repeat:no-repeat;padding-right:24px}
.row{display:grid;grid-template-columns:1fr 1fr;gap:8px}
.row3{display:grid;grid-template-columns:repeat(3,1fr);gap:8px}
.check{display:flex;align-items:center;gap:7px;font-size:10.5px;color:#bfa582;margin:6px 0;cursor:pointer}
.check input{accent-color:var(--teal);width:13px;height:13px}
.run{width:100%;justify-content:center;margin-top:5px}
.result{min-height:180px}
.placeholder{min-height:160px;display:grid;place-items:center;text-align:center;color:#8e7153;padding:24px}
.placeholder .big{font-size:26px;color:#5a4330;margin-bottom:6px}
.placeholder p{max-width:340px;font-size:10.5px;margin:0}
.cards{display:grid;grid-template-columns:repeat(4,minmax(0,1fr));gap:1px;margin-bottom:10px;border:1px solid var(--line);background:var(--line)}
.metric{background:#1a120c;padding:8px 9px}
.metric .k{font-size:8.5px;text-transform:uppercase;letter-spacing:.12em;color:#9c7c59}
.metric .v{font-size:17px;font-weight:700;margin-top:2px;font-variant-numeric:tabular-nums;color:#eedcbe}
.metric .v.good{color:var(--green)} .metric .v.warn{color:var(--amber)} .metric .v.bad{color:var(--rose)}
.metric .s{font-size:9.5px;color:#8e7153;margin-top:1px}
.verdict{border:1px solid var(--line2);border-left:3px solid var(--teal);background:#22180f;
  padding:8px 10px;font-size:11px;margin-bottom:10px;color:#d8c3a2}
.verdict b{color:var(--teal)}
.chartwrap{border:1px solid var(--line);background:#0d0906;padding:5px 4px 2px;margin-bottom:8px}
.chart{width:100%;height:auto;display:block;overflow:visible}
.chart text{fill:#947554;font-family:var(--mono);font-size:9px}
.chart .gridline{stroke:#3a2a1d;stroke-width:1;stroke-dasharray:2 5}
.chart .axis{stroke:#4a3524;stroke-width:1}
.chart .line{fill:none;stroke-width:1.15;stroke-linejoin:round;stroke-linecap:round;stroke-opacity:.92}
.chart .area{stroke:none;fill:rgba(255,154,77,.06)}
.chart .legend{font-family:var(--mono);font-size:9px;fill:#b39570}
.etable{width:100%;border-collapse:collapse;font-family:var(--mono);font-size:9px}
.etable th{text-align:left;color:#8f7052;font-weight:400;letter-spacing:.08em;text-transform:uppercase;border-bottom:1px solid #3c2b1d;padding:4px 5px}
.etable td{padding:5px;border-bottom:1px dashed #2f2218;color:#cbb18d}
.etable tr.primary td{color:#ffe2b0;background:#24180e}
.etable .pick{color:#ffb168;font-weight:700}
.etable .fallback{color:#8f7052}
.minicharts{display:grid;grid-template-columns:1fr 1fr;gap:8px}
.minicharts .chartwrap{margin:0}
.tagrow{display:flex;gap:4px;flex-wrap:wrap;margin-top:6px}
.tag{display:inline-flex;align-items:center;border:1px solid var(--line2);background:#22180f;color:#c7ab87;
  padding:2px 6px;font-size:9px;letter-spacing:.04em;text-transform:uppercase}
.tag.teal{border-color:#7b5b32;color:#f0c078;background:#2a1d10}
.tag.amber{border-color:#6e4a2c;color:#e0b15f;background:#2a1d10}
.warnlist{margin:8px 0 0;padding:0;list-style:none}
.warnlist li{font-size:10.5px;color:#d8b77e;background:#241a12;border-left:2px solid var(--amber);
  padding:5px 8px;margin-top:4px}
details{margin-top:8px;border-top:1px dashed #3c2b1d;padding-top:7px}
summary{cursor:pointer;color:#a38562;font-size:10px;letter-spacing:.06em}
pre{font-family:var(--mono);font-size:10px;line-height:1.55;background:#0a0705;border:1px solid #3c2b1d;
  padding:9px;overflow:auto;max-height:300px;color:#cdb18a;white-space:pre-wrap;overflow-wrap:anywhere}
.tabs{display:flex;gap:2px;border-bottom:1px solid var(--line2);padding:0 2px;margin-bottom:12px;flex-wrap:wrap}
.tab{background:transparent;border:1px solid transparent;border-bottom:0;color:#a38562;padding:7px 10px;
  font-size:10px;font-weight:700;letter-spacing:.06em;text-transform:uppercase}
.tab:hover{color:var(--text);background:#22180f}
.tab.active{color:var(--teal);border-color:var(--line2);background:#22180f}
.tabpane{display:none}.tabpane.active{display:block}
.split{display:grid;grid-template-columns:1fr 1fr;gap:12px}
.bigvalue{font-size:31px;letter-spacing:-.03em;color:var(--teal);line-height:1;
  font-variant-numeric:tabular-nums;margin:2px 0 4px}
.unit{font-size:10px;color:var(--muted);text-transform:uppercase;letter-spacing:.08em}
.kv{display:grid;grid-template-columns:1fr auto;gap:5px 12px;font-size:10.5px;margin-top:9px}
.kv .k{color:#9c7c59}.kv .v{color:#d8c3a2;text-align:right}
.coeffs{display:grid;grid-template-columns:repeat(4,1fr);gap:4px;margin-top:8px}
.coeff{border:1px solid var(--line);padding:5px 6px;background:#140e09}
.coeff b{display:block;font-size:8.5px;color:#947554;font-weight:400}
.coeff span{font-size:10px;color:#d8c3a2}
.table{width:100%;border-collapse:collapse;font-size:10.5px}
.table th{text-align:left;font-size:8.5px;letter-spacing:.12em;text-transform:uppercase;color:#9c7c59;
  border-bottom:1px solid var(--line2);padding:6px;font-weight:700}
.table td{border-bottom:1px solid #2e2116;padding:6px;color:#c7ab87;vertical-align:top}
.table tr:hover td{background:#22180f}
.table td.id{color:var(--teal);white-space:nowrap}
.table td.mut{color:#a38562}
.searchrow{display:flex;gap:8px;margin-bottom:10px;flex-wrap:wrap}
.searchrow input{flex:1;min-width:160px}.searchrow select{width:auto;min-width:140px}
.empty{padding:28px;text-align:center;color:#8e7153;font-size:10.5px}
.spinner{display:inline-block;width:10px;height:10px;border:2px solid rgba(0,0,0,.3);
  border-top-color:#160f0a;animation:spin .7s linear infinite;vertical-align:-1px;margin-right:6px}
@keyframes spin{to{transform:rotate(360deg)}}
.toast{position:fixed;right:22px;bottom:22px;z-index:10000;max-width:380px;border:1px solid #8a4a32;
  background:#2d150f;color:#ffd8c8;padding:9px 12px;font-size:10.5px;box-shadow:0 0 24px rgba(0,0,0,.7);display:none}
.toast.show{display:block}
.formula{font-size:9.5px;color:#b39570;background:#0d0906;border:1px dashed #3c2b1d;padding:6px 8px;margin-top:8px;line-height:1.6}
.note{font-size:9.5px;color:#8e7153;margin-top:7px;line-height:1.6}
.hero .scope{border:1px solid var(--line2);background:#0a0705;padding:7px;min-height:230px;display:flex;flex-direction:column}
.scopehead{display:flex;justify-content:space-between;font-size:8.5px;letter-spacing:.12em;color:#947554;
  border-bottom:1px dashed #3c2b1d;padding-bottom:5px;text-transform:uppercase}
.scope svg{width:100%;height:auto;flex:1}
.scopefoot{font-size:8.5px;color:#7a5c40;border-top:1px dashed #3c2b1d;padding-top:4px;letter-spacing:.08em}
.workflow{display:grid;grid-template-columns:272px minmax(0,1fr) 290px;gap:10px;align-items:start}
.wleft,.wright{display:flex;flex-direction:column;gap:10px;position:sticky;top:56px}
.wcenter{display:flex;flex-direction:column;gap:8px;min-width:0}
.step{border:1px solid var(--line);background:#18110c}
.step-head{display:flex;align-items:center;gap:9px;padding:8px 10px;border-bottom:1px solid var(--line2);background:#22180f}
.step-no{width:24px;height:24px;display:grid;place-items:center;border:1px solid var(--line2);font-size:10px;
  font-weight:700;color:var(--teal);flex:none}
.step-title{font-size:10.5px;letter-spacing:.1em;text-transform:uppercase;font-weight:700;color:#f0dfc2}
.step-status{margin-left:auto;font-size:8.5px;letter-spacing:.09em;text-transform:uppercase;color:#947554}
.step-body{padding:11px}
.step:not(.open) .step-body{display:none}
.step.locked{opacity:.55}
.step.done .step-no{background:var(--teal);color:#160f0a}
.step.done .step-status{color:var(--green)}
.facts{display:grid;grid-template-columns:repeat(3,1fr);gap:1px;background:var(--line);border:1px solid var(--line)}
.fact{background:#140e09;padding:7px 8px}
.fact .k{font-size:8px;letter-spacing:.1em;text-transform:uppercase;color:#947554}
.fact .v{font-size:11px;color:#e0c9a6;margin-top:2px}
.fact .s{font-size:9px;color:#8e7153}
.bar{height:6px;background:#2a1e14;border:1px solid #3c2b1d;margin:3px 0 7px}
.bar i{display:block;height:100%;background:var(--teal)}
.bar.warn i{background:var(--rose)}
.log{max-height:250px;overflow:auto;font-size:9.5px;line-height:1.55}
.logline{border-bottom:1px dotted #3c2b1d;padding:3px 0;color:#b39570}
.logline b{color:var(--teal);font-weight:400}
.logline.warn{color:#e0b15f}
.logline.err{color:#e08a6a}
.artlist{margin:0;padding:0;list-style:none;font-size:9.5px;color:#c7ab87}
.artlist li{border-left:2px solid var(--rose);padding:3px 7px;margin-bottom:4px;background:#22130e}
.varpick{max-height:300px;overflow:auto;border:1px solid var(--line);background:#0e0a07}
.varitem{display:flex;gap:7px;align-items:flex-start;padding:6px 7px;border-bottom:1px solid #2a1e14;cursor:pointer}
.varitem:hover{background:#22180f}
.varitem.active{background:#2b1e14;border-left:2px solid var(--teal)}
.varitem .code{font-size:8.5px;color:var(--teal);border:1px solid var(--line2);padding:1px 3px;flex:none;margin-top:1px}
.varitem .nm{display:block;font-size:10px;color:#d8c3a2}
.varitem .meta{display:block;font-size:8.5px;color:#8e7153;text-transform:uppercase;letter-spacing:.04em;margin-top:1px}
.sphere{font-size:8.5px;letter-spacing:.14em;text-transform:uppercase;color:#947554;background:#22180f;
  border-top:1px solid var(--line);border-bottom:1px solid var(--line);padding:4px 7px;margin-top:2px}
.cmdbar{display:flex;gap:6px}
.cmdbar input{flex:1}
.console{font-size:10px;line-height:1.55;background:#0a0705;border:1px solid #3c2b1d;padding:9px;
  min-height:260px;max-height:480px;overflow:auto;white-space:pre-wrap;color:#cdb18a}
.console .ok{color:var(--green)} .console .bad{color:#e08a6a} .console .cmd{color:var(--teal)}
.quick{display:flex;gap:4px;flex-wrap:wrap;margin-top:6px}
.quick button{font-size:9px;padding:3px 6px;background:#22180f;border:1px solid var(--line2);color:#c7ab87;text-transform:uppercase}
.quick button:hover{color:#f0c078;border-color:var(--teal)}
.flist{max-height:520px;overflow:auto}
.fitem{border-bottom:1px solid #2a1e14}
.fhead{display:flex;align-items:center;gap:7px;padding:5px 7px;cursor:pointer;font-size:10px;color:#d8c3a2}
.fhead:hover{background:#22180f}
.fhead .cnt{margin-left:auto;font-size:8.5px;color:#947554}
.fbody{display:none;padding:4px 7px 7px;background:#0e0a07}
.fitem.open .fbody{display:block}
.fsig{font-size:9px;color:#b39570;padding:2px 0;border-bottom:1px dotted #2a1e14;word-break:break-word}
.gridmap{display:grid;gap:1px;background:#2a1e14;border:1px solid #4a3524;padding:1px}
.gridcell{height:13px;min-width:2px;background:#2a1e14}
.dsgrid{display:grid;grid-template-columns:repeat(2,1fr);gap:8px}
.dscard{border:1px solid var(--line);background:#1a120c;padding:10px}
.dscard h3{font-size:10.5px;margin:0 0 3px;color:#f0dfc2;text-transform:uppercase;letter-spacing:.05em}
.dscard .id{font-size:8.5px;color:var(--teal);word-break:break-all}
.dscard p{font-size:9.5px;color:#a38562;margin:5px 0}
.dscard .meta{font-size:8.5px;color:#947554;text-transform:uppercase;letter-spacing:.05em}
.telemetry{font-size:8.5px;color:#8e7153;letter-spacing:.08em;text-transform:uppercase}
.statusline{display:flex;gap:12px;flex-wrap:wrap;border:1px solid var(--line);background:#0d0906;
  padding:5px 8px;font-size:8.5px;color:#947554;letter-spacing:.08em;text-transform:uppercase;margin-bottom:8px}
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
  <div class="navlabel">Mission</div>
  <button class="nav active" data-nav="workflow"><span class="ico">TD</span>Trend detective</button>
  <button class="nav" data-nav="nasa"><span class="ico">ND</span>NASA data</button>
  <div class="navlabel">Engine</div>
  <button class="nav" data-nav="formulas"><span class="ico">FX</span>Formula registry</button>
  <button class="nav" data-nav="console"><span class="ico">&gt;_</span>Command console</button>
  <button class="nav" data-nav="catalog"><span class="ico">CA</span>Variable catalog</button>
  <div class="navlabel">Laboratories</div>
  <button class="nav" data-nav="timeseries"><span class="ico">TS</span>Time series</button>
  <button class="nav" data-nav="co2"><span class="ico">CO</span>CO₂ curve fit</button>
  <button class="nav" data-nav="retrieval"><span class="ico">RS</span>Retrievals</button>
  <button class="nav" data-nav="physics"><span class="ico">PH</span>Physics</button>
  <button class="nav" data-nav="info"><span class="ico">MI</span>Information theory</button>
  <button class="nav" data-nav="uncertainty"><span class="ico">σ</span>Uncertainty</button>
  <button class="nav" data-nav="overview"><span class="ico">SY</span>System</button>
  <div class="sidefoot"><span class="dot"></span>LOCAL C++17 ENGINE<br><span id="foot-catalog">—</span> catalog variables<br><span id="foot-formulas">—</span> public functions</div>
</aside>
<main class="main">
  <header class="topbar">
    <div class="crumb"><span>TREND DETECTIVE</span><span>/</span><b id="crumb">Workflow</b></div>
    <div class="pill"><i></i>local engine online</div>
  </header>

  <!-- TREND DETECTIVE WORKFLOW -->
  <section class="page active" id="page-workflow">
    <div class="sectionhead" style="margin-top:0"><div><h2>Trend detective</h2><p>Change is always occurring in Earth's interconnected system, and a measured variable can trend one way in one region and the opposite way in another. Select a mission or model variable, see its time evolution, then read what is changing, where, how much, and whether it is significant. Nothing is tested until you confirm the data.</p></div><div class="telemetry">NASA / EARTH OBSERVATION / LOCAL</div></div>
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

  <!-- NASA DATA -->
  <section class="page" id="page-nasa">
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
  </section>

  <!-- FORMULA REGISTRY -->
  <section class="page" id="page-formulas">
    <div class="sectionhead" style="margin-top:0"><div><h2>Formula registry</h2><p>Every public C++ function declaration in the engine, indexed from the headers.</p></div><div class="telemetry" id="formula-count">—</div></div>
    <div class="searchrow"><input id="formula-search" type="text" placeholder="Search functions, modules, symbols…"><button class="btn" id="formula-refresh">Rescan headers</button></div>
    <div class="panel"><div class="flist" id="formula-list"><div class="empty">Loading registry…</div></div></div>
    <div class="note">The registry is generated from <span class="telemetry">src/*.hpp</span> at runtime. Use the Command Console to execute any command-line implementation; the Laboratories expose the retrieval and analysis paths interactively.</div>
  </section>

  <!-- COMMAND CONSOLE -->
  <section class="page" id="page-console">
    <div class="sectionhead" style="margin-top:0"><div><h2>Command console</h2><p>Every CLI formula and diagnostic, executed by the same binary.</p></div><div class="telemetry">STDIN → C++17</div></div>
    <div class="tool">
      <div class="panel"><div class="panelhead"><b>Transmit</b><span class="hint">climate &lt;command&gt;</span></div><div class="panelbody">
        <div class="field"><label>Command</label><div class="cmdbar"><input id="cmd-input" type="text" value="catalog --format json"><button class="btn primary" id="cmd-run">Run</button></div></div>
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
        <div class="eyebrow">NASA climate methods / dependency-free C++17</div>
        <h1>Every equation, <em>open to inspection.</em></h1>
        <p class="lede">A working bench for satellite retrievals, trend detection and climate diagnostics.
          No black boxes, no numerical libraries — the formulas from the notes, compiled and traceable.</p>
        <div class="actions">
          <button class="btn primary" data-nav="workflow">Open trend detective →</button>
          <button class="btn" data-nav="formulas">Browse all formulas</button>
        </div>
      </div>
      <div class="scope"><div class="scopehead"><span>SIGNAL MONITOR</span><span>LOCAL DEMO SST</span></div><div id="hero-chart" style="flex:1;display:grid;place-items:center;color:#5a4330;font-size:9px">AWAITING DATASET</div><div class="scopefoot">CH 01 / MONTHLY / DEG C / NO SMOOTHING</div></div>
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
            <input type="file" id="ts-file" accept=".csv,text/csv" style="margin-top:7px;font-size:11px;color:#8ba0ba">
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
            <input type="file" id="co2-file" accept=".csv,text/csv" style="margin-top:7px;font-size:11px;color:#8ba0ba">
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

  <!-- PHYSICS -->
  <section class="page" id="page-physics">
    <div class="sectionhead" style="margin-top:0"><div><h2>Physics calculators</h2><p>Energy balance, radiative forcing, Mie scattering and wind-driven transport</p></div></div>
    <div class="tool">
      <div class="panel"><div class="panelhead"><b>Calculator</b><span class="hint">SI unless marked</span></div><div class="panelbody">
        <div class="field"><label>Formula</label><select id="ph-formula">
          <option value="toa">TOA energy imbalance</option>
          <option value="forcing-co2">CO₂ radiative forcing</option>
          <option value="ecs">Equilibrium climate sensitivity</option>
          <option value="mie">Mie scattering · homogeneous sphere</option>
          <option value="wind">Bulk wind stress</option>
          <option value="ekman">Ekman transport</option>
          <option value="aod">Aerosol optical depth</option>
          <option value="equilibrium">Equilibrium temperature change</option>
          <option value="coriolis">Coriolis parameter</option>
          <option value="ekman-pumping">Ekman pumping</option>
          <option value="geostrophic">Geostrophic velocity</option>
          <option value="thermal-wind">Thermal wind shear</option>
          <option value="savi">SAVI vegetation index</option>
          <option value="evi">EVI vegetation index</option>
        </select></div>
        <div id="ph-fields"></div>
        <button class="btn primary run" id="ph-run">Calculate →</button>
        <div class="formula" id="ph-formula-desc"></div>
      </div></div>
      <div id="ph-output"><div class="panel"><div class="placeholder"><div><div class="big">ƒ</div><p>Pick a formula. The fields and the equation update with the selection.</p></div></div></div></div>
    </div>
  </section>

  <!-- CATALOG -->
  <section class="page" id="page-catalog">
    <div class="sectionhead" style="margin-top:0"><div><h2>Variable catalog</h2><p>Process and region tags drive the deterministic intent compiler</p></div></div>
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
<script>
const $ = s => document.querySelector(s);
const $$ = s => Array.from(document.querySelectorAll(s));
const state = {overview:null, catalog:[], physics:null, workflow:null, nasa:null, formulas:null, help:""};

function esc(s){return String(s ?? '').replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));}
function fmt(v,d=3){if(v===null||v===undefined||!isFinite(v))return '—';
  const a=Math.abs(v); if(a!==0&&(a<1e-3||a>=1e6))return v.toExponential(2);
  return Number(v).toFixed(d).replace(/\.?0+$/,'');}
function pval(v){if(!isFinite(v))return '—'; return v<0.0001?'< 0.0001':v.toFixed(4);}
function toast(msg){const t=$('#toast');t.textContent=msg;t.classList.add('show');clearTimeout(t._h);t._h=setTimeout(()=>t.classList.remove('show'),6000);}
async function api(path,body){
  const r=await fetch(path,{method:body?'POST':'GET',headers:body?{'Content-Type':'application/json'}:{},body:body?JSON.stringify(body):undefined});
  const text=await r.text(); let data;
  try{data=JSON.parse(text);}catch(e){throw new Error(text.slice(0,300)||'Invalid server response');}
  if(!r.ok)throw new Error(data.error||('HTTP '+r.status));
  return data;
}
function nav(name){
  $$('.page').forEach(p=>p.classList.toggle('active',p.id==='page-'+name));
  $$('.nav').forEach(b=>b.classList.toggle('active',b.dataset.nav===name));
  const label={workflow:'Trend detective',nasa:'NASA data',formulas:'Formula registry',console:'Command console',overview:'System',timeseries:'Time series',co2:'CO₂ curve fit',retrieval:'Retrievals',info:'Information theory',uncertainty:'Uncertainty',physics:'Physics',catalog:'Catalog'}[name]||name;
  $('#crumb').textContent=label; window.scrollTo({top:0,behavior:'smooth'});
}
document.addEventListener('click',e=>{const b=e.target.closest('[data-nav]');if(b){nav(b.dataset.nav);}});
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
    g+=`<path class="line" d="${d}" stroke="${s.color||'#f0d9b5'}"/>`;
  });
  let legend='';let lx=L;
  series.forEach(s=>{legend+=`<line x1="${lx}" y1="8" x2="${lx+10}" y2="8" stroke="${s.color||'#f0d9b5'}" stroke-width="2"/><text class="legend" x="${lx+14}" y="11">${esc(s.name)}</text>`;lx+=s.name.length*6+34;});
  return `<svg class="chart" viewBox="0 0 ${W} ${H}" role="img"><g>${g}</g><g>${legend}</g></svg>`;
}
function chartCard(title,series,opts={}){return `<div class="chartwrap"><div style="font-size:9.5px;color:#a38562;padding:2px 8px 0;text-transform:uppercase;letter-spacing:.08em">${esc(title)}</div>${chart(series,opts)}</div>`;}

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
    ${chartCard('Observed series and STL trend',[{name:'observed',values:d.series.observed,color:'#f0d9b5'},{name:'trend',values:d.components.trend,color:'#ff9a4d'}],{xlabels:labels,height:300})}
    <div class="minicharts">
      ${chartCard('Seasonal component',[{name:'seasonal',values:d.components.seasonal,color:'#5fc2c2'}],{xlabels:labels,height:210})}
      ${chartCard('Remainder',[{name:'remainder',values:d.components.remainder,color:'#d98a8a'}],{xlabels:labels,height:210})}
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
    ${chartCard('Observed and fitted curve',[{name:'observed',values:d.series.observed,color:'#f0d9b5'},{name:'fitted',values:d.series.fitted,color:'#ff9a4d'}],{xlabels:d.series.time.map(t=>fmt(t,1)),height:310})}
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
      <div class="bigvalue">${fmt(d.ksg,4)} <span style="font-size:13px;color:#8ba0ba">bits · KSG</span></div>
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
function renderPhysics(){
  const def=physics[$('#ph-formula').value];state.physics=def;
  $('#ph-formula-desc').textContent=def.desc;
  $('#ph-fields').innerHTML=def.fields.map(f=>`<div class="field"><label>${esc(f[1])}</label><input data-ph="${f[0]}" type="number" value="${f[2]}" step="any"></div>`).join('');
}
$('#ph-formula').addEventListener('change',renderPhysics);
$('#ph-run').addEventListener('click',async()=>{
  try{
    const body={formula:$('#ph-formula').value};
    $$('[data-ph]').forEach(i=>body[i.dataset.ph]=+i.value);
    const d=await api('/api/physics',body);
    $('#ph-output').innerHTML=`<div class="panel"><div class="panelhead"><b>${esc(d.title)}</b><span class="hint">${esc(d.formula)}</span></div><div class="panelbody">
      <div class="bigvalue">${fmt(d.primary.value,5)}</div><div class="unit">${esc(d.primary.unit)}</div>
      <div class="kv">${d.details.map(x=>`<span class="k">${esc(x.label)}</span><span class="v">${fmt(x.value,5)} ${esc(x.unit||'')}</span>`).join('')}</div>
      <div class="note">${esc(d.note||'')}</div>
    </div></div>`;
  }catch(err){toast(err.message);}
});

function renderCatalog(){
  const q=$('#cat-search').value.toLowerCase(),p=$('#cat-process').value,r=$('#cat-region').value;
  const list=state.catalog.filter(e=>(!p||e.processes.includes(p))&&(!r||e.regions.includes(r))&&
    (!q||(e.id+' '+e.name+' '+e.description+' '+e.retrieval).toLowerCase().includes(q)));
  if(!list.length){$('#cat-list').innerHTML='<div class="empty">No catalog variables match those filters.</div>';return;}
  $('#cat-list').innerHTML=`<table class="table"><thead><tr><th>Variable</th><th>Processes</th><th>Regions</th><th>Retrieval / analysis</th></tr></thead><tbody>
    ${list.map(e=>`<tr><td class="id">${esc(e.id)}<br><span style="font-family:var(--sans);color:#a9b9cf">${esc(e.name)}</span><br><span style="color:#6f829d">${esc(e.units)}</span></td>
    <td>${e.processes.map(x=>`<span class="tag">${esc(x)}</span>`).join(' ')}</td>
    <td class="mut">${e.regions.map(esc).join(', ')}</td>
    <td class="mut">${esc(e.retrieval)}<br><span style="color:#6f829d">${esc(e.analysis)}</span></td></tr>`).join('')}
  </tbody></table>`;
}
function loadCatalog(){
  api('/api/catalog').then(d=>{state.catalog=d;renderCatalog();renderWfVars();})
    .catch(e=>{$('#cat-list').innerHTML=`<div class="empty">${esc(e.message)}</div>`;});
}
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
  air_temperature:{file:'air_temperature_monthly.csv',value:'value',sphere:'ATMOSPHERE',product:'AIRS / MERRA-2 surface air temperature',pid:'AIRS L3 / M2T1NXSLV',provider:'GES DISC',res:'1 degree monthly',grid:false,
    artifacts:['Surface air temperature is not skin temperature; retrieval and reanalysis differ.','Station/reanalysis sampling is spatially uneven.','Demo series is synthetic; it is not an AIRS or MERRA-2 download.'],
    retrieval:'Gridded temperature retrieval / reanalysis variable',analysis:'STL decomposition + Mann-Kendall + Sen slope'},
  co2:{file:'co2_monthly.csv',value:'co2',sphere:'ATMOSPHERE',product:'OCO-2 Lite XCO2',pid:'OCO2_L2_Lite_FP',provider:'GES DISC',res:'point retrievals, monthly composite',grid:false,
    artifacts:['XCO2 is a column average, not a surface concentration.','OCO-2 has a narrow swath and sparse sampling.','This demo series is synthetic; use the CO2 curve-fit lab for the polyharmonic model.'],
    retrieval:'OCO-2 XCO2 retrieval',analysis:'Polynomial + harmonics (CO2 lab); STL + MK in this workflow'},
  ndvi:{file:'ndvi_monthly.csv',value:'value',sphere:'LAND / BIOSPHERE',product:'MODIS/Terra vegetation indices',pid:'MOD13A3',provider:'LP DAAC',res:'1 km monthly',grid:false,
    artifacts:['NDVI saturates in dense canopy.','Cloud and aerosol contamination persist in composites.','Land-cover change can mimic a trend.'],
    retrieval:'NDVI = (NIR - RED)/(NIR + RED)',analysis:'STL decomposition + Mann-Kendall + Sen slope'},
  precipitation:{file:'precipitation_monthly.csv',value:'value',sphere:'ATMOSPHERE',product:'GPCP monthly precipitation',pid:'GPCPMON',provider:'GES DISC',res:'0.5 degree monthly',grid:false,
    artifacts:['Satellite-gauge merging has regionally varying bias.','Monthly aggregation hides extreme daily events.','Demo series is synthetic.'],
    retrieval:'GPI cold-cloud index + GPCP inverse-variance merge',analysis:'STL decomposition + Mann-Kendall + Sen slope'},
  sea_ice_extent:{file:'sea_ice_extent_monthly.csv',value:'value',sphere:'CRYOSPHERE',product:'NSIDC sea ice concentration CDR',pid:'NSIDC-0051',provider:'NSIDC DAAC',res:'25 km / monthly extent',grid:false,
    artifacts:['Extent uses a 15% concentration threshold; it is not area.','Passive-microwave retrievals are weather-dependent.','Demo series is synthetic.'],
    retrieval:'NASA Team concentration + 15% extent rule',analysis:'STL decomposition + Mann-Kendall + Sen slope'},
  sea_ice_concentration:{file:'sea_ice_extent_monthly.csv',value:'value',sphere:'CRYOSPHERE',product:'NSIDC sea ice concentration CDR',pid:'NSIDC-0051',provider:'NSIDC DAAC',res:'25 km monthly',grid:false,
    artifacts:['Concentration is a fractional coverage, not thickness.','Melt ponds and snow cover affect the brightness temperatures.','The workflow file is an extent series, not a concentration field.'],
    retrieval:'NASA Team / Bootstrap concentration',analysis:'STL decomposition + Mann-Kendall + Sen slope'},
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
$('#wf-vars')&&$('#wf-vars').addEventListener('click',e=>{const item=e.target.closest('.varitem');if(item)selectWfVar(item.dataset.id);});
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
  let vi=meta.value?table.header.findIndex(h=>h.toLowerCase()===meta.value.toLowerCase()):-1;
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
  const file=meta.file||(entry.id+'_monthly.csv');
  if(!state.overview||!state.overview.files.includes(file)){
    wfLog(`No bundled demo file for ${entry.id}; requesting upload.`,'warn');
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
  return chartCard('Observed and STL trend',[{name:'observed',values:d.series.observed,color:'#f0d9b5'},{name:'trend',values:d.components.trend,color:'#ff9a4d'}],{xlabels:labels,height:250})+
    `<div class="minicharts">${chartCard('Seasonal',[{name:'seasonal',values:d.components.seasonal,color:'#5fc2c2'}],{xlabels:labels,height:210})}${chartCard('Remainder',[{name:'remainder',values:d.components.remainder,color:'#d98a8a'}],{xlabels:labels,height:210})}</div>`;
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
function renderWfReport(d){
  const e=state.workflow.entry,meta=state.workflow.meta,mk=d.trend.mann_kendall,sen=d.trend.sen;
  $('#wf-report').innerHTML=`<div class="statusline"><span>REPORT <b>GENERATED</b></span><span>VARIABLE <b>${esc(e.id)}</b></span><span>REGION <b>${esc($('#wf-region').value)}</b></span></div>
    <div class="verdict">${mk.significant?'The series shows a statistically significant '+esc(mk.direction)+' trend.':'No statistically significant monotonic trend was detected.'}</div>
    <div class="facts"><div class="fact"><div class="k">What</div><div class="v">${esc(e.name)}</div><div class="s">${esc(e.units||'')}</div></div><div class="fact"><div class="k">Where</div><div class="v">${esc($('#wf-region').value)}</div><div class="s">catalog region tag</div></div><div class="fact"><div class="k">How much</div><div class="v">${fmt(d.trend.sen_per_decade,3)} / decade</div><div class="s">Sen slope</div></div><div class="fact"><div class="k">Significance</div><div class="v">p ${pval(mk.p_value)}</div><div class="s">Mann-Kendall</div></div><div class="fact"><div class="k">Coverage</div><div class="v">${state.workflow.coverage.toFixed(1)}%</div><div class="s">requested months</div></div><div class="fact"><div class="k">Gaps</div><div class="v">${state.workflow.runs.length}</div><div class="s">runs of missing data</div></div></div>
    <div class="note">Formulas applied</div><div class="formula">${esc(meta.retrieval)}<br>STL: Y = T + S + R (Cleveland et al.)<br>Mann-Kendall S with tie correction + Hamed-Rao variance<br>Sen slope + rank-based 95% CI</div>
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
    wfLog('Data gate confirmed. Running STL + trend test.');
    const d=await api('/api/analyze',{csv:state.workflow.csv,time:state.workflow.timeName,value:state.workflow.valueName,period:12,alpha:0.05,units:e.units||'',hamed_rao:true});
    state.analysis=d;
    $('#wf-method').textContent=d.decomposed?'STL':'CLIMATOLOGY';
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
  }catch(e){out.textContent='ERROR: '+e.message;}
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
function renderFormulas(){
  const q=($('#formula-search').value||'').toLowerCase();
  if(!state.formulas)return;
  const mods=state.formulas.modules.map(m=>({...m,functions:m.functions.filter(f=>!q||f.toLowerCase().includes(q))})).filter(m=>m.functions.length);
  $('#formula-list').innerHTML=mods.map(m=>`<div class="fitem"><div class="fhead" data-fopen><span>${esc(m.file)}</span><span class="cnt">${m.functions.length}</span></div><div class="fbody">${m.functions.map(f=>`<div class="fsig">${esc(f)}</div>`).join('')}</div></div>`).join('')||'<div class="empty">No function signature matches that filter.</div>';
}
$('#formula-search')&&$('#formula-search').addEventListener('input',renderFormulas);
$('#formula-refresh')&&$('#formula-refresh').addEventListener('click',loadFormulas);
$('#formula-list')&&$('#formula-list').addEventListener('click',e=>{const h=e.target.closest('[data-fopen]');if(h)h.parentElement.classList.toggle('open');});

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
$('#cmd-run')&&$('#cmd-run').addEventListener('click',runConsole);
$('#cmd-input')&&$('#cmd-input').addEventListener('keydown',e=>{if(e.key==='Enter')runConsole();});
$$('.quick button').forEach(b=>b.addEventListener('click',()=>{$('#cmd-input').value=b.dataset.cmd;runConsole();}));

/* ---------------- BOOT ---------------- */
renderPhysics();
loadOverview();
loadNasa();
loadFormulas();
loadHelp();
</script>
</body>
</html>)CLIMATE_UI";
}

}  // namespace climate
