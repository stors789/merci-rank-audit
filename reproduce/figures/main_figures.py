"""Redraw column-width figures from accepted, locally selected aggregate rows."""
from pathlib import Path
import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from matplotlib.lines import Line2D

parser=argparse.ArgumentParser()
parser.add_argument('--data',type=Path,required=True)
parser.add_argument('--output',type=Path,required=True)
args=parser.parse_args()
DATA=args.data;OUT=args.output
OUT.mkdir(parents=True,exist_ok=True)
BLUE="#235e7a";ORANGE="#c36a38";TEAL="#258678";PURPLE="#715d96";GREY="#66747c";INK="#263744"
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":8.3,"axes.labelsize":8.3,
    "axes.titlesize":8.7,"xtick.labelsize":8.0,"ytick.labelsize":8.0,"legend.fontsize":8.0,
    "pdf.fonttype":42,"axes.spines.top":False,"axes.spines.right":False})
def read(name):return pd.read_csv(DATA/name,float_precision="round_trip")
def save(fig,name):
    fig.savefig(OUT/(name+".pdf"))
    fig.savefig(OUT/(name+".png"),dpi=300)
    plt.close(fig)

# A conceptual diagram, with analytic marginals and sampled orders separated.
fig,ax=plt.subplots(figsize=(3.45,2.95));ax.set_xlim(0,10);ax.set_ylim(0,10);ax.axis("off")
def box(x,y,w,h,lines,color="#eaf2f4"):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle="round,pad=0.08,rounding_size=0.15",
                 linewidth=.7,edgecolor=BLUE,facecolor=color))
    ax.text(x+w/2,y+h/2,lines,ha="center",va="center",fontsize=8.0,
            color=INK,linespacing=1.12)
def arrow(x1,y1,x2,y2):
    ax.add_patch(FancyArrowPatch((x1,y1),(x2,y2),arrowstyle="-|>",mutation_scale=7,
                                linewidth=.8,color=GREY))
box(.2,8.7,9.6,.95,"Fixed Neff blocks, RNA ranks and gate maxima")
box(.2,7.15,9.6,.95,"Legal tie orders and nested gate")
arrow(5,8.65,5,8.15)
box(.2,5.25,2.95,1.10,"Raw\nall orders")
box(3.52,5.25,2.95,1.10,"Default\nfailures empty","#fdf0e6")
box(6.85,5.25,2.95,1.10,"Passing-only\nconditional on G","#e8f4ef")
for x in (1.7,5,8.3):arrow(x,7.10,x,6.44)
box(.2,3.43,9.6,.96,"Exact core / union and full-order witnesses")
arrow(1.7,5.20,2.4,4.50);arrow(5,5.20,5,4.50);arrow(8.3,5.20,7.6,4.50)
box(.2,.77,4.45,1.55,"Analytic marginals\n(certain / single gate)","#eaf2f4")
box(5.35,.77,4.45,1.55,"Sampled full orders →\nnuclear outcomes","#e8f4ef")
arrow(2.4,3.39,2.4,2.38);arrow(7.6,3.39,7.6,2.38)
fig.subplots_adjust(left=.015,right=.985,bottom=.02,top=.98)
save(fig,"Figure1_workflow")

# Three equal-time examples, one panel per case, with common log axes.
d=read("figure2_equal_time.csv")
cases=[("P2","P2 · rare"),("GBML181","GBML181 · rare"),("bench2_split1_negative","Batch 2 negative · common")]
fig,axes=plt.subplots(3,1,figsize=(3.45,3.65),sharex=True,sharey=True)
styles=[("uniform",BLUE,"-","o"),("conditional_weighted",ORANGE,"--","s"),
        ("conditional_accepted",TEAL,":","^")]
for ax,(case,title) in zip(axes,cases):
    sub=d[d.case==case]
    for method,color,ls,marker in styles:
        z=sub[sub.method==method].sort_values("target_seconds")
        x=z.target_seconds.to_numpy(float);med=z.member_median.to_numpy(float)
        ax.plot(x,med,ls=ls,marker=marker,ms=3.2,lw=1.1,color=color)
        ax.fill_between(x,z.member_q25.to_numpy(float),z.member_q75.to_numpy(float),color=color,alpha=.11)
        if method=="uniform":
            for _,r in z[z.member_evaluable<20].iterrows():
                ax.annotate(f"{int(r.member_evaluable)}/20",(r.target_seconds,r.member_median),
                            xytext=(2,3),textcoords="offset points",fontsize=8,color=BLUE)
    ax.axhline(float(sub.reference_member_rms_MC_SE.iloc[0]),color=GREY,ls=":",lw=.9)
    ax.set_xscale("log");ax.set_yscale("log");ax.set_ylim(.0004,1)
    ax.set_yticks([.001,.01,.1,1])
    ax.set_title(title,loc="right",pad=1);ax.set_ylabel("RMSE");ax.grid(axis="y",alpha=.16)
axes[-1].set_xticks([.1,.5,2,10],labels=["0.1","0.5","2","10"])
axes[-1].set_xlabel("Sampling time (seconds)")
legend=[Line2D([],[],color=c,ls=ls,marker=m,ms=3.5,label=l) for (_,c,ls,m),l in zip(styles,["Uniform","Weighted","Accepted"])]
fig.legend(handles=legend,loc="lower center",ncol=3,frameon=False,bbox_to_anchor=(.53,-.005),columnspacing=.7,handlelength=1.4)
fig.subplots_adjust(left=.19,right=.97,top=.95,bottom=.19,hspace=.35)
save(fig,"Figure2_sampling")

# Two annotation contrasts with technical-split IQRs.
err=read("figure3_error_partitions.csv");pairs=read("figure3_paired_changes.csv")
stages=["frozen_evaluation","retrospective_second_experiment"]
fig=plt.figure(figsize=(3.45,3.4))
a=fig.add_axes([.25,.66,.72,.23]);b=fig.add_axes([.25,.18,.72,.23])
for i,stage in enumerate(stages):
    for off,metric,color,marker in [(-.12,"FP_core_fraction",BLUE,"o"),(.12,"FN_outside_fraction",ORANGE,"s")]:
        r=err[(err.stage==stage)&(err.metric==metric)].iloc[0];v=100*r["median"]
        a.errorbar(v,i+off,xerr=[[v-100*r.q25],[100*r.q75-v]],fmt=marker,ms=4,capsize=2,color=color,lw=1)
    for off,metric,color,marker in [(-.12,"precision_delta",TEAL,"o"),(.12,"sensitivity_delta",PURPLE,"s")]:
        r=pairs[(pairs.stage==stage)&(pairs.metric==metric)].iloc[0];v=100*r["median"]
        b.errorbar(v,i+off,xerr=[[v-100*r.q25],[100*r.q75-v]],fmt=marker,ms=4,capsize=2,color=color,lw=1)
for ax in (a,b):
    ax.set_yticks([0,1],["B1 frozen","B2 retro"]);ax.invert_yaxis();ax.grid(axis="x",alpha=.16)
a.set_xlim(70,101);a.set_xlabel("Fraction of original errors (%)")
fig.text(.25,.94,"● FP in core",color=BLUE,fontsize=8.0)
fig.text(.58,.94,"■ FN outside",color=ORANGE,fontsize=8.0)
b.set_xlim(-4.4,3.0);b.axvline(0,color=GREY,lw=.7);b.set_xlabel("Paired change (percentage points)")
fig.text(.25,.46,"● Precision",color=TEAL,fontsize=8.0)
fig.text(.58,.46,"■ Sensitivity",color=PURPLE,fontsize=8.0)
save(fig,"Figure3_annotation")

# Compact fixed-N curves and endpoint comparison (reference denominator in caption).
cond=read("figure4_core_curves.csv");tr=read("figure4_endpoint_transitions.csv")
fig=plt.figure(figsize=(3.45,3.35))
a=fig.add_axes([.25,.64,.72,.30]);b=fig.add_axes([.25,.19,.72,.22])
for stage,color,label in [("frozen_evaluation",BLUE,"B1 frozen"),("retrospective_second_experiment",ORANGE,"B2 retrospective")]:
    z=cond[cond.stage==stage].sort_values("level")
    a.plot(z.level,100*z["median"],marker="o",ms=3.2,lw=1.3,color=color,label=label)
    a.fill_between(z.level,100*z.q25.to_numpy(float),100*z.q75.to_numpy(float),color=color,alpha=.13)
a.set_ylim(0,30);a.set_xticks([0,25,50,75,100]);a.set_ylabel("Anchor in core (%)")
a.set_xlabel("CC replacement of slots (%)");a.grid(axis="y",alpha=.16);a.legend(frameon=False,loc="upper right",handlelength=1)
for i,stage in enumerate(stages):
    for off,design,color,marker in [(-.12,"fixedN",TEAL,"o"),(.12,"nested",PURPLE,"s")]:
        r=tr[(tr.stage==stage)&(tr.design==design)].iloc[0];v=100*r["median"]
        b.errorbar(v,i+off,xerr=[[v-100*r.q25],[100*r.q75-v]],fmt=marker,ms=4,capsize=2,color=color,lw=1)
b.set_yticks([0,1],["B1 frozen","B2 retro"]);b.invert_yaxis();b.set_xlim(0,20.5)
b.set_xlabel("Core → outside (% of reference)");b.grid(axis="x",alpha=.16)
fig.text(.29,.45,"● Fixed N",color=TEAL,fontsize=8.0)
fig.text(.59,.45,"■ Nested",color=PURPLE,fontsize=8.0)
save(fig,"Figure4_composition")

# Raw-label points, conditional quantiles, outcome-selected witness and MC/sign counts.
ref=read("figure5_reference_quantiles.csv")
sign=read("figure5_reference_signs.csv")
witness=read("figure5_selected_witness.csv")
fixed=pd.concat([read("figure5_fixed_original.csv"),read("figure5_fixed_alternative.csv")])
features=["HALLMARK_MTORC1_SIGNALING","HALLMARK_PROTEIN_SECRETION"]
fig,axes=plt.subplots(2,1,figsize=(3.45,4.2))
for ax,rep in zip(axes,["original","alternative"]):
    for i,feat in enumerate(features):
        z=ref[(ref.representation==rep)&(ref.feature==feat)].sort_values("quantile")
        q=z.reference_estimate.to_numpy(float)
        raw=fixed[(fixed.representation==rep)&(fixed.nuclear_feature==feat)].beta.iloc[0]
        selected=witness[(witness.representation==rep)&(witness.nuclear_feature==feat)].beta.iloc[0]
        ax.plot([q[0],q[2]],[i,i],color=BLUE,lw=3,solid_capstyle="round")
        ax.plot(q[1],i,"o",color=ORANGE,ms=4)
        ax.plot(raw,i,"D",color=INK,ms=4)
        ax.plot(selected,i,"x",color=PURPLE,ms=5,mew=1.2)
    ax.axvline(0,color=GREY,lw=.8,ls=":")
    ax.set_yticks([0,1],["mTORC1","Secretion"]);ax.invert_yaxis()
    ax.set_title(rep.capitalize()+" representation",loc="right",pad=1);ax.grid(axis="x",alpha=.15)
axes[0].set_xlim(-.012,.035);axes[0].set_xticks([-.01,0,.01,.02,.03])
axes[1].set_xlim(-.0022,.0092);axes[1].set_xticks([-.002,0,.002,.004,.006,.008])
axes[1].set_xlabel("β coefficient")
handles=[Line2D([],[],color=BLUE,lw=3,label="2.5–97.5%"),
         Line2D([],[],color=ORANGE,marker="o",ls="",ms=4,label="Median"),
         Line2D([],[],color=INK,marker="D",ls="",ms=4,label="Fixed raw"),
         Line2D([],[],color=PURPLE,marker="x",ls="",ms=5,label="Selected witness")]
fig.legend(handles=handles,loc="upper center",ncol=2,frameon=False,bbox_to_anchor=(.52,.35),
           columnspacing=.8,handlelength=1.4)
fig.text(.05,.155,"MC max (%) / negative counts in 196,549 orders",fontsize=8.0)
for row,rep,label in [(.105,"original","Orig"),(.055,"alternative","Alt")]:
    cells=[]
    for feat,short in zip(features,["mTOR","secretion"]):
        z=ref[(ref.representation==rep)&(ref.feature==feat)]
        max_mc=100*(z.reference_max_half_width/z.old_interval_width).max()
        negatives=int(sign[(sign.representation==rep)&(sign.feature==feat)].negative.iloc[0])
        cells.append(f"{short} {max_mc:.2f} / {negatives}")
    fig.text(.05,row,label+": "+"; ".join(cells),fontsize=8.0)
fig.subplots_adjust(left=.25,right=.94,top=.96,bottom=.46,hspace=.65)
save(fig,"Figure5_nuclear")

print('Redrew five main figures from supplied aggregate data.')
