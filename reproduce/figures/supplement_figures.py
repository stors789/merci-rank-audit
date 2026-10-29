"""Redraw the four accepted supplementary figures from portable aggregate tables."""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

p=argparse.ArgumentParser();p.add_argument("--data",type=Path,required=True);p.add_argument("--output",type=Path,required=True)
args=p.parse_args();D=args.data;O=args.output;O.mkdir(parents=True,exist_ok=True)
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":9.5,"axes.labelsize":9.5,
                     "xtick.labelsize":9,"ytick.labelsize":9,"pdf.fonttype":42,
                     "axes.spines.top":False,"axes.spines.right":False})
def save(fig,name):
    fig.savefig(O/(name+".pdf"));fig.savefig(O/(name+".png"),dpi=170);plt.close(fig)

# S1: seven real inputs; weighted and accepted share gate estimate.
d=pd.read_csv(D/"FigureS1_gate.csv")
cases=[("bench1_split21_mixed","B1 split 21 mixed"),("bench1_split21_negative","B1 split 21 negative"),
       ("bench2_split1_mixed","B2 split 1 mixed"),("bench2_split1_negative","B2 split 1 negative"),
       ("P2","P2"),("MC_only","MC only"),("GBML181","GBML181")]
fig,axes=plt.subplots(4,2,figsize=(7.0,9.0),sharex=True,sharey=True)
fig.subplots_adjust(left=.16,right=.98,top=.86,bottom=.13,wspace=.26,hspace=.78)
for ax,(case,label) in zip(axes.flat,cases):
    z=d[d.case==case];ax.set_xscale("log");ax.set_yscale("log");ax.set_xlim(.08,12);ax.set_ylim(1e-8,1)
    ax.set_xticks([.1,.5,2,10],["0.1","0.5","2","10"]);ax.grid(alpha=.22)
    ax.set_title(f"{label}\nP(G) ≈ {float(z.reference_P_G.iloc[0]):.4g}",fontsize=9)
    allzero=True
    for method,color,style in [("uniform","#2878a8","-"),("conditional_weighted","#d16b36","--")]:
        q=z[z.method==method].sort_values("target_seconds");x=q.target_seconds.to_numpy(float)
        y=q.gate_median.to_numpy(float);lo=q.gate_q25.to_numpy(float);hi=q.gate_q75.to_numpy(float)
        allzero &= np.isfinite(y).all() and np.all(y==0)
        y=np.where(y>0,y,np.nan);lo=np.where(lo>0,lo,np.nan);hi=np.where(hi>0,hi,np.nan)
        ax.plot(x,y,style,color=color,marker="o",ms=3.5,lw=1.5);ax.fill_between(x,lo,hi,color=color,alpha=.09)
    half=float(z.reference_gate_95_halfwidth.iloc[0])
    if half>0:ax.axhline(half,color="#666666",lw=1,ls="-.")
    if allzero:ax.text(.5,.5,"Gate error = 0\n20/20 at all checkpoints",transform=ax.transAxes,ha="center",va="center",fontsize=9)
    ax.tick_params(labelbottom=True)
axes.flat[-1].axis("off")
fig.legend(handles=[Line2D([],[],color="#2878a8",marker="o",label="Uniform"),
    Line2D([],[],color="#d16b36",ls="--",marker="o",label="Weighted"),
    Line2D([],[],color="#666666",ls="-.",label="Reference gate 95% half-width")],
    loc="upper center",bbox_to_anchor=(.5,.945),ncol=3,fontsize=9)
fig.suptitle("Equal-time gate error in seven fixed inputs",fontsize=11,y=.98)
fig.supylabel("Absolute gate-probability error",x=.025,fontsize=9.5)
fig.supxlabel("Target sampling time (s)",y=.08,fontsize=9.5)
fig.text(.5,.025,"20 technical repetitions; IQR bands; two-worker batch.",ha="center",fontsize=9)
save(fig,"FigureS1_gate")

# S2: fixed-reference states across three stages and both designs.
d=pd.read_csv(D/"FigureS2_cohorts.csv")
stages=["development_reused","frozen_evaluation","retrospective_second_experiment"]
labels=["B1 development 20","B1 frozen 80","B2 retrospective 20"];levels=[0,25,50,75,100]
fig,axes=plt.subplots(3,2,figsize=(7.0,8.0),sharex=True,sharey=True)
fig.subplots_adjust(left=.09,right=.98,top=.86,bottom=.13,hspace=.40,wspace=.16)
lines=[("raw","core_fraction","Raw core","#2878a8","-"),
       ("conditional","core_fraction","Passing-only core","#d16b36","--"),
       ("raw","projected_default_call_fraction","Parent-order projected default call","#5c8c50",":")]
for i,(stage,label) in enumerate(zip(stages,labels)):
 for j,design in enumerate(["nested","fixedN"]):
    ax=axes[i,j]
    for family,metric,lab,color,style in lines:
        q=d[(d.stage==stage)&(d.design==design)&(d.family==family)&(d.metric==metric)].set_index("level").reindex(levels)
        y=q["median"].to_numpy(float);lo=q.q25.to_numpy(float);hi=q.q75.to_numpy(float)
        ax.plot(levels,y,style,color=color,marker="o",ms=3.5,lw=1.5,label=lab);ax.fill_between(levels,lo,hi,color=color,alpha=.1)
    ax.set_title(label+" · "+("CC addition" if design=="nested" else "Fixed N replacement"),fontsize=10)
    ax.set_ylim(0,1);ax.set_xticks(levels);ax.grid(alpha=.22)
    if j==0:ax.set_ylabel("Fraction of fixed MC reference")
    if i==2:ax.set_xlabel("CC added (% available pool)" if design=="nested" else "CC replacement (% background slots)")
fig.legend(*axes[0,0].get_legend_handles_labels(),loc="upper center",bbox_to_anchor=(.5,.94),ncol=3,fontsize=9)
fig.suptitle("Candidate composition and fixed-reference label states",fontsize=11,y=.99)
fig.text(.5,.035,"Median and technical-split IQR; 100% fixed-N replacement is 50% CC of the full cohort.",ha="center",fontsize=8.5)
save(fig,"FigureS2_cohorts")

# S3: all four sequences and all prefix defined denominators.
summary=pd.read_csv(D/"FigureS3_precision.csv")
reference=pd.read_csv(D/"FigureS3_reference_quantiles.csv")
series=[("original","HALLMARK_MTORC1_SIGNALING"),("original","HALLMARK_PROTEIN_SECRETION"),
        ("alternative","HALLMARK_MTORC1_SIGNALING"),("alternative","HALLMARK_PROTEIN_SECRETION")]
colors={.025:"#2878a8",.5:"#d16b36",.975:"#648a39"};labels={.025:"2.5% quantile",.5:"Median",.975:"97.5% quantile"}
fig,axes=plt.subplots(2,2,figsize=(7.0,6.4));fig.subplots_adjust(left=.15,right=.98,top=.84,bottom=.18,hspace=.62,wspace=.27)
for ax,(rep,feat) in zip(axes.flat,series):
    part=summary[(summary.representation==rep)&(summary.feature==feat)]
    ref=reference[(reference.representation==rep)&(reference.feature==feat)]
    for q in [.025,.5,.975]:
        z=part[part["quantile"]==q].set_index("budget").reindex([100,1000,5000,20000])
        y=100*z.difference_over_old_width_median.to_numpy(float)
        lo=100*z.difference_over_old_width_q25.to_numpy(float);hi=100*z.difference_over_old_width_q75.to_numpy(float)
        ax.plot(range(4),y,"o-",color=colors[q],label=labels[q],lw=1.5,ms=4)
        ax.fill_between(range(4),lo,hi,color=colors[q],alpha=.1)
        rr=ref[ref["quantile"]==q]
        ax.axhline(100*rr.iloc[0].reference_max_half_width/rr.iloc[0].old_interval_width,color=colors[q],ls=":",lw=.8)
    ax.axhline(2.5,color="#666666",ls="--",lw=.8)
    ax.set_title(("mTORC1" if "MTORC1" in feat else "Protein secretion")+" · "+rep)
    rows=part[part["quantile"]==.5].set_index("budget").reindex([100,1000,5000,20000])
    ax.set_xticks(range(4),[f"{b}\n{int(r.comparison_defined)}/{int(r.planned)}" for b,r in rows.iterrows()])
    ax.set_xlim(-.15,3.15);ax.grid(alpha=.22)
fig.supylabel("Absolute difference / old interval width (%)",x=.025,fontsize=9.5)
fig.supxlabel("Proposal budget (defined / planned)",y=.08,fontsize=9.5)
fig.legend(*axes[0,0].get_legend_handles_labels(),loc="upper center",bbox_to_anchor=(.5,.94),ncol=3,fontsize=9)
fig.suptitle("Numerical precision across independent conditional streams",fontsize=11,y=.99)
fig.text(.5,.025,"20 streams; IQR bands; dotted lines: reference MC error.",ha="center",fontsize=8.5)
save(fig,"FigureS3_precision")

# S4: one-seed empirical union shortfall.
d=pd.read_csv(D/"FigureS4_union.csv")
cases=[("all_tied_certain","N8 · certain"),("all_tied_single","N8 · single gate"),
       ("all_tied_overlapping","N8 · overlapping gates"),("P2","P2 · rare multi-gate")]
fig,axes=plt.subplots(2,2,figsize=(7.0,5.7),sharex=True,sharey=True)
fig.subplots_adjust(left=.15,right=.97,top=.84,bottom=.18,hspace=.35,wspace=.2)
for ax,(case,label) in zip(axes.flat,cases):
 for method,color,style,marker in [("uniform","#2878a8","-","o"),("conditional","#d16b36","--","s")]:
    q=d[(d.case==case)&(d.method==method)].set_index("budget").reindex([100,1000,5000,20000])
    y=100*q.union_unseen.to_numpy(float)/q.exact_passing_union.to_numpy(float)
    ax.plot(range(4),y,style,marker=marker,ms=5,lw=1.5,color=color,
            markerfacecolor="white" if method=="conditional" else color,label=method.capitalize())
 ax.set_title(label);ax.set_xticks(range(4),["100","1,000","5,000","20,000"]);ax.set_ylim(-4,104)
 ax.grid(alpha=.22)
fig.supylabel("Exact union not yet seen (%)",x=.025,fontsize=9.5)
fig.supxlabel("Nested proposal budget",y=.10,fontsize=9.5)
fig.legend(*axes[0,0].get_legend_handles_labels(),loc="upper center",bbox_to_anchor=(.5,.94),ncol=2,fontsize=9)
fig.suptitle("Members not yet observed versus the exact passing-only union",fontsize=11,y=.99)
fig.text(.5,.035,"P2 accepted: uniform 0, 1, 3, 9; conditional 100, 981, 4,911, 19,661.\nZero acceptance gives undefined conditional frequencies.",ha="center",fontsize=8.5)
save(fig,"FigureS4_union")
print("Redrew all four supplementary figures from supplied aggregate rows.")
