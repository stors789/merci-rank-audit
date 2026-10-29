"""Redraw BIO-06 gate-probability sensitivity from its aggregate table."""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

p=argparse.ArgumentParser();p.add_argument("--data",type=Path,required=True);p.add_argument("--output",type=Path,required=True);a=p.parse_args()
data=pd.read_csv(a.data/"rarity_comparison.csv");a.output.mkdir(parents=True,exist_ok=True)
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":9.5,
                     "axes.spines.top":False,"axes.spines.right":False,"pdf.fonttype":42})
fig,axes=plt.subplots(1,2,figsize=(7.0,3.9),sharey=True)
for ax,case in zip(axes,["P2","GBML181"]):
 z=data[data.case==case].set_index("null_id").loc[["original"]+[f"null_{i:02d}" for i in range(1,11)]]
 x=np.arange(11);y=z.P_G.to_numpy(float);err=z.MC_95_halfwidth.to_numpy(float)
 ax.axhspan(y[0]-err[0],y[0]+err[0],color="#235e7a",alpha=.12);ax.axhline(y[0],color="#235e7a",lw=1,ls="--")
 ax.errorbar(x[1:],y[1:],yerr=err[1:],fmt="o",color="#c36a38",capsize=2,lw=1,label="New null seeds")
 ax.errorbar([0],[y[0]],yerr=[err[0]],fmt="s",color="#235e7a",capsize=2,lw=1,label="Original threshold")
 ax.set_yscale("log");ax.set_ylim(.0002,.002);ax.set_xticks(x,["Old"]+[str(i) for i in range(1,11)])
 ax.set_title(case);ax.set_xlabel("Fixed null-threshold seed");ax.grid(axis="y",alpha=.2)
axes[0].set_ylabel("P(gate passes)");axes[1].legend(loc="upper left",frameon=False,fontsize=9)
fig.text(.5,.035,"Original RNA order; 10,000 proposals; bars: numerical 95% MC error.",ha="center",fontsize=9)
fig.subplots_adjust(left=.15,right=.98,top=.89,bottom=.23,wspace=.20)
fig.savefig(a.output/"FigureS5_bio06.pdf");fig.savefig(a.output/"FigureS5_bio06.png",dpi=200)
print("Redrew BIO-06 figure from 22 supplied aggregate rows.")
