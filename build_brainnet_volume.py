# -*- coding: utf-8 -*-
"""为BrainNet Viewer生成统计体积.nii: AAL各脑区填网络的UPDRS-III标准化β。"""
import os; os.environ["MPLCONFIGDIR"]=r"G:/Temp/mpl"
import nibabel as nib, numpy as np, pandas as pd

AAL=r"D:/rest/Data/ExampleFiles/AAL90/aal.nii"   # AAL90, 标签1-90 (与模型一致)
BASE=r"F:/PD_brainage_data"; OUT=f"{BASE}/brainnet"; os.makedirs(OUT,exist_ok=True)

NET={"sns":[1,2,19,20,57,58,69,70],"fpn":[5,6,7,8,9,10,65,66],"drs":[3,4,59,60],
     "vnt":[11,12,13,14,63,64],"dft":[21,22,25,26,27,28,35,36,65,66,67,68,85,86],
     "slt":[29,30,31,32],"lng":[11,12,13,17,63],"adt":[79,80,81,82],
     "vsl":[43,44,45,46,47,48,49,50,51,52,53,54,55,56,89,90],
     "lmb":[15,16,23,24,33,34,39,40,83,84,87,88]}

img=nib.load(AAL); vol=np.round(img.get_fdata()).astype(int)
labs=np.unique(vol); print("AAL标签:",labs.min(),"-",labs.max(),"个数",len(labs))

def build(src_csv, outname, sig_only, absval=False):
    U=pd.read_csv(f"{BASE}/{src_csv}").set_index("r")
    beta={k:U.loc[k,"beta"] for k in NET}; fdr={k:U.loc[k,"fdr"] for k in NET}
    nets=[k for k in NET if (fdr[k]<0.05)] if sig_only else list(NET)
    out=np.zeros(vol.shape,np.float32)
    for L in range(1,91):
        bs=[beta[k] for k in nets if L in NET[k]]
        if bs:
            val=max(bs,key=abs)                  # 取绝对值最大的网络
            out[vol==L]=abs(val) if absval else val
    # 填补孔洞: 仅向背景膨胀2次(不覆盖已有脑区), 消除表面映射缝隙
    from scipy.ndimage import grey_dilation
    for _ in range(2):
        dil=grey_dilation(out,size=(3,3,3)); m=(out==0)&(dil!=0); out[m]=dil[m]
    nib.save(nib.Nifti1Image(out,img.affine,img.header),f"{OUT}/{outname}")
    print(f"  {outname}: 覆盖脑区={sum(1 for L in range(1,91) if any(L in NET[k] for k in nets))}, "
          f"β范围 {out[out!=0].min():.2f}~{out[out!=0].max():.2f}, 网络={nets}")

print("UPDRS-III:")
build("results_final_updrs.csv","updrs_beta_all.nii",False)   # 全部10网络
build("results_final_updrs.csv","updrs_beta_sig.nii",True)    # 仅4显著网络
print("MoCA (存为正的效应强度|β|, 便于BrainNet显示; 方向为负: 脑龄老->认知低):")
build("results_final_moca.csv","moca_beta_sig.nii",True,absval=True)
print("PD vs HC (存|β|, 方向为负: PD更年轻; 显著=感觉运动/边缘):")
build("results_final_pdvshc.csv","pdvshc_beta_sig.nii",True,absval=True)
print("\n输出目录:",OUT)
