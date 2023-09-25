# Exx-Symmetry
记一次完不成的约定，一个永不再问的问题。
> ——那等我解决了symmetry=1的问题再来问你
>
> ——那希望我到时能有别的回答...
>
> ——嗯，拭目以待


本仓库是尝试为ABACUS中的EXX添加symmetry==1的支持时的试验代码。

---

尝试用群论给出的波函数关系还原出DM(R)：
$$S(gk)=\ket{\phi_{gk}}\bra{\phi_{gk}}=S^T(gk)$$
$$S(g,k)\equiv\ket{\phi_{gk}}\bra{g\phi_{k}}$$
群论给出的波函数关系：
$$c^T_{gk}\ket{\phi_{gk}}=\pm c^T_k\ket{g\phi_k}$$
右乘行向量$\bra{\phi_{gk}}$得到$c^T_{gk}S^T(gk)=c^T_kS^T(g,k)$
两边转置得
$$S(g,k)c_k=S(gk)c_{gk}, c_k=S^{-1}(g,k)S(gk)c_{gk}$$


然而，实现出$S(g,k)$，发现$k_1\neq k_2$时是零矩阵...
$$S_{\mu\nu}(g,k)=\int{}d\mathbf{r}\sum_\mathbf{R_1}\phi_\mu(\mathbf{r}-\tau_\mu-\mathbf{R}_1)e^{-i\alpha\mathbf{k}\cdot\mathbf{R_1}}\sum_\mathbf{R_2}\phi_{\nu}(\mathbf{r}-\tau_{\nu}-\mathbf{R}_2)e^{i\mathbf{k}\cdot\mathbf{R}_2}$$
——$\ket{\phi_{gk}}$和$\ket{g\phi_k}$是正交的？！

那么每个k点的信息都不可或缺，symmetry=1无解。


## However,  something may still be useful: 
- calculate kstar
- `Symmetry_Basic::atom_ordering_new`: preserve index for atom-map
- doc: https://xmywuqhxb0.feishu.cn/wiki/A7ETwz0wSiOZILk8yEac6IWGnWg
---
<p align="center">
    <img src="docs/abacus-logo.svg">
</p>

<p align="center">
    <a href="https://github.com/deepmodeling/abacus-develop/actions/workflows/image.yml">
        <img src="https://github.com/deepmodeling/abacus-develop/actions/workflows/image.yml/badge.svg">
    </a>
    <a href="https://github.com/deepmodeling/abacus-develop/actions/workflows/test.yml">
        <img src="https://github.com/deepmodeling/abacus-develop/actions/workflows/test.yml/badge.svg">
    </a>
</p>

<a id="readme-top"></a>

# About ABACUS

ABACUS (Atomic-orbital Based Ab-initio Computation at UStc) is an open-source package based on density functional theory (DFT). The package utilizes both plane wave and numerical atomic basis sets with the usage of norm-conserving pseudopotentials to describe the interactions between nuclear ions and valence electrons. ABACUS supports LDA, GGA, meta-GGA, and hybrid functionals. Apart from single-point calculations, the package allows geometry optimizations and ab-initio molecular dynamics with various ensembles. The package also provides a variety of advanced functionalities for simulating materials, including the DFT+U, VdW corrections, and implicit solvation model, etc. In addition, ABACUS strives to provide a general infrastructure to facilitate the developments and applications of novel machine-learning-assisted DFT methods (DeePKS, DP-GEN, DeepH, etc.) in molecular and material simulations.

# Online Documentation
For detailed documentation, please refer to [our documentation website](https://abacus.deepmodeling.com/).
