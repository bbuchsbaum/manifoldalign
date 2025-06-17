```
RESEARCH ARTICLE
```
# Kernel Manifold Alignment for Domain

# Adaptation

```
Devis Tuia^1 *, Gustau Camps-Valls^2
1 MultiModal Remote Sensing, University of Zurich, Zurich, Switzerland, 2 Image Processing Laboratory,
Universitat of València, València, Spain
```
```
*devis.tuia@geo.uzh.ch
```
## Abstract

```
The wealth of sensory data coming from different modalities has opened numerous opportu-
nities for data analysis. The data are of increasing volume, complexity and dimensionality,
thus calling for new methodological innovations towards multimodal data processing. How-
ever, multimodal architectures must rely on models able to adapt to changes in the data dis-
tribution. Differences in the density functions can be due to changes in acquisition
conditions (pose, illumination), sensors characteristics (number of channels, resolution) or
different views (e.g. street level vs. aerial views of a same building). We call these different
acquisition modesdomains, and refer to the adaptation problem asdomain adaptation.In
this paper, instead of adapting the trained models themselves, we alternatively focus on
finding mappings of the data sources into a common, semantically meaningful, representa-
tion domain. This field ofmanifold alignmentextends traditional techniques in statistics
such as canonical correlation analysis (CCA) to deal with nonlinear adaptation and possibly
non-corresponding data pairs between the domains. We introduce a kernel method for man-
ifold alignment (KEMA) that can match an arbitrary number of data sources without needing
corresponding pairs, just few labeled examples in all domains. KEMA has interesting prop-
erties: 1) it generalizes other manifold alignment methods, 2) it can align manifolds of very
different complexities, performing a discriminative alignment preserving each manifold
inner structure, 3) it can define a domain-specific metric to cope with multimodal specifici-
ties, 4) it can align data spaces of different dimensionality, 5) it is robust to strong nonlinear
feature deformations, and 6) it is closed-form invertible, which allows transfer across-
domains and data synthesis. To authors’knowledge this is the first method addressing all
these important issues at once. We also present a reduced-rank version of KEMA for
computational efficiency, and discuss the generalization performance of KEMA under
Rademacher principles of stability. Aligning multimodal data with KEMA reports outstanding
benefits when used as a data pre-conditioner step in the standard data analysis processing
chain. KEMA exhibits very good performance over competing methods in synthetic con-
trolled examples, visual object recognition and recognition of facial expressions tasks.
KEMA is especially well-suited to deal with high-dimensional problems, such as images
and videos, and under complicated distortions, twists and warpings of the data manifolds. A
fully functional toolbox is available athttps://github.com/dtuia/KEMA.git.
```
```
OPENACCESS
```
Citation:Tuia D, Camps-Valls G (2016) Kernel
Manifold Alignment for Domain Adaptation. PLoS
ONE 11(2): e0148655. doi:10.1371/journal.
pone.

Editor:Zhaohong Deng, Jiangnan University, CHINA

Received:October 23, 2015

Accepted:January 21, 2016

Published:February 12, 2016

Copyright:© 2016 Tuia, Camps-Valls. This is an
open access article distributed under the terms of the
Creative Commons Attribution License, which permits
unrestricted use, distribution, and reproduction in any
medium, provided the original author and source are
credited.

Data Availability Statement:Access to all relevant
data is referenced within the paper. we provide fully
functional code URLs:http://isp.uv.es/code/KEMA.
htm https://github.com/dtuia/KEMA.git. These links
contain both the code and the data used in the
experiments.

Funding:DT was supported by Swiss National
Science Foundation, project PP00P2-150593 (http://
[http://www.snf.ch/en/Pages/default.aspx).](http://www.snf.ch/en/Pages/default.aspx).) GCV was
supported by Spanish Ministry of Economy and
Competitiveness, project TIN2012-38102-C03-
(http://www.mineco.gob.es/portal/site/mineco?lang_
choosen = en) and European research council,
project 647423 (http://erc.europa.eu/). The funders
had no role in study design, data collection and


### Introduction

```
Domain adaptation constitutes a field of high interest in pattern analysis and machine learn-
ing. Classification algorithms developed with data from one domain cannot be directly used
in another related domain, and hence adaptation of either the classifier or the data represen-
tation becomes strictly imperative [ 1 ]. For example, there is actually strong evidence that a
significant degradation in the performance of state-of-the-art image classifiers is due to test
domain shifts, such as changing image sensors and noise conditions [ 2 ], pose changes [ 3 ],
consumer vs. commercial video [ 4 ], and, more generally, datasets biased due to changing
acquisition procedures [ 5 ].
Adapting (modifying) the classifier for any new incoming situation requires either compu-
tationally demanding retraining, passive-aggressive strategies, online filtering, or sample-rele-
vance estimation and weighting. These approaches are algorithm-dependent, often resort to
heuristic parameters, require good estimates of sample relevance and information content. The
ever-evolving classifier is also very hard to analyze. Alternatively, one may also try to adapt the
domain representations to a single latent space, and then apply a unique single classifier in that
semantically meaningfulfeature space. In this paper, we focus on the latter pathway. Adapting
the representation space has been referred in the literature to asfeature representation transfer
[ 6 ]orfeature transformation learning[ 7 ].
```
### Related works

```
The literature of feature representation transfer can be divided into three families of adaptation
problems, depending on the availability of labels in the different domains. They are briefly
reviewed hereafter and their main properties are summarized inTable 1. We discuss on the
```
Table 1. Properties of feature representation transfer methods.

```
DA type Properties
```
Method Unsup. Semis. Sup. D 2 Unpaired dS 6 ¼dT Nonlinear

PCA [ 24 ] p p
KPCA [ 25 ] p pp
TCA [ 10 ] p pp
SSTCA [ 10 ] ppp
JDA [ 26 ] p ^y pp
CCA [ 22 ] ppp
kCCA [ 9 ] pppp
MA [ 20 ] p p
GM [ 27 ] p p
OT-lab [ 15 ] ppp
SGF [ 12 ] ppp p p
GFK [ 13 ] ppp p p
MMDT [ 18 ] pp
SSMA [ 23 ] pp p p
KEMA pp p p p

D: number of domains.
dS,dT: number of features in source and target.
^y: semilabels predicted by a classifier.
p: known corresponding samples, but no labels.

doi:10.1371/journal.pone.0148655.t

analysis, decision to publish, or preparation of the
manuscript.

Competing Interests:The authors have declared
that no competing interests exist.


main the type of domain adaptation method (supervised, unsupervised or semisupervised), the
capability to align several domains or possibly unpaired examples, and eventually of different
dimensionality, and the linear or nonlinear nature of the transformation.
Unsupervised adaptation. First attempts of unsupervised domain adaptation are found in
multiview analysis [ 8 ], and more precisely in canonical correlation analysis (CCA) and kernel
CCA (KCCA) [ 9 ]. Despite their good performance in general, they still require points in differ-
ent sources to be corresponding pairs, which is often hard to meet in real applications. Think,
for example, of exploiting images, text and video in Wikipedia for document categorization,
trying to align images with different geometrical resolutions containing similar (not necessarily
the same) objects, or comparing commercial product images with consumer snapshots of the
same product. These real applications seldom provide datasets with corresponding pairs and/
or features. Alternative methods seek for a set of projectors that minimize a measure of discrep-
ancy between the source and target data distributions, such as the Maximum Mean Discrep-
ancy (MMD) [ 10 ] or the recent geodesic distance between distributions [ 11 ]. However, to
compare distributions, the data are supposed to be represented by the same features in all
domains. The idea of exploiting geodesic distances along manifolds was also considered in
[ 12 ], where a finite set of intermediate transformed data distributions are sampled along the
geodesic flow (SGF) between the linear subspaces. The intermediate features are then used to
train the classifier. The idea was extended in [ 13 ], where a Geodesic Flow Kernel (GFK) was
constructed by considering the infinity of transformed subspaces along the geodesic path.
However, both SGF and GFK assume input data space of the same dimensionality.
Semi-supervised adaptation with labels in the source domain only. A second family of
methods exploits the wealth of unsupervised information along with the limited amount of
labeled data in the source domain to guide the adaptation. Actually, some of the above-men-
tioned methods can incorporate the information of labeled samples in the source domain: the
Transfer Component Analysis [ 10 ] becomes semi-supervised by maximizing the Hilbert-
Schmidt Independence Criterion (HSIC) [ 14 ] between a kernel on features and a kernel on
labels in the source domain, while SGF [ 12 ] and GFK [ 13 ] become semi-supervised if the eigen-
vectors of the source domain are found with a discriminative feature extractor such as partial
least squares (PLS). Another family of methods, collectively known as Optimal Transport (OT)
techniques, can use labeled samples in the source domain to maximize coherence in the trans-
portation plan of masses between source and target domains [ 15 ]. For this last method, the
transformation is defined such that the transformed source distribution has ideally the same
probability density as the target one, and simultaneously the labeled examples in the source
domain remain grouped together.
Supervised adaptation with labels in all domains. SGF and GFK can be also defined for
the case in which all the domains are labeled. Saenko et al. [ 2 ] learned transformations between
the domains as a dot product between the (linearly) transformed source samples. The method
was extended in [ 16 ] to domains of different dimensionality, and in [ 17 ] to problems with mul-
tiple domains. Alternative approaches try to align target and source features while simulta-
neously moving labeled examples to the correct side of a decision hyperplane (MMDT) [ 18 ].
Donahue et al. extended this reasoning by including Laplacian regularization [ 19 ]. A last family
of supervised methods is known asmanifold alignment, and aims at concurrently matching the
corresponding instances while preserving the topology of each input domain, generally using a
graph Laplacian [ 20 , 21 ]. Roughly speaking, aligning data manifolds reduces to finding projec-
tions to a common latent space where all datasets show similar statistical characteristics. Mani-
fold alignment (MA) is a new form of multivariate analysis that dates back to the work of
Hotelling in 1936 on canonical correlation analysis (CCA) [ 22 ], where projections try to corre-
late the data sources onto a common target domain. While appealing, these methods still


require specifying a small amount of cross-domain sample correspondences. The problem was
addressed in [ 23 ] by relaxing the constraint of paired correspondences with the constraint of
having the same class labels in all domains. The semi-supervised manifold alignment (SSMA)
method proposed in [ 23 ] projects data from different domains to a latent space where samples
belonging to the same class become closer, those of different classes are pushed far apart, and
the geometry of each domain is preserved. The method performs well in general and can deal
with multiple domains of different dimensionality. However, SSMA cannot cope with strong
nonlinear deformations and high-dimensional data problems.

### Contributions

This paper introduces a generalization of SSMA through kernelization for manifold alignment
and domain adaptation. The proposed Kernel Manifold Alignment (KEMA) has some remark-
able appealing properties:

1. KEMA generalizes other manifold alignment methods. Being a kernel method, KEMA
    reduces to SSMA [ 23 ] when using a linear kernel, thus allowing to deal with high-dimen-
    sional data efficiently in the dual form (Q-mode analysis): therefore KEMA can cope with
    input space of very large dimension, e.g. extracted by Fisher vectors or deep features. KEMA
    also generalizes other manifold alignment methods, e.g. [ 20 ] when used with a linear kernel
    and with sample correspondences instead of the class similarity matrices (see page 5);
2. KEMA goes beyond data rotations and can align manifolds of very different structure, per-
    forming a flexible discriminative alignment that preserves the manifold structure;
3. KEMA defines a domain-specific metric when using different kernel functions in the differ-
    ent domains. Contrarily to SSMA, KEMA can use different kernels in each domain, thus
    allowing to use the best descriptor for each data source at hand, e.g. when aligning text and
    images one could involve using (more appropriate) string or histogram kernels in the very
    same alignment procedure, or using the same kernel function with different hyperpara-
    meters in each domain;
4. As SSMA, KEMA can align data spaces of different dimensionality. This is an advantage
    with respect to other feature representation transfer approaches that require either sam-
    ple correspondences [ 9 , 12 , 15 , 20 ] or strict equivalence of the feature spaces across
    domains [ 2 , 10 , 25 ].
5. KEMA is robust to strong (nonlinear) deformations of the manifolds to be aligned, as the
    kernel compensates for problems in graph estimation and numerical problems. As noted
    above, the use of different metric stemming from different kernels reinforces the flexibility
    of the approach;
6. Mapping data between domains (and hence data synthesis) can be performed in closed-
    form, thus allowing to measure the quality of the alignment in physical units. Kernelization
    typically makes the new method not invertibleanalytically, and one commonly resorts to
    approximate methods for estimating pre-images [ 28 – 30 ]. For the case of KEMA, this is not
    straightforwad (see page 8). As an alternative, we propose a chain of transforms of different
    types as a simple, yet efficient way of performing the inversion accurately and in closed
    form.
    The reported theoretical advantages translate into outstanding convenience when working
with high-dimensional problems and strong distortions in the manifold structures, as illus-
trated on a large set of synthetic and real applications in the experimental section.


### Materials and Methods

In this section, we first recall the linear SSMA algorithm and then derive our proposed KEMA.
We discuss its theoretical properties, the stability bounds and propose a reduced rank algo-
rithm, as well as a closed-form inversion strategy.

### Semi-supervised manifold alignment

Semi-supervised learning consists in developing inference models that collectively incorporate
labeled and unlabeled data in the model definition. In semi-supervised learning (SSL) [ 31 ], the
algorithm is provided with some availablelabeledinformation in addition to theunlabeled
information, thus allowing to encode some knowledge about the geometry and the shape of the
dataset. There is an overwhelming amount of SSL methods in the literature, yet the vast major-
ity of algorithms try to encode the relations between labeled and unlabeled data through the
definition of an undirected graph, and more precisely through the graph Laplacian matrixL.
To defineL, let’s first define a graphG(V,E) with a set ofnnodes,V, connected by a set of
edges,E. The edge connecting nodesiandjhas an associated weight [ 31 ]. In this framework,
the nodes are the samples, and the edges represent the similarity among samples in the dataset.
A proper definition of the graph is the key to accurately introduce data structure in the model.
To understand how matrixLis constructed, two mathematical tools have to be introduced
[ 31 , 32 ]: First, theadjacencymatrixW, which contains the neighborhood relations between
samples. It has non-zero entries only between neighboring samples, which are generally found
byk-nearest neighbors or an-ball distance. Then, thedegreematrixD, which is a diagonal
matrix of sizen×ncontaining the number of connections to a node (degree). The Laplacian
matrixLis then defined asL=D−W. Intuitively,Lmeasures the variation (i.e. norm of deriv-
atives hence the name of Laplacian operator) of the decision function along the graph built
upon all (labeled and unlabeled) samples [ 31 ].
When it comes to manifold alignment, an interesting semisupervised approximation was
presented in [ 23 ]. Let us considerDdomainsXirepresenting similar classification problems.
The corresponding data matrices,Xi 2 Rdini,i=1,...,D, containniexamples (labeled,li, and

unlabeled,ui, withni=li+ui) of dimensiondi, andn¼

#### PD

i¼ 1 ni. The SSMA method [^23 ] maps
all the data to a latent spaceFsuch that samples belonging to the same class become closer,
those of different classes are pushed far apart, and the geometry of the data manifolds is pre-
served. Therefore, three entities have to be considered, leading to threen×nmatrices: 1) a sim-
ilarity matrixWsthat has componentsWsij¼ 1 ifxiandxjbelong to the same class, and 0

otherwise (including unlabeled); 2) a dissimilarity matrixWd, which has entriesWijd¼ 1 ifxi
andxjbelong to different classes, and 0 otherwise (including unlabeled); and 3) a similarity
matrix that represents the topology of a given domain,W, e.g. a radial basis function (RBF)
kernel or aknearest neighbors graph computed for each domain separately and joined in a
block-diagonal matrix. Since we are not interested in preserving geometrical similarity between
the domains (we are only interested in preserving their inner geometry), all the elements of the
off-diagonal blocks in the matrixWare zeros. On the contrary,WsandWdare defined between
the domains and therefore act as registration anchor points in the feature space. An illustrative
example of how SSMA works in given inFig 1. The three different entities lead to three differ-
ent graph Laplacians:Ls,Ld, andL, respectively. Then, the SSMA embedding must minimize a
joint cost function essentially given by the eigenvectors corresponding to the smallest non-zero
eigenvalues of the following generalized eigenvalue problem:

```
ZðLþmLsÞZ>V¼lZLdZ>V; ð 1 Þ
```

whereZis a block diagonal matrix containing the data matricesXi,Z= diag(X 1 ,,XD), andV
contains in the columns the eigenvectors organized in rows for the particular domain,V=[v 1 ,
v 2 ,...,vD]>, see details in [ 21 , 33 ]. The method allows to extract a maximum ofNf¼

#### PD

i¼ 1 di
features that serve for projecting the data to the common latent domain as follows:

```
PFðXiÞ¼v>iXi: ð 2 Þ
```
Advantageously, SSMA can easily project data between domainsjandi: first mapping the
data inXjto the latent domainF, and from there inverting back to the target domainXias fol-

lows:

```
PiðXjÞ¼ðvjvyiÞ>Xj; ð 3 Þ
```
where†represents the pseudo-inverse of the eigenvectors of the target domain. The operation
is depicted as:

Therefore, the method can be used for domain adaptation but also for data synthesis. This
property was pointed out in [ 23 ], and experimentally studied for image analysis in [ 34 ].

Fig 1. The idea behind semi-supervised manifold alignment.(A) Consider two data sources (red and
black small points) in a binary problem (labeled points in orange balls and blue squares). SSMA aligns the
dataset by (B) preserving their inner geometry and (C) registering the data clouds in the feature space using
labels. (D) After alignment the datasets live in a semantically meaningful space.

doi:10.1371/journal.pone.0148655.g


### Kernel manifold alignment

When using linear algorithms, a well-established theory and efficient methods are often avail-
able. Kernel methods exploit this fact by embedding the data setSdefined over the input or
attribute spaceXðSXÞinto a higher (possibly infinite) dimensional Hilbert spaceH,orfea-
ture space, and then they build a linear algorithm therein, resulting in an algorithm which is
nonlinear with respect to the input data space. The mapping function is denoted as:X!H.
Though linear algorithms will benefit from this mapping because of the higher dimensionality
of the feature space, the computational load would dramatically increase because we should
compute sample coordinates in that high dimensional space. This computation is avoided
through the use of the kernel trick by which, if an algorithm can be expressed with dot products
in the input space, its (nonlinear) kernel version only needs the dot products among mapped
samples. Kernel methods compute the similarity between training samplesS¼fxigni¼ 1 using
pair-wise inner products between mapped samples, and thus the so-called kernel matrixKij=K
(xi,xj)=hφ(xi),φ(xj)icontains all the necessary information to perform many classical linear
algorithms in the feature space.
Kernelization of SSMA. Kernelization of SSMA is apparently straightforward; one should
map the data to a Hilbert feature space and then replace all instances of dot products with ker-
nel functions. However, note that in the original formulation of SSMA, there areDdata sources
that need to be first mapped to a common feature space. For doing this, we need to defineD
different feature mappings to eventually different Hilbert feature spaces, and then ensure that
mapped data live in the same subspace in order to do linear operations therein withallmapped
data sources. This can be actually done by resorting to a property of Functional Analysis The-
ory [ 35 ], thedirect sum of Hilbert spaces.
Theorem 1 Direct sum of Hilbert spaces [ 35 ]:Given two Hilbert spaces,H 1 andH 2 ,the set
of pairs{x,y}withx 2 H 1 andy 2 H 2 is a Hilbert spaceHwith inner product
hfx 1 ;y 1 g;fx 2 ;y 2 gi ¼ hx 1 ;x 2 iH 2 þhy 1 ;y 2 iH 2 .This is called the direct sum of the spaces, and is
denoted asH¼H 1 H 2 .This property extends to afinite summation of D Hilbert spaces by
whichH¼Di¼ 1 Hiis a Hilbert space.
Now we have the necessary tools for kernelizing the SSMA algorithm. Let us first map theD
different datasets toDpossibly different Hilbert spacesHiof dimensionHi,
iðÞ:x 7 !iðxÞ2Hi,i=1,...,D. Now, by replacing all the samples with their mapped fea-
ture vectors, the problem becomes:

```
ΦðLþmLsÞΦ>U¼lΦLdΦ>U; ð 4 Þ
```
whereFis a block diagonal matrix containing the data matricesFi=[φi(x 1 ),...,φi(xni)]>and
Ucontains the eigenvectors organized in rows for the particular domain defined in Hilbert

spaceHi,U=[u 1 ,u 2 ,...,uH]>whereH¼

#### PD

iHi. Note that the eigenvectorsuiare of possibly
infinite dimension and cannot be explicitly computed. Instead, we resort to the definition ofD
corresponding Riesz representation theorems [ 36 ] so the eigenvectors can be expressed as a lin-
ear combination of mapped samples [ 37 ],ui=Fiαi, and in matrix notationU=FΛ. This
leads to the problem:

```
ΦðLþmLsÞΦ>ΦΛ¼lΦLdΦ>ΦΛ: ð 5 Þ
```
Now, by pre-multiplying both sides byF>and replacing the dot products with the correspond-
ing kernel matrices,Ki¼Φ>iΦi, we obtain thefinal solution:

```
KðLþmLsÞKΛ¼lKLdKΛ; ð 6 Þ
```

whereKis a block diagonal matrix containing the kernel matricesKi. Now the eigenproblem
becomes of sizen×ninstead ofd×d, and we can extract a maximum ofNf=nfeatures.
When a linear kernel is used for all the domains,Ki¼X>iXi, KEMA reduces to SSMA:

```
PFðXiÞ¼α>iX>iXi¼ðXiαiÞ>Xi¼v>iXi: ð 7 Þ
```
This dual formulation is advantageous when dealing with very high dimensional datasets,di
nifor which the SSMA problem is not well-conditioned. Operating inQ-mode endorses the
method with numerical stability and computational efficiency in current high-dimensional
problems, e.g. when using Fisher vectors or deep features for data representation. This type of
problems with much more dimensions than points are recurrent nowadays for example in the
fields of bioinformatics, chemometrics, and image and video processing. In this sense, even
KEMA with a linear kernel becomes a valid solution for these problems, as it has all the advan-
tages of CCA-like methods, but can also deal with unpaired data.
Projection to the latent space requires first mapping the dataXito its corresponding Hilbert
spaceHi, thus leading to the mapped dataFi, and then applying the projection vectorui
defined therein:

```
PFðXiÞ¼u>iΦi¼α>iΦ>iΦi¼α>iKi: ð 8 Þ
```
which can be depicted as:

Therefore, projection to the kernel latent space is possible through the use of dedicated repro-
ducing kernel functions.
In order to map data from domainXjto domainXiwith KEMA we would need to estimate
D−1 inverse mappings from the latent space to the corresponding target domainXi. Such
transformations are highly desirable in order to measure the accuracy of the alignment/adapta-
tion in meaningful physical units. In general, nevertheless, using kernel functions hampers the
invertibility of the transformation. One can show that if an exact pre-image exists, and if the
kernel can be written askðx;x^0 Þ¼ckðx>x^0 Þwith an invertible functionψk(), then one can
compute the pre-image analytically under mild assumptions. However, it is seldom the case
that exact pre-images exist, and one resorts to approximate methods such as those in [ 28 – 30 ].
In the case of KEMA, inversion from the latent space to the target domainXiis even harder,
and hampers the use of standard pre-imaging techniques. Standard pre-image methods in ker-
nel machines [ 28 – 30 ] typically assume a particular kernel method (e.g. kPCA) endorsed with a
particular kernel function (often the polynomial or the squared exponential). If other kernel
functions are used, the formulation should be derived again. Remember that our KEMA feature
vectors in the latent space were obtained using a complex (and supervised) function that con-
siders labeled and unlabeled samples from all available domains through the composition of
kernel functions and graph Laplacians. One could derive the equations for preimaging under
our eigenproblem setting,K^0 ≔Ks^1 KdwhereKs≔K(L+μLs)KandKd=KLdK, but this is
very complicated, data dependent, and sensitive because of the appearance of several hyper-
parameters. Another alternative could be performing a sort of multidimensional regression
(from the latent space toXi) in a similar way to the kernel dependency estimation (KDE)
method revised in [ 29 ], but the approach would be complicated (no guarantees about the exis-
tence of a kernel trying to reproduce the inverse mapping implicit inK^0 exist), computationally
demanding (many hyperparameters appear), and would not deliver a closed-form solution.


Here we propose a simple alternative solution to the mapping inversion: to use a linear ker-
nel for the latent-to-target transformationKi¼X>iXi, andKjforj 6 ¼iwith any desired form.
Following this intuition, projection of dataXjto the target domainibecomes:

```
PiðXjÞ¼ðuyiÞ>α>jKj¼ðαjðXiαiÞyÞ>Kj; ð 9 Þ
```
where for the target domain we usedui=Fiαi=Xiαi. We should note that the solution is not
unique sinceDdifferent inverse solutions can be obtained depending on the selected target
domain. Using different transforms to perform model inversion was also recently studied in
[ 38 ]: here, instead of using an alternate scheme, we perform direct inversion by chaining differ-
ent transforms, leading to an efficient closed-form solution. Such a simple idea yields impres-
sive results in practice (see the experimental section, page 14).

### Computational efficiency and stability of KEMA

One of the main shortcomings of KEMA is related to the computational cost since twon×n

kernel matrices are involved, beingn¼

#### PD

i¼ 1 ni. KEMA complexity scales quadratically withn
in terms of memory, and cubically with respect to the computation time. Also projection for
new data requires the evaluation ofnkernel functionsperexample, becoming computationally
expensive for largen. To alleviate this problem, we propose two alternatives to speed up
KEMA: a reduced-rank approximation (REKEMA) and a randomized features approximation
(rKEMA). We compare both approaches in CPU time, and for rKEMA we study the conver-
gence bound inℓ 2 -norm based on matrix Bernstein inequalities. Finally, we study the stability
of the obtained solution when solving a (regularized) generalized eigenproblem using afinite
number of samples based on Rademacher principles.

### Reduced rank approximation

The so-called reduced-rank Kernel Manifold Alignment (REKEMA) formulation imposes
reduced-rank solutions for the projection vectors,W=FrΛ, whereFris a subset of the train-
ing data containingrsamples (rn) andΛis the new argument for the maximization prob-
lem. PluggingWintoEq (5), and replacing the dot products with the corresponding kernels,
Krn¼Φ>rΦ, we obtain thefinal solution:

```
KrnðLþmLsÞKnrΛ¼lKrnLdKnrΛ; ð 10 Þ
```
whereKrnis a block diagonal matrix containing the kernel matricesKicomparing a reduced
set ofrrepresentative vectors andalltraining data points,n. REKEMA reports clear benefits
for obtaining the projection vectors (the eigenproblem becomes of sizer×rinstead ofn×n),
hence the computational cost becomesOðr^3 Þ,rn, compacting the solution (nowNf=rn
features), and in storage requirements (henceOðr^2 Þ). We want to highlight here that this is not
a simple subsampling, because the model considers correlations between all training data and
the reduced subset throughKrn. The selection of therpoints can be done in different ways and
degrees of sophistication: close to centroids provided by a pre-clustering stage, extremes of the
convex hull, sampling to minimize the reconstruction error or preserve information, form
compact basis in feature space, etc. While such strategies are crucial in low-to-moderate sam-
ple-size regimes, random selection offers an easy way to select therpoints and is the most
widely used strategy.Fig 2Ashows the evolution of the computational cost as a function of
(randomly selected)rsamples in a toy example of aligning two spirals (cf. experiment]1 in the
experiments section).


### Random features approximation

A recent alternative to reduced rank approximations exploits the classical Bochner’s theorem
in harmonic analysis, which has been recently introduced in the field of kernel methods [ 39 ].
The Bochner’s theorem states that a continuous kernelk(x,y)=k(x−y)onRdis positive defi-
nite (p.d.) if and only ifkis the Fourier transform of a non-negative measure. If a shift-invari-
ant kernelkis properly scaled, its Fourier transformp(w) is a proper probability distribution.
This property is used to approximate kernel functions and matrices with linear projections on
mrandom features as follows:

```
kðx;yÞ¼
```
#### Z

```
Rd
```
```
pðwÞejw>ðxyÞdw
```
```
Pm
```
i¼ (^11) mejw
>ixejw>iy

#### ¼

```
Pm
i¼ 1
```
#### 1

```
mcosðw
```
>ixþbiÞcosðw>iyþbiÞ¼h (^1) ffiffiffi
pmzðxÞ;
(^1) ffiffiffi
pmzðyÞi;
ð 11 Þ
wherep(w) is set to be the inverse Fourier transform ofkandbiUð 0 ; 2 pÞ[ 39 ]. Therefore, we
can randomly sample parameterswi 2 Rdfrom a data-independent distributionp(w) and con-
struct am-dimensional randomized feature mapz():X!Z, for dataX 2 Rndand
Z 2 Rnm, as follows:
w 1 ;...;wmpðwÞ;
zi ≔½cosðw>ix 1 þbiÞ;...;cosðw>ixnþbiÞ 2Rn;
zðXÞ ≔Z¼½z 1 zm2Rnm:
ð 12 Þ
For a collection ofndata points,fxigni¼ 1 , a kernel matrixK 2 Rnncan be approximated with
the explicitly mapped data,Z 2 Rnm,K^ZZ>. The Gaussian kernelkðx;yÞ¼expðkx
yk^22 =ð 2 s^2 ÞÞcan be approximated usingwiNð 0 ;I=s^2 Þ. For the case of KEMA, we have to
sample twice, hence obtain two sets of vectors and associated matricesZsandZd,toapproxi-
mate the similarity and dissimilarity kernel matrices,Ks≔KðLþmLsÞKZsZ>sand
Kd≔KLdKZdZ>d. The associated cost by using the random features approximation now
reduces toOðnm^2 Þ,seeFig 2B. It is also important to notice that solving the generalized eigen-
value problem in KEMA feature extraction with random features converges inℓ 2 -norm error
Fig 2. Average computational cost of REKEMA andrKEMA.CPU time [s], over 10 realizations as a
function ofrandmfor the (A) reduced rank KEMA (REKEMA) and (B) randomized KEMA (rKEMA) in black
lines. In both figures, the red line is the KEMA solution). We used synthetic example]1 (see experiments
section) withn= 1000 samples.
doi:10.1371/journal.pone.0148655.g


```
withOðm^1 =^2 Þand logarithmically in the number of samples when using an appropriate ran-
dom parameter sampling distribution [ 40 ](seetheAppendix).
```
### Stability of KEMA

```
The use of KEMA in practice raises, however, the important question of the amount of data
needed to provide an accurate empirical estimate, and how the quality of the solution differs
depending on the datasets. Such results have been previously derived for KPCA [ 41 ] and KPLS
[ 42 ] and here we extend them to our generalized eigenproblem setting. We focus on the con-
centration of sums of eigenvalues of the generalized KEMA eigenproblem solved using a finite
number of samples, where new points are projected into them-dimensional space spanned by
themeigenvectors corresponding to the largestmeigenvalues.
Following the notation in [ 41 ], we refer to the projection onto a subspaceUof the eigenvec-
tors of our eigenproblem asPU(φ(x)). We represent the projection onto the orthogonal comple-
ment ofUbyPU?(φ(x)). The norm of the orthogonal projection is also referred to as the
residual since it corresponds to the distance between the points and their projections.
Theorem 2 (Th. 1 and 2 in [ 41 ])Let us defineKs≔K(L+μLs)KandKd=KLdK.If we
perform KEMA in the feature space defined byK^ ≔Ks^1 Kd,then with probability greater than
1 −δover n random samples S,for all 1 rn,if we project data on the spaceU^r,the expected
squared residual is bounded by
```
```
Xn
```
```
j¼rþ 1
```
```
ljEkPU^r?k^2
```
```
hi
min 1 lr
```
#### 1

```
n
```
```
Xn
```
```
j¼lþ 1
```
```
^ljðSÞþ^1 þ
```
```
ffiffi
l
```
```
p
ffiffiffi
n
```
```
p
```
```
ffiffiffiffiffiffiffiffiffiffiffiffiffiffiffiffiffiffi
2
n
```
```
Xn
```
```
i¼ 1
```
```
Kii^2
```
```
"#s
þR^2
```
```
ffiffiffiffiffiffiffiffiffiffiffiffiffiffiffiffiffiffiffiffiffiffiffi
18
n
```
```
ln
```
```
2 n
d
```
```
s 
ð 13 Þ
```
```
and
```
```
Xr
```
```
j¼ 1
```
```
ljEkPU^rk^2
```
```
hi
max 1 lr
```
#### 1

```
n
```
```
Xl
```
```
j¼ 1
```
```
^ljðSÞ^1 þ
```
```
ffiffi
l
```
```
p
ffiffiffi
n
```
```
p
```
```
ffiffiffiffiffiffiffiffiffiffiffiffiffiffiffiffiffiffi
2
n
```
```
Xn
```
```
i¼ 1
```
```
Kii^2
```
```
"#s
R^2
```
```
ffiffiffiffiffiffiffiffiffiffiffiffiffiffiffiffiffiffiffiffiffiffiffiffiffiffiffiffiffiffiffiffiffiffiffi
19
n ln
```
```
2 ðnþ 1 Þ
d
```
```
s 
;ð 14 Þ
```
```
where the support of the distribution is in a ball of radius R in the feature space andλiare^liare
the process and empirical eigenvalues, respectively.
Theorem 3 (Regularized KEMA)The previous theorem holds only when the inverseKs^1
exists. Otherwise, we typically resort to matrix conditioning via regularization. Among the many
possibilities in problem conditioning, the standard direct Tikhonov-Arnoldi approach helps solv-
ing the generalized eigenproblem on a shifted and inverted matrix, which damps the eigenvalues.
Now we aim to bound a well-conditioned matrixK^0 ≔(Ks+γKd)−^1 Kd,whereγ> 0 is the regu-
larization parameter. It is easy to show that its estimated eigenvalues,^yiare related to the unre-
gularized ones asl^j¼y^j=ð 1 g^yjÞ.Therefore, with probability greater than 1 −δover n
random samples S,for all 1 rn,if we project data on the spaceU^r,the expected squared
residual is bounded by
```
```
Xn
```
```
j¼rþ 1
```
```
ljEkPU^r?k^2
```
```
hi
min 1 lr
```
#### 1

```
n
```
```
Xn
```
```
j¼lþ 1
```
```
^yjðSÞ
1 gy^jðSÞ
```
```
þ
```
```
1 þ
```
```
ffiffi
l
```
```
p
ffiffiffi
n
```
```
p
```
```
ffiffiffiffiffiffiffiffiffiffiffiffiffiffiffiffiffiffi
2
n
```
```
Xn
```
```
i¼ 1
```
```
K^0 ii^2
```
```
"#s
þR^2
```
```
ffiffiffiffiffiffiffiffiffiffiffiffiffiffiffiffiffiffiffiffiffiffiffi
18
n
```
```
ln
```
```
2 n
d
```
```
s 
ð 15 Þ
```
```
and
```
Xr

```
j¼ 1
```
```
ljEkPU^rk^2
```
```
hi
max 1 lr
```
#### 1

```
n
```
```
Xl
```
```
j¼ 1
```
```
y^jðSÞ
1 g^yjðSÞ
```
#### 

```
1 þ
```
```
ffiffi
l
```
```
p
ffiffiffi
n
```
```
p
```
```
ffiffiffiffiffiffiffiffiffiffiffiffiffiffiffiffiffiffi
2
n
```
```
Xn
```
```
i¼ 1
```
```
Kii^02
```
```
"#s
R^2
```
```
ffiffiffiffiffiffiffiffiffiffiffiffiffiffiffiffiffiffiffiffiffiffiffiffiffiffiffiffiffiffiffiffiffiffiffi
19
n
```
```
ln
```
```
2 ðnþ 1 Þ
d
```
```
s 
;ð 16 Þ
```

where the support of the distribution is in a ball of radius R in the feature space,θiand^yiare the
process and empirical eigenvalues.
In either case, the lower bound confirms that a good representation of the data can be
achieved by using the firstffiffiffiffiffiffiffi reigenvectors if the empirical eigenvalues quickly decrease before

```
l=n
```
p
becomes large, while the upper bound suggests that a good approximation is achievable

for values ofrwhere

```
ffiffiffiffiffiffiffi
r=n
```
p
is small. These results can be used as a benchmark to test different
approaches or to select among possible candidate kernels. Also, note that depending on how
much non-diagonal isK^ (orK^0 ), i.e. how large are the manifold mis-alignments, the KEMA
bounds may be tighter than those of KPCA. With an appropriate estimation of the manifold
structures via the graph Laplacians and tuning of the kernel parameters, the performance of
KEMA will be at least asfitted as that of KPCA. Note that when intense regularization is
needed, the trace of the squaredK^0 can be upper bounded byn^1 g 2 and then the expected squared

residuals are mainly governed bynandγ.

### Results and discussion

We analyze the behavior of KEMA in a series of artificial datasets of controlled level of distor-
tion and mis-alignment, and on real domain adaptation problems of visual object recognition
from multi-source commercial databases and recognition of multi-subject facial expressions.

### Toy examples with controlled distortions and manifold mis-alignments

Setup. the first set of experiments considers a series of toy examples composed of two
domains with data matricesX 1 andX 2 , which are spirals with three classes (see the two first
columns ofFig 3). Each dataset is visualized by

- domain(first column ofFig 3): the first domain is characterized by a red square marker and
    the second by a blue cross. With this plot, we see if the domains are misaligned, irrespectively
    of the classes.
- class(second column ofFig 3): in this case, both domains are characterized by the class col-
    ors (red, green and blue circles). With this plot we see if the classes are aligned, irrespectively
    of the domain.
       Then, a series of deformations are applied to the second domain: scaling, rotation, inversion
of the order of the classes, the shape of the domain (spiral or line) or the data dimensionality
(seeTable 2). These experiments are designed to study the flexibility of KEMA to handle align-
ment problems of increasing complexity and between data of different dimensionality (Ex. #2).
The last experiment (#6) considers the same setting of Exp. #1, but adds 50 features of Gaussian
noise to the two informative features.
    For each experiment, 60 labeled pixelsperclass were sampled in each domain, as well as
1000 unlabeled samples that were randomly selected. Classification performance was assessed
on 1000 held-out samples from each domain. The toy classification results can be reproduced
using the MATLAB toolbox available athttps://github.com/dtuia/KEMA.git. Theσbandwidth
parameter of the RBF kernel was set in each domain as half of the median distance between all
the samples in the domain, thus enforcing a domain-specific metric in each domain.
    Latent space and domain adaptation. Fig 3illustrates the projections obtained by KEMA
when using a linear and an RBF kernel (lengthscale was set as the average distance between
labeled samples). Looking at the alignment results, we observe that the linear KEMAlinaligns
effectively the domains only in experiments #1 and #4, which are basically scalings and rota-
tions of the data. However, it fails on experiments #2, #3 and #5, where the manifolds have


undergone stronger deformations. The use of a nonlinear kernel (KEMARBF) allows much
more flexible solution, performing a discriminative transform plus alignment in all experi-
ments. In Experiment #6, even though the two discriminative dimensions (out of 52) are the
same as in Exp. #1, only KEMARBFcan align the data effectively, since KEMAlinis strongly
affected by the noise and returns a non-discriminative alignment for the eigenvectors corre-
sponding to the smallest eigenvalues.
Classification performances. Fig 4reports the classification errors obtained by a linear
discriminant analysis (LDA,Fig 4A) and the nearest neighbor classifier (1-NN,Fig 4B). For
each classifier, classification errors are reported for the samples from the source domain (left
inset) and the target domain (right inset). LDA is used to show the ability of projecting the

Fig 3. Illustration of linear and kernel manifold alignment on the toy experiments.(A) data in the original
domains (X1 is designated with red squares, X2 is designated with blue crosses) andperclass (red, green
and blue circles, respectively), data projected (B) with the linear and (C) the RBF kernels.

doi:10.1371/journal.pone.0148655.g


```
domains in a joint discriminative latent space, where even the simplest linear classifier can be
successful. 1-NN is used to show the increase in performance that can be obtained by using a
nonlinear, yet simple, classifier on top of the projected data.
When using a linear model (LDA), a large improvement of KEMARBFover KEMAlin(thus
over SSMA) is observed. In experiment #1, even if the alignment is correct (Fig 3), the linear
classifier trained on the projections of KEMAlincannot resolve the classification of the two
domains, while KEMARBFsolution provides a latent space where both domains can be classi-
fied correctly. Experiment #2 shows a different picture: the baseline error (green line inFig 4)
is much smaller in the source domain, since the dataset in 3D is linearly separable. Even if the
classification of this first domain (red square inFig 3) is correct for all methods, classification
after SSMA/KEMAlinprojection of the second domain (blue x inFig 3) is poor, since their pro-
jection in the latent space does not“unfold”the blue spiral. KEMARBFprovides the best result.
For experiment #3, the same trend as in experiment #2 is observed. Experiments #4 and #
deal with reversed classes (the brown class is the top one in the source domain and the bottom
one in the target domain). In both experiments, we observe a very accurate baseline (both
domains are linearly separable in their own input spaces), but only KEMARBFprovides the cor-
rect match in a low-dimensional latent space (2 dimensions), including a discriminative V-
shaped projection leading to nearly 0% errors on average; KEMAlinrequires 5 dimensions to
achieve a correct manifold alignment and a classification as accurate as the baseline (that still
includes misclassifications in the linear classifier). The missclassifications can be explained by
the projected space (3rd and 4th columns inFig 3), where classes are aligned at best, but no
real matching of the two data clouds is performed. The last experiment (#6) deals with noisy
data, where only two out of the 52 dimensions are discriminative: KEMARBFfinds the two first
eigenvectors that align the data accurately (classification errors close to 0% in both domains),
while KEMAlinshows a much noisier alignment that, due to the rigidity of a linear transform,
leads to about 20% misclassification in both domains.
When using the nonlinear 1-NN, both the KEMARBFand KEMAlinperform similarly.
KEMARBFstill leads to correct classification with close to zero errors in all cases, thus confirm-
ing that the latent space projects samples of the same class close. KEMAlinleads to correct clas-
sification in almost all the cases, since the 1-NN can cope with multimodal class distributions
and nonlinear patterns in the latent space. KEMAlinstill fails in Exp #3, where the projection of
the source domain (red circle inFig 3) stretches over the target domain, and in Exp. # 6, where
the latent space is not discriminative and harms the performance of the 1-NN.
Alignment with REKEMA. We now consider the reduced-rank approximation of KEMA
proposed. We used the data in the experiment #1 above.Fig 5illustrates the solutions of the stan-
dard SSMA (or KEMAlin), and for REKEMA using a varying rate of samples. We also give the
classification accuracies of a SVM (with both a linear and an RBF kernel) in the projected latent
```
Table 2. Specification of the toy examples.

Exp. Dimension Deformations Noisy

```
STShape ofS Scaling Rotation Classflip dimensions
#1 2 2 Spiral p -- 0
#2 3 2 Spiral - - - 0
#3 3 3 Line - - - 0
#4 3 3 Spiral - pp 0
#5 3 3 Spiral p - p 0
#6 52 52 Spiral p -- 50
```
doi:10.1371/journal.pone.0148655.t


```
space. Samples were randomly chosen and the sigma parameter for the RBF kernel in KEMARBF
was fixed to the average distance between all used labeled samples. We can observe that SSMA
successfully aligns the two domains, but we still need to resort to nonlinear classification to
achieve good results. REKEMA, on the contrary, essentially does two operations simultaneously:
aligns the manifolds and increases class separability. Excessive sparsification leads to poor results.
Virtually no difference between the full and the reduced-rank solutions are obtained for small
```
Fig 4. Classification performances on the toy examples.Error rates as a function of the extracted features (Nf) when predicting data for the first (left inset)
or the second (right inset) domain. In all plots KEMALinis in blue, KEMARBFin red, SSMA in cyan and the Baseline in green. Panel (A) shows the LDA results,
panel (B) the 1-NN.

doi:10.1371/journal.pone.0148655.g


values ofr: just 10% of examples are actually needed to saturate accuracies. The proposed
rKEMA showed similar behaviour but results are omitted for the sake of simplicity.
Invertibility of the projections. Fig 6shows the results of invertibility of SSMA and
KEMA (usingEq (9)) on the previous toy examples (we excluded Exp. # 6 to avoid synthesizing
data with 50 noisy dimensions). We use a linear kernel for the inversion part (latent-to-source)

Fig 5. Linear and kernel manifold alignment on the scaled interwined spirals toy experiment (Exp. #
inFig 3).REKEMA is compared to SSMA for different rates of training samples (we usedli= 100 andui=
per class for both domains).

doi:10.1371/journal.pone.0148655.g

Fig 6. Domain inversion with SSMA and KEMA.For each panel, the left inset represents domains: the red
squares are samples in the source domain, while the blue crosses are target domain samples projected onto
the source domain. The right inset represents the three classes (red, green and blue circles). Each plot shows
the result of a single run, and the averagedℓ 2 -norm reconstruction error over 10 runs.

doi:10.1371/journal.pone.0148655.g


and use for the direct part (target-to-latent space) an RBF kernel. All results are shown in the
source domain space. All the other settings (# labeled and unlabeled,μ, graphs) are kept as in
the experiments shown inFig 3. The reconstruction error, averaged on 10 runs, is also reported:
KEMARBF!linis capable of inverting the projections and is always as accurate as the SSMA
method in the simplest cases (#1, #4). For the cases related to higher levels of deformation,
KEMA is either as accurate as SSMA (#3, where the inversion is basically a projection on a
line) or significantly better: for experiment #2, where the two domain are strongly deformed,
and experiment #5, where we deal with both scaling and inverted classes, only KEMARBF!lin
can achieve satisfying inversion, as it“unfolds”the target domain and then only needs a rota-
tion to match the distribution in the source domain.

### Visual object recognition in multi-modal datasets

We here evaluate KEMA on visual object recognition tasks by using theOffice-Caltechdataset
introduced in [ 2 ]. We consider the four domains Webcam (W), Caltech (C), Amazon (A) and
DSLR (D), and selected the 10 common classes in the four datasets following [ 13 ]. By doing so,
the domains contain 295 (Webcam), 1123 (Caltech), 958 (Amazon) and 157 (DSLR) images,
respectively. The features were extracted in two ways

- SURF features, as described in [ 2 ]: we use a 800-dimensional normalized histogram of visual
    words obtained from a codebook constructed from a subset of the Amazon dataset on points
    of interest detected by the Speeded Up Robust Features (SURF) method. The features are
    included in the in the MATLAB package onhttps://github.com/dtuia/KEMA.git. Alterna-
    tively, they can be downloaded from their original repository onhttps://www.eecs.berkeley.
    edu/jhoffman/domainadapt/.
- Deep features from DeCAF [ 43 ]: these features are extracted as the sparse activations of the
    fully connected 7th layer of of a convolutional network trained on imageNet and then fine
    tuned on the visual recognition tasks considered here. It forms a 4096-dimensional vector.
    The features are included in the MATLAB package onhttps://github.com/dtuia/KEMA.git.

Experimental setup. We compare our proposed KEMA with the following unsupervised
and semi-supervised domain adaptation methods: GFK [ 13 ], OT-lab [ 15 ] and JDA [ 26 ]. We
used the same experimental setting as [ 13 ], in order to compare with these unsupervised
domain adaptation methods. For all methods, we used 20 labeled pixelsperclass in the source
domain for theC,AandWdomains and 8 samplesperclass for theDdomain. After alignment,
an ordinary 1-NN classifier was trained with the labeled samples. The same labeled samples in
the source domain were used to define the PLS eigenvectors for GFK and OT-lab. For all the
methods using labeled samples in the target domain (including KEMA), we used 3 labeled sam-
ples in target domain to define the projections.
We used a sensible kernel for this problem in KEMA: the (fast) histogram intersection ker-
nel [ 44 ]. Using aχ 2 kernel resulted in similar performances. We usedu= 300 unlabeled sam-
ples to compute the graph Laplacians, for which ak-NN graph withk= 21 was used.
Numerical results. The projections obtained by KEMA in the visual object recognition
experiments remain discriminative, as shown byFig 7, where projections on the first three
dimensions of the latent space are reported for theA!W(top) andC!A(bottom) using the
SURF features. The numerical results obtained in all the eight problems are reported in
Table 3: KEMA outperforms the unsupervised GFK and, in most of the cases, improves the
results obtained by the semi-supervised methods using labels in the source domain only.
KEMA provides the most accurate results in 5 out of the 8 settings. KEMA is as accurate as the
state of the art, but with the advantage of handling naturally domains of different


```
Fig 7. Example of the three first dimensions of the latent space.(A) illustrates theA!Wexperiment. (B)
illustrates theC!Aexperiment. Left: by domain (red circles are the source samples, blue crosses are the
target samples), right: by class (each color represents a different class).
doi:10.1371/journal.pone.0148655.g
```
Table 3. 1-NN classification accuracy in the visual object recognition study using the SURF features.

```
Train on source No adapt. Unsup. GFK [ 13 ]DA
Labels:S
```
```
Labels:S,TKEMAKint Train on target No adapt.
```
```
OT-lab [ 15 ] JDA [ 26 ]
```
lS 02020 20
lT 00 ^y 3
C!A 21.4±3.7 35.3±3.2 43.5±2.1 40.7±4.0 47.1±3.0 35.4±2.
C!D 12.3±2.8 35.6±5.0 41.8±2.8 40.0±4.0 61.5±2.8 65.1±1.
A!C 35.3±0.5 32.9±2.5 35.2±0.8 34.0±3.1 29.5±3.0 28.4±1.
A!W 31.0±0.7 32.0±3.4 38.4±5.4 36.0±5.1 65.4±2.7 63.5±2.
W!C 21.7±0.4 27.7±2.4 35.5±0.9 31.8±1.9 32.9±3.3 28.4±1.
W!A 27.0±1.5 33.3±2.1 40.0±1.0 31.5±4.7 44.9±4.5 35.4±2.
D!A 19.0±2.2 33.0±1.3 34.9±1.3 32.9±2.9 44.2±3.1 35.4±2.
D!W 37.4±3.0 69.7±3.8 84.2±1.0 80.0±4.1 64.1±2.9 63.5±2.
Mean 22.3 37.4 44.2 40.9 48.7 44.

C: Caltech,A: Amazon,D: DSLR,W: Webcam.
ldomain: number of labels per class.
^y: predicted labels.

doi:10.1371/journal.pone.0148655.t


```
dimensionality, and not requiring semilabeled examples (^yin the Table) to align the domains
as JDA. The results obtained when using the deep DeCAF features are reported inTable 4:a
strong improvement in performance is observed for all methods. This general increase was
expected, since the deep features in DeCAF are naturally suited for domain adaptation (they
are extracted withfine tuning on this specific dataset): but nonetheless, even if the boost in per-
formance is visible for all the methods (including the case without adaptation), KEMA
improves performances even further and leads to the best average results. Looking at the single
experiments, KEMA performs most often on a tie with OT-lab [ 15 ]. Summing up, KEMA
leads to results as accurate as the state of art, but is much more versatile, since it allows to han-
dle unpaired data, works with datasets of different dimensionality, and has a significantly
smaller computational load (see alsoTable 1for a taxonomical comparison of the properties of
the different methods).
```
### Recognition of facial expressions in multi-subject databases

```
This experiment deals with the task of recognizing facial expressions. We used the dataset in
[ 45 ], where 185 photos of three subjects depicting three facial expressions (happy, neutral and
shocked) are available. The features are included in the MATLAB package onhttps://github.
com/dtuia/KEMA.git. Alternatively, they can be downloaded from their original repository on
http://www.cc.gatech.edu/lsong/code.html. Each image is 217 × 308 pixels and we take each
pixel as one dimension for classification (66836 dimensional problem). Each pair {subject,
expression} has around 20 repetitions.
Experimental setup. Different subjects represent the domains and we align them with
respect to the three expression classes. We used only three labeled examplesperclass and sub-
ject, and held out 70% of the data for testing and used the remaining 30% (55 samples) for the
extraction of the labeled samples. The examples which have not been selected as labeled points
are used as unlabeled data. The three domains are aligned simultaneously into a common latent
space, and then all classifications are run therein for all subjects. Below, we report the results
```
Table 4. 1-NN classification accuracy in the visual object recognition study using the DeCAF fully connected layer fc7.

```
fc7 Train on source No adapt. Unsup. GFK [ 13 ]DA
labels:S
```
```
labels:S,TKEMAKint Train on target No adapt.
```
```
OT-lab [ 15 ] JDA [ 26 ]
```
lS 02020 20
lT 00 y^ 3
C!A 84.5±1.5 87.8±2.1 92.1±1.3 89.6±2.0 91.5±1.5 84.4±3.
C!D 73.1±4.9 83.5±3.6 85.4±6.0 85.0±4.9 93.6±3.1 92.2±1.
A!C 72.0±1.7 80.2±1.9 87.2±1.2 82.6±2.9 80.3±3.4 66.3±3.
A!W 61.3±3.4 78.0±4.8 84.5±2.4 83.0±4.6 92.7±2.5 88.1±3.
W!C 68.9±3.0 75.1±2.5 83.7±1.5 79.8±2.0 82.1±2.3 66.3±3.
W!A 73.5±2.7 81.2±2.2 91.9±1.4 90.9±1.2 91.6±1.3 84.4±3.
D!A 74.6±3.9 85.4±2.1 92.9±1.1 91.9±0.8 90.3±1.1 84.4±3.
D!W 93.8±1.5 96.7±1.9 94.1±3.4 97.0±1.5 91.0±3.5 88.1±3.
Mean 75.2 83.5 88.9 87.49 89.1 81.

C: Caltech,A: Amazon,D: DSLR,W: Webcam.
ldomain: number of labels per class.
^y: predicted labels.

doi:10.1371/journal.pone.0148655.t


obtained by using a LDA classifier trained in that common latent space. We consider three
experimental settings:

- Single resolution: all images are considered at their maximal resolution accounting for the
    three domains. Each domain is therefore a 66836-dimensional dataset. SSMA could not han-
    dle these data, since it would involve a 200508-dimensional eigendecomposition.
- Multiresolution, factor 2: the resolution of one of the domains (Subject #1) is downgraded by
    a factor two. 154 × 109, leading to a 16786-dimensional domain. The alignment problem in
    the primal would then be 16786 + (2 × 66836) = 150458-dimensional. With this experiment,
    we aim at showing the capability of KEMA to handle data of different dimensionality.
- Multiresolution, factor 4: the resolution of one of the domains (Subject #1) is downgraded by
    a factor four. 62 × 44, leading to a 2728-dimensional domain. The alignment problem in the
    primal would then be 136400-dimensional.

Numerical results. Average results over ten realizations are given inFig 8: since it works
directly in the dual, KEMA can effectively cast the three-domains problem into a low dimen-
sional space. In the single resolution case (Fig 8B) all domains are classified with less than 5%
error. This shows an additional advantage of KEMA with respect to SSMA in high dimensional
spaces: SSMA would have required to solve a 200508-dimensional eigenproblem, while KEMA
solves only a 55-dimensional problem. Subject #1 seems to be the most difficult to align with

Fig 8. Results of the classification of facial expressions (top: error rates, middle: predicted
expressions; bottom: subjects).(A) single resolution experiment; (B) multiresolution experiment with a
factor-two reduction for the images of subject 1; (C) multiresolution experiment with a factor-four reduction for
the images of subject 1.

doi:10.1371/journal.pone.0148655.g


the two others, difficulty that is also reflected in the higher classification errors. Actually, sub-
ject #1 shows little variations in his facial traits from one expression to the other compared to
the other subjects (see Fig 3 in [ 45 ]).
In the multi-resolution cases, similar error rates are observed for subjects #2 and #3, even
though the images of subject #1 were of coarse resolution. The reduced resolution of the images
of subject #1 made the expression recognition harder, but error rates lower than 20% are still
achieved by using KEMA. By looking at the projections (second and third rows ofFig 8), those
of the multiresolution experiment with a factor 2 reduction ((B) panel) are very similar to
those in the single resolution experiment ((A) panel).

### Conclusions

We introduced a kernel method for semi-supervised manifold alignment. We want to stress
that this particular kernelization goes beyond the standard academic exercise as the method
addresses many problems in the literature of domain adaptation and manifold learning. The
so-called KEMA can actually align an arbitrary number of domains of different dimensionality
without needing corresponding pairs, just few labeled examples in all domains. We also showed
that KEMA generalizes SSMA when using a linear kernel, which allows us to deal with high-
dimensional data efficiently in the dual form. Working in the dual can be computationally
costly because of the construction of the graph Laplacians and the size of the involved kernel
matrices. Regarding the Laplacians, they can be computed just once and off-line, while regard-
ing the size of the kernels, we introduced a reduced-ranked version that allows to work with a
fraction of the samples while maintaining the accuracy of the representation. Advantageously,
KEMA can align manifolds of very different structures and dimensionality, performing a dis-
criminative transform along with the alignment. We have also provided a simple yet effective
way to map data between domains as an alternative to standard pre-imaging techniques in the
kernel methods literature. This is an important feature that allows synthesis applications, but
more remarkably allows to study and characterize the distortion of the manifolds in physically
meaningful units. To the authors’knowledge this is the first method in addressing all these
important issues at once. All these features were illustrated through toy examples of increasing
complexity (including data of different dimensionality, noise, warps and strong nonlinearities)
and real problems in computer vision, and face recognition, thus showing the versatility of the
method and its interest for numerous application domains. It does not escape our attention
that KEMA may become a standard multivariate method for data preprocessing in general
applications where multisensor, multimodal, sensory data is acquired.

### Acknowledgments

This work has been supported by the Swiss National Science Foundation under project
PP00P2-150593, the Spanish Ministry of Economy and Competitiveness (MINECO) under
project TIN2012-38102-C03-01 (LIFE-VISION), and the ERC Consolidator Grant (ERC-CoG)
entitled SEDAL with grant agreement]647423.

### Appendix: Convergence bounds for rKEMA

In this Appendix, we study some theoretical properties of the proposedrandomizedKEMA
(rKEMA, page 9) to provide some guarantees of its convergence to KEMA. The solution of
KEMA is the eigensystem of the matrixKd^1 Ks—or alternatively its Tikhonov-regularized

problem (Kd+γI)−^1 Ks. The matrices are now approximated byK^s¼ZsZ>sandK^d¼ZdZ>d,
seeEq (11). Our aim is to give a bound on the approximation error to a product of these two


```
matrices from products of their approximations through random projection matrices. First we
recall the Hermitian Matrix Bernstein theorem, which is then used to derive the bound on
rKEMA.
Theorem 4 (Matrix Bernstein, [ 46 ])LetZ 1 ,...,Zmbe independent n×n random matrices.
Assume thatE½Zi¼ 0 and that the norm of the error matrices is boundedkZikR.Define the
variance parameters^2 ≔maxfk
```
#### P

```
iE½Z
```
```
>
iZik;k
```
#### P

```
iE½ZiZ
```
```
>
ikg.Then, for all t0,
```
#### P

## X

```
i
```
```
Zit
```
(^)
(^)
(^)
(^)

#### !

```
 2 nexp
```
```
t^2
3 s^2 þ 2 Rt
```
#### 

```
and E
```
## X

```
i
```
```
Zi
```
(^)
(^)
(^)
(^)

ffiffiffiffiffiffiffiffiffiffiffiffiffiffiffiffiffiffiffi
3 s^2 logðnÞ
p
þRlogðnÞ:ð 17 Þ
Theorem 5Given two kernel matricesKdandKs,we aim to solve the eigensystem ofKd^1 Ks.
Let us define the corresponding kernel approximationsK^d,K^susing md,msrandom features as
inEq (12),respectively, and m: = min(md,ms).Then, theℓ 2 approximation error bound can be
bounded as
EkK^d^1 K^sKd^1 Ksk
ffiffiffiffiffiffiffiffiffiffiffiffiffiffiffiffiffiffiffiffiffi
3 n^4 logðnÞ
m
r
þ
2 n^2 logðnÞ
m : ð^18 Þ
Proof 1For the sake of simplicity, let us renameD^¼K^d^1 andD¼Kd^1 .We follow a simi-
lar derivation to[ 47 ]for randomized nonlinear CCA. The total error matrix can be decomposed
as a sum of individual error terms,E¼
Pms
i¼ 1 Ei,which are defined asEi¼

#### 1

```
msD^K^
```
```
ðiÞ
s DKs
```
#### 

#### .

```
Now recall that the md+msrandom features are sampled i.i.d. and that the data matrices for
each domain are constant. Therefore, the random matricesfD^ð^1 Þ;...;D^ðmdÞ;K^ðs^1 Þ;...;K^ðsmsÞg
are i.i.d. random variables. Hence, their expectations factorize,E½Ei¼m^1 s E½D^KsDKs
```
#### 

#### ,

```
where we usedE½K^ðsiÞ¼Ks.The deviation of the individual error matrices from their expecta-
tions isZi¼EiE½Ei¼m^1 sD^K^ðsiÞE½D^Ks
```
#### 

```
.Now we can apply Hölder’s condition twice
after using the triangle inequality on the norm, and Jensen’s inequality on the expected values
and obtain a bound of the error matrices,R:
```
```
kZik¼
```
#### 1

```
msk
D^K^ðsiÞE½D^KskBðBþkKskÞ
ms ; ð^19 Þ
```
```
where B is a bound on the norm of the randomized feature map,kzk^2 B.The variance is
defined ass^2 ≔maxfk
```
```
Pms
i¼ 1 E½ZiZ
```
```
>
ik;k
```
```
Pms
i¼ 1 E½Z
```
```
>
iZikg.Let us expand the individual terms
in the (first) summand:
```
```
Z>iZi¼
```
#### 1

```
m^2 s
```
K^ðsiÞD^ (^2) K^ðsiÞþKsE½D^ (^2) KsK^ðsiÞD^E½D^KsE½D^KsD^K^ðsiÞ; ð 20 Þ
and now taking the norm of the expectation, and using Jensen’s inequality, we obtain
kEZ>iZi

#### 	

```
kB
```
(^2) kKsk^2
m^2 ,which is the same forkE½ZiZ
>
ik,and therefore the worst-case estimate of
the variance iss^2 B
(^2) kKsk^2
ms .The bound can be readily obtained by appealing to the matrix Bern-
stein inequality (Theorem 4) and using the fact that random features and kernel evaluations are
upper-bounded by 1, and thus both B andkKkare upper-bounded by n.


```
Theorem 6Equivalently, when we define the corresponding bound for a Tikhonov-regular-
```
ized problem as(Kd+γI)−^1 Ks,and its approximation asðK^dþgIÞ^1 K^s,the bound reduces to

```
EkðK^dþgIÞ^1 K^sðKdþgIÞ^1 Ksk^1
g
```
```
ffiffiffiffiffiffiffiffiffiffiffiffiffiffiffiffiffiffiffiffiffi
3 n^2 logðnÞ
m
```
```
r
þ^2 nlogðnÞ
m
```
#### !

```
; ð 21 Þ
```
whereγ> 0 is a regularization parameter.
Proof 2The demonstration is trivial by following the same rationale and derivations in Theo-

rem 5, and simply boundingkðK^dþgIÞ^1 k^2  1 =g.Interestingly, the bound is exactly the same
as that of the randomized nonlinear CCA in[ 47 ]for the case of paired examples in the domains
and no graph Laplacian terms.
Fig 9shows the absolute error committed by doing an approximation with random features
of the corresponding kernels for rKEMA, along with the derived theoretical bound. We analyze
the issue as a function ofm(here for the sake of simplicity we usedmd=ms=m), and the num-
ber of samplesn. The curves are the result of 300 realizations. The reported results match the
previous bound: we observe a logarithmical trend as a function ofm(linear in the log-scale,
Oðm^1 =^2 Þ), andnlog(n) for the case of training examples, as expected.

### Author Contributions

Conceived and designed the experiments: DT GCV. Performed the experiments: DT. Analyzed
the data: DT. Wrote the paper: DT GCV.

### References

1. Quiñonero-Candela J, Sugiyama M, Schwaighofer A, Lawrence ND. Dataset shift in machine learning.
    Neural information processing series. Cambridge, Mass., London: MIT Press; 2009.
2. Saenko K, Kulis B, Fritz M, Darrell T. Adapting visual category models to new domains. In: Proc.
    ECCV. Berlin, Heidelberg: Springer-Verlag; 2010. p. 213–226.
3. Farhadi A, Tabrizi MK. Learning to Recognize Activities from the Wrong View Point. In: Proc. ECCV.
    Berlin, Heidelberg: Springer-Verlag; 2008. p. 154–166.
4. Duan L, Xu D, Tsang IW, Luo J. Visual Event Recognition in Videos by Learning from Web Data. IEEE
    Trans Pattern Anal Mach Intell. 2012; 34(9):1667–1680. doi:10.1109/TPAMI.2011.265PMID:
    22201057
5. Torralba A, Efros AA. Unbiased look at dataset bias. In: Proc. CVPR. Colorado Springs, CO; 2011. p.
    1521 – 1528.

Fig 9.Error of the approximation,Error≔kK^d^1 K^sKd^1 KskFas a function of the number of samples n (left)
and number of random features m (right) used in the approximation of the eigensystem equation of KEMA.

doi:10.1371/journal.pone.0148655.g009


6. Pan SJ, Yang Q. A Survey on Transfer Learning. IEEE Transactions on Knowledge and Data Engineer-
    ing. 2010 October; 22(10):1345–1359. doi:10.1109/TKDE.2009.191
7. Patel VM, Gopalan R, Li R, Chellappa R. Visual Domain Adaptation: A survey of recent advances.
    IEEE Signal Proc Mag. 2015 May; 32(3):53–69. doi:10.1109/MSP.2014.2347059
8. Jacobs DW, Daume H, Kumar A, Sharma A. Generalized Multiview Analysis: A discriminative latent
    space. In: Proc. CVPR. Providence, RH; 2012. p. 2160–2167.
9. Lai PL, Fyfe C. Kernel and Nonlinear Canonical Correlation Analysis. In: Int. J. Neural Sys.; 2000. p.
    365 – 377. doi:10.1142/S012906570000034X
10. Pan SJ, Yang Q. Domain adaptation via transfer component analysis. IEEE Trans Neural Networks.
2011; 22:199–210. doi:10.1109/TNN.2010.2091281PMID: 21095864
11. Baktashmotlagh M, Harandi MT, Lovell BC, Salzmann M. Domain adaptation on the statistical manifold.
In: Proc. CVPR. Columbus, OH; 2014. p. 2481–2488.
12. Gopalan R, Li R, Chellappa R. Domain adaptation for object recognition: An unsupervised approach.
In: Proc. ICCV. Barcelona, Spain; 2011. p. 999–1006.
13. Gong B, Shi Y, Sha F, Grauman K. Geodesic flow kernel for unsupervised domain adaptation. In: Proc.
CVPR. Providence, RH: IEEE; 2012. p. 2066–2073.
14. Gretton A, Bousquet O, Smola AJ, Schölkopf B. Measuring statistical dependence with Hilbert-Schmidt
norms. In: Jain S, Lee WS, editors. Proc. Algorithmic Learn. Theory; 2005. p. 63–77.
15. Courty N, Flamary R, Tuia D. Domain adaptation with regularized optimal transport. In: Proc. ECML.
Nancy, France; 2014. p. 274–289.
16. Kulis B, Saenko K, Darrell T. What you saw is not what you get: domain adaptation using asymmetric
kernel transforms. In: Proc. CVPR. Colorado Springs, CO; 2011. p. 1785–1792.
17. Jhuo IH, Liu D, Lee DT, Chang SF. Robust visual domain adaptation with low-rank reconstruction. In:
Proc. CVPR. Providence, RH; 2012. p. 2168–2175.
18. Hoffman J, Rodner E, Donahue J, Saenko K, Darrell T. Efficient Learning of Domain Invariant Image
Representations. In: Proc. ICLR. Scottsdale, AZ; 2013.
19. Donahue J, Hoffman J, Rodner E, Saenko K, Darrell T. Semi-supervised Domain Adaptation with
Instance Constraints. In: CVPR; 2013. p. 668–675.
20. Ham J, Lee D, Saul L. Semisupervised alignment of manifolds. In: Cowell RG, Ghahramani Z, editors.
Proc. AISTATS. London, UK; 2005. p. 120–127.
21. Wang C, Krafft P, Mahadevan S. Manifold alignment. In: Ma Y, Fu Y, editors. Manifold Learning: Theory
and Applications. CRC Press; 2011.
22. Hotelling H. Relations Between Two Sets of Variates. Biometrika. 1936 Dec; 28(3/4):321–377. doi:10.
1093/biomet/28.3-4.321
23. Wang C, Mahadevan S. Heterogeneous domain adaptation using manifold alignment. In: IJCAI. Barce-
lona, Spain; 2011. p. 1541–1546.
24. Jolliffe IT. Principal Component Analysis. New York: Springer; 1986.
25. Schölkopf B, Smola AJ, Müller KR. Nonlinear component analysis as a kernel Eigenvalue problem.
Neural Comput. 1998; 10:1299–1319. doi:10.1162/089976698300017467
26. Long M, Wang J, Ding G, Sun J, Yu PS. Transfer Feature Learning with Joint Distribution Adaptation.
In: ICCV; 2013. p. 2200–2207.
27. Tuia D, Muñoz-Marí J, Gómez-Chova L, Malo J. Graph matching for adaptation in remote sensing.
IEEE Trans Geosci Remote Sens. 2013; 51(1):329–341. doi:10.1109/TGRS.2012.2200045
28. Mika S, Schölkopf B, Smola A, Müller KR, Scholz M, Rätsch G. Kernel PCA and De-Noising in Feature
Spaces. In: NIPS 11. MIT Press; 1999. p. 536–542.
29. Bakιr G, Weston J, Schölkopf B. Learning to find Pre-images. In: Proc. NIPS; 2003.
30. Kwok JT, Tsang IW. The Pre-Image Problem in Kernel Methods. IEEE Trans Neural Networks. 2004;
15(6):1517–1525. doi:10.1109/TNN.2004.837781PMID: 15565778
31. Chapelle O, Schölkopf B, Zien A. Semi-Supervised Learning. 1st ed. Cambridge, MA and London,
England: MIT Press; 2006.
32. Camps-Valls G, Bandos T, Zhou D. Semi-supervised Graph-based Hyperspectral Image Classification.
IEEE Transactions on Geoscience and Remote Sensing. 45(10):2044–3054.
33. Tuia D, Volpi M, Trolliet M, Camps-Valls G. Semisupervised Manifold Alignment of Multimodal Remote
Sensing Images. IEEE Trans Geosci Remote Sens. 2014; 52(12):7708–7720. doi:10.1109/TGRS.
2014.2317499


34. Tuia D, Trolliet M, Volpi M. Multisource alignment of image manifolds. In: IEEE International Geosci-
    ence and Remote Sensing Symposium, IGARSS. Melbourne, Australia; 2013.
35. Reed M, Simon B. I: Functional Analysis, Volume 1 (Methods of Modern Mathematical Physics) ( vol 1).
    1st ed. Academic Press; 1981.
36. Riesz F, Nagy BS. Functional Analysis. Frederick Ungar Publishing Co.; 1955.
37. Yan S, Xu D, Zhang B, Zhang HJ, Yang Q, Lin S. Graph Embedding and Extensions: A General Frame-
    work for Dimensionality Reduction. IEEE Trans Patt Anal Mach Intell. 2007; 29(1):40–51. doi:10.1109/
    TPAMI.2007.250598
38. Zhu F, P H, Kallas M. Kernel nonnegative matrix factorization without the pre-image problem. In:
    Machine Learning for Signal Processing. Reims,France; 2014.
39. Rahimi A, Recht B. Random Features for Large-Scale Kernel Machines. In: Neural Information Pro-
    cessing Systems; 2007.
40. Jones LK. Annals of Statistics. A simple lemma on greedy approximation in Hilbert space and conver-
    gence rates for projection pursuit regression and neural network training. 1992; 20:608–613.
41. Shawe-Taylor J, Williams CKI, Cristianini N, Kandola J. On the eigenspectrum of the Gram matrix and
    the generalization error of kernel-PCA. IEEE Trans Info Theory. 2005; 51(7):2510–2522. doi:10.1109/
    TIT.2005.850052
42. Dhanjal C, Gunn SR, Shawe-Taylor J. Efficient Sparse Kernel Feature Extraction Based on Partial
    Least Squares. IEEE Trans Pattern Anal Mach Intell. 2009; 31(8):1347–1361. doi:10.1109/TPAMI.
    2008.171PMID: 19542571
43. Donahue J, Jia Y, Vinyals O, Hoffman J, Zhang N, Tzeng E, et al. DeCAF: A Deep Convolutional Acti-
    vation Feature for Generic Visual Recognition. In: Proceedings of The 31st International Conference on
    Machine Learning; 2014. p. 647–655.
44. Maji S, Berg AC, Malik J. Classification using intersection kernel support vector machines is efficient.
    In: 2008 IEEE Computer Society Conference on Computer Vision and Pattern Recognition (CVPR
    2008), 24–26 June 2008, Anchorage, Alaska, USA; 2008.
45. Song L, Smola A, Gretton A, Borgwardt KM. A Dependence Maximization View of Clustering. In: Proc.
    ICML. Corvallis, OR; 2007. p. 815–822.
46. Mackey L, Jordan MI, Chen RY, Farrell B, Tropp JA. Matrix Concentration Inequalities via the Method
    of Exchangeable Pairs. Annals of Probability. 2014;.
47. Lopez-Paz D, Sra S, Smola AJ, Ghahramani Z, Schölkopf B. Randomized Nonlinear Component Anal-
    ysis. In: Proceedings of the 31 st International Conference on Machine Learning. Beijing, China; 2014.
    p. 1–9.


