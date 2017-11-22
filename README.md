# ReVEAL
This repository provides Matlab codes to recover accelerated phase contrast MRI data using <a href= "https://www.ncbi.nlm.nih.gov/pubmed/26444911">ReVEAL</a>.  Codes are provided to recover planar images with 1-3 velocity endoding dimensions, and 4D flow data.  Two example data sets are provided 

<h4> Description </h4>
To regularize the ill-posed inverse problem, we present a novel inversion algorithm that applies regularization based on structure unique to phase contrast MRI data. Adopting an empirical Bayes approach, spatial and temporal redundancies are exploited via sparsity in the wavelet domain, and the voxel-wise magnitudeand phase structure across encodings is captured in a conditional mixture prior that applies regularizing constraints based on the presence of flow.

<h4> Installation  </h4>
To install the codes, simply download the ReVEAL repository and add the folder and its subfolders to the Matlab path.  This repository also requires the GAMPLAB software package available at <a href="http://sourceforge.net/projects/gampmatlab/files/">http://sourceforge.net/projects/gampmatlab/files/</a>.  To generate additional VISTA sampling patterns, download the <a href="https://github.com/OSU-CMR/VISTA">VISTA</a> repository. 

<h4> Publications </h4>
<p>
1) <a href = "https://www.ncbi.nlm.nih.gov/pubmed/26444911">A Bayesian model for highly accelerated phase-contrast MRI</a><br/>
2)<br/>
3) <a href = "https://www.ncbi.nlm.nih.gov/pubmed/28270219">Quantification of aortic stenosis diagnostic parameters: comparison of fast 3 direction and 1 direction phase contrast CMR and transthoracic echocardiography</a>
</p>
