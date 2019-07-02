# robust-autoreg readme

This code for the robust autoregressive occupancy model of Hall et al. 2019, "Validating dispersal distances inferred from occupancy models with genetic parentage assignments." We developed a novel rigorous auto-logistic occupancy model in order to create and test the performance of connectivity metrics that relate colonization of sites in a metapopulation to measures of site connectivity (or it’s inverse, isolation) that account for two forms of uncertainty: (1) the occupancy status of sites that have not been surveyed, and (2) the chance that sites where the species of interest was not detected could have been occupied by the species. We validated the performance of these metrics with genetic parent-offspring data. 

Our model is a modified version of MacKenzie et al. (2003)’s single-species multi-season occupancy model, which does not assume dynamic equilibrium and analyzes binary detection histories with revisits (i.e., 11 00 01 for a patch where the species was detected on both visits in the season one, not detected in the season two, and detected on the second visit only in season three). By assuming closure to changes in occupancy within a season and no false positive detections, these revisits allow for the estimation of probability of detection per visit j (ρ<sub>i,t,j</sub>), and detection-corrected estimates of patch occupancy in the first year (ψ<sub>i,1</sub>) and patch colonization (γ<sub>i,t</sub>) and extinction (ε<sub>i,t</sub>) in subsequent years (for patch i in year t). These four parameters can be used to derive the conditional occupancy probability in each year (ψ<sup>c</sup><sub>i,t</sub>, the probability that a patch was occupied conditional on its detection history), which can then act as an auto-covariate to model dispersal’s influence on metapopulation dynamics, for example, colonization:

logit(γ<sub>i,t</sub>) = β<sub>0</sub> + β<sub>1</sub>×X<sub>i,t</sub> + ⋯ + β<sub>n</sub>*∑<sub>p≠i</sub>(w<sub>i,p</sub>ψ<sup>c</sup><sub>p,t-1</sub>)

where w<sub>i,p</sub> is a weighting function between focal patch i and other patches p (in this study, a buffer radius or incidence function; see below). This auto-covariate builds on those implemented in previous auto-logistic occupancy models (Yackulic et al. 2012, Eaton et al. 2014), but is more general by allowing for weights that are specific to both the focal patch and colonist source patch (i.e., the incidence function model or asymmetric resistance surfaces). Our model also improves on these previous efforts via a more rigorous calculation of conditional occupancy probabilities. 

Unconditional occupancy in seasons after the first can be calculated as a derived parameter as the probability that the site was previous occupied and did not go extinct, or was previously unoccupied and was colonized:

ψ<sub>i,t</sub> = ψ<sub>i,t-1</sub> (1-ε<sub>i,t</sub>) + (1-ψ<sub>i,t-1</sub>)γ<sub>i,t</sub>

The probability of unconditional occupancy in each season can thus be estimated on the basis of site-level and annual covariates that inform ψ<sub>i,1</sub>, γ<sub>i,t</sub>, ε<sub>i,t</sub>, which are estimated from a comprehensive analysis of sites that were surveyed. When combined with p<sub>i,t,j</sub> these four parameters can be used to derive the conditional occupancy probability in each season, c.ψ<sub>i,t</sub>, the probability that a surveyed patch was occupied conditional on its detection history. For sites that were surveyed and a species was detected, ψ<sup>c</sup><sub>i,t</sub> = 1 because the model assumes no false positive detections.

The calculation of conditional occupancy for a single-species single-season occupancy model (Royal & Dorazio 2008, page 107; MacKenzie et al. 2006, pages 97-98), using Bayes Law, is equal to the probability that the site was occupied and not detected divided by the probability of the species not being detected (which includes the former probability plus the probability the species was not detected because the site was not occupied):  

ψ<sup>c</sup><sub>i,t</sub> = Pr(species present & not detected) / Pr(species not detected)					

which MacKenzie et al. (2006) solve for single-season occupancy models (where K is the total number of visits in that season):

ψ<sup>c</sup><sub>i,t</sub> = (ψ<sub>i,t</sub>∏<sub>j = 1 to K</sub>(1-p<sub>i,t,j</sub>)) / ((1-ψ<sub>i,t</sub>) + ψ<sub>i,t</sub>∏<sub>j=1 to K</sub>(1-p<sub>i,t,j</sub>))			

Previously developed autoregressive occupancy models (Yackulic et al. 2012, Eaton et al. 2014) took a computationally efficient approach to estimating conditional occupancy by using only information about the occupancy state in the current season and previous season. They used the above equation to calculate conditional occupancy in the first season, and used the first equation above to calculate ψ<sub>i,t</sub> when calculating ψ<sup>c</sup><sub>i,t</sub> in all subsequent seasons. For example, given the detection history y = <sub>11 00 00</sub>, their conditional probability of occupancy in the second season would be (omitting subscripts):

ψ<sup>c</sup><sub>i,2</sub> = (ψ(1-ε) + (1-ψ)γ)(1-p)<sup>2</sup>) / (1-(ψ(1-ε) + (1-ψ)γ) + (ψ(1-ε) + (1-ψ)γ)(1-p)<sup>2</sup>)		
			     
While this estimation method is unbiased for unsurveyed sites and the single-season occupancy models it was derived from, it is biased for multi-season models because the occupancy probability at time t + 1 is informed not only by our knowledge of the occupancy probability in the previous time step (t), but also in the next time steps (t + 2, 3, … T). For example, in a metapopulation with low colonization probability the likelihood that a species that was not detected in its year two of three was actually present is higher with detection history 11 00 11 than with detection history 11 00 00. Knowing that the bird was detected in the next year increases the likelihood that they were present while undetected, rather than a rare re-colonization event, while in the latter history, additional consecutive non-detections increase our confidence that it was a true extinction. Thus, while colonization and extinction operate as first-order Markov processes, conditional occupancy probability is a function of the entire detection history such that. We expand this to include additional survey seasons, changing our equation to a more general:

ψ<sup>c</sup><sub>i,t</sub> = Pr(occupancy histories with patch i occupied in season t | detection history) / Pr(all occupancy histories for patch i | detection history) 

We implemented this equation to include a site’s entire detection history to inform the estimation of ψ<sup>c</sup><sub>i,t</sub>. For the above example, and where the patch’s “occupancy history” is a vector of the unknown true occupancy states in each season, z, our formulation would yield for the second season in the above example:

ψ<sup>c</sup><sub>i,2</sub> = Pr((z=1 1 0|y) +  Pr(z=1 1 1|y)) / (Pr(z=1 0 0|y)  + Pr(z=1 0 1|y) + Pr(z=1 1 0|y) + Pr(z=1 1 1|y))

   = (ψp<sup>2</sup>(1-ε)(1-p)<sup>2</sup>ε + ψp<sup>2</sup>(1-ε)<sup>2</sup>(1-p)<sup>4</sup>) / (ψp<sup>2</sup>ε(1-γ) + ψ<sup>2</sup>εγ(1-p)<sup>2</sup> + p<sup>2</sup>(1-ε)(1-p)<sup>2</sup>ε + ψ<sup>2</sup>(1-ε)<sup>2</sup>(1-p)<sup>4</sup>)

omitting subscripts and occupancy histories with zero probability (i.e., those with the species absent in seasons when it was detected). For seasons when the species was detected, the numerator and denominator become the same and thus ψ<sup>c</sup><sub>i,t</sub> = 1. For seasons when a patch was unsurveyed, p can simply be omitted from that season in the function, as in MacKenzie et al (2003). For sites that were never surveyed, occupancy cannot be conditioned on visits and instead unconditional occupancy (the first equation above) is used.

Finally, these estimates of conditional occupancy can be used for calculating the connectivity auto-covariates on occupancy, colonization, and extinction rates. This code includes two options, the first of which is the buffer radius measure (BRM) that sums the area of all occupied sites within a radius of dispersal (Moilanen & Nieminen 2002):

BRM<sub>i,t</sub> = ∑<sub>s≠i</sub>(ψ<sup>c</sup><sub>s,t-1</sub>A<sub>s</sub>B<sub>s</sub>); if d<sub>i,s</sub> ≤ r then B = 1, else B = 0

where A<sub>s</sub> is the area of the colonist source patch s, d<sub>i,s</sub> is the distance between the focal patch centroid and colonist source patch centroid, and r is a user-specified constant buffer radius. The second autoregressive measure was the incidence function measure (IFM), a function representing the likelihood of successful dispersal decaying continuously over distance (Moilanen & Nieminen 2002):

IFM<sub>i,t</sub> = ∑<sub>s≠i</sub>(Ψ<sub>s,t-1</sub>A<sub>s</sub>e^(-∝d<sub>i,s</sub>))

where and 1/∝ is the mean dispersal distance assuming a negative exponential kernel. These conditional occupancy estimates can then be fit into auto-covariates and models iteratively re-fitted until convergence is achieved (Augustin et al. 1996).

This model is estimated via pseudo-maximum likelihood in a combination of [Program PRESENCE](https://www.mbr-pwrc.usgs.gov/software/presence.html) (Hines 2006) for fitting multi-species occupancy models and R (this code) for calculating auto-covariates.


## References

Augustin, N.H., Mugglestone, M.A. & Buckland, S.T. 1996. An autologistic model for the spatial distribution of wildlife. The Journal of Applied Ecology 33:339-347

Eaton, M. J., P. T. Hughes, J. E. Hines, and J. D. Nichols. 2014. Testing metapopulation concepts: effects of patch characteristics and neighborhood occupancy on the dynamics of an endangered lagomorph. Oikos 123:662–676.

Hines, J.E. 2013. PRESENCE - Software to estimate patch occupancy and related parameters. U.S.G.S. Patuxent Wildlife Research Center, Laurel, MD, USA.

MacKenzie, D. I., J. D. Nichols, J. E. Hines, M. G. Knutson, and A. B. Franklin. 2003. Estimating site occupancy, colonization, and local extinction when a species is detected imperfectly. Ecology 84:2200–2207.

MacKenzie, D. I., J. D. Nichols, J. A. Royle, K. H. Pollock, L. L. Bailey, and J. E. Hines. 2006. Occupancy estimation and modeling: Inferring patterns and dynamics of species occurrence. Academic Press, Burlington, MA.

Moilanen, A., & Nieminen, M. 2002. Simple connectivity measures in spatial ecology. Ecology 83:1131–1145.

Yackulic, C. B., J. Reid, R. Davis, J. E. Hines, J. D. Nichols, and E. Forsman. 2012. Neighborhood and habitat effects on vital rates: expansion of the Barred Owl in the Oregon Coast Ranges. Ecology 93:1953–1966.
