## Experiment 1: Getting to grips with time-dependency in Influenza

Setup is:

- Create ``50`` Influenza (A?) alignments (HA?) with ``50`` sequences each, with varying [time spans](https://github.com/sebastianduchene/virus_time_dep_rates/blob/master/utils/glossary.md);

- Estimate evolutionary rates with [Tempest](http://tree.bio.ed.ac.uk/software/tempest/) and [BEAST](http://beast.bio.ed.ac.uk/), producing a figure similar to figure 4 [here](http://biorxiv.org/content/early/2016/08/15/069492), possibly with credibility intervals for the BEAST estimates.

- For the BEAST analyses, use three different priors: CTMC reference prior, vague Gamma and uniform, generate a figure similar to the above for each;

- For each dataset, sample from each of three priors, generating  the same figure, but for the prior only.
