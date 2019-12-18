FROM rocker/tidyverse:3.6.1

RUN apt install libffi-dev
RUN Rscript -e 'install.packages("synapser", repos=c("http://ran.synapse.org", "http://cran.fhcrc.org"))'
RUN Rscript -e 'install.packages("optparse")'

# This retriggers the cache when building with --no-cache
ARG CACHE_DATE=not_a_date

RUN git clone https://github.com/Sage-Bionetworks/beat-pd.git /beat-pd
