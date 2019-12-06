FROM rocker/tidyverse:3.6.1

RUN apt install libffi-dev
RUN Rscript -e 'install.packages("synapser", repos=c("http://ran.synapse.org", "http://cran.fhcrc.org"))'
RUN Rscript -e 'install.packages("optparse")'
RUN git clone https://github.com/philerooski/pddb2.git /root/beat-pd
RUN cd /root/beat-pd && git checkout phil-changes

COPY /root/beat-pd/scoring_code/BEAT-PD_Scoring_Code.R /usr/local/bin/
CMD BEAT-PD_Scoring_Code.R
