FROM rocker/tidyverse:3.6.1

RUN apt install libffi-dev
RUN Rscript -e 'install.packages("synapser", repos=c("http://ran.synapse.org", "http://cran.fhcrc.org"))'
RUN git clone https://github.com/philerooski/pddb2.git /root/beat-pd
RUN cd /root/beat-pd && git checkout phil-changes

CMD Rscript /root/beat-pd/scoring_code/BEAT-PD_Scoring_Code.R
