FROM continuumio/miniconda
MAINTAINER Nick Swainston <nickaswainston@gmail.com>

RUN apt-get -y install ttf-dejavu

COPY conda.yml .
RUN \
   conda env update -n root -f conda.yml \
&& conda clean -a

#grab my accel_sift from my repo
RUN mkdir /tmp/build && \
    cd /tmp/build && \
    git clone https://github.com/NickSwainston/mwa_search.git && \
    chmod ugo+x mwa_search/scripts/ACCEL_sift.py && \
    mv mwa_search/scripts/ACCEL_sift.py /usr/bin/

RUN git clone https://github.com/scottransom/presto.git && \
    chmod ugo+x presto/lib/python/*py && \
    mv presto/lib/python/*py /usr/bin/
    
