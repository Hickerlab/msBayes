FROM r-base:latest
RUN apt-get update
RUN apt-get install -y libgsl0-dev gsl-bin make
RUN mkdir /msbayes
COPY * /msbayes/
COPY src/*.r /usr/bin/
WORKDIR /msbayes
RUN make
RUN make install
RUN cp /root/bin/* /usr/bin/
RUN mkdir /workspace
WORKDIR /workspace
RUN Rscript /msbayes/r/dependencies.r
CMD ["msbayes.pl", "-h"]