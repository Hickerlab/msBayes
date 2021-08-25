FROM r-base:latest
RUN apt-get update
RUN apt-get install -y libgsl0-dev gsl-bin make
RUN mkdir /msbayes
COPY * /msbayes/
WORKDIR /msbayes
RUN Rscript r/dependencies.r
RUN make
RUN make install
RUN cp /root/bin/* /usr/bin/
RUN cp /root/lib/* /usr/lib/
RUN mkdir /workspace
WORKDIR /workspace
CMD ["msbayes.pl", "-h"]