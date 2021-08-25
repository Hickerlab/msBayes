FROM r-base:latest
RUN apt-get update
RUN apt-get install -y libgsl0-dev gsl-bin make
RUN mkdir /msbayes
COPY * /msbayes/
WORKDIR /msbayes
RUN make
RUN make install
RUN cp /root/bin/* /usr/bin/
CMD ["msbayes.pl", "-h"]