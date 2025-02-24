################################################################################
# HOWTO use this file for local builds                                         #
################################################################################
# To build Khoca, run this command from your khoca/ directory:                 #
# $ docker build -f Dockerfile --tag khoca:<your_tag> .                        #
# To run Khoca:                                                                #
# $ docker run -it khoca:<your_tag>                                            #
################################################################################
FROM ubuntu:latest
LABEL maintainer="Sebastian Oehms <seb.oehms@gmail.com>"
# Set sane defaults for common environment variables.
ENV LC_ALL C.UTF-8
ENV LANG C.UTF-8
ENV SHELL /bin/bash
# install prerequisites
RUN apt-get -qq update \
    && apt-get -qq install -y --no-install-recommends git g++ make libgmp-dev pari-gp2c python3 python3-dev cython python-is-python3\
    && apt-get -qq clean \
    && rm -r /var/lib/apt/lists/* \
    && git clone https://github.com/soehms/khoca.git \
    && cd khoca/ \
    && make
WORKDIR "/khoca/"
COPY ./docker/entrypoint.sh /usr/local/bin/khoca-entrypoint
ENTRYPOINT ["/usr/local/bin/khoca-entrypoint"]
CMD ["bash"]
