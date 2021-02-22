FROM debian:stable
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update
RUN apt-get install -q --yes less sudo vim screen
RUN apt-get install -q --yes git
RUN apt-get install --yes gcc
RUN apt-get install --yes gfortran
RUN apt-get install --yes screen git make autoconf automake libtool;
RUN apt-get install --yes openmpi-bin
RUN apt-get install --yes libopenmpi-dev
COPY . psblas3
WORKDIR psblas3
RUN bash autogen.sh
RUN ./configure FCOPT=-O0\ -pipe CCOPT=-O0\ -pipe
RUN make
