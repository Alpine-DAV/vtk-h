FROM alpinedav/ascent-ci:ubuntu-16-cuda-10.1-devel
 
COPY Build.sh /
COPY Test_Both.sh /
COPY config.yaml   /

RUN apt-get update 

RUN apt-get -y install git python gfortran openmpi-bin libopenmpi-dev zlib1g-dev curl
    
