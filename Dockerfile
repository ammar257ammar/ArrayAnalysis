FROM        rocker/shiny:3.6.3
MAINTAINER  ammar257ammar@gmail.com

USER root

RUN apt-get update && apt-get install -y \
    sudo \
    pandoc \
    pandoc-citeproc \
    libcurl4-gnutls-dev \
    libcairo2-dev \
    libxt-dev \
    libssl-dev \
    libssh2-1-dev

WORKDIR /shiny 

#COPY ./install_packages.R  ./ 

#RUN Rscript ./install_packages.R

COPY app/server.R ./
COPY app/ui.R ./

COPY ./shiny-server.conf  /etc/shiny-server/shiny-server.conf

USER root
RUN usermod -a -G staff root

EXPOSE 7730
EXPOSE 3838

CMD ["/usr/bin/shiny-server.sh"]
