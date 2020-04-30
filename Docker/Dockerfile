FROM ubuntu:bionic

ENV MAKEFLAGS=-j4
ENV DEBIAN_FRONTEND=noninteractive
ENV PATH=/usr/meme/bin:/usr/meme/libexec/meme-5.1.1:$PATH
ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8

## Initial package install to get GNUPG
## Add CRAN Repos for R 3.6
RUN apt-get update && \
	apt-get install -yq gnupg ca-certificates && \
	echo "deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/" >> /etc/apt/sources.list && \
	apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 && \
	apt-get update && \
	apt-get install -yq \
	imagemagick python3 python3-pip build-essential \
	wget bzip2 tar zlib1g-dev libexpat1-dev ghostscript \
	aptitude libcurl4-openssl-dev libxml2-dev \
	r-base r-base-dev

## Install TFEA and Python3 dependencies
RUN pip3 install Cython numpy && pip3 install tfea

## Install DESeq and DESeq2, as well as dependencies
RUN R --slave -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager')" && \
	R --slave -e "BiocManager::install(c('GenomeInfoDb','XML','annotate','latticeExtra','Hmisc','DESeq','DESeq2'))"

## Install bedtools
RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.29.2/bedtools.static.binary && \
	mv bedtools.static.binary /usr/bin/bedtools && \
	chmod a+x /usr/bin/bedtools

## Install samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.10/samtools-1.10.tar.bz2 && \
	bzip2 -d samtools-1.10.tar.bz2 && \
	tar -xf samtools-1.10.tar && \
	cd samtools-1.10 && \
	./configure --prefix=/ && \
	make && \
	make install

# Install MEME Perl Dependencies
RUN cpan \
	XML::Simple HTML::Treebuilder JSON File::Which Log::Log4perl \
	Math::CDF XML::Compile::SOAP11 XML::Compile::WSDL11

## Install MEME Suite
RUN wget http://meme-suite.org/meme-software/5.1.1/meme-5.1.1.tar.gz && \
	tar zxf meme-5.1.1.tar.gz && \
	cd meme-5.1.1 && \
	./configure --prefix=/usr/meme/ --with-url=http://meme-suite.org/ --enable-build-libxml2 --enable-build-libxslt && \
	make && \
	make install && \
	perl scripts/dependencies.pl

CMD	["TFEA"]

# TFEA --test-install
# TFEA --test-full
