FROM memesuite/memesuite:5.5.2 as meme

FROM rocker/r-ubuntu:22.04
COPY --from=meme /opt/meme/libexec/meme-5.5.2 /bin
COPY --from=meme /opt/meme/bin /bin
COPY --from=meme /opt/meme/share/meme-5.5.2/ /opt/meme/share/meme-5.5.2/

ENV MAKEFLAGS=-j4
ENV DEBIAN_FRONTEND=noninteractive
ENV PATH=/usr/meme/bin:/usr/meme/libexec/meme-5.1.1:$PATH
ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8

# Initial package install to get GNUPG
# Add CRAN Repos for R 3.6
RUN apt-get update && \
	apt-get install -yq gnupg ca-certificates && \
	echo "deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/" >> /etc/apt/sources.list && \
	apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 && \
	apt-get update && \
	apt-get upgrade -y &&\
	apt-get install -yq \
	imagemagick python3 python3-pip build-essential \
	wget bzip2 tar zlib1g-dev libexpat1-dev ghostscript \
	aptitude libcurl4-openssl-dev libxml2-dev libopenmpi-dev
	# r-base r-base-dev

# Install DESeq and DESeq2, as well as dependencies
RUN R --slave -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager')" && \
	R --slave -e "BiocManager::install(c('GenomeInfoDb','XML','annotate','latticeExtra','Hmisc','DESeq','DESeq2'))"

# Install bedtools
RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.29.2/bedtools.static.binary && \
	mv bedtools.static.binary /usr/bin/bedtools && \
	chmod a+x /usr/bin/bedtools

# Install samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.10/samtools-1.10.tar.bz2 && \
	bzip2 -d samtools-1.10.tar.bz2 && \
	tar -xf samtools-1.10.tar && \
	cd samtools-1.10 && \
	./configure --prefix=/ && \
	make && \
	make install

# Fix ImageMagick policy to remove errors on conversion
RUN sed -i 's/rights="none" pattern="PDF"/rights="read|write" pattern="PDF"/g' /etc/ImageMagick-6/policy.xml && \
	mv /etc/ImageMagick-6/policy.xml /etc/ImageMagick-6/policy.xmlout

# Install TFEA and Python3 dependencies
# Do this last to cache as much as possible
ADD . /app
RUN pip3 install Cython numpy && pip3 install /app

# Actual run command
CMD	["TFEA"]

# TFEA --test-install
# TFEA --test-full
