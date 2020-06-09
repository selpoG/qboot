FROM gcc:latest
WORKDIR /srv

# apt
ENV DEBCONF_NOWARNINGS yes
RUN apt-get update \
	&& apt-get install -y --no-install-recommends \
		cmake ca-certificates lzip

# build gmp
RUN curl https://gmplib.org/download/gmp/gmp-6.2.0.tar.lz | tar --lzip -xf - \
	&& cd gmp-6.2.0 \
	&& ./configure --enable-cxx \
	&& make -j4 \
	&& make install \
	&& ldconfig \
	&& cd .. \
	&& rm -r gmp-6.2.0

# build mpfr
RUN curl https://www.mpfr.org/mpfr-current/mpfr-4.0.2.tar.gz | tar zxf - \
	&& cd mpfr-4.0.2 \
	&& ./configure \
	&& make -j4 \
	&& make install \
	&& ldconfig \
	&& cd .. \
	&& rm -r mpfr-4.0.2

# build qboot
RUN git clone https://github.com/selpoG/qboot \
	&& mkdir qboot/build \
	&& cd qboot/build \
	&& cmake .. -DCMAKE_BUILD_TYPE=Debug \
	&& make -j4 \
	&& make install \
	&& cmake . -DCMAKE_BUILD_TYPE=Release \
	&& make -j4 \
	&& make install \
	&& cd /srv \
	&& rm -r qboot

WORKDIR /
CMD [ "/bin/bash" ]
