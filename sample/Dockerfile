FROM selpo/qboot

WORKDIR /srv/sample
ADD CMakeLists.txt main.cpp ./

RUN mkdir build \
	&& cd build \
	&& cmake .. -DCMAKE_BUILD_TYPE=Debug \
	&& make -j4

CMD [ "/bin/bash" ]
