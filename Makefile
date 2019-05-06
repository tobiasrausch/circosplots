SHELL := /bin/bash

# Targets
TARGETS = .circos circos.png
PBASE=$(shell pwd)

all:   	$(TARGETS)

.circos:
	wget 'http://circos.ca/distribution/circos-0.69-6.tgz' && tar -xzf circos-0.69-6.tgz && rm circos-0.69-6.tgz && touch .circos

%.png: %.conf
	${PBASE}/circos-0.69-6/bin/circos -conf $<

clean:
	rm -rf $(TARGETS)
