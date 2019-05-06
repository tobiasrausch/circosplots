SHELL := /bin/bash

# Targets
TARGETS = circos.png
PBASE=$(shell pwd)

all: .circos $(TARGETS)

.circos:
	wget 'http://circos.ca/distribution/circos-0.69-6.tgz' && tar -xzf circos-0.69-6.tgz && rm circos-0.69-6.tgz && touch .circos

%.png: %.conf
	${PBASE}/circos-0.69-6/bin/circos -conf $<

clean:
	rm -rf $(TARGETS)

distclean:
	rm -rf $(TARGETS) .circos circos-0.69-6/

.PHONY: distclean clean all
