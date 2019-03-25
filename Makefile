###
###     Makefile for JM H.264/AVC encoder/decoder
###
###             generated for UNIX/LINUX environments
###             by Limin Wang 
###

SUBDIRS := lencod ldecod rtpdump

.PHONY: default all clean tags depend $(SUBDIRS)

default: all

all: $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -C $@

clean tags depend:
	for i in $(SUBDIRS); do make -C $$i $@; done

