###
###     Makefile for H.264/AVC RTP dump tool
###
###             generated for UNIX/LINUX environments
###             by H. Schwarz, Limin Wang
###



NAME=   rtpdump

### include debug information: 1=yes, 0=no
DBG?= 0

DEPEND= dependencies

BINDIR= ../bin
INCDIR= .
SRCDIR= .
OBJDIR= obj

CC=     $(shell which gcc)

LIBS=   -lm
FLAGS=  -ffloat-store -Wall -I$(INCDIR) -D __USE_LARGEFILE64 -D _FILE_OFFSET_BITS=64

ifeq ($(DBG),1)
SUFFIX= .dbg
FLAGS+= -g
else
SUFFIX=
FLAGS+= -O2
endif

OBJSUF= .o$(SUFFIX)

SRC=    $(wildcard $(SRCDIR)/*.cpp) 
OBJ=    $(SRC:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o$(SUFFIX)) 
BIN=    $(BINDIR)/$(NAME)$(SUFFIX).exe

.PHONY: default distclean clean tags depend

default: objdir_mk depend bin

clean:
	@echo remove all objects
	@rm -rf $(OBJDIR)

distclean: clean
	@rm -f $(DEPEND) tags
	@rm -f $(BIN)

tags:
	@echo update tag table
	@ctags *.h *.cpp

bin:    $(OBJ)
	@echo
	@echo 'creating binary "$(BIN)"'
	@$(CXX) -o $(BIN) $(OBJ) $(LIBS)
	@echo '... done'
	@echo

depend:
	@echo
	@echo 'checking dependencies'
	@$(SHELL) -ec '$(CXX) -MM $(CFLAGS) -I$(INCDIR) $(SRC)                  \
         | sed '\''s@\(.*\)\.o[ :]@$(OBJDIR)/\1.o$(SUFFIX):@g'\''               \
         >$(DEPEND)'
	@echo

$(OBJDIR)/%.o$(SUFFIX): $(SRCDIR)/%.cpp
	@echo 'compiling object file "$@" ...'
	@$(CC) -c -o $@ $(FLAGS) $<

objdir_mk:
	@echo 'Creating $(OBJDIR) ...'
	@mkdir -p $(OBJDIR)


-include $(DEPEND)

