# ------------------------------------------------
# Generic Makefile
# ------------------------------------------------

# project name (generate executable with this name)
TARGET    = calibration

CC        = g++
# compiling flags here
CFLAGS    = -g -fopenmp -ggdb -Wall -I$(INCDIR) $(LIBFLAGS)

LINKER    = g++
# linking flags here
LFLAGS    = -ggdb -fopenmp -Wall -I. -lm -lconfig++ \
            -L$(LIBDIR) -l$(LIB)

# change these to proper directories where each file should be
SRCDIR    = src
INCDIR    = inc
OBJDIR    = obj
BINDIR    = bin

# Library flags and folders
LIBFLAGS  = -fpermissive -Ilib/include
LIBDIR    = lib
LIB       = ga

SOURCES   := $(wildcard $(SRCDIR)/*.cpp)
INCLUDES  := $(wildcard $(INCDIR)/*.h)
OBJECTS   := $(SOURCES:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o)
rm        = rm -f


$(BINDIR)/$(TARGET).exe: $(OBJECTS)
	@$(LINKER) $(OBJECTS) $(LFLAGS) -o $@
	@echo "Linking complete!"

$(OBJECTS): $(OBJDIR)/%.o : $(SRCDIR)/%.cpp
	@$(CC) $(CFLAGS) -c $< -o $@
	@echo "Compiled "$<" successfully!"

.PHONY: clean
clean:
	@$(rm) $(OBJECTS)
	@echo "Cleanup complete!"

.PHONY: remove
remove: clean
	@$(rm) $(BINDIR)/$(TARGET).exe
	@echo "Executable removed!"
