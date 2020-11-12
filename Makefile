# ------------------------------------------------
# Generic Makefile
# ------------------------------------------------

# project name (generate executable with this name)
TARGET    = calibration

CC        = g++
# compiling flags here
CFLAGS    = -g -fopenmp -ggdb -w -Wall -I$(INCDIR) $(LIBFLAGS)

LINKER    = g++
# linking flags here
LFLAGS    = -ggdb -fopenmp -Wall -I. -lm \
            -L$(LIBDIR) $(LIB)

# change these to proper directories where each file should be
SRCDIR    = src
INCDIR    = inc
OBJDIR    = obj
BINDIR    = bin

# Library flags and folders
LIBFLAGS  = -fpermissive -Istatic_lib/include
LIBDIR    = static_lib/lib
LIB       = -lga -lconfig++

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
